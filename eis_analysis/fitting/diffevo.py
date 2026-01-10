"""
Differential Evolution global optimization for circuit fitting.

Clean design: No logging in core functions, all diagnostics returned as data.

Strategy options:
    1 = 'randtobest1bin' (default) - balanced exploration/exploitation
    2 = 'best1bin' - faster convergence, may miss global minimum
    3 = 'rand1bin' - better exploration, slower convergence
"""

import numpy as np
import logging
import warnings
from typing import Tuple, List, Optional
from numpy.typing import NDArray
from dataclasses import dataclass, field
from scipy.optimize import differential_evolution, least_squares, OptimizeWarning

from .circuit import FitResult, FitDiagnostics, Circuit
from .bounds import generate_simple_bounds
from .covariance import compute_covariance_matrix
from .diagnostics import compute_weights, compute_fit_metrics
from .jacobian import make_jacobian_function

logger = logging.getLogger(__name__)

# Strategy mapping: number -> scipy strategy name
DE_STRATEGIES = {
    1: 'randtobest1bin',
    2: 'best1bin',
    3: 'rand1bin',
}


class _DECostFunction:
    """Picklable cost function for differential evolution with workers > 1."""

    def __init__(self, circuit: Circuit, frequencies: NDArray[np.float64],
                 Z: NDArray[np.complex128], weights: NDArray[np.float64],
                 fixed_params: list = None, full_initial_guess: list = None):
        self.circuit = circuit
        self.frequencies = frequencies
        self.Z = Z
        self.weights = weights
        self.fixed_params = fixed_params
        self.full_initial_guess = full_initial_guess

    def _reconstruct_params(self, free_params):
        if self.fixed_params is None or not any(self.fixed_params):
            return list(free_params)
        full, idx = [], 0
        for i, is_fixed in enumerate(self.fixed_params):
            if is_fixed:
                full.append(self.full_initial_guess[i])
            else:
                full.append(free_params[idx])
                idx += 1
        return full

    def __call__(self, params):
        full_params = self._reconstruct_params(params)
        Z_pred = self.circuit.impedance(self.frequencies, full_params)
        residuals_real = (self.Z.real - Z_pred.real) * self.weights
        residuals_imag = (self.Z.imag - Z_pred.imag) * self.weights
        return np.sum(residuals_real**2 + residuals_imag**2)


@dataclass
class DiffEvoDiagnostics:
    """Diagnostics from differential evolution optimization."""
    # DE settings
    strategy: str
    popsize: int
    maxiter: int
    tol: float
    workers: int
    weighting: str
    jacobian_type: str

    # DE results
    de_converged: bool
    de_iterations: int
    de_evaluations: int
    de_error: float

    # Refinement results
    refined_error: float
    refinement_improved: bool
    total_evaluations: int

    # Fixed params info
    n_fixed_params: int = 0
    fixed_param_indices: List[int] = field(default_factory=list)

    # Warnings
    warnings: List[str] = field(default_factory=list)


@dataclass
class DiffEvoResult:
    """
    Result from differential evolution optimization.

    Attributes
    ----------
    best_result : FitResult
        Best fitting result after least_squares refinement
    de_result : scipy.optimize.OptimizeResult
        Raw result from differential_evolution
    de_error : float
        Fit error after DE (before refinement) [%]
    final_error : float
        Fit error after least_squares refinement [%]
    n_evaluations : int
        Number of function evaluations
    strategy : str
        DE strategy used
    improvement : float
        Relative improvement from DE to refined [%]
    diagnostics : DiffEvoDiagnostics
        Detailed diagnostics
    """
    best_result: FitResult
    de_result: any
    de_error: float
    final_error: float
    n_evaluations: int
    strategy: str
    improvement: float
    diagnostics: Optional[DiffEvoDiagnostics] = None


def fit_circuit_diffevo(
    circuit: Circuit,
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    strategy: int = 1,
    popsize: int = 15,
    maxiter: int = 1000,
    tol: float = 0.01,
    workers: int = 1,
    weighting: str = 'modulus',
    verbose: bool = True,  # Kept for backward compatibility, ignored
    use_analytic_jacobian: bool = True
) -> Tuple[DiffEvoResult, NDArray[np.complex128], any]:
    """
    Fit circuit using Differential Evolution global optimization.

    Parameters
    ----------
    circuit : Circuit
        Circuit object with initial parameter guesses
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance data [Ohm]
    strategy : int, optional
        DE strategy: 1='randtobest1bin', 2='best1bin', 3='rand1bin'
    popsize : int, optional
        Population size multiplier (default: 15)
    maxiter : int, optional
        Maximum number of generations (default: 1000)
    tol : float, optional
        Relative tolerance for convergence (default: 0.01)
    workers : int, optional
        Number of parallel workers (default: 1)
    weighting : str, optional
        Weighting scheme
    verbose : bool, optional
        Ignored (kept for backward compatibility)
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian for refinement (default: True)

    Returns
    -------
    diffevo_result : DiffEvoResult
        Differential evolution result with all diagnostics
    Z_fit : ndarray
        Best fit impedance
    fig : matplotlib.figure.Figure
        Nyquist plot with best fit
    """
    strategy_name = DE_STRATEGIES.get(strategy, 'randtobest1bin')
    diag_warnings = []

    # Get initial guess from circuit definition
    initial_guess_full = list(circuit.get_all_params())

    # Get parameter labels and bounds
    param_labels = None
    if hasattr(circuit, 'get_param_labels'):
        param_labels = circuit.get_param_labels()
        lower_bounds_full, upper_bounds_full = generate_simple_bounds(param_labels)
    else:
        n_params_full = len(initial_guess_full)
        lower_bounds_full = [1e-15] * n_params_full
        upper_bounds_full = [1e15] * n_params_full

    # Get fixed parameters
    fixed_params = None
    fixed_param_indices = []
    if hasattr(circuit, 'get_all_fixed_params'):
        fixed_params = circuit.get_all_fixed_params()
        fixed_param_indices = [i for i, f in enumerate(fixed_params) if f]

    # Filter to free parameters only
    if fixed_params is not None and any(fixed_params):
        initial_guess = [v for v, f in zip(initial_guess_full, fixed_params) if not f]
        lower_bounds = [lb for lb, f in zip(lower_bounds_full, fixed_params) if not f]
        upper_bounds = [ub for ub, f in zip(upper_bounds_full, fixed_params) if not f]
    else:
        initial_guess = initial_guess_full
        lower_bounds = lower_bounds_full
        upper_bounds = upper_bounds_full

    initial_guess = np.array(initial_guess)
    bounds = list(zip(lower_bounds, upper_bounds))
    n_params = len(bounds)

    # Clip initial guess to bounds
    initial_guess = np.clip(initial_guess, lower_bounds, upper_bounds)

    # Precompute weights
    weights = compute_weights(Z, weighting)

    # Cost function for DE
    cost_function = _DECostFunction(
        circuit, frequencies, Z, weights,
        fixed_params=fixed_params,
        full_initial_guess=initial_guess_full
    )

    # Helper to reconstruct full params from free params
    def reconstruct_params(free_params):
        if fixed_params is None or not any(fixed_params):
            return list(free_params)
        full, idx = [], 0
        for i, is_fixed in enumerate(fixed_params):
            if is_fixed:
                full.append(initial_guess_full[i])
            else:
                full.append(free_params[idx])
                idx += 1
        return full

    # Residual function for least_squares
    def residual_function(params):
        full_params = reconstruct_params(params)
        Z_pred = circuit.impedance(frequencies, full_params)
        return np.concatenate([
            (Z.real - Z_pred.real) * weights,
            (Z.imag - Z_pred.imag) * weights
        ])

    # Step 1: Run Differential Evolution
    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            de_result = differential_evolution(
                cost_function,
                bounds,
                x0=initial_guess,
                strategy=strategy_name,
                popsize=popsize,
                maxiter=maxiter,
                tol=tol,
                workers=workers,
                polish=False,
                seed=None,
                disp=False,
                updating='deferred' if workers != 1 else 'immediate',
            )
    except Exception as e:
        raise RuntimeError(f"DE optimization failed: {e}") from e

    # Compute DE error
    de_params_full = reconstruct_params(de_result.x)
    Z_fit_de = circuit.impedance(frequencies, de_params_full)
    de_error_rel, _, _ = compute_fit_metrics(Z, Z_fit_de, weighting)

    # Step 2: Refine with least_squares
    jacobian_type = 'analytic'
    if use_analytic_jacobian:
        try:
            jac_func = make_jacobian_function(
                circuit, frequencies, weights,
                fixed_params=fixed_params,
                full_initial_guess=initial_guess_full
            )
        except NotImplementedError:
            jac_func = '2-point'
            jacobian_type = 'numeric'
    else:
        jac_func = '2-point'
        jacobian_type = 'numeric'

    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always", OptimizeWarning)

            x_scale = np.maximum(np.abs(de_result.x), 1e-10)
            ls_result = least_squares(
                residual_function,
                de_result.x,
                jac=jac_func,
                bounds=(lower_bounds, upper_bounds),
                method='trf',
                x_scale=x_scale,
                ftol=1e-10,
                xtol=1e-10,
                gtol=1e-10,
                max_nfev=5000,
            )
    except Exception as e:
        diag_warnings.append(f"Refinement failed: {e}, using DE result")
        ls_result = de_result
        ls_result.x = de_result.x

    # Reconstruct full params
    ls_params_free = np.array(ls_result.x)
    ls_params_full = np.array(reconstruct_params(ls_params_free))

    # Compute refined fit error
    Z_fit_ls = circuit.impedance(frequencies, list(ls_params_full))
    ls_error_rel, _, _ = compute_fit_metrics(Z, Z_fit_ls, weighting)

    improvement = (de_error_rel - ls_error_rel) / de_error_rel * 100 if de_error_rel > 0 else 0

    # Choose better result
    if ls_error_rel <= de_error_rel:
        params_opt_free = ls_params_free
        params_opt = ls_params_full
        Z_fit = Z_fit_ls
        fit_error_rel = ls_error_rel
        used_refinement = True
    else:
        params_opt_free = np.array(de_result.x)
        params_opt = np.array(de_params_full)
        Z_fit = Z_fit_de
        fit_error_rel = de_error_rel
        used_refinement = False
        diag_warnings.append("Refinement worsened fit, using DE result")

    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

    # Step 3: Compute covariance
    final_residuals = residual_function(params_opt_free)
    n_params_full = len(params_opt)

    if hasattr(ls_result, 'jac') and ls_result.jac is not None:
        cov_result = compute_covariance_matrix(
            jacobian=ls_result.jac,
            residuals=final_residuals,
            n_params=n_params_full,
            fixed_params=fixed_params
        )
        params_stderr = cov_result.stderr
        cov = cov_result.cov
        condition_number = cov_result.condition_number
        is_well_conditioned = cov_result.is_well_conditioned
    else:
        params_stderr = np.full(n_params_full, np.inf)
        cov = None
        condition_number = np.inf
        is_well_conditioned = False

    # Step 4: Update circuit with fitted parameters
    if hasattr(circuit, 'update_params'):
        circuit.update_params(list(params_opt))

    # Create indexed param labels
    if param_labels is not None:
        label_counts = {}
        param_labels_indexed = []
        for label in param_labels:
            if label in label_counts:
                label_counts[label] += 1
            else:
                label_counts[label] = 0
            param_labels_indexed.append(f"{label}{label_counts[label]}")
    else:
        param_labels_indexed = None

    # Build FitDiagnostics
    fit_diagnostics = FitDiagnostics(
        optimizer_status=ls_result.status if hasattr(ls_result, 'status') else -1,
        optimizer_message=ls_result.message if hasattr(ls_result, 'message') else 'DE only',
        optimizer_success=ls_result.success if hasattr(ls_result, 'success') else de_result.success,
        n_function_evals=de_result.nfev + (ls_result.nfev if hasattr(ls_result, 'nfev') else 0),
        jacobian_type=jacobian_type,
        condition_number=condition_number,
        covariance_rank=cov_result.rank if cov_result else 0,
        covariance_warning=cov_result.warning_message if cov_result else None,
        warnings=diag_warnings
    )

    # Step 5: Create FitResult
    n_data = len(frequencies)
    fit_result = FitResult(
        circuit=circuit,
        params_opt=params_opt,
        params_stderr=params_stderr,
        fit_error_rel=fit_error_rel,
        fit_error_abs=fit_error_abs,
        quality=quality,
        condition_number=condition_number,
        is_well_conditioned=is_well_conditioned,
        cov=cov,
        diagnostics=fit_diagnostics,
        param_labels=param_labels_indexed,
        _n_data=n_data
    )

    # Build DiffEvoDiagnostics
    de_diagnostics = DiffEvoDiagnostics(
        strategy=strategy_name,
        popsize=popsize,
        maxiter=maxiter,
        tol=tol,
        workers=workers,
        weighting=weighting,
        jacobian_type=jacobian_type,
        de_converged=de_result.success,
        de_iterations=de_result.nit,
        de_evaluations=de_result.nfev,
        de_error=de_error_rel,
        refined_error=ls_error_rel,
        refinement_improved=used_refinement,
        total_evaluations=de_result.nfev + (ls_result.nfev if hasattr(ls_result, 'nfev') else 0),
        n_fixed_params=len(fixed_param_indices),
        fixed_param_indices=fixed_param_indices,
        warnings=diag_warnings
    )

    # Step 6: Create visualization
    from ..visualization.plots import plot_circuit_fit

    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = circuit.impedance(freq_plot, list(params_opt))
    fig = plot_circuit_fit(frequencies, Z, Z_fit_plot, circuit, Z_fit_at_data=Z_fit)

    # Create result object
    diffevo_result = DiffEvoResult(
        best_result=fit_result,
        de_result=de_result,
        de_error=de_error_rel,
        final_error=fit_error_rel,
        n_evaluations=de_diagnostics.total_evaluations,
        strategy=strategy_name,
        improvement=improvement,
        diagnostics=de_diagnostics
    )

    return diffevo_result, Z_fit, fig


__all__ = [
    'DiffEvoResult',
    'DiffEvoDiagnostics',
    'DE_STRATEGIES',
    'fit_circuit_diffevo',
]
