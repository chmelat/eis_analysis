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
from typing import Tuple, List, Optional, Any
from numpy.typing import NDArray
from dataclasses import dataclass, field
from scipy.optimize import differential_evolution, least_squares, OptimizeWarning

from .circuit import FitResult, FitDiagnostics, Circuit
from .bounds import generate_simple_bounds, build_bound_status, log_scale_ci_mask
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
                 fixed_params: Optional[List[bool]] = None,
                 full_initial_guess: Optional[List[float]] = None):
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

    # Optimized objective values (weighted SSR, S = sum w^2 |dZ|^2).
    # These drive the DE-vs-refinement selection; the *_error fields above
    # are the human-readable weighted mean relative error (%) used for display.
    de_cost: float = 0.0
    refined_cost: float = 0.0

    # Fixed params info
    n_fixed_params: int = 0
    fixed_param_indices: List[int] = field(default_factory=list)

    # Initial guess passed to DE (full parameter vector, including fixed).
    # Captured before the optimizer runs so it survives circuit.update_params()
    # at the end of fit_circuit_diffevo.
    initial_guess: List[float] = field(default_factory=list)

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
    de_result: Any
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
    use_analytic_jacobian: bool = True,
    seed: Optional[int] = None
) -> Tuple[DiffEvoResult, NDArray[np.complex128], Any]:
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
    seed : int, optional
        Seed for differential_evolution's random generator. Default None
        (non-deterministic). Set an int for reproducible runs (e.g. tests).

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
                seed=seed,
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

    refinement_ran = True
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
        refinement_ran = False
        diag_warnings.append(f"Refinement failed: {e}, using DE result")
        ls_result = de_result

    # Reconstruct full params
    ls_params_free = np.array(ls_result.x)
    ls_params_full = np.array(reconstruct_params(ls_params_free))

    # Compute refined fit error
    Z_fit_ls = circuit.impedance(frequencies, list(ls_params_full))
    ls_error_rel, _, _ = compute_fit_metrics(Z, Z_fit_ls, weighting)

    # Selection and improvement use the *optimized* objective (weighted SSR,
    # S = sum w^2 |dZ|^2), not the weighted mean relative error: DE and
    # least_squares both minimize S, so choosing on a different metric could
    # discard a genuinely better refined fit. cost_function returns S directly.
    de_cost = float(cost_function(de_result.x))
    ls_cost = float(cost_function(ls_result.x))

    improvement = (de_cost - ls_cost) / de_cost * 100 if de_cost > 0 else 0

    # Choose better result. A failed refinement must not masquerade as a
    # successful one: ls_result then aliases de_result (equal costs), so gate
    # on refinement_ran as well (audit 2026-07-02 finding 2.3).
    if refinement_ran and ls_cost <= de_cost:
        params_opt_free = ls_params_free
        params_opt = ls_params_full
        Z_fit = Z_fit_ls
        used_refinement = True
    else:
        params_opt_free = np.array(de_result.x)
        params_opt = np.array(de_params_full)
        Z_fit = Z_fit_de
        used_refinement = False
        if refinement_ran:
            diag_warnings.append("Refinement worsened fit, using DE result")

    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

    # Step 3: Compute covariance
    # Both the residuals and the Jacobian must be evaluated at the *chosen*
    # point (params_opt_free). When the DE result is kept, ls_result.jac is the
    # Jacobian at the LS point, not the chosen one, so it must not be reused.
    final_residuals = residual_function(params_opt_free)
    n_params_full = len(params_opt)

    if jacobian_type == 'analytic':
        jac_at_opt = jac_func(params_opt_free)
    elif used_refinement and getattr(ls_result, 'jac', None) is not None:
        # Numeric Jacobian; LS point coincides with the chosen point.
        jac_at_opt = ls_result.jac
    else:
        jac_at_opt = None

    cov_result = None
    if jac_at_opt is not None:
        cov_result = compute_covariance_matrix(
            jacobian=jac_at_opt,
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

    # Per-parameter bound status and derived warnings — same contract as
    # fit_equivalent_circuit (full-space indices, classify_bound_status
    # criterion, labels when available).
    bound_status = build_bound_status(
        params_opt, lower_bounds_full, upper_bounds_full, fixed_params
    )
    bounds_warnings = []
    params_at_bounds = []
    for i, status in enumerate(bound_status):
        if status not in ('lower', 'upper'):
            continue
        params_at_bounds.append(i)
        name = param_labels_indexed[i] if param_labels_indexed is not None else str(i)
        bound_val = lower_bounds_full[i] if status == 'lower' else upper_bounds_full[i]
        bounds_warnings.append(
            f"Parameter {name} = {params_opt[i]:.3e} near {status} "
            f"bound {bound_val:.1e}"
        )

    # Build FitDiagnostics
    # Optimizer metadata belongs to least_squares only when it actually ran;
    # on failure ls_result aliases de_result, whose message/nfev would be
    # misattributed (and DE evaluations double-counted).
    fit_diagnostics = FitDiagnostics(
        optimizer_status=ls_result.status if refinement_ran else -1,
        optimizer_message=ls_result.message if refinement_ran else 'DE only (refinement failed)',
        optimizer_success=ls_result.success if refinement_ran else de_result.success,
        n_function_evals=de_result.nfev + (ls_result.nfev if refinement_ran else 0),
        jacobian_type=jacobian_type,
        condition_number=condition_number,
        covariance_rank=cov_result.rank if cov_result else 0,
        covariance_warning=cov_result.warning_message if cov_result else None,
        params_at_bounds=params_at_bounds,
        bounds_warnings=bounds_warnings,
        warnings=diag_warnings
    )

    # Step 5: Create FitResult
    # When cov_result is None, stderr is inf so the CI is +/-inf regardless of dof.
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
        bound_status=bound_status,
        _dof=cov_result.dof if cov_result is not None else 0,
        _ci_log_scale=log_scale_ci_mask(lower_bounds_full, upper_bounds_full)
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
        de_cost=de_cost,
        refined_cost=ls_cost,
        refinement_improved=used_refinement,
        total_evaluations=de_result.nfev + (ls_result.nfev if refinement_ran else 0),
        n_fixed_params=len(fixed_param_indices),
        fixed_param_indices=fixed_param_indices,
        initial_guess=list(initial_guess_full),
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
