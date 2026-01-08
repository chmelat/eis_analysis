"""
Differential Evolution global optimization for circuit fitting.

Provides global optimization using scipy's differential_evolution,
followed by least_squares refinement for accurate parameter estimation
and covariance computation.

Strategy options:
    1 = 'randtobest1bin' (default) - balanced exploration/exploitation
    2 = 'best1bin' - faster convergence, may miss global minimum
    3 = 'rand1bin' - better exploration, slower convergence
"""

import numpy as np
import logging
import warnings
from typing import Tuple
from numpy.typing import NDArray
from dataclasses import dataclass
from scipy.optimize import differential_evolution, least_squares, OptimizeWarning

from .circuit import FitResult, Circuit
from .bounds import generate_simple_bounds
from .covariance import compute_covariance_matrix
from .diagnostics import compute_weights, compute_fit_metrics, log_fit_results
from .jacobian import make_jacobian_function

logger = logging.getLogger(__name__)

# Strategy mapping: number -> scipy strategy name
DE_STRATEGIES = {
    1: 'randtobest1bin',
    2: 'best1bin',
    3: 'rand1bin',
}


class _DECostFunction:
    """
    Picklable cost function for differential evolution with workers > 1.

    Must be defined at module level to be picklable by multiprocessing.
    Supports fixed parameters - only free parameters are optimized.
    """
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
        """Reconstruct full parameter list from free parameters."""
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


def _log_separator(length: int = 50, char: str = "=") -> None:
    """Log a separator line for visual clarity."""
    logger.info(char * length)


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
    """
    best_result: FitResult
    de_result: any  # scipy.optimize.OptimizeResult
    de_error: float
    final_error: float
    n_evaluations: int
    strategy: str
    improvement: float


def fit_circuit_diffevo(
    circuit: Circuit,
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    strategy: int = 1,
    popsize: int = 15,
    maxiter: int = 1000,
    tol: float = 0.01,
    workers: int = 1,
    weighting: str = 'proportional',
    verbose: bool = True,
    use_analytic_jacobian: bool = True
) -> Tuple[DiffEvoResult, NDArray[np.complex128], any]:
    """
    Fit circuit using Differential Evolution global optimization.

    Algorithm:
    1. Run differential_evolution to find global minimum region
    2. Refine with least_squares for accurate parameters and covariance
    3. Return best result with error estimates

    Parameters
    ----------
    circuit : Circuit
        Circuit object with initial parameter guesses (used for bounds)
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance data [Ohm]
    strategy : int, optional
        DE strategy: 1='randtobest1bin' (default), 2='best1bin', 3='rand1bin'
    popsize : int, optional
        Population size multiplier (default: 15, actual size = popsize * n_params)
    maxiter : int, optional
        Maximum number of generations (default: 1000)
    tol : float, optional
        Relative tolerance for convergence (default: 0.01)
    workers : int, optional
        Number of parallel workers (default: 1, use -1 for all CPUs)
    weighting : str, optional
        Weighting scheme: 'uniform', 'sqrt', 'proportional', 'modulus'
    verbose : bool, optional
        Log progress (default: True)
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian for least_squares refinement (default: True).
        More accurate and faster than numeric approximation. Set to False
        only if using unsupported custom elements.

    Returns
    -------
    diffevo_result : DiffEvoResult
        Differential evolution result with best fit and statistics
    Z_fit : ndarray
        Best fit impedance
    fig : matplotlib.figure.Figure
        Nyquist plot with best fit

    Notes
    -----
    DE is a stochastic global optimizer that:
    - Does not require gradients
    - Handles bounds naturally
    - Can escape local minima

    The least_squares refinement after DE provides:
    - More accurate parameters (quadratic convergence)
    - Covariance matrix for uncertainty estimation
    - Better numerical stability

    Example
    -------
    >>> circuit = R(100) - (R(5000) | C(1e-6))
    >>> result, Z_fit, fig = fit_circuit_diffevo(circuit, freq, Z, strategy=1)
    >>> print(f"Final error: {result.final_error:.2f}%")
    """
    strategy_name = DE_STRATEGIES.get(strategy, 'randtobest1bin')

    if verbose:
        weighting_labels = {
            'uniform': 'uniform (w=1)',
            'sqrt': 'sqrt (w=1/sqrt|Z|)',
            'proportional': 'proportional (w=1/|Z|)',
            'modulus': 'modulus (w=1/|Z|^2)'
        }
        _log_separator()
        logger.info("Differential Evolution optimization")
        _log_separator()
        logger.info(f"  Strategy: {strategy_name} (option {strategy})")
        logger.info(f"  Population: {popsize} * n_params")
        logger.info(f"  Max iterations: {maxiter}")
        logger.info(f"  Tolerance: {tol}")
        logger.info(f"  Workers: {workers}")
        logger.info(f"  Weighting: {weighting_labels.get(weighting, weighting)}")

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
    if hasattr(circuit, 'get_all_fixed_params'):
        fixed_params = circuit.get_all_fixed_params()
        n_fixed = sum(fixed_params)
        if n_fixed > 0 and verbose:
            logger.info(f"  Fixed parameters: {n_fixed} of {len(fixed_params)}")
            for i, (is_fixed, val) in enumerate(zip(fixed_params, initial_guess_full)):
                if is_fixed:
                    label = param_labels[i] if param_labels else f"p{i}"
                    logger.info(f"    {label} = {val:.4g} (fixed)")

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

    # Clip initial guess to bounds (in case user provided out-of-bounds values)
    initial_guess = np.clip(initial_guess, lower_bounds, upper_bounds)

    if verbose:
        logger.info(f"  Parameters: {n_params} (free)")
        logger.info(f"  Initial guess: {initial_guess}")
        logger.info("")

    # Precompute weights
    weights = compute_weights(Z, weighting)

    # Cost function for DE - use picklable class for multiprocessing support
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

    # Residual function for least_squares (returns vector)
    def residual_function(params):
        full_params = reconstruct_params(params)
        Z_pred = circuit.impedance(frequencies, full_params)
        return np.concatenate([
            (Z.real - Z_pred.real) * weights,
            (Z.imag - Z_pred.imag) * weights
        ])

    # Step 1: Run Differential Evolution
    if verbose:
        logger.info("Running differential evolution...")

    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            de_result = differential_evolution(
                cost_function,
                bounds,
                x0=initial_guess,  # Seed one population member with initial guess
                strategy=strategy_name,
                popsize=popsize,
                maxiter=maxiter,
                tol=tol,
                workers=workers,
                polish=False,  # We do our own refinement
                seed=None,  # Random seed for reproducibility if needed
                disp=False,
                updating='deferred' if workers != 1 else 'immediate',
            )
    except Exception as e:
        logger.error(f"Differential evolution failed: {e}")
        raise RuntimeError(f"DE optimization failed: {e}") from e

    # Compute DE error (reconstruct full params)
    de_params_full = reconstruct_params(de_result.x)
    Z_fit_de = circuit.impedance(frequencies, de_params_full)
    de_error_rel, _, _ = compute_fit_metrics(Z, Z_fit_de, weighting)

    if verbose:
        logger.info(f"  DE converged: {de_result.success}")
        logger.info(f"  DE iterations: {de_result.nit}")
        logger.info(f"  DE evaluations: {de_result.nfev}")
        logger.info(f"  DE error: {de_error_rel:.3f}%")
        logger.info("")

    # Step 2: Refine with least_squares
    if verbose:
        jac_type = "analytic" if use_analytic_jacobian else "numeric"
        logger.info(f"Refining with least_squares ({jac_type} Jacobian)...")

    # Prepare Jacobian function if using analytic
    if use_analytic_jacobian:
        try:
            jac_func = make_jacobian_function(
                circuit, frequencies, weights,
                fixed_params=fixed_params,
                full_initial_guess=initial_guess_full
            )
        except NotImplementedError as e:
            logger.warning(f"Analytic Jacobian failed: {e}, falling back to numeric")
            jac_func = '2-point'
            use_analytic_jacobian = False
    else:
        jac_func = '2-point'

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
        logger.warning(f"least_squares refinement failed: {e}, using DE result")
        ls_result = de_result
        ls_result.x = de_result.x

    # Reconstruct full params (free params from optimizer + fixed params)
    ls_params_free = np.array(ls_result.x)
    ls_params_full = np.array(reconstruct_params(ls_params_free))

    # Compute refined fit error
    Z_fit_ls = circuit.impedance(frequencies, list(ls_params_full))
    ls_error_rel, _, _ = compute_fit_metrics(Z, Z_fit_ls, weighting)

    improvement = (de_error_rel - ls_error_rel) / de_error_rel * 100 if de_error_rel > 0 else 0

    # Choose better result: DE or refined
    if ls_error_rel <= de_error_rel:
        # Refinement improved or matched - use refined result
        params_opt_free = ls_params_free
        params_opt = ls_params_full
        Z_fit = Z_fit_ls
        fit_error_rel = ls_error_rel
        used_refinement = True
    else:
        # Refinement worsened - keep DE result
        params_opt_free = np.array(de_result.x)
        params_opt = np.array(de_params_full)
        Z_fit = Z_fit_de
        fit_error_rel = de_error_rel
        used_refinement = False

    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

    if verbose:
        logger.info(f"  Refined error: {ls_error_rel:.3f}%")
        logger.info(f"  Improvement: {improvement:+.1f}%")
        if not used_refinement:
            logger.info("  Refinement worsened fit, using DE result")
        logger.info("")

    # Step 3: Compute covariance
    # Compute residuals at final parameters (use free params for residual_function)
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
        # Fallback: compute Jacobian numerically
        try:
            from scipy.optimize import approx_fprime
            eps = 1e-8
            jac = np.zeros((2 * len(Z), n_params))
            for i in range(n_params):
                def f_i(x, idx=i):
                    p = params_opt_free.copy()
                    p[idx] = x
                    return residual_function(p)
                jac[:, i] = approx_fprime([params_opt_free[i]], lambda x: residual_function(
                    np.array([params_opt_free[j] if j != i else x[0] for j in range(n_params)])
                ), eps).flatten()[:2*len(Z)]
            cov_result = compute_covariance_matrix(
                jacobian=jac,
                residuals=final_residuals,
                n_params=n_params_full,
                fixed_params=fixed_params
            )
            params_stderr = cov_result.stderr
            cov = cov_result.cov
            condition_number = cov_result.condition_number
            is_well_conditioned = cov_result.is_well_conditioned
        except Exception:
            params_stderr = np.full(n_params_full, np.inf)
            cov = None
            condition_number = np.inf
            is_well_conditioned = False

    # Step 4: Update circuit with fitted parameters
    if hasattr(circuit, 'update_params'):
        circuit.update_params(list(params_opt))

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
        _n_data=n_data
    )

    # Create indexed param labels for logging
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

    # Step 5: Log results
    if verbose:
        _log_separator()
        logger.info("Differential Evolution results")
        _log_separator()
        logger.info(f"  Strategy: {strategy_name}")
        logger.info(f"  Total evaluations: {de_result.nfev + (ls_result.nfev if hasattr(ls_result, 'nfev') else 0)}")
        logger.info(f"  DE error: {de_error_rel:.3f}% -> Refined: {ls_error_rel:.3f}%")
        if not used_refinement:
            logger.info("  Using DE result (refinement worsened fit)")
        logger.info("")

        ci_low, ci_high = fit_result.params_ci_95
        log_fit_results(
            params_opt, params_stderr, ci_low, ci_high,
            fit_error_rel, fit_error_abs, quality, param_labels_indexed
        )
        _log_separator()

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
        n_evaluations=de_result.nfev + (ls_result.nfev if hasattr(ls_result, 'nfev') else 0),
        strategy=strategy_name,
        improvement=improvement
    )

    return diffevo_result, Z_fit, fig


__all__ = [
    'DiffEvoResult',
    'DE_STRATEGIES',
    'fit_circuit_diffevo',
]
