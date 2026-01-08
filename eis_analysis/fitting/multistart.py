"""
Multi-start optimization for circuit fitting.

Provides adaptive multi-start optimization that uses covariance information
from initial fits to intelligently generate perturbations for subsequent
optimization attempts.

The key insight is that parameters with large uncertainty (high stderr)
have a "flat" cost function landscape, so exploring different values
is more likely to find better minima. Parameters with low uncertainty
are well-determined by the data and should be perturbed less.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List
from numpy.typing import NDArray
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed

from .circuit import fit_equivalent_circuit, FitResult, Circuit
from .bounds import generate_simple_bounds
from .diagnostics import log_fit_results

logger = logging.getLogger(__name__)


def _log_separator(length: int = 50, char: str = "=") -> None:
    """Log a separator line for visual clarity."""
    logger.info(char * length)


@dataclass
class MultistartResult:
    """
    Result from multi-start optimization.

    Attributes
    ----------
    best_result : FitResult
        Best fitting result (lowest error)
    all_results : list of FitResult
        All fitting results from all starts
    n_starts : int
        Number of optimization starts performed
    n_successful : int
        Number of successful optimizations
    improvement : float
        Relative improvement over initial fit [%]
    """
    best_result: FitResult
    all_results: List[FitResult]
    n_starts: int
    n_successful: int
    improvement: float


def perturb_from_covariance(
    params: NDArray[np.float64],
    cov: NDArray[np.float64],
    scale: float = 2.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None
) -> NDArray[np.float64]:
    """
    Generate correlated perturbation using Cholesky decomposition of covariance.

    Parameters with high covariance (uncertainty) are perturbed more,
    and correlations between parameters are preserved.

    Parameters
    ----------
    params : ndarray
        Current parameter values
    cov : ndarray
        Covariance matrix from previous fit
    scale : float, optional
        Scale factor for perturbation magnitude (default: 2.0 = 2 sigma)
    bounds : tuple of (lower, upper), optional
        Parameter bounds for clipping

    Returns
    -------
    perturbed : ndarray
        Perturbed parameters

    Notes
    -----
    Uses Cholesky decomposition: cov = L @ L.T
    Then: perturbation = scale * L @ z, where z ~ N(0, I)

    This preserves parameter correlations. For example, if R and C in a
    Voigt element are anticorrelated (tau = RC is well-determined but
    individual values are not), the perturbation will move them together
    in a way that approximately preserves tau.
    """
    n_params = len(params)

    try:
        # Add small regularization for numerical stability
        cov_reg = cov + 1e-10 * np.eye(n_params)

        # Cholesky decomposition
        L = np.linalg.cholesky(cov_reg)

        # Generate correlated random vector
        z = np.random.randn(n_params)
        perturbation = scale * L @ z

        perturbed = params + perturbation

    except np.linalg.LinAlgError:
        # Fallback to diagonal (uncorrelated) perturbation
        logger.debug("Cholesky failed, using diagonal perturbation")
        stderr = np.sqrt(np.abs(np.diag(cov)))
        perturbation = scale * stderr * np.random.randn(n_params)
        perturbed = params + perturbation

    # Apply bounds if provided
    if bounds is not None:
        lower, upper = bounds
        perturbed = np.clip(perturbed, lower, upper)

    # Ensure positive values for typical EIS parameters
    perturbed = np.maximum(perturbed, 1e-15)

    return perturbed


def perturb_from_stderr(
    params: NDArray[np.float64],
    stderr: NDArray[np.float64],
    scale: float = 2.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None
) -> NDArray[np.float64]:
    """
    Generate uncorrelated perturbation scaled by standard errors.

    Simpler alternative to covariance-based perturbation when full
    covariance matrix is not available or is ill-conditioned.

    Parameters
    ----------
    params : ndarray
        Current parameter values
    stderr : ndarray
        Standard errors of parameters
    scale : float, optional
        Scale factor (default: 2.0 = 2 sigma)
    bounds : tuple of (lower, upper), optional
        Parameter bounds for clipping

    Returns
    -------
    perturbed : ndarray
        Perturbed parameters
    """
    # Handle infinite or zero stderr
    stderr_safe = np.where(
        np.isfinite(stderr) & (stderr > 0),
        stderr,
        np.abs(params) * 0.1  # Fallback: 10% of parameter value
    )

    perturbation = scale * stderr_safe * np.random.randn(len(params))
    perturbed = params + perturbation

    if bounds is not None:
        lower, upper = bounds
        perturbed = np.clip(perturbed, lower, upper)

    perturbed = np.maximum(perturbed, 1e-15)

    return perturbed


def perturb_log_uniform(
    params: NDArray[np.float64],
    factor: float = 3.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None
) -> NDArray[np.float64]:
    """
    Generate log-uniform perturbation (multiplicative).

    Each parameter is multiplied by a random factor in [1/factor, factor].
    Useful for parameters spanning many orders of magnitude.

    Parameters
    ----------
    params : ndarray
        Current parameter values
    factor : float, optional
        Maximum multiplicative factor (default: 3.0)
    bounds : tuple of (lower, upper), optional
        Parameter bounds for clipping

    Returns
    -------
    perturbed : ndarray
        Perturbed parameters
    """
    log_factor = np.log(factor)
    multipliers = np.exp(np.random.uniform(-log_factor, log_factor, len(params)))
    perturbed = params * multipliers

    if bounds is not None:
        lower, upper = bounds
        perturbed = np.clip(perturbed, lower, upper)

    perturbed = np.maximum(perturbed, 1e-15)

    return perturbed


def fit_circuit_multistart(
    circuit: Circuit,
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_restarts: int = 10,
    scale: float = 2.0,
    weighting: str = 'modulus',
    parallel: bool = False,
    max_workers: int = 4,
    verbose: bool = True,
    use_analytic_jacobian: bool = True
) -> Tuple[MultistartResult, NDArray[np.complex128], any]:
    """
    Fit circuit using adaptive multi-start optimization.

    Algorithm:
    1. Perform initial fit from circuit's initial guess
    2. Use covariance matrix to generate intelligent perturbations
    3. Run additional optimizations from perturbed starting points
    4. Return best result

    Parameters
    ----------
    circuit : Circuit
        Circuit object with initial parameter guesses
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance data [Ohm]
    n_restarts : int, optional
        Total number of optimization starts (default: 10)
    scale : float, optional
        Perturbation scale in units of sigma (default: 2.0)
    weighting : str, optional
        Weighting scheme: 'uniform', 'sqrt', 'proportional', 'modulus'
    parallel : bool, optional
        Use parallel execution (default: False)
    max_workers : int, optional
        Maximum parallel workers (default: 4)
    verbose : bool, optional
        Log progress (default: True)
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian for least_squares optimization (default: True).
        More accurate and faster than numeric approximation.

    Returns
    -------
    multistart_result : MultistartResult
        Multi-start optimization result with best fit and statistics
    Z_fit : ndarray
        Best fit impedance
    fig : matplotlib.figure.Figure
        Nyquist plot with best fit

    Notes
    -----
    The adaptive strategy uses information from the covariance matrix:
    - Parameters with high uncertainty are perturbed more
    - Correlated parameters are perturbed together
    - This is more efficient than random sampling

    Example
    -------
    >>> circuit = R(100) - (R(5000) | C(1e-6))
    >>> result, Z_fit, fig = fit_circuit_multistart(circuit, freq, Z, n_restarts=10)
    >>> print(f"Best error: {result.best_result.fit_error_rel:.2f}%")
    >>> print(f"Improvement: {result.improvement:.1f}%")
    """
    if verbose:
        weighting_labels = {
            'uniform': 'uniform (w=1)',
            'sqrt': 'sqrt (w=1/sqrt|Z|)',
            'proportional': 'proportional (w=1/|Z|)',
            'modulus': 'modulus (w=1/|Z|^2)'
        }
        _log_separator()
        logger.info("Multi-start optimization")
        _log_separator()
        jac_type = "analytic" if use_analytic_jacobian else "numeric"
        logger.info(f"  Restarts: {n_restarts}, scale: {scale} sigma")
        logger.info(f"  Jacobian: {jac_type}")
        logger.info(f"  Weighting: {weighting_labels.get(weighting, weighting)}")

    all_results = []
    all_errors = []  # For condensed progress display
    n_successful = 0

    # Step 1: Initial fit (no plot, no verbose - we'll show results at the end)
    try:
        result0, Z_fit0, _ = fit_equivalent_circuit(
            frequencies, Z, circuit, weighting=weighting, plot=False, verbose=False,
            use_analytic_jacobian=use_analytic_jacobian
        )
        all_results.append(result0)
        all_errors.append(result0.fit_error_rel)
        n_successful += 1
        initial_error = result0.fit_error_rel

    except Exception as e:
        logger.error(f"Initial fit failed: {e}")
        raise RuntimeError("Multi-start failed: initial fit unsuccessful") from e

    # Extract bounds from circuit
    param_labels = circuit.get_param_labels() if hasattr(circuit, 'get_param_labels') else None
    if param_labels is not None:
        lower_bounds, upper_bounds = generate_simple_bounds(param_labels)
        lower_bounds = np.array(lower_bounds)
        upper_bounds = np.array(upper_bounds)
    else:
        # Fallback: use wide bounds
        n_params = len(result0.params_opt)
        lower_bounds = np.full(n_params, 1e-15)
        upper_bounds = np.full(n_params, 1e15)
    bounds = (lower_bounds, upper_bounds)

    # Step 2: Generate perturbations and run additional fits
    def run_single_fit(start_idx: int, initial_params: NDArray) -> Optional[FitResult]:
        """Run single optimization from perturbed start (no plot, no verbose)."""
        try:
            # Use initial_guess parameter to override starting point
            # plot=False, verbose=False to reduce output during multi-start
            result, _, _ = fit_equivalent_circuit(
                frequencies, Z, circuit,
                weighting=weighting,
                initial_guess=list(initial_params),
                plot=False,
                verbose=False,
                use_analytic_jacobian=use_analytic_jacobian
            )
            return result
        except Exception as e:
            logger.debug(f"Start {start_idx} failed: {e}")
            return None

    # Generate all perturbations
    perturbations = []
    for i in range(1, n_restarts):
        if result0.cov is not None and result0.is_well_conditioned:
            # Use covariance-based perturbation
            perturbed = perturb_from_covariance(
                result0.params_opt, result0.cov, scale=scale, bounds=bounds
            )
        elif not np.any(np.isinf(result0.params_stderr)):
            # Fallback to stderr-based perturbation
            perturbed = perturb_from_stderr(
                result0.params_opt, result0.params_stderr, scale=scale, bounds=bounds
            )
        else:
            # Last resort: log-uniform perturbation
            perturbed = perturb_log_uniform(
                result0.params_opt, factor=3.0, bounds=bounds
            )
        perturbations.append((i + 1, perturbed))

    # Run fits (parallel or sequential)
    if parallel and n_restarts > 2:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(run_single_fit, idx, params): idx
                for idx, params in perturbations
            }

            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    all_results.append(result)
                    all_errors.append(result.fit_error_rel)
                    n_successful += 1
                else:
                    all_errors.append(None)
    else:
        for idx, params in perturbations:
            result = run_single_fit(idx, params)
            if result is not None:
                all_results.append(result)
                all_errors.append(result.fit_error_rel)
                n_successful += 1
            else:
                all_errors.append(None)

    # Step 3: Find best result
    if not all_results:
        raise RuntimeError("Multi-start failed: no successful fits")

    best_result = min(all_results, key=lambda r: r.fit_error_rel)
    best_error = best_result.fit_error_rel
    improvement = (initial_error - best_error) / initial_error * 100 if initial_error > 0 else 0

    # Find which start gave the best result
    best_idx = all_errors.index(best_error) + 1 if best_error in all_errors else 1

    if verbose:
        # Condensed progress display
        logger.info("")
        logger.info("Progress (error %):")
        error_strs = []
        for i, err in enumerate(all_errors):
            if err is None:
                error_strs.append("fail")
            elif i + 1 == best_idx:
                error_strs.append(f"[{err:.2f}]")  # Highlight best
            else:
                error_strs.append(f"{err:.2f}")
        # Display in rows of 10
        for row_start in range(0, len(error_strs), 10):
            row = error_strs[row_start:row_start + 10]
            logger.info(f"  {row_start + 1:2d}-{row_start + len(row):2d}: {', '.join(row)}")

        # Summary
        _log_separator()
        logger.info("Multi-start summary")
        _log_separator()
        logger.info(f"  Successful fits: {n_successful}/{n_restarts}")
        logger.info(f"  Initial error:   {initial_error:.3f}%")
        logger.info(f"  Best error:      {best_error:.3f}% (start #{best_idx})")
        logger.info(f"  Improvement:     {improvement:+.1f}%")

        # Detailed results of best fit
        _log_separator()
        logger.info("Best fit results")
        _log_separator()

        # Get parameter labels
        param_labels = None
        if hasattr(circuit, 'get_param_labels'):
            param_labels_raw = circuit.get_param_labels()
            label_counts = {}
            param_labels = []
            for label in param_labels_raw:
                if label in label_counts:
                    label_counts[label] += 1
                else:
                    label_counts[label] = 0
                param_labels.append(f"{label}{label_counts[label]}")

        # Use shared log function for consistent output
        ci_low, ci_high = best_result.params_ci_95
        log_fit_results(
            best_result.params_opt,
            best_result.params_stderr,
            ci_low, ci_high,
            best_result.fit_error_rel,
            best_result.fit_error_abs,
            best_result.quality,
            param_labels
        )
        _log_separator()

    # Create visualization for best result
    from ..visualization.plots import plot_circuit_fit

    Z_fit_best = best_result.circuit.impedance(frequencies, list(best_result.params_opt))

    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = best_result.circuit.impedance(freq_plot, list(best_result.params_opt))
    # Pass Z_fit_best (at original frequencies) for residuals calculation
    fig = plot_circuit_fit(frequencies, Z, Z_fit_plot, best_result.circuit, Z_fit_at_data=Z_fit_best)

    multistart_result = MultistartResult(
        best_result=best_result,
        all_results=all_results,
        n_starts=n_restarts,
        n_successful=n_successful,
        improvement=improvement
    )

    return multistart_result, Z_fit_best, fig


__all__ = [
    'MultistartResult',
    'perturb_from_covariance',
    'perturb_from_stderr',
    'perturb_log_uniform',
    'fit_circuit_multistart',
]
