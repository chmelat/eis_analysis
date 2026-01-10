"""
Multi-start optimization for circuit fitting.

Clean design: No logging in core functions, all diagnostics returned as data.

Provides adaptive multi-start optimization that uses covariance information
from initial fits to intelligently generate perturbations.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List
from numpy.typing import NDArray
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed

from .circuit import fit_equivalent_circuit, FitResult, Circuit
from .bounds import generate_simple_bounds

logger = logging.getLogger(__name__)


@dataclass
class MultistartDiagnostics:
    """Diagnostics from multi-start optimization."""
    n_restarts: int
    n_successful: int
    scale: float
    weighting: str
    jacobian_type: str
    parallel: bool

    # Error progression
    initial_error: float
    best_error: float
    best_start_index: int
    all_errors: List[Optional[float]]

    # Perturbation method used
    perturbation_method: str  # 'covariance', 'stderr', 'log_uniform'

    warnings: List[str] = field(default_factory=list)
    failed_errors: List[str] = field(default_factory=list)  # Exception messages from failed fits


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
    diagnostics : MultistartDiagnostics
        Detailed diagnostics
    """
    best_result: FitResult
    all_results: List[FitResult]
    n_starts: int
    n_successful: int
    improvement: float
    diagnostics: Optional[MultistartDiagnostics] = None


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
    """
    n_params = len(params)

    try:
        cov_reg = cov + 1e-10 * np.eye(n_params)
        L = np.linalg.cholesky(cov_reg)
        z = np.random.randn(n_params)
        perturbation = scale * L @ z
        perturbed = params + perturbation

    except np.linalg.LinAlgError:
        stderr = np.sqrt(np.abs(np.diag(cov)))
        perturbation = scale * stderr * np.random.randn(n_params)
        perturbed = params + perturbation

    if bounds is not None:
        lower, upper = bounds
        perturbed = np.clip(perturbed, lower, upper)

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
    """
    stderr_safe = np.where(
        np.isfinite(stderr) & (stderr > 0),
        stderr,
        np.abs(params) * 0.1
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
    verbose: bool = True,  # Kept for backward compatibility, ignored
    use_analytic_jacobian: bool = True
) -> Tuple[MultistartResult, NDArray[np.complex128], any]:
    """
    Fit circuit using adaptive multi-start optimization.

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
        Weighting scheme
    parallel : bool, optional
        Use parallel execution (default: False)
    max_workers : int, optional
        Maximum parallel workers (default: 4)
    verbose : bool, optional
        Ignored (kept for backward compatibility)
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian (default: True)

    Returns
    -------
    multistart_result : MultistartResult
        Multi-start optimization result with all diagnostics
    Z_fit : ndarray
        Best fit impedance
    fig : matplotlib.figure.Figure
        Nyquist plot with best fit
    """
    all_results = []
    all_errors = []
    n_successful = 0
    diag_warnings = []
    failed_errors = []  # Track exceptions from failed fits
    perturbation_method = 'covariance'

    jacobian_type = 'analytic' if use_analytic_jacobian else 'numeric'

    # Step 1: Initial fit
    try:
        result0, Z_fit0, _ = fit_equivalent_circuit(
            frequencies, Z, circuit, weighting=weighting, plot=False,
            use_analytic_jacobian=use_analytic_jacobian
        )
        all_results.append(result0)
        all_errors.append(result0.fit_error_rel)
        n_successful += 1
        initial_error = result0.fit_error_rel

    except Exception as e:
        raise RuntimeError(f"Multi-start failed: initial fit unsuccessful: {e}") from e

    # Extract bounds from circuit
    param_labels = circuit.get_param_labels() if hasattr(circuit, 'get_param_labels') else None
    if param_labels is not None:
        lower_bounds, upper_bounds = generate_simple_bounds(param_labels)
        lower_bounds = np.array(lower_bounds)
        upper_bounds = np.array(upper_bounds)
    else:
        n_params = len(result0.params_opt)
        lower_bounds = np.full(n_params, 1e-15)
        upper_bounds = np.full(n_params, 1e15)
    bounds = (lower_bounds, upper_bounds)

    # Step 2: Generate perturbations and run additional fits
    def run_single_fit(start_idx: int, initial_params: NDArray) -> Optional[FitResult]:
        try:
            result, _, _ = fit_equivalent_circuit(
                frequencies, Z, circuit,
                weighting=weighting,
                initial_guess=list(initial_params),
                plot=False,
                use_analytic_jacobian=use_analytic_jacobian
            )
            return result
        except Exception as e:
            logger.debug(f"Multistart fit #{start_idx} failed: {e}")
            failed_errors.append(f"Start #{start_idx}: {e}")
            return None

    # Generate all perturbations
    perturbations = []
    for i in range(1, n_restarts):
        if result0.cov is not None and result0.is_well_conditioned:
            perturbed = perturb_from_covariance(
                result0.params_opt, result0.cov, scale=scale, bounds=bounds
            )
            perturbation_method = 'covariance'
        elif not np.any(np.isinf(result0.params_stderr)):
            perturbed = perturb_from_stderr(
                result0.params_opt, result0.params_stderr, scale=scale, bounds=bounds
            )
            perturbation_method = 'stderr'
        else:
            perturbed = perturb_log_uniform(
                result0.params_opt, factor=3.0, bounds=bounds
            )
            perturbation_method = 'log_uniform'
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

    # Build diagnostics
    diagnostics = MultistartDiagnostics(
        n_restarts=n_restarts,
        n_successful=n_successful,
        scale=scale,
        weighting=weighting,
        jacobian_type=jacobian_type,
        parallel=parallel,
        initial_error=initial_error,
        best_error=best_error,
        best_start_index=best_idx,
        all_errors=all_errors,
        perturbation_method=perturbation_method,
        warnings=diag_warnings,
        failed_errors=failed_errors
    )

    # Create visualization for best result
    from ..visualization.plots import plot_circuit_fit

    Z_fit_best = best_result.circuit.impedance(frequencies, list(best_result.params_opt))

    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = best_result.circuit.impedance(freq_plot, list(best_result.params_opt))
    fig = plot_circuit_fit(frequencies, Z, Z_fit_plot, best_result.circuit, Z_fit_at_data=Z_fit_best)

    multistart_result = MultistartResult(
        best_result=best_result,
        all_results=all_results,
        n_starts=n_restarts,
        n_successful=n_successful,
        improvement=improvement,
        diagnostics=diagnostics
    )

    return multistart_result, Z_fit_best, fig


__all__ = [
    'MultistartResult',
    'MultistartDiagnostics',
    'perturb_from_covariance',
    'perturb_from_stderr',
    'perturb_log_uniform',
    'fit_circuit_multistart',
]
