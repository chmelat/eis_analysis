"""
Multi-start optimization for circuit fitting.

Clean design: No logging in core functions, all diagnostics returned as data.

Provides adaptive multi-start optimization that uses covariance information
from initial fits to intelligently generate perturbations.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List, Union
from numpy.typing import NDArray
from dataclasses import dataclass, field
from concurrent.futures import ThreadPoolExecutor, as_completed
from copy import deepcopy

from .circuit import fit_equivalent_circuit, FitResult, Circuit
from .bounds import generate_simple_bounds
from .diagnostics import compute_information_criteria

logger = logging.getLogger(__name__)

# Anything np.random.default_rng() accepts: None (fresh entropy), an int seed,
# or a Generator (used as is, so one stream can feed several calls)
Seed = Union[None, int, np.random.Generator]


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
        Relative improvement of the weighted SSR over the initial fit [%] -
        the quantity the best start is selected on, so it is never negative
    diagnostics : MultistartDiagnostics
        Detailed diagnostics
    """
    best_result: FitResult
    all_results: List[FitResult]
    n_starts: int
    n_successful: int
    improvement: float
    diagnostics: Optional[MultistartDiagnostics] = None


def _clip_to_bounds(perturbed, bounds):
    """Clip a perturbed start into `bounds`, or to >= 0 without them.

    Every parameter of the library is non-negative, so 0 is the physical
    floor. The former absolute floor of 1e-15 lifted G (bounds 0..1e4 S)
    and alpha_CC (0..0.9) off a valid 0 and pushed a G below 1e-15 S - an
    oxide's conductance - up to it. A zero start is fine: the optimizer
    gives it a step scale only (optimizer.ZERO_START_SCALE).
    """
    lower, upper = bounds if bounds is not None else (0.0, np.inf)
    return np.clip(perturbed, lower, upper)


def perturb_from_covariance(
    params: NDArray[np.float64],
    cov: NDArray[np.float64],
    scale: float = 2.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None,
    rng: Seed = None
) -> NDArray[np.float64]:
    """
    Generate correlated perturbation using Cholesky decomposition of covariance.

    Parameters with high covariance (uncertainty) are perturbed more,
    and correlations between parameters are preserved.
    """
    rng = np.random.default_rng(rng)
    n_params = len(params)

    try:
        # Sample in the correlation matrix, then scale by the standard errors.
        # The regularization that makes Cholesky work (fixed parameters leave
        # zero rows) must be relative: an absolute 1e-10 added to cov swamped
        # the variance of a capacitor (C ~ 1e-7 F -> var ~ 1e-20) and spread
        # its perturbation ~2e4 times wider than 2*stderr, onto the bounds.
        d = np.sqrt(np.clip(np.diag(cov), 0.0, None))  # 0 for fixed params
        safe = np.where(d > 0, d, 1.0)
        corr = cov / np.outer(safe, safe)
        L = np.linalg.cholesky(corr + 1e-10 * np.eye(n_params))
        perturbation = scale * d * (L @ rng.standard_normal(n_params))
        perturbed = params + perturbation

    except np.linalg.LinAlgError:
        stderr = np.sqrt(np.abs(np.diag(cov)))
        perturbation = scale * stderr * rng.standard_normal(n_params)
        perturbed = params + perturbation

    return _clip_to_bounds(perturbed, bounds)


def perturb_from_stderr(
    params: NDArray[np.float64],
    stderr: NDArray[np.float64],
    scale: float = 2.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None,
    rng: Seed = None
) -> NDArray[np.float64]:
    """
    Generate uncorrelated perturbation scaled by standard errors.
    """
    rng = np.random.default_rng(rng)
    stderr_safe = np.where(
        np.isfinite(stderr) & (stderr > 0),
        stderr,
        np.abs(params) * 0.1
    )

    perturbation = scale * stderr_safe * rng.standard_normal(len(params))
    perturbed = params + perturbation

    return _clip_to_bounds(perturbed, bounds)


def perturb_log_uniform(
    params: NDArray[np.float64],
    factor: float = 3.0,
    bounds: Optional[Tuple[NDArray, NDArray]] = None,
    rng: Seed = None
) -> NDArray[np.float64]:
    """
    Generate log-uniform perturbation (multiplicative).
    """
    rng = np.random.default_rng(rng)
    log_factor = np.log(factor)
    multipliers = np.exp(rng.uniform(-log_factor, log_factor, len(params)))
    perturbed = params * multipliers

    return _clip_to_bounds(perturbed, bounds)


def fit_circuit_multistart(
    circuit: Circuit,
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_restarts: int = 10,
    scale: float = 2.0,
    weighting: str = 'modulus',
    parallel: bool = False,
    max_workers: int = 4,
    use_analytic_jacobian: bool = True,
    rng: Seed = None
) -> Tuple[MultistartResult, NDArray[np.complex128]]:
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
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian (default: True)
    rng : None, int or numpy.random.Generator, optional
        Source of the restart perturbations. None draws fresh entropy; an int
        seed or a Generator makes the result reproducible. The global
        np.random state is never touched.

    Returns
    -------
    multistart_result : MultistartResult
        Multi-start optimization result with all diagnostics
        (plot: `plot_circuit_fit(frequencies, Z, multistart_result.best_result)`)
    Z_fit : ndarray
        Best fit impedance
    """
    all_results: List[FitResult] = []
    all_errors: List[Optional[float]] = []
    result_indices: List[int] = []  # start_idx aligned with all_results
    n_successful = 0
    diag_warnings: List[str] = []
    failed_errors = []  # Track exceptions from failed fits
    perturbation_method = 'covariance'

    jacobian_type = 'analytic' if use_analytic_jacobian else 'numeric'
    rng = np.random.default_rng(rng)

    # Step 1: Initial fit (on a copy too, so every result owns its circuit -
    # see run_single_fit; the caller's circuit is synced to the best fit at the end)
    try:
        result0, _ = fit_equivalent_circuit(
            frequencies, Z, deepcopy(circuit), weighting=weighting,
            use_analytic_jacobian=use_analytic_jacobian
        )
        all_results.append(result0)
        all_errors.append(result0.fit_error_rel)
        result_indices.append(1)  # initial fit is restart #1
        n_successful += 1
        initial_error = result0.fit_error_rel

    except Exception as e:
        raise RuntimeError(f"Multi-start failed: initial fit unsuccessful: {e}") from e

    # Extract bounds from circuit
    lower_bounds, upper_bounds = generate_simple_bounds(circuit.get_param_labels())
    bounds = (np.array(lower_bounds), np.array(upper_bounds))

    # Step 2: Generate perturbations and run additional fits
    def run_single_fit(start_idx: int, initial_params: NDArray) -> Optional[FitResult]:
        try:
            # Each restart fits its OWN copy of the circuit: fit_equivalent_circuit()
            # writes its parameters into the circuit object it is given, so a shared
            # circuit would be overwritten by every restart in turn (and, in parallel
            # mode, concurrently by several threads at once - a data race, since
            # impedance() reads the same parameters it is being written into).
            result, _ = fit_equivalent_circuit(
                frequencies, Z, deepcopy(circuit),
                weighting=weighting,
                initial_guess=list(initial_params),
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
                result0.params_opt, result0.cov, scale=scale, bounds=bounds, rng=rng
            )
            perturbation_method = 'covariance'
        elif not np.any(np.isinf(result0.params_stderr)):
            perturbed = perturb_from_stderr(
                result0.params_opt, result0.params_stderr, scale=scale, bounds=bounds, rng=rng
            )
            perturbation_method = 'stderr'
        else:
            perturbed = perturb_log_uniform(
                result0.params_opt, factor=3.0, bounds=bounds, rng=rng
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
                idx = futures[future]
                result = future.result()
                if result is not None:
                    all_results.append(result)
                    all_errors.append(result.fit_error_rel)
                    result_indices.append(idx)
                    n_successful += 1
                else:
                    all_errors.append(None)
    else:
        for idx, params in perturbations:
            result = run_single_fit(idx, params)
            if result is not None:
                all_results.append(result)
                all_errors.append(result.fit_error_rel)
                result_indices.append(idx)
                n_successful += 1
            else:
                all_errors.append(None)

    # Step 3: Find best result
    if not all_results:
        raise RuntimeError("Multi-start failed: no successful fits")

    # Select on the optimized objective (weighted RSS), not fit_error_rel: every
    # start minimizes RSS, and choosing on another metric could return a point
    # that is not the least-squares minimum - its covariance and AIC/BIC would
    # then describe the wrong fit (same reasoning as in diffevo.py).
    rss = [compute_information_criteria(
               Z, r.circuit.impedance(frequencies, list(r.params_opt)),
               weighting, r.n_free_params)[0]
           for r in all_results]
    best_pos = int(np.argmin(rss))
    best_result = all_results[best_pos]
    best_error = best_result.fit_error_rel
    # On the selection criterion: the initial fit is rss[0], and best_error
    # (fit_error_rel) can exceed initial_error when the two metrics disagree.
    improvement = (rss[0] - rss[best_pos]) / rss[0] * 100 if rss[0] > 0 else 0

    # Which start produced the best result. result_indices is aligned with
    # all_results, so this is robust to completion-order shuffling that
    # ThreadPoolExecutor/as_completed introduces in parallel mode.
    best_idx = result_indices[best_pos]

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

    # best_result.circuit already holds the best parameters (every fit ran on
    # its own circuit copy). The circuit the caller passed in was left at its
    # initial guesses, so sync it to the best fit as the single-fit and DE
    # paths do - consumers that read parameters from the circuit tree
    # (e.g. oxide capacitance/permittivity extraction) rely on it.
    circuit.update_params(list(best_result.params_opt))

    Z_fit_best = best_result.circuit.impedance(frequencies, list(best_result.params_opt))

    multistart_result = MultistartResult(
        best_result=best_result,
        all_results=all_results,
        n_starts=n_restarts,
        n_successful=n_successful,
        improvement=improvement,
        diagnostics=diagnostics
    )

    return multistart_result, Z_fit_best


__all__ = [
    'MultistartResult',
    'MultistartDiagnostics',
    'perturb_from_covariance',
    'perturb_from_stderr',
    'perturb_log_uniform',
    'fit_circuit_multistart',
]
