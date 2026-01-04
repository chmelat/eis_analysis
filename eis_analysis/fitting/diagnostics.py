"""
Fit diagnostics and quality assessment for circuit fitting.

Provides parameter diagnostics (uncertainty, correlation, extreme values)
and fit quality metrics computation.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import List, Tuple, Optional
from numpy.typing import NDArray

from .config import FIT_QUALITY_EXCELLENT_ERROR, FIT_QUALITY_GOOD_ERROR

logger = logging.getLogger(__name__)


def compute_weights(Z: NDArray[np.complex128], weighting: str) -> NDArray[np.float64]:
    """
    Compute weights based on weighting type.

    Parameters
    ----------
    Z : ndarray of complex
        Impedance data
    weighting : str
        Type of weighting: 'uniform', 'sqrt', 'proportional', or 'square'

    Returns
    -------
    weights : ndarray of float
        Normalized weights (mean = 1)
    """
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    if weighting == 'uniform':
        weights = np.ones_like(Z_mag)
    elif weighting == 'sqrt':
        weights = 1.0 / np.sqrt(Z_mag_safe)
    elif weighting == 'proportional':
        weights = 1.0 / Z_mag_safe
    elif weighting == 'square':
        weights = Z_mag ** 2
    else:
        weights = np.ones_like(Z_mag)

    return weights / np.mean(weights)


def check_parameter_diagnostics(
    params_opt: NDArray[np.float64],
    params_stderr: NDArray[np.float64],
    cov: Optional[NDArray[np.float64]],
    param_labels: Optional[List[str]]
) -> None:
    """
    Check for high uncertainties, correlations, and extreme values.

    Parameters
    ----------
    params_opt : ndarray
        Optimized parameters
    params_stderr : ndarray
        Standard errors of parameters
    cov : ndarray or None
        Covariance matrix
    param_labels : list of str or None
        Parameter labels for logging
    """
    # Check high uncertainties
    if not np.all(np.isinf(params_stderr)):
        relative_errors = params_stderr / np.maximum(np.abs(params_opt), 1e-15) * 100
        high_uncertainty_mask = relative_errors > 100

        if np.any(high_uncertainty_mask):
            logger.warning("=" * 50)
            logger.warning("WARNING: Very high uncertainty for some parameters!")
            for i, (param, stderr, rel_err) in enumerate(zip(params_opt, params_stderr, relative_errors)):
                if high_uncertainty_mask[i]:
                    logger.warning(f"  Parameter {i}: {param:.3e} +/- {stderr:.3e} ({rel_err:.1f}%)")
            logger.warning("=" * 50)

    # Check correlations
    if cov is not None:
        n_params = len(params_stderr)
        high_corr_pairs = []

        for i in range(n_params):
            for j in range(i + 1, n_params):
                if params_stderr[i] == 0 or params_stderr[j] == 0:
                    continue
                if np.isinf(params_stderr[i]) or np.isinf(params_stderr[j]):
                    continue

                corr_ij = cov[i, j] / (params_stderr[i] * params_stderr[j])
                if abs(corr_ij) > 0.95:
                    high_corr_pairs.append((i, j, corr_ij))

        if high_corr_pairs:
            logger.warning("=" * 50)
            logger.warning("WARNING: High correlation between parameters!")
            for i, j, corr_ij in high_corr_pairs:
                if param_labels and i < len(param_labels) and j < len(param_labels):
                    logger.warning(f"    {param_labels[i]} <-> {param_labels[j]}: corr = {corr_ij:+.3f}")
                else:
                    logger.warning(f"    p[{i}] <-> p[{j}]: corr = {corr_ij:+.3f}")
            logger.warning("  Recommendation: Try simpler circuit or fix a parameter")
            logger.warning("=" * 50)

    # Check extreme values
    extreme_small = np.abs(params_opt) < 1e-15
    extreme_large = np.abs(params_opt) > 1e15

    if np.any(extreme_small) or np.any(extreme_large):
        logger.warning("=" * 50)
        logger.warning("WARNING: Extreme parameter values!")
        if np.any(extreme_small):
            logger.warning(f"  Very small: {params_opt[extreme_small]}")
        if np.any(extreme_large):
            logger.warning(f"  Very large: {params_opt[extreme_large]}")
        logger.warning("=" * 50)


def compute_fit_metrics(
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128],
    weighting: str
) -> Tuple[float, float, str]:
    """
    Compute fit error metrics and quality assessment.

    Parameters
    ----------
    Z : ndarray of complex
        Measured impedance data
    Z_fit : ndarray of complex
        Fitted impedance
    weighting : str
        Weighting type used in fitting

    Returns
    -------
    fit_error_rel : float
        Weighted relative error [%]
    fit_error_abs : float
        Mean absolute error [Ohm]
    quality : str
        Quality assessment: 'excellent', 'good', 'acceptable', 'poor'
    """
    weights = compute_weights(Z, weighting)
    Z_mag_safe = np.maximum(np.abs(Z), 1e-15)
    relative_errors = np.abs(Z - Z_fit) / Z_mag_safe

    fit_error_rel = np.sum(weights * relative_errors) / np.sum(weights) * 100
    fit_error_abs = np.mean(np.abs(Z - Z_fit))

    # Log unweighted vs weighted difference if significant
    fit_error_rel_unweighted = np.mean(relative_errors) * 100
    if abs(fit_error_rel_unweighted - fit_error_rel) > 10:
        logger.info(f"  Note: Unweighted error {fit_error_rel_unweighted:.1f}%, weighted {fit_error_rel:.2f}%")

    # Quality assessment
    if fit_error_rel < FIT_QUALITY_EXCELLENT_ERROR:
        quality = 'excellent'
    elif fit_error_rel < FIT_QUALITY_GOOD_ERROR:
        quality = 'good'
    elif fit_error_rel < FIT_QUALITY_GOOD_ERROR * 2:
        quality = 'acceptable'
    else:
        quality = 'poor'

    return fit_error_rel, fit_error_abs, quality


def log_fit_results(
    params_opt: NDArray[np.float64],
    params_stderr: NDArray[np.float64],
    ci_low: NDArray[np.float64],
    ci_high: NDArray[np.float64],
    fit_error_rel: float,
    fit_error_abs: float,
    quality: str,
    param_labels: Optional[List[str]]
) -> None:
    """
    Log fit results to console.

    Parameters
    ----------
    params_opt : ndarray
        Optimized parameters
    params_stderr : ndarray
        Standard errors
    ci_low : ndarray
        Lower 95% CI bounds
    ci_high : ndarray
        Upper 95% CI bounds
    fit_error_rel : float
        Relative fit error [%]
    fit_error_abs : float
        Absolute fit error [Ohm]
    quality : str
        Fit quality assessment
    param_labels : list of str or None
        Parameter labels
    """
    logger.info("")
    logger.info("Fit results:")
    logger.info("  Parameters:")

    for i, (param, stderr) in enumerate(zip(params_opt, params_stderr)):
        label = param_labels[i] if param_labels and i < len(param_labels) else f"p[{i}]"
        if np.isinf(ci_low[i]) or np.isinf(ci_high[i]):
            logger.info(f"    {label:5s} = {param:.6e} +/- {stderr:.6e}")
        else:
            logger.info(f"    {label:5s} = {param:.2e} +/- {stderr:.2e}  [95% CI: {ci_low[i]:.2e}, {ci_high[i]:.2e}]")

    logger.info(f"  Fit error: {fit_error_rel:.2f}% (rel), {fit_error_abs:.2f} Ohm (abs)")

    if quality == 'excellent':
        logger.info(f"  Quality: Excellent (<{FIT_QUALITY_EXCELLENT_ERROR}%)")
    elif quality == 'good':
        logger.info(f"  Quality: Good (<{FIT_QUALITY_GOOD_ERROR}%)")
    elif quality == 'acceptable':
        logger.warning("  Quality: Acceptable (consider checking the model)")
    else:
        logger.warning("  Quality: POOR! Model does not fit the data")


__all__ = [
    'compute_weights',
    'check_parameter_diagnostics',
    'compute_fit_metrics',
    'log_fit_results',
]
