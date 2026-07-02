"""
Fit diagnostics and quality assessment for circuit fitting.

Provides weight computation and fit quality metrics.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import Tuple
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
        Type of weighting: 'uniform', 'sqrt', 'modulus', or 'proportional'
        - 'modulus': w = 1/|Z| (DEFAULT, Lin-KK standard)
        - 'proportional': w = 1/|Z|^2 (strong low-Z emphasis)

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
    elif weighting == 'modulus':
        weights = 1.0 / Z_mag_safe
    elif weighting == 'proportional':
        weights = 1.0 / (Z_mag_safe ** 2)
    else:
        logger.warning(f"Unknown weighting '{weighting}', using uniform weights")
        weights = np.ones_like(Z_mag)

    return weights / np.mean(weights)


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
        Weighting-consistent relative error [%]:
        ``sum(w_i * |Z_i - Z_fit_i|) / sum(w_i * |Z_i|) * 100``. The weight is
        applied once (to both residual and magnitude), so it is not
        double-counted with the 1/|Z| of a relative error. For modulus
        weighting this equals the mean relative error ``mean(|dZ|/|Z|)``.
    fit_error_abs : float
        Mean absolute error [Ohm]
    quality : str
        Quality assessment: 'excellent', 'good', 'acceptable', 'poor'
    """
    weights = compute_weights(Z, weighting)
    Z_mag_safe = np.maximum(np.abs(Z), 1e-15)
    abs_errors = np.abs(Z - Z_fit)
    relative_errors = abs_errors / Z_mag_safe

    # Weighting-consistent relative error: the weight is applied once, to both
    # the residual and the magnitude, so it is not double-counted with the
    # 1/|Z| that already defines a relative error. For modulus weighting
    # (w = 1/|Z|) this reduces to the mean relative error mean(|dZ|/|Z|).
    fit_error_rel = np.sum(weights * abs_errors) / np.sum(weights * Z_mag_safe) * 100
    fit_error_abs = np.mean(abs_errors)

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


__all__ = [
    'compute_weights',
    'compute_fit_metrics',
]
