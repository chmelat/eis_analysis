"""
Fit diagnostics and quality assessment for circuit fitting.

Provides weight computation and fit quality metrics.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import Optional, Sequence, Tuple
from numpy.typing import NDArray

from .config import FIT_QUALITY_EXCELLENT_ERROR, FIT_QUALITY_GOOD_ERROR
from .jacobian import circuit_jacobian

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


def compute_information_criteria(
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128],
    weighting: str,
    n_free_params: int
) -> Tuple[float, float, float]:
    """
    Compute AIC and BIC for a fitted circuit model.

    Adding an element to a circuit almost always lowers the residual, so the
    fit error alone cannot tell an improvement apart from fitting the noise.
    Both criteria answer that by charging for each free parameter: the model
    with the lowest value is the one the data actually supports.

    Only *differences* between models carry meaning - the absolute value
    depends on an additive constant that is dropped here. The conventional
    reading of a difference is: below 2 the models are indistinguishable,
    4-7 is a noticeable difference, above 10 is decisive.

    BIC charges ``ln(n)`` per parameter against AIC's 2, so for a typical
    spectrum (n = 160 residuals, ln(n) = 5.1) it penalises complexity about
    2.5x harder and prefers simpler circuits. Reporting both is deliberate:
    where they disagree, the data does not settle the extra element.

    The small-sample correction AICc adds ``2k(k+1)/(n-k-1)``, which stays
    under 0.4 for a typical spectrum (n = 160, k = 5). It is not applied.

    Parameters
    ----------
    Z : ndarray of complex
        Measured impedance [Ohm]
    Z_fit : ndarray of complex
        Fitted impedance [Ohm]
    weighting : str
        Weighting used for the fit; see compute_weights
    n_free_params : int
        Number of freely optimized parameters (fixed parameters excluded)

    Returns
    -------
    rss : float
        Weighted residual sum of squares - the quantity the optimizer
        minimizes
    aic : float
        Akaike information criterion, ``n*ln(RSS/n) + 2k``
    bic : float
        Bayesian information criterion, ``n*ln(RSS/n) + k*ln(n)``

    Notes
    -----
    Values are comparable ONLY between models fitted to the same data with
    the same weighting. Change the frequency range or the weighting between
    two candidates and the comparison becomes meaningless, with nothing to
    signal it.

    The real and imaginary parts are separate residuals, so ``n = 2*len(Z)``.
    """
    k = int(n_free_params)
    # A fitted model always has at least one free parameter, so k <= 0 means
    # the caller passed a FitResult whose n_free_params was never populated.
    # Scoring it would silently charge no complexity penalty at all, handing
    # that model every comparison it enters.
    if k <= 0:
        raise ValueError(
            f"n_free_params must be positive, got {k}. A FitResult built "
            "outside the standard optimizers may not populate it."
        )

    weights = compute_weights(Z, weighting)
    residuals_re = weights * (Z.real - Z_fit.real)
    residuals_im = weights * (Z.imag - Z_fit.imag)
    rss = float(np.sum(residuals_re ** 2) + np.sum(residuals_im ** 2))

    n = 2 * len(Z)

    # A perfect fit (RSS = 0) makes ln(RSS/n) diverge to -inf. That is the
    # correct limit - such a model wins every comparison - but -inf poisons
    # the delta arithmetic that follows, so clamp to the smallest positive
    # normal instead and let the ranking stay finite.
    rss_safe = max(rss, np.finfo(float).tiny)

    log_likelihood_term = n * np.log(rss_safe / n)
    aic = log_likelihood_term + 2.0 * k
    bic = log_likelihood_term + k * np.log(n)

    return rss, float(aic), float(bic)


def compute_significance(
    circuit,
    frequencies: NDArray[np.float64],
    params: Sequence[float]
) -> Optional[NDArray[np.float64]]:
    """
    Sensitivity of the network impedance to each parameter (Zahner significance).

    S_i = max_n |d ln|Z_n| / d ln P_i|

    One number per parameter, answering a different question than the standard
    error: not "how precisely is this parameter determined?" but "does this
    parameter matter in this frequency window at all?". A large standard error
    conflates a parameter that is irrelevant with one that is merely correlated
    with another; the significance separates them.

    Parameters
    ----------
    circuit : CircuitElement or CompositeCircuit
        The fitted circuit
    frequencies : ndarray of float
        Measurement frequencies [Hz] - the significance is a property of the
        model *in the measured window*, not of the model in the abstract
    params : sequence of float
        All circuit parameters (fixed ones included), as in FitResult.params_opt

    Returns
    -------
    significance : ndarray of float or None
        One value per parameter, aligned with `params`. None if the circuit
        contains an element with no analytic derivative.

    Notes
    -----
    Z is the impedance of the *whole network*, not of the individual element,
    and P is a scalar parameter, not an element - a CPE contributes two.
    S ~ 1 means the parameter dominates the impedance somewhere in the window,
    below SIGNIFICANCE_NEGLIGIBLE the element may be dropped, and above 1 is
    normal for a parameter entering non-linearly.

    One deviation from the source: it maximises the signed quantity, this
    maximises the absolute value.

    Full derivation, the bound on a linear element, the edge cases and that
    deviation: doc/WEIGHTING_AND_STATISTICS.md section 3.7.

    References
    ----------
    Zahner Analysis manual (11/2023), section 2.2.2 "Significance".
    """
    # circuit_jacobian is the right source: unweighted, un-negated dZ/dp with a
    # column for every parameter. The optimizer's Jacobian is none of those.
    try:
        Z, dZ = circuit_jacobian(circuit, frequencies, list(params))
    except NotImplementedError:
        return None

    # d|Z|/dp = Re(conj(Z) * dZ/dp) / |Z|, so the whole expression is
    # Re(conj(Z) * dZ/dp) * P / |Z|^2. Clamp |Z|^2 the way compute_weights
    # clamps |Z|: a data point at exactly zero impedance is not physical, but
    # it must not turn the diagnostic into a division by zero.
    Z_mag2 = np.maximum(np.abs(Z) ** 2, 1e-30)
    relative = (Z.conj()[:, np.newaxis] * dZ).real * np.asarray(params, dtype=float)
    return np.max(np.abs(relative / Z_mag2[:, np.newaxis]), axis=0)


__all__ = [
    'compute_weights',
    'compute_fit_metrics',
    'compute_information_criteria',
    'compute_significance',
]
