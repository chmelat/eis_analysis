"""
Per-point outlier flagging from KK / Z-HIT residuals.

Both validation methods already compute a residual for every measured point;
this module turns those arrays into a list of individual suspicious points.

The detection is purely diagnostic - nothing is removed or down-weighted.
That is deliberate: in EIS most low-frequency deviations are sample drift
(a real property of the measurement), not bad points, and discarding them
would hide a non-stationary experiment rather than fix it.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# A flagged point must also stand out from the method's OWN reconstruction
# error, not just from the absolute threshold. A residual is the sum of the
# data error and the method's error, and the latter is not the same for both:
# Z-HIT rebuilds |Z| from the numerically differentiated phase, so on noisy
# data its residual baseline reaches ~3% on spectra containing no defect at
# all, while a Lin-KK chain fits an order of magnitude tighter. Without this
# factor the list fills up with Z-HIT artifacts.
#
# Calibration (25 noise realizations of the 1%-noise synthetic spectrum,
# 1750 defect-free points, vs. detection of a single injected spike):
#
#   factor   false alarms   5% spike   8% spike   10% spike
#     none      8.2%          8/25      24/25       25/25
#      3        3.4%          7/25      24/25       25/25
#      4        1.9%          5/25      23/25       25/25
#      5        1.2%          2/25      12/25       23/25
#
# 4 is the knee: it removes three quarters of the false alarms at no real cost
# in sensitivity, while 5 starts missing half of the 8% spikes. On measured
# spectra (EISPOT-w2-10/11/13, median residual 0.3-0.8%) the factor changes
# nothing at all - there the absolute threshold stays binding. It only bites
# where a method cannot reconstruct the spectrum better than the threshold.
#
# Rejected alternative: scaling by the noise_estimate the validations already
# report. It is derived from pseudo-chi^2, so it absorbs genuine misfit and
# raises the bar exactly on the spectra that have real problems (on w2-13 it
# lands at 29%, flagging 1 point instead of 17).
_METHOD_BASELINE_FACTOR = 4.0


@dataclass
class OutlierPoint:
    """
    A single point whose residual exceeds the flagging threshold.

    Attributes
    ----------
    index : int
        Position in the input frequency array.
    frequency : float
        Frequency of the point [Hz]
    residual_kk : float or None
        Relative deviation from the Lin-KK fit [%]. None if KK did not
        contribute (not run, failed, or excluded by the global guard).
    residual_zhit : float or None
        Relative deviation from the Z-HIT reconstruction [%]. None if Z-HIT
        did not contribute.
    worst : float
        Larger of the contributing residuals [%]
    methods : str
        Which methods flagged the point: "KK", "Z-HIT" or "KK+Z-HIT"
    """
    index: int
    frequency: float
    residual_kk: Optional[float]
    residual_zhit: Optional[float]
    worst: float
    methods: str


@dataclass
class OutlierReport:
    """
    Result of per-point outlier flagging.

    Attributes
    ----------
    points : list of OutlierPoint
        Flagged points, ordered by descending `worst`.
    skipped : list of str
        Methods excluded by the global guard (their own mean residual already
        exceeds the threshold, so per-point flagging is not meaningful).
    max_residual : float
        Threshold used [%]
    """
    points: List[OutlierPoint] = field(default_factory=list)
    skipped: List[str] = field(default_factory=list)
    max_residual: float = 5.0


def residual_percent(
    residuals_real: NDArray[np.float64],
    residuals_imag: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Relative deviation from a reconstruction, per point.

    Parameters
    ----------
    residuals_real, residuals_imag : ndarray
        Real and imaginary residuals normalized by |Z| (fractions), as
        produced by both KK and Z-HIT validation.

    Returns
    -------
    ndarray
        |Z - Z_fit| / |Z| in percent.

    Notes
    -----
    One definition serves both methods. For Z-HIT the reconstruction keeps
    the measured phase (Z_fit = |Z_recon| * exp(i*phi)), so the two residual
    components share a phase factor and this reduces exactly to the magnitude
    residual |residuals_mag|. For Lin-KK it is the full complex deviation
    from the Voigt fit.
    """
    return np.hypot(residuals_real, residuals_imag) * 100.0


def _usable_residuals(
    result,
    n_points: int,
    max_residual: float,
    label: str,
    skipped: List[str]
) -> Optional[NDArray[np.float64]]:
    """
    Per-point residuals [%] of one validation method, or None if unusable.

    Returns None when the method did not run, failed, returned arrays that do
    not pair with the frequency grid, or was excluded by the global guard.
    """
    if result is None:
        return None

    res_real = getattr(result, 'residuals_real', None)
    res_imag = getattr(result, 'residuals_imag', None)
    if res_real is None or res_imag is None:
        return None
    if len(res_real) != n_points or len(res_imag) != n_points:
        # Failed run (Z-HIT returns empty arrays on reconstruction failure)
        # or a caller pairing residuals with the wrong frequency grid.
        logger.debug(f"{label}: residual length != {n_points}, skipping outlier check")
        return None

    r = residual_percent(res_real, res_imag)

    # Global guard: a method whose own mean residual already exceeds the
    # threshold is describing a spectrum that fails as a whole, not a few bad
    # points. Z-HIT in particular carries a free global offset (the anchor at
    # ref_freq), and a shifted offset lifts every residual at once - without
    # this guard the entire spectrum would be listed point by point.
    mean_r = np.nanmean(r) if np.any(np.isfinite(r)) else np.nan
    if not np.isfinite(mean_r):
        logger.debug(f"{label}: residuals not finite, skipping outlier check")
        return None
    if mean_r > max_residual:
        skipped.append(label)
        return None

    return r


def _flagged(
    r: Optional[NDArray[np.float64]],
    max_residual: float,
    n_points: int
) -> NDArray[np.bool_]:
    """
    Per-point flag mask for one method: over the absolute threshold AND over
    the method's own baseline (see _METHOD_BASELINE_FACTOR).

    Non-finite residuals never flag a point (NaN comparisons are False).
    """
    if r is None:
        return np.zeros(n_points, dtype=bool)
    baseline = _METHOD_BASELINE_FACTOR * np.nanmedian(r)
    return (r > max_residual) & (r > baseline)


def find_outliers(
    frequencies: NDArray[np.float64],
    kk_result=None,
    zhit_result=None,
    max_residual: float = 5.0
) -> OutlierReport:
    """
    Flag individual points whose KK or Z-HIT residual exceeds a threshold.

    Parameters
    ----------
    frequencies : ndarray
        Measured frequencies [Hz]
    kk_result : KKResult, optional
        Result of Kramers-Kronig validation. None if it did not run.
    zhit_result : ZHITResult, optional
        Result of Z-HIT validation. None if it did not run.
    max_residual : float, optional
        Flagging threshold [%] (default: 5.0, matching the pass/fail
        threshold used for the mean residual).

    Returns
    -------
    OutlierReport
        Flagged points (worst first) and the methods excluded by the guard.

    Notes
    -----
    A point is flagged when EITHER method exceeds the threshold. The two are
    complementary: Lin-KK fits both impedance components and is sensitive to
    phase errors, while Z-HIT reconstructs the magnitude from the phase and
    is sensitive to magnitude errors and drift.

    The threshold is necessary but not sufficient: the point must also exceed
    the method's own residual baseline (see _METHOD_BASELINE_FACTOR), which
    keeps a method that reconstructs the spectrum poorly from filling the list
    with its own artifacts.

    Non-finite residuals never flag a point (NaN comparisons are False).
    """
    report = OutlierReport(max_residual=max_residual)
    n = len(frequencies)
    if n == 0:
        return report

    r_kk = _usable_residuals(kk_result, n, max_residual, "KK", report.skipped)
    r_zhit = _usable_residuals(zhit_result, n, max_residual, "Z-HIT", report.skipped)
    if r_kk is None and r_zhit is None:
        return report

    bad_kk = _flagged(r_kk, max_residual, n)
    bad_zhit = _flagged(r_zhit, max_residual, n)

    # `worst` counts only the methods that actually flagged the point, so a
    # method that stayed under the threshold cannot raise it.
    worst_kk = np.where(bad_kk, r_kk, 0.0) if r_kk is not None else np.zeros(n)
    worst_zhit = np.where(bad_zhit, r_zhit, 0.0) if r_zhit is not None else np.zeros(n)

    for i in np.flatnonzero(bad_kk | bad_zhit):
        methods = [m for m, bad in (("KK", bad_kk[i]), ("Z-HIT", bad_zhit[i])) if bad]
        report.points.append(OutlierPoint(
            index=int(i),
            frequency=float(frequencies[i]),
            residual_kk=float(r_kk[i]) if r_kk is not None else None,
            residual_zhit=float(r_zhit[i]) if r_zhit is not None else None,
            worst=float(max(worst_kk[i], worst_zhit[i])),
            methods="+".join(methods)
        ))

    report.points.sort(key=lambda p: p.worst, reverse=True)
    return report
