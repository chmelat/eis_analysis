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

# Points below this multiple of the lowest analysed frequency (the lowest
# decade) are where a flagged deviation is often sample drift rather than a bad
# point. The CLI says so under the per-point list and when --exclude-flagged
# leaves such points out; one constant, so both name the same points.
LOWEST_DECADE_FACTOR = 10.0

# A Lin-KK point far above the median residual is flagged even under the
# absolute threshold: at 0.1 % noise a point 2-5 % off is 20-50 sigma and pulls
# a fit (cold agent T13: 5 such points among 81, R1/R2 of two close arcs
# 142/38 instead of 100/80) while the threshold lets all of them through.
# Calibration (stress-test spectra, 6 families x 40, own noise 0.1-3 %, 6 % of
# the points given an extra error; detection by error size in units of the
# noise sigma, and spectra without bad points that get a false flag):
#
#   factor   10-20x   20-40x   >40x    false, 0.1 % / 0.3 % / 1 % noise
#     4       0.85     0.93    0.97         0.30 / 0.14 / 0.07
#     6       0.58     0.82    0.92         0.11 / 0.04 / 0.00
#     8       0.31     0.67    0.84         0.08 / 0.01 / 0.00
#    10       0.13     0.50    0.73         0.05 / 0.01 / 0.00
#
# 8: the false flags drop to the 0.1 % level, and those sit in the lowest
# fifth of the points (37 of 52 in the lowest 5 %), where the Voigt chain
# misses the end of a diffusion or oxide spectrum; the CLI note on deviations
# in the lowest decade covers most of them. Off on a sparse fit
# (KKResult.sparse_fit, 9 % of these spectra): no loss of detection.
# The chain absorbs part of each bad point, so detection stays below what the
# factor alone would give. Lin-KK only: the Z-HIT residual is a smoothed,
# integrated quantity whose scatter this was not measured on.
_ISOLATED_POINT_FACTOR = 8.0


@dataclass
class OutlierPoint:
    """
    A single point whose residual exceeds the flagging threshold.

    Attributes
    ----------
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
        exceeds the threshold, so per-point flagging is not meaningful), and
        Lin-KK on a capacitive low-frequency end without the series C
        (KKResult.capacitive_end: its residuals describe the model) when the
        rerun with the series C is missing; with it, its residuals are used.
    """
    points: List[OutlierPoint] = field(default_factory=list)
    skipped: List[str] = field(default_factory=list)

    def frequencies_flagged_by(self, mode: str = 'kk') -> List[float]:
        """
        Frequencies [Hz] of the flagged points, ascending, by flagging method.

        Parameters
        ----------
        mode : {'kk', 'both', 'any'}
            'kk': flagged by Lin-KK (alone or with Z-HIT) - the default for
            excluding points, because Z-HIT smears one bad point into its
            neighbours (it rebuilds |Z| from the integrated phase; cold agent
            T14: 3 of the 5 listed points were neighbours flagged by Z-HIT
            alone). 'both': flagged by both methods. 'any': every listed point.
        """
        accept = {'kk': ("KK", "KK+Z-HIT"), 'both': ("KK+Z-HIT",),
                  'any': ("KK", "Z-HIT", "KK+Z-HIT")}
        if mode not in accept:
            raise ValueError(f"mode must be one of {sorted(accept)}, got '{mode}'")
        return sorted(p.frequency for p in self.points if p.methods in accept[mode])


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
    if getattr(result, 'capacitive_end', False):
        # The residuals describe the Voigt chain's misfit of the end; the
        # rerun with the series C describes the data (and locates drift)
        result = getattr(result, 'with_series_c', None)
        if result is None:
            skipped.append(label)
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

    if not np.any(np.isfinite(r)):
        logger.debug(f"{label}: residuals not finite, skipping outlier check")
        return None

    # The global guard (see find_outliers) takes the median, not the mean: a
    # handful of genuine outliers pulls a mean over the threshold and would
    # switch the report off exactly when it has something to say - 5 points at
    # 40% among 30 at 0.5% give a mean of 6.1% but a median of 0.5%. The median
    # also keeps one non-finite residual (Z-HIT divides by |Z| without a floor,
    # unlike KK) from disabling the whole method.
    if np.nanmedian(r) > max_residual:
        skipped.append(label)
        return None

    return r


def _flagged(
    r: Optional[NDArray[np.float64]],
    max_residual: float,
    n_points: int,
    relative: bool = False
) -> NDArray[np.bool_]:
    """
    Per-point flag mask for one method: over the absolute threshold AND over
    the method's own baseline (see _METHOD_BASELINE_FACTOR). With `relative`,
    also over _ISOLATED_POINT_FACTOR x the median alone.

    Non-finite residuals never flag a point (NaN comparisons are False).
    """
    if r is None:
        return np.zeros(n_points, dtype=bool)
    median = np.nanmedian(r)
    flagged = (r > max_residual) & (r > _METHOD_BASELINE_FACTOR * median)
    if relative:
        flagged |= r > _ISOLATED_POINT_FACTOR * median
    return flagged


def _relative_applies(kk_result) -> bool:
    """
    Does the relative criterion apply to this Lin-KK run? Not when its
    median residual is below the noise floor (noise-free synthetic data) or
    the chain has more elements than half the points (KKResult.sparse_fit):
    then the median says nothing about the noise. The run is the one whose
    residuals are used: with_series_c on a capacitive end. A stand-in without
    sparse_fit (not a KKResult) gets the absolute criterion only.
    """
    from ..fitting.config import NOISE_FLOOR_PERCENT
    if getattr(kk_result, 'capacitive_end', False):
        kk_result = kk_result.with_series_c
    if kk_result is None or getattr(kk_result, 'sparse_fit', True):
        return False
    return bool(kk_result.noise_estimate_median >= NOISE_FLOOR_PERCENT)


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
        Absolute flagging threshold [%] (default: 5.0, the same number as
        KK_RESIDUAL_THRESHOLD). Note this metric combines both residual
        components, so it runs up to sqrt(2) larger than the per-point
        max(|res_real|, |res_imag|) that `KKResult.is_valid` counts - the two
        are not the same test.

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

    A Lin-KK point is also flagged below the threshold when its residual is
    over _ISOLATED_POINT_FACTOR (8) x the median Lin-KK residual: on a
    low-noise spectrum a point a few percent off is many sigma away and pulls
    a fit. Not when the median is below the noise floor
    (NOISE_FLOOR_PERCENT, noise-free synthetic data), not on a sparse fit
    (KKResult.sparse_fit), not for Z-HIT.

    A method whose MEDIAN residual is already over the threshold is dropped
    into `skipped` instead: over half its points would be flagged, which is a
    global failure rather than a list of local ones.

    Non-finite residuals never flag a point (NaN comparisons are False).
    """
    report = OutlierReport()
    n = len(frequencies)
    if n == 0:
        return report

    r_kk = _usable_residuals(kk_result, n, max_residual, "KK", report.skipped)
    r_zhit = _usable_residuals(zhit_result, n, max_residual, "Z-HIT", report.skipped)
    if r_kk is None and r_zhit is None:
        return report

    bad_kk = _flagged(r_kk, max_residual, n, relative=_relative_applies(kk_result))
    bad_zhit = _flagged(r_zhit, max_residual, n)

    # `worst` counts only the methods that actually flagged the point.
    worst = np.maximum(np.where(bad_kk, r_kk, 0.0) if r_kk is not None else 0.0,
                       np.where(bad_zhit, r_zhit, 0.0) if r_zhit is not None else 0.0)

    for i in np.flatnonzero(bad_kk | bad_zhit):
        methods = [m for m, bad in (("KK", bad_kk[i]), ("Z-HIT", bad_zhit[i])) if bad]
        report.points.append(OutlierPoint(
            frequency=float(frequencies[i]),
            residual_kk=float(r_kk[i]) if r_kk is not None else None,
            residual_zhit=float(r_zhit[i]) if r_zhit is not None else None,
            worst=float(worst[i]),
            methods="+".join(methods)
        ))

    report.points.sort(key=lambda p: p.worst, reverse=True)
    return report
