"""
Local CPE exponent n(f): where and how a spectrum departs from one CPE.

For a layer whose admittance is Y = G + Q(jw)^n + jwC, the real part of the
admittance (series R_inf and L removed) is G + Q w^n cos(n pi/2): C does not
enter it, and where the CPE dominates its log-log slope is n. A single CPE
therefore shows a constant slope; any variation with frequency is a second
process or a different kind of element. The map needs no circuit fit, only
R_inf and L.

CONTRACT

Callers. The CLI with --local-exponent (cli/handlers/local_exponent.py)
passes R_inf, L and R_inf_range from estimate_rinf: finite, R_inf and L
>= 0, the range holding R_inf. Its fallback R_inf can be Re(Z) of a
measured point (hf_bound), so with L = 0 and Im(Z) exactly 0 there,
Z - R_inf is exactly 0 at that point. ValueError is printed as "Local
exponent failed". The data come from the loaders: finite, f > 0, any
order, duplicate frequencies possible. plot_local_exponent reads
frequencies, n, n_uncertainty, valid and span. Public API
(eis_analysis.local_exponent): any caller, nothing guaranteed beyond the
input rules below. From other modules it uses validate_eis_data
(ValueError for unequal lengths, no points, f <= 0, NaN or inf in f or Z,
Z = 0 at a point, although Y = 1/(Z - R_inf - j*omega*L) would be finite
there)
and RINF_REL_STDERR_MAX.

local_exponent(frequencies, Z, R_inf, L=0, R_inf_range=None)
-> LocalExponentResult.
  Input: ValueError for what validate_eis_data rejects, for fewer than
  LOCAL_EXPONENT_MIN_POINTS points, for R_inf or L not finite, and for
  R_inf_range not finite or not holding R_inf. The caller passes 1-D
  arrays, R_inf and L as real numbers and R_inf_range as a pair of them;
  other types are not checked. Any finite R_inf and L are subtracted as
  given, negative ones too. Default range: R_inf -+
  LOCAL_EXPONENT_RINF_REL x |R_inf| (zero width for R_inf = 0: only the
  regression stderr counts).
  Y = 1 / (Z - R_inf - j*omega*L), omega = 2 pi f. A point has an exponent
  only where Re Y is finite and > 0 (not where Z - R_inf - j*omega*L = 0).
  n: at a point with an exponent, the least-squares slope of ln Re Y
  against ln omega over the points with an exponent within
  |log10 f - log10 f_i| <= LOCAL_EXPONENT_WINDOW_DECADES / 2 (inclusive,
  as computed in floating point). NaN at a point without an exponent, and
  where such points in the window have fewer than
  LOCAL_EXPONENT_MIN_POINTS distinct frequencies (repeated sweeps count
  once).
  n_uncertainty: hypot of the slope's stderr (residual variance over
  k - 2 for k window points) and the sensitivity to R_inf, the larger
  |n(R) - n| for R at the two ends of R_inf_range. NaN where n is NaN at
  R_inf or at either end of the range.
  valid: n finite and n_uncertainty <= LOCAL_EXPONENT_UNCERTAINTY_MAX.
  frequencies, n, n_uncertainty, valid: one entry per input point, sorted
  by ascending frequency; points at the same frequency in input order.
  n_min, f_n_min, n_max, f_n_max: extremes of n over the valid points and
  their frequencies, at a tie the lowest frequency. All NaN when no point
  is valid, with a warning that n is determined nowhere. A span above
  LOCAL_EXPONENT_SPAN_WARN gives a warning naming both extremes and their
  frequencies. No other warnings.
  R_inf, L, R_inf_range: the values used, as floats; window_decades and
  uncertainty_max: the constants.

What must hold:
  - The inputs are not changed.
  - Raises only ValueError, and only for the input cases above; no
    RuntimeWarning.
  - n_uncertainty >= 0 where finite; valid implies n finite.
  - The uncertainty covers the ends of R_inf_range: |n(R) - n| <=
    n_uncertainty at a valid point for R at either end. Between them n(R)
    is not evaluated; it can leave that band where a window point's Re Z
    lies inside the range, since its Re Y goes to 0 as R approaches Re Z
    (constructed case: n 0.80 +- 0.014, 0.69 at R 0.1 % below Re Z; none
    in random spectra through estimate_rinf).
  - Unit invariance: Z, R_inf, L and R_inf_range scaled by c > 0 leave n,
    n_uncertainty and valid unchanged (to rounding). Input frequencies
    scaled by c with L scaled by 1/c leave them unchanged too and scale
    result.frequencies by c, except where a point lies on a window edge
    to rounding (an exact log grid at
    10 points per decade): whether it is in depends on the rounding of
    log10 and changes with c. Measured frequencies, rounded by the
    instrument, are not affected.
  - The order of the points changes nothing but the order of points at
    the same frequency (and n by rounding, ~1e-16, where duplicates
    reorder a window's sums).
  - One CPE in parallel with C, with the exact R_inf and L subtracted:
    n equals the CPE exponent at every point with a window (to rounding).
  - The result is deterministic.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..utils.validation import validate_eis_data
from ..rinf_estimation.estimate import RINF_REL_STDERR_MAX

# Width of the sliding log-log regression. One decade is ~10 points at the
# usual 10 points/decade; on the M136 ZrO2 spectra 0.6 decade gave the same
# n(f) with up to twice the scatter, so narrower resolves nothing more.
LOCAL_EXPONENT_WINDOW_DECADES = 1.0

# Fewest distinct frequencies a regression window may have: two for the line,
# two more so the residual scatter (its stderr) means something. Also the
# fewest points the input may have.
LOCAL_EXPONENT_MIN_POINTS = 4

# A point counts as determined when its uncertainty (regression stderr and the
# R_inf sensitivity below, combined) is at most this. 0.02 keeps the span
# warning (LOCAL_EXPONENT_SPAN_WARN = 0.1) at 5x the largest point uncertainty,
# so noise alone cannot raise it.
LOCAL_EXPONENT_UNCERTAINTY_MAX = 0.02

# Relative shift of R_inf used to test each point's sensitivity to it. The
# formal stderr of R_inf understates the error under model mismatch (0.1-0.3 %
# on the M136 spectra, where n still drifted near f_max), so the threshold at
# which R_inf counts as determined at all is used instead.
LOCAL_EXPONENT_RINF_REL = RINF_REL_STDERR_MAX

# Span of n over the determined points above which one CPE cannot describe the
# dispersion. Typical fitted CPE exponents carry stderr ~0.005; the M136 ZrO2
# spectra, where one CPE leaves 4.5-5 % residuals, span ~0.15.
LOCAL_EXPONENT_SPAN_WARN = 0.1


@dataclass
class LocalExponentResult:
    """Local exponent n(f) of the real part of the admittance.

    Arrays are sorted by ascending frequency. `n` is NaN where Re Y is not
    finite and > 0 or the window holds too few distinct frequencies; `valid`
    marks the points determined to
    LOCAL_EXPONENT_UNCERTAINTY_MAX, and the summary values use only those
    (NaN if there are none).
    """
    frequencies: NDArray[np.float64]  # [Hz]
    n: NDArray[np.float64]
    n_uncertainty: NDArray[np.float64]
    valid: NDArray[np.bool_]
    R_inf: float  # subtracted series resistance [Ohm]
    R_inf_range: Tuple[float, float]  # where the true R_inf may lie; n_uncertainty covers its ends [Ohm]
    L: float  # subtracted series inductance [H]
    window_decades: float
    uncertainty_max: float  # threshold behind `valid`
    n_min: float
    f_n_min: float  # [Hz]
    n_max: float
    f_n_max: float  # [Hz]
    warnings: List[str] = field(default_factory=list)

    @property
    def span(self) -> float:
        """n_max - n_min over the determined points."""
        return self.n_max - self.n_min


def _sliding_slope(log_f: NDArray[np.float64], x: NDArray[np.float64],
                   y: NDArray[np.float64], window: float):
    """Slope of y against x and its stderr, over log_f +- window/2."""
    slope = np.full(len(x), np.nan)
    stderr = np.full(len(x), np.nan)
    for i in range(len(x)):
        if not np.isfinite(y[i]):  # Re Y <= 0 here: no exponent to read
            continue
        m = (np.abs(log_f - log_f[i]) <= window / 2) & np.isfinite(y)
        xm, ym = x[m], y[m]
        # Frequencies, not points: a repeated sweep adds no spread in x. One
        # repeated frequency gave n = 0 +- 0, valid (its mean rounded off by
        # ~1e-16), two with two sweeps a line with no lack-of-fit term
        if len(np.unique(xm)) < LOCAL_EXPONENT_MIN_POINTS:
            continue
        k = len(xm)
        dx = xm - xm.mean()
        b = np.sum(dx * (ym - ym.mean())) / np.sum(dx**2)
        resid = ym - ym.mean() - b * dx
        slope[i] = b
        stderr[i] = np.sqrt(np.sum(resid**2) / (k - 2) / np.sum(dx**2))
    return slope, stderr


def local_exponent(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    R_inf: float,
    L: float = 0.0,
    R_inf_range: Optional[Tuple[float, float]] = None,
) -> LocalExponentResult:
    """
    Map the local CPE exponent n(f) = d ln Re Y / d ln w.

    Y = 1/(Z - R_inf - jwL). The slope is a sliding linear regression over
    LOCAL_EXPONENT_WINDOW_DECADES. Each point's uncertainty combines the
    regression stderr with the change of n when R_inf moves to either end of
    `R_inf_range`, which is what limits it near f_max, where Re Z
    approaches R_inf.

    Parameters
    ----------
    frequencies : ndarray
        Frequencies [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    R_inf : float
        Series resistance to subtract [Ohm], e.g. `estimate_rinf(...).R_inf`
    L : float, optional
        Series inductance to subtract [H] (default 0)
    R_inf_range : (float, float), optional
        Interval the true R_inf lies in [Ohm]; must contain R_inf. Default:
        R_inf +-LOCAL_EXPONENT_RINF_REL. Pass `estimate_rinf(...).R_inf_range`:
        when R_inf is only an upper bound there, it is (0, bound + noise).
        The bound can be 100x R_s or more on an oxide, and +-5 % of it marked
        points as determined to 0.02 that were 0.2 off (stress test,
        oxide/119).

    Returns
    -------
    LocalExponentResult

    Raises
    ------
    ValueError
        On invalid data, or fewer points than one regression window needs.

    Notes
    -----
    Read against the expected shape: one CPE in parallel with C gives a flat
    n(f); a DC conductance G pulls n towards 0 at low frequency; a blocked
    transport channel pulls it up towards its own exponent there. Spectra
    without a blocking layer (several closed arcs) have Re Y nearly flat,
    n(f) near 0, and the map says little about them.
    """
    frequencies = np.asarray(frequencies, dtype=float)
    Z = np.asarray(Z, dtype=complex)
    validate_eis_data(frequencies, Z, context="local_exponent")
    if len(frequencies) < LOCAL_EXPONENT_MIN_POINTS:
        raise ValueError(f"local_exponent: need at least {LOCAL_EXPONENT_MIN_POINTS} "
                         f"points, got {len(frequencies)}")
    if not np.isfinite(R_inf) or not np.isfinite(L):
        raise ValueError("local_exponent: R_inf and L must be finite")
    if R_inf_range is None:
        shift = LOCAL_EXPONENT_RINF_REL * abs(R_inf)
        R_inf_range = (R_inf - shift, R_inf + shift)
    if not (np.all(np.isfinite(R_inf_range)) and R_inf_range[0] <= R_inf <= R_inf_range[1]):
        raise ValueError(f"local_exponent: R_inf_range {R_inf_range} must be finite "
                         f"and contain R_inf = {R_inf}")

    order = np.argsort(frequencies, kind='stable')
    f, Z = frequencies[order], Z[order]
    omega = 2 * np.pi * f
    log_f, x = np.log10(f), np.log(omega)

    def slopes(r_inf: float):
        # Z - R_inf = 0 at a point (estimate_rinf's fallback with Im = 0
        # there) has no exponent, and ends in NaN as it should
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            re_Y = (1 / (Z - r_inf - 1j * omega * L)).real
            y = np.where(re_Y > 0, np.log(np.abs(re_Y)), np.nan)
        return _sliding_slope(log_f, x, y, LOCAL_EXPONENT_WINDOW_DECADES)

    n, stderr = slopes(R_inf)
    n_down, _ = slopes(R_inf_range[0])
    n_up, _ = slopes(R_inf_range[1])
    # np.maximum, not fmax: n undefined at either end (Re Y <= 0 near f_max)
    # must leave the point undetermined, not judged by the other side
    sensitivity = np.maximum(np.abs(n_up - n), np.abs(n_down - n))
    uncertainty = np.hypot(stderr, sensitivity)
    valid = np.isfinite(n) & (uncertainty <= LOCAL_EXPONENT_UNCERTAINTY_MAX)

    warnings: List[str] = []
    if valid.any():
        i_min = int(np.nanargmin(np.where(valid, n, np.nan)))
        i_max = int(np.nanargmax(np.where(valid, n, np.nan)))
        n_min, f_n_min = float(n[i_min]), float(f[i_min])
        n_max, f_n_max = float(n[i_max]), float(f[i_max])
        if n_max - n_min > LOCAL_EXPONENT_SPAN_WARN:
            warnings.append(
                f"n varies from {n_min:.2f} at {f_n_min:.3g} Hz to {n_max:.2f} at "
                f"{f_n_max:.3g} Hz (span {n_max - n_min:.2f}); one CPE has a constant "
                "n. A fall towards 0 at low frequency is a DC conductance, a rise "
                "there a blocked channel; otherwise the dispersion holds more than "
                "one process")
    else:
        n_min = f_n_min = n_max = f_n_max = float('nan')
        warnings.append(
            f"n is determined to {LOCAL_EXPONENT_UNCERTAINTY_MAX} at no frequency: "
            "Re Y is too small or too noisy, or the spectrum is dominated by R_inf")

    return LocalExponentResult(
        frequencies=f, n=n, n_uncertainty=uncertainty, valid=valid,
        R_inf=float(R_inf), R_inf_range=(float(R_inf_range[0]), float(R_inf_range[1])),
        L=float(L),
        window_decades=LOCAL_EXPONENT_WINDOW_DECADES,
        uncertainty_max=LOCAL_EXPONENT_UNCERTAINTY_MAX,
        n_min=n_min, f_n_min=f_n_min, n_max=n_max, f_n_max=f_n_max,
        warnings=warnings,
    )


__all__ = ['local_exponent', 'LocalExponentResult']
