"""
Linearity check from the total harmonic distortion an instrument records.

Gamry writes THD of the current and voltage per point when its THD option is
on (`LoadResult.current_thd` / `voltage_thd`). It is the only direct
instrument diagnostic of linearity; KK and Z-HIT see a nonlinear response only
indirectly, and not always. This module only summarizes it against a
threshold - it does not try to tell nonlinearity from noise or interference
(drift alone puts energy into the 2nd harmonic, so the harmonic pattern does
not decide it reliably).
"""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

# For an odd nonlinearity i = g1*v + g3*v^3 driven by v = A*cos(wt), the
# cubic term adds 3/4*g3*A^3 to the fundamental and 1/4*g3*A^3 to the 3rd
# harmonic: the error of the fundamental is 3x the 3rd-harmonic ratio. A THD
# of 1 % therefore bounds the |Z| error from nonlinearity to about 3 %, below
# the 5 % residual threshold of the KK validation. An even nonlinearity does
# not shift the fundamental to first order, so the bound is conservative.
THD_THRESHOLD = 0.01

# |Z| error per unit of THD in the worst (purely cubic) case, see above.
THD_TO_Z_ERROR = 3.0


@dataclass
class THDChannel:
    """
    THD summary of one channel (current or voltage). All values are fractions.

    Attributes
    ----------
    median : float
        Median THD over the points with a value
    maximum : float
        Largest THD
    f_at_max : float
        Frequency of the largest THD [Hz]
    n_valid : int
        Number of points with a value (finite and > 0) - what `n_above`
        counts among
    n_above : int
        Number of points above the threshold
    f_above_min, f_above_max : float or None
        Frequency span of the points above the threshold [Hz]; None if none is
    """
    median: float
    maximum: float
    f_at_max: float
    n_valid: int
    n_above: int
    f_above_min: Optional[float] = None
    f_above_max: Optional[float] = None


@dataclass
class THDResult:
    """
    Result of `thd_check`.

    Attributes
    ----------
    current, voltage : THDChannel or None
        Summary per channel; None if the data has no usable value for it.
        In potentiostatic EIS the current THD is the response of the sample,
        the voltage THD the purity of the excitation (galvanostatic: reversed).
    threshold : float
        THD threshold used (fraction)
    n_points : int
        Number of points of the spectrum
    warnings : list of str
        Caveats about the summary: a channel present but without any value,
        or missing at some points
    """
    current: Optional[THDChannel]
    voltage: Optional[THDChannel]
    threshold: float
    n_points: int
    warnings: List[str] = field(default_factory=list)


def _summarize(frequencies: NDArray[np.float64], thd: NDArray[np.float64],
               threshold: float) -> Optional[THDChannel]:
    # No measurement gives exactly zero distortion - the noise floor alone is
    # above it - so 0 (or a negative placeholder) means "not measured". Counting
    # it would turn an empty column into a clean bill of linearity.
    valid = np.isfinite(thd) & (thd > 0)
    if not valid.any():
        return None
    f, values = frequencies[valid], thd[valid]
    i_max = int(np.argmax(values))
    above = values > threshold
    return THDChannel(
        median=float(np.median(values)),
        maximum=float(values[i_max]),
        f_at_max=float(f[i_max]),
        n_valid=int(valid.sum()),
        n_above=int(above.sum()),
        f_above_min=float(f[above].min()) if above.any() else None,
        f_above_max=float(f[above].max()) if above.any() else None,
    )


def thd_check(
    frequencies: NDArray[np.float64],
    current_thd: Optional[NDArray[np.float64]],
    voltage_thd: Optional[NDArray[np.float64]],
    threshold: float = THD_THRESHOLD,
) -> Optional[THDResult]:
    """
    Summarize per-point THD against a linearity threshold.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz]
    current_thd, voltage_thd : ndarray of float or None
        THD per point as a fraction, aligned with `frequencies`; NaN or a
        value <= 0 marks a missing one. None when the data carries no such
        column.
    threshold : float
        THD above which a point is counted (fraction, default 1 %)

    Returns
    -------
    THDResult or None
        None when neither channel is given (nothing to report).
    """
    if current_thd is None and voltage_thd is None:
        return None

    frequencies = np.asarray(frequencies, dtype=np.float64)
    warnings: List[str] = []
    channels: List[Optional[THDChannel]] = []
    for label, thd in (('current', current_thd), ('voltage', voltage_thd)):
        if thd is None:
            channels.append(None)
            continue
        thd = np.asarray(thd, dtype=np.float64)
        if len(thd) != len(frequencies):
            raise ValueError(f"{label} THD has {len(thd)} values for {len(frequencies)} frequencies")
        summary = _summarize(frequencies, thd, threshold)
        if summary is None:
            warnings.append(f"{label.capitalize()} THD column present but holds no value")
        elif summary.n_valid < len(frequencies):
            warnings.append(f"{label.capitalize()} THD missing at "
                            f"{len(frequencies) - summary.n_valid}/{len(frequencies)} points")
        channels.append(summary)

    return THDResult(current=channels[0], voltage=channels[1], threshold=threshold,
                     n_points=len(frequencies), warnings=warnings)
