"""
The spectrum a loader returns, and the checks every loader applies to it.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
from numpy.typing import NDArray

# Validation constants
MIN_DATA_POINTS = 10  # Minimum number of data points for analysis
MIN_FREQUENCY_RANGE = 10  # Minimum ratio f_max/f_min


@dataclass
class LoadResult:
    """
    A spectrum as it came out of a file.

    The only points removed are the high-frequency run with Re(Z) < 0, a lead
    artifact; `warnings` says how many (see _drop_negative_real_hf).

    Caveats about the data land in `warnings` rather than on the console:
    the loader has no idea whether it runs under the CLI, in a notebook or
    in a batch script, and the caveat qualifies the returned spectrum the
    way an uncertainty qualifies a measurement. Failures of the operation
    itself (unreadable file, missing section) still raise or log, since
    there is no result for them to qualify.

    Attributes
    ----------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ohm]
    filename : str
        Path the data was read from
    metadata : dict or None
        DTA header metadata; None for formats that carry none (CSV)
    warnings : list of str
        Caveats about the data, in the order they were found
    current_thd, voltage_thd : ndarray of float or None
        Total harmonic distortion of the current and voltage per point, from
        the Gamry ``Ithd``/``Vthd`` columns (written when the THD option is
        on), aligned with `frequencies`; NaN where a cell is empty. None when
        the file has no such column. A fraction, not percent: it equals
        sqrt(sum_{n=2..10} |H_n|^2) / |H_1| of the harmonic columns exactly
        (checked on example/EISPOT-test1.DTA), although the header gives
        only '#' as the unit.
    """
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    filename: str
    metadata: Optional[Dict[str, Any]] = None
    warnings: List[str] = field(default_factory=list)
    current_thd: Optional[NDArray[np.float64]] = None
    voltage_thd: Optional[NDArray[np.float64]] = None

    def keep_points(self, mask: NDArray[np.bool_]) -> None:
        """Keep only the points where `mask` is True, in every per-point field.

        The one place that knows which fields are per point, so a column
        added later cannot be left misaligned with `frequencies`.
        """
        self.frequencies, self.Z = self.frequencies[mask], self.Z[mask]
        if self.current_thd is not None:
            self.current_thd = self.current_thd[mask]
        if self.voltage_thd is not None:
            self.voltage_thd = self.voltage_thd[mask]


def _drop_negative_real_hf(frequencies: NDArray[np.float64], Z: NDArray[np.complex128],
                           warnings: List[str]) -> NDArray[np.bool_]:
    """
    Mask of the points to keep: all but the high-frequency run with Re(Z) < 0,
    which is noted in `warnings`. A mask rather than the trimmed arrays, so
    `LoadResult.keep_points` trims every per-point column the same way.

    No passive system has Re(Z) < 0. At the top of a sweep it is a lead
    artifact (cable inductance resonating with stray capacitance) that breaks
    Lin-KK, the DRT and the circuit fit, and turns the HF R_inf negative.
    Only the contiguous run from the highest frequency down is removed: a
    negative differential resistance (passivation, oscillating systems under
    DC bias) gives Re(Z) < 0 at low frequencies as real physics, so points
    elsewhere are kept and only noted. A run interrupted by a positive point is
    not bridged: Re(Z) flipping sign there is within noise of zero, not a
    clear artifact.

    Raises
    ------
    ValueError
        If every point has Re(Z) < 0
    """
    # Descending; stable, so a duplicated top frequency resolves in file order
    order = np.argsort(-frequencies, kind='stable')
    negative = Z.real[order] < 0
    if negative.all():
        raise ValueError("All points have Re(Z) < 0 - check the sign convention of the data")
    n_hf = int(np.argmin(negative))

    mask = np.ones(len(frequencies), dtype=bool)
    mask[order[:n_hf]] = False

    if n_hf:
        dropped = frequencies[~mask]
        warnings.append(f"Dropped {n_hf} high-frequency point(s) with Re(Z) < 0 "
                        f"({dropped.min():.2e} - {dropped.max():.2e} Hz): not possible "
                        f"for a passive system, typically a lead artifact")

    remaining = mask & (Z.real < 0)
    if remaining.any():
        f_neg = frequencies[remaining]
        warnings.append(f"{int(remaining.sum())} point(s) below the HF end have Re(Z) < 0 "
                        f"({f_neg.min():.2e} - {f_neg.max():.2e} Hz), kept: a negative "
                        f"resistance or a measurement problem")

    return mask


def _check_spectrum(frequencies: NDArray[np.float64], warnings: List[str]) -> None:
    """
    Checks every loader applies to the spectrum it read, noting caveats in `warnings`.

    Raises
    ------
    ValueError
        If there are fewer than MIN_DATA_POINTS points
    """
    if len(frequencies) < MIN_DATA_POINTS:
        raise ValueError(f"Dataset must have at least {MIN_DATA_POINTS} points, got {len(frequencies)}")

    freq_range = frequencies.max() / frequencies.min()
    if freq_range < MIN_FREQUENCY_RANGE:
        warnings.append(
            f"Small frequency range: {freq_range:.1f}x "
            f"(recommended >{MIN_FREQUENCY_RANGE}x); DRT analysis may have poor resolution")

    # A sweep is strictly monotonic in file order. Several sweeps in one file
    # break that even when the instrument logged slightly different measured
    # frequencies, which exact-equality (np.unique) would miss. Z-HIT takes
    # its phase derivative over a minimum step (MIN_DERIVATIVE_STEP), so such
    # points no longer break it; repeated sweeps that disagree show up as its
    # residuals.
    steps = np.diff(frequencies)
    n_against = int(min(np.sum(steps >= 0), np.sum(steps <= 0)))
    if n_against:
        warnings.append(
            f"Dataset contains duplicate or out-of-order frequencies ({n_against} step(s) "
            f"against the sweep direction) - several sweeps in one file?")

