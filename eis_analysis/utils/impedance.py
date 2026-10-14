"""
Impedance data utility functions.

This module contains common operations on impedance data, following the
DRY (Don't Repeat Yourself) principle from CLAUDE.md.

Functions
---------
calculate_rpol : Calculate polarization resistance R_pol = R_dc - R_inf
sort_by_frequency : Sort impedance data by frequency in ascending order
"""

import numpy as np
from numpy.typing import NDArray
from typing import Optional, Tuple


# HF median: up to 5 highest points, but no more than 10 % of the spectrum,
# so short spectra do not reach down into the arc.
HF_MEDIAN_MAX_POINTS = 5

# Re(Z) at the top capacitive point (hf_bound) is an upper bound of R_inf for
# a passive system. On a flat end the HF median lies within the noise of it;
# one more than 5 % above both of the two top capacitive points reaches into
# an arc still open there (cold-agent T10: median 17.5 Ohm, top points 11.6
# and 14.1, R_s 5). Two points, not one: a single low glitch at f_max must
# not replace the median it was taken to resist. 5 % is the
# RINF_REL_STDERR_MAX margin: at 1 % noise the median of 5 points and one
# point differ by ~1.2 % (1 sigma).
HF_MEDIAN_ABOVE_BOUND_MAX = 0.05


def calculate_rpol(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    n_avg: int = 5
) -> Tuple[float, float, float]:
    """
    Calculate polarization resistance R_pol = R_dc - R_inf.

    R_dc is the median of the n_avg lowest-frequency points, R_inf that of
    the highest ones (hf_r_inf: capped at the top capacitive point where an
    arc is still open there).

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (N points)
    Z : ndarray of complex
        Complex impedance [Ω] (N points)
    n_avg : int, optional
        Number of low-frequency points for R_dc (default: 5), limited to
        max(1, len(frequencies) // 10); the top uses HF_MEDIAN_MAX_POINTS

    Returns
    -------
    R_pol : float
        Polarization resistance R_pol = R_dc - R_inf [Ω]
    R_inf : float
        High-frequency resistance (series resistance) [Ω]
    R_dc : float
        Low-frequency resistance (total resistance) [Ω]

    Notes
    -----
    Uses median instead of mean for robustness against:
    - Noise in high-frequency region
    - Deviations in low-frequency region (non-steady state)
    - Measurement outliers

    Physical meaning:
    - R_inf: Ohmic resistance of electrolyte + contacts + wiring
    - R_dc: Total resistance including electrode polarization
    - R_pol: Resistance of polarization processes (charge transfer, diffusion, etc.)

    Examples
    --------
    >>> freq = np.logspace(0, 5, 100)  # 1 Hz - 100 kHz
    >>> Z = 10 + 100/(1 + 1j*2*np.pi*freq*0.001)  # R0 + RC circuit
    >>> R_pol, R_inf, R_dc = calculate_rpol(freq, Z)
    >>> print(f"R_inf = {R_inf:.1f} Ω")
    R_inf = 10.0 Ω
    >>> print(f"R_pol = {R_pol:.1f} Ω")
    R_pol = 100.0 Ω

    See Also
    --------
    sort_by_frequency : Sort data by frequency (recommended before calling)
    """
    # Automatically limit n_avg to be reasonable for small datasets
    n_avg = min(n_avg, max(1, len(frequencies) // 10))

    # Indices of the n_avg lowest frequencies
    low_freq_indices = np.argsort(frequencies)[:n_avg:]

    # Median of Re(Z) at each end; at the top capped where an arc is still
    # open (hf_r_inf), as the DRT does, so both see the same R_inf
    R_inf, _ = hf_r_inf(frequencies, Z)
    R_dc = np.median(Z.real[low_freq_indices])

    # Polarization resistance
    R_pol = R_dc - R_inf

    return R_pol, R_inf, R_dc


def sort_by_frequency(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating]
) -> Tuple[NDArray[np.float64], NDArray[np.complexfloating]]:
    """
    Sort impedance data by frequency in ascending order.

    Instruments save spectra in either direction; callers that assume
    ascending order (edge-point selection in R_inf estimation, R_pol
    from sorted data) normalize through this function first.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (N points), may be unsorted
    Z : ndarray of complex
        Complex impedance [Ω] (N points)

    Returns
    -------
    frequencies_sorted, Z_sorted : ndarray
        Both arrays sorted by ascending frequency
    """
    sort_idx = np.argsort(frequencies)
    return frequencies[sort_idx], Z[sort_idx]


def hf_median(frequencies: NDArray, Z: NDArray) -> Tuple[float, int]:
    """Median of Re(Z) over the highest frequencies.

    Returns ``(R_inf, n_points)``; see HF_MEDIAN_MAX_POINTS.
    """
    n = min(HF_MEDIAN_MAX_POINTS, max(1, len(frequencies) // 10))
    idx = np.argsort(frequencies)[-n:]
    return float(np.median(Z.real[idx])), n


def hf_bound(frequencies: NDArray, Z: NDArray) -> Tuple[float, float, int]:
    """Re(Z) at the highest frequency with Im(Z) <= 0, f_max if there is none.

    Returns ``(R, f, index)``. Every point of the R-L-(R|Q) model is an upper bound
    of R_s; this one is the first below an inductive top, where lead
    artifacts can pull Re(Z) below R_s (redoxED flow cell: 0.168 Ohm at the
    Im = 0 crossing, 0.003 and -0.064 Ohm above it). On a capacitive top it
    is f_max, the tightest bound.
    """
    order = np.argsort(frequencies)[::-1]
    capacitive = np.flatnonzero(Z.imag[order] <= 0)
    i = order[capacitive[0]] if capacitive.size else order[0]
    return float(Z.real[i]), float(frequencies[i]), int(i)


def hf_r_inf(frequencies: NDArray, Z: NDArray) -> Tuple[float, Optional[str]]:
    """R_inf from the top of the spectrum: the HF median, unless it overshoots.

    A median more than HF_MEDIAN_ABOVE_BOUND_MAX (of the median) above the two
    highest points with Im(Z) <= 0 reaches into an arc still open there; then
    hf_bound, floored at 0, is returned - still an upper bound of R_inf.

    Returns ``(R_inf, warning)``, warning None when the median is used.
    """
    median, _ = hf_median(frequencies, Z)
    order = np.argsort(frequencies)[::-1]
    top = order[Z.imag[order] <= 0][:2]
    margin = HF_MEDIAN_ABOVE_BOUND_MAX * abs(median)
    if len(top) < 2 or not all(median - Z.real[i] > margin for i in top):
        return median, None
    bound, f_bound, _ = hf_bound(frequencies, Z)
    where = ("f_max" if f_bound == float(np.max(frequencies))
             else f"{f_bound:.3g} Hz, the highest point with Im(Z) <= 0")
    return max(bound, 0.0), (
        f"HF median R_inf = {median:.4g} Ohm lies above Re(Z) = {bound:.4g} Ohm at "
        f"{where}, an upper bound of R_inf: an arc is still open at the top of the "
        f"spectrum. Using the bound, still an upper bound; a fit with the series "
        f"elements (or --ri-fit) can determine R_s")
