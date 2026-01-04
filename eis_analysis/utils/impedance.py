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
from typing import Tuple


def calculate_rpol(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_avg: int = 5
) -> Tuple[float, float, float]:
    """
    Calculate polarization resistance R_pol = R_dc - R_inf.

    Uses median of n_avg edge points for robustness against outliers.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (N points)
    Z : ndarray of complex
        Complex impedance [Ω] (N points)
    n_avg : int, optional
        Number of points for averaging at each end of spectrum (default: 5)
        Automatically limited to max(1, len(frequencies) // 10)

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

    # Find indices of n_avg highest and lowest frequencies
    high_freq_indices = np.argsort(frequencies)[-n_avg:]
    low_freq_indices = np.argsort(frequencies)[:n_avg:]

    # Use median of real part of impedance
    R_inf = np.median(Z.real[high_freq_indices])
    R_dc = np.median(Z.real[low_freq_indices])

    # Polarization resistance
    R_pol = R_dc - R_inf

    return R_pol, R_inf, R_dc


def sort_by_frequency(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128]
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Sort impedance data by frequency in ascending order.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (N points), may be unsorted
    Z : ndarray of complex
        Complex impedance [Ω] (N points)

    Returns
    -------
    frequencies_sorted : ndarray of float
        Frequencies sorted in ascending order [Hz]
    Z_sorted : ndarray of complex
        Impedance sorted by frequency [Ω]

    Notes
    -----
    Some EIS instruments save data in descending order (from high frequencies),
    others in ascending order. This function ensures consistent ordering (ascending)
    for further processing.

    Ascending order is preferred because:
    - Corresponds to time axis (low f = slow processes)
    - Consistent with DRT calculations (integrals over log(f))
    - Easier selection of edge points (first = low freq, last = high freq)

    Examples
    --------
    >>> freq = np.array([1000, 100, 10, 1])  # Descending
    >>> Z = np.array([10+1j, 10+5j, 10+20j, 10+50j])
    >>> freq_sorted, Z_sorted = sort_by_frequency(freq, Z)
    >>> print(freq_sorted)
    [1 10 100 1000]

    See Also
    --------
    calculate_rpol : Calculate R_pol (prefers sorted data for speed)
    """
    sort_idx = np.argsort(frequencies)
    return frequencies[sort_idx], Z[sort_idx]
