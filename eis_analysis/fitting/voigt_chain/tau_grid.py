"""
Tau grid generation for Voigt chain fitting.

This module provides functions to generate logarithmically-spaced
time constants for Voigt elements.
"""

import numpy as np
import logging
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def generate_tau_grid(
    frequencies: NDArray[np.float64],
    n_per_decade: int = 3,
    extend_decades: float = 0.0
) -> NDArray[np.float64]:
    """
    Generate logarithmically-spaced time constants tau for Voigt elements.

    The grid covers the measured frequency range (optionally extended by
    `extend_decades` toward lower frequencies) with `n_per_decade` points
    per decade.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz] (N points)
    n_per_decade : int, optional
        Number of tau values per decade (default: 3)
        Higher values -> better resolution, more parameters
    extend_decades : float, optional
        Extend tau range toward lower frequencies (default: 0.0 for Lin-KK compatibility)
        Set to >0 to capture slow processes beyond measured range

    Returns
    -------
    tau : ndarray of float
        Time constants tau = 1/(2*pi*f) [s], logarithmically spaced

    Notes
    -----
    Time constant tau corresponds to frequency f via:
        tau = 1 / (2*pi*f)
        f = 1 / (2*pi*tau)

    For Lin-KK compatibility, extend_decades should be 0.0 so that tau range
    exactly matches the measured frequency range.

    Examples
    --------
    >>> freq = np.logspace(4, -1, 30)  # 10 kHz to 0.1 Hz
    >>> tau = generate_tau_grid(freq, n_per_decade=3)
    >>> # tau will span exactly the frequency range with 3 points per decade
    """
    # Frequency range
    f_min = frequencies.min()
    f_max = frequencies.max()

    # Extend toward LOWER frequencies (higher tau) if requested
    # tau = 1/(2*pi*f), so lower f -> higher tau
    if extend_decades > 0:
        f_min_extended = f_min / (10 ** extend_decades)
    else:
        f_min_extended = f_min  # Lin-KK compatible: exact frequency range

    # Convert to tau range
    # tau_max corresponds to f_min_extended (slowest process)
    # tau_min corresponds to f_max (fastest process)
    tau_min = 1.0 / (2 * np.pi * f_max)
    tau_max = 1.0 / (2 * np.pi * f_min_extended)

    # Number of decades in tau space
    n_decades = np.log10(tau_max / tau_min)

    # Total number of tau points
    n_tau = int(np.ceil(n_decades * n_per_decade)) + 1

    # Generate logarithmic grid
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    logger.debug(f"Generated {len(tau)} tau values:")
    logger.debug(f"  tau range: [{tau.min():.3e}, {tau.max():.3e}] s")
    logger.debug(f"  f range (equiv): [{1/(2*np.pi*tau.max()):.3e}, {1/(2*np.pi*tau.min()):.3e}] Hz")
    logger.debug(f"  Measured f range: [{f_min:.3e}, {f_max:.3e}] Hz")
    logger.debug(f"  Extension: {extend_decades} decades toward lower f")

    return tau


def generate_tau_grid_fixed_M(
    frequencies: NDArray[np.float64],
    M: int,
    extend_decades: float = 0.0
) -> NDArray[np.float64]:
    """
    Generate logarithmically-spaced time constants tau with exact count M.

    This follows the Lin-KK approach from Schonleber et al. (2014):
        tau_1 = 1/omega_max ; tau_M = 1/omega_min
        tau_k = 10^(log(tau_min) + (k-1)/(M-1) * log(tau_max/tau_min))

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz] (N points)
    M : int
        Exact number of time constants to generate
    extend_decades : float, optional
        Extend tau range toward lower frequencies (default: 0.0 for Lin-KK compatibility)

    Returns
    -------
    tau : ndarray of float
        Time constants tau [s], exactly M points, logarithmically spaced

    Notes
    -----
    For Lin-KK compatibility, extend_decades should be 0.0 so that tau range
    exactly matches the measured frequency range.

    Examples
    --------
    >>> freq = np.logspace(4, -1, 30)  # 10 kHz to 0.1 Hz
    >>> tau = generate_tau_grid_fixed_M(freq, M=5, extend_decades=0)
    >>> len(tau)
    5
    """
    # Frequency range
    f_min = frequencies.min()
    f_max = frequencies.max()

    # Extend toward LOWER frequencies (higher tau) if requested
    if extend_decades > 0:
        f_min_extended = f_min / (10 ** extend_decades)
    else:
        f_min_extended = f_min

    # Convert to tau range (Lin-KK standard: tau = 1/omega, not 1/(2*pi*f))
    # But for consistency with Voigt elements, we use tau = 1/(2*pi*f)
    tau_min = 1.0 / (2 * np.pi * f_max)
    tau_max = 1.0 / (2 * np.pi * f_min_extended)

    # Generate exactly M points logarithmically
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), M)

    logger.debug(f"Generated {len(tau)} tau values (fixed M={M}):")
    logger.debug(f"  tau range: [{tau.min():.3e}, {tau.max():.3e}] s")
    logger.debug(f"  f range (equiv): [{1/(2*np.pi*tau.max()):.3e}, {1/(2*np.pi*tau.min()):.3e}] Hz")

    return tau


__all__ = ['generate_tau_grid', 'generate_tau_grid_fixed_M']
