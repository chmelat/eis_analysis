"""
Validation functions for Voigt chain fitting.

This module provides input validation for EIS data and tau grids.
"""

import numpy as np
from numpy.typing import NDArray


def validate_eis_data(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    context: str = "EIS data"
) -> None:
    """
    Validate EIS input data.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance array [Ohm]
    context : str
        Context string for error messages

    Raises
    ------
    ValueError
        If data is invalid
    """
    # Check array lengths match
    if len(frequencies) != len(Z):
        raise ValueError(
            f"{context}: frequencies and Z must have same length "
            f"(got {len(frequencies)} vs {len(Z)})"
        )

    # Check for empty arrays
    if len(frequencies) == 0:
        raise ValueError(f"{context}: empty data arrays")

    # Check frequencies are positive
    if np.any(frequencies <= 0):
        n_invalid = np.sum(frequencies <= 0)
        raise ValueError(
            f"{context}: frequencies must be positive "
            f"(found {n_invalid} non-positive values)"
        )

    # Check for NaN/Inf in frequencies
    if not np.all(np.isfinite(frequencies)):
        n_invalid = np.sum(~np.isfinite(frequencies))
        raise ValueError(
            f"{context}: frequencies contain {n_invalid} NaN/Inf values"
        )

    # Check for NaN/Inf in impedance
    if not np.all(np.isfinite(Z)):
        n_invalid = np.sum(~np.isfinite(Z))
        raise ValueError(
            f"{context}: impedance Z contains {n_invalid} NaN/Inf values"
        )


def validate_tau(
    tau: NDArray[np.float64],
    context: str = "tau grid"
) -> None:
    """
    Validate time constant array.

    Parameters
    ----------
    tau : ndarray
        Time constants [s]
    context : str
        Context string for error messages

    Raises
    ------
    ValueError
        If tau values are invalid
    """
    if len(tau) == 0:
        raise ValueError(f"{context}: empty tau array")

    if np.any(tau <= 0):
        n_invalid = np.sum(tau <= 0)
        raise ValueError(
            f"{context}: tau values must be positive "
            f"(found {n_invalid} non-positive values)"
        )

    if not np.all(np.isfinite(tau)):
        n_invalid = np.sum(~np.isfinite(tau))
        raise ValueError(
            f"{context}: tau contains {n_invalid} NaN/Inf values"
        )


__all__ = ['validate_eis_data', 'validate_tau']
