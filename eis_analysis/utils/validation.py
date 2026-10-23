"""
Input validation shared by Lin-KK, Z-HIT, the DRT, the Voigt chain fit and
the local exponent map.
"""

import numpy as np
from numpy.typing import NDArray


def validate_eis_data(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
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
        If data is invalid: unequal lengths, no points, f <= 0, NaN or inf
        in f or Z, or |Z|^2 below the smallest normal float (|Z| < ~1.5e-154
        Ohm, 1/|Z|^2 would overflow).

    Notes
    -----
    Known limit: the analyses are scale-invariant only for |Z| within about
    1e-150..1e150 Ohm (measured impedances span 1e-6..1e13). Near the
    bounds sums of 1/|Z|^2 can still overflow, and compute_weights floors
    |Z| at 1e-15, so modulus weights turn uniform below 1e-15 Ohm.
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
        n_invalid = int(np.sum(frequencies <= 0))
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

    # Weights (modulus: 1/|Z|), residuals and pseudo chi^2 are relative to |Z|,
    # and the fits divide by |Z|^2: below the smallest normal float (|Z| <
    # ~1.5e-154 Ohm) it is 0 or subnormal, and 1/|Z|^2 overflows to inf
    tiny = np.finfo(float).tiny
    n_small = int(np.count_nonzero(Z.real**2 + Z.imag**2 < tiny))
    if n_small:
        raise ValueError(f"{context}: |Z| is 0 or below {np.sqrt(tiny):.1e} Ohm at {n_small} "
                         f"point(s), but weights and residuals are relative to |Z|")
