"""
Validation of the tau grid for Voigt chain fitting (the EIS data:
utils.validation).
"""

import numpy as np
from numpy.typing import NDArray


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
        n_invalid = int(np.sum(tau <= 0))
        raise ValueError(
            f"{context}: tau values must be positive "
            f"(found {n_invalid} non-positive values)"
        )

    if not np.all(np.isfinite(tau)):
        n_invalid = int(np.sum(~np.isfinite(tau)))
        raise ValueError(
            f"{context}: tau contains {n_invalid} NaN/Inf values"
        )


__all__ = ['validate_tau']
