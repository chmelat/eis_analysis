"""
High-frequency data selection for R_inf estimation.

Provides selection of highest frequency decade for R-L-K() fitting.
"""

import numpy as np
import logging
from typing import Tuple
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def select_highest_decade(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    verbose: bool = True
) -> Tuple[NDArray[np.float64], NDArray[np.complex128], str, list]:
    """
    Select all points from the highest frequency decade for R-L-K() fit.

    The decade is defined as f_max/10 to f_max, where f_max is the
    highest frequency in the dataset.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies sorted ascending [Hz]
    Z : ndarray of complex
        Impedance corresponding to frequencies [Ω]
    verbose : bool
        Log selection info (default: True)

    Returns
    -------
    f_high : ndarray of float
        Selected frequencies [Hz]
    Z_high : ndarray of complex
        Selected impedance [Ω]
    method_str : str
        Method description (for diagnostics)
    warnings : list of str
        List of warnings (if any issues)

    Examples
    --------
    >>> f_high, Z_high, method, warnings = select_highest_decade(freq, Z)
    >>> print(method)
    'highest_decade_11points'
    """
    warnings = []

    if len(frequencies) == 0:
        raise ValueError("Empty frequency array")

    f_max = frequencies[-1]
    f_min_decade = f_max / 10.0  # Exactly 1 decade

    # Find points in highest decade
    mask_decade = frequencies >= f_min_decade
    n_in_decade = np.sum(mask_decade)

    if n_in_decade == 0:
        # No points in decade - use all data (fallback)
        msg = (
            f"No points in highest decade ({f_min_decade/1e6:.3f} - {f_max/1e6:.3f} MHz). "
            f"Using all data."
        )
        warnings.append(msg)
        logger.warning(msg)
        f_high = frequencies
        Z_high = Z
        method_str = f"fallback_all_{len(frequencies)}points"

    else:
        # Use all points in decade
        f_high = frequencies[mask_decade]
        Z_high = Z[mask_decade]

        method_str = f"highest_decade_{n_in_decade}points"

        if verbose:
            logger.info("Selecting highest frequency decade")
            logger.info(f"  Decade range: {f_min_decade/1e6:.3f} - {f_max/1e6:.3f} MHz")
            logger.info(f"  Points in decade: {n_in_decade}")
            logger.info(f"  Frequency range: {f_high.min()/1e6:.3f} - {f_high.max()/1e6:.3f} MHz")

    if verbose:
        n_inductive = np.sum(Z_high.imag > 0)
        n_capacitive = np.sum(Z_high.imag < 0)
        logger.info(f"  Inductive points (Im(Z)>0): {n_inductive}/{len(f_high)}")
        logger.info(f"  Capacitive points (Im(Z)<0): {n_capacitive}/{len(f_high)}")

    return f_high, Z_high, method_str, warnings


__all__ = ['select_highest_decade']
