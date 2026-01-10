"""
High-frequency data selection for R_inf estimation.

Provides selection of highest frequency decade for R-L-K() fitting.
Clean design: No logging, all info returned as data.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List
from numpy.typing import NDArray


@dataclass
class DataSelectionResult:
    """Result of high-frequency data selection."""
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    method: str
    n_points: int
    freq_range: tuple  # (f_min, f_max)
    decade_range: tuple  # (f_min_decade, f_max)
    n_inductive: int
    n_capacitive: int
    warnings: List[str] = field(default_factory=list)


def select_highest_decade(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128]
) -> DataSelectionResult:
    """
    Select all points from the highest frequency decade for R-L-K() fit.

    The decade is defined as f_max/10 to f_max, where f_max is the
    highest frequency in the dataset.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies sorted ascending [Hz]
    Z : ndarray of complex
        Impedance corresponding to frequencies [Ohm]

    Returns
    -------
    DataSelectionResult
        Selected data with diagnostics

    Examples
    --------
    >>> result = select_highest_decade(freq, Z)
    >>> print(result.method)
    'highest_decade_11points'
    """
    warnings = []

    if len(frequencies) == 0:
        raise ValueError("Empty frequency array")

    f_max = frequencies[-1]
    f_min_decade = f_max / 10.0

    # Find points in highest decade
    mask_decade = frequencies >= f_min_decade
    n_in_decade = int(np.sum(mask_decade))

    if n_in_decade == 0:
        # No points in decade - use all data (fallback)
        warnings.append(
            f"No points in highest decade ({f_min_decade/1e6:.3f} - {f_max/1e6:.3f} MHz), "
            f"using all data"
        )
        f_high = frequencies
        Z_high = Z
        method = f"fallback_all_{len(frequencies)}points"
    else:
        f_high = frequencies[mask_decade]
        Z_high = Z[mask_decade]
        method = f"highest_decade_{n_in_decade}points"

    n_inductive = int(np.sum(Z_high.imag > 0))
    n_capacitive = int(np.sum(Z_high.imag < 0))

    return DataSelectionResult(
        frequencies=f_high,
        Z=Z_high,
        method=method,
        n_points=len(f_high),
        freq_range=(float(f_high.min()), float(f_high.max())),
        decade_range=(float(f_min_decade), float(f_max)),
        n_inductive=n_inductive,
        n_capacitive=n_capacitive,
        warnings=warnings
    )


__all__ = ['select_highest_decade', 'DataSelectionResult']
