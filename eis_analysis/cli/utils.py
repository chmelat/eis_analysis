"""
Utility functions and dataclasses for EIS CLI.

Contains:
- Exception classes
- Data containers (dataclasses)
- Helper functions (save_figure, parse_circuit_expression)
"""

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..fitting import R, C, Q, L, W, Wo, K

logger = logging.getLogger(__name__)


# =============================================================================
# Exceptions
# =============================================================================

class EISAnalysisError(Exception):
    """Base exception for EIS analysis errors."""
    pass


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class LoadedData:
    """
    Container for loaded EIS data.

    Attributes
    ----------
    frequencies : ndarray
        Frequency values [Hz]
    Z : ndarray
        Complex impedance values [Ohm]
    title : str
        Data title (filename or "Synthetic data")
    metadata : dict or None
        Metadata from DTA file
    ocv_data : dict or None
        OCV curve data from DTA file
    """
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    title: str
    metadata: Optional[dict]
    ocv_data: Optional[dict] = None


@dataclass
class CLIDRTResult:
    """
    Container for DRT analysis results in CLI context.

    This is a simplified wrapper for CLI use. For full DRT results,
    use eis_analysis.drt.DRTResult directly.

    Attributes
    ----------
    tau : ndarray or None
        Time constants [s]
    gamma : ndarray or None
        DRT spectrum [Ohm]
    peaks_gmm : list or None
        GMM-detected peaks
    fig_drt : Figure or None
        DRT plot figure
    fig_rinf : Figure or None
        R_inf estimation figure
    """
    tau: Optional[NDArray]
    gamma: Optional[NDArray]
    peaks_gmm: Optional[list]
    fig_drt: Optional[plt.Figure]
    fig_rinf: Optional[plt.Figure]


# =============================================================================
# Helper Functions
# =============================================================================

def save_figure(
    fig: Optional[plt.Figure],
    prefix: Optional[str],
    suffix: str,
    fmt: str = 'png'
) -> None:
    """
    Save figure to file if fig and prefix are provided.

    Parameters
    ----------
    fig : Figure or None
        Matplotlib figure to save
    prefix : str or None
        File prefix (from --save argument)
    suffix : str
        File suffix (e.g., 'nyquist_bode', 'kk', 'drt')
    fmt : str
        Output format: 'png', 'pdf', 'svg', 'eps' (default: 'png')
    """
    if fig is None or prefix is None:
        return

    filepath = f"{prefix}_{suffix}.{fmt}"
    try:
        # For bitmap formats use dpi, for vector formats use bbox_inches
        if fmt == 'png':
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
        else:
            # Vector formats (pdf, svg, eps)
            fig.savefig(filepath, bbox_inches='tight')
        logger.info(f"Saved: {filepath}")
    except Exception as e:
        logger.error(f"Error saving figure: {e}")


def parse_circuit_expression(expr: str):
    """
    Parse circuit expression string into circuit object.

    Parameters
    ----------
    expr : str
        Circuit expression, e.g. "R(100) - (R(5000) | C(1e-6))"

    Returns
    -------
    circuit : Circuit object
        Parsed circuit

    Raises
    ------
    ValueError
        If expression cannot be parsed

    Examples
    --------
    >>> circuit = parse_circuit_expression("R(100) - (R(5000) | C(1e-6))")
    >>> print(circuit.get_all_params())
    [100, 5000, 1e-6]

    Notes
    -----
    Uses eval() with restricted namespace for safety. Only circuit element
    classes (R, C, Q, L, W, Wo, K) are available in the evaluation context.
    """
    # Safe namespace for eval - only circuit elements
    safe_namespace = {
        'R': R,
        'C': C,
        'Q': Q,
        'L': L,
        'W': W,
        'Wo': Wo,
        'K': K,
    }

    try:
        circuit = eval(expr, {"__builtins__": {}}, safe_namespace)
        return circuit
    except Exception as e:
        raise ValueError(f"Cannot parse circuit expression '{expr}': {e}")
