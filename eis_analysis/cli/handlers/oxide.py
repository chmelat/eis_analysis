"""
Oxide layer analysis handler for the EIS CLI.

- run_oxide_analysis: oxide thickness estimation from capacitance
"""

import argparse
import logging
from typing import Optional

from numpy.typing import NDArray

from ...analysis import analyze_oxide_layer
from ...fitting import FitResult

logger = logging.getLogger(__name__)


def run_oxide_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    fitted_result: Optional[FitResult],
    metadata: Optional[dict]
) -> None:
    """
    Run oxide layer analysis.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: analyze_oxide, epsilon_r, area)
    fitted_result : FitResult or None
        Circuit fitting result
    metadata : dict or None
        DTA file metadata
    """
    if not args.analyze_oxide:
        return

    # Use area from metadata if available and not explicitly specified
    area_to_use = args.area
    if metadata is not None and metadata.get('area') is not None:
        if args.area == 1.0:  # Default value was not changed
            area_to_use = metadata['area']
            logger.info(f"Using area from DTA metadata: {area_to_use:.4f} cm^2")
        else:
            logger.info(f"Using explicitly specified area: {area_to_use:.4f} cm^2 "
                        f"(metadata: {metadata['area']:.4f} cm^2)")

    analyze_oxide_layer(
        frequencies, Z,
        epsilon_r=args.epsilon_r,
        area_cm2=area_to_use,
        fit_result=fitted_result
    )
