"""
Oxide layer analysis handler for the EIS CLI.

- run_oxide_analysis: oxide thickness estimation from capacitance, or the
  inverse (permittivity from a known thickness) when --thickness is given
"""

import argparse
import logging
from typing import Optional

from numpy.typing import NDArray

from ...analysis import analyze_oxide_layer, estimate_permittivity
from ...analysis.config import DEFAULT_EPSILON_R
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

    With --thickness the analysis runs in reverse: the thickness becomes
    the input and the relative permittivity the estimated quantity.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: analyze_oxide, epsilon_r, thickness, area)
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

    if args.thickness is not None:
        if args.epsilon_r is not None:
            logger.warning(f"--epsilon-r {args.epsilon_r:g} ignored: --thickness "
                           f"given, permittivity is the estimated quantity")
        estimate_permittivity(
            frequencies, Z,
            thickness_nm=args.thickness,
            area_cm2=area_to_use,
            fit_result=fitted_result
        )
    else:
        analyze_oxide_layer(
            frequencies, Z,
            epsilon_r=args.epsilon_r if args.epsilon_r is not None
                      else DEFAULT_EPSILON_R,
            area_cm2=area_to_use,
            fit_result=fitted_result
        )
