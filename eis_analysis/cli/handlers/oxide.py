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
from ...analysis.config import DEFAULT_AREA_CM2, DEFAULT_EPSILON_R
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

    # An explicit --area always wins; metadata only fills in when the flag
    # was omitted (args.area is None), so that --area 1.0 is honored like
    # any other value rather than being mistaken for the default
    area_from_metadata = metadata.get('area') if metadata is not None else None
    if args.area is not None:
        area_to_use = args.area
        if area_from_metadata is not None:
            logger.info(f"Using explicitly specified area: {area_to_use:.4f} cm^2 "
                        f"(metadata: {area_from_metadata:.4f} cm^2)")
    elif area_from_metadata is not None:
        area_to_use = area_from_metadata
        logger.info(f"Using area from DTA metadata: {area_to_use:.4f} cm^2")
    else:
        area_to_use = DEFAULT_AREA_CM2

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
