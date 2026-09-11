"""
Oxide layer analysis handler for the EIS CLI.

- run_oxide_analysis: oxide thickness estimation from capacitance, or the
  inverse (permittivity from a known thickness) when --thickness is given
"""

import argparse
import logging
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ..logging import log_separator
from ...analysis import (
    analyze_oxide_layer,
    estimate_permittivity,
    OxideAnalysisResult,
)
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
        oxide = estimate_permittivity(
            frequencies, Z,
            thickness_nm=args.thickness,
            area_cm2=area_to_use,
            fit_result=fitted_result
        )
        _print_oxide_section(oxide, "Permittivity estimation from known thickness")
    else:
        oxide = analyze_oxide_layer(
            frequencies, Z,
            epsilon_r=args.epsilon_r if args.epsilon_r is not None
                      else DEFAULT_EPSILON_R,
            area_cm2=area_to_use,
            fit_result=fitted_result
        )
        _print_oxide_section(oxide, "Oxide layer analysis")


def _print_candidates(oxide: OxideAnalysisResult) -> None:
    """List every capacitive element the circuit offered, and which one won.

    The dominant-element choice is a heuristic - the largest R may equally
    belong to a charge-transfer process - so the alternatives are printed
    beside it rather than left implicit.
    """
    logger.info(f"Found {len(oxide.candidates)} capacitive element(s):")
    for i, e in enumerate(oxide.candidates, 1):
        R_str = (f"R = {e['R']:.1f} Ω" if e['R'] is not None
                 else "no parallel R")
        if e['type'] == 'Q':
            logger.info(f"  [{i}] Q: {R_str}, "
                        f"Q = {e['Q']:.3e}, n = {e['n']:.3f}")
        elif e['type'] == 'CC':
            logger.info(f"  [{i}] CC: C_inf = {e['C_inf']:.3e} F, "
                        f"ΔC = {e['dC']:.3e} F "
                        f"(C_s = {e['C_static']:.3e} F), "
                        f"tau = {e['tau']:.2e} s, alpha = {e['alpha']:.3f}")
        elif e['type'] == 'DQ':
            logger.info(f"  [{i}] DQ: R_pol = {e['R']:.3e} Ω, "
                        f"C_eff = {e['C']:.3e} F, n = {e['n_power']:.3f}, "
                        f"tau = {e['tau_min']:.2e}..{e['tau']:.2e} s "
                        f"(U = {e['U']:.2f})")
        else:
            tau_str = (f", tau = {e['tau']:.2e} s" if e['tau'] is not None
                       else "")
            logger.info(f"  [{i}] {e['type']}: {R_str}, "
                        f"C = {e['C']:.3e} F{tau_str}")

    params = oxide.element_params
    if oxide.element_type == 'CC':
        logger.info(f"Dominant element: CC with C = {params['C']:.3e} F")
    elif oxide.element_R is not None:
        logger.info(f"Dominant element: {oxide.element_type} with "
                    f"R = {oxide.element_R:.1f} Ω")
    else:
        size = params.get('C', params.get('Q'))
        logger.info(f"Dominant element: {oxide.element_type} with "
                    f"C = {size:.3e} F (no parallel resistance to rank by)")
    if oxide.selection_reason:
        logger.info(oxide.selection_reason)


def _print_oxide_section(oxide: Optional[OxideAnalysisResult],
                         title: str) -> None:
    """Print the oxide analysis section from the result."""
    log_separator()
    logger.info(title)
    log_separator()

    if oxide is None:
        logger.warning("Could not extract capacitance from data")
        log_separator()
        return

    if oxide.mode == 'circuit':
        _print_candidates(oxide)
        logger.info("")
        logger.info("Results:")
        logger.info(f"  Element type:       {oxide.element_type}")
        if oxide.element_R is not None:
            logger.info(f"  Resistance:         {oxide.element_R:.1f} Ω")
        elif oxide.element_type == 'CC':
            logger.info("  Resistance:         n/a (blocking dielectric)")
        else:
            logger.info("  Resistance:         n/a (no parallel R "
                        "in the circuit)")

        params = oxide.element_params
        cc_suffix = ""
        if oxide.element_type == 'CC':
            logger.info(f"  Broadening alpha:   {params['alpha']:.3f}")
            logger.info(f"  C_inf / ΔC:         {params['C_inf']:.3e} F / "
                        f"{params['dC']:.3e} F")
            cc_suffix = ("  (C_inf, high-frequency limit)"
                         if params['C_regime'] == 'high_frequency'
                         else "  (static, C_inf + ΔC)")
        logger.info(f"  Capacitance:        {oxide.capacitance:.3e} F{cc_suffix}")
        if oxide.capacitance_brug is not None:
            logger.info(f"  C (Brug, 2D):       {oxide.capacitance_brug:.3e} F "
                        f"(comparison; primary value is Hsu-Mansfeld, 3D)")
        logger.info(f"  Specific cap.:      "
                    f"{oxide.capacitance_specific * 1e6:.2f} µF/cm²")
        if oxide.element_tau is not None and oxide.element_tau > 0:
            logger.info(f"  Time constant:      {oxide.element_tau:.3e} s")
            logger.info(f"  Char. frequency:    "
                        f"{1 / (2 * np.pi * oxide.element_tau):.2e} Hz")
        else:
            # A capacitance with no parallel resistance has no RC time
            # constant - there is nothing to discharge through
            logger.info("  Time constant:      n/a (no parallel R)")
    else:
        logger.info("Mode: High-frequency estimate (simplified)")
        if oxide.n_hf_points is not None:
            logger.info(f"  Median over {oxide.n_hf_points} point(s) "
                        f"in the top frequency decade")
        logger.info(f"  Capacitance:        {oxide.capacitance:.3e} F")
        logger.info(f"  Specific cap.:      "
                    f"{oxide.capacitance_specific * 1e6:.2f} µF/cm²")

    if oxide.permittivity is not None:
        logger.info(f"  Known thickness:    {oxide.thickness_nm:.1f} nm")
        logger.info(f"  Permittivity ε_r:   {oxide.permittivity:.1f}")
        if oxide.permittivity_brug is not None:
            # .3g, not .1f: when n is far from 1 the two CPE models diverge by
            # orders of magnitude and a fixed-point format would print "0.0"
            logger.info(f"  ε_r (Brug):         {oxide.permittivity_brug:.3g} "
                        f"(2D model, for comparison)")
        logger.info(f"  (area={oxide.area_cm2} cm²)")
    else:
        logger.info(f"  Oxide thickness:    {oxide.thickness_nm:.1f} nm")
        if oxide.thickness_brug_nm is not None:
            logger.info(f"  Thickness (Brug):   {oxide.thickness_brug_nm:.1f} nm "
                        f"(2D model, for comparison)")
        logger.info(f"  (assuming ε_r={oxide.epsilon_r}, "
                    f"area={oxide.area_cm2} cm²)")

    for warning in oxide.warnings:
        logger.warning(warning)
    log_separator()
