"""
Oxide layer thickness estimation from EIS data.

Simplified implementation that finds the dominant Voigt or Q element
and estimates oxide thickness from its capacitance.
"""

import numpy as np
import logging
from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from numpy.typing import NDArray

from ..fitting.circuit import FitResult
from ..fitting.circuit_elements import R, C, Q, K
from ..fitting.circuit_builder import Series, Parallel
from .config import (
    EPSILON_0,
    CPE_N_RELIABLE_MIN,
    HF_ESTIMATE_DECADE_FACTOR,
    HF_C_SPREAD_MAX_RATIO,
)

logger = logging.getLogger(__name__)


@dataclass
class OxideAnalysisResult:
    """Result of oxide layer analysis."""
    capacitance: float          # Effective capacitance [F]
    capacitance_specific: float # Specific capacitance [F/cm²]
    thickness_nm: float         # Oxide thickness [nm]
    element_type: str           # 'C', 'K', 'Q', or 'estimate' (HF fallback)
    element_R: Optional[float]  # Associated resistance [Ω]
    element_tau: Optional[float] # Time constant [s]
    element_params: Dict[str, float]  # All element parameters
    # Brug (2D) comparison values - set only for a dominant Q element
    # when a series resistance is present in the circuit
    capacitance_brug: Optional[float] = None          # Brug C_eff [F]
    capacitance_specific_brug: Optional[float] = None # Brug C_eff/area [F/cm²]
    thickness_brug_nm: Optional[float] = None         # Thickness from Brug C [nm]


def _find_parallel_rc_elements(circuit) -> List[Dict[str, Any]]:
    """
    Find all parallel R-C, R-Q combinations and K elements in circuit.

    Returns list of dicts with keys: 'type', 'R', 'C' or 'Q'/'n', 'tau'
    """
    results = []

    def traverse(node):
        if isinstance(node, Parallel):
            # Check if this parallel contains R with C/Q
            elements = node.elements
            R_elem = None
            cap_elem = None

            for elem in elements:
                if isinstance(elem, R):
                    if R_elem is not None:
                        logger.warning("Multiple R elements in one parallel "
                                       "combination - using the last one")
                    R_elem = elem
                elif isinstance(elem, (C, Q)):
                    if cap_elem is not None:
                        logger.warning("Multiple C/Q elements in one parallel "
                                       "combination - using the last one")
                    cap_elem = elem

            if R_elem is not None and cap_elem is not None:
                R_val = R_elem.params[0]
                if isinstance(cap_elem, C):
                    C_val = cap_elem.params[0]
                    results.append({
                        'type': 'C',
                        'R': R_val,
                        'C': C_val,
                        'tau': R_val * C_val
                    })
                elif isinstance(cap_elem, Q):
                    Q_val = cap_elem.params[0]
                    n_val = cap_elem.params[1]
                    results.append({
                        'type': 'Q',
                        'R': R_val,
                        'Q': Q_val,
                        'n': n_val,
                        'tau': None  # Will be computed later
                    })

            # Continue traversing children
            for elem in elements:
                traverse(elem)

        elif isinstance(node, Series):
            for elem in node.elements:
                traverse(elem)

        elif isinstance(node, K):
            # K element directly provides R and tau
            R_val = node.params[0]
            tau_val = node.params[1]
            if R_val <= 0:
                # C = tau/R is undefined; a non-positive R would be dropped
                # by the dominant-element filter anyway
                logger.warning(f"K element with non-positive R = {R_val:g} Ω - skipping")
                return
            C_val = tau_val / R_val
            results.append({
                'type': 'K',
                'R': R_val,
                'C': C_val,
                'tau': tau_val
            })

    traverse(circuit)
    return results


def _estimate_cpe_capacitance(Q_val: float, n: float, R_val: float) -> float:
    """
    Estimate effective capacitance of Q (CPE) element.

    Hsu-Mansfeld formula (requires the parallel resistance R):
        C_eff = (R × Q)^(1/n) / R    (via τ = (R × Q)^(1/n))

    Assumes a normal (3D, through-layer) distribution of time
    constants — appropriate for oxide layers. For a surface (2D)
    distribution the Brug (1984) formula would apply instead,
    which also involves the series resistance:
    C = Q^(1/n) × (1/Rs + 1/Rct)^((n-1)/n).

    Reference: Hsu & Mansfeld, Corrosion 57, 747 (2001).
    """
    C_eff = (R_val * Q_val) ** (1.0 / n) / R_val
    logger.debug(f"Q C_eff (Hsu-Mansfeld): {C_eff:.3e} F")
    return C_eff


def _estimate_cpe_capacitance_brug(
    Q_val: float, n: float, R_ct: float, R_s: float
) -> float:
    """
    Estimate effective capacitance of a Q (CPE) element by the Brug formula:

        C_eff = Q^(1/n) × (1/Rs + 1/Rct)^((n-1)/n)

    Assumes a surface (2D, lateral) distribution of time constants.
    Reported alongside the Hsu-Mansfeld (3D) value as a comparison;
    the spread between the two brackets the model uncertainty of C_eff.

    Reference: Brug et al., J. Electroanal. Chem. 176, 275 (1984).
    """
    C_eff = Q_val ** (1.0 / n) * (1.0 / R_s + 1.0 / R_ct) ** ((n - 1.0) / n)
    logger.debug(f"Q C_eff (Brug): {C_eff:.3e} F")
    return C_eff


def _find_series_resistance(circuit) -> Optional[float]:
    """
    Sum of R elements on the series path of the circuit (outside any
    parallel combination) — the ohmic/electrolyte resistance Rs needed
    by the Brug formula.

    Returns None if no such element exists or the sum is not positive.
    """
    total = 0.0
    found = False

    def traverse(node):
        nonlocal total, found
        if isinstance(node, R):
            total += node.params[0]
            found = True
        elif isinstance(node, Series):
            for elem in node.elements:
                traverse(elem)
        # Parallel, K, C, Q: not part of the series path

    traverse(circuit)
    return total if found and total > 0 else None


def _extract_capacitance(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    area_cm2: float,
    fit_result: Optional[FitResult]
) -> Optional[Dict[str, Any]]:
    """
    Extract effective capacitance of the dominant capacitive element.

    Shared core of analyze_oxide_layer() and estimate_permittivity():
    element selection, capacitance estimation, and related logging.
    Deliberately does NOT compute or log thickness — each caller logs
    only the quantity it actually derives from the capacitance.

    Returns
    -------
    extracted : dict or None
        Keys: 'C_eff' [F], 'C_specific' [F/cm²], 'element_type',
        'element_R', 'element_tau', 'element_params'.
        None if capacitance could not be extracted.
    """
    # === Mode 1: From fitted circuit (preferred) ===
    if fit_result is not None:
        circuit = fit_result.circuit

        # Find all parallel R-C/Q combinations and K elements
        elements = _find_parallel_rc_elements(circuit)

        if not elements:
            logger.warning("No Voigt (R||C), K, or R||Q elements found in circuit")
            logger.warning("Falling back to high-frequency estimate...")
            fit_result = None
        else:
            # List all candidates so the dominant-element choice can be
            # verified (largest R may also be a charge-transfer process)
            logger.info(f"Found {len(elements)} capacitive element(s):")
            for i, e in enumerate(elements, 1):
                if e['type'] == 'Q':
                    logger.info(f"  [{i}] Q: R = {e['R']:.1f} Ω, "
                                f"Q = {e['Q']:.3e}, n = {e['n']:.3f}")
                else:
                    logger.info(f"  [{i}] {e['type']}: R = {e['R']:.1f} Ω, "
                                f"C = {e['C']:.3e} F, tau = {e['tau']:.2e} s")

            # Find element with largest R (dominant barrier = compact oxide)
            elements_with_R = [e for e in elements if e.get('R') is not None and e['R'] > 0]

            if not elements_with_R:
                logger.warning("No elements with valid resistance found")
                fit_result = None
            else:
                dominant = max(elements_with_R, key=lambda e: e['R'])

                logger.info(f"Dominant element: {dominant['type']} with R = {dominant['R']:.1f} Ω")
                logger.info("Selection assumes the largest-R element is the compact "
                            "oxide barrier (verify: a charge-transfer process can "
                            "also have the largest R)")

                # Get capacitance
                C_eff_brug = None
                if dominant['type'] in ('C', 'K'):
                    C_eff = dominant['C']
                    tau = dominant['tau']
                else:  # Q
                    if dominant['n'] < CPE_N_RELIABLE_MIN:
                        logger.warning(
                            f"CPE exponent n = {dominant['n']:.3f} < "
                            f"{CPE_N_RELIABLE_MIN}: effective capacitance is not "
                            "well-defined; thickness estimate may be unreliable")
                    C_eff = _estimate_cpe_capacitance(
                        dominant['Q'], dominant['n'], dominant['R']
                    )
                    # Estimate tau from R and C_eff
                    tau = dominant['R'] * C_eff
                    dominant['tau'] = tau

                    # Brug (2D) comparison estimate - needs series resistance
                    R_s = _find_series_resistance(circuit)
                    if R_s is not None:
                        C_eff_brug = _estimate_cpe_capacitance_brug(
                            dominant['Q'], dominant['n'], dominant['R'], R_s
                        )
                    else:
                        logger.info("No series R element in circuit - "
                                    "Brug (2D) estimate not available")

                C_specific = C_eff / area_cm2
                C_specific_brug = (C_eff_brug / area_cm2
                                   if C_eff_brug is not None else None)

                # Log element info (thickness/permittivity logged by caller)
                logger.info("")
                logger.info("Results:")
                logger.info(f"  Element type:       {dominant['type']}")
                logger.info(f"  Resistance:         {dominant['R']:.1f} Ω")
                logger.info(f"  Capacitance:        {C_eff:.3e} F")
                if C_eff_brug is not None:
                    logger.info(f"  C (Brug, 2D):       {C_eff_brug:.3e} F "
                                f"(comparison; primary value is Hsu-Mansfeld, 3D)")
                logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")
                logger.info(f"  Time constant:      {tau:.3e} s")
                logger.info(f"  Char. frequency:    {1/(2*np.pi*tau):.2e} Hz")

                return {
                    'C_eff': C_eff,
                    'C_specific': C_specific,
                    'C_eff_brug': C_eff_brug,
                    'C_specific_brug': C_specific_brug,
                    'element_type': dominant['type'],
                    'element_R': dominant['R'],
                    'element_tau': tau,
                    'element_params': dict(dominant)
                }

    # === Mode 2: Fallback - high-frequency estimate ===
    logger.info("Mode: High-frequency estimate (simplified)")
    logger.warning("For better accuracy, provide fitted circuit via fit_result")
    logger.warning("For multilayer (series) systems the high-frequency estimate "
                   "yields the series combination of layer capacitances")

    # Estimate C from imaginary impedance, C = -1 / (ω × Z''), as the
    # median over capacitive points in the top frequency decade
    high_freq_idx = np.argmax(frequencies)
    f_max = frequencies[high_freq_idx]
    decade_mask = frequencies >= f_max / HF_ESTIMATE_DECADE_FACTOR
    capacitive_mask = decade_mask & (Z.imag < -1e-10)

    if np.any(capacitive_mask):
        omega = 2 * np.pi * frequencies[capacitive_mask]
        C_values = -1 / (omega * Z.imag[capacitive_mask])
        C_estimate = float(np.median(C_values))
        logger.info(f"  Median over {C_values.size} point(s) in the top frequency decade")

        # C_i is frequency-independent only when the capacitance dominates
        # (ωRC ≫ 1); a large spread means that assumption does not hold
        spread = float(np.max(C_values) / np.min(C_values))
        if spread > HF_C_SPREAD_MAX_RATIO:
            logger.warning(
                f"C estimates vary by factor {spread:.2f} across the top "
                f"frequency decade (ωRC ≫ 1 may not hold); "
                "estimate may be unreliable")
    else:
        # No capacitive point in the top decade: fall back to the single
        # highest-frequency point (original pre-0.16.16 behavior)
        Z_imag_hf = Z[high_freq_idx].imag

        if abs(Z_imag_hf) < 1e-10:
            logger.error("Imaginary impedance too small at high frequency")
            return None

        if Z_imag_hf > 0:
            logger.warning("Positive imaginary impedance (inductive) - result may be invalid")

        omega_hf = 2 * np.pi * f_max
        C_estimate = -1 / (omega_hf * Z_imag_hf)

    C_specific = C_estimate / area_cm2

    logger.info(f"  Capacitance:        {C_estimate:.3e} F")
    logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")

    return {
        'C_eff': C_estimate,
        'C_specific': C_specific,
        'C_eff_brug': None,
        'C_specific_brug': None,
        'element_type': 'estimate',
        'element_R': None,
        'element_tau': None,
        'element_params': {}
    }


def analyze_oxide_layer(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    epsilon_r: float = 22.0,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[OxideAnalysisResult]:
    """
    Estimate oxide layer thickness from dominant capacitive element.

    Finds the Voigt (R||C), K, or R||Q element with the largest resistance
    (dominant barrier) and calculates oxide thickness from its capacitance.

    Parameters
    ----------
    frequencies : ndarray
        Measurement frequencies [Hz]
    Z : ndarray
        Complex impedance [Ω]
    epsilon_r : float, optional
        Relative permittivity of oxide (default: 22 for ZrO₂)
    area_cm2 : float, optional
        Electrode area [cm²] (default: 1.0)
    fit_result : FitResult, optional
        Result from fit_equivalent_circuit(). If None, uses simple
        high-frequency estimate (less accurate).

    Returns
    -------
    result : OxideAnalysisResult or None
        Analysis result with capacitance and thickness, or None if failed.

    Notes
    -----
    Thickness formula (parallel plate capacitor model):
        d = ε₀ × εᵣ / C_specific

    For Q elements, effective capacitance is estimated using the
    Hsu-Mansfeld formula: C_eff = (R × Q)^(1/n) / R
    (assumes a normal/3D distribution of time constants).
    When the circuit also contains a series resistance, the Brug (1984)
    formula (surface/2D distribution) is evaluated as well and reported
    in capacitance_brug / thickness_brug_nm for comparison; the spread
    between the two estimates brackets the model uncertainty.
    See doc/OXIDE_ANALYSIS_GUIDE.md for the 2D vs 3D discussion.

    Examples
    --------
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
    >>> oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
    >>> print(f"Thickness: {oxide.thickness_nm:.1f} nm")
    """
    logger.info("=" * 50)
    logger.info("Oxide layer analysis")
    logger.info("=" * 50)

    extracted = _extract_capacitance(frequencies, Z, area_cm2, fit_result)
    if extracted is None:
        return None

    # Calculate thickness (parallel plate capacitor model)
    C_specific = extracted['C_specific']
    d_cm = EPSILON_0 * epsilon_r / C_specific
    d_nm = d_cm * 1e7

    # Brug (2D) comparison thickness, when available
    C_specific_brug = extracted['C_specific_brug']
    d_brug_nm = None
    if C_specific_brug is not None:
        d_brug_nm = EPSILON_0 * epsilon_r / C_specific_brug * 1e7

    logger.info(f"  Oxide thickness:    {d_nm:.1f} nm")
    if d_brug_nm is not None:
        logger.info(f"  Thickness (Brug):   {d_brug_nm:.1f} nm "
                    f"(2D model, for comparison)")
    logger.info(f"  (assuming ε_r={epsilon_r}, area={area_cm2} cm²)")
    logger.info("=" * 50)

    return OxideAnalysisResult(
        capacitance=extracted['C_eff'],
        capacitance_specific=C_specific,
        thickness_nm=d_nm,
        element_type=extracted['element_type'],
        element_R=extracted['element_R'],
        element_tau=extracted['element_tau'],
        element_params=extracted['element_params'],
        capacitance_brug=extracted['C_eff_brug'],
        capacitance_specific_brug=C_specific_brug,
        thickness_brug_nm=d_brug_nm
    )


def estimate_permittivity(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    thickness_nm: float,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[float]:
    """
    Estimate relative permittivity from known oxide thickness.

    Inverse of analyze_oxide_layer(): given thickness, calculates epsilon_r.

    Parameters
    ----------
    frequencies : ndarray
        Measurement frequencies [Hz]
    Z : ndarray
        Complex impedance [Ω]
    thickness_nm : float
        Known oxide layer thickness [nm]
    area_cm2 : float, optional
        Electrode area [cm²] (default: 1.0)
    fit_result : FitResult, optional
        Result from fit_equivalent_circuit(). If None, uses simple
        high-frequency estimate (less accurate).

    Returns
    -------
    epsilon_r : float or None
        Estimated relative permittivity, or None if failed.

    Notes
    -----
    Formula (from parallel plate capacitor model):
        ε_r = d × C_specific / ε₀

    Examples
    --------
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
    >>> eps_r = estimate_permittivity(freq, Z, thickness_nm=20, fit_result=result)
    >>> print(f"Permittivity: {eps_r:.1f}")
    """
    logger.info("=" * 50)
    logger.info("Permittivity estimation from known thickness")
    logger.info("=" * 50)

    # Get capacitance using the same element-selection logic as
    # analyze_oxide_layer (no thickness is computed or logged here)
    extracted = _extract_capacitance(frequencies, Z, area_cm2, fit_result)

    if extracted is None:
        logger.error("Could not extract capacitance from data")
        return None

    # Calculate permittivity from thickness and capacitance
    # d = ε₀ × εᵣ / C_specific  =>  εᵣ = d × C_specific / ε₀
    d_cm = thickness_nm * 1e-7  # nm -> cm
    C_specific = extracted['C_specific']
    epsilon_r = d_cm * C_specific / EPSILON_0

    logger.info(f"  Known thickness:    {thickness_nm:.1f} nm")
    logger.info(f"  Permittivity ε_r:   {epsilon_r:.1f}")
    logger.info(f"  (area={area_cm2} cm²)")
    logger.info("=" * 50)

    return epsilon_r


__all__ = ['analyze_oxide_layer', 'estimate_permittivity', 'OxideAnalysisResult']
