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
from .config import EPSILON_0

logger = logging.getLogger(__name__)


@dataclass
class OxideAnalysisResult:
    """Result of oxide layer analysis."""
    capacitance: float          # Effective capacitance [F]
    capacitance_specific: float # Specific capacitance [F/cm²]
    thickness_nm: float         # Oxide thickness [nm]
    element_type: str           # 'C', 'K', or 'Q'
    element_R: Optional[float]  # Associated resistance [Ω]
    element_tau: Optional[float] # Time constant [s]
    element_params: Dict[str, float]  # All element parameters


def _find_parallel_rc_elements(circuit) -> List[Dict[str, Any]]:
    """
    Find all parallel R-C, R-Q combinations and K elements in circuit.

    Returns list of dicts with keys: 'type', 'R', 'C' or 'Q'/'n', 'tau'
    """
    results = []

    def traverse(node, parent_is_parallel=False, sibling_R=None):
        if isinstance(node, Parallel):
            # Check if this parallel contains R with C/Q
            elements = node.elements
            R_elem = None
            cap_elem = None

            for elem in elements:
                if isinstance(elem, R):
                    R_elem = elem
                elif isinstance(elem, (C, Q)):
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
                traverse(elem, parent_is_parallel=True)

        elif isinstance(node, Series):
            for elem in node.elements:
                traverse(elem, parent_is_parallel=False)

        elif isinstance(node, K):
            # K element directly provides R and tau
            R_val = node.params[0]
            tau_val = node.params[1]
            C_val = tau_val / R_val
            results.append({
                'type': 'K',
                'R': R_val,
                'C': C_val,
                'tau': tau_val
            })

    traverse(circuit)
    return results


def _estimate_cpe_capacitance(
    Q: float,
    n: float,
    R: Optional[float],
    frequencies: Optional[NDArray] = None,
    Z: Optional[NDArray] = None
) -> float:
    """
    Estimate effective capacitance of Q (CPE) element.

    Method 1 (preferred): Brug formula when R is known
        C_eff = (R × Q)^(1/n) / R

    Method 2 (fallback): From Z'' maximum frequency
        C_eff = Q × ω_max^(n-1)

    Method 3 (last resort): Simple approximation at 1 kHz
        C_eff ≈ Q × (2π × 1000)^(n-1)
    """
    # Method 1: Brug formula (most accurate when R is known)
    if R is not None and R > 0:
        C_eff = (R * Q) ** (1.0 / n) / R
        logger.debug(f"Q C_eff (Brug): {C_eff:.3e} F")
        return C_eff

    # Method 2: From Z'' maximum
    if frequencies is not None and Z is not None:
        Z_imag = -Z.imag  # -Z'' (positive for capacitive)
        if np.any(Z_imag > 0):
            max_idx = np.argmax(Z_imag)
            omega_max = 2 * np.pi * frequencies[max_idx]
            C_eff = Q * omega_max ** (n - 1)
            logger.debug(f"Q C_eff (from Z'' max at {frequencies[max_idx]:.1f} Hz): {C_eff:.3e} F")
            return C_eff

    # Method 3: Approximation at 1 kHz
    omega_ref = 2 * np.pi * 1000  # 1 kHz reference
    C_eff = Q * omega_ref ** (n - 1)
    logger.debug(f"Q C_eff (1 kHz approx): {C_eff:.3e} F")
    return C_eff


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

    For Q elements, effective capacitance is estimated using:
        - Brug formula: C_eff = (R × Q)^(1/n) / R  (when R is known)
        - From Z'' maximum frequency (fallback)

    Examples
    --------
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
    >>> oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
    >>> print(f"Thickness: {oxide.thickness_nm:.1f} nm")
    """
    logger.info("=" * 50)
    logger.info("Oxide layer analysis")
    logger.info("=" * 50)

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
            logger.info(f"Found {len(elements)} capacitive element(s)")

            # Find element with largest R (dominant barrier = compact oxide)
            elements_with_R = [e for e in elements if e.get('R') is not None and e['R'] > 0]

            if not elements_with_R:
                logger.warning("No elements with valid resistance found")
                fit_result = None
            else:
                dominant = max(elements_with_R, key=lambda e: e['R'])

                logger.info(f"Dominant element: {dominant['type']} with R = {dominant['R']:.1f} Ω")

                # Get capacitance
                if dominant['type'] in ('C', 'K'):
                    C_eff = dominant['C']
                    tau = dominant['tau']
                else:  # Q
                    C_eff = _estimate_cpe_capacitance(
                        dominant['Q'], dominant['n'], dominant['R'],
                        frequencies, Z
                    )
                    # Estimate tau from R and C_eff
                    tau = dominant['R'] * C_eff
                    dominant['tau'] = tau

                # Calculate thickness
                C_specific = C_eff / area_cm2
                d_cm = EPSILON_0 * epsilon_r / C_specific
                d_nm = d_cm * 1e7

                # Log results
                logger.info("")
                logger.info("Results:")
                logger.info(f"  Element type:       {dominant['type']}")
                logger.info(f"  Resistance:         {dominant['R']:.1f} Ω")
                logger.info(f"  Capacitance:        {C_eff:.3e} F")
                logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")
                logger.info(f"  Time constant:      {tau:.3e} s")
                logger.info(f"  Char. frequency:    {1/(2*np.pi*tau):.2e} Hz")
                logger.info(f"  Oxide thickness:    {d_nm:.1f} nm")
                logger.info(f"  (assuming ε_r={epsilon_r}, area={area_cm2} cm²)")
                logger.info("=" * 50)

                return OxideAnalysisResult(
                    capacitance=C_eff,
                    capacitance_specific=C_specific,
                    thickness_nm=d_nm,
                    element_type=dominant['type'],
                    element_R=dominant['R'],
                    element_tau=tau,
                    element_params=dominant
                )

    # === Mode 2: Fallback - high-frequency estimate ===
    logger.info("Mode: High-frequency estimate (simplified)")
    logger.warning("For better accuracy, provide fitted circuit via fit_result")

    # Estimate C from imaginary impedance at highest frequency
    # C = -1 / (ω × Z_imag)
    high_freq_idx = np.argmax(frequencies)
    omega_hf = 2 * np.pi * frequencies[high_freq_idx]
    Z_imag_hf = Z[high_freq_idx].imag

    if abs(Z_imag_hf) < 1e-10:
        logger.error("Imaginary impedance too small at high frequency")
        return None

    if Z_imag_hf > 0:
        logger.warning("Positive imaginary impedance (inductive) - result may be invalid")

    C_estimate = -1 / (omega_hf * Z_imag_hf)
    C_specific = C_estimate / area_cm2
    d_cm = EPSILON_0 * epsilon_r / C_specific
    d_nm = d_cm * 1e7

    logger.info(f"  Capacitance:        {C_estimate:.3e} F")
    logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")
    logger.info(f"  Oxide thickness:    {d_nm:.1f} nm")
    logger.info(f"  (assuming ε_r={epsilon_r}, area={area_cm2} cm²)")
    logger.info("=" * 50)

    return OxideAnalysisResult(
        capacitance=C_estimate,
        capacitance_specific=C_specific,
        thickness_nm=d_nm,
        element_type='estimate',
        element_R=None,
        element_tau=None,
        element_params={}
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

    # Get capacitance using same logic as analyze_oxide_layer
    # but with dummy epsilon_r (we'll calculate the real one)
    oxide_result = analyze_oxide_layer(
        frequencies, Z,
        epsilon_r=1.0,  # dummy value
        area_cm2=area_cm2,
        fit_result=fit_result
    )

    if oxide_result is None:
        logger.error("Could not extract capacitance from data")
        return None

    # Calculate permittivity from thickness and capacitance
    # d = ε₀ × εᵣ / C_specific  =>  εᵣ = d × C_specific / ε₀
    d_cm = thickness_nm * 1e-7  # nm -> cm
    C_specific = oxide_result.capacitance_specific
    epsilon_r = d_cm * C_specific / EPSILON_0

    logger.info("")
    logger.info("Results:")
    logger.info(f"  Known thickness:    {thickness_nm:.1f} nm")
    logger.info(f"  Capacitance:        {oxide_result.capacitance:.3e} F")
    logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")
    logger.info(f"  Permittivity ε_r:   {epsilon_r:.1f}")
    logger.info(f"  (area={area_cm2} cm²)")
    logger.info("=" * 50)

    return epsilon_r


__all__ = ['analyze_oxide_layer', 'estimate_permittivity', 'OxideAnalysisResult']
