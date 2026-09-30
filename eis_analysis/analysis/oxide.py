"""
Oxide layer thickness estimation from EIS data.

Picks the dielectric element of a fitted circuit (C, Q, K, CC, DQ or YG
with admittance ~ omega^n, n near 1; the largest parallel resistance wins),
reads its capacitance - a Q through Hsu-Mansfeld, with Brug for comparison -
and turns it into a thickness or a permittivity by the parallel-plate model.
Without a fitted circuit the capacitance is a high-frequency spectral estimate.

Element search and selection live in oxide_elements, notes on what the
measured window supports and the spectral estimate in oxide_window.
"""

import numpy as np
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field
from numpy.typing import NDArray

from ..fitting.circuit import FitResult
from .config import (
    EPSILON_0,
    DEFAULT_EPSILON_R,
    CPE_N_RELIABLE_MIN,
    BRUG_HM_DIVERGENCE_MAX,
)
from .oxide_elements import (
    _find_capacitive_elements, _select_dielectric_element, _q_capacitance,
    _add_penetration_depth)
from .oxide_window import _window_notes, _hf_capacitance_estimate
from .oxide_power_law import _power_law

@dataclass
class OxideAnalysisResult:
    """
    Result of oxide layer analysis.

    Besides the numbers, this carries how they were arrived at: every
    capacitive candidate the circuit offered, why one of them was picked,
    and whether the capacitance came from the fit at all (`mode`). Without
    those a reader cannot check the dominant-element choice, which is a
    heuristic - the largest-R element may equally be a charge-transfer
    process.
    """
    capacitance: float          # Effective capacitance [F]
    capacitance_specific: float # Specific capacitance [F/cm²]
    thickness_nm: float         # Oxide thickness [nm]
    element_type: str           # 'C', 'K', 'Q', 'CC', 'DQ', 'YG', or 'estimate' (HF fallback)
    element_R: Optional[float]  # Parallel resistance [Ω] (None: no DC path beside it)
    element_tau: Optional[float] # Time constant [s]
    element_params: Dict[str, Any]  # All element parameters (+ 'type', flags)
    # Brug (2D) comparison values - set only for a dominant Q element
    # when a series resistance is present in the circuit
    capacitance_brug: Optional[float] = None          # Brug C_eff [F]
    capacitance_specific_brug: Optional[float] = None # Brug C_eff/area [F/cm²]
    thickness_brug_nm: Optional[float] = None         # Thickness from Brug C [nm]
    # Inverse mode - set only by estimate_permittivity(), where the
    # thickness is the input and the permittivity the derived quantity
    permittivity: Optional[float] = None       # ε_r from known thickness
    permittivity_brug: Optional[float] = None  # ε_r from Brug C_eff
    # Power-law model (Hirschorn-Orazem) comparison values - set only for a
    # dominant Q element when rho_delta_ohm_cm was given
    thickness_pl_nm: Optional[float] = None    # d from the power-law model [nm]
    permittivity_pl: Optional[float] = None    # ε_r from the power-law model
    rho_delta_ohm_cm: Optional[float] = None   # rho_delta it assumed [Ω·cm]

    # How the capacitance was obtained
    mode: str = 'circuit'       # 'circuit' (from the fit) or 'hf_estimate'
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    selection_reason: str = ''  # why the dominant element was picked
    n_hf_points: Optional[int] = None  # points behind the HF median estimate
    epsilon_r: Optional[float] = None  # value assumed by analyze_oxide_layer
    area_cm2: float = 1.0
    warnings: List[str] = field(default_factory=list)


def _extract_capacitance(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    area_cm2: float,
    fit_result: Optional[FitResult],
    warnings: List[str]
) -> Optional[Dict[str, Any]]:
    """
    Extract effective capacitance of the dominant capacitive element.

    Shared core of analyze_oxide_layer() and estimate_permittivity():
    element selection and capacitance estimation. Deliberately does NOT
    compute thickness — each caller derives only its own quantity from the
    capacitance. Caveats are appended to `warnings`.

    Returns
    -------
    extracted : dict or None
        Keys: 'C_eff' [F], 'C_specific' [F/cm²], 'C_eff_brug' [F],
        'C_specific_brug' [F/cm²], 'element_type', 'element_R',
        'element_tau', 'element_params'.
        The two Brug (2D) keys are always present but None unless the
        dominant element is a Q and the circuit has a series resistance
        of at least BRUG_RS_MIN_OHM (below that the fit has not identified
        R_s and Brug's R_s^((1-n)/n) scaling makes the value meaningless);
        in fallback mode 'element_R' and 'element_tau' are None too,
        and 'element_R' is None for an element with no DC path beside
        it (a blocking dielectric).
        None if capacitance could not be extracted.
    """
    # === Mode 1: From fitted circuit (preferred) ===
    if fit_result is not None:
        circuit = fit_result.circuit

        # Every C, Q, K, CC, DQ and YG in the circuit, parallel R or not
        elements = _find_capacitive_elements(circuit, frequencies, warnings)

        if not elements:
            warnings.append("No capacitive element (C, Q, K, CC, DQ, YG) found in "
                            "circuit - falling back to high-frequency estimate")
        else:
            # Which element carries the dielectric response? The criterion is
            # physical - admittance rising as omega^n with n close to 1 - not
            # the element type and not its position in the expression. C and K
            # are n = 1 by construction, a CC is n = 1 in both of its limits,
            # and a CPE qualifies only when its exponent is near-ideal: a CPE
            # at n ~ 0.6 describes transport or a distribution of resistivity,
            # not a dielectric.
            usable = []
            for e in elements:
                if e['type'] == 'Q' and not (e['R'] is not None and e['R'] > 0):
                    # Both Hsu-Mansfeld and Brug need the parallel resistance;
                    # without it a CPE cannot be converted to a capacitance
                    warnings.append(
                        f"CPE with n = {e['n']:.3f} has no parallel resistance - "
                        "neither Hsu-Mansfeld nor Brug can convert it to a "
                        "capacitance, so it is not a candidate.")
                else:
                    usable.append(e)

            dielectric = [e for e in usable if e['n'] >= CPE_N_RELIABLE_MIN]
            non_ideal = [e for e in usable if e['n'] < CPE_N_RELIABLE_MIN]

            if dielectric:
                dominant, selection_reason = _select_dielectric_element(
                    dielectric, warnings)
            elif non_ideal:
                # Nothing in the circuit behaves as a dielectric. Better than
                # the spectral fallback - these are at least fitted parameters
                # - but the result is not a dielectric capacitance.
                warnings.append(
                    f"No dielectric element in circuit: the only capacitive "
                    f"element(s) are CPEs with n < {CPE_N_RELIABLE_MIN} "
                    f"(largest n = {max(e['n'] for e in non_ideal):.3f}). Such a "
                    "CPE describes transport or a distribution of resistivity, "
                    "not a dielectric, so the capacitance below - and any "
                    "thickness or permittivity from it - has no dielectric "
                    "meaning. Reported only because it is still a fitted "
                    "parameter, unlike the high-frequency estimate.")
                dominant, selection_reason = _select_dielectric_element(
                    non_ideal, warnings)
            else:
                warnings.append("No convertible capacitive element found - "
                                "falling back to high-frequency estimate")
                dominant = None

            if dominant is not None:

                # Get capacitance
                C_eff_brug = None
                if dominant['type'] in ('C', 'K', 'CC', 'DQ', 'YG'):
                    # For CC this is exact, not an effective capacitance: both
                    # C_s = C_inf + dC (the omega -> 0 limit of C*(omega), for
                    # any alpha) and C_inf (the omega -> inf limit) are model
                    # limits, not fits. No Hsu-Mansfeld / Brug choice arises;
                    # which limit the data support was decided in
                    # _cc_capacitance_regime.
                    C_eff = dominant['C']
                    tau = dominant['tau']
                else:  # Q; one with n < CPE_N_RELIABLE_MIN was warned about above
                    C_eff, C_eff_brug = _q_capacitance(dominant, circuit, warnings)
                    # Estimate tau from R and C_eff
                    tau = dominant['R'] * C_eff

                C_specific = C_eff / area_cm2
                C_specific_brug = (C_eff_brug / area_cm2
                                   if C_eff_brug is not None else None)

                warnings.extend(_window_notes(dominant, frequencies))
                if C_eff_brug is not None:
                    # Ratio is exactly (1 + R_ct/R_s)^((1-n)/n) and is always
                    # >= 1 for n <= 1; a large value means the pair no longer
                    # brackets C_eff, it just reflects how well R_s is known.
                    divergence = C_eff / C_eff_brug
                    if divergence > BRUG_HM_DIVERGENCE_MAX:
                        warnings.append(
                            f"Hsu-Mansfeld and Brug C_eff differ by {divergence:.0f}x "
                            f"(> {BRUG_HM_DIVERGENCE_MAX:.0f}x): the two CPE models do "
                            "not bracket a single C_eff here. The ratio is "
                            "(1 + R_ct/R_s)^((1-n)/n), so it is driven by the large "
                            "R_ct/R_s and by n far from 1 - treat both values, and "
                            "any thickness or permittivity derived from them, as "
                            "order-of-magnitude estimates only.")

                return {
                    'C_eff': C_eff,
                    'C_specific': C_specific,
                    'C_eff_brug': C_eff_brug,
                    'C_specific_brug': C_specific_brug,
                    'element_type': dominant['type'],
                    'element_R': dominant['R'],
                    'element_tau': tau,
                    # A copy: tau and delta_nm must not leak into `candidates`
                    'element_params': dict(dominant, tau=tau),
                    'mode': 'circuit',
                    'candidates': elements,
                    'selection_reason': selection_reason,
                    'n_hf_points': None,
                }

    # === Mode 2: Fallback - high-frequency estimate ===
    estimate = _hf_capacitance_estimate(frequencies, Z, warnings)
    if estimate is None:
        return None
    C_estimate, n_hf_points = estimate

    C_specific = C_estimate / area_cm2

    return {
        'C_eff': C_estimate,
        'C_specific': C_specific,
        'C_eff_brug': None,
        'C_specific_brug': None,
        'element_type': 'estimate',
        'element_R': None,
        'element_tau': None,
        'element_params': {},
        'mode': 'hf_estimate',
        'candidates': [],
        'selection_reason': '',
        'n_hf_points': n_hf_points,
    }


def _validate_inputs(frequencies, Z, **positive: Optional[float]):
    """Reject data and scalars that would give a crash or a signed thickness.

    Each keyword must be a finite number > 0 (area, epsilon_r, thickness):
    a zero divides by zero, a negative one flips the sign of the result.
    None means an optional input that was not given and is skipped.
    Returns the data as arrays.
    """
    for name, value in positive.items():
        if value is not None and not (np.isfinite(value) and value > 0):
            raise ValueError(f"{name} must be a finite number > 0, got {value}")
    frequencies = np.asarray(frequencies, dtype=float)
    Z = np.asarray(Z, dtype=complex)
    if frequencies.shape != Z.shape:
        raise ValueError(f"frequencies and Z differ in shape: "
                         f"{frequencies.shape} vs {Z.shape}")
    if frequencies.size == 0:
        raise ValueError("Empty data: no frequency points")
    if not np.all(np.isfinite(frequencies) & np.isfinite(Z) & (frequencies > 0)):
        raise ValueError("Data contain NaN/Inf or a frequency <= 0")
    return frequencies, Z


def analyze_oxide_layer(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    epsilon_r: float = DEFAULT_EPSILON_R,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None,
    rho_delta_ohm_cm: Optional[float] = None
) -> Optional[OxideAnalysisResult]:
    """
    Estimate oxide layer thickness from dominant capacitive element.

    Collects every capacitive element in the circuit (C, Q, K, CC, DQ, YG), keeps
    those that behave as a dielectric (admittance ~ omega^n with n near 1),
    and reports the dominant one: the largest parallel resistance - the
    compact barrier - whatever its type, with a blocking element (no DC path
    beside it) counting as infinite. The parallel resistance is the DC
    resistance of the sibling branches, so R_ct in Q | (R_ct - W) is found.
    A parallel resistance is required only to convert a Q. The thickness
    follows from the selected element's capacitance.

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
    rho_delta_ohm_cm : float, optional
        Film resistivity at the electrolyte interface [Ω·cm] for the
        power-law model; the spectrum does not determine it. If None, no
        power-law value is computed.

    Returns
    -------
    result : OxideAnalysisResult or None
        Analysis result with capacitance and thickness, or None if failed.

    Raises
    ------
    ValueError
        If epsilon_r, area_cm2 or rho_delta_ohm_cm is not a finite number > 0, or the data
        are empty, differ in shape, or contain NaN/Inf or f <= 0.

    Notes
    -----
    Thickness formula (parallel plate capacitor model):
        d = ε₀ × εᵣ / C_specific

    A Q is converted by Hsu-Mansfeld, C_eff = (R × Q)^(1/n) / R, as the
    primary value. Two comparison values follow when their inputs exist:
    Brug (1984), with a series resistance, and the power-law model
    (Hirschorn-Orazem 2010), with rho_delta_ohm_cm. Their spread is the
    model uncertainty; see doc/OXIDE_ANALYSIS_GUIDE.md.

    Examples
    --------
    >>> result, Z_fit = fit_equivalent_circuit(freq, Z, circuit)
    >>> oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
    >>> print(f"Thickness: {oxide.thickness_nm:.1f} nm")
    """
    frequencies, Z = _validate_inputs(frequencies, Z, epsilon_r=epsilon_r, area_cm2=area_cm2,
                                      rho_delta_ohm_cm=rho_delta_ohm_cm)
    warnings: List[str] = []
    extracted = _extract_capacitance(frequencies, Z, area_cm2, fit_result,
                                     warnings)
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

    d_pl_cm = _power_law(extracted['element_params'], area_cm2, rho_delta_ohm_cm,
                         frequencies.max(), warnings, eps_r=epsilon_r)

    _add_penetration_depth(extracted['element_params'], d_nm)

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
        thickness_brug_nm=d_brug_nm,
        thickness_pl_nm=d_pl_cm * 1e7 if d_pl_cm is not None else None,
        rho_delta_ohm_cm=rho_delta_ohm_cm,
        mode=extracted['mode'],
        candidates=extracted['candidates'],
        selection_reason=extracted['selection_reason'],
        n_hf_points=extracted['n_hf_points'],
        epsilon_r=epsilon_r,
        area_cm2=area_cm2,
        warnings=warnings
    )


def estimate_permittivity(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    thickness_nm: float,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None,
    rho_delta_ohm_cm: Optional[float] = None
) -> Optional[OxideAnalysisResult]:
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
    rho_delta_ohm_cm : float, optional
        Film resistivity at the electrolyte interface [Ω·cm] for the
        power-law model; the spectrum does not determine it. If None, no
        power-law value is computed.

    Returns
    -------
    result : OxideAnalysisResult or None
        Analysis result with capacitance and permittivity, or None if
        failed. The estimate is in `permittivity`; `thickness_nm` holds
        the thickness that was passed in, since here it is the input.

    Raises
    ------
    ValueError
        If thickness_nm, area_cm2 or rho_delta_ohm_cm is not a finite number > 0, or the data
        are empty, differ in shape, or contain NaN/Inf or f <= 0.

    Notes
    -----
    Formula (from parallel plate capacitor model):
        ε_r = d × C_specific / ε₀

    A Q is converted as in analyze_oxide_layer(), with the Brug and
    power-law comparison values in permittivity_brug / permittivity_pl.

    Examples
    --------
    >>> result, Z_fit = fit_equivalent_circuit(freq, Z, circuit)
    >>> oxide = estimate_permittivity(freq, Z, thickness_nm=20, fit_result=result)
    >>> print(f"Permittivity: {oxide.permittivity:.1f}")
    """
    frequencies, Z = _validate_inputs(frequencies, Z, thickness_nm=thickness_nm,
                                      area_cm2=area_cm2, rho_delta_ohm_cm=rho_delta_ohm_cm)
    warnings: List[str] = []

    # Get capacitance using the same element-selection logic as
    # analyze_oxide_layer (no thickness is computed here)
    extracted = _extract_capacitance(frequencies, Z, area_cm2, fit_result,
                                     warnings)

    if extracted is None:
        return None

    # Calculate permittivity from thickness and capacitance
    # d = ε₀ × εᵣ / C_specific  =>  εᵣ = d × C_specific / ε₀
    d_cm = thickness_nm * 1e-7  # nm -> cm
    C_specific = extracted['C_specific']
    epsilon_r = d_cm * C_specific / EPSILON_0

    # Brug (2D) comparison permittivity, when available
    C_specific_brug = extracted['C_specific_brug']
    eps_r_brug = None
    if C_specific_brug is not None:
        eps_r_brug = d_cm * C_specific_brug / EPSILON_0

    eps_r_pl = _power_law(extracted['element_params'], area_cm2, rho_delta_ohm_cm,
                          frequencies.max(), warnings, d_cm=d_cm)

    _add_penetration_depth(extracted['element_params'], thickness_nm)

    return OxideAnalysisResult(
        capacitance=extracted['C_eff'],
        capacitance_specific=C_specific,
        thickness_nm=thickness_nm,
        element_type=extracted['element_type'],
        element_R=extracted['element_R'],
        element_tau=extracted['element_tau'],
        element_params=extracted['element_params'],
        capacitance_brug=extracted['C_eff_brug'],
        capacitance_specific_brug=C_specific_brug,
        permittivity=epsilon_r,
        permittivity_brug=eps_r_brug,
        permittivity_pl=eps_r_pl,
        rho_delta_ohm_cm=rho_delta_ohm_cm,
        mode=extracted['mode'],
        candidates=extracted['candidates'],
        selection_reason=extracted['selection_reason'],
        n_hf_points=extracted['n_hf_points'],
        area_cm2=area_cm2,
        warnings=warnings
    )


__all__ = ['analyze_oxide_layer', 'estimate_permittivity', 'OxideAnalysisResult']
