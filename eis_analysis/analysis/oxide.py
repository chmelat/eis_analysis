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

from ..fitting.bounds import PARAMETER_BOUNDS, classify_bound_status
from ..fitting.circuit import FitResult
from ..fitting.circuit_elements import R, C, G, Q, K, CC
from ..fitting.circuit_builder import Series, Parallel
from .config import (
    EPSILON_0,
    DEFAULT_EPSILON_R,
    CPE_N_RELIABLE_MIN,
    HF_ESTIMATE_DECADE_FACTOR,
    HF_C_SPREAD_MAX_RATIO,
    BRUG_RS_MIN_OHM,
    BRUG_HM_DIVERGENCE_MAX,
    CC_WINDOW_EDGE_MARGIN_DECADES,
)

logger = logging.getLogger(__name__)


@dataclass
class OxideAnalysisResult:
    """Result of oxide layer analysis."""
    capacitance: float          # Effective capacitance [F]
    capacitance_specific: float # Specific capacitance [F/cm²]
    thickness_nm: float         # Oxide thickness [nm]
    element_type: str           # 'C', 'K', 'Q', 'CC', or 'estimate' (HF fallback)
    element_R: Optional[float]  # Associated resistance [Ω] (None for CC)
    element_tau: Optional[float] # Time constant [s]
    element_params: Dict[str, float]  # All element parameters
    # Brug (2D) comparison values - set only for a dominant Q element
    # when a series resistance is present in the circuit
    capacitance_brug: Optional[float] = None          # Brug C_eff [F]
    capacitance_specific_brug: Optional[float] = None # Brug C_eff/area [F/cm²]
    thickness_brug_nm: Optional[float] = None         # Thickness from Brug C [nm]
    # Inverse mode - set only by estimate_permittivity(), where the
    # thickness is the input and the permittivity the derived quantity
    permittivity: Optional[float] = None       # ε_r from known thickness
    permittivity_brug: Optional[float] = None  # ε_r from Brug C_eff


def _cc_capacitance_regime(tau: float, frequencies: NDArray[np.float64]) -> str:
    """
    Which limit of C*(omega) the measured frequency window actually determines.

    A Cole-Cole element disperses around f_char = 1/(2*pi*tau):

        omega*tau << 1  ->  C*(omega) -> C_s = C_inf + dC   (static limit)
        omega*tau >> 1  ->  C*(omega) -> C_inf              (high-frequency limit)

    Returns
    -------
    'high_frequency'
        f_char lies below the lowest measured frequency: every measured point
        sits at omega*tau >> 1, so the data constrain only C_inf. dC is then
        an extrapolation to DC through a region with no measurements and C_s
        must not be reported as if it had been measured.
    'static'
        f_char lies inside the window (the relaxation is traced, C_s is the
        exact omega -> 0 limit) or above it (the whole window sits at
        omega*tau << 1, so C_s is measured - though only as a sum, the
        C_inf/dC split being unidentified). Both cases report C_s.

    A non-positive tau or an empty frequency array leaves no window test to
    make; the static limit is the unconditional pre-0.25.2 behavior.
    """
    if tau <= 0 or frequencies.size == 0:
        return 'static'
    f_char = 1.0 / (2.0 * np.pi * tau)
    return 'high_frequency' if f_char < float(np.min(frequencies)) else 'static'


def _log_cc_capacitance_choice(
    cc: Dict[str, Any],
    frequencies: NDArray[np.float64]
) -> None:
    """
    Explain which Cole-Cole capacitance was reported, and why.

    Two independent checks, both worth logging:

    1. Where f_char = 1/(2*pi*tau) sits relative to the measured window.
       This is what selects the reported value in `_cc_capacitance_regime`;
       outside the window (either side) at least one of C_inf / dC is not
       determined by the data, and anything derived from it - thickness,
       permittivity - inherits that.
    2. Whether tau itself landed on a fitting bound. A parameter on its bound
       always means the data did not determine it, so the test is made
       explicitly rather than inferred from f_char. `classify_bound_status`
       is the project-wide definition of "at a bound", and the bounds cannot
       be overridden by a caller (`generate_simple_bounds` derives them from
       the parameter labels), so PARAMETER_BOUNDS is authoritative here.
    """
    tau = cc['tau']
    f_char = 1.0 / (2.0 * np.pi * tau) if tau > 0 else None
    f_min = float(np.min(frequencies)) if frequencies.size else 0.0
    f_max = float(np.max(frequencies)) if frequencies.size else 0.0

    if cc['C_regime'] == 'high_frequency' and f_char is not None:
        logger.warning(
            f"Cole-Cole relaxation lies BELOW the measured window: "
            f"f_char = 1/(2*pi*tau) = {f_char:.2e} Hz vs f_min = {f_min:.2e} Hz "
            f"({np.log10(f_min / f_char):.1f} decades below it). Every measured "
            f"point sits at omega*tau >> 1, where C*(omega) -> C_inf, so "
            f"dC = {cc['dC']:.3e} F is an extrapolation to DC through a region "
            f"with no data. Reporting C_inf = {cc['C_inf']:.3e} F instead of "
            f"C_s = {cc['C_static']:.3e} F - the thickness/permittivity below is "
            f"the high-frequency value. Extend the sweep to lower frequencies "
            f"to determine dC.")
    elif f_char is not None and frequencies.size > 0 and f_min > 0 and f_char > f_max:
        logger.warning(
            f"Cole-Cole relaxation lies ABOVE the measured window: "
            f"f_char = {f_char:.2e} Hz vs f_max = {f_max:.2e} Hz. The whole "
            f"window sits at omega*tau << 1, so the reported "
            f"C_s = {cc['C_static']:.3e} F is what the data determine - but only "
            f"as a sum: the split into C_inf = {cc['C_inf']:.3e} F and "
            f"dC = {cc['dC']:.3e} F is not identified. Check their confidence "
            f"intervals before reading either value on its own.")
    elif (f_char is not None and frequencies.size > 0 and f_min > 0
          and min(np.log10(f_char / f_min),
                  np.log10(f_max / f_char)) < CC_WINDOW_EDGE_MARGIN_DECADES):
        logger.warning(
            f"Cole-Cole relaxation sits within "
            f"{CC_WINDOW_EDGE_MARGIN_DECADES:g} decade(s) of the edge of the "
            f"measured window (f_char = {f_char:.2e} Hz, window "
            f"{f_min:.2e} .. {f_max:.2e} Hz). Only the tail of the dispersion is "
            f"traced, so the C_inf/dC split rests on the last few points of the "
            f"sweep; C_s = {cc['C_static']:.3e} F is reported but is only "
            f"marginally determined.")

    if not cc['tau_fixed']:
        tau_lo, tau_hi = PARAMETER_BOUNDS['τ_CC']
        status = classify_bound_status(tau, tau_lo, tau_hi)
        if status:
            logger.warning(
                f"Cole-Cole tau = {tau:.2e} s sits at its {status} fitting bound "
                f"({tau_lo:.0e} .. {tau_hi:.0e} s): the data did not determine "
                f"it, so the relaxation strength dC and everything derived from "
                f"it are unconstrained. Either extend the frequency range or fix "
                f"tau to an independently known value.")


def _find_capacitive_elements(
    circuit,
    frequencies: NDArray[np.float64]
) -> List[Dict[str, Any]]:
    """
    Find every capacitive element in the circuit: C, Q, K and CC.

    A capacitive element is a candidate on its own account, whether or not it
    shares a parallel combination with a resistance. Requiring a parallel R -
    a Voigt element - used to hide a perfectly well fitted capacitance: in
    `L - R0 - (Q|C)` the `C` has no resistance beside it, so the whole
    analysis fell through to the high-frequency spectral estimate even though
    `C` was a fitted parameter with a 0.7 % confidence interval.

    A parallel resistance, where one exists, is still recorded: it is what the
    Hsu-Mansfeld/Brug conversion of a `Q` needs, what `tau = R*C` needs, and
    what the largest-R barrier heuristic compares. `R` = None simply means the
    element has no DC path beside it.

    Returns list of dicts. Common keys: 'type', 'R' (may be None), 'tau' (may
    be None), and 'n' - the exponent of the element's admittance, Y ~ omega^n,
    which is how the dielectric element is identified downstream. A CC entry
    additionally carries the capacitance the measured window determines (see
    `_cc_capacitance_regime`), which is why `frequencies` is needed here.
    """
    results = []

    def parallel_resistance(node) -> Optional[float]:
        """Resistance of the R/G elements that are direct children of a Parallel.

        G is a resistance parametrized as its own reciprocal, so it defines
        the parallel resistance exactly as R does. G = 0 is an open branch:
        R = 1/G is infinite, which is reported as None (no DC path), the
        same convention the blocking-dielectric branches use.
        """
        r_elements = [e for e in node.elements if isinstance(e, (R, G))]
        if not r_elements:
            return None
        if len(r_elements) > 1:
            logger.warning("Multiple R/G elements in one parallel "
                           "combination - using the last one")
        last = r_elements[-1]
        if isinstance(last, G):
            return 1 / last.G if last.G > 0 else None
        return last.params[0]

    def traverse(node, R_parallel: Optional[float]) -> None:
        if isinstance(node, Parallel):
            # An inner Parallel with no R of its own stays under the enclosing
            # one's resistance: R | (C | Q) is just R | C | Q rewritten.
            own_R = parallel_resistance(node)
            R_here = own_R if own_R is not None else R_parallel
            for elem in node.elements:
                traverse(elem, R_here)

        elif isinstance(node, Series):
            # A series boundary ends the scope of any enclosing parallel R:
            # in R | (C - R2) the capacitor is in series with R2, not parallel
            # to R.
            for elem in node.elements:
                traverse(elem, None)

        elif isinstance(node, C):
            C_val = node.params[0]
            results.append({
                'type': 'C',
                'R': R_parallel,
                'C': C_val,
                'n': 1.0,           # an ideal capacitor, by construction
                'tau': R_parallel * C_val if R_parallel else None,
            })

        elif isinstance(node, Q):
            results.append({
                'type': 'Q',
                'R': R_parallel,
                'Q': node.params[0],
                'n': node.params[1],
                'tau': None,        # needs C_eff first; computed later
            })

        elif isinstance(node, K):
            # K element directly provides R and tau
            R_val = node.params[0]
            tau_val = node.params[1]
            if R_val <= 0:
                # C = tau/R is undefined; a non-positive R would be dropped
                # by the dominant-element filter anyway
                logger.warning(f"K element with non-positive R = {R_val:g} Ω - skipping")
                return
            results.append({
                'type': 'K',
                'R': R_val,
                'C': tau_val / R_val,
                'n': 1.0,           # K is a Voigt R||C reparameterised
                'tau': tau_val,
            })

        elif isinstance(node, CC):
            C_inf_val, dC_val = node.params[0], node.params[1]
            tau_val, alpha_val = node.params[2], node.params[3]
            # The static (fully relaxed) capacitance is the one that pairs
            # with a static permittivity such as eps_r = 22 for ZrO2 - but
            # only when the data reach omega*tau << 1. With the relaxation
            # below the measured window the fit determines C_inf alone, and
            # C_s = C_inf + dC is an extrapolation to DC (see
            # _cc_capacitance_regime and _log_cc_capacitance_choice).
            regime = _cc_capacitance_regime(tau_val, frequencies)
            results.append({
                'type': 'CC',
                'R': None,          # blocking dielectric - no DC path
                'C': C_inf_val if regime == 'high_frequency' else C_inf_val + dC_val,
                'n': 1.0,           # Y = j*omega*C*(omega) in both limits
                'tau': tau_val,
                'C_inf': C_inf_val,
                'dC': dC_val,
                'C_static': C_inf_val + dC_val,
                'C_regime': regime,
                # A tau pinned by the user (CC(..., tau="1e4")) is a choice,
                # not an undetermined fit parameter - it must not raise the
                # "parameter sits at its bound" warning.
                'tau_fixed': bool(node.fixed_params[2]),
                'alpha': alpha_val,
            })

    traverse(circuit, None)
    return results


def _element_size(element: Dict[str, Any]) -> float:
    """Capacitive magnitude used to rank elements the R heuristic cannot separate.

    The fitted capacitance for C, K and CC; the CPE coefficient for Q, which
    is only ever compared against other Q elements (a tier holds one kind).
    """
    return float(element['C'] if 'C' in element else element['Q'])


def _select_dielectric_element(
    candidates: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """
    Pick the element that carries the dielectric response.

    Caller has already applied the physical criterion (admittance ~ omega^n
    with n near 1); this only resolves which of several qualifying elements to
    report, in order of how directly each one's capacitance is determined:

    1. `CC` - the general dielectric relaxation model. A plain `C` is its
       degenerate case (dC = 0), so if both somehow appear the general model
       wins; preferring the simpler element over the more general one would be
       backwards.
    2. `C` and `K` - the capacitance is a fitted parameter, exact.
    3. `Q` - a near-ideal CPE, whose capacitance still needs the
       Hsu-Mansfeld/Brug model on top of the fit.

    Within a tier the largest-R element is taken to be the compact barrier -
    the long-standing heuristic, which only distinguishes an oxide barrier
    from a charge-transfer process. Elements with no parallel resistance have
    no R to compare and are ranked by capacitance instead.
    """
    cc = [e for e in candidates if e['type'] == 'CC']
    if cc:
        if len(cc) > 1:
            logger.warning(
                f"{len(cc)} Cole-Cole elements in circuit - using "
                "the one with the largest static capacitance. With several "
                "dielectric relaxations the layer assignment is yours to make.")
        dominant = max(cc, key=lambda e: e['C'])
        logger.info(f"Dominant element: CC with C = {dominant['C']:.3e} F")
        logger.info("Selected because the circuit models the dielectric "
                    "relaxation explicitly; a plain C is its degenerate case "
                    "(ΔC = 0), so the general model wins")
        return dominant

    exact = [e for e in candidates if e['type'] in ('C', 'K')]
    tier = exact if exact else candidates
    with_R = [e for e in tier if e['R'] is not None and e['R'] > 0]

    if with_R:
        # Rank by R, then by size: elements sharing one parallel resistance sit
        # in the same combination, where the heuristic has nothing to say and
        # the larger capacitance is the dominant contribution
        dominant = max(with_R, key=lambda e: (e['R'], _element_size(e)))
        tied = [e for e in with_R if e['R'] == dominant['R']]
        if len(tied) > 1:
            logger.warning(
                f"{len(tied)} capacitive elements share one parallel resistance "
                f"(R = {dominant['R']:.1f} Ω) - using the largest, "
                f"{_element_size(dominant):.3e}. They are one parallel "
                "combination, so their individual values are not separately "
                "identifiable from the spectrum.")
        if len(tier) > len(with_R):
            logger.warning(
                f"{len(tier) - len(with_R)} capacitive element(s) have no "
                "parallel resistance and cannot be ranked by the largest-R "
                "heuristic - they were not considered for the dominant "
                "element. Name the layer explicitly if one of them is it.")
        logger.info(f"Dominant element: {dominant['type']} with "
                    f"R = {dominant['R']:.1f} Ω")
        logger.info("Selection assumes the largest-R element is the compact "
                    "oxide barrier (verify: a charge-transfer process can "
                    "also have the largest R)")
        return dominant

    # No resistance anywhere in the tier: rank by capacitance instead
    dominant = max(tier, key=_element_size)
    if len(tier) > 1:
        logger.warning(
            f"{len(tier)} capacitive elements with no parallel resistance - "
            "using the largest. Their individual values are not separately "
            "identifiable from the spectrum, so check the fit before relying "
            "on the split.")
    logger.info(f"Dominant element: {dominant['type']} with "
                f"C = {_element_size(dominant):.3e} F "
                f"(no parallel resistance to rank by)")
    return dominant


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
        Keys: 'C_eff' [F], 'C_specific' [F/cm²], 'C_eff_brug' [F],
        'C_specific_brug' [F/cm²], 'element_type', 'element_R',
        'element_tau', 'element_params'.
        The two Brug (2D) keys are always present but None unless the
        dominant element is a Q and the circuit has a series resistance
        of at least BRUG_RS_MIN_OHM (below that the fit has not identified
        R_s and Brug's R_s^((1-n)/n) scaling makes the value meaningless);
        in fallback mode 'element_R' and 'element_tau' are None too,
        and 'element_R' is None for a CC (a blocking dielectric has
        no parallel resistance).
        None if capacitance could not be extracted.
    """
    # === Mode 1: From fitted circuit (preferred) ===
    if fit_result is not None:
        circuit = fit_result.circuit

        # Every C, Q, K and CC in the circuit, parallel resistance or not
        elements = _find_capacitive_elements(circuit, frequencies)

        if not elements:
            logger.warning("No capacitive element (C, Q, K, CC) found in circuit")
            logger.warning("Falling back to high-frequency estimate...")
            fit_result = None
        else:
            # List all candidates so the dominant-element choice can be
            # verified (largest R may also be a charge-transfer process)
            logger.info(f"Found {len(elements)} capacitive element(s):")
            for i, e in enumerate(elements, 1):
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
                else:
                    tau_str = (f", tau = {e['tau']:.2e} s" if e['tau'] is not None
                               else "")
                    logger.info(f"  [{i}] {e['type']}: {R_str}, "
                                f"C = {e['C']:.3e} F{tau_str}")

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
                    logger.warning(
                        f"CPE with n = {e['n']:.3f} has no parallel resistance - "
                        "neither Hsu-Mansfeld nor Brug can convert it to a "
                        "capacitance, so it is not a candidate.")
                else:
                    usable.append(e)

            dielectric = [e for e in usable if e['n'] >= CPE_N_RELIABLE_MIN]
            non_ideal = [e for e in usable if e['n'] < CPE_N_RELIABLE_MIN]

            if dielectric:
                dominant = _select_dielectric_element(dielectric)
            elif non_ideal:
                # Nothing in the circuit behaves as a dielectric. Better than
                # the spectral fallback - these are at least fitted parameters
                # - but the result is not a dielectric capacitance.
                logger.warning(
                    f"No dielectric element in circuit: the only capacitive "
                    f"element(s) are CPEs with n < {CPE_N_RELIABLE_MIN} "
                    f"(largest n = {max(e['n'] for e in non_ideal):.3f}). Such a "
                    "CPE describes transport or a distribution of resistivity, "
                    "not a dielectric, so the capacitance below - and any "
                    "thickness or permittivity from it - has no dielectric "
                    "meaning. Reported only because it is still a fitted "
                    "parameter, unlike the high-frequency estimate.")
                dominant = _select_dielectric_element(non_ideal)
            else:
                logger.warning("No convertible capacitive element found")
                logger.warning("Falling back to high-frequency estimate...")
                fit_result = None
                dominant = None

            if dominant is not None:

                # Get capacitance
                C_eff_brug = None
                if dominant['type'] in ('C', 'K', 'CC'):
                    # For CC this is exact, not an effective capacitance: both
                    # C_s = C_inf + dC (the omega -> 0 limit of C*(omega), for
                    # any alpha) and C_inf (the omega -> inf limit) are model
                    # limits, not fits. No Hsu-Mansfeld / Brug choice arises;
                    # which limit the data support was decided in
                    # _cc_capacitance_regime.
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
                    if R_s is None:
                        logger.info("No series R element in circuit - "
                                    "Brug (2D) estimate not available")
                    elif R_s < BRUG_RS_MIN_OHM:
                        # A CPE with n < 1 mimics a series resistance at high
                        # frequency, so R_s is often unidentifiable and the fit
                        # drives it to the optimizer floor. Brug's
                        # C ~ R_s^((1-n)/n) would then be arbitrarily small.
                        logger.warning(
                            f"Series resistance R_s = {R_s:.3e} Ohm < "
                            f"{BRUG_RS_MIN_OHM:.3g} Ohm: the fit did not identify it "
                            "(a CPE with n < 1 mimics a series resistance at high "
                            "frequency, so R_s collapses to its lower bound). "
                            "Brug (2D) estimate suppressed - it would scale as "
                            "R_s^((1-n)/n) and be meaningless. Check the fitted "
                            "R_s against Re(Z) at the highest measured frequency.")
                    else:
                        C_eff_brug = _estimate_cpe_capacitance_brug(
                            dominant['Q'], dominant['n'], dominant['R'], R_s
                        )

                C_specific = C_eff / area_cm2
                C_specific_brug = (C_eff_brug / area_cm2
                                   if C_eff_brug is not None else None)

                # Log element info (thickness/permittivity logged by caller)
                logger.info("")
                logger.info("Results:")
                logger.info(f"  Element type:       {dominant['type']}")
                if dominant['R'] is not None:
                    logger.info(f"  Resistance:         {dominant['R']:.1f} Ω")
                elif dominant['type'] == 'CC':
                    logger.info("  Resistance:         n/a (blocking dielectric)")
                else:
                    logger.info("  Resistance:         n/a (no parallel R "
                                "in the circuit)")
                if dominant['type'] == 'CC':
                    logger.info(f"  Broadening alpha:   {dominant['alpha']:.3f}")
                    logger.info(f"  C_inf / ΔC:         {dominant['C_inf']:.3e} F / "
                                f"{dominant['dC']:.3e} F")
                cc_suffix = ""
                if dominant['type'] == 'CC':
                    cc_suffix = ("  (C_inf, high-frequency limit)"
                                 if dominant['C_regime'] == 'high_frequency'
                                 else "  (static, C_inf + ΔC)")
                logger.info(f"  Capacitance:        {C_eff:.3e} F{cc_suffix}")
                if dominant['type'] == 'CC':
                    _log_cc_capacitance_choice(dominant, frequencies)
                if C_eff_brug is not None:
                    logger.info(f"  C (Brug, 2D):       {C_eff_brug:.3e} F "
                                f"(comparison; primary value is Hsu-Mansfeld, 3D)")
                    # Ratio is exactly (1 + R_ct/R_s)^((1-n)/n) and is always
                    # >= 1 for n <= 1; a large value means the pair no longer
                    # brackets C_eff, it just reflects how well R_s is known.
                    divergence = C_eff / C_eff_brug
                    if divergence > BRUG_HM_DIVERGENCE_MAX:
                        logger.warning(
                            f"Hsu-Mansfeld and Brug C_eff differ by {divergence:.0f}x "
                            f"(> {BRUG_HM_DIVERGENCE_MAX:.0f}x): the two CPE models do "
                            "not bracket a single C_eff here. The ratio is "
                            "(1 + R_ct/R_s)^((1-n)/n), so it is driven by the large "
                            "R_ct/R_s and by n far from 1 - treat both values, and "
                            "any thickness or permittivity derived from them, as "
                            "order-of-magnitude estimates only.")
                logger.info(f"  Specific cap.:      {C_specific * 1e6:.2f} µF/cm²")
                if tau is not None and tau > 0:
                    logger.info(f"  Time constant:      {tau:.3e} s")
                    logger.info(f"  Char. frequency:    {1/(2*np.pi*tau):.2e} Hz")
                else:
                    # A capacitance with no parallel resistance has no RC time
                    # constant - there is nothing to discharge through
                    logger.info("  Time constant:      n/a (no parallel R)")

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
    # Stated as plainly as a parameter sitting on its bound: the number below
    # is a spectral guess, not a fitted quantity, and nothing downstream
    # distinguishes the two once they are printed side by side.
    logger.warning(
        "NOT FROM THE FIT: the capacitance below is estimated directly from "
        "the spectrum (median of C = -1/(omega*Z'') over the top frequency "
        "decade), because the circuit offered no capacitive element to read it "
        "from. It carries no confidence interval and the thickness or "
        "permittivity derived from it is an order-of-magnitude figure - treat "
        "it as such even when it happens to land near the expected value.")
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
    epsilon_r: float = DEFAULT_EPSILON_R,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[OxideAnalysisResult]:
    """
    Estimate oxide layer thickness from dominant capacitive element.

    Collects every capacitive element in the circuit (C, Q, K, CC), keeps
    those that behave as a dielectric (admittance ~ omega^n with n near 1),
    and reports the dominant one: the general dielectric model first (CC),
    then the exactly fitted capacitances (C, K), then a near-ideal CPE (Q),
    and within a tier the largest parallel resistance - the compact barrier.
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

    Returns
    -------
    result : OxideAnalysisResult or None
        Analysis result with capacitance and permittivity, or None if
        failed. The estimate is in `permittivity`; `thickness_nm` holds
        the thickness that was passed in, since here it is the input.

    Notes
    -----
    Formula (from parallel plate capacitor model):
        ε_r = d × C_specific / ε₀

    For Q elements the capacitance conversion is the same as in
    analyze_oxide_layer(): Hsu-Mansfeld (normal/3D distribution) as the
    primary value, and the Brug (1984) formula (surface/2D distribution)
    reported in permittivity_brug for comparison when the circuit also
    contains a series resistance. The spread between the two brackets
    the model uncertainty. See doc/OXIDE_ANALYSIS_GUIDE.md.

    Examples
    --------
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
    >>> oxide = estimate_permittivity(freq, Z, thickness_nm=20, fit_result=result)
    >>> print(f"Permittivity: {oxide.permittivity:.1f}")
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

    # Brug (2D) comparison permittivity, when available
    C_specific_brug = extracted['C_specific_brug']
    eps_r_brug = None
    if C_specific_brug is not None:
        eps_r_brug = d_cm * C_specific_brug / EPSILON_0

    logger.info(f"  Known thickness:    {thickness_nm:.1f} nm")
    logger.info(f"  Permittivity ε_r:   {epsilon_r:.1f}")
    if eps_r_brug is not None:
        # .3g, not .1f: when n is far from 1 the two CPE models diverge by
        # orders of magnitude and a fixed-point format would print "0.0"
        logger.info(f"  ε_r (Brug):         {eps_r_brug:.3g} "
                    f"(2D model, for comparison)")
    logger.info(f"  (area={area_cm2} cm²)")
    logger.info("=" * 50)

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
        permittivity_brug=eps_r_brug
    )


__all__ = ['analyze_oxide_layer', 'estimate_permittivity', 'OxideAnalysisResult']
