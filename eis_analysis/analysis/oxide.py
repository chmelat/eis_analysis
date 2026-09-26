"""
Oxide layer thickness estimation from EIS data.

Simplified implementation that finds the dominant Voigt or Q element
and estimates oxide thickness from its capacitance.
"""

import numpy as np
import logging
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass, field
from numpy.typing import NDArray

from ..fitting.bounds import PARAMETER_BOUNDS, classify_bound_status
from ..fitting.circuit import FitResult
from ..fitting.circuit_elements import R, C, G, Q, K, CC, DQ, YG
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

    # How the capacitance was obtained
    mode: str = 'circuit'       # 'circuit' (from the fit) or 'hf_estimate'
    candidates: List[Dict[str, Any]] = field(default_factory=list)
    selection_reason: str = ''  # why the dominant element was picked
    n_hf_points: Optional[int] = None  # points behind the HF median estimate
    epsilon_r: Optional[float] = None  # value assumed by analyze_oxide_layer
    area_cm2: float = 1.0
    warnings: List[str] = field(default_factory=list)


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


def _plateau_notes(
    label: str,
    tau_cap: float,
    frequencies: NDArray[np.float64],
    capacitance: str,
    leans_on: str
) -> List[str]:
    """Say whether an element's capacitive plateau was measured or extrapolated.

    Every element whose capacitance is the omega -> inf limit of a fitted
    model - DQ above 1/tau_min, YG above 1/(2*pi*tau) - is exact only where
    the window reaches that plateau. If f_cap lies above the window the
    capacitance, and every thickness or permittivity from it, is an
    extrapolation; if it sits within the last decade, only a sliver was
    measured and the value leans on whichever parameter shapes the approach.
    Same two tests, and the same edge margin, as the Cole-Cole notes.

    Attached to the dominant element only: one that lost the selection
    contributed nothing to the reported thickness, and warning about its
    capacitance would point the reader at the wrong element.

    Parameters
    ----------
    label : str
        Element symbol, opening each message.
    tau_cap : float
        Time constant of the plateau's corner [s]; the element is capacitive
        above 1/(2*pi*tau_cap). Taken rather than the frequency so that a
        non-positive tau is rejected before the division, not after it.
    frequencies : ndarray
        The measured window.
    capacitance : str
        How the element names its capacitance, e.g. 'C_eff' or 'C'.
    leans_on : str
        What the edge case leaves the value resting on, as a sentence tail.
    """
    if tau_cap <= 0 or frequencies.size == 0:
        return []

    f_cap = 1.0 / (2.0 * np.pi * tau_cap)
    f_max = float(np.max(frequencies))

    if f_cap > f_max:
        return [f"{label}: the capacitive plateau starts at {f_cap:.3g} Hz, above "
                f"the highest measured frequency {f_max:.3g} Hz - {capacitance} "
                "is an extrapolation, and so is any thickness derived from it. "
                "Extend the sweep upwards to measure it."]
    if np.log10(f_max / f_cap) < CC_WINDOW_EDGE_MARGIN_DECADES:
        return [f"{label}: the capacitive plateau starts at {f_cap:.3g} Hz, less "
                f"than {CC_WINDOW_EDGE_MARGIN_DECADES:.0f} decade below the "
                f"highest measured frequency {f_max:.3g} Hz - only its edge is "
                f"measured, so {capacitance} {leans_on}"]
    return []


def _yg_R_dc_note(
    yg: Dict[str, Any],
    frequencies: NDArray[np.float64]
) -> List[str]:
    """State what YG's R_dc is, for a window that does not reach it.

    Unlike the plateau notes this is not a warning. f_R lies e^(1/p) below
    the capacitive corner, which for any p <~ 0.1 is dozens of decades under
    any real sweep, so the note describes what the number is - a value of the
    model - and is phrased that way on purpose.
    """
    f_R = yg['f_R']
    if frequencies.size == 0 or not (0 < f_R < float(np.min(frequencies))):
        return []

    f_min = float(np.min(frequencies))
    return [f"YG: R_dc = {yg['R_dc']:.3g} Ohm is the omega -> 0 limit of the "
            f"model, reached below {f_R:.3g} Hz - "
            f"{np.log10(f_min / f_R):.0f} decades under the lowest measured "
            f"frequency {f_min:.3g} Hz. Expected for an exponential "
            "conductivity profile, and the reason it is not used to rank "
            "elements; read it as a model value, not a measured resistance."]


def _cc_capacitance_notes(
    cc: Dict[str, Any],
    frequencies: NDArray[np.float64]
) -> List[str]:
    """
    Explain which Cole-Cole capacitance was reported, and why.

    Two independent checks, both worth reporting:

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
    notes: List[str] = []
    tau = cc['tau']
    f_char = 1.0 / (2.0 * np.pi * tau) if tau > 0 else None
    f_min = float(np.min(frequencies)) if frequencies.size else 0.0
    f_max = float(np.max(frequencies)) if frequencies.size else 0.0

    if cc['C_regime'] == 'high_frequency' and f_char is not None:
        notes.append(
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
        notes.append(
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
        notes.append(
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
            notes.append(
                f"Cole-Cole tau = {tau:.2e} s sits at its {status} fitting bound "
                f"({tau_lo:.0e} .. {tau_hi:.0e} s): the data did not determine "
                f"it, so the relaxation strength dC and everything derived from "
                f"it are unconstrained. Either extend the frequency range or fix "
                f"tau to an independently known value.")

    return notes


def _find_capacitive_elements(
    circuit,
    frequencies: NDArray[np.float64],
    warnings: List[str]
) -> List[Dict[str, Any]]:
    """
    Find every capacitive element in the circuit: C, Q, K, CC, DQ and YG.

    A capacitive element is a candidate on its own account, whether or not it
    shares a parallel combination with a resistance. Requiring a parallel R -
    a Voigt element - used to hide a perfectly well fitted capacitance: in
    `L - R0 - (Q|C)` the `C` has no resistance beside it, so the whole
    analysis fell through to the high-frequency spectral estimate even though
    `C` was a fitted parameter with a 0.7 % confidence interval.

    A parallel resistance, where one exists, is still recorded: it is what the
    Hsu-Mansfeld/Brug conversion of a `Q` needs, what `tau = R*C` needs, and
    what the largest-R barrier heuristic compares. `R` = None means the
    fitted model gives the element no DC path beside it.

    Returns list of dicts. Common keys: 'type', 'R' (may be None), 'tau' (may
    be None), and 'n' - the exponent of the element's admittance, Y ~ omega^n,
    which is how the dielectric element is identified downstream. A CC entry
    additionally carries the capacitance the measured window determines (see
    `_cc_capacitance_regime`), which is why `frequencies` is needed here.
    """
    results = []

    def traverse(node, R_parallel: Optional[float]) -> None:
        if isinstance(node, Parallel):
            for i, elem in enumerate(node.elements):
                # The resistance beside an element: its sibling branches
                # (_branch_resistance) and whatever encloses this Parallel,
                # combined in parallel - R_ct in Q | (R_ct - W), R_ox in
                # Q_ox | (R_ox - (Q_dl | R_ct)), R1 || R2 in R1 | (R2 | C).
                siblings = node.elements[:i] + node.elements[i + 1:]
                enclosing = [R_parallel] if R_parallel is not None else []
                R_here = _parallel_combination(
                    [_branch_resistance(s) for s in siblings] + enclosing)
                traverse(elem, R_here if R_here < np.inf else None)

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
                warnings.append(f"K element with non-positive R = {R_val:g} Ω - skipping")
                return
            results.append({
                'type': 'K',
                'R': R_val,
                'C': tau_val / R_val,
                'n': 1.0,           # K is a Voigt R||C reparameterised
                'tau': tau_val,
            })

        elif isinstance(node, DQ):
            # A truncated CPE brings its own DC path and its own capacitive
            # plateau: R_pol and C_eff are limits of the fitted distribution,
            # exact within the model, so no Hsu-Mansfeld / Brug estimate is
            # needed on top.
            # Whether C_eff is measured or extrapolated is decided later, on
            # the winner alone (_plateau_notes).
            R_dc = node.R_pol
            if R_parallel is not None and R_parallel > 0:
                # In R | DQ both paths reach DC. The reported resistance and
                # the largest-R barrier heuristic must see the combination;
                # R_pol alone can be four orders too large.
                R_dc = R_parallel * R_dc / (R_parallel + R_dc)
            results.append({
                'type': 'DQ',
                'R': R_dc,          # R_pol, shunted by an enclosing parallel R
                'C': node.C_eff,
                # 1.0, like CC and for the same reason: above 1/tau_min the
                # element *is* a capacitor, so C_eff is a limit of the model
                # and not a Hsu-Mansfeld conversion whose reliability decays
                # with the exponent. The power law is reported as n_power;
                # what qualifies the capacitance here is whether the plateau
                # is inside the measured window, which is checked above.
                'n': 1.0,
                'n_power': node.n,
                # One number cannot stand for a distribution: tau is the
                # slow end, which sets the low-frequency arc, and tau_min
                # travels beside it to give the range.
                'tau': node.tau_max,
                'tau_min': node.tau_min,
                'U': node.U,
            })

        elif isinstance(node, YG):
            # The film capacitance is a fitted parameter and the exact
            # omega -> inf limit of the element, so it needs no Hsu-Mansfeld /
            # Brug conversion. Whether that limit is inside the window is checked
            # later, on the winner alone (_plateau_notes).
            results.append({
                'type': 'YG',
                # NOT node.R_dc. The DC resistance of an exponential
                # conductivity profile is e^(1/p) times below the capacitive
                # corner - 2.7e45 Ohm for Zahner's own p = 0.01 example - so
                # it is an extrapolation of the model, not a fitted arc, and
                # must not be reported as the element's resistance. It
                # travels as R_dc instead, for reporting only. This is a
                # deliberate departure from DQ, whose R_pol is a measurable
                # arc. The ranking is the same either way: a YG with no
                # parallel R counts as blocking (infinite R) in the
                # largest-R barrier heuristic, which is where R_dc would
                # have put it too.
                'R': R_parallel,
                'C': node.C,
                # 1.0, like CC and DQ: above 1/(2*pi*tau) the element *is* a
                # capacitor, so C is a limit of the model rather than a
                # conversion whose reliability decays with an exponent.
                'n': 1.0,
                'tau': node.tau,
                'p': node.p,
                'R_dc': node.R_dc,
                'f_R': node.dc_corner_freq,
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
                # CC has no DC path of its own; a leakage branch beside it
                # (R_leak | CC) is what the largest-R heuristic must compare,
                # or a side relaxation would count as blocking and win
                'R': R_parallel,
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
    """Capacitance used to rank elements the R heuristic cannot separate.

    The fitted capacitance for C, K, CC, DQ and YG; the Hsu-Mansfeld C_eff for
    Q, so that a Q and a C sharing one resistance compare in the same unit.
    A Q always has its R here - one without it is not a candidate.
    """
    if 'C' in element:
        return float(element['C'])
    return _estimate_cpe_capacitance(element['Q'], element['n'], element['R'])


def _barrier_resistance(element: Dict[str, Any]) -> float:
    """Parallel resistance for the largest-R heuristic; no DC path is infinite."""
    return element['R'] if element['R'] is not None else np.inf


def _select_dielectric_element(
    candidates: List[Dict[str, Any]],
    warnings: List[str]
) -> Tuple[Dict[str, Any], str]:
    """
    Pick the element that carries the dielectric response.

    Caller has already applied the physical criterion (admittance ~ omega^n
    with n near 1); this only resolves which of several qualifying elements to
    report. The element with the largest parallel resistance is taken to be
    the compact barrier - the long-standing heuristic, which distinguishes an
    oxide barrier from a charge-transfer process. An element with no DC path
    beside it is blocking, i.e. its resistance is infinite, so it wins.

    The element type plays no part here. It decides only how the capacitance
    is read (exact for C, K, CC, DQ and YG; Hsu-Mansfeld/Brug for Q), not
    which element is the barrier: ranking by type let a 50 Ohm | 1 nF side
    arc beat a 10 MOhm | Q (n = 0.95) barrier.

    Elements sharing the largest resistance are not separable by it and are
    ranked by capacitance.

    Returns the element and a one-line statement of why it won, which the
    caller reports: the choice is a heuristic and has to be checkable.
    """
    dominant = max(candidates,
                   key=lambda e: (_barrier_resistance(e), _element_size(e)))
    R_max = _barrier_resistance(dominant)
    tied = [e for e in candidates if _barrier_resistance(e) == R_max]

    blocking = R_max == np.inf
    if len(tied) > 1:
        shared = ("with no DC path beside them" if blocking else
                  f"share one parallel resistance (R = {R_max:.1f} Ω)")
        warnings.append(
            f"{len(tied)} capacitive elements {shared} - using the largest, "
            f"{_element_size(dominant):.3e} F. Their individual values are not "
            "separately identifiable from the spectrum, so check the fit "
            "before relying on the split.")
    if blocking:
        return dominant, ("Selected as the barrier: no DC path beside it "
                          "(blocking), so its resistance exceeds any other "
                          "candidate's")
    return dominant, ("Selection assumes the largest-R element is the "
                      "compact oxide barrier (verify: a charge-transfer "
                      "process can also have the largest R)")


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


def _parallel_combination(resistances: List[float]) -> float:
    """Parallel combination [Ohm]; inf when every branch blocks, 0 if one shorts."""
    G_total = sum(1 / r if r > 0 else np.inf for r in resistances)
    return 1 / G_total if G_total > 0 else np.inf


def _branch_resistance(node, across: bool = True) -> float:
    """Resistance one branch of a Parallel puts beside its siblings [Ohm].

    inf when the branch blocks DC: a C, Q, CC or YG in series means no DC
    path, so the Debye branch R_rel - C_rel leaves C_geo blocking.

    A relaxation (Parallel, K, DQ) directly across the element (`across`) is
    part of the element's own parallel combination and counts with its DC
    resistance: C sees R in C | (R | C2), R in C | K(R, tau). One in series
    inside the branch is a separate arc, shorted by its own capacitance at
    this element's frequency: Q_ox in Q_ox | (R_ox - (Q_dl | R_ct)) sees R_ox,
    not R_ox + R_ct - unless it blocks DC, which blocks the branch.
    """
    if isinstance(node, (C, Q, CC, YG)):
        return np.inf
    if isinstance(node, R):
        return node.params[0]
    if isinstance(node, G):
        return 1 / node.G if node.G > 0 else np.inf
    if isinstance(node, Series):
        return sum(_branch_resistance(e, across=False) for e in node.elements)
    if isinstance(node, (K, DQ, Parallel)):
        if isinstance(node, K):
            R_dc = node.params[0]
        elif isinstance(node, DQ):
            R_dc = node.R_pol
        else:
            R_dc = _parallel_combination(
                [_branch_resistance(e) for e in node.elements])
        return R_dc if across or R_dc == np.inf else 0.0
    # ponytail: every other element (W, Wo, L, ...) counts as a short. Exact
    # for L; for a Warburg it keeps R_ct in the Randles branch R_ct - W, where
    # the DC-exact infinity would hide the arc. Ceiling: an element that
    # blocks DC or adds resistance and is not listed above is misread - add
    # it above if that matters.
    return 0.0


def _find_series_resistance(circuit) -> Optional[float]:
    """
    Sum of R and G elements on the series path of the circuit (outside any
    parallel combination) — the ohmic/electrolyte resistance Rs needed
    by the Brug formula. G contributes 1/G; G = 0 is an open series branch
    and contributes nothing.

    Returns None if no such element exists or the sum is not positive.
    """
    total = 0.0
    found = False

    def traverse(node):
        nonlocal total, found
        if isinstance(node, R):
            total += node.params[0]
            found = True
        elif isinstance(node, G) and node.G > 0:
            total += 1 / node.G
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
            fit_result = None
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
                fit_result = None
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
                else:  # Q
                    if dominant['n'] < CPE_N_RELIABLE_MIN:
                        warnings.append(
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
                        warnings.append("No series R element in circuit - "
                                        "Brug (2D) estimate not available")
                    elif R_s < BRUG_RS_MIN_OHM:
                        # A CPE with n < 1 mimics a series resistance at high
                        # frequency, so R_s is often unidentifiable and the fit
                        # drives it to the optimizer floor. Brug's
                        # C ~ R_s^((1-n)/n) would then be arbitrarily small.
                        warnings.append(
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

                if dominant['type'] == 'CC':
                    warnings.extend(_cc_capacitance_notes(dominant, frequencies))
                if dominant['type'] == 'DQ':
                    warnings.extend(_plateau_notes(
                        'DQ', dominant['tau_min'], frequencies, 'C_eff',
                        'leans on the fitted power law. Check the confidence '
                        'interval on tau_DQ.'))
                if dominant['type'] == 'YG':
                    warnings.extend(_plateau_notes(
                        'YG', dominant['tau'], frequencies, 'C',
                        'leans on the fitted p. Check the confidence '
                        'interval on p_YG.'))
                    warnings.extend(_yg_R_dc_note(dominant, frequencies))
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
                    'element_params': dict(dominant),
                    'mode': 'circuit',
                    'candidates': elements,
                    'selection_reason': selection_reason,
                    'n_hf_points': None,
                }

    # === Mode 2: Fallback - high-frequency estimate ===
    # Stated as plainly as a parameter sitting on its bound: the number below
    # is a spectral guess, not a fitted quantity, and nothing downstream
    # distinguishes the two once they are printed side by side.
    warnings.append(
        "NOT FROM THE FIT: the capacitance below is estimated directly from "
        "the spectrum (median of C = -1/(omega*Z'') over the top frequency "
        "decade), because the circuit offered no capacitive element to read it "
        "from. It carries no confidence interval and the thickness or "
        "permittivity derived from it is an order-of-magnitude figure - treat "
        "it as such even when it happens to land near the expected value.")
    warnings.append("For better accuracy, provide fitted circuit via fit_result")
    warnings.append("For multilayer (series) systems the high-frequency estimate "
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
        n_hf_points = int(C_values.size)

        # C_i is frequency-independent only when the capacitance dominates
        # (ωRC ≫ 1); a large spread means that assumption does not hold
        spread = float(np.max(C_values) / np.min(C_values))
        if spread > HF_C_SPREAD_MAX_RATIO:
            warnings.append(
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
            warnings.append("Positive imaginary impedance (inductive) - "
                            "result may be invalid")

        n_hf_points = None
        omega_hf = 2 * np.pi * f_max
        C_estimate = -1 / (omega_hf * Z_imag_hf)

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


def _add_penetration_depth(
    element_params: Dict[str, Any],
    thickness_nm: float
) -> None:
    """Turn Young-Göhr's p into the penetration depth it stands for, in nm.

    p = delta/d is the ratio the fit returns; the depth itself needs the
    thickness, which only exists once the capacitance has been through the
    parallel-plate model. Mutates in place, on both the thickness and the
    permittivity paths - they compute d differently but need the same number.

    This is the quantity the element exists to deliver: how far the
    conductivity reaches into the film. The CPE route has no equivalent.
    """
    if 'p' in element_params:
        element_params['delta_nm'] = element_params['p'] * thickness_nm


def analyze_oxide_layer(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    epsilon_r: float = DEFAULT_EPSILON_R,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[OxideAnalysisResult]:
    """
    Estimate oxide layer thickness from dominant capacitive element.

    Collects every capacitive element in the circuit (C, Q, K, CC, DQ), keeps
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
    warnings: List[str] = []

    # Get capacitance using the same element-selection logic as
    # analyze_oxide_layer (no thickness is computed here)
    extracted = _extract_capacitance(frequencies, Z, area_cm2, fit_result,
                                     warnings)

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
        mode=extracted['mode'],
        candidates=extracted['candidates'],
        selection_reason=extracted['selection_reason'],
        n_hf_points=extracted['n_hf_points'],
        area_cm2=area_cm2,
        warnings=warnings
    )


__all__ = ['analyze_oxide_layer', 'estimate_permittivity', 'OxideAnalysisResult']
