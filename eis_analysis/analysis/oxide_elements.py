"""
Capacitive elements of a fitted circuit, for oxide analysis.

Finds every capacitive element with the resistance beside it, picks the
dielectric one, and converts a CPE (Q) to an effective capacitance.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting.circuit_elements import R, C, G, Q, K, CC, DQ, YG
from ..fitting.circuit_builder import Series, Parallel
from .config import BRUG_RS_MIN_OHM
from .oxide_window import _cc_capacitance_regime

logger = logging.getLogger(__name__)


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
            # _cc_capacitance_regime and _cc_capacitance_notes).
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

    Commonly used for a normal (through-layer) distribution of time
    constants, but not exact for one: Hirschorn et al. (2010) showed it
    does not in general return the film capacitance there, and derived the
    power-law model (oxide_power_law) for that case. For a surface
    distribution the Brug (1984) formula applies instead.

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
    Reported alongside the Hsu-Mansfeld value as a comparison;
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
    # ponytail: every other element counts as a short. Exact for L; W, Ws, Wo
    # are the semi-infinite W at the arc's frequency, so R_ct - W keeps R_ct
    # (their DC limits, inf or R_W, belong to the slower diffusion). Ceiling:
    # an element that blocks DC or adds resistance at the arc's frequency and
    # is not listed above is misread - add it above if that matters.
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


def _q_capacitance(
    q: Dict[str, Any],
    circuit,
    warnings: List[str]
) -> Tuple[float, Optional[float]]:
    """
    Effective capacitance of a Q candidate: Hsu-Mansfeld C_eff, and the
    Brug (2D) comparison value, or None when the circuit has no identified
    series resistance. Caveats are appended to `warnings`.
    """
    C_eff = _estimate_cpe_capacitance(q['Q'], q['n'], q['R'])

    # Brug (2D) comparison estimate - needs series resistance
    C_eff_brug = None
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
        C_eff_brug = _estimate_cpe_capacitance_brug(q['Q'], q['n'], q['R'], R_s)
    return C_eff, C_eff_brug


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
