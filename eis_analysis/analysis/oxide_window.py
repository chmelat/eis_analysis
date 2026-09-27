"""
What the measured frequency window supports, for oxide analysis.

Notes on whether a dominant element's capacitance was measured or
extrapolated (Cole-Cole regime, DQ/YG plateaus, YG R_dc), and the
high-frequency spectral estimate used when the circuit offers no
capacitive element.
"""

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting.bounds import PARAMETER_BOUNDS, classify_bound_status
from .config import (
    HF_ESTIMATE_DECADE_FACTOR,
    HF_C_SPREAD_MAX_RATIO,
    HF_ZIMAG_MIN_REL,
    CC_WINDOW_EDGE_MARGIN_DECADES,
)

logger = logging.getLogger(__name__)


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




def _window_notes(
    dominant: Dict[str, Any],
    frequencies: NDArray[np.float64]
) -> List[str]:
    """Whether the window measured the dominant element's capacitance.

    Only CC, DQ and YG have a capacitance that is a limit of the model and
    so can lie outside the window; the other types get no note.
    """
    if dominant['type'] == 'CC':
        return _cc_capacitance_notes(dominant, frequencies)
    if dominant['type'] == 'DQ':
        return _plateau_notes(
            'DQ', dominant['tau_min'], frequencies, 'C_eff',
            'leans on the fitted power law. Check the confidence '
            'interval on tau_DQ.')
    if dominant['type'] == 'YG':
        return (_plateau_notes(
                    'YG', dominant['tau'], frequencies, 'C',
                    'leans on the fitted p. Check the confidence '
                    'interval on p_YG.')
                + _yg_R_dc_note(dominant, frequencies))
    return []

def _hf_capacitance_estimate(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    warnings: List[str]
) -> Optional[Tuple[float, Optional[int]]]:
    """
    Capacitance read directly from the spectrum, for a circuit with none.

    Returns (C [F], number of points behind the median, or None for the
    single-point fallback), or None if no capacitance could be read.
    Caveats are appended to `warnings`.
    """
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
    capacitive_mask = decade_mask & (Z.imag < -HF_ZIMAG_MIN_REL * np.abs(Z))

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

        if abs(Z_imag_hf) < HF_ZIMAG_MIN_REL * abs(Z[high_freq_idx]):
            logger.warning("Imaginary impedance too small at high frequency")
            return None

        if Z_imag_hf > 0:
            warnings.append("Positive imaginary impedance (inductive) - "
                            "result may be invalid")

        n_hf_points = None
        omega_hf = 2 * np.pi * f_max
        C_estimate = -1 / (omega_hf * Z_imag_hf)

    return C_estimate, n_hf_points
