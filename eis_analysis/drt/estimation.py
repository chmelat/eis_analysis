"""
Scalar-quantity estimation for DRT analysis.

Quantities derived from the impedance data and the recovered gamma(tau):
high-frequency resistance R_inf, polarization resistance R_pol, per-peak
resistances, and the effective-bins shape metric.
"""

import numpy as np
import logging
from typing import Any, Dict, Optional, List, Tuple, Union
from numpy.typing import NDArray

from .results import RinfEstimate
from ..fitting.config import DRT_INDUCTANCE_DECADES, DRT_PEAK_EDGE_DECADES
from ..rinf_estimation import hf_median

logger = logging.getLogger(__name__)


def _rpol_from_gamma(gamma: NDArray, d_ln_tau: float) -> float:
    """Polarization resistance R_pol = sum(gamma) * d_ln_tau (rectangle rule).

    This matches the DRT kernel, which integrates with the rectangle rule:
    A_re -> d_ln_tau as omega -> 0, so the model's own DC limit is
    Z'(0) - R_inf = sum_m gamma_m * d_ln_tau. Using the same quadrature here
    (rather than trapz) keeps R_pol consistent with the reconstructed model
    (audit finding F10).
    """
    return float(np.sum(gamma) * d_ln_tau)


def _estimate_peak_resistance(tau: NDArray, gamma: NDArray,
                               peak_indices: NDArray) -> List[float]:
    """
    Estimate resistance for each peak by integrating gamma over a partition
    of the tau axis.

    The tau axis is split at the valleys (gamma minima) between consecutive
    peaks into disjoint half-open segments, so each grid point is assigned to
    exactly one peak. Each segment is integrated with the rectangle rule
    (consistent with the DRT kernel, F10), so sum(R_i) equals the total R_pol
    over the spanned range exactly — unlike per-peak threshold windows, which
    double-count the overlap region of adjacent peaks.
    """
    if len(peak_indices) == 0:
        return []

    ln_tau = np.log(tau)
    d_ln_tau = float(np.mean(np.diff(ln_tau)))
    peaks = np.sort(peak_indices)

    # Partition boundaries: start of the array plus the valley (argmin) between
    # each pair of consecutive peaks, plus the end. Segments are half-open
    # [bounds[j], bounds[j+1]) so every grid point belongs to exactly one peak.
    bounds = [0]
    for j in range(len(peaks) - 1):
        lo, hi = int(peaks[j]), int(peaks[j + 1])
        valley = lo + int(np.argmin(gamma[lo:hi + 1]))
        bounds.append(valley)
    bounds.append(len(gamma))

    resistances = []
    for j in range(len(peaks)):
        left, right = bounds[j], bounds[j + 1]
        R_peak = _rpol_from_gamma(gamma[left:right], d_ln_tau) if right > left else 0.0
        resistances.append(R_peak)

    return resistances


def _flag_boundary_peaks(tau_window: Tuple[float, float],
                         peaks: List[Dict[str, Any]], tau_key: str) -> None:
    """
    Annotate each peak with its distance to the nearer end of the measured window.

    A peak close to either end has one flank that no measurement constrains:
    its position is poorly localized and its ``R_estimate`` uncertain. The
    distance is signed - negative for a peak past the window, which only an
    extended grid allows - and such a peak also gets ``outside_window``. Sets
    ``edge_distance_decades``, ``boundary_sensitive`` and ``outside_window``
    on every peak in place, on both the scipy dicts (keyed ``tau``) and the
    GMM dicts (keyed ``tau_center``).
    """
    log_lo, log_hi = np.log10(tau_window[0]), np.log10(tau_window[1])

    for peak in peaks:
        log_tau = float(np.log10(peak[tau_key]))
        distance = float(min(log_tau - log_lo, log_hi - log_tau))
        peak['edge_distance_decades'] = distance
        peak['boundary_sensitive'] = bool(distance < DRT_PEAK_EDGE_DECADES)
        peak['outside_window'] = bool(distance < 0)


def _edge_pile_up_fractions(gamma: NDArray, d_ln_tau: float,
                            R_pol: float) -> Tuple[float, float]:
    """Share of R_pol in the falling run at the (fast, slow) end, see _edge_pile_up."""
    if R_pol <= 1e-10:
        return 0.0, 0.0

    def falling_run_mass(g: NDArray) -> float:
        """Mass of the descending run starting at g[0]; 0.0 if g rises."""
        i = 0
        while i + 1 < len(g) and g[i + 1] < g[i]:
            i += 1
        return float(np.sum(g[:i + 1])) if i > 0 else 0.0

    return (falling_run_mass(gamma) * d_ln_tau / R_pol,
            falling_run_mass(gamma[::-1]) * d_ln_tau / R_pol)


def _edge_pile_up(gamma: NDArray, d_ln_tau: float,
                  R_pol: float) -> Tuple[float, Optional[str]]:
    """
    Share of R_pol sitting in a falling run at an end of the tau grid.

    Non-negative NNLS has nowhere to put response whose time constant lies
    outside the measured window, so it piles gamma up against the first or
    last bin instead. The signature is a gamma that is already falling as it
    leaves the boundary: a relaxation inside the window makes gamma *rise*
    from the edge towards its maximum, so a descending run at the edge means
    the maximum lies outside it.

    Measuring the whole run rather than the outermost bin keeps the number
    independent of ``n_tau``: the pile-up lobe has a width in decades, and a
    finer grid merely spreads the same mass over more bins.

    Returns ``(fraction, end)`` with ``end`` the loaded side ('low' for the
    fast end, 'high' for the slow end), or ``(0.0, None)`` when neither end
    is loaded or R_pol is too small for the ratio to mean anything.
    """
    low, high = _edge_pile_up_fractions(gamma, d_ln_tau, R_pol)
    if max(low, high) <= 0.0:
        return 0.0, None
    return (low, 'low') if low >= high else (high, 'high')


def _extrapolated_fraction(tau: NDArray, gamma: NDArray, d_ln_tau: float,
                           R_pol: float, tau_max: float) -> float:
    """Share of R_pol placed past the slow end of the measured window (tau > tau_max)."""
    if R_pol <= 1e-10:
        return 0.0
    # Tolerance: the window's own last grid point equals tau_max up to rounding.
    beyond = tau > tau_max * (1 + 1e-9)
    return float(np.sum(gamma[beyond]) * d_ln_tau / R_pol)


def _inductance_choice(frequencies: NDArray, Z: NDArray,
                       inductance: Union[bool, str]) -> Tuple[bool, Optional[str]]:
    """Whether to model a series L, and for 'auto' why (DRT_INDUCTANCE_DECADES)."""
    if inductance != 'auto':
        return bool(inductance), None
    top = frequencies >= frequencies.max() / 10**DRT_INDUCTANCE_DECADES
    n_inductive = int(np.sum(Z.imag[top] > 0))
    if n_inductive:
        return True, f"{n_inductive} point(s) with Im(Z) > 0 in the top decade"
    return False, "no point with Im(Z) > 0 in the top decade"


def _lf_rc_ratio(frequencies: NDArray, Z: NDArray) -> float:
    """
    r = (-dZ'/d ln omega) / (-Z'') over the four lowest frequencies.

    Separates a relaxation just past the window (r ~ 0.2-1) from a capacitive
    end (r -> 0), see DRT_LF_RC_RATIO_MIN. Four points: enough for a slope fit
    to average single-point noise, few enough to stay at the low-frequency
    end. Returns inf when -Z'' is not positive there - no capacitive response,
    so nothing for the guard to reject.
    """
    idx = np.argsort(frequencies)[:4]
    minus_z_imag = -float(np.mean(Z.imag[idx]))
    if len(idx) < 2 or minus_z_imag <= 0:
        return float('inf')
    slope = np.polyfit(np.log(2 * np.pi * frequencies[idx]), Z.real[idx], 1)[0]
    return float(-slope / minus_z_imag)


def _effective_bins(gamma: NDArray) -> float:
    """Participation ratio N_eff = (sum gamma)^2 / sum(gamma^2).

    ~1 for a single-bin spike, grows to tens for a smooth distribution.
    Used to flag DRT too sparse/spiky for peak-shape analysis (audit F3).
    """
    s = float(np.sum(gamma))
    denom = float(np.sum(gamma**2))
    return (s * s) / denom if denom > 0 else 0.0


def _estimate_r_inf(frequencies: NDArray, Z: NDArray,
                    r_inf_preset: Optional[float] = None) -> RinfEstimate:
    """R_inf from the caller (preset) or the HF median, with the median kept for comparison."""
    R_inf_median, n_median = hf_median(frequencies, Z)
    if r_inf_preset is not None:
        return RinfEstimate(R_inf=r_inf_preset, method='preset',
                            R_inf_median=R_inf_median)
    return RinfEstimate(R_inf=R_inf_median, method='median',
                        R_inf_median=R_inf_median, n_points_used=n_median)
