"""
Scalar-quantity estimation for DRT analysis.

Quantities derived from the impedance data and the recovered gamma(tau):
high-frequency resistance R_inf, polarization resistance R_pol, per-peak
resistances, and the effective-bins shape metric.
"""

import numpy as np
import logging
from typing import Any, Dict, Optional, List, Tuple
from numpy.typing import NDArray

from .results import RinfEstimate
from ..fitting.config import DRT_PEAK_EDGE_DECADES
from ..rinf_estimation import estimate_rinf_with_inductance

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


def _flag_boundary_peaks(tau: NDArray, peaks: List[Dict[str, Any]],
                        tau_key: str) -> None:
    """
    Annotate each peak with its distance to the nearer end of the tau grid.

    The grid spans exactly the measured window, so a peak close to either end
    has one flank that no measurement constrains: its position is poorly
    localized and the basin integral behind ``R_estimate`` is truncated by the
    end of the array. Sets ``edge_distance_decades`` and the
    ``boundary_sensitive`` flag on every peak in place, on both the scipy
    dicts (keyed ``tau``) and the GMM dicts (keyed ``tau_center``).
    """
    log_lo, log_hi = np.log10(tau[0]), np.log10(tau[-1])

    for peak in peaks:
        log_tau = float(np.log10(peak[tau_key]))
        distance = float(min(log_tau - log_lo, log_hi - log_tau))
        peak['edge_distance_decades'] = distance
        peak['boundary_sensitive'] = bool(distance < DRT_PEAK_EDGE_DECADES)


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
    if R_pol <= 1e-10:
        return 0.0, None

    def falling_run_mass(g: NDArray) -> float:
        """Mass of the descending run starting at g[0]; 0.0 if g rises."""
        i = 0
        while i + 1 < len(g) and g[i + 1] < g[i]:
            i += 1
        return float(np.sum(g[:i + 1])) if i > 0 else 0.0

    low = falling_run_mass(gamma) * d_ln_tau / R_pol
    high = falling_run_mass(gamma[::-1]) * d_ln_tau / R_pol

    if max(low, high) <= 0.0:
        return 0.0, None
    return (low, 'low') if low >= high else (high, 'high')


def _effective_bins(gamma: NDArray) -> float:
    """Participation ratio N_eff = (sum gamma)^2 / sum(gamma^2).

    ~1 for a single-bin spike, grows to tens for a smooth distribution.
    Used to flag DRT too sparse/spiky for peak-shape analysis (audit F3).
    """
    s = float(np.sum(gamma))
    denom = float(np.sum(gamma**2))
    return (s * s) / denom if denom > 0 else 0.0


def _estimate_r_inf(frequencies: NDArray, Z: NDArray,
                    use_rl_fit: bool = False,
                    r_inf_preset: Optional[float] = None) -> RinfEstimate:
    """
    Estimate high-frequency resistance R_inf.

    Returns structured RinfEstimate with all diagnostics.
    """
    n_avg = min(5, max(1, len(frequencies) // 10))
    high_freq_indices = np.argsort(frequencies)[-n_avg:]
    R_inf_median = float(np.median(Z.real[high_freq_indices]))

    if r_inf_preset is not None:
        return RinfEstimate(
            R_inf=r_inf_preset,
            method='preset',
            R_inf_median=R_inf_median
        )

    if use_rl_fit:
        try:
            fit, fig_rl = estimate_rinf_with_inductance(frequencies, Z, plot=True)

            warnings = list(fit.warnings)
            if fit.L_nH > 500:
                warnings.append(f"High inductance L = {fit.L_nH:.1f} nH detected")

            return RinfEstimate(
                R_inf=fit.R_inf,
                method='rl_fit',
                R_inf_median=R_inf_median,
                figure=fig_rl,
                behavior=fit.behavior,
                n_points_used=fit.n_points_used,
                R_squared=fit.R_squared,
                L_nH=fit.L_nH,
                warnings=warnings
            )
        except Exception as e:
            logger.debug(f"R_inf fit failed: {e}, using median fallback")

    # Median method (default)
    return RinfEstimate(
        R_inf=R_inf_median,
        method='median',
        R_inf_median=R_inf_median,
        n_points_used=n_avg
    )
