"""
Scalar-quantity estimation for DRT analysis.

Quantities derived from the impedance data and the recovered gamma(tau):
high-frequency resistance R_inf, polarization resistance R_pol, per-peak
resistances, and the effective-bins shape metric.
"""

import numpy as np
import logging
from typing import Optional, List
from numpy.typing import NDArray

from .results import RinfEstimate
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
                    use_voigt_fit: bool = False,
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

    if use_rl_fit or use_voigt_fit:
        try:
            R_inf, L_fit, circuit_rl, diag, fig_rl = estimate_rinf_with_inductance(
                frequencies, Z, use_voigt_fit=use_voigt_fit, verbose=False, plot=True
            )

            warnings = []
            if diag.get('warnings'):
                warnings.extend(diag['warnings'])
            if diag.get('L_nH', 0) > 500:
                warnings.append(f"High inductance L = {diag['L_nH']:.1f} nH detected")

            return RinfEstimate(
                R_inf=R_inf,
                method='voigt_fit' if use_voigt_fit else 'rl_fit',
                R_inf_median=R_inf_median,
                figure=fig_rl,
                behavior=diag.get('behavior'),
                n_points_used=diag.get('n_points_used'),
                R_squared=diag.get('r_squared'),
                L_nH=diag.get('L_nH'),
                R_ct=diag.get('R_ct'),
                C_nF=diag.get('C_nF'),
                f_characteristic=diag.get('f_characteristic'),
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
