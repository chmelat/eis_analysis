"""
Core DRT (Distribution of Relaxation Times) calculation.

Clean design: No logging in core functions, all diagnostics returned as data.
CLI layer is responsible for user output.

This module is the orchestrator: it keeps the public ``calculate_drt`` entry
point and ties together the pipeline stages, which live in sibling modules
(``results``, ``estimation``, ``linear_system``, ``plotting``, ``peaks``). Those
symbols are re-exported below so they remain importable from ``drt.core``.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List, Dict
from numpy.typing import NDArray
from scipy.signal import find_peaks

from .results import (
    RinfEstimate,
    LambdaSelection,
    NNLSSolution,
    DRTDiagnostics,
    DRTMatrices,
    DRTResult,
    LambdaProbePoint,
    PeakStability,
    StabilityDiagnostics,
)
from .estimation import (
    _rpol_from_gamma,
    _estimate_peak_resistance,
    _edge_pile_up,
    _effective_bins,
    _estimate_r_inf,
    _flag_boundary_peaks,
)
from .linear_system import (
    _validate_frequencies,
    _build_drt_matrices,
    _select_lambda,
    _solve_nnls,
)
from .plotting import _create_visualization
from .peaks import gmm_peak_detection
from .stability import probe_lambda_stability
from ..fitting.config import (DRT_PEAK_HEIGHT_THRESHOLD, DRT_MIN_EFFECTIVE_BINS,
                             DRT_EDGE_BIN_RPOL_FRACTION, DRT_PEAK_EDGE_DECADES,
                             GMM_N_COMPONENTS_RANGE)

logger = logging.getLogger(__name__)

# Public API re-exported from the pipeline submodules so it stays importable
# from ``drt.core`` (and via ``drt/__init__.py``) after the split.
__all__ = [
    'calculate_drt',
    'DRTResult',
    'DRTDiagnostics',
    'RinfEstimate',
    'LambdaSelection',
    'NNLSSolution',
    'DRTMatrices',
    'LambdaProbePoint',
    'PeakStability',
    'StabilityDiagnostics',
]



# =============================================================================
# Peak Detection
# =============================================================================

def _detect_peaks(tau: NDArray, gamma: NDArray,
                  peak_method: str,
                  gmm_bic_threshold: float = 10.0,
                  n_data: Optional[int] = None
                  ) -> Tuple[Optional[List[Dict]], Optional[List[float]], Optional[List[Dict]]]:
    """
    Detect peaks in DRT spectrum.

    n_data: počet skutečných měření (frekvencí) pro penalizaci BIC v GMM.

    Returns:
        (gmm_peaks, bic_scores, scipy_peaks)
    """
    use_gmm = (peak_method == 'gmm')

    # Always calculate scipy peaks for diagnostics
    peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD)
    peak_resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)

    scipy_peaks = []
    for i, idx in enumerate(peaks_idx):
        R_peak = peak_resistances[i] if i < len(peak_resistances) else 0.0
        scipy_peaks.append({
            'index': int(idx),
            'tau': float(tau[idx]),
            'frequency': float(1/(2 * np.pi * tau[idx])),
            'R_estimate': float(R_peak)
        })

    _flag_boundary_peaks(tau, scipy_peaks, 'tau')

    if use_gmm:
        peaks_result, gmm_model, bic_scores = gmm_peak_detection(
            tau, gamma, n_components_range=GMM_N_COMPONENTS_RANGE,
            bic_threshold=gmm_bic_threshold, n_data=n_data
        )

        if len(peaks_result) == 0 or gmm_model is None:
            # GMM failed, scipy_peaks available as fallback
            return None, None, scipy_peaks

        _flag_boundary_peaks(tau, peaks_result, 'tau_center')

        return peaks_result, bic_scores, scipy_peaks

    return None, None, scipy_peaks


# =============================================================================
# Main Function
# =============================================================================

def calculate_drt(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_tau: int = 100,
    lambda_reg: Optional[float] = None,
    auto_lambda: bool = False,
    normalize_rpol: bool = False,
    use_rl_fit: bool = False,
    peak_method: str = 'scipy',
    r_inf_preset: Optional[float] = None,
    gmm_bic_threshold: float = 10.0,
    lambda_probe: bool = False
) -> DRTResult:
    """
    Calculate DRT (Distribution of Relaxation Times) using Tikhonov regularization.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    n_tau : int
        Number of tau points (default: 100)
    lambda_reg : float, optional
        Regularization parameter
    auto_lambda : bool
        Auto-select lambda using hybrid GCV + L-curve
    normalize_rpol : bool
        Normalize gamma by R_pol
    use_rl_fit : bool
        Use R-L-K fit for R_inf estimation instead of the HF median
    peak_method : str
        Peak detection method ('scipy' or 'gmm')
    r_inf_preset : float, optional
        Preset R_inf value
    gmm_bic_threshold : float
        BIC threshold for GMM peak detection
    lambda_probe : bool
        Re-solve the DRT at lambdas around the selected one and report
        per-peak stability (see drt.stability)

    Returns
    -------
    DRTResult
        Complete analysis result with all diagnostics
    """
    f_min, f_max = float(frequencies.min()), float(frequencies.max())
    freq_range_ratio = f_max / f_min

    # === Step 1: Validate frequencies ===
    try:
        _validate_frequencies(frequencies)
    except ValueError as e:
        logger.error(str(e))
        return DRTResult()

    # === Step 2: R_inf Estimation ===
    rinf_est = _estimate_r_inf(
        frequencies, Z,
        use_rl_fit=use_rl_fit,
        r_inf_preset=r_inf_preset
    )
    R_inf = rinf_est.R_inf

    # === Step 3: Build Matrices ===
    matrices = _build_drt_matrices(frequencies, Z, R_inf, n_tau)

    # === Step 4: Select Lambda ===
    lambda_sel = _select_lambda(
        matrices.A, matrices.b, matrices.L,
        lambda_reg, auto_lambda
    )

    # === Step 5: Solve NNLS ===
    nnls_result = _solve_nnls(
        matrices.A, matrices.b, matrices.L,
        lambda_sel.lambda_value, n_tau, Z
    )

    if not nnls_result.success:
        return DRTResult(
            R_inf=R_inf,
            figure_rinf=rinf_est.figure,
            diagnostics=DRTDiagnostics(
                freq_min=f_min, freq_max=f_max,
                freq_range_ratio=freq_range_ratio,
                n_points=len(frequencies),
                n_tau=n_tau,
                condition_number=matrices.condition_number,
                d_ln_tau=matrices.d_ln_tau,
                rinf=rinf_est,
                lambda_sel=lambda_sel,
                nnls=nnls_result,
                R_pol_from_data=0.0,
                R_pol_from_gamma=0.0,
                normalized=False,
                reconstruction_error_rel=0.0,
                peak_method=peak_method,
                n_peaks=0
            )
        )

    # success=True guarantees a valid gamma (see _solve_nnls); assert narrows
    # the Optional for the type checker.
    gamma = nnls_result.gamma
    assert gamma is not None

    # === Step 6: R_pol Calculation & Normalization ===
    n_avg = min(5, max(1, len(frequencies) // 10))
    low_freq_indices = np.argsort(frequencies)[:n_avg]
    R_dc = float(np.median(Z.real[low_freq_indices]))
    R_pol_from_data = R_dc - R_inf
    R_pol_from_gamma = _rpol_from_gamma(gamma, matrices.d_ln_tau)

    gamma_original = None
    normalized = False
    if normalize_rpol:
        if R_pol_from_gamma > 1e-10:
            gamma_original = gamma.copy()
            gamma = gamma / R_pol_from_gamma
            normalized = True

    # Physical (unnormalized) gamma for everything downstream that must stay
    # in Ohm — peak R_estimates, reconstruction, shape metrics — even when the
    # returned gamma is normalized by R_pol (audit 2026-07-02 finding 2.2).
    gamma_physical = gamma_original if normalized else gamma
    assert gamma_physical is not None  # set whenever normalized; narrows Optional

    # === Step 7: Peak Detection ===
    peaks_result, bic_scores, scipy_peaks = _detect_peaks(
        matrices.tau, gamma_physical, peak_method, gmm_bic_threshold,
        n_data=len(frequencies)
    )

    # The peaks the run reports: GMM components when GMM ran and succeeded,
    # otherwise the scipy maxima. Everything downstream counts this set.
    reported_peaks = peaks_result or scipy_peaks or []
    n_peaks = len(reported_peaks)

    # === Step 8: Reconstruction & Error ===
    Z_reconstructed = R_inf + (matrices.A_re + 1j * matrices.A_im) @ gamma_physical
    rel_error = float(np.mean(np.abs(Z - Z_reconstructed) / np.abs(Z)) * 100)

    # Add warning for high reconstruction error
    if rel_error > 10.0:
        nnls_result.warnings.append(
            f"High reconstruction error ({rel_error:.1f}%) - DRT model may not be suitable"
        )
    elif rel_error > 5.0:
        nnls_result.warnings.append(
            f"Elevated reconstruction error ({rel_error:.1f}%)"
        )

    # === Step 8b: Window-edge diagnostics ===
    # Pile-up at an end of the tau grid: NNLS has nowhere to put response
    # whose time constant lies outside the measured window, so it heaps gamma
    # up against the boundary instead.
    edge_pile_up_fraction, edge_pile_up_end = _edge_pile_up(
        gamma_physical, matrices.d_ln_tau, R_pol_from_gamma
    )
    if edge_pile_up_fraction > DRT_EDGE_BIN_RPOL_FRACTION:
        nnls_result.warnings.append(
            f"{edge_pile_up_fraction*100:.0f}% of R_pol is heaped against the "
            f"{'fast' if edge_pile_up_end == 'low' else 'slow'} end of the tau "
            f"window - likely response from outside it (series inductance, "
            f"unresolved tail, a process slower than the lowest frequency)"
        )
        # That mass is not left unattributed: the scipy basin partition runs to
        # the end of the array, so the outermost peak on the loaded side
        # absorbs it, while the GMM path divides R_pol by component weight, so
        # every component carries a share. Mark whichever applies - a distance
        # test cannot catch this, the contaminated peak can sit decades away.
        if reported_peaks:
            contaminated = (reported_peaks if peaks_result else
                            [reported_peaks[0] if edge_pile_up_end == 'low'
                             else reported_peaks[-1]])
            for peak in contaminated:
                peak['edge_contaminated'] = True
            nnls_result.warnings.append(
                f"R_estimate of {len(contaminated)} of {len(reported_peaks)} "
                f"peaks includes that out-of-window mass and is inflated"
            )

    # Peaks too close to the window edge to be localized or integrated fully.
    n_boundary_peaks = sum(1 for p in reported_peaks
                           if p.get('boundary_sensitive'))
    if n_boundary_peaks:
        nnls_result.warnings.append(
            f"{n_boundary_peaks} of {len(reported_peaks)} peaks lie within "
            f"{DRT_PEAK_EDGE_DECADES} decade of the measured tau window edge - "
            f"position and R_estimate are only partly supported by data"
        )

    # === Step 8c: Shape-quality diagnostics (F3) ===
    # Warn if the DRT is too sparse/spiky for peak-shape analysis, or if
    # auto-lambda hit the search-range edge (regularization too low). Advisory
    # only - gamma and detected peaks are unchanged.
    n_eff = _effective_bins(gamma_physical)
    if n_eff < DRT_MIN_EFFECTIVE_BINS:
        nnls_result.warnings.append(
            f"DRT is sparse/spiky (effective bins {n_eff:.1f} < "
            f"{DRT_MIN_EFFECTIVE_BINS:.0f}); peak-shape analysis may be "
            f"unreliable - consider a higher lambda"
        )
    if lambda_sel.lambda_at_edge or lambda_sel.corner_at_edge:
        nnls_result.warnings.append(
            f"Auto-lambda at search-range edge (lambda="
            f"{lambda_sel.lambda_value:.2e}); regularization may be too low "
            f"for reliable DRT shape"
        )

    # === Step 8d: Lambda-probe peak stability (opt-in) ===
    # Track the reported peaks across lambdas around the selected one; peaks
    # that vanish or drift under a modest lambda change are likely
    # regularization artifacts. Reference peaks and probe run on the physical
    # gamma [Ohm].
    stability = None
    probe_curves = None
    if lambda_probe:
        if peaks_result:
            reference_peaks = [(p['tau_center'], p['R_estimate'])
                               for p in peaks_result]
        else:
            reference_peaks = [(p['tau'], p['R_estimate'])
                               for p in (scipy_peaks or [])]
        stability = probe_lambda_stability(
            matrices, lambda_sel.lambda_value, reference_peaks, Z, R_inf,
            n_tau
        )
        # Overlay curves for the DRT figure; match the normalization of the
        # displayed gamma.
        scale = R_pol_from_gamma if normalized else 1.0
        probe_curves = [
            (p.lambda_value, p.gamma / scale)
            for p in stability.probe_points
            if p.success and p.gamma is not None
        ]

    # === Step 9: Visualization ===
    fig = _create_visualization(
        matrices.tau, gamma, gamma_original,
        Z, Z_reconstructed,
        lambda_sel.lambda_value,
        normalized, peak_method,
        peaks_result, bic_scores,
        probe_curves=probe_curves
    )

    # === Build diagnostics ===
    diagnostics = DRTDiagnostics(
        freq_min=f_min,
        freq_max=f_max,
        freq_range_ratio=freq_range_ratio,
        n_points=len(frequencies),
        n_tau=n_tau,
        condition_number=nnls_result.condition_number,
        d_ln_tau=matrices.d_ln_tau,
        rinf=rinf_est,
        lambda_sel=lambda_sel,
        nnls=nnls_result,
        R_pol_from_data=R_pol_from_data,
        R_pol_from_gamma=R_pol_from_gamma,
        normalized=normalized,
        reconstruction_error_rel=rel_error,
        peak_method=peak_method,
        n_peaks=n_peaks,
        scipy_peaks=scipy_peaks,
        n_effective_bins=n_eff,
        edge_pile_up_fraction=edge_pile_up_fraction,
        edge_pile_up_end=edge_pile_up_end,
        n_boundary_peaks=n_boundary_peaks,
        stability=stability
    )

    return DRTResult(
        tau=matrices.tau,
        gamma=gamma,
        gamma_original=gamma_original,
        figure=fig,
        figure_rinf=rinf_est.figure,
        peaks=peaks_result,
        bic_scores=bic_scores,
        R_inf=R_inf,
        R_pol=R_pol_from_gamma,
        lambda_used=lambda_sel.lambda_value,
        reconstruction_error=rel_error,
        diagnostics=diagnostics
    )
