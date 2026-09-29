"""
Core DRT (Distribution of Relaxation Times) calculation.

Clean design: No logging in core functions, all diagnostics returned as data.
CLI layer is responsible for user output.

This module is the orchestrator: it keeps the public ``calculate_drt`` entry
point and ties together the pipeline stages, which live in sibling modules
(``results``, ``estimation``, ``linear_system``, ``extension``, ``plotting``,
``peaks``). Those
symbols are re-exported below so they remain importable from ``drt.core``.
"""

import numpy as np
import logging
from typing import Tuple, Optional, List, Dict, Union
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
    _extrapolated_fraction,
    _inductance_choice,
)
from .linear_system import _reconstruct, _validate_frequencies
from .plotting import _create_visualization
from .peaks import gmm_peak_detection
from .stability import probe_lambda_stability
from .extension import _solve_with_extension
from ..fitting.config import (DRT_PEAK_HEIGHT_THRESHOLD, DRT_MIN_EFFECTIVE_BINS,
                             DRT_EDGE_BIN_RPOL_FRACTION, DRT_PEAK_EDGE_DECADES,
                             DRT_EXTRAPOLATED_RPOL_FRACTION, GMM_N_COMPONENTS_RANGE)

logger = logging.getLogger(__name__)

# Data weightings understood by fitting.diagnostics.compute_weights.
DRT_WEIGHTINGS = ('uniform', 'sqrt', 'modulus', 'proportional')

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
                  n_data: Optional[int] = None,
                  *, tau_window: Tuple[float, float]
                  ) -> Tuple[Optional[List[Dict]], Optional[List[float]], Optional[List[Dict]]]:
    """
    Detect peaks in DRT spectrum.

    n_data: počet skutečných měření (frekvencí) pro penalizaci BIC v GMM.
    tau_window: měřené okno pro okrajové příznaky píků.

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

    _flag_boundary_peaks(tau_window, scipy_peaks, 'tau')

    if use_gmm:
        peaks_result, gmm_model, bic_scores = gmm_peak_detection(
            tau, gamma, n_components_range=GMM_N_COMPONENTS_RANGE,
            bic_threshold=gmm_bic_threshold, n_data=n_data
        )

        if len(peaks_result) == 0 or gmm_model is None:
            # GMM failed, scipy_peaks available as fallback
            return None, None, scipy_peaks

        _flag_boundary_peaks(tau_window, peaks_result, 'tau_center')

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
    peak_method: str = 'scipy',
    r_inf_preset: Optional[float] = None,
    gmm_bic_threshold: float = 10.0,
    lambda_probe: bool = False,
    weighting: str = 'sqrt',
    tau_extend_decades: Union[float, str] = 0.0,
    inductance: Union[bool, str] = 'auto'
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
    peak_method : str
        Peak detection method ('scipy' or 'gmm')
    r_inf_preset : float, optional
        Preset R_inf value, e.g. from `estimate_rinf`; the HF median otherwise
    gmm_bic_threshold : float
        BIC threshold for GMM peak detection
    lambda_probe : bool
        Re-solve the DRT at lambdas around the selected one and report
        per-peak stability (see drt.stability)
    weighting : str
        Weighting of the data term from the measured |Z|: 'sqrt' (1/sqrt|Z|,
        default), 'uniform' (unweighted, before v0.38), 'modulus' (1/|Z|) or
        'proportional' (1/|Z|^2). Choice rationale: README, DRT analysis.
    tau_extend_decades : float or 'auto'
        Extend the tau grid this many decades past the slow end of the
        measured window, at the same log spacing (``n_tau`` points stay on
        the window). 0 (default) keeps the grid on the window. 'auto' extends
        only when that resolves a slow-end pile-up, see
        DRT_TAU_EXTEND_STEPS and DRT_LF_RC_RATIO_MIN; the choice is reported
        in ``diagnostics.tau_extend_note``. Mass past the window is an
        extrapolation, reported as ``R_pol_extrapolated_fraction``.
    inductance : bool or 'auto'
        Model a series inductance j*omega*L next to R_inf (unregularized,
        L >= 0, returned as ``L_series``). Without it the DRT cannot produce
        Im(Z) > 0 and an inductive high-frequency end deforms gamma. 'auto'
        (default) adds it only when the top decade has a point with
        Im(Z) > 0 (DRT_INDUCTANCE_DECADES); the choice is reported in
        ``diagnostics.inductance_note``.

    Returns
    -------
    DRTResult
        Complete analysis result with all diagnostics
    """
    if weighting not in DRT_WEIGHTINGS:
        # compute_weights would fall back to uniform with only a log line,
        # while diagnostics reported the misspelled name as if it applied.
        raise ValueError(f"Unknown DRT weighting '{weighting}', "
                         f"expected one of {DRT_WEIGHTINGS}")
    if tau_extend_decades != 'auto' and not (
            isinstance(tau_extend_decades, (int, float))
            and np.isfinite(tau_extend_decades) and tau_extend_decades >= 0):
        raise ValueError(f"tau_extend_decades must be a finite number >= 0 or 'auto', "
                         f"got {tau_extend_decades!r}")
    if inductance not in (True, False, 'auto'):
        raise ValueError(f"inductance must be True, False or 'auto', got {inductance!r}")

    f_min, f_max = float(frequencies.min()), float(frequencies.max())
    freq_range_ratio = f_max / f_min

    # === Step 1: Validate input ===
    try:
        _validate_frequencies(frequencies)
        # A single NaN would turn every data weight into NaN (compute_weights
        # divides by their mean) and crash the SVD in _build_drt_matrices.
        if not np.all(np.isfinite(Z)):
            raise ValueError("Impedance contains NaN or Inf values - "
                             "remove invalid points before DRT analysis")
    except ValueError as e:
        logger.error(str(e))
        return DRTResult()

    # === Step 2: R_inf Estimation ===
    rinf_est = _estimate_r_inf(frequencies, Z, r_inf_preset=r_inf_preset)
    R_inf = rinf_est.R_inf
    inductance_used, inductance_note = _inductance_choice(frequencies, Z, inductance)

    # === Steps 3-5: Build matrices, select lambda, solve NNLS ===
    matrices, lambda_sel, nnls_result, tau_extend, tau_extend_note = \
        _solve_with_extension(frequencies, Z, R_inf, n_tau, weighting,
                              tau_extend_decades, inductance_used,
                              lambda_reg, auto_lambda)
    n_grid = len(matrices.tau)
    tau_window = matrices.tau_window

    if not nnls_result.success:
        return DRTResult(
            R_inf=R_inf,
            diagnostics=DRTDiagnostics(
                freq_min=f_min, freq_max=f_max,
                freq_range_ratio=freq_range_ratio,
                n_points=len(frequencies),
                n_tau=n_grid,
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
                n_peaks=0,
                weighting=weighting,
                tau_extend_decades=tau_extend,
                tau_extend_note=tau_extend_note,
                inductance_used=inductance_used,
                inductance_note=inductance_note
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
        n_data=len(frequencies), tau_window=tau_window
    )

    # The peaks the run reports: GMM components when GMM ran and succeeded,
    # otherwise the scipy maxima. Everything downstream counts this set.
    reported_peaks = peaks_result or scipy_peaks or []
    n_peaks = len(reported_peaks)

    # === Step 8: Reconstruction & Error ===
    L_series = nnls_result.L_series
    Z_reconstructed = _reconstruct(matrices, gamma_physical, L_series, R_inf)
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
    # Same criterion as 'auto': a forced L on a non-inductive top decade is
    # fitted to noise or model error, not to a measured inductance.
    if L_series > 0 and not _inductance_choice(frequencies, Z, 'auto')[0]:
        nnls_result.warnings.append(
            f"L = {L_series * 1e9:.3g} nH fitted to data without inductive points "
            f"in the top decade - likely absorbs model error"
        )

    # === Step 8b: Window-edge diagnostics ===
    # Pile-up at an end of the tau grid: NNLS has nowhere to put response
    # whose time constant lies outside the grid, so it heaps gamma up against
    # the boundary instead.
    edge_pile_up_fraction, edge_pile_up_end = _edge_pile_up(
        gamma_physical, matrices.d_ln_tau, R_pol_from_gamma
    )
    if edge_pile_up_fraction > DRT_EDGE_BIN_RPOL_FRACTION:
        nnls_result.warnings.append(
            f"{edge_pile_up_fraction*100:.0f}% of R_pol is heaped against the "
            f"{'fast' if edge_pile_up_end == 'low' else 'slow'} end of the tau "
            f"grid - likely response from beyond it (series inductance, "
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
            f"{DRT_PEAK_EDGE_DECADES} decade of the measured tau window edge "
            f"or past it - position and R_estimate are only partly supported "
            f"by data"
        )

    # Mass the extended grid placed past the window is extrapolated, not measured.
    R_pol_extrapolated_fraction = _extrapolated_fraction(
        matrices.tau, gamma_physical, matrices.d_ln_tau, R_pol_from_gamma,
        tau_window[1]
    )
    if R_pol_extrapolated_fraction > DRT_EXTRAPOLATED_RPOL_FRACTION:
        nnls_result.warnings.append(
            f"{R_pol_extrapolated_fraction*100:.0f}% of R_pol lies past the "
            f"measured window (tau > {tau_window[1]:.3g} s) - extrapolated "
            f"from its high-frequency flank, not measured"
        )

    # === Step 8c: Shape-quality diagnostics (F3) ===
    # Warn if the DRT is too sparse/spiky for peak-shape analysis, or if
    # auto-lambda hit the search-range edge (either end: GCV can pin at the
    # top on very noisy data and then supplies lambda). Advisory
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
            f"{lambda_sel.lambda_value:.2e}); the optimum may lie outside the "
            f"searched range - DRT shape may be unreliable"
        )
    # A corner below GCV (flagged by the hybrid search) is the unusual
    # direction: the two criteria genuinely disagree and GCV's lambda was
    # taken without the L-curve's support. A corner above GCV is the expected
    # NNLS effect (GCV underestimates lambda) and not warned about, although
    # it exceeds a decade on 5 of 18 two-ZARC synthetics with 0.1-2 % noise.
    if lambda_sel.corner_below_gcv and lambda_sel.lambda_gcv and lambda_sel.lambda_lcurve:
        decades = np.log10(lambda_sel.lambda_gcv / lambda_sel.lambda_lcurve)
        nnls_result.warnings.append(
            f"L-curve corner (lambda={lambda_sel.lambda_lcurve:.2e}) lies "
            f"{decades:.1f} decades below GCV (lambda="
            f"{lambda_sel.lambda_gcv:.2e}); the larger (GCV) was used"
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
            matrices, lambda_sel.lambda_value, reference_peaks, Z, R_inf
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
        probe_curves=probe_curves,
        tau_window=tau_window
    )

    # === Build diagnostics ===
    diagnostics = DRTDiagnostics(
        freq_min=f_min,
        freq_max=f_max,
        freq_range_ratio=freq_range_ratio,
        n_points=len(frequencies),
        n_tau=n_grid,
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
        stability=stability,
        weighting=weighting,
        tau_extend_decades=tau_extend,
        tau_extend_note=tau_extend_note,
        R_pol_extrapolated_fraction=R_pol_extrapolated_fraction,
        inductance_used=inductance_used,
        inductance_note=inductance_note
    )

    return DRTResult(
        tau=matrices.tau,
        gamma=gamma,
        gamma_original=gamma_original,
        figure=fig,
        peaks=peaks_result,
        bic_scores=bic_scores,
        R_inf=R_inf,
        L_series=L_series,
        R_pol=R_pol_from_gamma,
        lambda_used=lambda_sel.lambda_value,
        reconstruction_error=rel_error,
        diagnostics=diagnostics
    )
