"""
DRT analysis handlers for the EIS CLI.

- run_drt_analysis: DRT computation + diagnostics logging
- run_voigt_analysis: Voigt element analysis from the DRT spectrum
"""

import argparse
import logging
from math import isfinite
from typing import List, Optional

from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import save_figure
from ...drt import calculate_drt, DRTResult
from ...fitting import analyze_voigt_elements, VoigtSuggestion
from ...fitting.config import GMM_N_COMPONENTS_RANGE
from ...visualization import plot_drt

logger = logging.getLogger(__name__)


# =============================================================================
# DRT Analysis
# =============================================================================

def _log_rinf_estimation(rinf) -> None:
    """Report the R_inf estimate this DRT run made for itself."""
    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    logger.info("Method: Median of HF points (--ri-fit for an R-L-(R|Q) fit)")
    logger.info(f"R_inf = {rinf.R_inf:.3f} Ohm ({rinf.n_points_used} HF points)")


def _log_lambda_value(lambda_sel) -> None:
    """
    Report the selected lambda, and for the hybrid search both of its stages.

    The solver computes lambda mid-run, before any of these section headers
    exist, so it only records the numbers - reporting them is done here, where
    they belong to the DRT section. Run with -v for the search itself.
    """
    both_stages = lambda_sel.lambda_gcv and lambda_sel.lambda_lcurve
    if lambda_sel.method != 'hybrid' or not both_stages:
        logger.info(f"  lambda = {lambda_sel.lambda_value:.2e}")
        return

    ratio = lambda_sel.lambda_lcurve / lambda_sel.lambda_gcv
    # A corner more than a decade below GCV is in the DRT warnings, printed below.
    winner = 'L-curve' if lambda_sel.hybrid_stage == 'lcurve' else 'GCV'
    logger.info(f"  lambda = {lambda_sel.lambda_value:.2e}  "
                f"(GCV {lambda_sel.lambda_gcv:.2e}, L-curve corner "
                f"{lambda_sel.lambda_lcurve:.2e}, ratio {ratio:.2f}; larger: {winner})")


def _edge_marker(peak: dict) -> str:
    """Mark a peak the measured window leaves unsupported, or inflates."""
    marks = []
    if peak.get('outside_window'):
        marks.append(f"past window: {-peak['edge_distance_decades']:.2f} dec, extrapolated")
    elif peak.get('boundary_sensitive'):
        marks.append(f"edge: {peak['edge_distance_decades']:.2f} dec from window")
    if peak.get('edge_contaminated'):
        marks.append("R inflated by out-of-window pile-up")
    return f"  [{'; '.join(marks)}]" if marks else ""


def _log_gmm_selection(bic_scores: List[float], n_components: int) -> None:
    """
    Report which GMM model BIC settled on, and why it may not be the obvious one.

    `bic_scores` covers GMM_N_COMPONENTS_RANGE in order, and the chosen model is
    the number of peaks returned - so its score, the raw minimum and the edge
    test all read straight off the list. Run with -v for the search itself.
    """
    lo, hi = GMM_N_COMPONENTS_RANGE
    # A model that failed to fit is stored as inf, so it is no baseline to
    # measure an improvement against. The improvement is measured against the
    # smallest model, so quoting it for that model itself would just print zero.
    gain = bic_scores[0] - bic_scores[n_components - lo]
    gain_str = (f" (BIC improvement {gain:.1f} over {lo})"
                if n_components != lo and isfinite(gain) else "")
    logger.info(f"  BIC selection: {n_components} of {lo}-{hi} components{gain_str}")

    n_failed = sum(1 for bic in bic_scores if not isfinite(bic))
    if n_failed:
        logger.warning(f"  {n_failed} of {len(bic_scores)} candidate models failed "
                       f"to fit and were skipped (-v for the errors)")

    # Early stopping deliberately prefers the simpler model unless the extra
    # component earns more than gmm_bic_threshold - so the raw BIC minimum can
    # sit higher without that being an error.
    n_bic_min = bic_scores.index(min(bic_scores)) + lo
    if n_bic_min != n_components:
        logger.info(f"  Early stop (Occam): raw BIC minimum would be "
                    f"{n_bic_min} components, see --gmm-bic-threshold")

    # Only the upper bound is actionable: nothing lies below one component.
    if n_components == hi:
        logger.warning(f"  Optimum at the upper bound of the {lo}-{hi} range - the "
                       f"true component count may lie above it")


def _log_drt_diagnostics(result: DRTResult) -> None:
    """Log DRT analysis results from diagnostics."""
    diag = result.diagnostics
    if diag is None:
        return

    # R_inf estimation. A preset value did not come from this stage - it was
    # measured by --ri-fit, which already reported it in full, or handed in by
    # a caller of calculate_drt(). Repeating the section would restate someone
    # else's number under the heading of an estimation that never ran; the DRT
    # section states the value it uses either way.
    rinf = diag.rinf
    if rinf.method != 'preset':
        _log_rinf_estimation(rinf)

    # DRT Analysis
    log_separator()
    logger.info("DRT Analysis")
    log_separator()
    # The comparison is what places a caller-supplied R_inf against the data,
    # and the preset path carries no other diagnostic worth reporting.
    note = ""
    if rinf.method == 'preset' and rinf.R_inf_median:
        diff_pct = (rinf.R_inf - rinf.R_inf_median) / rinf.R_inf_median * 100
        note = f" (preset; HF median = {rinf.R_inf_median:.3f} Ohm, {diff_pct:+.1f}%)"
    logger.info(f"Using R_inf = {rinf.R_inf:.3f} Ohm{note}")
    logger.info(f"Weighting: {diag.weighting}")
    if diag.inductance_used:
        note = f" (auto: {diag.inductance_note})" if diag.inductance_note else ""
        logger.info(f"Series inductance: L = {result.L_series * 1e9:.1f} nH{note}")
    if diag.tau_extend_decades > 0 or diag.tau_extend_note:
        note = f" (auto: {diag.tau_extend_note})" if diag.tau_extend_note else ""
        logger.info(f"Tau grid extension: {diag.tau_extend_decades:.1f} decade past "
                    f"the measured window{note}")

    # Lambda selection
    lambda_sel = diag.lambda_sel
    lambda_method_names = {
        'user': 'User-specified',
        'default': 'Default',
        'gcv': 'GCV (L-curve correction failed)',
        'hybrid': 'Hybrid GCV + L-curve',
        'fallback': 'Fallback (GCV failed)'
    }
    logger.info(f"Lambda: {lambda_method_names.get(lambda_sel.method, lambda_sel.method)}")
    _log_lambda_value(lambda_sel)
    if diag.n_effective_bins is not None:
        logger.info(f"  DRT effective bins (N_eff): {diag.n_effective_bins:.1f}")

    # Matrix condition
    if diag.condition_number > 1e15:
        logger.warning(f"Matrix A is ill-conditioned ({diag.condition_number:.2e})")
    elif diag.condition_number > 1e12:
        logger.info(f"Matrix A has high condition number ({diag.condition_number:.2e})")

    # R_pol
    logger.info(f"R_pol (from data) = {diag.R_pol_from_data:.2f} Ohm")
    logger.info(f"R_pol (from DRT integral) = {diag.R_pol_from_gamma:.2f} Ohm")
    if diag.R_pol_extrapolated_fraction > 0:
        logger.info(f"  of which past the measured window: "
                    f"{diag.R_pol_extrapolated_fraction*100:.1f}%")
    if diag.normalized:
        logger.info("gamma(tau) normalized by R_pol")

    # Reconstruction error
    logger.info(f"Mean relative reconstruction error: {diag.reconstruction_error_rel:.1f}%")

    # NNLS warnings
    for warning in diag.nnls.warnings:
        logger.warning(f"  {warning}")

    # Peak detection
    log_separator()
    logger.info("Peak detection in DRT spectrum")
    log_separator()
    method_str = "GMM" if diag.peak_method == 'gmm' else "scipy.signal.find_peaks"
    logger.info(f"Method: {method_str}")
    if result.bic_scores:
        _log_gmm_selection(result.bic_scores, diag.n_peaks)
    logger.info(f"Found {diag.n_peaks} peaks")

    # Print the same peaks that n_peaks counts. With GMM, n_peaks reflects the
    # GMM components (result.peaks), which may differ from the raw scipy peaks
    # kept for diagnostics (GMM merges nearby maxima via BIC). Listing
    # scipy_peaks here would contradict the reported count.
    if diag.peak_method == 'gmm' and result.peaks:
        for i, peak in enumerate(result.peaks):
            logger.info(f"  Peak {i+1}: tau = {peak['tau_center']:.2e} s "
                        f"(f = {peak['f_center']:.2e} Hz), R ~ {peak['R_estimate']:.2f} Ohm, "
                        f"width = {peak['log_tau_std']:.2f} dec, "
                        f"weight = {peak['weight']:.3f}"
                        f"{_edge_marker(peak)}")
    elif diag.scipy_peaks:
        for i, peak in enumerate(diag.scipy_peaks):
            logger.info(f"  Peak {i+1}: tau = {peak['tau']:.2e} s "
                        f"(f = {peak['frequency']:.2e} Hz), R ~ {peak['R_estimate']:.2f} Ohm"
                        f"{_edge_marker(peak)}")

    # Lambda-probe peak stability
    if diag.stability is not None:
        stability = diag.stability
        log_separator()
        logger.info("Peak stability (lambda probe)")
        log_separator()
        logger.info(f"Reference lambda* = {stability.lambda_star:.2e}")
        logger.info(f"Probed span: {stability.span_decades:.2f} decades of lambda"
                    + (f" ({stability.n_clipped} probes clipped at the bounds)"
                       if stability.n_clipped else ""))

        for point in stability.probe_points:
            if point.success:
                logger.info(f"  lambda = {point.lambda_value:.2e}: "
                            f"{len(point.peaks)} peaks, "
                            f"gamma_max = {point.gamma_max:.3g} Ohm, "
                            f"reconstruction error = {point.reconstruction_error_rel:.1f}%")
            else:
                logger.warning(f"  lambda = {point.lambda_value:.2e}: "
                               f"solver failed ({point.error})")

        for i, peak_stab in enumerate(stability.peak_stability):
            line = (f"  Peak {i+1}: tau = {peak_stab.tau_ref:.2e} s  "
                    f"persistence {peak_stab.persistence}/{peak_stab.n_probes}  "
                    f"drift {peak_stab.tau_drift_decades:.2f} dec  "
                    f"R var {peak_stab.R_variation_rel*100:.0f}%  "
                    f"{peak_stab.verdict.upper()}")
            if peak_stab.verdict == 'artifact':
                logger.warning(line)
            else:
                logger.info(line)

        for warning in stability.warnings:
            logger.warning(f"  {warning}")

    log_separator()


def run_drt_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    R_inf_computed: Optional[float],
    peak_method: str
) -> DRTResult:
    """
    Run DRT analysis.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_drt, lambda_reg, n_tau, normalize_rpol,
                       gmm_bic_threshold, lambda_probe, drt_weighting, drt_inductance,
                       tau_extend, save, format)
    R_inf_computed : float or None
        Pre-computed R_inf from --ri-fit
    peak_method : str
        Peak detection method ('scipy' or 'gmm')

    Returns
    -------
    DRTResult
        Container with tau, gamma, peaks and diagnostics
    """
    if args.no_drt:
        return DRTResult()

    use_auto_lambda = args.lambda_reg is None

    result = calculate_drt(
        frequencies, Z,
        n_tau=args.n_tau,
        lambda_reg=args.lambda_reg,
        auto_lambda=use_auto_lambda,
        normalize_rpol=args.normalize_rpol,
        peak_method=peak_method,
        r_inf_preset=R_inf_computed,
        gmm_bic_threshold=args.gmm_bic_threshold,
        lambda_probe=args.lambda_probe,
        weighting=args.drt_weighting,
        tau_extend_decades=args.tau_extend,
        inductance=('auto' if args.drt_inductance == 'auto'
                    else args.drt_inductance == 'on')
    )

    # Log diagnostics
    _log_drt_diagnostics(result)

    # A plotting error must not take the computed DRT, and the stages
    # after it, down with it.
    if result.success:
        try:
            save_figure(plot_drt(Z, result), args.save, 'drt', args.format)
        except Exception as e:
            logger.warning(f"DRT figure failed: {e}")
            logger.debug(f"Traceback: {e}", exc_info=True)

    return result


# =============================================================================
# Voigt Element Analysis
# =============================================================================

def run_voigt_analysis(
    drt_result: DRTResult,
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> None:
    """
    Run Voigt element analysis from DRT results.

    Parameters
    ----------
    drt_result : DRTResult
        DRT analysis results
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_drt, no_voigt_info)
    """
    if args.no_drt or args.no_voigt_info:
        return
    if drt_result.tau is None or drt_result.gamma is None:
        return

    # With --normalize-rpol, drt_result.gamma is gamma/R_pol; the R and C
    # estimates need the unnormalized gamma [Ohm] (audit 2026-07-02, 2.2).
    gamma_ohm = (drt_result.gamma_original
                 if drt_result.gamma_original is not None
                 else drt_result.gamma)

    try:
        suggestion = analyze_voigt_elements(
            drt_result.tau, gamma_ohm, frequencies, Z,
            peaks_gmm=drt_result.peaks
        )
        _log_voigt_report(suggestion)

    except Exception as e:
        logger.warning(f"Voigt element analysis failed: {e}")
        logger.debug(f"Traceback: {e}", exc_info=True)


def _log_voigt_report(suggestion: VoigtSuggestion) -> None:
    """Report the Voigt elements read off the DRT, with how they were chosen."""
    log_separator()
    logger.info("Voigt elements (R||C) from DRT")
    log_separator()
    logger.info(f"Peaks: {suggestion.n_peaks_raw} detected ({suggestion.method}), "
                f"{suggestion.n_peaks_valid} valid")
    for note in suggestion.excluded_peaks:
        logger.warning(note)

    # Never empty: without DRT peaks there is one element from the -Z'' maximum
    elements = suggestion.elements
    logger.info("")
    logger.info("  ID | tau [s]    | f [Hz]     | R [Ohm]   | C [F]      | Warnings")
    logger.info("  " + "-" * 72)
    for elem in elements:
        warnings_str = ", ".join(elem.warnings) if elem.warnings else "-"
        if len(warnings_str) > 20:
            warnings_str = warnings_str[:17] + "..."
        logger.info(f"  {elem.id:2d} | "
                    f"{elem.tau:10.2e} | "
                    f"{elem.freq:10.2e} | "
                    f"{elem.R:9.1f} | "
                    f"{elem.C:10.2e} | "
                    f"{warnings_str}")
    logger.info("")

    # A ratio outside 0.5-2 is already one of suggestion.warnings
    logger.info("Consistency validation:")
    logger.info(f"  Sum R_i (from peaks): {suggestion.total_R:9.1f} Ohm")
    logger.info(f"  R_pol (from data):    {suggestion.R_pol:9.1f} Ohm")
    if suggestion.ratio == float('inf'):
        logger.info("  Ratio:                INF (R_pol = 0)")
    else:
        logger.info(f"  Ratio:                {suggestion.ratio:9.2f}")
    logger.info("")

    logger.info(f"Analysis quality: {suggestion.quality.upper()}")
    for warning in suggestion.warnings:
        logger.warning(warning)
    circuit = " - ".join(["R(R_inf)"] + [f"(R(R{i}) | C(C{i}))"
                                          for i in range(1, len(elements) + 1)])
    logger.info(f"Suggested circuit: {circuit}")
    log_separator()
