"""
DRT analysis handlers for the EIS CLI.

- run_drt_analysis: DRT computation + diagnostics logging
- run_voigt_analysis: Voigt element analysis from the DRT spectrum
"""

import argparse
import logging
from typing import Optional

from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import save_figure
from ...drt import calculate_drt, DRTResult
from ...fitting import analyze_voigt_elements, format_voigt_report

logger = logging.getLogger(__name__)


# =============================================================================
# DRT Analysis
# =============================================================================

def _log_rinf_estimation(rinf) -> None:
    """Report the R_inf estimate this DRT run made for itself."""
    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    method_names = {
        'median': 'Median of HF points',
        'rl_fit': 'R-L fit (auto-detection)'
    }
    logger.info(f"Method: {method_names.get(rinf.method, rinf.method)}")

    if rinf.behavior:
        logger.info(f"  Detected: {rinf.behavior.capitalize()} behavior")

    if rinf.n_points_used:
        logger.info(f"R_inf = {rinf.R_inf:.3f} Ohm ({rinf.n_points_used} HF points)")
    else:
        logger.info(f"R_inf = {rinf.R_inf:.3f} Ohm")

    if rinf.R_squared and rinf.R_squared > 0:
        logger.info(f"  Quality: R^2 = {rinf.R_squared:.4f}")
    if rinf.L_nH and rinf.L_nH > 0:
        logger.info(f"  Inductance: L = {rinf.L_nH:.2f} nH")

    if rinf.R_inf_median and rinf.method != 'median':
        diff_abs = rinf.R_inf - rinf.R_inf_median
        diff_pct = (diff_abs / rinf.R_inf_median * 100) if rinf.R_inf_median != 0 else 0
        logger.info(f"  Comparison: median = {rinf.R_inf_median:.3f} Ohm "
                    f"(diff: {diff_abs:+.3f} Ohm, {diff_pct:+.1f}%)")

    for warning in rinf.warnings:
        logger.warning(f"  {warning}")


def _preset_rinf_note(rinf) -> str:
    """
    How a preset R_inf compares with the HF median, appended to "Using R_inf".

    The preset path carries only these two numbers, and the comparison is the
    one worth keeping: it is what places a caller-supplied R_inf against the
    data. Empty string for a value this stage estimated itself.
    """
    if rinf.method != 'preset' or not rinf.R_inf_median:
        return ""

    diff_pct = (rinf.R_inf - rinf.R_inf_median) / rinf.R_inf_median * 100
    return (f" (preset; HF median = {rinf.R_inf_median:.3f} Ohm, "
            f"{diff_pct:+.1f}%)")


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
    logger.info(f"  lambda = {lambda_sel.lambda_value:.2e}  "
                f"(GCV {lambda_sel.lambda_gcv:.2e} -> L-curve corner "
                f"{lambda_sel.lambda_lcurve:.2e}, ratio {ratio:.2f})")

    # A corner within a decade of the GCV guess is the expected outcome; the
    # other two stages mean the two criteria disagreed and are worth flagging.
    if lambda_sel.hybrid_stage == 'lcurve_correction':
        logger.warning(f"  L-curve raised GCV's lambda (ratio {ratio:.1f}) - GCV "
                       f"underestimates lambda under the NNLS non-negativity constraint")
    elif lambda_sel.hybrid_stage == 'geometric_mean':
        logger.warning(f"  GCV and L-curve disagree (ratio {ratio:.2f}); "
                       f"geometric mean used")


def _log_gmm_selection(gmm) -> None:
    """
    Report which GMM model BIC settled on, and why it may not be the obvious one.

    Like the lambda search, model selection runs inside calculate_drt(), before
    this section exists - so the solver only records the outcome. Run with -v
    for the per-model BIC scores.
    """
    lo, hi = gmm.n_components_range
    # The improvement is measured against the smallest model, so quoting it for
    # the smallest model itself would just print zero.
    gain = (f" (BIC improvement {gmm.bic_improvement:.1f} over {lo})"
            if gmm.n_components != lo else "")
    logger.info(f"  BIC selection: {gmm.n_components} of {lo}-{hi} components{gain}")

    # Early stopping deliberately prefers the simpler model unless the extra
    # component earns more than gmm_bic_threshold - so the raw BIC minimum can
    # sit higher without that being an error.
    if gmm.n_components_bic_min != gmm.n_components:
        logger.info(f"  Early stop (Occam): raw BIC minimum would be "
                    f"{gmm.n_components_bic_min} components, see --gmm-bic-threshold")

    if gmm.at_range_edge:
        logger.warning(f"  Optimum at the edge of the {lo}-{hi} range - the true "
                       f"component count may lie outside it")


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
    logger.info(f"Using R_inf = {rinf.R_inf:.3f} Ohm{_preset_rinf_note(rinf)}")


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
    if diag.gmm is not None:
        _log_gmm_selection(diag.gmm)
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
                        f"weight = {peak['weight']:.3f}")
    elif diag.scipy_peaks:
        for i, peak in enumerate(diag.scipy_peaks):
            logger.info(f"  Peak {i+1}: tau = {peak['tau']:.2e} s "
                        f"(f = {peak['frequency']:.2e} Hz), R ~ {peak['R_estimate']:.2f} Ohm")

    # Lambda-probe peak stability
    if diag.stability is not None:
        stability = diag.stability
        log_separator()
        logger.info("Peak stability (lambda probe)")
        log_separator()
        logger.info(f"Reference lambda* = {stability.lambda_star:.2e}")

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
        CLI arguments (uses: no_drt, lambda_reg, n_tau, normalize_rpol, ri_fit,
                       gmm_bic_threshold, lambda_probe, save, format)
    R_inf_computed : float or None
        Pre-computed R_inf from --ri-fit
    peak_method : str
        Peak detection method ('scipy' or 'gmm')

    Returns
    -------
    DRTResult
        Container with tau, gamma, peaks, and figures
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
        use_rl_fit=args.ri_fit,
        r_inf_preset=R_inf_computed,
        gmm_bic_threshold=args.gmm_bic_threshold,
        lambda_probe=args.lambda_probe
    )

    # Log diagnostics
    _log_drt_diagnostics(result)

    save_figure(result.figure, args.save, 'drt', args.format)

    # Save R_inf figure from DRT only if not already saved via --ri-fit
    if not args.ri_fit:
        save_figure(result.figure_rinf, args.save, 'ri_fit', args.format)

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
        voigt_info = analyze_voigt_elements(
            drt_result.tau, gamma_ohm, frequencies, Z,
            peaks_gmm=drt_result.peaks
        )
        report = format_voigt_report(voigt_info)
        logger.info(report)

    except Exception as e:
        logger.warning(f"Voigt element analysis failed: {e}")
        logger.debug(f"Traceback: {e}", exc_info=True)
