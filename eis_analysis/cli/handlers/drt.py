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

def _log_drt_diagnostics(result: DRTResult) -> None:
    """Log DRT analysis results from diagnostics."""
    diag = result.diagnostics
    if diag is None:
        return

    # R_inf estimation
    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    rinf = diag.rinf
    method_names = {
        'preset': 'Preset value',
        'median': 'Median of HF points',
        'rl_fit': 'R-L fit (auto-detection)',
        'voigt_fit': 'Voigt fit'
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
    if rinf.R_ct and rinf.C_nF:
        logger.info(f"  Voigt params: R_ct = {rinf.R_ct:.2f} Ohm, C = {rinf.C_nF:.2f} nF")
        if rinf.f_characteristic:
            logger.info(f"  Characteristic freq: f_char = {rinf.f_characteristic/1e6:.3f} MHz")

    if rinf.R_inf_median and rinf.method != 'median':
        diff_abs = rinf.R_inf - rinf.R_inf_median
        diff_pct = (diff_abs / rinf.R_inf_median * 100) if rinf.R_inf_median != 0 else 0
        logger.info(f"  Comparison: median = {rinf.R_inf_median:.3f} Ohm "
                    f"(diff: {diff_abs:+.3f} Ohm, {diff_pct:+.1f}%)")

    for warning in rinf.warnings:
        logger.warning(f"  {warning}")

    # DRT Analysis
    log_separator()
    logger.info("DRT Analysis")
    log_separator()
    logger.info(f"Using R_inf = {rinf.R_inf:.3f} Ohm")

    # Lambda selection
    lambda_sel = diag.lambda_sel
    lambda_method_names = {
        'user': 'User-specified',
        'default': 'Default',
        'gcv': 'GCV (automatic)',
        'hybrid': 'Hybrid GCV + L-curve',
        'fallback': 'Fallback (GCV failed)'
    }
    logger.info(f"Lambda: {lambda_method_names.get(lambda_sel.method, lambda_sel.method)}")
    if lambda_sel.method == 'hybrid' and lambda_sel.lambda_gcv:
        logger.info(f"  L-curve correction: lambda_gcv={lambda_sel.lambda_gcv:.2e} -> "
                    f"lambda={lambda_sel.lambda_value:.2e}")
    else:
        logger.info(f"  lambda = {lambda_sel.lambda_value:.2e}")
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
    logger.info(f"Found {diag.n_peaks} peaks")

    if diag.scipy_peaks:
        for i, peak in enumerate(diag.scipy_peaks):
            logger.info(f"  Peak {i+1}: tau = {peak['tau']:.2e} s "
                        f"(f = {peak['frequency']:.2e} Hz), R ~ {peak['R_estimate']:.2f} Ohm")

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
                       gmm_bic_threshold, save, format)
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
        use_rl_fit=False,
        use_voigt_fit=args.ri_fit,
        r_inf_preset=R_inf_computed,
        gmm_bic_threshold=args.gmm_bic_threshold
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

    try:
        voigt_info = analyze_voigt_elements(
            drt_result.tau, drt_result.gamma, frequencies, Z,
            peaks_gmm=drt_result.peaks
        )
        report = format_voigt_report(voigt_info)
        logger.info(report)

    except Exception as e:
        logger.warning(f"Voigt element analysis failed: {e}")
        logger.debug(f"Traceback: {e}", exc_info=True)
