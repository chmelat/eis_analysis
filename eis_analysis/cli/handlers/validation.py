"""
Data validation handlers for the EIS CLI.

- run_kk_validation: Kramers-Kronig validation
- run_zhit_validation: Z-HIT validation
- report_outliers: per-point suspicious-point report from both methods
"""

import argparse
import logging
from typing import Optional

from numpy.typing import NDArray

from ..utils import save_figure
from ...validation import (
    kramers_kronig_validation,
    zhit_validation,
    find_outliers,
    KKResult,
    ZHITResult,
)
from ...validation.zhit import _quality_label

logger = logging.getLogger(__name__)


# =============================================================================
# Kramers-Kronig Validation
# =============================================================================

def run_kk_validation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[KKResult]:
    """
    Run Kramers-Kronig validation.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_kk, mu_threshold, auto_extend, extend_decades_max,
        kk_series_c, save, format)

    Returns
    -------
    result : KKResult or None
        Validation result (including `figure` and per-point residuals),
        or None if KK was skipped or failed.
    """
    if args.no_kk:
        return None

    logger.info("=" * 60)
    logger.info("Kramers-Kronig validation")
    logger.info("=" * 60)

    result = kramers_kronig_validation(
        frequencies, Z,
        mu_threshold=args.mu_threshold,
        auto_extend_decades=args.auto_extend,
        extend_decades_range=(0.0, args.extend_decades_max),
        include_C=args.kk_series_c
    )
    if not result.success:
        logger.warning(f"KK validation failed: {result.error}")
        return None

    # Summary (format consistent with Z-HIT validation)
    logger.info(f"KK: M={result.M}, mu={result.mu:.4f} "
                f"(Lin-KK stop, threshold {args.mu_threshold}), "
                f"extend_decades={result.extend_decades:.2f}")
    logger.info(f"  Mean |res_real|: {result.mean_residual_real:.2f}%")
    logger.info(f"  Mean |res_imag|: {result.mean_residual_imag:.2f}%")
    logger.info(f"  Pseudo chi^2: {result.pseudo_chisqr:.2e}")
    logger.info(f"  Estimated noise (upper bound): {result.noise_estimate:.2f}%")
    if result.capacitance is not None:
        logger.info(f"  Series C: {result.capacitance:.2e} F")

    mean_abs_residual = max(result.mean_residual_real, result.mean_residual_imag)
    log_fn = logger.info if result.is_valid else logger.warning
    log_fn(f"Data quality: {_quality_label(mean_abs_residual)} "
           f"(max mean |res|={mean_abs_residual:.2f}%, threshold=5.0%)")

    # Signature of a missing series term (L or C): the real part fits well
    # while imaginary residuals dominate - typical for blocking/2-electrode
    # cells whose series capacitance the Voigt chain cannot represent.
    if (not result.is_valid and not args.kk_series_c
            and result.mean_residual_imag > 5.0
            and result.mean_residual_real < 1.0):
        logger.info("Hint: imag residuals dominate while the real fit is good - "
                    "try --kk-series-c (blocking/2-electrode behavior)")

    save_figure(result.figure, args.save, 'kk', args.format)
    return result


# =============================================================================
# Z-HIT Validation
# =============================================================================

def run_zhit_validation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[ZHITResult]:
    """
    Run Z-HIT validation.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_zhit, zhit_optimize_offset, save, format)

    Returns
    -------
    result : ZHITResult or None
        Validation result (including `figure` and per-point residuals),
        or None if Z-HIT was skipped.
    """
    if args.no_zhit:
        return None

    result = zhit_validation(
        frequencies, Z,
        optimize_offset=args.zhit_optimize_offset
    )

    save_figure(result.figure, args.save, 'zhit', args.format)
    return result


# =============================================================================
# Per-point outlier report
# =============================================================================

def report_outliers(
    frequencies: NDArray,
    kk_result: Optional[KKResult],
    zhit_result: Optional[ZHITResult],
    args: argparse.Namespace
) -> None:
    """
    Report individual points whose KK or Z-HIT residual exceeds the threshold.

    Silent when nothing is flagged. Flagged points are also marked in the
    residual panel of both validation figures.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz] the validations ran on
    kk_result : KKResult or None
        Result from run_kk_validation
    zhit_result : ZHITResult or None
        Result from run_zhit_validation
    args : argparse.Namespace
        CLI arguments (uses: max_residual)
    """
    report = find_outliers(
        frequencies, kk_result, zhit_result, max_residual=args.max_residual
    )

    for method in report.skipped:
        logger.info(f"{method}: mean residual exceeds {args.max_residual:.1f}% - "
                    f"the spectrum fails as a whole, per-point flagging skipped")

    if not report.points:
        return

    logger.warning(f"Suspicious points ({len(report.points)}, residual > "
                   f"{args.max_residual:.1f}%, see --max-residual):")
    logger.warning(f"  {'f [Hz]':>10}  {'KK [%]':>8}  {'Z-HIT [%]':>10}   flagged by")
    for p in report.points:
        kk = f"{p.residual_kk:8.2f}" if p.residual_kk is not None else f"{'-':>8}"
        zhit = f"{p.residual_zhit:10.2f}" if p.residual_zhit is not None else f"{'-':>10}"
        logger.warning(f"  {p.frequency:10.3e}  {kk}  {zhit}   {p.methods}")

    # Deviations in the lowest decade are usually sample drift (a real,
    # non-stationary measurement) rather than bad points, and deleting them
    # would hide the problem instead of fixing it.
    f_min = float(min(frequencies))
    if any(p.frequency < f_min * 10 for p in report.points):
        logger.warning("  Note: deviations at the lowest frequencies are often "
                       "sample drift, not bad points")

    _mark_outliers(kk_result, report, 'kk', args)
    _mark_outliers(zhit_result, report, 'zhit', args)


def _mark_outliers(result, report, suffix: str, args: argparse.Namespace) -> None:
    """
    Mark flagged frequencies in a validation figure's residual panel.

    Done here rather than inside kramers_kronig_validation / zhit_validation so
    the computation core stays independent of the CLI threshold. The figure was
    already written to disk by the validation handler, so with --save it is
    re-saved to pick up the markers; that only happens on spectra that actually
    have flagged points.
    """
    fig = getattr(result, 'figure', None)
    if fig is None or len(fig.axes) < 2:
        return

    # Both validation figures use a 1x2 layout with the residuals on the right.
    ax = fig.axes[1]
    for i, p in enumerate(report.points):
        ax.axvline(x=p.frequency, color='red', linestyle='-', alpha=0.25,
                   linewidth=3, zorder=0,
                   label='Flagged point' if i == 0 else None)
    ax.legend()
    save_figure(fig, args.save, suffix, args.format)
