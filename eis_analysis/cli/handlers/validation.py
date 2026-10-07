"""
Data validation handlers for the EIS CLI.

- run_kk_validation: Kramers-Kronig validation
- run_zhit_validation: Z-HIT validation
- report_thd: linearity from the THD the instrument recorded (Gamry THD option)
- apply_zhit_reconstruction: --fit-on, Z-HIT reconstruction as a data correction
- report_outliers: per-point suspicious-point report from both methods
- plot_validation: KK and Z-HIT figures, with the flagged points marked
"""

import argparse
import logging
from dataclasses import replace
from typing import Optional

from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import EISAnalysisError, LoadedData, save_figure
from ...validation import (
    kramers_kronig_validation,
    zhit_validation,
    find_outliers,
    thd_check,
    KKResult,
    OutlierReport,
    THDResult,
    ZHITResult,
)
from ...visualization import plot_kk_validation, plot_zhit_validation, plot_thd
from ...validation.kramers_kronig import (KK_MAX_FRACTION_ABOVE, KK_RESIDUAL_THRESHOLD,
                                          low_frequency_slope)
from ...validation.zhit import _quality_label
from ...validation.thd import THD_TO_Z_ERROR

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
        kk_series_c)

    Returns
    -------
    result : KKResult or None
        Validation result with per-point residuals, or None if KK was
        skipped or failed. The figure is drawn by plot_validation.
    """
    if args.no_kk:
        return None

    log_separator()
    logger.info("Kramers-Kronig validation")
    log_separator()

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
    logger.info(f"KK: M={result.M} (from M={result.M_lower}, chi^2 plateau), "
                f"mu={result.mu:.4f} "
                f"(Lin-KK stop, threshold {args.mu_threshold}), "
                f"extend_decades={result.extend_decades:.2f}")
    logger.info(f"  Mean |res_real|: {result.mean_residual_real:.2f}%")
    logger.info(f"  Mean |res_imag|: {result.mean_residual_imag:.2f}%")
    logger.info(f"  Pseudo chi^2: {result.pseudo_chisqr:.2e}")
    logger.info(f"  Estimated noise (upper bound): {result.noise_estimate:.2f}%")
    if result.capacitance is not None:
        logger.info(f"  Series C: {result.capacitance:.2e} F")
    for warning in result.warnings:
        logger.warning(warning)

    mean_abs_residual = max(result.mean_residual_real, result.mean_residual_imag)
    # The label is graded on the mean, the verdict on points above the line;
    # the verdict wins. A local violation fails with an "excellent" mean, and
    # a few wild edge points pass with a mean past "poor", which means invalid.
    label = _quality_label(mean_abs_residual)
    if not result.is_valid:
        label = "poor"
    elif label == "poor":
        label = "marginal (check for drift/nonlinearity)"
    log_fn = logger.info if result.is_valid else logger.warning
    log_fn(f"Data quality: {label} "
           f"({result.n_above_threshold}/{len(frequencies)} points above "
           f"{KK_RESIDUAL_THRESHOLD}%, allowed {KK_MAX_FRACTION_ABOVE:.0%}; "
           f"max mean |res|={mean_abs_residual:.2f}%)")

    # A capacitive low-frequency end (-Z'' still rising) is what the series C
    # fixes: blocking/2-electrode cells, Warburg diffusion, an arc continuing
    # past f_min. The test is the shape of the data, not of the residuals: a
    # Warburg fails on a few edge points with a 1.5 % mean, and on a closing
    # end the series C would only absorb drift (doc/KK_INTUITION.md).
    if not result.is_valid and not args.kk_series_c:
        slope = low_frequency_slope(frequencies, Z)
        if slope < 0:
            logger.info(f"Hint: the low-frequency end is capacitive (-Z'' still "
                        f"rising, slope {slope:.2f} per decade) - try --kk-series-c. "
                        f"If the test still fails with it, the series C was not the cause.")

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
        CLI arguments (uses: no_zhit)

    Returns
    -------
    result : ZHITResult or None
        Validation result with per-point residuals, or None if Z-HIT was
        skipped. The figure is drawn by plot_validation.
    """
    if args.no_zhit:
        return None

    log_separator()
    logger.info("Z-HIT validation")
    log_separator()

    result = zhit_validation(frequencies, Z)
    if not result.success:
        # zhit_validation already logged why; the empty result still goes back
        # so the outlier report can skip Z-HIT rather than mistake it for data.
        return result

    # Summary (format consistent with KK validation)
    logger.info("Z-HIT: second order, offset = median over the spectrum")
    logger.info(f"  Mean |res_real|: {result.mean_residual_real:.2f}%")
    logger.info(f"  Mean |res_imag|: {result.mean_residual_imag:.2f}%")
    logger.info(f"  Pseudo chi^2: {result.pseudo_chisqr:.2e}")
    logger.info(f"  Estimated noise (upper bound): {result.noise_estimate:.2f}%")

    log_fn = logger.info if result.is_valid else logger.warning
    log_fn(f"Data quality: {result.quality_label} "
           f"(mean |res_mag|={result.mean_residual_mag:.2f}%, "
           f"threshold={result.quality_threshold:.1f}%)")

    return result


# =============================================================================
# THD recorded by the instrument
# =============================================================================

def report_thd(data: LoadedData, args: argparse.Namespace) -> Optional[THDResult]:
    """
    Summarize the per-point THD the instrument recorded and plot it against
    frequency.

    Silent when the data carries no THD (CSV, Gamry without the THD option).

    Parameters
    ----------
    data : LoadedData
        Full spectrum, with `current_thd` / `voltage_thd` from the loader
    args : argparse.Namespace
        CLI arguments (uses: save, format)

    Returns
    -------
    THDResult or None
        None when there is no THD to report
    """
    result = thd_check(data.frequencies, data.current_thd, data.voltage_thd)
    if result is None:
        return None

    log_separator()
    logger.info("THD (Gamry harmonic analysis)")
    log_separator()

    for warning in result.warnings:
        logger.warning(warning)

    limit = result.threshold * 100
    present = [(label, ch) for label, ch in (('Current', result.current), ('Voltage', result.voltage))
               if ch is not None]
    if not present:
        return result  # nothing was measured: no verdict, no empty figure

    # Channel-neutral cause: which channel is the sample's response depends on
    # the control mode (potentiostatic: current; galvanostatic: voltage), the
    # other one measures the purity of the excitation.
    for label, ch in present:
        logger.info(f"  {label}: median {ch.median * 100:.2f} %, "
                    f"max {ch.maximum * 100:.2f} % at {ch.f_at_max:.2e} Hz")
    for label, ch in present:
        if ch.n_above:
            logger.warning(f"{label} THD above {limit:.2g} % at {ch.n_above}/{ch.n_valid} "
                           f"points ({ch.f_above_min:.2e} - {ch.f_above_max:.2e} Hz): "
                           f"nonlinear response, distorted excitation, or noise at low signal")
    if not any(ch.n_above for _, ch in present):
        logger.info(f"All points with a THD value below {limit:.2g} % (|Z| error from "
                    f"nonlinearity below ~{limit * THD_TO_Z_ERROR:.2g} %)")

    fig = plot_thd(data.frequencies, data.current_thd, data.voltage_thd, result.threshold)
    save_figure(fig, args.save, 'thd', args.format)
    return result


# =============================================================================
# Z-HIT reconstruction as a data correction (--fit-on)
# =============================================================================

def apply_zhit_reconstruction(
    data: LoadedData,
    zhit_result: Optional[ZHITResult],
    args: argparse.Namespace
) -> LoadedData:
    """
    Attach or substitute the Z-HIT reconstruction according to --fit-on.

    Z-HIT is not only a validator ("was the system stationary?") but also a
    correction: where the low-frequency modulus drifts during the measurement
    while the phase stays sound - a coating taking up water, say - |Z| can be
    reconstructed from the phase and the circuit fitted against that instead.

    Must be called on the FULL spectrum, before frequency filtering, because
    that is where the reconstruction was computed: Z-HIT integrates the phase
    over log-omega, so a truncated range is a different reconstruction.

    Parameters
    ----------
    data : LoadedData
        Loaded spectrum, unfiltered
    zhit_result : ZHITResult or None
        Result of run_zhit_validation
    args : argparse.Namespace
        CLI arguments (uses: fit_on)

    Returns
    -------
    LoadedData
        Unchanged for --fit-on original; with `Z_zhit` attached for `zhit`;
        with `Z` itself replaced for `all`.

    Raises
    ------
    EISAnalysisError
        If the reconstruction was asked for but is not available.
    """
    if args.fit_on == 'original':
        return data

    if zhit_result is None or not zhit_result.success:
        # Falling back to the original would quietly fit something other than
        # what was asked for, and the fit report has no way to say so.
        raise EISAnalysisError(
            f"--fit-on {args.fit_on} needs the Z-HIT reconstruction, "
            "but Z-HIT produced no result for this spectrum"
        )

    log_separator()
    logger.info(f"Z-HIT reconstruction (--fit-on {args.fit_on})")
    log_separator()

    # How far the data moved is the whole point of the switch; without it the
    # user cannot tell whether the correction did anything.
    logger.info(f"|Z| replaced by the reconstruction from the phase "
                f"(mean shift {zhit_result.mean_residual_mag:.2f}%, max "
                f"{abs(zhit_result.residuals_mag).max():.2f}%)")
    logger.info("Applied to: " + ("every stage below" if args.fit_on == 'all'
                                  else "the circuit fit only"))

    if args.fit_on == 'all':
        # Truncation error peaks at relaxations, not edges; R_inf/DRT read it
        # wherever one sits (doc/ZHIT_REVIEW.md 1.3).
        logger.warning("R_inf and the DRT now read the Z-HIT reconstruction, "
                       "which deviates from |Z| by up to ~3% around sharp "
                       "relaxations even on clean data")
        # Title flows into visualize_data, so the Nyquist/Bode plot says which
        # curve it is showing.
        return replace(data, Z=zhit_result.Z_fit,
                       title=f"{data.title} (Z-HIT)")

    return replace(data, Z_zhit=zhit_result.Z_fit)


# =============================================================================
# Per-point outlier report
# =============================================================================

def report_outliers(
    frequencies: NDArray,
    kk_result: Optional[KKResult],
    zhit_result: Optional[ZHITResult],
    args: argparse.Namespace
) -> OutlierReport:
    """
    Report individual points whose KK or Z-HIT residual exceeds the threshold.

    Silent when nothing is flagged. The report goes back to the caller, so
    plot_validation can mark the flagged points.

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

    Returns
    -------
    OutlierReport
        The flagged points, empty when there are none
    """
    report = find_outliers(
        frequencies, kk_result, zhit_result, max_residual=args.max_residual
    )

    if not report.skipped and not report.points:
        return report

    # Own section header: the table draws on BOTH validations (see the
    # `flagged by` column), so printing it bare right after the Z-HIT block
    # made it read as part of Z-HIT.
    log_separator()
    logger.info("Per-point residual check")
    log_separator()

    for method in report.skipped:
        logger.info(f"{method}: over half the points exceed {args.max_residual:.1f}% - "
                    f"the spectrum fails as a whole, per-point flagging skipped")

    if not report.points:
        return report

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

    return report


def plot_validation(
    frequencies: NDArray,
    Z: NDArray,
    kk_result: Optional[KKResult],
    zhit_result: Optional[ZHITResult],
    report: OutlierReport,
    args: argparse.Namespace
) -> None:
    """
    Draw and save the KK and Z-HIT figures, each with its own flagged points.

    Drawn after report_outliers, so the figure is written once with the
    markers. Only the points THIS method flagged are marked: a band on the KK
    panel at a frequency where the KK residual is 0.3% would contradict the
    table's `flagged by` column.

    Parameters
    ----------
    frequencies, Z : ndarray
        Spectrum the validations ran on
    kk_result, zhit_result : KKResult, ZHITResult or None
        Results of the validation handlers; failed ones are skipped
    report : OutlierReport
        Result of report_outliers
    args : argparse.Namespace
        CLI arguments (uses: save, format)
    """
    for result, plot, method, suffix in ((kk_result, plot_kk_validation, 'KK', 'kk'),
                                         (zhit_result, plot_zhit_validation, 'Z-HIT', 'zhit')):
        if result is None or not result.success:
            continue
        flagged = [p.frequency for p in report.points if method in p.methods.split('+')]
        fig = plot(frequencies, Z, result, flagged_frequencies=flagged)
        save_figure(fig, args.save, suffix, args.format)
