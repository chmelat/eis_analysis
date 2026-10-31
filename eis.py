#!/usr/bin/env python3
"""
EIS Analysis with DRT
=====================

CLI tool for electrochemical impedance spectroscopy analysis.

Version: Imported from eis_analysis.version (single source of truth)

Features:
- Kramers-Kronig validation (Lin-KK method)
- DRT analysis with automatic lambda selection (GCV + L-curve)
- GMM peak detection in DRT spectrum (--peak-method gmm)
- Equivalent circuit fitting with operator syntax
- Oxide layer thickness estimation
- R_inf from an R-L-(R|Q) fit over the top two decades (--rinf-plot to see it)
- Fitting on the Z-HIT reconstruction (--fit-on zhit)

Usage:
    eis                             # synthetic data demo
    eis data.DTA                    # Gamry file (auto-lambda default)
    eis data.DTA --lambda 1e-6      # manual regularization
    eis data.DTA --peak-method gmm  # GMM peak detection
    eis data.DTA --rinf-plot        # plot the R_inf window fit
    eis data.DTA --voigt-chain --analyze-oxide  # oxide analysis

    eis --help                      # help

This is a thin wrapper around the CLI modules in eis_analysis.cli.
All logic has been refactored into separate modules for maintainability.
"""

import logging
import os

import matplotlib.pyplot as plt

from eis_analysis import get_version_string, visualize_data, visualize_ocv
from eis_analysis.cli import (
    # Logging
    setup_logging,
    log_separator,
    # Parser
    parse_arguments,
    # Data handling
    load_eis_data,
    filter_by_frequency,
    frequency_mask,
    # Handlers
    run_kk_validation,
    kk_noise_for_fit,
    run_zhit_validation,
    apply_zhit_reconstruction,
    report_outliers,
    select_excluded,
    flagged_in_fit,
    report_thd,
    plot_validation,
    run_drt_analysis,
    run_voigt_analysis,
    run_circuit_fitting,
    run_oxide_analysis,
    run_local_exponent,
    # Utils
    EISAnalysisError,
    save_figure,
)

logger = logging.getLogger(__name__)


def _run_analysis(args) -> None:
    """
    Run the full analysis pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments
    """
    peak_method = args.peak_method

    log_separator()
    logger.info(f"EIS Analysis ({get_version_string()})")
    log_separator()

    # Load the full measured spectrum
    data = load_eis_data(args)
    # For the fit JSON: which impedance of a drift-corrected Gamry file was read
    args.drift_corrected = data.drift_corrected

    # Validation runs on the FULL spectrum, before any frequency filtering.
    # KK is an integral relation over all frequencies, so truncating the range
    # (e.g. via --f-min) corrupts the imaginary-part reconstruction and yields
    # spurious residuals. The --f-min/--f-max filter applies only to the
    # analysis stages below (R_inf, DRT, circuit fit).
    kk_result = run_kk_validation(data.frequencies, data.Z, args)
    zhit_result = run_zhit_validation(data.frequencies, data.Z, args)
    report_thd(data, args)

    report = report_outliers(data.frequencies, kk_result, zhit_result, args)
    plot_validation(data.frequencies, data.Z, kk_result, zhit_result, report, args)

    # --fit-on: the Z-HIT reconstruction as a drift correction. Attached here
    # because this is the last point where the full spectrum still exists -
    # the reconstruction was computed on it and stays aligned with it.
    data = apply_zhit_reconstruction(data, zhit_result, args)

    # --exclude-flagged: the points the per-point check lists leave the
    # analysis region like the range outside --f-min/--f-max. Kept on args
    # (documented in run_circuit_fitting) for the fit export, which records
    # them; threading them through every fitting path would cost more.
    excluded = select_excluded(data.frequencies, report, kk_result, zhit_result, args)
    args.excluded_frequencies = excluded
    # Without --exclude-flagged: the flagged points the fit keeps, for its note
    args.flagged_in_fit = flagged_in_fit(data.frequencies, report, args)

    # Filter to the analysis region. The Lin-KK noise ran on the full
    # spectrum; the fit is compared with it on the points it is fitted to.
    mask = frequency_mask(data.frequencies, args, excluded)
    data = filter_by_frequency(data, args, excluded)
    kk_noise = kk_noise_for_fit(kk_result, mask, args.fit_on)
    kk_note = None
    if kk_result is not None and kk_result.capacitive_end:
        # The noise of the rerun with the series C (kk_noise_for_fit); the
        # note names it, or says why there is none
        kk_note = ("Lin-KK with the series C" if kk_result.with_series_c is not None else
                   "Lin-KK failed on a capacitive low-frequency end; rerun with --kk-series-c")

    # Visualization
    fig_vis = visualize_data(data.frequencies, data.Z, data.title)
    save_figure(fig_vis, args.save, 'nyquist_bode', args.format)

    # OCV visualization (for DTA files with OCV data)
    if args.ocv and data.ocv_data is not None:
        fig_ocv = visualize_ocv(data.ocv_data, data.title)
        save_figure(fig_ocv, args.save, 'ocv', args.format)

    # Local CPE exponent n(f)
    run_local_exponent(data.frequencies, data.Z, args)

    # DRT analysis
    drt_result = run_drt_analysis(data.frequencies, data.Z, args, peak_method)

    # Voigt element analysis
    run_voigt_analysis(drt_result, data.frequencies, data.Z, args)

    # Circuit fitting
    fitted_result, _ = run_circuit_fitting(data.frequencies, data.Z_for_fit, args,
                                           noise_percent=kk_noise, noise_note=kk_note)

    # Oxide layer analysis
    run_oxide_analysis(
        data.frequencies, data.Z, args, fitted_result, data.metadata
    )

    # Display plots
    if not args.no_show:
        plt.show()

    log_separator()
    logger.info("Analysis complete")
    log_separator()


def main():
    """Main CLI entry point."""
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args)

    # Every --save output goes under one prefix; a missing directory failed
    # each of them separately and the run still exited 0
    save_dir = os.path.dirname(args.save) if args.save else ''
    if save_dir:
        try:
            os.makedirs(save_dir, exist_ok=True)
        except OSError as e:
            logger.error(f"Cannot create the --save directory {save_dir}: {e}")
            exit(1)

    try:
        _run_analysis(args)
    except EISAnalysisError as e:
        logger.error(str(e))
        exit(1)
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        exit(130)


if __name__ == "__main__":
    main()
