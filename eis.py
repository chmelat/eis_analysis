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
- Robust R_inf estimation (--ri-fit)

Usage:
    eis                             # synthetic data demo
    eis data.DTA                    # Gamry file (auto-lambda default)
    eis data.DTA --lambda 0.5       # manual regularization
    eis data.DTA --peak-method gmm  # GMM peak detection
    eis data.DTA --ri-fit -v        # robust R_inf estimation
    eis data.DTA --voigt-chain --analyze-oxide  # oxide analysis

    eis --help                      # help

This is a thin wrapper around the CLI modules in eis_analysis.cli.
All logic has been refactored into separate modules for maintainability.
"""

import logging

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
    # Handlers
    run_kk_validation,
    run_zhit_validation,
    run_rinf_estimation,
    run_drt_analysis,
    run_voigt_analysis,
    run_circuit_fitting,
    run_oxide_analysis,
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
    # Auto-enable GMM if --classify-terms is used
    peak_method = args.peak_method
    if args.classify_terms and peak_method != 'gmm':
        logger.info("--classify-terms requires GMM peak detection")
        logger.info("Automatically enabling --peak-method gmm")
        peak_method = 'gmm'

    log_separator(60)
    logger.info(f"EIS Analysis ({get_version_string()})")
    log_separator(60)

    # Load and filter data
    data = load_eis_data(args)
    data = filter_by_frequency(data, args)

    # Visualization
    fig_vis = visualize_data(data.frequencies, data.Z, data.title)
    save_figure(fig_vis, args.save, 'nyquist_bode', args.format)

    # OCV visualization (for DTA files with OCV data)
    if args.ocv and data.ocv_data is not None:
        fig_ocv = visualize_ocv(data.ocv_data, data.title)
        save_figure(fig_ocv, args.save, 'ocv', args.format)

    # Kramers-Kronig validation
    run_kk_validation(data.frequencies, data.Z, args)

    # Z-HIT validation (optional)
    run_zhit_validation(data.frequencies, data.Z, args)

    # R_inf estimation
    R_inf_computed, _ = run_rinf_estimation(data.frequencies, data.Z, args)

    # DRT analysis
    drt_result = run_drt_analysis(
        data.frequencies, data.Z, args, R_inf_computed, peak_method
    )

    # Voigt element analysis
    run_voigt_analysis(drt_result, data.frequencies, data.Z, args)

    # Circuit fitting
    fitted_result, _ = run_circuit_fitting(data.frequencies, data.Z, args)

    # Oxide layer analysis
    run_oxide_analysis(
        data.frequencies, data.Z, args, fitted_result, data.metadata
    )

    # Display plots
    if not args.no_show:
        plt.show()

    log_separator(60)
    logger.info("Analysis complete")
    log_separator(60)


def main():
    """Main CLI entry point."""
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args)

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
