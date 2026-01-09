"""
Data loading and filtering for EIS CLI.

Contains:
- load_eis_data: Load from file or generate synthetic data
- filter_by_frequency: Apply frequency range filter
"""

import argparse
import logging
import os

import numpy as np

from .utils import EISAnalysisError, LoadedData
from ..io import (
    load_data,
    load_csv_data,
    parse_dta_metadata,
    parse_ocv_curve,
    log_metadata,
    generate_synthetic_data,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Synthetic Data Configuration
# =============================================================================
# Default parameters for synthetic data demo
# These values represent a typical two-RC circuit with CPE behavior

SYNTHETIC_DATA_PARAMS = {
    'Rs': 10,           # Series resistance [Ohm]
    'R0': 1e5,          # First RC parallel resistance [Ohm]
    'Q0': (1e-6, 0.6),  # First CPE: (Q [F*s^(n-1)], n [-])
                        # n=0.6: moderate distribution of relaxation times
    'R1': 8e5,          # Second RC parallel resistance [Ohm]
    'Q1': (3e-5, 0.43), # Second CPE: (Q, n)
                        # n=0.43: wider distribution, near Warburg behavior
    'noise': 0.01,      # 1% noise level for realistic data
}


# =============================================================================
# Data Loading
# =============================================================================

def load_eis_data(args: argparse.Namespace) -> LoadedData:
    """
    Load EIS data from file or generate synthetic data.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments. Uses:
        - args.input: Input file path (None for synthetic)

    Returns
    -------
    LoadedData
        Container with frequencies, Z, title, metadata, and ocv_data

    Raises
    ------
    EISAnalysisError
        If file does not exist or cannot be parsed
    """
    metadata = None
    ocv_data = None

    if args.input is None:
        # Synthetic data
        frequencies, Z = generate_synthetic_data(**SYNTHETIC_DATA_PARAMS)
        title = "Synthetic data"
    else:
        # File input
        if not os.path.exists(args.input):
            raise EISAnalysisError(f"File '{args.input}' does not exist!")

        ext = os.path.splitext(args.input)[1].lower()
        try:
            if ext == '.dta':
                frequencies, Z = load_data(args.input)
                metadata = parse_dta_metadata(args.input)
                log_metadata(metadata)
                # Load OCV data if available
                ocv_data = parse_ocv_curve(args.input)
                if ocv_data is not None:
                    logger.info(f"OCV data: {len(ocv_data['time'])} points, "
                                f"duration {ocv_data['time'][-1]/60:.1f} min")
            elif ext == '.csv':
                frequencies, Z = load_csv_data(args.input)
            else:
                raise EISAnalysisError(
                    f"Unsupported format '{ext}'. Supported: .DTA (Gamry), .csv"
                )
        except EISAnalysisError:
            raise
        except Exception as e:
            raise EISAnalysisError(f"Error loading file: {e}") from e

        title = os.path.basename(args.input)

    return LoadedData(
        frequencies=frequencies,
        Z=Z,
        title=title,
        metadata=metadata,
        ocv_data=ocv_data
    )


# =============================================================================
# Data Filtering
# =============================================================================

def filter_by_frequency(
    data: LoadedData,
    args: argparse.Namespace
) -> LoadedData:
    """
    Filter data by frequency range.

    Parameters
    ----------
    data : LoadedData
        Loaded EIS data
    args : argparse.Namespace
        Parsed command line arguments. Uses:
        - args.f_min: Minimum frequency [Hz] or None
        - args.f_max: Maximum frequency [Hz] or None

    Returns
    -------
    LoadedData
        Filtered data (or original if no filtering needed)

    Raises
    ------
    EISAnalysisError
        If no data remains after filtering
    """
    if args.f_min is None and args.f_max is None:
        return data

    frequencies = data.frequencies
    Z = data.Z
    original_count = len(frequencies)

    mask = np.ones(len(frequencies), dtype=bool)

    if args.f_min is not None:
        mask &= (frequencies >= args.f_min)
        logger.info(f"Applying f_min = {args.f_min} Hz")

    if args.f_max is not None:
        mask &= (frequencies <= args.f_max)
        logger.info(f"Applying f_max = {args.f_max} Hz")

    frequencies = frequencies[mask]
    Z = Z[mask]

    filtered_count = len(frequencies)
    removed_count = original_count - filtered_count
    logger.info(f"Frequency filter: {original_count} -> {filtered_count} points "
                f"(removed: {removed_count})")

    if filtered_count == 0:
        raise EISAnalysisError(
            "No data remaining after filtering! Check --f-min and --f-max."
        )

    return LoadedData(
        frequencies=frequencies,
        Z=Z,
        title=data.title,
        metadata=data.metadata,
        ocv_data=data.ocv_data
    )
