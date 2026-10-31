"""
Data loading and filtering for EIS CLI.

Contains:
- load_eis_data: Load from file or generate synthetic data
- frequency_mask: The --f-min/--f-max selection as a boolean mask
- filter_by_frequency: Apply frequency range filter
"""

import argparse
import logging
import os
from typing import Any, Dict, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from .logging import log_separator
from .utils import EISAnalysisError, LoadedData
from ..io import (
    LoadResult,
    load_data,
    load_csv_data,
    parse_ocv_curve,
    generate_synthetic_data,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Synthetic Data Configuration
# =============================================================================
# Default parameters for synthetic data demo
# These values represent a typical two-RC circuit with CPE behavior

SYNTHETIC_DATA_PARAMS: Dict[str, Any] = {
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

def _log_synthetic_params() -> None:
    """
    Announce the circuit the demo data is generated from.

    Reported here rather than by generate_synthetic_data(), which is a library
    function that should return data and print nothing - and SYNTHETIC_DATA_PARAMS
    is this caller's own choice, so this is where it is known.
    """
    Y0_0, n0 = SYNTHETIC_DATA_PARAMS['Q0']
    Y0_1, n1 = SYNTHETIC_DATA_PARAMS['Q1']

    log_separator()
    logger.info("Generating synthetic data")
    log_separator()
    logger.info("Circuit: Rs - (R0||Q0) - (R1||Q1)")
    logger.info(f"Rs = {SYNTHETIC_DATA_PARAMS['Rs']} Ω")
    logger.info(f"R0 = {SYNTHETIC_DATA_PARAMS['R0']:.2e} Ω, Q0 = ({Y0_0:.2e} S·s^n, n={n0})")
    logger.info(f"R1 = {SYNTHETIC_DATA_PARAMS['R1']:.2e} Ω, Q1 = ({Y0_1:.2e} S·s^n, n={n1})")


def print_metadata(metadata: Dict[str, Any]) -> None:
    """
    Print DTA file metadata as a CLI section.

    Parameters
    ----------
    metadata : dict
        Metadata dictionary from parse_dta_metadata()
    """
    log_separator()
    logger.info("DTA file metadata")
    log_separator()

    # Sample identification
    if metadata.get('title'):
        logger.info(f"Sample: {metadata['title']}")
    if metadata.get('date') or metadata.get('time'):
        date_str = metadata.get('date', '?')
        time_str = metadata.get('time', '?')
        logger.info(f"Measurement date: {date_str} {time_str}")

    # Notes
    if metadata.get('notes'):
        logger.info("Notes:")
        for note in metadata['notes']:
            logger.info(f"  - {note}")

    # EIS parameters
    logger.info("")
    logger.info("Measurement parameters:")

    if metadata.get('area') is not None:
        logger.info(f"  Sample area: {metadata['area']:.4f} cm²")

    if metadata.get('vdc') is not None:
        logger.info(f"  DC voltage: {metadata['vdc']:.4f} V")

    if metadata.get('vac') is not None:
        logger.info(f"  AC voltage: {metadata['vac']:.2f} mV rms")

    if metadata.get('freq_init') is not None and metadata.get('freq_final') is not None:
        logger.info(f"  Frequency range: {metadata['freq_final']:.2e} - {metadata['freq_init']:.2e} Hz")

    if metadata.get('pts_per_dec') is not None:
        logger.info(f"  Points per decade: {metadata['pts_per_dec']:.0f}")

    if metadata.get('pstat'):
        logger.info(f"  Potentiostat: {metadata['pstat']}")

    log_separator()


def _print_load_summary(result: LoadResult) -> None:
    """Print the caveats and the one-line summary for a loaded spectrum."""
    for warning in result.warnings:
        logger.warning(warning)
    f = result.frequencies
    logger.info(f"Loaded {len(f)} points from {result.filename}")
    if result.drift_corrected:
        logger.info("Drift-corrected impedance (Gamry ZrealDrCor/ZimagDrCor); "
                    "--raw-impedance reads the measured Zreal/Zimag")
    elif result.drift_corrected is not None:
        logger.info("Measured Zreal/Zimag; the file also holds Gamry's "
                    "drift-corrected impedance")
    logger.info(f"Frequency range: {f.min():.2e} - {f.max():.2e} Hz")


def load_eis_data(args: argparse.Namespace) -> LoadedData:
    """
    Load EIS data from file or generate synthetic data.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments. Uses:
        - args.input: Input file path (None for synthetic)
        - args.raw_impedance (optional, default False): for a Gamry file
          with drift correction, the measured instead of the corrected Z

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
    current_thd = voltage_thd = None
    drift_corrected = None
    raw_impedance = getattr(args, 'raw_impedance', False)
    ignored_raw = "--raw-impedance applies to Gamry .DTA files only; ignored"

    if args.input is None:
        # Synthetic data
        _log_synthetic_params()
        frequencies, Z = generate_synthetic_data(**SYNTHETIC_DATA_PARAMS)
        title = "Synthetic data"
        if raw_impedance:
            logger.warning(ignored_raw)
    else:
        # File input
        if not os.path.exists(args.input):
            raise EISAnalysisError(f"File '{args.input}' does not exist!")

        ext = os.path.splitext(args.input)[1].lower()
        if raw_impedance and ext != '.dta':
            logger.warning(ignored_raw)
        try:
            if ext == '.dta':
                loaded = load_data(args.input, drift_corrected=not raw_impedance)
                frequencies, Z = loaded.frequencies, loaded.Z
                metadata = loaded.metadata
                current_thd, voltage_thd = loaded.current_thd, loaded.voltage_thd
                drift_corrected = loaded.drift_corrected
                _print_load_summary(loaded)
                if raw_impedance and drift_corrected is None:
                    logger.info("The file has no drift-corrected impedance "
                                "(ZrealDrCor/ZimagDrCor): --raw-impedance has no effect")
                if metadata is not None:
                    print_metadata(metadata)
                # Load OCV data if available
                ocv_data = parse_ocv_curve(args.input)
                if ocv_data is not None:
                    logger.info(f"OCV data: {len(ocv_data['time'])} points, "
                                f"duration {ocv_data['time'][-1]/60:.1f} min")
                    vf = ocv_data['Vf']
                    ocv_mV = vf[-1] * 1000
                    mean_mV = float(np.mean(vf)) * 1000
                    drift_mV = abs(vf[-1] - vf[0]) * 1000
                    logger.info(f"  OCV = {ocv_mV:.1f} mV "
                                f"(mean {mean_mV:.1f} mV, drift {drift_mV:.2f} mV)")
            elif ext == '.csv':
                loaded = load_csv_data(args.input)
                frequencies, Z = loaded.frequencies, loaded.Z
                _print_load_summary(loaded)
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
        ocv_data=ocv_data,
        current_thd=current_thd,
        voltage_thd=voltage_thd,
        drift_corrected=drift_corrected
    )


# =============================================================================
# Data Filtering
# =============================================================================

def frequency_mask(frequencies: NDArray[np.float64],
                   args: argparse.Namespace,
                   excluded: Optional[Sequence[float]] = None) -> NDArray[np.bool_]:
    """
    The --f-min/--f-max selection, less the `excluded` points
    (--exclude-flagged), as a mask over `frequencies`.

    The one definition of the analysis range: filter_by_frequency applies it,
    and the pipeline uses it to restrict full-spectrum results (Lin-KK noise)
    to the same points.
    """
    mask: NDArray[np.bool_] = np.ones(len(frequencies), dtype=bool)
    if args.f_min is not None:
        mask &= (frequencies >= args.f_min)
    if args.f_max is not None:
        mask &= (frequencies <= args.f_max)
    if excluded:
        # The excluded frequencies are copied from this very array
        # (OutlierReport), so exact equality is safe
        mask &= ~np.isin(frequencies, excluded)
    return mask


def filter_by_frequency(
    data: LoadedData,
    args: argparse.Namespace,
    excluded: Optional[Sequence[float]] = None
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
    excluded : sequence of float, optional
        Frequencies [Hz] left out (--exclude-flagged), from select_excluded,
        which also logs them. Every point at such a frequency goes - a
        spectrum measuring one frequency twice loses both copies.

    Returns
    -------
    LoadedData
        Filtered data (or original if no filtering needed)

    Raises
    ------
    EISAnalysisError
        If no data remains after filtering
    """
    if args.f_min is None and args.f_max is None and not excluded:
        return data

    # Own section header so the filter output is not visually attached to the
    # preceding Z-HIT validation block. The filter applies to all analysis
    # stages below (visualization, R_inf, DRT, circuit fit), not to validation.
    log_separator()
    logger.info("Frequency filtering (analysis range)")
    log_separator()

    frequencies = data.frequencies
    Z = data.Z
    original_count = len(frequencies)

    mask = frequency_mask(frequencies, args, excluded)
    if args.f_min is not None:
        logger.info(f"Applying f_min = {args.f_min} Hz")
    if args.f_max is not None:
        logger.info(f"Applying f_max = {args.f_max} Hz")

    frequencies = frequencies[mask]
    Z = Z[mask]
    # The Z-HIT reconstruction was computed on the full spectrum and is aligned
    # element-wise with Z, so it has to survive the same mask - otherwise the
    # fit would pair reconstructed points with the wrong frequencies.
    Z_zhit = data.Z_zhit[mask] if data.Z_zhit is not None else None
    current_thd = data.current_thd[mask] if data.current_thd is not None else None
    voltage_thd = data.voltage_thd[mask] if data.voltage_thd is not None else None

    filtered_count = len(frequencies)
    removed_count = original_count - filtered_count
    logger.info(f"Frequency filter: {original_count} -> {filtered_count} points "
                f"(removed: {removed_count})")

    if filtered_count == 0:
        raise EISAnalysisError(
            "No data remaining after filtering! Check --f-min, --f-max and "
            "--exclude-flagged."
        )

    logger.info(f"Analysis frequency range: {frequencies.min():.2e} - "
                f"{frequencies.max():.2e} Hz")

    return LoadedData(
        frequencies=frequencies,
        Z=Z,
        title=data.title,
        metadata=data.metadata,
        ocv_data=data.ocv_data,
        Z_zhit=Z_zhit,
        current_thd=current_thd,
        voltage_thd=voltage_thd,
        drift_corrected=data.drift_corrected
    )
