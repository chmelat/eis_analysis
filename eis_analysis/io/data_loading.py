"""
Data loading functions for EIS analysis.

This module provides functions to load impedance data from various file formats.
Native Gamry DTA parser - no external dependencies.
"""

import numpy as np
import logging
from typing import Tuple, Dict, Optional, Any, List
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Validation constants
MIN_DATA_POINTS = 10  # Minimum number of data points for analysis
MIN_FREQUENCY_RANGE = 10  # Minimum ratio f_max/f_min


def read_gamry_native(filename: str) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Native parser for Gamry .DTA files.

    Parses EIS data from ZCURVE section of Gamry potentiostat output files.
    No external dependencies required.

    Parameters
    ----------
    filename : str
        Path to .DTA file

    Returns
    -------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ω]

    Raises
    ------
    ValueError
        If file cannot be parsed or contains no valid data

    Notes
    -----
    Gamry DTA format:
    - Tab/whitespace-separated values
    - European decimal format (comma as separator)
    - ZCURVE section contains EIS data
    - Columns: Pt, Time, Freq, Zreal, Zimag, ...
    - EXPERIMENTABORTED marks interrupted experiments
    """
    frequencies: List[float] = []
    z_real: List[float] = []
    z_imag: List[float] = []

    # Read entire file (needed to detect EXPERIMENTABORTED)
    try:
        # Try ISO-8859-1 first (Gamry default), fallback to utf-8
        try:
            with open(filename, 'r', encoding='ISO-8859-1') as f:
                lines = f.readlines()
        except UnicodeDecodeError:
            with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
    except FileNotFoundError:
        raise ValueError(f"File not found: {filename}")
    except Exception as e:
        raise ValueError(f"Error reading file {filename}: {e}")

    # Find ZCURVE section and optional EXPERIMENTABORTED
    start_line = None
    end_line = None

    for i, line in enumerate(lines):
        if 'ZCURVE' in line:
            start_line = i
        if 'EXPERIMENTABORTED' in line:
            end_line = i

    if start_line is None:
        raise ValueError(f"No ZCURVE section found in {filename}")

    # Extract data lines (skip ZCURVE header + column names + units = 3 lines)
    if end_line is not None:
        raw_data = lines[start_line + 3:end_line]
    else:
        raw_data = lines[start_line + 3:]

    # Parse data lines
    for line in raw_data:
        line = line.strip()
        if not line:
            continue

        # Convert European decimal format and split by whitespace
        line = line.replace(',', '.')
        parts = line.split()

        # Need at least 5 columns: Pt, Time, Freq, Zreal, Zimag
        if len(parts) < 5:
            continue

        try:
            freq = float(parts[2])
            zr = float(parts[3])
            zi = float(parts[4])

            # Validate values
            if freq > 0 and np.isfinite(freq) and np.isfinite(zr) and np.isfinite(zi):
                frequencies.append(freq)
                z_real.append(zr)
                z_imag.append(zi)
        except (ValueError, IndexError):
            # Skip malformed lines
            continue

    if len(frequencies) == 0:
        raise ValueError(f"No valid EIS data found in {filename}. "
                        f"Check if file contains ZCURVE section.")

    # Convert to numpy arrays
    freq_array = np.array(frequencies, dtype=np.float64)
    Z = np.array(z_real, dtype=np.float64) + 1j * np.array(z_imag, dtype=np.float64)

    logger.debug(f"Parsed {len(freq_array)} data points from {filename}")

    return freq_array, Z


def parse_ocv_curve(filename: str) -> Optional[Dict[str, NDArray]]:
    """
    Parse OCVCURVE (Open Circuit Voltage) data from Gamry .DTA file.

    Parameters
    ----------
    filename : str
        Path to .DTA file

    Returns
    -------
    dict or None
        Dictionary with keys 'time', 'Vf', 'Vm' (numpy arrays), or None if not found
        - time: time in seconds [s]
        - Vf: filtered voltage [V]
        - Vm: measured voltage [V]
    """
    try:
        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except Exception as e:
        logger.warning(f"Error reading file {filename}: {e}")
        return None

    # Find OCVCURVE section
    start_line = None
    for i, line in enumerate(lines):
        if line.startswith('OCVCURVE'):
            start_line = i
            break

    if start_line is None:
        logger.debug(f"No OCVCURVE section found in {filename}")
        return None

    # Parse number of points from header: OCVCURVE<tab>TABLE<tab>N
    try:
        header_parts = lines[start_line].strip().split('\t')
        n_points = int(header_parts[2]) if len(header_parts) > 2 else 0
    except (ValueError, IndexError):
        n_points = 0

    # Data starts after header + column names + units (3 lines)
    data_start = start_line + 3

    time_data: List[float] = []
    vf_data: List[float] = []
    vm_data: List[float] = []

    for i in range(data_start, min(data_start + n_points, len(lines))):
        line = lines[i].strip()
        if not line:
            continue

        # Stop if we hit another section
        if not line[0].isdigit() and not line.startswith('\t'):
            break

        # Convert European decimal format
        line = line.replace(',', '.')
        parts = line.split()

        # Columns: Pt, T, Vf, Vm, Ach, Over, Temp
        if len(parts) >= 4:
            try:
                t = float(parts[1])
                vf = float(parts[2])
                vm = float(parts[3])
                time_data.append(t)
                vf_data.append(vf)
                vm_data.append(vm)
            except (ValueError, IndexError):
                continue

    if len(time_data) == 0:
        return None

    return {
        'time': np.array(time_data, dtype=np.float64),
        'Vf': np.array(vf_data, dtype=np.float64),
        'Vm': np.array(vm_data, dtype=np.float64),
    }


def parse_dta_metadata(filename: str) -> Dict[str, Any]:
    """
    Parse metadata from Gamry .DTA file.

    Extracts useful information such as:
    - AREA: sample area [cm²]
    - VDC: DC voltage [V]
    - VAC: AC voltage [mV rms]
    - FREQINIT: initial frequency [Hz]
    - FREQFINAL: final frequency [Hz]
    - PTSPERDEC: points per decade
    - DATE/TIME: measurement date and time
    - NOTES: experiment notes
    - PSTAT: potentiostat model

    Parameters
    ----------
    filename : str
        Path to .DTA file

    Returns
    -------
    dict
        Dictionary with metadata (empty values if field is missing)
    """
    metadata = {
        'area': None,
        'vdc': None,
        'vac': None,
        'freq_init': None,
        'freq_final': None,
        'pts_per_dec': None,
        'date': None,
        'time': None,
        'notes': [],
        'pstat': None,
        'title': None,
    }

    try:
        with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i]
            parts = line.strip().split('\t')

            if len(parts) < 2:
                i += 1
                continue

            key = parts[0]

            if key == 'AREA' and parts[1] == 'QUANT':
                try:
                    metadata['area'] = float(parts[2].replace(',', '.'))
                except (ValueError, IndexError):
                    pass

            elif key == 'VDC' and parts[1] == 'POTEN':
                try:
                    metadata['vdc'] = float(parts[2].replace(',', '.'))
                except (ValueError, IndexError):
                    pass

            elif key == 'VAC' and parts[1] == 'QUANT':
                try:
                    metadata['vac'] = float(parts[2].replace(',', '.'))  # mV rms
                except (ValueError, IndexError):
                    pass

            elif key == 'FREQINIT' and parts[1] == 'QUANT':
                try:
                    metadata['freq_init'] = float(parts[2].replace(',', '.'))
                except (ValueError, IndexError):
                    pass

            elif key == 'FREQFINAL' and parts[1] == 'QUANT':
                try:
                    metadata['freq_final'] = float(parts[2].replace(',', '.'))
                except (ValueError, IndexError):
                    pass

            elif key == 'PTSPERDEC' and parts[1] == 'QUANT':
                try:
                    metadata['pts_per_dec'] = float(parts[2].replace(',', '.'))
                except (ValueError, IndexError):
                    pass

            elif key == 'DATE' and parts[1] == 'LABEL':
                metadata['date'] = parts[2] if len(parts) > 2 else None

            elif key == 'TIME' and parts[1] == 'LABEL':
                metadata['time'] = parts[2] if len(parts) > 2 else None

            elif key == 'NOTES' and parts[1] == 'NOTES':
                # NOTES format: NOTES<tab>NOTES<tab>N<tab>label
                # Followed by N lines starting with tab
                try:
                    n_lines = int(parts[2])
                    for j in range(n_lines):
                        if i + 1 + j < len(lines):
                            note_line = lines[i + 1 + j].strip()
                            if note_line:
                                metadata['notes'].append(note_line)
                    i += n_lines  # Skip the note lines
                except (ValueError, IndexError):
                    pass

            elif key == 'PSTAT' and parts[1] == 'PSTAT':
                metadata['pstat'] = parts[2] if len(parts) > 2 else None

            elif key == 'TITLE' and parts[1] == 'LABEL':
                metadata['title'] = parts[2] if len(parts) > 2 else None

            # Stop at data section
            elif key == 'ZCURVE':
                break

            i += 1

    except Exception as e:
        logger.warning(f"Error parsing metadata from {filename}: {e}")

    return metadata


def log_metadata(metadata: Dict[str, Any]) -> None:
    """
    Log metadata in a readable format.

    Parameters
    ----------
    metadata : dict
        Metadata dictionary from parse_dta_metadata()
    """
    logger.info("="*60)
    logger.info("DTA file metadata")
    logger.info("="*60)

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

    logger.info("="*60)


def load_data(filename: str) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Load data from Gamry .DTA file.

    Uses native parser (no dependency on impedance.py).
    For CSV files use load_csv_data().

    Parameters
    ----------
    filename : str
        Path to .DTA file

    Returns
    -------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ω]

    Raises
    ------
    ValueError
        If data is invalid (empty, NaN, negative frequencies)
    """
    frequencies, Z = read_gamry_native(filename)

    # Data validation
    if len(frequencies) == 0 or len(Z) == 0:
        raise ValueError("File contains no data")

    if len(frequencies) != len(Z):
        raise ValueError(f"Array length mismatch: {len(frequencies)} frequencies, {len(Z)} impedances")

    if np.any(frequencies <= 0):
        raise ValueError("Frequencies must be positive")

    if np.any(~np.isfinite(frequencies)) or np.any(~np.isfinite(Z)):
        raise ValueError("Data contains NaN or Inf values")

    # Edge case: minimum number of points
    if len(frequencies) < MIN_DATA_POINTS:
        raise ValueError(f"Dataset must have at least {MIN_DATA_POINTS} points, got {len(frequencies)}")

    # Edge case: frequency range
    freq_range = frequencies.max() / frequencies.min()
    if freq_range < MIN_FREQUENCY_RANGE:
        logger.warning(f"Small frequency range: {freq_range:.1f}x (recommended >{MIN_FREQUENCY_RANGE}x)")
        logger.warning("DRT analysis may have poor resolution")

    # Edge case: duplicate frequencies
    if len(np.unique(frequencies)) != len(frequencies):
        logger.warning("Dataset contains duplicate frequencies")

    logger.info(f"Loaded {len(frequencies)} points from {filename}")
    logger.info(f"Frequency range: {frequencies.min():.2e} - {frequencies.max():.2e} Hz")
    return frequencies, Z


def _detect_delimiter(header_line: str) -> str:
    """
    Auto-detect CSV delimiter from header line.

    Tries common delimiters: tab, semicolon, comma.
    Returns the one that produces most columns.
    """
    delimiters = ['\t', ';', ',']
    best_delimiter = ','
    max_columns = 0

    for delim in delimiters:
        columns = len(header_line.split(delim))
        if columns > max_columns:
            max_columns = columns
            best_delimiter = delim

    return best_delimiter


def _find_column_index(headers: List[str], patterns: List[str]) -> Optional[int]:
    """
    Find column index matching any of the patterns (case-insensitive).

    Parameters
    ----------
    headers : list of str
        Column header names
    patterns : list of str
        Patterns to match (case-insensitive)

    Returns
    -------
    int or None
        Column index if found, None otherwise
    """
    headers_lower = [h.lower().strip() for h in headers]

    for pattern in patterns:
        pattern_lower = pattern.lower()
        for i, header in enumerate(headers_lower):
            if pattern_lower in header or header in pattern_lower:
                return i

    return None


def load_csv_data(
    filename: str,
    delimiter: Optional[str] = None
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Load EIS data from CSV file with auto-detection of columns and delimiter.

    Automatically detects:
    - Delimiter: comma, semicolon, or tab
    - Columns: frequency, Z_real, Z_imag by header names
    - Comments: lines starting with '#' are ignored

    Supported column names (case-insensitive):
    - Frequency: 'freq', 'frequency', 'f', 'hz'
    - Z real: 'zreal', 'z_real', 'z\'', 're', 'real', 'z.real'
    - Z imag: 'zimag', 'z_imag', 'z\'\'', 'im', 'imag', 'z.imag'

    Parameters
    ----------
    filename : str
        Path to CSV file
    delimiter : str, optional
        Column delimiter. If None, auto-detected from header.

    Returns
    -------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ω]

    Raises
    ------
    ValueError
        If file cannot be parsed or required columns not found

    Examples
    --------
    Supported CSV formats:

    Format 1 (comma-separated):
        frequency,Z_real,Z_imag
        100000,10.5,-5.2

    Format 2 (semicolon, European):
        freq;Zreal;Zimag
        100000;10,5;-5,2

    Format 3 (tab-separated):
        f	Re(Z)	Im(Z)
        100000	10.5	-5.2

    Format 4 (with comments):
        # EIS data exported from instrument
        # Units: Hz, Ω, Ω
        frequency[Hz],Z_real[Ω],Z_imag[Ω]
        100000,10.5,-5.2
    """
    # Read file
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(filename, 'r', encoding='ISO-8859-1') as f:
            lines = f.readlines()

    # Skip comment lines (starting with #) to find header
    header_idx = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and not stripped.startswith('#'):
            header_idx = i
            break

    if len(lines) - header_idx < 2:
        raise ValueError(f"CSV file {filename} must have header and at least one data row")

    # Detect delimiter from header
    header_line = lines[header_idx].strip()
    if delimiter is None:
        delimiter = _detect_delimiter(header_line)

    logger.debug(f"CSV delimiter: '{repr(delimiter)}'")

    # Parse header
    headers = header_line.split(delimiter)
    logger.debug(f"CSV headers: {headers}")

    # Find column indices
    freq_patterns = ['freq', 'frequency', 'f', 'hz']
    zreal_patterns = ['zreal', 'z_real', "z'", 're(z)', 'real', 'z.real', 're']
    zimag_patterns = ['zimag', 'z_imag', "z''", 'im(z)', 'imag', 'z.imag', 'im']

    freq_col = _find_column_index(headers, freq_patterns)
    zreal_col = _find_column_index(headers, zreal_patterns)
    zimag_col = _find_column_index(headers, zimag_patterns)

    # Fallback to positional if headers not found
    if freq_col is None or zreal_col is None or zimag_col is None:
        logger.warning("Could not detect columns from headers, using positional (0, 1, 2)")
        freq_col, zreal_col, zimag_col = 0, 1, 2

    logger.debug(f"Column indices: freq={freq_col}, zreal={zreal_col}, zimag={zimag_col}")

    # Parse data rows
    frequencies: List[float] = []
    z_real: List[float] = []
    z_imag: List[float] = []

    for line_num, line in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Handle European decimal format (comma -> dot) for semicolon-delimited
        if delimiter == ';':
            line = line.replace(',', '.')

        parts = line.split(delimiter)

        try:
            freq = float(parts[freq_col].replace(',', '.'))
            zr = float(parts[zreal_col].replace(',', '.'))
            zi = float(parts[zimag_col].replace(',', '.'))

            if freq > 0 and np.isfinite(freq) and np.isfinite(zr) and np.isfinite(zi):
                frequencies.append(freq)
                z_real.append(zr)
                z_imag.append(zi)
        except (ValueError, IndexError) as e:
            logger.debug(f"Skipping line {line_num}: {e}")
            continue

    if len(frequencies) == 0:
        raise ValueError(f"No valid data found in {filename}")

    freq_array = np.array(frequencies, dtype=np.float64)
    Z = np.array(z_real, dtype=np.float64) + 1j * np.array(z_imag, dtype=np.float64)

    # Validation
    if len(freq_array) < MIN_DATA_POINTS:
        raise ValueError(f"Dataset must have at least {MIN_DATA_POINTS} points, got {len(freq_array)}")

    freq_range = freq_array.max() / freq_array.min()
    if freq_range < MIN_FREQUENCY_RANGE:
        logger.warning(f"Small frequency range: {freq_range:.1f}x (recommended >{MIN_FREQUENCY_RANGE}x)")

    logger.info(f"Loaded {len(freq_array)} points from {filename}")
    logger.info(f"Frequency range: {freq_array.min():.2e} - {freq_array.max():.2e} Hz")

    return freq_array, Z
