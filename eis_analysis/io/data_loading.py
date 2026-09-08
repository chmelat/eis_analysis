"""
Data loading functions for EIS analysis.

This module provides functions to load impedance data from various file formats.
Native Gamry DTA parser - no external dependencies.
"""

import math
import unicodedata
import numpy as np
import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, Any, List
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Validation constants
MIN_DATA_POINTS = 10  # Minimum number of data points for analysis
MIN_FREQUENCY_RANGE = 10  # Minimum ratio f_max/f_min


@dataclass
class LoadResult:
    """
    A spectrum as it came out of a file.

    Caveats about the data land in `warnings` rather than on the console:
    the loader has no idea whether it runs under the CLI, in a notebook or
    in a batch script, and the caveat qualifies the returned spectrum the
    way an uncertainty qualifies a measurement. Failures of the operation
    itself (unreadable file, missing section) still raise or log, since
    there is no result for them to qualify.

    Attributes
    ----------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ohm]
    filename : str
        Path the data was read from
    metadata : dict or None
        DTA header metadata; None for formats that carry none (CSV)
    warnings : list of str
        Caveats about the data, in the order they were found
    """
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    filename: str
    metadata: Optional[Dict[str, Any]] = None
    warnings: List[str] = field(default_factory=list)


def _read_dta_lines(filename: str) -> List[str]:
    """
    Read a Gamry .DTA file as ASCII text lines. Raises OSError like open().

    Notes
    -----
    Gamry writes the Windows code page of the machine that recorded the run,
    cp1250 on the Czech systems these files come from. Reading with
    ``errors='ignore'`` instead silently deleted the diacritics from operator
    notes, so the decode is explicit and lossless.

    Folding to ASCII is deliberate: it makes the result identical whichever
    encoding was right, so a wrong guess degrades to unaccented text instead
    of mojibake. No number contains a non-ASCII character, so measurements
    are untouched.
    """
    with open(filename, 'rb') as f:
        raw = f.read()

    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError:
        # cp1250 leaves 5 bytes undefined; 'replace' maps them to U+FFFD,
        # which the fold below drops. So this can never raise.
        text = raw.decode('cp1250', errors='replace')

    # NFKD splits an accented letter into a base letter plus a combining
    # accent, so dropping the non-ASCII remainder leaves the letter behind.
    return unicodedata.normalize('NFKD', text).encode('ascii', 'ignore').decode('ascii').splitlines()


def read_gamry_native(filename: str) -> LoadResult:
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
    LoadResult
        Spectrum plus any caveat about how it was read (an unnamed ZCURVE
        header falls back to the standard column order and says so).

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
        lines = _read_dta_lines(filename)
    except FileNotFoundError:
        raise ValueError(f"File not found: {filename}")
    except OSError as e:
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

    # An abort marker ahead of the sweep means the run stopped while settling,
    # so the ZCURVE table is written but empty. Falling through would report a
    # missing ZCURVE section and send the reader after the wrong thing.
    if end_line is not None and end_line < start_line:
        raise ValueError(f"Experiment in {filename} was aborted before the impedance "
                         f"sweep started; the ZCURVE section contains no data")

    # Follow the header row rather than assuming columns 3 to 5.
    warnings: List[str] = []
    header = lines[start_line + 1].split() if start_line + 1 < len(lines) else []
    try:
        col_freq, col_zreal, col_zimag = (header.index(n) for n in ('Freq', 'Zreal', 'Zimag'))
    except ValueError:
        warnings.append(f"ZCURVE header in {filename} does not name Freq/Zreal/Zimag "
                        f"({header or 'header row missing'}); assuming standard order")
        col_freq, col_zreal, col_zimag = 2, 3, 4

    # A row must be long enough to hold the rightmost column we actually read.
    min_columns = max(col_freq, col_zreal, col_zimag) + 1

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

        if len(parts) < min_columns:
            continue

        try:
            freq = float(parts[col_freq])
            zr = float(parts[col_zreal])
            zi = float(parts[col_zimag])

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

    return LoadResult(freq_array, Z, filename, warnings=warnings)


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
        lines = _read_dta_lines(filename)
    except OSError as e:
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

        # Stop if we hit another section; data rows open with the point index.
        if not line[0].isdigit():
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
    metadata: Dict[str, Any] = {
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
        lines = _read_dta_lines(filename)

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

    except OSError as e:
        # Only IO errors are swallowed here (open/readlines). Per-field guards
        # inside the loop handle malformed data; a genuine parsing bug now
        # surfaces instead of being silently turned into partial metadata.
        logger.warning(f"Error reading file {filename}: {e}")

    return metadata


def expected_points(metadata: Dict[str, Any]) -> Optional[int]:
    """
    Points the sweep should contain, from parse_dta_metadata() output.

    Gamry steps logarithmically from FREQINIT down to FREQFINAL at PTSPERDEC
    points per decade, so the requested sweep length follows from the header
    alone. A run that was stopped early comes up short, which makes this a
    cheap integrity check on any file handed to us. Returns None when the
    header lacks usable sweep parameters.

    Notes
    -----
    The result is exact only up to the endpoint: Gamry stops at the first
    point at or below FREQFINAL, so a complete sweep may overshoot by one.
    Only a shortfall indicates a truncated run.
    """
    f_init = metadata.get('freq_init')
    f_final = metadata.get('freq_final')
    per_decade = metadata.get('pts_per_dec')

    if not f_init or not f_final or not per_decade or f_init <= f_final:
        return None

    return round(math.log10(f_init / f_final) * per_decade) + 1


def load_data(filename: str) -> LoadResult:
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
    LoadResult
        Spectrum, DTA header metadata, and any caveat about the data.

    Raises
    ------
    ValueError
        If data is invalid (empty, NaN, negative frequencies)
    """
    result = read_gamry_native(filename)
    frequencies, Z = result.frequencies, result.Z

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
        result.warnings.append(
            f"Small frequency range: {freq_range:.1f}x "
            f"(recommended >{MIN_FREQUENCY_RANGE}x); DRT analysis may have poor resolution")

    # Edge case: duplicate frequencies
    if len(np.unique(frequencies)) != len(frequencies):
        result.warnings.append("Dataset contains duplicate frequencies")

    # Edge case: sweep stopped before reaching the requested final frequency.
    # Only a shortfall is reported - see expected_points() on the overshoot.
    # The metadata is kept on the result so that callers need not parse the
    # file a second time for it.
    result.metadata = parse_dta_metadata(filename)
    n_expected = expected_points(result.metadata)
    if n_expected is not None and len(frequencies) < n_expected:
        result.warnings.append(
            f"Sweep may be truncated: {len(frequencies)} points, header implies "
            f"{n_expected} (lowest measured {frequencies.min():.2e} Hz, "
            f"header FREQFINAL {result.metadata['freq_final']:.2e} Hz)")

    return result


def _detect_delimiter(header_line: str) -> str:
    """
    Auto-detect CSV delimiter from header line.

    Returns whichever of comma, tab, semicolon occurs most often. Comma is
    listed first so it wins the all-zero case (single-column header).
    """
    return max(',', '\t', ';', key=header_line.count)


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
) -> LoadResult:
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
    LoadResult
        Spectrum plus any caveat about the data (metadata is None; CSV
        carries no header of its own).

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
    warnings: List[str] = []
    if freq_col is None or zreal_col is None or zimag_col is None:
        warnings.append("Could not detect columns from headers, using positional (0, 1, 2)")
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
        warnings.append(f"Small frequency range: {freq_range:.1f}x "
                        f"(recommended >{MIN_FREQUENCY_RANGE}x)")

    return LoadResult(freq_array, Z, filename, warnings=warnings)
