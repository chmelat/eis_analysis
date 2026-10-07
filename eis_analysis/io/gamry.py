"""
Loading EIS data from Gamry .DTA files: the ZCURVE sweep, the OCVCURVE
and the header metadata. Native parser - no external dependencies.
"""

import logging
import math
import unicodedata
from typing import Any, Dict, List, Optional

import numpy as np
from numpy.typing import NDArray

from .spectrum import LoadResult, _check_spectrum, _drop_negative_real_hf

logger = logging.getLogger(__name__)


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
    - Tab-separated values, rows open with a tab
    - European decimal format (comma as separator)
    - ZCURVE section contains EIS data
    - Columns: Pt, Time, Freq, Zreal, Zimag, ...
    - EXPERIMENTABORTED marks interrupted experiments
    """
    frequencies: List[float] = []
    z_real: List[float] = []
    z_imag: List[float] = []
    thd: Dict[str, List[float]] = {'Ithd': [], 'Vthd': []}

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

    # Follow the header row rather than assuming columns. Tab split keeps empty
    # cells as fields; rows open with a tab, so standard order starts at 3.
    warnings: List[str] = []
    header = lines[start_line + 1].split('\t') if start_line + 1 < len(lines) else []
    try:
        col_freq, col_zreal, col_zimag = (header.index(n) for n in ('Freq', 'Zreal', 'Zimag'))
    except ValueError:
        warnings.append(f"ZCURVE header in {filename} does not name Freq/Zreal/Zimag "
                        f"({header or 'header row missing'}); assuming standard order")
        col_freq, col_zreal, col_zimag = 3, 4, 5

    # A row must be long enough to hold the rightmost column we actually read.
    # THD is optional per row: a short or empty cell gives NaN, not a lost point.
    min_columns = max(col_freq, col_zreal, col_zimag) + 1
    col_thd = {name: header.index(name) for name in thd if name in header}

    # Extract data lines (skip ZCURVE header + column names + units = 3 lines)
    if end_line is not None:
        raw_data = lines[start_line + 3:end_line]
    else:
        raw_data = lines[start_line + 3:]

    # Parse data lines
    for line in raw_data:
        if not line.strip():
            continue

        # Convert European decimal format; a comma is never a separator here
        parts = line.replace(',', '.').split('\t')

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
                for name, col in col_thd.items():
                    thd[name].append(_float_or_nan(parts, col))
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

    current_thd, voltage_thd = (np.array(thd[name], dtype=np.float64) if name in col_thd else None
                                for name in ('Ithd', 'Vthd'))
    return LoadResult(freq_array, Z, filename, warnings=warnings,
                      current_thd=current_thd, voltage_thd=voltage_thd)


def _float_or_nan(parts: List[str], col: int) -> float:
    """Cell `col` as a finite float, NaN if it is missing, empty or not finite."""
    try:
        value = float(parts[col])
    except (ValueError, IndexError):
        return math.nan
    return value if math.isfinite(value) else math.nan


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

    # Columns by name, as for ZCURVE; tab split keeps empty cells in place.
    header = lines[start_line + 1].split('\t') if start_line + 1 < len(lines) else []
    try:
        col_t, col_vf, col_vm = (header.index(n) for n in ('T', 'Vf', 'Vm'))
    except ValueError:
        # OCV is auxiliary, so no curve beats a guessed column assignment.
        logger.warning(f"OCVCURVE header in {filename} does not name T/Vf/Vm "
                       f"({header or 'header row missing'}); OCV curve skipped")
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
        parts = lines[i].replace(',', '.').split('\t')

        try:
            t = float(parts[col_t])
            vf = float(parts[col_vf])
            vm = float(parts[col_vm])
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

    n_measured = len(frequencies)
    result.keep_points(_drop_negative_real_hf(frequencies, Z, result.warnings))
    frequencies = result.frequencies
    _check_spectrum(frequencies, result.warnings)

    # Edge case: sweep stopped before reaching the requested final frequency.
    # Only a shortfall is reported - see expected_points() on the overshoot.
    # The metadata is kept on the result so that callers need not parse the
    # file a second time for it.
    result.metadata = parse_dta_metadata(filename)
    n_expected = expected_points(result.metadata)
    # Counted before the Re(Z) < 0 drop: those points were measured.
    if n_expected is not None and n_measured < n_expected:
        result.warnings.append(
            f"Sweep may be truncated: {n_measured} points, header implies "
            f"{n_expected} (lowest measured {frequencies.min():.2e} Hz, "
            f"header FREQFINAL {result.metadata['freq_final']:.2e} Hz)")

    return result

