"""
Data loading functions for EIS analysis.

This module provides functions to load impedance data from various file formats.
Native Gamry DTA parser - no external dependencies.
"""

import math
import re
import unicodedata
import numpy as np
import logging
from dataclasses import dataclass, field
from typing import Dict, Optional, Any, List, Tuple
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Validation constants
MIN_DATA_POINTS = 10  # Minimum number of data points for analysis
MIN_FREQUENCY_RANGE = 10  # Minimum ratio f_max/f_min


@dataclass
class LoadResult:
    """
    A spectrum as it came out of a file.

    The only points removed are the high-frequency run with Re(Z) < 0, a lead
    artifact; `warnings` says how many (see _drop_negative_real_hf).

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
    current_thd, voltage_thd : ndarray of float or None
        Total harmonic distortion of the current and voltage per point, from
        the Gamry ``Ithd``/``Vthd`` columns (written when the THD option is
        on), aligned with `frequencies`; NaN where a cell is empty. None when
        the file has no such column. A fraction, not percent: it equals
        sqrt(sum_{n=2..10} |H_n|^2) / |H_1| of the harmonic columns exactly
        (checked on example/EISPOT-test1.DTA), although the header gives
        only '#' as the unit.
    """
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    filename: str
    metadata: Optional[Dict[str, Any]] = None
    warnings: List[str] = field(default_factory=list)
    current_thd: Optional[NDArray[np.float64]] = None
    voltage_thd: Optional[NDArray[np.float64]] = None

    def keep_points(self, mask: NDArray[np.bool_]) -> None:
        """Keep only the points where `mask` is True, in every per-point field.

        The one place that knows which fields are per point, so a column
        added later cannot be left misaligned with `frequencies`.
        """
        self.frequencies, self.Z = self.frequencies[mask], self.Z[mask]
        if self.current_thd is not None:
            self.current_thd = self.current_thd[mask]
        if self.voltage_thd is not None:
            self.voltage_thd = self.voltage_thd[mask]


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


def _drop_negative_real_hf(frequencies: NDArray[np.float64], Z: NDArray[np.complex128],
                           warnings: List[str]) -> NDArray[np.bool_]:
    """
    Mask of the points to keep: all but the high-frequency run with Re(Z) < 0,
    which is noted in `warnings`. A mask rather than the trimmed arrays, so
    `LoadResult.keep_points` trims every per-point column the same way.

    No passive system has Re(Z) < 0. At the top of a sweep it is a lead
    artifact (cable inductance resonating with stray capacitance) that breaks
    Lin-KK, the DRT and the circuit fit, and turns the HF R_inf negative.
    Only the contiguous run from the highest frequency down is removed: a
    negative differential resistance (passivation, oscillating systems under
    DC bias) gives Re(Z) < 0 at low frequencies as real physics, so points
    elsewhere are kept and only noted. A run interrupted by a positive point is
    not bridged: Re(Z) flipping sign there is within noise of zero, not a
    clear artifact.

    Raises
    ------
    ValueError
        If every point has Re(Z) < 0
    """
    # Descending; stable, so a duplicated top frequency resolves in file order
    order = np.argsort(-frequencies, kind='stable')
    negative = Z.real[order] < 0
    if negative.all():
        raise ValueError("All points have Re(Z) < 0 - check the sign convention of the data")
    n_hf = int(np.argmin(negative))

    mask = np.ones(len(frequencies), dtype=bool)
    mask[order[:n_hf]] = False

    if n_hf:
        dropped = frequencies[~mask]
        warnings.append(f"Dropped {n_hf} high-frequency point(s) with Re(Z) < 0 "
                        f"({dropped.min():.2e} - {dropped.max():.2e} Hz): not possible "
                        f"for a passive system, typically a lead artifact")

    remaining = mask & (Z.real < 0)
    if remaining.any():
        f_neg = frequencies[remaining]
        warnings.append(f"{int(remaining.sum())} point(s) below the HF end have Re(Z) < 0 "
                        f"({f_neg.min():.2e} - {f_neg.max():.2e} Hz), kept: a negative "
                        f"resistance or a measurement problem")

    return mask


def _check_spectrum(frequencies: NDArray[np.float64], warnings: List[str]) -> None:
    """
    Checks every loader applies to the spectrum it read, noting caveats in `warnings`.

    Raises
    ------
    ValueError
        If there are fewer than MIN_DATA_POINTS points
    """
    if len(frequencies) < MIN_DATA_POINTS:
        raise ValueError(f"Dataset must have at least {MIN_DATA_POINTS} points, got {len(frequencies)}")

    freq_range = frequencies.max() / frequencies.min()
    if freq_range < MIN_FREQUENCY_RANGE:
        warnings.append(
            f"Small frequency range: {freq_range:.1f}x "
            f"(recommended >{MIN_FREQUENCY_RANGE}x); DRT analysis may have poor resolution")

    # A sweep is strictly monotonic in file order. Several sweeps in one file
    # break that even when the instrument logged slightly different measured
    # frequencies, which exact-equality (np.unique) would miss. Z-HIT takes
    # its phase derivative over a minimum step (MIN_DERIVATIVE_STEP), so such
    # points no longer break it; repeated sweeps that disagree show up as its
    # residuals.
    steps = np.diff(frequencies)
    n_against = int(min(np.sum(steps >= 0), np.sum(steps <= 0)))
    if n_against:
        warnings.append(
            f"Dataset contains duplicate or out-of-order frequencies ({n_against} step(s) "
            f"against the sweep direction) - several sweeps in one file?")


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


def _sweep_segments(frequencies: NDArray[np.float64]) -> List[Tuple[int, int]]:
    """
    Half-open [start, stop) runs of the points between steps strictly against
    the sweep direction (the direction most steps take; equal frequencies
    are no step). One run for a single sweep; a second sweep starts with a
    step back to the first sweep's start frequency.
    """
    steps = np.diff(frequencies)
    direction = 1.0 if np.sum(steps > 0) >= np.sum(steps < 0) else -1.0
    bounds = [0] + [int(i) + 1 for i in np.flatnonzero(direction * steps < 0)] + [len(frequencies)]
    return list(zip(bounds[:-1], bounds[1:]))


def _detect_delimiter(header_line: str) -> str:
    """
    Auto-detect CSV delimiter from header line.

    A line of numbers (a file without a header) takes the first of semicolon,
    tab, comma and whitespace that splits it into at least three numbers:
    counting characters would pick the decimal commas of "0,1;2,5;-0,3".
    A header line returns whichever of comma, tab, semicolon occurs most
    often, or ' ' (runs of whitespace, see _split) when it has none of them
    but a space. Comma is listed first so it wins a single-column header.
    """
    for candidate in (';', '\t', ',', ' '):
        fields = _split(header_line, candidate)
        if len(fields) >= 3 and all(_is_number(x) for x in fields):
            return candidate
    if not any(d in header_line for d in ',\t;') and ' ' in header_line.strip():
        return ' '
    return max(',', '\t', ';', key=header_line.count)


def _is_number(field: str) -> bool:
    """Whether a CSV field reads as a number (decimal comma allowed)."""
    try:
        float(field.strip().replace(',', '.'))
    except ValueError:
        return False
    return True


def _split(line: str, delimiter: str) -> List[str]:
    """Fields of a line; ' ' splits at runs of whitespace (aligned columns)."""
    return line.split() if delimiter == ' ' else line.split(delimiter)


# Header words per quantity. A header is split into words (_header_words), so
# 're' in 'freq' or 'im' in 'time' never count, and units or labels such as
# 'ohms', 'hz' in brackets or ZView's '(a)'/'(b)' are simply other words
# (the prefix of a unit is read separately, see _unit_factor).
_FREQ_WORDS = {'freq', 'frequency', 'hz'}  # plus 'f' as the first word, see below
_ZREAL_WORDS = {'re', 'real', 'zreal', 'zre', 'zr', 'rez', "z'"}
_ZIMAG_WORDS = {'im', 'imag', 'imaginary', 'zimag', 'zim', 'zi', 'imz', "z''"}
# Polar form: |Z| (written |Z| or abs(Z); 'modulus' stays the electric modulus
# M below) and its phase, in degrees unless a word 'rad' is in the header.
_ZMOD_WORDS = {'zmod', 'mod', 'magnitude', 'abs'}
_ZPHASE_WORDS = {'phase', 'phz', 'zphz', 'phi', 'theta', 'arg', 'angle'}
# A word naming another quantity or a derived column rules the header out:
# Re(Y) is admittance, Re(M) modulus, Re(C) capacitance, "Z' err" an error
# bar, "Zreal fit" a model curve next to the data.
_OTHER_WORDS = {'y', 'm', 'c', 'admittance', 'modulus', 'capacitance', 'permittivity',
                'conductivity', 'eps', 'epsilon', 'sigma', 'err', 'error', 'std', 'stdev',
                'fit', 'fitted', 'sim', 'calc', 'model'}
# Typographic minus, primes and curly quotes (Origin, Word) as ASCII; Z" for Z''
_TYPED = (('−', '-'), ('″', "''"), ('′', "'"), ('”', "''"),
          ('“', "''"), ('’', "'"), ('‘', "'"), ('"', "''"))


def _header_words(header: str) -> Tuple[List[str], List[str], bool]:
    """
    Words of a column header, outside brackets and in all, and a leading minus.

    "-Im(Z)/Ohm" -> (['im', 'ohm'], ['im', 'z', 'ohm'], True). Words keep
    trailing primes, so Z' and Z'' stay apart.
    """
    name = header.strip().lower()
    if len(name) >= 2 and name[0] == name[-1] == '"':  # "freq" quoted; Z" is not
        name = name[1:-1]
    for typed, plain in _TYPED:
        name = name.replace(typed, plain)
    name = name.replace('|z|', 'zmod').strip()  # '|' is no word character
    negated = name.startswith('-')
    name = name.lstrip('- ')
    words = r"[^\W_]+'*"
    outside = re.sub(r'\([^)]*\)|\[[^\]]*\]', ' ', name)
    return re.findall(words, outside), re.findall(words, name), negated


# Unit prefixes, read case-sensitively as in SI: lowercased, the MOhm of an
# oxide would be the mOhm of a battery. K is the common misspelling of k,
# u stands in for µ in ASCII-only exports.
_PREFIXES = {'G': 1e9, 'M': 1e6, 'k': 1e3, 'K': 1e3, '': 1.0,
             'm': 1e-3, 'µ': 1e-6, 'μ': 1e-6, 'u': 1e-6}
# A unit with its prefix, neither preceded nor followed by a letter: Gamry's
# Zphz holds no Hz, 'Ohmic' no Ohm, while '_', '/', '(' and '[' are no letters.
_UNIT_PATTERNS = {unit: re.compile(rf"(?<![^\W\d_])([GMkKmµμu]?)(?:{names})(?![^\W\d_])")
                  for unit, names in (('Hz', 'Hz|hz|HZ'), ('Ohm', 'Ohms?|ohms?|OHMS?|Ω|Ω'))}


def _unit_factor(header: str, unit: str, filename: str) -> Tuple[float, str]:
    """
    Factor from a column's unit to `unit` ('Hz' or 'Ohm'), and the unit as written.

    (1.0, '') when the header names no such unit. A frequency header with the
    word rad is an angular frequency, divided by 2 pi. Two different prefixes
    in one header raise: either guess could be off by orders of magnitude.
    """
    found = {m.group(1): m.group(0) for m in _UNIT_PATTERNS[unit].finditer(header)}
    if len(found) > 1:
        raise ValueError(f"CSV header {header.strip()!r} of {filename} gives several units: "
                         f"{', '.join(found.values())}")
    if found:
        ((prefix, written),) = found.items()
        return _PREFIXES[prefix], written
    if unit == 'Hz' and 'rad' in _header_words(header)[1]:
        return 1 / (2 * np.pi), 'rad/s'
    return 1.0, ''


def _classify_header(header: str) -> Tuple[Optional[str], bool]:
    """
    Quantity a column header names ('frequency', 'Z_real', 'Z_imag', 'Z_mod',
    'Z_phase' or None) and whether it holds the negative (-Im(Z), EC-Lab;
    -phase; -Z').

    Only words outside brackets name the quantity: 'Re(Z)' is Re, 'C (F)' is
    not a frequency. Words anywhere rule it out: 'Re(Y)' is not Re(Z).
    A lone 'f' counts only as the first word ('f_Hz'), not as a unit ('Cs_F').
    """
    outside, every, negated = _header_words(header)
    if not outside or _OTHER_WORDS.intersection(every):
        return None, False
    words = set(outside)
    hits = [q for q, names in (('frequency', _FREQ_WORDS), ('Z_real', _ZREAL_WORDS),
                               ('Z_imag', _ZIMAG_WORDS), ('Z_mod', _ZMOD_WORDS),
                               ('Z_phase', _ZPHASE_WORDS)) if words & names]
    if outside[0] == 'f' and 'frequency' not in hits:
        hits.append('frequency')
    # A component or the phase can be stored negated (-Im(Z), -phase, -Z');
    # a negative frequency or |Z| cannot
    if len(hits) != 1 or (negated and hits[0] in ('frequency', 'Z_mod')):
        return None, False
    return hits[0], negated


def _detect_columns(headers: List[str], filename: str
                    ) -> Optional[Tuple[int, int, int, float, float, Optional[float]]]:
    """
    Columns of frequency and the two impedance components, with their signs.

    Returns
    -------
    tuple or None
        (freq_col, a_col, b_col, sign_a, sign_b, phase_scale): a, b are Re(Z)
        and Im(Z) with phase_scale None, or |Z| and the phase with
        phase_scale the factor to radians (pi/180 for degrees, 1 for a header
        with the word 'rad'). The signs turn a stored -Im(Z), -Z' or -phase
        back. None when no header is recognised (unknown names or no
        header), leaving column order to the caller.

    Raises
    ------
    ValueError
        If neither frequency, Re(Z), Im(Z) nor frequency, |Z|, phase are
        named once each: guessing the rest would load one column as another.
        Re/Im take precedence, so the |Z| and phase columns many exports
        carry next to them are ignored.
    """
    classes = [_classify_header(h) for h in headers]
    found: Dict[str, List[int]] = {q: [] for q in
                                   ('frequency', 'Z_real', 'Z_imag', 'Z_mod', 'Z_phase')}
    for i, (quantity, _) in enumerate(classes):
        if quantity is not None:
            found[quantity].append(i)
    if not any(found.values()):
        return None

    def once(*quantities: str) -> bool:
        return all(len(found[q]) == 1 for q in quantities)

    def sign(col: int) -> float:
        return -1.0 if classes[col][1] else 1.0

    if once('frequency', 'Z_real', 'Z_imag'):
        (f_col,), (a_col,), (b_col,) = found['frequency'], found['Z_real'], found['Z_imag']
        return f_col, a_col, b_col, sign(a_col), sign(b_col), None
    if once('frequency', 'Z_mod', 'Z_phase') and not (found['Z_real'] or found['Z_imag']):
        (f_col,), (a_col,), (b_col,) = found['frequency'], found['Z_mod'], found['Z_phase']
        radians = 'rad' in _header_words(headers[b_col])[1]
        return f_col, a_col, b_col, 1.0, sign(b_col), 1.0 if radians else np.pi / 180
    read_as = ', '.join(f"{h.strip()!r} -> {q or '-'}" for h, (q, _) in zip(headers, classes))
    raise ValueError(
        f"CSV header of {filename} must name frequency, Z_real and Z_imag (or frequency, "
        f"|Z| and phase) once each; columns read as: {read_as}. "
        f"Rename them, e.g. frequency, Z_real, Z_imag")


def load_csv_data(
    filename: str,
    delimiter: Optional[str] = None
) -> LoadResult:
    """
    Load EIS data from CSV file with auto-detection of columns and delimiter.

    Automatically detects:
    - Delimiter: comma, semicolon, or tab; whitespace when the header has
      none of them (aligned columns - their names then must not contain spaces)
    - Columns: frequency, Z_real, Z_imag by header names, or frequency, |Z|
      and phase (polar form)
    - Comments: lines starting with '#' are ignored

    Column names (case-insensitive) are split into words at spaces,
    punctuation and brackets; units and labels are just other words, so
    'Frequency (Hz)', "Z' (Ohms)", 'Z Real', "Z'(a)" (ZView) all work.
    A unit prefix G, M, k, m or µ (u) on Hz or Ohm/Ω converts the column to
    Hz and Ohm ('Freq (kHz)', "Z' (MΩ)", 'Re(Z)/mOhm'), case-sensitively as
    in SI (M mega, m milli; K counts as k), noted in the warnings; a
    frequency in rad/s is divided by 2 pi. Names:
    - Frequency: a word freq, frequency or hz, or f as the first word
    - Z real: a word re, real, zreal, zre, zr, rez or z'
    - Z imag: a word im, imag, imaginary, zimag, zim, zi, imz or z'' (or z");
      with a leading minus (-Im(Z), - Z'', EC-Lab) the column holds -Im(Z)
      and is negated; likewise -Z' for Re(Z)
    - |Z|: |Z|, a word zmod, mod, magnitude or abs; phase: a word phase, phz,
      zphz, phi, theta, arg or angle, in degrees unless its header has the
      word rad, negated with a leading minus. Only used when Re(Z) and Im(Z)
      are not both named (exports often carry |Z| and phase next to them).
    Words inside brackets only rule a column out: a word for another
    quantity or a derived column (y, m, c, admittance, modulus, err, std,
    fit, ...) anywhere means it is not Z, so Re(Y) is not read as Re(Z).
    Typographic minus, primes and curly quotes count as - ' ''.
    Without any recognised name the columns are taken in order (0, 1, 2);
    a first line of numbers is no header but the first data row.

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
        If file cannot be parsed, the header names only some of the three
        quantities or one of them twice, or the file holds several sweeps
        (runs of at least MIN_DATA_POINTS separated by a step back against
        the sweep direction): they would be read as one spectrum

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
        with open(filename, 'r', encoding='utf-8-sig') as f:  # -sig: Excel "CSV UTF-8" BOM
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
    headers = _split(header_line, delimiter)
    logger.debug(f"CSV headers: {headers}")

    warnings: List[str] = []
    # A first line of numbers is the first data row of a file without a
    # header, not a header: taken as one, the first point was lost
    headerless = len(headers) >= 3 and all(_is_number(h) for h in headers)
    columns = None if headerless else _detect_columns(headers, filename)
    if columns is None:
        warnings.append("No header row, columns taken in order (frequency, Z_real, Z_imag)"
                        if headerless else
                        "Could not detect columns from headers, using positional (0, 1, 2)")
        columns = (0, 1, 2, 1.0, 1.0, None)
    freq_col, a_col, b_col, sign_a, sign_b, phase_scale = columns

    # Unit prefixes scale the columns; a polar b is the phase, see phase_scale
    roles = [(freq_col, 'Hz', 'frequency'),
             (a_col, 'Ohm', 'Z_real' if phase_scale is None else '|Z|')]
    if phase_scale is None:
        roles.append((b_col, 'Ohm', 'Z_imag'))
    factors = [1.0, sign_a, sign_b]
    conversions = []
    for k, (col, unit, name) in enumerate(roles):
        if headerless or col >= len(headers):
            break
        factor, written = _unit_factor(headers[col], unit, filename)
        factors[k] *= factor
        if factor != 1.0:
            conversions.append(f"{name} {written} -> {unit}")
    if conversions:
        warnings.append("Units from the header: " + ", ".join(conversions))
    f_factor, a_factor, b_factor = factors

    logger.debug(f"Column indices: freq={freq_col}, a={a_col}, b={b_col} "
                 f"(signs {sign_a:+.0f} {sign_b:+.0f}, "
                 f"{'polar' if phase_scale is not None else 'Re/Im'})")

    # Parse data rows
    frequencies: List[float] = []
    impedances: List[complex] = []
    line_nums: List[int] = []

    first_data = header_idx if headerless else header_idx + 1
    for line_num, line in enumerate(lines[first_data:], start=first_data + 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Handle European decimal format (comma -> dot) for semicolon-delimited
        if delimiter == ';':
            line = line.replace(',', '.')

        parts = _split(line, delimiter)

        try:
            freq = f_factor * float(parts[freq_col].replace(',', '.'))
            a = a_factor * float(parts[a_col].replace(',', '.'))
            b = b_factor * float(parts[b_col].replace(',', '.'))

            if freq > 0 and np.isfinite(freq) and np.isfinite(a) and np.isfinite(b):
                frequencies.append(freq)
                impedances.append(complex(a, b) if phase_scale is None
                                  else a * np.exp(1j * b * phase_scale))
                line_nums.append(line_num)
        except (ValueError, IndexError) as e:
            logger.debug(f"Skipping line {line_num}: {e}")
            continue

    if len(frequencies) == 0:
        raise ValueError(f"No valid data found in {filename}")

    freq_array = np.array(frequencies, dtype=np.float64)
    sweeps = [(start, stop) for start, stop in _sweep_segments(freq_array)
              if stop - start >= MIN_DATA_POINTS]
    if len(sweeps) > 1:
        lines_str = ', '.join(f"{line_nums[a]}-{line_nums[b - 1]}" for a, b in sweeps)
        raise ValueError(f"{filename} holds {len(sweeps)} sweeps (lines {lines_str}); "
                         f"split the file and load one sweep at a time")
    Z = np.array(impedances, dtype=np.complex128)
    result = LoadResult(freq_array, Z, filename, warnings=warnings)
    result.keep_points(_drop_negative_real_hf(freq_array, Z, warnings))
    _check_spectrum(result.frequencies, warnings)
    return result
