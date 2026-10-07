"""
Loading EIS data from CSV files: delimiter, header names and units.
"""

import logging
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .spectrum import (MIN_DATA_POINTS, LoadResult, _check_spectrum,
                       _drop_negative_real_hf)

logger = logging.getLogger(__name__)


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
_FREQ_WORDS = {'freq', 'frequency', 'hz', 'khz', 'mhz', 'ghz'}  # plus 'f' as the first word, see below
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
                  for unit, names in (('Hz', 'Hz|hz|HZ'),
                                       ('Ohm', 'Ohms?|ohms?|OHMS?|\u03a9|\u2126'))}  # omega, ohm sign


def _unit_factor(header: str, unit: str, filename: str) -> Tuple[float, str]:
    """
    Factor from a column's unit to `unit` ('Hz' or 'Ohm'), and the unit as written.

    (1.0, '') when the header names no such unit. A frequency header with the
    word rad is an angular frequency, divided by 2 pi. Two different prefixes
    in one header raise, as does m or M on a unit written in one case ('mhz',
    'MOHM'): either guess could be off by orders of magnitude.
    """
    found: Dict[float, str] = {}  # by factor: kOhm and KOhm are one unit
    for m in _UNIT_PATTERNS[unit].finditer(header):
        prefix, written = m.group(1), m.group(0)
        # 'mhz', 'MOHM': a header in one case cannot tell milli from mega
        if prefix and prefix in 'mM' and written[1:].isascii() and (
                written.islower() or written.isupper()):
            raise ValueError(f"CSV header {header.strip()!r} of {filename}: {written!r} "
                             f"may be milli or mega; write m{unit} or M{unit}")
        found[_PREFIXES[prefix]] = written
    if len(found) > 1:
        raise ValueError(f"CSV header {header.strip()!r} of {filename} gives several units: "
                         f"{', '.join(found.values())}")
    if found:
        ((factor, written),) = found.items()
        return factor, written
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
    frequency in rad/s is divided by 2 pi. m or M on a unit in one case
    ('mhz', 'MOHM') is ambiguous and raises. Names:
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
