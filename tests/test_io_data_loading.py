#!/usr/bin/env python3
"""Unit tests for the IO module (eis_analysis/io/data_loading.py).

Covers the public DTA/CSV parsing surface (audit finding 2.4 — the parser
had zero direct unit tests despite tight coupling to user file formats):

- read_gamry_native / load_data  (Gamry .DTA parser + validation)
- load_csv_data                  (CSV parser, delimiter/column auto-detect)
- parse_dta_metadata             (metadata block)
- parse_ocv_curve                (OCVCURVE section)

Synthetic in-memory fixtures (written to tmp_path) exercise edge cases
deterministically; smoke tests against the real export
example/EISPOT-test1.DTA (OCVCURVE + ZCURVE + full metadata, European
decimals) anchor the parsers to a current real-world Gamry file.
"""

import logging
import os

import numpy as np
import pytest

from eis_analysis.io.data_loading import (
    read_gamry_native,
    load_data,
    load_csv_data,
    parse_dta_metadata,
    parse_ocv_curve,
    expected_points,
    MIN_DATA_POINTS,
)

EXAMPLE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example")
REAL_DTA = os.path.join(EXAMPLE_DIR, "EISPOT-test1.DTA")


# ---------------------------------------------------------------------------
# Fixtures / builders
# ---------------------------------------------------------------------------

def _fmt(x, decimal="."):
    """Round-trippable float formatting; ',' emulates the European Gamry format."""
    s = repr(float(x))
    return s.replace(".", ",") if decimal == "," else s


def _rows(n, fmin=0.1, fmax=1e5):
    """n (freq, Zreal, Zimag) tuples, frequency log-spaced descending."""
    freqs = np.logspace(np.log10(fmax), np.log10(fmin), n)
    return [(float(fr), 100.0 + i, -(10.0 + i)) for i, fr in enumerate(freqs)]


def _make_dta(rows, decimal=",", aborted_after=None, sweep=None):
    """Minimal valid Gamry .DTA with a ZCURVE section built from rows.

    ``sweep`` is an optional ``(f_init, f_final, pts_per_dec)`` triple written
    into the header, so expected_points() has something to work with.
    """
    lines = [
        "TAG\tEISPOT",
        "TITLE\tLABEL\tTest\tT",
    ]
    if sweep is not None:
        f_init, f_final, per_dec = sweep
        lines += [
            f"FREQINIT\tQUANT\t{_fmt(f_init, decimal)}\tInitial Fre&q. (Hz)",
            f"FREQFINAL\tQUANT\t{_fmt(f_final, decimal)}\tFinal Fre&q. (Hz)",
            f"PTSPERDEC\tQUANT\t{_fmt(per_dec, decimal)}\tPoints/&decade",
        ]
    lines += [
        "ZCURVE\tTABLE",
        "\tPt\tTime\tFreq\tZreal\tZimag",
        "\t#\ts\tHz\tohm\tohm",
    ]
    for i, (fr, zr, zi) in enumerate(rows):
        lines.append(f"\t{i}\t{i}\t{_fmt(fr, decimal)}\t{_fmt(zr, decimal)}\t{_fmt(zi, decimal)}")
        if aborted_after is not None and i + 1 == aborted_after:
            lines.append("EXPERIMENTABORTED")
    return "\n".join(lines) + "\n"


def _make_csv(rows, delimiter=",", decimal=".", headers=("frequency", "Z_real", "Z_imag"), comments=()):
    lines = list(comments)
    lines.append(delimiter.join(headers))
    for fr, zr, zi in rows:
        lines.append(delimiter.join([_fmt(fr, decimal), _fmt(zr, decimal), _fmt(zi, decimal)]))
    return "\n".join(lines) + "\n"


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def _write_encoded(tmp_path, name, text, encoding):
    """Write a file in the code page a Gamry machine would have used."""
    p = tmp_path / name
    p.write_bytes(text.encode(encoding))
    return str(p)


# ---------------------------------------------------------------------------
# A) read_gamry_native / load_data
# ---------------------------------------------------------------------------

def test_load_data_happy_path(tmp_path):
    rows = _rows(12)
    f, Z = load_data(_write(tmp_path, "ok.DTA", _make_dta(rows)))
    assert len(f) == 12 and len(Z) == 12
    assert np.iscomplexobj(Z)
    assert np.allclose(f, [r[0] for r in rows])
    assert np.isclose(Z[0].real, 100.0) and np.isclose(Z[0].imag, -10.0)


def test_european_decimal_parsed(tmp_path):
    rows = [(1000.5, 12.25, -3.75)] + _rows(11)
    f, Z = read_gamry_native(_write(tmp_path, "eu.DTA", _make_dta(rows, decimal=",")))
    assert np.isclose(f[0], 1000.5)
    assert np.isclose(Z[0].real, 12.25) and np.isclose(Z[0].imag, -3.75)


def test_experimentaborted_truncates(tmp_path):
    f, Z = load_data(_write(tmp_path, "abort.DTA", _make_dta(_rows(20), aborted_after=12)))
    assert len(f) == 12


def test_missing_zcurve_raises(tmp_path):
    text = "TAG\tEISPOT\nTITLE\tLABEL\tNo curve\tT\n"
    with pytest.raises(ValueError, match="ZCURVE"):
        read_gamry_native(_write(tmp_path, "noz.DTA", text))


def test_file_not_found_raises():
    with pytest.raises(ValueError, match="not found"):
        read_gamry_native("/nonexistent/path/missing.DTA")


def test_malformed_lines_skipped(tmp_path):
    text = "\n".join([
        "ZCURVE\tTABLE",
        "\tPt\tTime\tFreq\tZreal\tZimag",
        "\t#\ts\tHz\tohm\tohm",
        "\t0\t0\t1000.0\t100.0\t-10.0",   # valid
        "\t1\t1\tgarbage\t100.0\t-10.0",  # non-numeric freq -> skipped
        "\t2\t2\t900.0",                  # too few columns -> skipped
        "",                              # blank -> skipped
        "\t3\t3\t800.0\t101.0\t-11.0",   # valid
        "\t4\t4\t700.0\t102.0\t-12.0",   # valid
    ]) + "\n"
    f, Z = read_gamry_native(_write(tmp_path, "mal.DTA", text))
    assert len(f) == 3


def test_nonpositive_and_nonfinite_filtered(tmp_path):
    text = "\n".join([
        "ZCURVE\tTABLE",
        "\tPt\tTime\tFreq\tZreal\tZimag",
        "\t#\ts\tHz\tohm\tohm",
        "\t0\t0\t1000.0\t100.0\t-10.0",   # valid
        "\t1\t1\t-5.0\t100.0\t-10.0",     # negative freq -> filtered
        "\t2\t2\t0.0\t100.0\t-10.0",      # zero freq -> filtered
        "\t3\t3\t900.0\tinf\t-10.0",      # non-finite Zreal -> filtered
        "\t4\t4\t800.0\t100.0\t-11.0",    # valid
    ]) + "\n"
    f, Z = read_gamry_native(_write(tmp_path, "nf.DTA", text))
    assert len(f) == 2
    assert np.all(f > 0) and np.all(np.isfinite(Z))


def test_load_data_too_few_points_raises(tmp_path):
    with pytest.raises(ValueError, match="at least"):
        load_data(_write(tmp_path, "few.DTA", _make_dta(_rows(5))))


def test_load_data_duplicate_freq_warns(tmp_path, caplog):
    rows = _rows(12)
    rows[1] = (rows[0][0], rows[1][1], rows[1][2])  # duplicate first frequency
    with caplog.at_level(logging.WARNING, logger="eis_analysis.io.data_loading"):
        f, Z = load_data(_write(tmp_path, "dup.DTA", _make_dta(rows)))
    assert len(f) == 12
    assert any("duplicate" in r.message.lower() for r in caplog.records)


def test_expected_points_from_header():
    """3 decades at 10 pts/decade is 31 points, endpoints included."""
    assert expected_points({'freq_init': 1e5, 'freq_final': 1e2,
                            'pts_per_dec': 10.0}) == 31


@pytest.mark.parametrize("metadata", [
    {},                                                             # non-EIS file
    {'freq_init': 1e5, 'freq_final': None, 'pts_per_dec': 10.0},    # partial header
    {'freq_init': 1e2, 'freq_final': 1e5, 'pts_per_dec': 10.0},     # inverted range
])
def test_expected_points_returns_none_without_usable_header(metadata):
    assert expected_points(metadata) is None


def test_truncated_sweep_warns(tmp_path, caplog):
    """A run stopped above FREQFINAL is short of the header count."""
    rows = _rows(20, fmin=1e2, fmax=1e5)  # header asks for 31
    path = _write(tmp_path, "short.DTA",
                  _make_dta(rows, sweep=(1e5, 1e2, 10.0)))
    with caplog.at_level(logging.WARNING, logger="eis_analysis.io.data_loading"):
        f, _ = load_data(path)
    assert len(f) == 20
    assert any("truncated" in r.message.lower() for r in caplog.records)


@pytest.mark.parametrize("n", [31, 32])  # exact, and Gamry's one-point overshoot
def test_complete_sweep_does_not_warn(tmp_path, caplog, n):
    """A full sweep, and the endpoint overshoot, stay silent."""
    path = _write(tmp_path, f"full{n}.DTA",
                  _make_dta(_rows(n, fmin=1e2, fmax=1e5), sweep=(1e5, 1e2, 10.0)))
    with caplog.at_level(logging.WARNING, logger="eis_analysis.io.data_loading"):
        load_data(path)
    assert not any("truncated" in r.message.lower() for r in caplog.records)


# ---------------------------------------------------------------------------
# B) load_csv_data
# ---------------------------------------------------------------------------

def test_csv_comma_standard(tmp_path):
    rows = _rows(12)
    f, Z = load_csv_data(_write(tmp_path, "c.csv", _make_csv(rows)))
    assert len(f) == 12
    assert np.allclose(f, [r[0] for r in rows])
    assert np.isclose(Z[0].real, 100.0)


def test_csv_semicolon_european(tmp_path):
    text = _make_csv(_rows(12), delimiter=";", decimal=",", headers=("freq", "Zreal", "Zimag"))
    f, Z = load_csv_data(_write(tmp_path, "eu.csv", text))
    assert len(f) == 12
    assert np.isclose(Z[0].real, 100.0) and np.isclose(Z[0].imag, -10.0)


def test_csv_tab_delimited(tmp_path):
    text = _make_csv(_rows(12), delimiter="\t", headers=("f", "Re(Z)", "Im(Z)"))
    f, Z = load_csv_data(_write(tmp_path, "t.csv", text))
    assert len(f) == 12


def test_csv_comment_lines_skipped(tmp_path):
    text = _make_csv(_rows(12), comments=("# exported data", "# units: Hz, Ohm"))
    f, Z = load_csv_data(_write(tmp_path, "cm.csv", text))
    assert len(f) == 12


def test_csv_header_autodetect(tmp_path):
    text = _make_csv(_rows(12), headers=("Frequency [Hz]", "Re(Z)", "Im(Z)"))
    f, Z = load_csv_data(_write(tmp_path, "ad.csv", text))
    assert len(f) == 12
    assert np.isclose(Z[0].real, 100.0) and np.isclose(Z[0].imag, -10.0)


def test_csv_positional_fallback(tmp_path, caplog):
    # Headers that match no known pattern -> fall back to columns 0,1,2.
    text = _make_csv(_rows(12), headers=("alpha", "beta", "gamma"))
    with caplog.at_level(logging.WARNING, logger="eis_analysis.io.data_loading"):
        f, Z = load_csv_data(_write(tmp_path, "pos.csv", text))
    assert len(f) == 12
    assert any("positional" in r.message.lower() for r in caplog.records)


def test_csv_header_only_raises(tmp_path):
    with pytest.raises(ValueError, match="header"):
        load_csv_data(_write(tmp_path, "ho.csv", "frequency,Z_real,Z_imag\n"))


def test_csv_too_few_points_raises(tmp_path):
    with pytest.raises(ValueError, match="at least"):
        load_csv_data(_write(tmp_path, "fewc.csv", _make_csv(_rows(5))))


# ---------------------------------------------------------------------------
# C) parse_dta_metadata
# ---------------------------------------------------------------------------

def _full_metadata_dta():
    return "\n".join([
        "TAG\tEISPOT",
        "TITLE\tLABEL\tMy Sample\tTest &Identifier",
        "DATE\tLABEL\t01.01.2026\tDate",
        "TIME\tLABEL\t12:34:56\tTime",
        "NOTES\tNOTES\t2\t&Notes",
        "\tFirst note",
        "\tSecond note",
        "PSTAT\tPSTAT\tREF600\tPotentiostat",
        "VDC\tPOTEN\t0,5\tT\tDC Voltage",
        "FREQINIT\tQUANT\t1,00000E+005\tInit",
        "FREQFINAL\tQUANT\t1,00000E-002\tFinal",
        "PTSPERDEC\tQUANT\t1,00000E+001\tPts",
        "VAC\tQUANT\t1,00000E+001\tAC",
        "AREA\tQUANT\t2,5\tArea",
        "ZCURVE\tTABLE",
    ]) + "\n"


def test_metadata_full(tmp_path):
    m = parse_dta_metadata(_write(tmp_path, "meta.DTA", _full_metadata_dta()))
    assert m["title"] == "My Sample"
    assert m["date"] == "01.01.2026" and m["time"] == "12:34:56"
    assert m["notes"] == ["First note", "Second note"]
    assert m["pstat"] == "REF600"
    assert np.isclose(m["vdc"], 0.5)
    assert np.isclose(m["freq_init"], 1e5) and np.isclose(m["freq_final"], 1e-2)
    assert np.isclose(m["pts_per_dec"], 10.0)
    assert np.isclose(m["vac"], 10.0) and np.isclose(m["area"], 2.5)


def test_metadata_missing_fields_none(tmp_path):
    m = parse_dta_metadata(_write(tmp_path, "min.DTA", "TAG\tEISPOT\nZCURVE\tTABLE\n"))
    for key in ("area", "vdc", "vac", "freq_init", "freq_final", "pts_per_dec", "date", "time", "pstat", "title"):
        assert m[key] is None
    assert m["notes"] == []


def test_metadata_malformed_value_no_crash(tmp_path):
    text = "TAG\tEISPOT\nAREA\tQUANT\tNOTANUMBER\tArea\nZCURVE\tTABLE\n"
    m = parse_dta_metadata(_write(tmp_path, "bad.DTA", text))
    assert m["area"] is None  # unparseable value left as default, no exception


def test_metadata_missing_file_returns_defaults():
    # Narrowed except (OSError) still handles a missing file gracefully.
    m = parse_dta_metadata("/nonexistent/path/missing.DTA")
    assert m["area"] is None and m["notes"] == []


# Operator notes are typed on a Czech Windows machine, so the file arrives in
# cp1250. Reading it as UTF-8 with errors='ignore' used to delete every
# accented character without a trace.
ACCENTED_DTA = "\n".join([
    "TAG\tEISPOT",
    "TITLE\tLABEL\tM\u011b\u0159en\u00ed vzorku\tTest &Identifier",
    "NOTES\tNOTES\t2\t&Notes",
    "\t380 \u00b0C, 252 bar",
    "\tautokl\u00e1v, elektroda \u010d. 2",
    "ZCURVE\tTABLE",
]) + "\n"


@pytest.mark.parametrize("encoding", ["cp1250", "utf-8"])
def test_metadata_diacritics_folded_not_dropped(tmp_path, encoding):
    """Accents are folded to ASCII; the letters underneath must survive."""
    path = _write_encoded(tmp_path, f"acc_{encoding}.DTA", ACCENTED_DTA, encoding)
    m = parse_dta_metadata(path)
    assert m["title"] == "Mereni vzorku"
    assert m["notes"] == ["380 C, 252 bar", "autoklav, elektroda c. 2"]


def test_metadata_undecodable_bytes_do_not_raise(tmp_path):
    """Bytes that are valid in neither UTF-8 nor cp1250 fall back to latin-1."""
    path = tmp_path / "raw.DTA"
    path.write_bytes(b"TAG\tEISPOT\nTITLE\tLABEL\tA\x81B\tT\nZCURVE\tTABLE\n")
    assert parse_dta_metadata(str(path))["title"] == "AB"


def test_ocv_reads_cp1250_file(tmp_path):
    """The OCV parser shares the reader, so it must tolerate the same file."""
    text = ACCENTED_DTA.replace(
        "ZCURVE\tTABLE",
        "OCVCURVE\tTABLE\t2\n\tPt\tT\tVf\tVm\n\t#\ts\tV\tV\n"
        "\t0\t2,5\t-0,0739\t-0,0739\n\t1\t5,0\t-0,0738\t-0,0738",
    )
    ocv = parse_ocv_curve(_write_encoded(tmp_path, "ocv1250.DTA", text, "cp1250"))
    assert ocv is not None and len(ocv["time"]) == 2


def test_zcurve_data_unaffected_by_folding(tmp_path):
    """Folding touches text only; the numeric columns must round-trip."""
    rows = _rows(15)
    text = _make_dta(rows).replace("TITLE\tLABEL\tTest\tT",
                                   "TITLE\tLABEL\tVzorek \u010d. 1 p\u0159i 380 \u00b0C\tT")
    f, Z = read_gamry_native(_write_encoded(tmp_path, "num.DTA", text, "cp1250"))
    assert len(f) == 15
    assert np.isclose(f[0], rows[0][0]) and np.isclose(Z[0].real, rows[0][1])


# ---------------------------------------------------------------------------
# D) parse_ocv_curve
# ---------------------------------------------------------------------------

def _make_ocv_dta(n=5):
    lines = ["TAG\tEISPOT", f"OCVCURVE\tTABLE\t{n}", "\tPt\tT\tVf\tVm\tAch", "\t#\ts\tV\tV\tbits"]
    for i in range(n):
        t, vf, vm = 0.25 * (i + 1), -0.44 - i * 0.001, -0.44 - i * 0.002
        lines.append(f"\t{i}\t{_fmt(t, ',')}\t{_fmt(vf, ',')}\t{_fmt(vm, ',')}\t0,0004")
    return "\n".join(lines) + "\n"


def test_ocv_happy_path(tmp_path):
    ocv = parse_ocv_curve(_write(tmp_path, "ocv.DTA", _make_ocv_dta(5)))
    assert ocv is not None
    assert set(ocv) == {"time", "Vf", "Vm"}
    assert len(ocv["time"]) == 5 and len(ocv["Vf"]) == 5 and len(ocv["Vm"]) == 5
    assert np.isclose(ocv["time"][0], 0.25)


def test_ocv_missing_section_returns_none(tmp_path):
    assert parse_ocv_curve(_write(tmp_path, "noocv.DTA", _make_dta(_rows(12)))) is None


def test_ocv_missing_file_returns_none():
    # Narrowed except (OSError) still handles a missing file gracefully.
    assert parse_ocv_curve("/nonexistent/path/missing.DTA") is None


# ---------------------------------------------------------------------------
# E) Smoke tests against the real Gamry export
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.exists(REAL_DTA), reason="example/EISPOT-test1.DTA missing")
def test_smoke_load_real_dta():
    f, Z = load_data(REAL_DTA)
    assert len(f) >= MIN_DATA_POINTS
    assert np.all(f > 0) and np.all(np.isfinite(Z))


@pytest.mark.skipif(not os.path.exists(REAL_DTA), reason="example/EISPOT-test1.DTA missing")
def test_smoke_metadata_real_dta():
    m = parse_dta_metadata(REAL_DTA)
    assert np.isclose(m["area"], 1.0)
    assert np.isclose(m["vac"], 10.0)
    assert np.isclose(m["freq_init"], 1e5)
    assert np.isclose(m["freq_final"], 1e-3)
    assert np.isclose(m["pts_per_dec"], 10.0)
    assert m["title"] == "Potentiostatic EIS"
    assert m["pstat"] == "REF620-51061"
    assert m["date"] == "15.5.2026"


@pytest.mark.skipif(not os.path.exists(REAL_DTA), reason="example/EISPOT-test1.DTA missing")
def test_smoke_ocv_real_dta():
    ocv = parse_ocv_curve(REAL_DTA)
    assert ocv is not None
    assert len(ocv["time"]) == 1200


@pytest.mark.parametrize("name", ["example_eis_data.csv", "example_eis_data_eu.csv"])
def test_smoke_load_example_csv(name):
    path = os.path.join(EXAMPLE_DIR, name)
    if not os.path.exists(path):
        pytest.skip(f"{name} missing")
    f, Z = load_csv_data(path)
    assert len(f) >= MIN_DATA_POINTS
    assert np.all(f > 0)
