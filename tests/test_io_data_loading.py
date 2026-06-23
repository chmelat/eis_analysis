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


def _make_dta(rows, decimal=",", aborted_after=None):
    """Minimal valid Gamry .DTA with a ZCURVE section built from rows."""
    lines = [
        "TAG\tEISPOT",
        "TITLE\tLABEL\tTest\tT",
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
