#!/usr/bin/env python3
"""Unit tests for the THD linearity check (validation/thd.py).

Covers:
- thd_check: threshold counting, frequency span, NaN, missing channels
- report_thd: no verdict and no figure when the THD columns hold nothing
- report_thd: the CLI section and the THD figure appear for a Gamry file
  recorded with THD (example/EISPOT-test1.DTA); nothing for one without
"""

import argparse
import logging
import os
from types import SimpleNamespace

import numpy as np
import pytest

# Suppress matplotlib GUI (report_thd draws a figure)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.cli import load_eis_data, report_thd
from eis_analysis.validation import thd_check, THD_THRESHOLD

EXAMPLE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example")

F = np.array([1e3, 1e2, 1e1, 1e0, 1e-1])


def test_none_without_thd():
    assert thd_check(F, None, None) is None


@pytest.mark.parametrize("delta, n_above", [(+1e-6, 1), (-1e-6, 0)])
def test_threshold_is_strict(delta, n_above):
    thd = np.full(len(F), 1e-3)
    thd[2] = THD_THRESHOLD + delta
    ch = thd_check(F, thd, None).current
    assert ch.n_above == n_above
    assert ch.f_above_min == (1e1 if n_above else None)


def test_summary_and_span():
    current = np.array([0.001, 0.02, 0.002, 0.015, np.nan])
    result = thd_check(F, current, None)
    ch = result.current
    assert result.voltage is None and result.n_points == 5
    assert ch.n_above == 2
    assert (ch.f_above_min, ch.f_above_max) == (1e0, 1e2)
    assert ch.maximum == 0.02 and ch.f_at_max == 1e2
    assert ch.median == pytest.approx(np.median([0.001, 0.02, 0.002, 0.015]))  # NaN ignored


def test_all_nan_channel_warns():
    result = thd_check(F, np.full(len(F), np.nan), np.full(len(F), 1e-3))
    assert result.current is None and result.voltage is not None
    assert any("Current THD column present but holds no value" in w for w in result.warnings)


def test_zero_thd_is_missing_not_linear():
    """A column of zeros is an unmeasured column, not a perfectly linear one."""
    result = thd_check(F, np.zeros(len(F)), None)
    assert result.current is None
    assert any("holds no value" in w for w in result.warnings)


def test_plot_thd_accepts_lists_and_skips_empty_channel():
    from eis_analysis.visualization import plot_thd
    fig = plot_thd(list(F), [0.001, 0.02, 0.002, 0.015, 0.003], [np.nan] * 5, 0.01)
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    plt.close(fig)
    assert labels == ['Current', 'Threshold 1 %']


def test_report_thd_bound_not_rounded_to_zero(caplog, monkeypatch):
    """A sub-percent threshold must not print a ~0 % error bound."""
    from eis_analysis.cli.handlers import validation
    monkeypatch.setattr(validation, "thd_check",
                        lambda f, i, v: thd_check(f, i, v, threshold=0.0015))
    data = SimpleNamespace(frequencies=F, current_thd=np.full(len(F), 1e-4), voltage_thd=None)
    with caplog.at_level(logging.INFO, logger="eis_analysis.cli.handlers.validation"):
        report_thd(data, argparse.Namespace(save=None, format="png"))
    plt.close('all')
    assert "below 0.15 %" in caplog.text and "below ~0.45 %" in caplog.text


def test_partly_nan_counts_valid_points():
    """n_above is counted among the points with a value, and the gap is reported."""
    current = np.array([0.02, np.nan, np.nan, 0.001, np.nan])
    result = thd_check(F, current, None)
    assert (result.current.n_above, result.current.n_valid) == (1, 2)
    assert any("Current THD missing at 3/5 points" in w for w in result.warnings)


def test_report_thd_without_values_gives_no_verdict(caplog, tmp_path):
    """THD columns that hold nothing: warnings only - no 'All points below', no figure."""
    data = SimpleNamespace(frequencies=F, current_thd=np.full(len(F), np.nan),
                           voltage_thd=np.full(len(F), np.nan))
    args = argparse.Namespace(save=str(tmp_path / "out"), format="png")
    with caplog.at_level(logging.INFO, logger="eis_analysis.cli.handlers.validation"):
        report_thd(data, args)
    assert "holds no value" in caplog.text
    assert "All points" not in caplog.text
    assert not (tmp_path / "out_thd.png").exists()


def test_length_mismatch_raises():
    with pytest.raises(ValueError, match="current THD"):
        thd_check(F, np.zeros(3), None)


@pytest.mark.parametrize("name, expected", [
    ("EISPOT-test1.DTA", "Current THD above 1 % at 3/72 points"),
    ("real_gamry_example.DTA", None),
])
def test_report_thd_on_real_files(caplog, tmp_path, name, expected):
    """Section and figure only for a file with THD columns; values checked by hand on the file."""
    data = load_eis_data(argparse.Namespace(input=os.path.join(EXAMPLE_DIR, name)))
    args = argparse.Namespace(save=str(tmp_path / "out"), format="png")
    with caplog.at_level(logging.INFO, logger="eis_analysis.cli.handlers.validation"):
        result = report_thd(data, args)
    plt.close('all')
    figure = tmp_path / "out_thd.png"
    if expected is None:
        assert result is None and "THD" not in caplog.text
        assert not figure.exists()
    else:
        assert expected in caplog.text
        assert result.current.maximum == pytest.approx(0.0216618)
        assert figure.exists()
