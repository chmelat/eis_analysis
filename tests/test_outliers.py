#!/usr/bin/env python3
"""Unit tests for per-point outlier flagging (validation/outliers.py).

Covers:
- residual_percent: the shared metric, and its equivalence to the Z-HIT
  magnitude residual
- find_outliers: absolute threshold, the method-baseline factor, the global
  guard, method attribution, ordering
- Edge cases: no results, failed runs, mismatched lengths, NaN, empty input
- report_outliers: the CLI-facing report is silent when nothing is flagged
- One end-to-end check through the real KK / Z-HIT validation
"""

import argparse
from types import SimpleNamespace

import numpy as np
import pytest

# Suppress matplotlib GUI (the validation functions build figures)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.validation import (
    find_outliers,
    kramers_kronig_validation,
    zhit_validation,
    KKResult,
    ZHITResult,
)
from eis_analysis.validation.outliers import (
    residual_percent,
    _METHOD_BASELINE_FACTOR,
)


# =============================================================================
# Helpers
# =============================================================================

def fake_result(residuals_percent):
    """Minimal stand-in carrying per-point residuals of a given size [%].

    The residual is split evenly between the real and imaginary component so
    that residual_percent() reproduces the requested magnitude.
    """
    r = np.asarray(residuals_percent, dtype=float) / 100.0
    component = r / np.sqrt(2.0)
    return SimpleNamespace(residuals_real=component, residuals_imag=component)


def voigt_impedance(frequencies, Rs, voigt):
    """KK-compliant impedance: Rs + sum_i R_i / (1 + j*omega*tau_i)."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, Rs, dtype=complex)
    for R, tau in voigt:
        Z += R / (1 + 1j * omega * tau)
    return Z


# =============================================================================
# residual_percent
# =============================================================================

def test_residual_percent_is_relative_deviation():
    # res_real = 0.03, res_imag = 0.04 -> hypot = 0.05 -> 5%
    r = residual_percent(np.array([0.03]), np.array([0.04]))
    assert r[0] == pytest.approx(5.0)


def test_residual_percent_zero_for_perfect_reconstruction():
    zeros = np.zeros(5)
    assert np.all(residual_percent(zeros, zeros) == 0.0)


def test_residual_percent_matches_zhit_magnitude_residual():
    """For Z-HIT the reconstruction keeps the measured phase, so the metric
    must reduce exactly to |residuals_mag|."""
    freq = np.logspace(-1, 4, 40)
    Z = voigt_impedance(freq, Rs=10.0, voigt=[(100.0, 1e-2), (500.0, 1e0)])
    result = zhit_validation(freq, Z)
    r = residual_percent(result.residuals_real, result.residuals_imag)
    assert np.allclose(r, np.abs(result.residuals_mag), rtol=1e-9, atol=1e-9)
    plt.close(result.figure)


# =============================================================================
# find_outliers: core criterion
# =============================================================================

def test_clean_spectrum_flags_nothing():
    freq = np.logspace(-2, 5, 30)
    kk = fake_result(np.full(30, 1.0))
    report = find_outliers(freq, kk, None, max_residual=5.0)
    assert report.points == []
    assert report.skipped == []


def test_single_bad_point_is_flagged():
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 1.0)
    r[7] = 20.0
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert len(report.points) == 1
    p = report.points[0]
    assert p.frequency == pytest.approx(freq[7])
    assert p.residual_kk == pytest.approx(20.0, rel=1e-6)
    assert p.residual_zhit is None
    assert p.methods == "KK"


def test_threshold_is_strict_and_respected():
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 1.0)
    r[3] = 5.0        # exactly at the threshold -> not flagged
    r[4] = 5.001      # just over -> flagged
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert [p.frequency for p in report.points] == pytest.approx([freq[4]])


def test_lower_threshold_flags_more_points():
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 0.5)
    r[[5, 10, 15]] = [3.0, 6.0, 9.0]
    n_at = lambda t: len(find_outliers(freq, fake_result(r), None, max_residual=t).points)
    assert n_at(8.0) == 1
    assert n_at(5.0) == 2
    assert n_at(2.0) == 3


def test_method_baseline_factor_suppresses_a_noisy_reconstruction():
    """A point over the absolute threshold but within the method's own
    residual scatter must not be flagged."""
    freq = np.logspace(-2, 5, 30)
    # Baseline 4%: a 6% point is over the 5% threshold but only 1.5x the
    # median, i.e. ordinary scatter for this reconstruction.
    r = np.full(30, 4.0)
    r[9] = 6.0
    assert find_outliers(freq, fake_result(r), None, max_residual=5.0).points == []

    # Same baseline, but now far above it -> genuine outlier.
    r[9] = 4.0 * _METHOD_BASELINE_FACTOR + 1.0
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert [p.frequency for p in report.points] == pytest.approx([freq[9]])


def test_baseline_factor_does_not_bite_on_a_well_reconstructed_spectrum():
    """With a low baseline the absolute threshold alone decides."""
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 0.5)
    r[2] = 5.5
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert [p.frequency for p in report.points] == pytest.approx([freq[2]])


# =============================================================================
# find_outliers: two methods
# =============================================================================

def test_worst_of_both_methods_and_attribution():
    freq = np.logspace(-2, 5, 30)
    rk = np.full(30, 0.5)
    rz = np.full(30, 0.5)
    rk[1] = 12.0               # KK only
    rz[2] = 9.0                # Z-HIT only
    rk[3], rz[3] = 7.0, 11.0   # both
    report = find_outliers(freq, fake_result(rk), fake_result(rz), max_residual=5.0)

    by_freq = {round(p.frequency, 9): p for p in report.points}
    assert set(by_freq) == {round(freq[i], 9) for i in (1, 2, 3)}
    assert by_freq[round(freq[1], 9)].methods == "KK"
    assert by_freq[round(freq[2], 9)].methods == "Z-HIT"
    assert by_freq[round(freq[3], 9)].methods == "KK+Z-HIT"
    # `worst` takes the larger of the contributing residuals
    assert by_freq[round(freq[3], 9)].worst == pytest.approx(11.0, rel=1e-6)
    # Both residuals are reported even when only one method flagged the point
    assert by_freq[round(freq[1], 9)].residual_zhit == pytest.approx(0.5, rel=1e-6)


def test_points_are_ordered_worst_first():
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 0.5)
    r[[4, 8, 12]] = [7.0, 30.0, 15.0]
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert [p.frequency for p in report.points] == pytest.approx(
        [freq[8], freq[12], freq[4]]
    )
    assert [p.worst for p in report.points] == sorted(
        [p.worst for p in report.points], reverse=True
    )


def test_global_guard_excludes_a_method_that_fails_as_a_whole():
    """A method whose mean residual already exceeds the threshold describes a
    globally bad spectrum, not a few bad points."""
    freq = np.logspace(-2, 5, 30)
    rk = np.full(30, 0.5)
    rk[6] = 9.0
    rz = np.full(30, 12.0)      # mean well over the threshold
    report = find_outliers(freq, fake_result(rk), fake_result(rz), max_residual=5.0)

    assert report.skipped == ["Z-HIT"]
    # Z-HIT contributes neither flags nor residual columns
    assert [p.frequency for p in report.points] == pytest.approx([freq[6]])
    assert report.points[0].residual_zhit is None
    assert report.points[0].methods == "KK"


def test_global_guard_prevents_flagging_the_whole_spectrum():
    freq = np.logspace(-2, 5, 30)
    rng = np.random.default_rng(0)
    r = 20.0 + rng.normal(0, 2.0, 30)     # everything far over the threshold
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert report.points == []
    assert report.skipped == ["KK"]


# =============================================================================
# Edge cases
# =============================================================================

def test_no_results_at_all():
    report = find_outliers(np.logspace(-2, 5, 10), None, None)
    assert report.points == []
    assert report.skipped == []


def test_empty_frequency_array():
    assert find_outliers(np.array([]), None, None).points == []


def test_failed_kk_result_is_ignored():
    """KKResult(error=...) carries no residuals."""
    report = find_outliers(np.logspace(-2, 5, 10), KKResult(error="boom"), None)
    assert report.points == []


def test_empty_zhit_result_is_ignored():
    """Z-HIT returns empty arrays when the reconstruction fails."""
    empty = np.array([])
    zhit = ZHITResult(
        Z_mag_reconstructed=empty, Z_fit=np.array([], dtype=complex),
        residuals_mag=empty, residuals_real=empty, residuals_imag=empty,
        pseudo_chisqr=0.0, noise_estimate=0.0, quality=0.0, ref_freq=1.0,
        figure=None
    )
    report = find_outliers(np.logspace(-2, 5, 10), None, zhit)
    assert report.points == []


def test_mismatched_residual_length_is_ignored():
    report = find_outliers(np.logspace(-2, 5, 10), fake_result(np.full(7, 30.0)), None)
    assert report.points == []


def test_nan_residuals_never_flag():
    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 1.0)
    r[5] = np.nan
    r[6] = 20.0
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert [p.frequency for p in report.points] == pytest.approx([freq[6]])


def test_all_nan_residuals_are_ignored():
    freq = np.logspace(-2, 5, 30)
    report = find_outliers(freq, fake_result(np.full(30, np.nan)), None)
    assert report.points == []
    assert report.skipped == []


# =============================================================================
# End-to-end through the real validations
# =============================================================================

def test_injected_spike_is_the_worst_point_end_to_end():
    """A KK-compliant spectrum with one corrupted point: both validations must
    put that point at the top of the list."""
    freq = np.logspace(-1, 4, 50)
    Z = voigt_impedance(freq, Rs=10.0, voigt=[(100.0, 1e-2), (500.0, 1e0)])
    spike_index = 25
    Z_bad = Z.copy()
    Z_bad[spike_index] *= 1.20

    kk = kramers_kronig_validation(freq, Z_bad)
    zhit = zhit_validation(freq, Z_bad)
    report = find_outliers(freq, kk, zhit, max_residual=5.0)

    assert report.points, "a 20% spike must be flagged"
    assert report.points[0].frequency == pytest.approx(freq[spike_index])
    assert "KK" in report.points[0].methods

    plt.close(kk.figure)
    plt.close(zhit.figure)


def test_clean_compliant_spectrum_is_silent_end_to_end():
    freq = np.logspace(-1, 4, 50)
    Z = voigt_impedance(freq, Rs=10.0, voigt=[(100.0, 1e-2), (500.0, 1e0)])
    kk = kramers_kronig_validation(freq, Z)
    zhit = zhit_validation(freq, Z)
    report = find_outliers(freq, kk, zhit, max_residual=5.0)
    assert report.points == [], f"clean data flagged {len(report.points)} points"
    plt.close(kk.figure)
    plt.close(zhit.figure)


# =============================================================================
# CLI report
# =============================================================================

def test_report_is_silent_when_nothing_is_flagged(caplog):
    from eis_analysis.cli.handlers.validation import report_outliers

    freq = np.logspace(-2, 5, 30)
    args = argparse.Namespace(max_residual=5.0)
    with caplog.at_level('INFO'):
        report_outliers(freq, fake_result(np.full(30, 1.0)), None, args)
    assert "Suspicious points" not in caplog.text


def test_report_lists_the_flagged_frequency(caplog):
    from eis_analysis.cli.handlers.validation import report_outliers

    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 0.5)
    r[11] = 25.0
    args = argparse.Namespace(max_residual=5.0)
    with caplog.at_level('WARNING'):
        report_outliers(freq, fake_result(r), None, args)
    assert "Suspicious points (1" in caplog.text
    assert f"{freq[11]:10.3e}".strip() in caplog.text


# =============================================================================
# Regressions: the global guard must not silence the report
# =============================================================================

def test_a_few_real_outliers_do_not_trip_the_global_guard():
    """Regression: the guard used the mean, which a handful of genuine
    outliers pulls over the threshold - switching the report off exactly when
    it has something to say."""
    freq = np.logspace(-2, 5, 35)
    r = np.full(35, 0.5)
    r[:5] = 40.0
    assert np.mean(r) > 5.0, "the mean must be over the threshold for this to bite"

    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert report.skipped == []
    assert len(report.points) == 5


def test_guard_trips_when_over_half_the_points_are_bad():
    freq = np.logspace(-2, 5, 35)
    r = np.full(35, 0.5)
    r[:20] = 40.0                      # majority bad -> median over threshold
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)
    assert report.skipped == ["KK"]
    assert report.points == []


def test_single_infinite_residual_does_not_disable_the_method():
    """Regression: Z-HIT divides by |Z| with no floor, so a single inf used to
    poison the mean and drop the whole method silently."""
    freq = np.logspace(-2, 5, 35)
    r = np.full(35, 0.5)
    r[3] = np.inf
    r[10] = 30.0
    report = find_outliers(freq, fake_result(r), None, max_residual=5.0)

    assert report.skipped == []
    assert [p.frequency for p in report.points] == pytest.approx([freq[3], freq[10]])


# =============================================================================
# Figure marking
# =============================================================================

def test_each_figure_is_marked_only_with_its_own_method_s_points():
    """A band on the KK panel at a frequency KK considers fine would
    contradict the table's `flagged by` column."""
    from eis_analysis.cli.handlers.validation import report_outliers

    freq = np.logspace(-2, 5, 30)
    rk = np.full(30, 0.5)
    rk[4] = 12.0
    rz = np.full(30, 0.5)
    rz[9] = 12.0

    kk_fig, kk_axes = plt.subplots(1, 2)
    zhit_fig, zhit_axes = plt.subplots(1, 2)
    kk = fake_result(rk)
    kk.figure = kk_fig
    zhit = fake_result(rz)
    zhit.figure = zhit_fig

    args = argparse.Namespace(max_residual=5.0, save=None, format='png')
    report_outliers(freq, kk, zhit, args)

    # One axvline each: KK's own point on the KK panel, Z-HIT's on its own.
    assert len(kk_axes[1].lines) == 1
    assert len(zhit_axes[1].lines) == 1
    assert kk_axes[1].lines[0].get_xdata()[0] == pytest.approx(freq[4])
    assert zhit_axes[1].lines[0].get_xdata()[0] == pytest.approx(freq[9])

    plt.close(kk_fig)
    plt.close(zhit_fig)


def test_report_only_reads_the_documented_args(tmp_path):
    """The whole marking path runs, including the re-save with --save."""
    from eis_analysis.cli.handlers.validation import report_outliers

    freq = np.logspace(-2, 5, 30)
    r = np.full(30, 0.5)
    r[6] = 15.0
    result = fake_result(r)
    result.figure, _ = plt.subplots(1, 2)

    prefix = str(tmp_path / "out")
    args = argparse.Namespace(max_residual=5.0, save=prefix, format='png')
    report_outliers(freq, result, None, args)

    assert (tmp_path / "out_kk.png").exists()
    plt.close(result.figure)
