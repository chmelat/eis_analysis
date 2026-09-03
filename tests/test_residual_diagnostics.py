"""Tests for residual shape diagnostics (fitting/residual_diagnostics.py).

The module answers one question - are these residuals noise, or is the model
missing something - so the tests are built around five residual shapes whose
answer is known by construction:

    white noise        not systematic
    ripple             systematic; amplitude recovered, no trend
    single bump        systematic; no trend, structure without a trend
    monotone drift     systematic; the trend carries it, the structure is small
    drift + ripple     systematic; both are reported, neither hides the other

Plus the degenerate inputs that would otherwise divide by zero, and one
end-to-end check on a real under-parametrized circuit fit.
"""

import numpy as np
import pytest

from eis_analysis.fitting import C, R, fit_equivalent_circuit
from eis_analysis.fitting.residual_diagnostics import (
    MIN_PERIODOGRAM_POWER,
    analyze_residuals,
    residual_structure,
    lag1_autocorrelation,
    runs_test,
)

FREQ = np.logspace(-2, 5, 80)
LOG_FREQ = np.log10(FREQ)
FLAT_Z = np.full(len(FREQ), 100 + 0j)


def _with_residual(shape):
    """analyze_residuals arguments whose real-part residual is exactly `shape`."""
    return FREQ, FLAT_Z, FLAT_Z - shape.astype(complex)


def _noise(seed, scale=0.5):
    return np.random.default_rng(seed).standard_normal(len(FREQ)) * scale


# --- the five shapes ---

def test_white_noise_is_not_systematic():
    d = analyze_residuals(*_with_residual(_noise(3, scale=1.0)), weighting='uniform')
    assert not d.is_systematic
    assert abs(d.real.lag1_autocorr) < 0.3
    assert d.real.runs_p > 0.05
    assert d.real.power < MIN_PERIODOGRAM_POWER   # a peak, but not structure
    assert d.real.slope_p > 0.05                  # and no trend either


def test_ripple_is_systematic_and_its_amplitude_is_recovered():
    ripple = 3 * np.sin(2 * np.pi * LOG_FREQ / 1.5) + _noise(3)
    d = analyze_residuals(*_with_residual(ripple), weighting='uniform')

    assert d.is_systematic
    assert d.real.lag1_autocorr > 0.8
    assert d.real.runs_p < 1e-6
    assert d.real.amplitude == pytest.approx(3.0, rel=0.15)
    assert d.real.power > 0.5
    # the structure carries it, not the trend
    assert abs(d.real.slope) * d.window_decades < 0.5 * d.real.amplitude


def test_single_bump_is_systematic_and_reads_as_structure_not_trend():
    """A symmetric bump has no net slope, so all of it must land on the structure line."""
    bump = 3 * np.exp(-((LOG_FREQ - 1.5) / 1.2) ** 2) + _noise(3)
    d = analyze_residuals(*_with_residual(bump), weighting='uniform')

    assert d.is_systematic
    assert d.real.runs_p < 1e-6
    assert d.real.slope_p > 0.05                        # no trend to report
    assert d.real.power > MIN_PERIODOGRAM_POWER
    assert d.real.amplitude > abs(d.real.slope) * d.window_decades


def test_monotone_drift_is_carried_by_the_trend():
    drift = 3 * np.tanh((LOG_FREQ - 1.5) / 2) + _noise(3)
    d = analyze_residuals(*_with_residual(drift), weighting='uniform')

    assert d.is_systematic
    assert d.real.lag1_autocorr > 0.8
    assert d.real.slope_p < 1e-6
    # What is left after the line is a fraction of what the line itself does.
    assert abs(d.real.slope) * d.window_decades > 5 * d.real.amplitude


def test_a_drift_and_a_ripple_are_both_reported():
    """The reason the verdict is two numbers: a strong trend used to hide the wave.

    A 3.7-decade ripple in a 7-decade window sits above the half-window mark
    that the old binary rule used to call a trend, so the drift won and the
    ripple was never mentioned - even though both are real by construction.
    """
    both = (2 * np.tanh((LOG_FREQ - 1.5) / 3)
            + 1.0 * np.sin(2 * np.pi * LOG_FREQ / 3.7) + _noise(3, scale=0.2))
    d = analyze_residuals(*_with_residual(both), weighting='uniform')

    assert d.is_systematic
    assert d.real.slope_p < 1e-6                        # the trend is there
    assert d.real.amplitude == pytest.approx(1.0, rel=0.3)
    assert d.real.power > MIN_PERIODOGRAM_POWER         # ...and so is the ripple


# --- ordering: the statistics read neighbours, so the sort is load-bearing ---

def test_result_is_independent_of_input_ordering():
    ripple = (3 * np.sin(2 * np.pi * LOG_FREQ / 1.5)).astype(complex)
    Z, Z_fit = FLAT_Z, FLAT_Z - ripple
    perm = np.random.default_rng(0).permutation(len(FREQ))

    ordered = analyze_residuals(FREQ, Z, Z_fit, weighting='uniform')
    shuffled = analyze_residuals(FREQ[perm], Z[perm], Z_fit[perm], weighting='uniform')

    assert shuffled.real.lag1_autocorr == pytest.approx(ordered.real.lag1_autocorr)
    assert shuffled.real.runs == ordered.real.runs
    assert shuffled.real.amplitude == pytest.approx(ordered.real.amplitude)


def test_descending_frequency_matches_ascending():
    """Instruments write high-to-low; that must not change the answer."""
    ripple = (3 * np.sin(2 * np.pi * LOG_FREQ / 1.5)).astype(complex)
    up = analyze_residuals(FREQ, FLAT_Z, FLAT_Z - ripple, weighting='uniform')
    down = analyze_residuals(FREQ[::-1], FLAT_Z, (FLAT_Z - ripple)[::-1],
                             weighting='uniform')
    assert down.real.runs == up.real.runs
    assert down.real.lag1_autocorr == pytest.approx(up.real.lag1_autocorr)


# --- degenerate inputs ---

def test_too_few_points_reports_a_warning_instead_of_failing():
    f = np.array([1.0, 10.0, 100.0])
    Z = np.full(3, 100 + 0j)
    d = analyze_residuals(f, Z, Z * 1.01)
    assert d.real is None and d.imag is None
    assert d.warnings and "at least 4 points" in d.warnings[0]


def test_non_positive_frequency_reports_a_warning():
    f = np.array([0.0, 1.0, 10.0, 100.0])
    Z = np.full(4, 100 + 0j)
    d = analyze_residuals(f, Z, Z * 1.01)
    assert d.real is None
    assert d.warnings and "positive frequencies" in d.warnings[0]


def test_perfect_fit_has_no_variance_to_test():
    """Zero residuals must not divide by zero anywhere."""
    d = analyze_residuals(FREQ, FLAT_Z, FLAT_Z.copy())
    assert d.real.lag1_autocorr == 0.0
    assert not d.is_systematic
    assert np.isnan(d.real.runs_z)            # every residual is at the median
    assert np.isnan(d.real.slope_p)           # no variance to regress
    assert np.isnan(d.real.amplitude)


def test_runs_test_undefined_when_every_crossing_is_on_one_side():
    """Reachable with ties: five values at the median, three above, none below."""
    runs, expected, z, p = runs_test(np.array([0.0] * 5 + [5.0] * 3))
    assert np.isnan(z) and np.isnan(p)


def test_runs_test_undefined_below_three_points():
    runs, expected, z, p = runs_test(np.array([1.0, -1.0]))
    assert np.isnan(z) and np.isnan(p)


def test_lag1_of_constant_series_is_zero():
    assert lag1_autocorrelation(np.full(10, 5.0)) == 0.0


def test_residual_structure_needs_a_window():
    """A window narrower than the shortest period asked about yields nothing."""
    power, amplitude = residual_structure(np.linspace(0.0, 0.2, 10), _noise(1)[:10])
    assert np.isnan(power) and np.isnan(amplitude)


def test_amplitude_cannot_exceed_the_residual_it_came_from():
    """The bound the power form gives for free: A <= sqrt(2 * variance).

    A least-squares sinusoid at the peak does not have this bound - near the
    Nyquist floor its basis is ill-conditioned, and on this very series it
    returned 11.3 for residuals spanning 6.
    """
    n = 10
    f = np.logspace(-2, 4, n)
    x = np.log10(f)
    drift = 3 * np.tanh((x - 1.0) / 2) + np.random.default_rng(3).standard_normal(n) * 0.5
    flat = np.full(n, 100 + 0j)
    d = analyze_residuals(f, flat, flat - drift.astype(complex), weighting='uniform')

    assert d.real.amplitude <= np.ptp(drift)


def test_a_known_slope_is_recovered():
    d = analyze_residuals(*_with_residual(0.5 * LOG_FREQ - 2 + _noise(5, scale=0.1)),
                          weighting='uniform')
    assert d.real.slope == pytest.approx(0.5, rel=0.05)
    assert d.real.slope_p < 1e-20


# --- known-value checks against the closed form ---

def test_runs_test_matches_the_analytic_expectation():
    """Perfectly alternating signs: the maximum possible number of runs."""
    alternating = np.array([1.0, -1.0] * 20)
    runs, expected, z, p = runs_test(alternating)
    assert runs == 40                          # every step is a crossing
    assert expected == pytest.approx(2 * 20 * 20 / 40 + 1)
    assert z > 0 and p < 1e-6                  # too many crossings is also


def test_lag1_of_alternating_series_is_minus_one():
    assert lag1_autocorrelation(np.array([1.0, -1.0] * 20)) == pytest.approx(
        -39 / 40, abs=0.02)


# --- end to end on a real fit ---

def _two_branch_data(seed=7):
    rng = np.random.default_rng(seed)
    truth = R(10) - (R(500) | C(1e-5)) - (R(2000) | C(1e-3))
    Z0 = truth.impedance(FREQ, truth.get_all_params())
    return Z0 + 0.01 * np.abs(Z0) * (rng.standard_normal(len(FREQ))
                                     + 1j * rng.standard_normal(len(FREQ)))


def test_correct_circuit_passes_and_under_parametrized_one_fails():
    Z = _two_branch_data()

    good, Z_good, _ = fit_equivalent_circuit(
        FREQ, Z, R(10) - (R(500) | C(1e-5)) - (R(2000) | C(1e-3)), plot=False)
    bad, Z_bad, _ = fit_equivalent_circuit(
        FREQ, Z, R(10) - (R(2000) | C(1e-3)), plot=False)

    d_good = analyze_residuals(FREQ, Z, Z_good, 'modulus')
    d_bad = analyze_residuals(FREQ, Z, Z_bad, 'modulus')

    assert not d_good.is_systematic
    assert abs(d_good.real.lag1_autocorr) < 0.3
    assert d_good.real.runs_p > 0.05

    assert d_bad.is_systematic
    assert d_bad.real.lag1_autocorr > 0.8
    assert d_bad.real.runs_p < 1e-3


def test_n_eff_collapses_when_residuals_are_correlated():
    """The number that must NOT reach AIC/BIC - this is why."""
    Z = _two_branch_data()
    _, Z_bad, _ = fit_equivalent_circuit(
        FREQ, Z, R(10) - (R(2000) | C(1e-3)), plot=False)
    d = analyze_residuals(FREQ, Z, Z_bad, 'modulus')

    # rho1 ~ 0.99 drives n_eff towards 1, where ln(n_eff) = 0 would remove the
    # BIC complexity penalty entirely.
    assert d.real.n_eff < 0.05 * len(FREQ)


# --- regressions from the v0.30.0 review ---

def test_sparse_drift_is_reported_on_the_trend_line_not_the_structure_line():
    """A 10-point spectrum is where a drift most easily leaks into the structure.

    On this sampling (6 decades, spacing 0.67) the periodogram used to read a
    pure monotone drift as a 0.62-decade ripple at power 0.98 for every seed,
    and the CLI then prescribed adding a branch - the opposite of the fix.
    Fitting the line first is what removed that; the check is that the drift
    still lands on the trend for every seed, sparse sampling and all.
    """
    n = 10
    f = np.logspace(-2, 4, n)
    x = np.log10(f)
    flat = np.full(n, 100 + 0j)

    for seed in range(6):
        drift = 3 * np.tanh((x - 1.0) / 2) + np.random.default_rng(seed).standard_normal(n) * 0.5
        d = analyze_residuals(f, flat, flat - drift.astype(complex), weighting='uniform')
        assert abs(d.real.slope) * d.window_decades > 5 * d.real.amplitude, f"seed {seed}"


# --- the shape real fits actually have (v0.31.1) ---

def test_structure_survives_a_spacing_that_spreads():
    """Measured residuals are not periodic - the reason no period is reported.

    On a reference oxide fit the extrema spread out towards high frequency
    (successive spacings of 1.6, 1.5 and 2.5 decades) while their amplitude
    decayed. A chirp reproduces that: no single period describes it, so the
    period the periodogram used to print was an amplitude-weighted average of
    local spacings, and it differed between Re and Im on the same fit. Power
    and amplitude make no such claim, and must still see the structure.
    """
    # Phase quadratic in log f -> spacing grows; amplitude decays with it.
    phase = 2 * np.pi * (LOG_FREQ - LOG_FREQ[0]) ** 2 / 25
    chirp = 3 * np.exp(-(LOG_FREQ - LOG_FREQ[0]) / 6) * np.sin(phase) + _noise(3, scale=0.2)
    d = analyze_residuals(*_with_residual(chirp), weighting='uniform')

    assert d.is_systematic                              # the runs test does not care
    assert d.real.power > MIN_PERIODOGRAM_POWER         # nor does the power collapse
    # The envelope runs from 3.0 down to 0.9 across the window, and the one
    # amplitude reported lands inside that range instead of claiming either end.
    assert 0.9 < d.real.amplitude < 3.0


def test_a_noise_part_reports_numbers_that_disown_themselves(caplog):
    """Both parts are always printed, so the noise one must read as noise.

    Nothing is suppressed any more - the reader tells the parts apart from the
    p-value and the power, which is why those two numbers are on the line.
    """
    from eis_analysis.cli.handlers.fitting import _log_residual_diagnostics

    noise_part = 0.3 * np.sin(2 * np.pi * LOG_FREQ / 1.5) + _noise(11, scale=1.5)
    drift_part = 30 * np.tanh((LOG_FREQ - 1.5) / 2)
    residual = noise_part + 1j * drift_part
    d = analyze_residuals(FREQ, FLAT_Z, FLAT_Z - residual, weighting='uniform')

    assert d.is_systematic                      # driven by the imaginary part
    assert not d.real.is_systematic
    assert d.real.power < MIN_PERIODOGRAM_POWER
    assert d.real.slope_p > 0.05

    with caplog.at_level('WARNING'):
        _log_residual_diagnostics(d)
    printed = caplog.text
    assert 'trend:' in printed and 'structure:' in printed


def test_cli_wrapper_survives_a_broken_diagnostic():
    """A diagnostic failure must never cost the user a converged fit."""
    from eis_analysis.cli.handlers import fitting as handler

    broken = np.array([1.0, 2.0])               # length mismatch -> raises
    assert handler._residual_diagnostics(FREQ, FLAT_Z, broken, 'modulus') is None
