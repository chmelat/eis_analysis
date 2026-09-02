"""Tests for residual shape diagnostics (fitting/residual_diagnostics.py).

The module answers one question - are these residuals noise, or is the model
missing something - so the tests are built around four residual shapes whose
answer is known by construction:

    white noise        not systematic
    ripple             systematic, and the period is recoverable
    single bump        systematic, but a trend rather than a period
    monotone drift     systematic, trend

Plus the degenerate inputs that would otherwise divide by zero, and one
end-to-end check on a real under-parametrized circuit fit.
"""

import numpy as np
import pytest

from eis_analysis.fitting import C, R, fit_equivalent_circuit
from eis_analysis.fitting.residual_diagnostics import (
    MAX_PERIOD_FRACTION,
    analyze_residuals,
    dominant_period,
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


# --- the four shapes ---

def test_white_noise_is_not_systematic():
    d = analyze_residuals(*_with_residual(_noise(3, scale=1.0)), weighting='uniform')
    assert not d.is_systematic
    assert abs(d.real.lag1_autocorr) < 0.3
    assert d.real.runs_p > 0.05
    assert d.real.period_decades is None      # nothing to name


def test_ripple_is_systematic_and_its_period_is_recovered():
    ripple = 3 * np.sin(2 * np.pi * LOG_FREQ / 1.5) + _noise(3)
    d = analyze_residuals(*_with_residual(ripple), weighting='uniform')

    assert d.is_systematic
    assert d.real.lag1_autocorr > 0.8
    assert d.real.runs_p < 1e-6
    assert d.real.period_decades == pytest.approx(1.5, abs=0.15)
    assert d.real.power > 0.5


def test_single_bump_is_systematic_but_reports_no_period():
    """A missing relaxation is a bump, not a wave - naming a period would lie."""
    bump = 3 * np.exp(-((LOG_FREQ - 1.5) / 1.2) ** 2) + _noise(3)
    d = analyze_residuals(*_with_residual(bump), weighting='uniform')

    assert d.is_systematic
    assert d.real.runs_p < 1e-6
    assert d.real.period_decades is None
    # ...because the peak sits at the window width, not at a real wavelength
    assert d.real.period_over_window > MAX_PERIOD_FRACTION


def test_monotone_drift_is_systematic_but_reports_no_period():
    drift = 3 * np.tanh((LOG_FREQ - 1.5) / 2) + _noise(3)
    d = analyze_residuals(*_with_residual(drift), weighting='uniform')

    assert d.is_systematic
    assert d.real.lag1_autocorr > 0.8
    assert d.real.period_decades is None
    assert d.real.period_over_window > MAX_PERIOD_FRACTION


# --- ordering: the statistics read neighbours, so the sort is load-bearing ---

def test_result_is_independent_of_input_ordering():
    ripple = (3 * np.sin(2 * np.pi * LOG_FREQ / 1.5)).astype(complex)
    Z, Z_fit = FLAT_Z, FLAT_Z - ripple
    perm = np.random.default_rng(0).permutation(len(FREQ))

    ordered = analyze_residuals(FREQ, Z, Z_fit, weighting='uniform')
    shuffled = analyze_residuals(FREQ[perm], Z[perm], Z_fit[perm], weighting='uniform')

    assert shuffled.real.lag1_autocorr == pytest.approx(ordered.real.lag1_autocorr)
    assert shuffled.real.runs == ordered.real.runs
    assert shuffled.real.period_decades == pytest.approx(ordered.real.period_decades)


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


def test_runs_test_undefined_when_every_crossing_is_on_one_side():
    """Reachable with ties: five values at the median, three above, none below."""
    runs, expected, z, p = runs_test(np.array([0.0] * 5 + [5.0] * 3))
    assert np.isnan(z) and np.isnan(p)


def test_runs_test_undefined_below_three_points():
    runs, expected, z, p = runs_test(np.array([1.0, -1.0]))
    assert np.isnan(z) and np.isnan(p)


def test_lag1_of_constant_series_is_zero():
    assert lag1_autocorrelation(np.full(10, 5.0)) == 0.0


def test_dominant_period_needs_a_window():
    """A window narrower than the shortest period asked about yields nothing."""
    period, power = dominant_period(np.linspace(0.0, 0.2, 10), _noise(1)[:10])
    assert np.isnan(period) and np.isnan(power)


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

def test_sparse_spectrum_does_not_alias_a_drift_into_a_ripple():
    """A period below the Nyquist limit of the sampling can only be an alias.

    Before the floor, a 10-point spectrum over 6 decades (spacing 0.67) read a
    pure monotone drift as a 0.62-decade ripple at power 0.98 for every seed -
    and the CLI then prescribed adding a branch, the opposite of the fix.
    """
    n = 10
    f = np.logspace(-2, 4, n)
    x = np.log10(f)
    flat = np.full(n, 100 + 0j)

    for seed in range(6):
        drift = 3 * np.tanh((x - 1.0) / 2) + np.random.default_rng(seed).standard_normal(n) * 0.5
        d = analyze_residuals(f, flat, flat - drift.astype(complex), weighting='uniform')
        assert d.real.period_decades is None, f"aliased at seed {seed}"


def test_nyquist_floor_still_admits_a_resolvable_ripple():
    """The floor must not silence periods the sampling can actually resolve."""
    ripple = 3 * np.sin(2 * np.pi * LOG_FREQ / 1.5) + _noise(3)
    d = analyze_residuals(*_with_residual(ripple), weighting='uniform')
    assert d.real.period_decades == pytest.approx(1.5, abs=0.15)


def test_shape_verdict_ignores_a_part_that_is_not_systematic():
    """Only a part that failed the runs test may name the shape.

    A noise part still has a periodogram peak somewhere; letting it speak
    would prescribe 'add a branch' while the genuinely broken part is asking
    for a different element type.
    """
    from eis_analysis.cli.handlers.fitting import _log_residual_diagnostics

    noise_part = 0.3 * np.sin(2 * np.pi * LOG_FREQ / 1.5) + _noise(11, scale=1.5)
    drift_part = 30 * np.tanh((LOG_FREQ - 1.5) / 2)
    residual = noise_part + 1j * drift_part
    d = analyze_residuals(FREQ, FLAT_Z, FLAT_Z - residual, weighting='uniform')

    assert d.is_systematic                      # driven by the imaginary part
    periods = [s.period_decades for s in (d.real, d.imag)
               if s.is_systematic and s.period_decades is not None]
    assert periods == [], "a non-systematic part must not name the shape"
    _log_residual_diagnostics(d)                # must not raise


def test_cli_wrapper_survives_a_broken_diagnostic():
    """A diagnostic failure must never cost the user a converged fit."""
    from eis_analysis.cli.handlers import fitting as handler

    broken = np.array([1.0, 2.0])               # length mismatch -> raises
    assert handler._residual_diagnostics(FREQ, FLAT_Z, broken, 'modulus') is None
