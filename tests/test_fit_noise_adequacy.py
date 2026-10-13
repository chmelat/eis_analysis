#!/usr/bin/env python3
"""Tests for the fit-to-noise comparison (fitting/diagnostics.py,
compare_fit_to_noise) and its CLI wiring.

The claim under test: a circuit fit and Lin-KK, converted to noise by the
same estimator, land at a ratio ~1 for a correct model and well above
FIT_NOISE_RATIO_WARN for a model that misses an arc - at any noise level.
"""

import argparse
import logging

import numpy as np
import pytest

from eis_analysis.cli import run_circuit_fitting
from eis_analysis.cli.handlers.fit_noise import format_noise_ratio
from eis_analysis.cli.handlers.model_comparison import score_candidates
from eis_analysis.fitting import (compare_fit_to_noise, fit_equivalent_circuit,
                                  NoiseAdequacyResult)
from eis_analysis.fitting.circuit_elements import Q, R
from eis_analysis.fitting.config import FIT_NOISE_RATIO_WARN, NOISE_FLOOR_PERCENT
from eis_analysis.validation import kramers_kronig_validation

FREQ = np.logspace(5, -2, 71)


def _true_circuit():
    return R(50) - (R(1000) | Q(1e-6, 0.9)) - (R(3000) | Q(1e-4, 0.8))


def _synthesize(noise_percent, seed=1):
    """Two-arc spectrum with relative Gaussian noise on Re and Im."""
    circuit = _true_circuit()
    Z = circuit.impedance(FREQ, circuit.get_all_params())
    rng = np.random.default_rng(seed)
    sigma = noise_percent / 100 * np.abs(Z)
    return Z + sigma * (rng.standard_normal(len(Z)) + 1j * rng.standard_normal(len(Z)))


def _fit(circuit, Z):
    _, Z_fit = fit_equivalent_circuit(FREQ, Z, circuit, weighting='modulus')
    return Z_fit


# =============================================================================
# compare_fit_to_noise
# =============================================================================

@pytest.mark.parametrize('noise', [0.1, 0.5, 2.0])
@pytest.mark.parametrize('seed', [1, 2])
def test_correct_model_is_adequate(noise, seed):
    Z = _synthesize(noise, seed)
    kk = kramers_kronig_validation(FREQ, Z)
    a = compare_fit_to_noise(Z, _fit(_true_circuit(), Z), kk.noise_estimate)

    assert a.adequate
    assert a.ratio < 1.5
    assert a.warnings == []


@pytest.mark.parametrize('noise', [0.1, 0.5, 2.0])
def test_missing_arc_is_inadequate(noise):
    Z = _synthesize(noise)
    kk = kramers_kronig_validation(FREQ, Z)
    a = compare_fit_to_noise(Z, _fit(R(50) - (R(4000) | Q(1e-5, 0.7)), Z),
                             kk.noise_estimate)

    assert not a.adequate
    assert a.ratio > FIT_NOISE_RATIO_WARN
    assert len(a.warnings) == 1 and 'misses structure' in a.warnings[0]


def test_noise_floor_keeps_noiseless_ratio_bounded():
    """Noiseless data: Lin-KK noise ~0.001 %; the floor caps the divisor."""
    Z = _synthesize(0.0)
    Z_fit = Z * 1.001  # 0.1 % per component -> fit noise 0.1 %
    a = compare_fit_to_noise(Z, Z_fit, 1e-6)

    assert a.noise_percent == 1e-6
    assert a.ratio == pytest.approx(a.fit_noise_percent / NOISE_FLOOR_PERCENT)
    assert a.ratio == pytest.approx(np.sqrt(0.5) * 0.1 / NOISE_FLOOR_PERCENT, rel=1e-6)


def test_ratio_is_fit_noise_over_noise():
    """Re and Im each off by 1 % of |Z| -> pseudo chi^2 = 2N * 1e-4 -> 1 %."""
    Z = np.full(10, 3.0 + 4.0j)
    Z_fit = Z + 0.05 * (1 + 1j)  # 0.05 / |Z| = 1 %
    a = compare_fit_to_noise(Z, Z_fit, 0.5)

    assert a.fit_noise_percent == pytest.approx(1.0)
    assert a.ratio == pytest.approx(2.0)
    assert a.adequate  # the limit itself is still adequate


def test_non_finite_is_reported_not_flagged():
    Z = _synthesize(0.5)
    Z_fit = Z.copy()
    Z_fit[3] = np.nan
    for a in (compare_fit_to_noise(Z, Z_fit, 0.5),
              compare_fit_to_noise(Z, Z, np.nan)):
        assert a.adequate
        assert np.isnan(a.ratio)
        assert 'not compared' in a.warnings[0]


@pytest.mark.parametrize('Z, Z_fit', [
    (np.array([], dtype=complex), np.array([], dtype=complex)),
    (np.ones(3, dtype=complex), np.ones(2, dtype=complex)),
])
def test_bad_shapes_raise(Z, Z_fit):
    with pytest.raises(ValueError):
        compare_fit_to_noise(Z, Z_fit, 0.5)


# =============================================================================
# KKResult.noise_estimate_on
# =============================================================================

def test_kk_noise_full_range_equals_estimate():
    Z = _synthesize(0.5)
    kk = kramers_kronig_validation(FREQ, Z)
    mask = np.ones(len(FREQ), dtype=bool)
    assert kk.noise_estimate_on(mask) == pytest.approx(kk.noise_estimate)


def test_kk_noise_subset_uses_only_the_fit_points():
    Z = _synthesize(0.5)
    kk = kramers_kronig_validation(FREQ, Z)
    mask = FREQ >= 1.0
    res2 = kk.residuals_real[mask] ** 2 + kk.residuals_imag[mask] ** 2
    expected = 100 * np.sqrt(np.sum(res2) / (2 * mask.sum()))

    assert kk.noise_estimate_on(mask) == pytest.approx(expected)


@pytest.mark.parametrize('mask', [np.ones(5, dtype=bool), np.zeros(len(FREQ), dtype=bool)])
def test_kk_noise_bad_mask_raises(mask):
    kk = kramers_kronig_validation(FREQ, _synthesize(0.5))
    with pytest.raises(ValueError):
        kk.noise_estimate_on(mask)


# =============================================================================
# CLI wiring
# =============================================================================

def _args(circuit):
    return argparse.Namespace(
        no_fit=False, input='x', circuit=[circuit], voigt_chain=False,
        weighting='modulus', optimizer='single', numeric_jacobian=False,
        fit_on='original', save=None, format='png', no_show=True)


def test_cli_reports_the_ratio_and_warns(caplog):
    Z = _synthesize(0.5)
    kk = kramers_kronig_validation(FREQ, Z)
    with caplog.at_level(logging.INFO, logger='eis_analysis.cli.handlers'):
        run_circuit_fitting(FREQ, Z, _args('R(50)-(R(4000)|Q(1e-5,0.7))'),
                            Z_measured=Z, noise_percent=kk.noise_estimate)

    assert 'Fit vs. noise' in caplog.text
    assert 'misses structure' in caplog.text


def test_cli_without_noise_stays_silent(caplog):
    Z = _synthesize(0.5)
    with caplog.at_level(logging.INFO, logger='eis_analysis.cli.handlers'):
        result, _ = run_circuit_fitting(FREQ, Z, _args('R(50)-(R(4000)|Q(1e-5,0.7))'))

    assert result is not None
    assert 'Fit vs. noise' not in caplog.text


def test_table_survives_a_failing_comparison():
    """Z_measured of the wrong length makes compare_fit_to_noise raise; the
    ranking of converged fits must not go down with it."""
    Z = _synthesize(0.5)
    result, _ = fit_equivalent_circuit(FREQ, Z, _true_circuit(), weighting='modulus')
    scores = score_candidates(FREQ, Z, [('a', result)], 'modulus',
                              Z_measured=Z[:-1], noise_percent=0.5)

    assert scores[0].result is result
    assert scores[0].noise is None


@pytest.mark.parametrize('ratio, adequate, cell', [
    (0.8, True, '  0.8 '),
    (9.3, False, '  9.3*'),
    (219.4, False, '219.4*'),
    (1000.0, False, '1e+03*'),
    (20080.0, False, '2e+04*'),
    (np.nan, True, '     -'),
])
def test_ratio_cell_keeps_its_width(ratio, adequate, cell):
    a = NoiseAdequacyResult(1.0, 1.0, ratio, adequate)
    assert format_noise_ratio(a) == cell
    assert len(format_noise_ratio(a)) == 6


def test_ratio_cell_without_comparison():
    assert format_noise_ratio(None) == '     -'
