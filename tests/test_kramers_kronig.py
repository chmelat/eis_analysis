#!/usr/bin/env python3
"""Unit tests for the Kramers-Kronig validation core (validation/kramers_kronig.py).

Covers the public API:
- Pure helpers: compute_pseudo_chisqr, estimate_noise_percent,
  reconstruct_impedance
- Native Lin-KK fitting: lin_kk_native
- High-level wrapper + figure: kramers_kronig_validation
- Tau-range optimization: find_optimal_extend_decades
- Result dataclasses: KKResult / LinKKResult contracts

Includes a regression for the v0.13.10 auto-extend fix (auto_extend_decades
on by default removes the spurious imaginary-part residual on data with a
strong capacitive tail).
"""

import os

import numpy as np
import pytest

# Suppress matplotlib GUI (kramers_kronig_validation builds a figure)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.io import load_data
from eis_analysis.validation.kramers_kronig import (
    compute_pseudo_chisqr,
    estimate_noise_percent,
    reconstruct_impedance,
    lin_kk_native,
    kramers_kronig_validation,
    find_optimal_extend_decades,
    KKResult,
)

EXAMPLE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example"
)
REAL_DTA = os.path.join(EXAMPLE_DIR, "EISPOT-test1.DTA")


# =============================================================================
# Helpers
# =============================================================================

def voigt_impedance(frequencies, Rs, voigt):
    """KK-compliant impedance: Rs + sum_i R_i / (1 + j*omega*tau_i).

    `voigt` is a list of (R, tau) tuples. By construction this is causal,
    linear and stable, so it must pass Kramers-Kronig validation.
    """
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, Rs, dtype=complex)
    for R, tau in voigt:
        Z += R / (1 + 1j * omega * tau)
    return Z


# =============================================================================
# compute_pseudo_chisqr (Boukamp 1995)
# =============================================================================

def test_pseudo_chisqr_perfect_fit_is_zero():
    Z = np.array([3 + 4j, 5 + 12j, 8 - 6j])
    assert compute_pseudo_chisqr(Z, Z) == 0.0


def test_pseudo_chisqr_known_value():
    # |Z|^2 = 25 and 169; fit = 0 -> weighted residual = 1 + 1 = 2.0
    Z_exp = np.array([3 + 4j, 5 + 12j])
    Z_fit = np.zeros_like(Z_exp)
    assert compute_pseudo_chisqr(Z_exp, Z_fit) == pytest.approx(2.0)


def test_pseudo_chisqr_increases_with_deviation():
    Z = np.array([10 + 10j, 20 + 5j])
    small = compute_pseudo_chisqr(Z, Z * 1.01)
    large = compute_pseudo_chisqr(Z, Z * 1.10)
    assert 0 < small < large


# =============================================================================
# estimate_noise_percent (Yrjana & Bobacka 2024)
# =============================================================================

def test_estimate_noise_exact_formula():
    # sqrt(chi2 * 5000 / n): chi2=2, n=2 -> sqrt(5000)
    assert estimate_noise_percent(2.0, 2) == pytest.approx(np.sqrt(5000.0))


def test_estimate_noise_zero_chisqr():
    assert estimate_noise_percent(0.0, 50) == 0.0


def test_estimate_noise_monotonic():
    # Increases with chi-squared, decreases with point count.
    assert estimate_noise_percent(0.01, 50) < estimate_noise_percent(0.04, 50)
    assert estimate_noise_percent(0.01, 100) < estimate_noise_percent(0.01, 50)


# =============================================================================
# reconstruct_impedance
# =============================================================================

def test_reconstruct_single_voigt_matches_analytic():
    f = np.logspace(-1, 5, 40)
    omega = 2 * np.pi * f
    Rs, R, tau = 10.0, 50.0, 1e-3
    # include_L=True ignores the trailing element value, L from L_value=0
    elements = np.array([Rs, R, 0.0])
    Z = reconstruct_impedance(f, elements, np.array([tau]), L_value=0.0)
    Z_analytic = Rs + R / (1 + 1j * omega * tau)
    np.testing.assert_allclose(Z, Z_analytic, rtol=1e-12, atol=1e-12)


def test_reconstruct_adds_inductance_term():
    f = np.logspace(0, 5, 30)
    omega = 2 * np.pi * f
    Rs, R, tau, L = 5.0, 20.0, 1e-2, 1e-6
    elements = np.array([Rs, R, 0.0])
    Z = reconstruct_impedance(f, elements, np.array([tau]), L_value=L)
    Z_expected = Rs + R / (1 + 1j * omega * tau) + 1j * omega * L
    np.testing.assert_allclose(Z, Z_expected, rtol=1e-12, atol=1e-12)


def test_reconstruct_l_value_none_has_no_inductive_term():
    f = np.logspace(0, 5, 20)
    elements = np.array([5.0, 20.0, 0.0])
    tau = np.array([1e-2])
    Z = reconstruct_impedance(f, elements, tau, L_value=None, include_L=True)
    Z_no_L = reconstruct_impedance(f, elements, tau, L_value=0.0, include_L=True)
    np.testing.assert_allclose(Z, Z_no_L, rtol=1e-12, atol=1e-12)


def test_reconstruct_include_L_false_treats_all_as_resistors():
    f = np.logspace(-1, 5, 30)
    omega = 2 * np.pi * f
    Rs = 10.0
    elements = np.array([Rs, 50.0, 30.0])  # no inductance slot
    tau = np.array([1e-3, 1e-1])
    Z = reconstruct_impedance(f, elements, tau, L_value=None, include_L=False)
    Z_expected = (Rs
                  + 50.0 / (1 + 1j * omega * 1e-3)
                  + 30.0 / (1 + 1j * omega * 1e-1))
    np.testing.assert_allclose(Z, Z_expected, rtol=1e-12, atol=1e-12)


def test_reconstruct_low_frequency_limit_is_sum_of_resistances():
    f = np.array([1e-8])  # omega*tau -> 0
    Rs = 10.0
    elements = np.array([Rs, 50.0, 30.0])
    tau = np.array([1e-3, 1e-1])
    Z = reconstruct_impedance(f, elements, tau, L_value=None, include_L=False)
    assert Z[0].real == pytest.approx(Rs + 50.0 + 30.0, rel=1e-6)
    assert Z[0].imag == pytest.approx(0.0, abs=1e-3)


# =============================================================================
# KKResult dataclass contract
# =============================================================================

def test_kkresult_empty_defaults():
    r = KKResult()
    assert r.success is False
    assert r.mean_residual_real == float('inf')
    assert r.mean_residual_imag == float('inf')
    assert r.is_valid is False


def test_kkresult_valid_when_residuals_small():
    r = KKResult(
        Z_fit=np.ones(3, dtype=complex),
        residuals_real=np.array([0.01, 0.02, 0.01]),   # ~1-2%
        residuals_imag=np.array([0.01, 0.01, 0.02]),
    )
    assert r.success is True
    assert r.mean_residual_real == pytest.approx(4.0 / 3.0)
    assert r.is_valid is True


def test_kkresult_invalid_when_residual_exceeds_threshold():
    r = KKResult(
        Z_fit=np.ones(2, dtype=complex),
        residuals_real=np.array([0.01, 0.01]),
        residuals_imag=np.array([0.20, 0.20]),  # 20% >> 5%
    )
    assert r.success is True
    assert r.is_valid is False


def test_kkresult_error_keeps_success_false():
    r = KKResult(error="fitting failed")
    assert r.success is False
    assert r.is_valid is False


# =============================================================================
# lin_kk_native
# =============================================================================

def test_lin_kk_native_compliant_data_low_residuals():
    # Two-RC KK-compliant spectrum; mu_threshold=0.7 admits enough Voigt
    # elements to fit it well, so residuals should be comfortably < 5%.
    f = np.logspace(-1, 5, 60)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    r = lin_kk_native(f, Z, mu_threshold=0.7)
    assert r.is_valid is True
    assert r.mean_residual_real < 1.0
    assert r.mean_residual_imag < 1.0


def test_lin_kk_native_output_shapes():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    r = lin_kk_native(f, Z)
    assert r.Z_fit.shape == f.shape
    assert r.residuals_real.shape == f.shape
    assert r.residuals_imag.shape == f.shape
    assert r.tau.shape[0] == r.M
    # elements = [R_s, R_1..R_M, L_slot]
    assert r.elements.shape[0] == r.M + 2
    assert np.all(np.isfinite(r.Z_fit))


def test_lin_kk_native_mu_in_unit_interval():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    r = lin_kk_native(f, Z)
    assert 0.0 < r.mu <= 1.0


def test_lin_kk_native_fits_inductance_when_requested():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3)])
    r = lin_kk_native(f, Z, include_L=True)
    assert r.inductance is not None


# =============================================================================
# kramers_kronig_validation (high-level wrapper)
# =============================================================================

def test_kk_validation_success_and_returns_figure():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    result = kramers_kronig_validation(f, Z)
    try:
        assert result.success is True
        assert isinstance(result.figure, plt.Figure)
        assert result.residuals_real.shape == f.shape
        assert result.residuals_imag.shape == f.shape
        assert result.pseudo_chisqr >= 0.0
        assert result.noise_estimate >= 0.0
    finally:
        plt.close('all')


def test_kk_validation_compliant_data_is_valid():
    f = np.logspace(-1, 5, 60)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    result = kramers_kronig_validation(f, Z, mu_threshold=0.7)
    try:
        assert result.is_valid is True
        assert result.mean_residual_real < 1.0
        assert result.mean_residual_imag < 1.0
    finally:
        plt.close('all')


# =============================================================================
# find_optimal_extend_decades
# =============================================================================

def test_find_optimal_extend_within_range():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    ext, chi2, tau, elements, L, C = find_optimal_extend_decades(
        f, Z, M=7, search_range=(0.0, 1.0), n_evaluations=6
    )
    assert 0.0 <= ext <= 1.0
    assert chi2 >= 0.0
    assert tau.shape[0] == 7
    assert np.isfinite(elements).all()


# =============================================================================
# Regression: auto-extend default (v0.13.10)
# =============================================================================

@pytest.mark.skipif(not os.path.exists(REAL_DTA),
                    reason="example/EISPOT-test1.DTA missing")
def test_auto_extend_reduces_imag_residual_on_capacitive_tail():
    # On data with a strong low-frequency capacitive tail, the bounded tau
    # grid (auto_extend_decades=False) cannot reconstruct Z'' near the edge,
    # giving a large imaginary residual. Extending the grid fixes it.
    f, Z = load_data(REAL_DTA)
    off = lin_kk_native(f, Z, auto_extend_decades=False)
    on = lin_kk_native(f, Z, auto_extend_decades=True,
                       extend_decades_range=(0.0, 1.0))
    assert off.mean_residual_imag > 10.0          # spurious without extension
    assert on.mean_residual_imag < 5.0            # compliant with extension
    assert on.mean_residual_imag < off.mean_residual_imag / 5
    assert on.extend_decades > 0.0


@pytest.mark.skipif(not os.path.exists(REAL_DTA),
                    reason="example/EISPOT-test1.DTA missing")
def test_kk_validation_default_enables_auto_extend():
    # Wrapper default (auto_extend_decades=True) must not show the spurious
    # imaginary residual on the example spectrum.
    f, Z = load_data(REAL_DTA)
    result = kramers_kronig_validation(f, Z)
    try:
        assert result.success is True
        assert result.mean_residual_imag < 5.0
    finally:
        plt.close('all')


# =============================================================================
# Series capacitance (include_C, Schonleber add_cap)
# =============================================================================

def blocking_impedance(frequencies, Rs, voigt, C):
    """KK-compliant impedance with a blocking series capacitance."""
    omega = 2 * np.pi * frequencies
    return voigt_impedance(frequencies, Rs, voigt) + 1.0 / (1j * omega * C)


def test_reconstruct_adds_capacitance_term():
    f = np.logspace(-2, 5, 30)
    omega = 2 * np.pi * f
    Rs, R, tau, C = 5.0, 20.0, 1e-2, 1e-4
    elements = np.array([Rs, R, 0.0])
    Z = reconstruct_impedance(f, elements, np.array([tau]), L_value=0.0, C_value=C)
    Z_expected = Rs + R / (1 + 1j * omega * tau) + 1.0 / (1j * omega * C)
    np.testing.assert_allclose(Z, Z_expected, rtol=1e-12, atol=1e-12)


def test_reconstruct_c_value_none_has_no_capacitive_term():
    f = np.logspace(0, 5, 20)
    elements = np.array([5.0, 20.0, 0.0])
    tau = np.array([1e-2])
    Z_none = reconstruct_impedance(f, elements, tau, L_value=None, C_value=None)
    Z_default = reconstruct_impedance(f, elements, tau, L_value=None)
    np.testing.assert_allclose(Z_none, Z_default, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("fit_type", ['real', 'imag', 'complex'])
def test_estimate_R_linear_recovers_series_C(fit_type):
    # A series C has zero real part, so only include_C can represent it;
    # on clean data the fitted C must match the true value.
    from eis_analysis.fitting.voigt_chain import estimate_R_linear

    f = np.logspace(-2, 4, 50)
    Rs_true, C_true = 10.0, 1e-4
    voigt_true = [(50.0, 1e-3), (30.0, 1e-1)]
    Z = blocking_impedance(f, Rs_true, voigt_true, C_true)

    elements, _, L_value, C_value = estimate_R_linear(
        f, Z, np.array([1e-3, 1e-1]),
        include_Rs=True, include_L=True, include_C=True,
        fit_type=fit_type, allow_negative=True, weighting='modulus'
    )

    assert C_value is not None
    assert abs(C_value - C_true) / C_true < 0.02, \
        f"C = {C_value:.3e} (true {C_true:.3e}), fit_type={fit_type}"
    # C stays out of the elements array: [R_s, R_1, R_2, L]
    assert elements.shape[0] == 4


def test_lin_kk_native_series_cap_reduces_lf_residuals():
    # Regression (M136119, 2026-07-09): a blocking series capacitance makes
    # the imaginary residuals grow toward low frequencies (real fit stays
    # good) because the Voigt chain cannot represent a diverging Z''.
    # include_C must fix this; without it the data look non-compliant.
    # C=0.2 keeps 1/(omega*C) comparable to the Voigt part at f_min, and
    # mu_threshold=0.7 admits enough elements (as in the compliant-data test).
    f = np.logspace(-2, 5, 60)
    C_true = 0.2
    Z = blocking_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)], C=C_true)

    off = lin_kk_native(f, Z, mu_threshold=0.7, auto_extend_decades=True,
                        extend_decades_range=(0.0, 1.0))
    on = lin_kk_native(f, Z, mu_threshold=0.7, include_C=True,
                       auto_extend_decades=True, extend_decades_range=(0.0, 1.0))

    assert off.mean_residual_imag > 5.0           # spurious without series C
    assert off.mean_residual_real < 1.0           # ...while the real fit is good
    assert on.mean_residual_imag < 0.5            # compliant with series C
    assert on.capacitance == pytest.approx(C_true, rel=0.02)


def test_lin_kk_native_include_C_default_off():
    f = np.logspace(-1, 5, 50)
    Z = voigt_impedance(f, 10.0, [(50.0, 1e-3), (30.0, 1e-1)])
    r = lin_kk_native(f, Z)
    assert r.capacitance is None


# =============================================================================
# mu semantics in the CLI log (regression, KK audit 2026-07-03, K2)
# =============================================================================

def test_kk_cli_log_explains_mu(caplog):
    """Regression (audit K2): the CLI log must present mu as the Lin-KK stop
    value (with its threshold), not as a bare number that reads like a
    failed quality metric. On normal termination mu < threshold by
    construction."""
    import argparse
    import logging

    from eis_analysis.cli.handlers.validation import run_kk_validation

    freq = np.logspace(4, -1, 40)
    Z = voigt_impedance(freq, 100.0, [(5000.0, 5e-3)])
    args = argparse.Namespace(
        no_kk=False, mu_threshold=0.85, auto_extend=True,
        extend_decades_max=1.0, kk_series_c=False, save=None, format='png',
    )

    with caplog.at_level(logging.INFO,
                         logger='eis_analysis.cli.handlers.validation'):
        fig = run_kk_validation(freq, Z, args)

    try:
        assert fig is not None
        assert 'Lin-KK stop, threshold 0.85' in caplog.text
        # Pin the semantics the message describes: normal termination
        # ends below the threshold.
        mu_line = [ln for ln in caplog.text.splitlines() if 'Lin-KK stop' in ln][0]
        mu_value = float(mu_line.split('mu=')[1].split()[0])
        assert mu_value < 0.85
    finally:
        plt.close('all')
