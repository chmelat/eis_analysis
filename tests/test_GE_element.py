#!/usr/bin/env python3
"""Test GE element (Gerischer for reaction-diffusion processes)."""

import numpy as np
import pytest
from eis_analysis.fitting import R, GE, fit_equivalent_circuit


@pytest.fixture
def freq():
    """Test frequencies: 100 kHz to 0.01 Hz."""
    return np.logspace(5, -2, 60)


@pytest.fixture
def g_element_params():
    """GE element parameters."""
    sigma = 100.0
    tau = 1e-3
    return sigma, tau


def test_g_impedance_matches_analytical(freq, g_element_params):
    """Test GE element impedance matches analytical formula."""
    sigma, tau = g_element_params

    g_elem = GE(sigma, tau)
    Z_G = g_elem.impedance(freq, [sigma, tau])

    omega = 2 * np.pi * freq
    Z_expected = sigma / np.sqrt(1 + 1j * omega * tau)

    max_diff = np.max(np.abs(Z_G - Z_expected))
    assert max_diff < 1e-10, f"Impedance differs from analytical: {max_diff}"


def test_g_limiting_behavior(g_element_params):
    """Test GE element limiting behavior at low/high frequency."""
    sigma, tau = g_element_params
    g_elem = GE(sigma, tau)

    # Low frequency: Z -> sigma (real)
    Z_low = g_elem.impedance(np.array([1e-6]), [sigma, tau])
    assert abs(Z_low[0].real - sigma) / sigma < 0.01, "Low freq limit error"
    assert abs(Z_low[0].imag) < 0.01 * sigma, "Low freq should be real"

    # High frequency: |Z| -> 0
    Z_high = g_elem.impedance(np.array([1e8]), [sigma, tau])
    assert np.abs(Z_high[0]) < 0.01 * sigma, "High freq limit should be ~0"


def test_g_element_properties(g_element_params):
    """Test GE element computed properties."""
    sigma, tau = g_element_params
    g_elem = GE(sigma, tau)

    assert g_elem.sigma == sigma
    assert g_elem.tau == tau

    expected_fc = 1 / (2 * np.pi * tau)
    assert abs(g_elem.characteristic_freq - expected_fc) < 1e-10


def test_g_element_fitting(freq):
    """Test circuit fitting with GE element."""
    R_s_true, sigma_true, tau_true = 15.0, 120.0, 5e-4

    omega = 2 * np.pi * freq
    Z_true = R_s_true + sigma_true / np.sqrt(1 + 1j * omega * tau_true)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    circuit = R(12) - GE(100, 3e-4)
    result, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, weighting='modulus')

    assert result.fit_error_rel < 2.0, f"Fit error too high: {result.fit_error_rel:.2f}%"

    R_s_fit, sigma_fit, tau_fit = result.params_opt
    assert abs(R_s_fit - R_s_true) / R_s_true < 0.15, "R_s recovery error > 15%"
    assert abs(sigma_fit - sigma_true) / sigma_true < 0.15, "sigma recovery error > 15%"
    assert abs(tau_fit - tau_true) / tau_true < 0.15, "tau recovery error > 15%"


def test_g_fixed_parameters():
    """Test fixed parameter handling."""
    # One fixed
    g_fixed = GE("100", 1e-3)
    assert g_fixed.fixed_params[0]
    assert not g_fixed.fixed_params[1]

    # Both fixed
    g_both = GE("100", "1e-3")
    assert g_both.fixed_params == [True, True]
