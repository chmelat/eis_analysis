#!/usr/bin/env python3
"""Test G element (Gerischer for reaction-diffusion processes)."""

import numpy as np
import pytest
from eis_analysis.fitting import R, G, fit_equivalent_circuit
from eis_analysis.fitting.jacobian import element_jacobian


@pytest.fixture
def freq():
    """Test frequencies: 100 kHz to 0.01 Hz."""
    return np.logspace(5, -2, 60)


@pytest.fixture
def g_element_params():
    """G element parameters."""
    sigma = 100.0
    tau = 1e-3
    return sigma, tau


def test_g_impedance_matches_analytical(freq, g_element_params):
    """Test G element impedance matches analytical formula."""
    sigma, tau = g_element_params

    g_elem = G(sigma, tau)
    Z_G = g_elem.impedance(freq, [sigma, tau])

    omega = 2 * np.pi * freq
    Z_expected = sigma / np.sqrt(1 + 1j * omega * tau)

    max_diff = np.max(np.abs(Z_G - Z_expected))
    assert max_diff < 1e-10, f"Impedance differs from analytical: {max_diff}"


def test_g_limiting_behavior(g_element_params):
    """Test G element limiting behavior at low/high frequency."""
    sigma, tau = g_element_params
    g_elem = G(sigma, tau)

    # Low frequency: Z -> sigma (real)
    Z_low = g_elem.impedance(np.array([1e-6]), [sigma, tau])
    assert abs(Z_low[0].real - sigma) / sigma < 0.01, "Low freq limit error"
    assert abs(Z_low[0].imag) < 0.01 * sigma, "Low freq should be real"

    # High frequency: |Z| -> 0
    Z_high = g_elem.impedance(np.array([1e8]), [sigma, tau])
    assert np.abs(Z_high[0]) < 0.01 * sigma, "High freq limit should be ~0"


def test_g_element_properties(g_element_params):
    """Test G element computed properties."""
    sigma, tau = g_element_params
    g_elem = G(sigma, tau)

    assert g_elem.sigma == sigma
    assert g_elem.tau == tau

    expected_fc = 1 / (2 * np.pi * tau)
    assert abs(g_elem.characteristic_freq - expected_fc) < 1e-10


def test_g_series_circuit(freq, g_element_params):
    """Test series circuit R - G."""
    sigma, tau = g_element_params
    R_s = 10.0

    circuit = R(R_s) - G(sigma, tau)
    Z_series = circuit.impedance(freq, circuit.get_all_params())

    g_elem = G(sigma, tau)
    Z_G = g_elem.impedance(freq, [sigma, tau])
    Z_expected = R_s + Z_G

    max_diff = np.max(np.abs(Z_series - Z_expected))
    assert max_diff < 1e-10, f"Series circuit error: {max_diff}"


def test_g_parallel_circuit(freq, g_element_params):
    """Test parallel circuit R | G."""
    sigma, tau = g_element_params
    R_p = 200.0

    circuit = R(R_p) | G(sigma, tau)
    Z_par = circuit.impedance(freq, circuit.get_all_params())

    g_elem = G(sigma, tau)
    Z_G = g_elem.impedance(freq, [sigma, tau])
    Z_expected = 1 / (1/R_p + 1/Z_G)

    max_diff = np.max(np.abs(Z_par - Z_expected))
    assert max_diff < 1e-10, f"Parallel circuit error: {max_diff}"


def test_g_jacobian(freq, g_element_params):
    """Test analytical Jacobian matches numerical."""
    sigma, tau = g_element_params
    g_elem = G(sigma, tau)

    Z_jac, dZ_jac = element_jacobian(g_elem, freq, [sigma, tau])

    # Numerical Jacobian
    eps = 1e-7
    dZ_dsigma_num = (g_elem.impedance(freq, [sigma + eps, tau]) -
                    g_elem.impedance(freq, [sigma - eps, tau])) / (2 * eps)
    dZ_dtau_num = (g_elem.impedance(freq, [sigma, tau + eps]) -
                  g_elem.impedance(freq, [sigma, tau - eps])) / (2 * eps)

    rel_err_sigma = np.max(np.abs(dZ_jac[:, 0] - dZ_dsigma_num)) / np.max(np.abs(dZ_jac[:, 0]))
    rel_err_tau = np.max(np.abs(dZ_jac[:, 1] - dZ_dtau_num)) / np.max(np.abs(dZ_jac[:, 1]))

    assert rel_err_sigma < 1e-5, f"dZ/dsigma error: {rel_err_sigma}"
    assert rel_err_tau < 1e-5, f"dZ/dtau error: {rel_err_tau}"


def test_g_element_fitting(freq):
    """Test circuit fitting with G element."""
    R_s_true, sigma_true, tau_true = 15.0, 120.0, 5e-4

    omega = 2 * np.pi * freq
    Z_true = R_s_true + sigma_true / np.sqrt(1 + 1j * omega * tau_true)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    circuit = R(12) - G(100, 3e-4)
    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, weighting='modulus', plot=False)

    assert result.fit_error_rel < 2.0, f"Fit error too high: {result.fit_error_rel:.2f}%"

    R_s_fit, sigma_fit, tau_fit = result.params_opt
    assert abs(R_s_fit - R_s_true) / R_s_true < 0.15, "R_s recovery error > 15%"
    assert abs(sigma_fit - sigma_true) / sigma_true < 0.15, "sigma recovery error > 15%"
    assert abs(tau_fit - tau_true) / tau_true < 0.15, "tau recovery error > 15%"


def test_g_fixed_parameters():
    """Test fixed parameter handling."""
    # One fixed
    g_fixed = G("100", 1e-3)
    assert g_fixed.fixed_params[0] == True
    assert g_fixed.fixed_params[1] == False

    # Both fixed
    g_both = G("100", "1e-3")
    assert g_both.fixed_params == [True, True]
