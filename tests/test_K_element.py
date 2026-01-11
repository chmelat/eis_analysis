#!/usr/bin/env python3
"""Test K element (Voigt with tau parametrization)."""

import numpy as np
import pytest
from eis_analysis.fitting import R, C, K, fit_equivalent_circuit


@pytest.fixture
def freq():
    """Test frequencies: 10 kHz to 0.1 Hz."""
    return np.logspace(4, -1, 50)


@pytest.fixture
def k_element_params():
    """K element parameters: R=1000, tau=1e-4."""
    R_val = 1000.0
    tau_val = 1e-4
    C_val = tau_val / R_val  # equivalent capacitance
    return R_val, tau_val, C_val


def test_k_equivalent_to_parallel_rc(freq, k_element_params):
    """Test that K(R, tau) is equivalent to (R || C) where C = tau/R."""
    R_val, tau_val, C_val = k_element_params

    k_elem = K(R_val, tau_val)
    Z_K = k_elem.impedance(freq, [R_val, tau_val])

    rc_circuit = R(R_val) | C(C_val)
    Z_RC = rc_circuit.impedance(freq, [R_val, C_val])

    max_diff = np.max(np.abs(Z_K - Z_RC))
    assert max_diff < 1e-10, f"K and R||C differ by {max_diff}"


def test_k_element_properties(k_element_params):
    """Test K element computed properties."""
    R_val, tau_val, C_val = k_element_params

    k_elem = K(R_val, tau_val)

    assert abs(k_elem.capacitance - C_val) < 1e-15, "Capacitance mismatch"

    expected_fc = 1 / (2 * np.pi * tau_val)
    assert abs(k_elem.characteristic_freq - expected_fc) < 1e-10, "Characteristic freq mismatch"


def test_k_to_rc_conversion(freq, k_element_params):
    """Test K.to_RC() conversion produces equivalent impedance."""
    R_val, tau_val, C_val = k_element_params

    k_elem = K(R_val, tau_val)
    rc_from_k = k_elem.to_RC()

    Z_K = k_elem.impedance(freq, [R_val, tau_val])
    Z_RC = rc_from_k.impedance(freq, [R_val, C_val])

    max_diff = np.max(np.abs(Z_K - Z_RC))
    assert max_diff < 1e-10, f"to_RC() conversion error: {max_diff}"


def test_k_series_chain(freq):
    """Test series of K elements (Voigt chain)."""
    circuit = R(100) - K(500, 1e-4) - K(2000, 1e-3)
    params = circuit.get_all_params()

    assert len(params) == 5, "Should have 5 params: R_s, R1, tau1, R2, tau2"

    Z = circuit.impedance(freq, params)
    assert len(Z) == len(freq), "Impedance array length mismatch"
    assert np.all(np.isfinite(Z)), "Impedance contains NaN or Inf"


def test_k_element_fitting(freq):
    """Test circuit fitting with K element."""
    # True circuit: R_s + K(R, tau)
    R_s_true, R_true, tau_true = 100.0, 1000.0, 2e-4
    omega = 2 * np.pi * freq
    Z_true = R_s_true + R_true / (1 + 1j * omega * tau_true)

    # Add 1% noise
    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    # Fit
    circuit = R(90) - K(900, 1.5e-4)
    result, Z_fit, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, weighting='modulus', plot=False)

    assert result.fit_error_rel < 2.0, f"Fit error too high: {result.fit_error_rel:.2f}%"

    # Check parameter recovery
    R_s_fit, R_fit, tau_fit = result.params_opt
    assert abs(R_s_fit - R_s_true) / R_s_true < 0.1, "R_s recovery error > 10%"
    assert abs(R_fit - R_true) / R_true < 0.1, "R recovery error > 10%"
    assert abs(tau_fit - tau_true) / tau_true < 0.1, "tau recovery error > 10%"
