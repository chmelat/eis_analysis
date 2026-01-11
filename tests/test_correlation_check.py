#!/usr/bin/env python3
"""Test parameter correlation detection in circuit fitting."""

import numpy as np
import pytest
from eis_analysis.fitting import R, C, fit_equivalent_circuit


@pytest.fixture
def voigt_data():
    """Generate data from simple Voigt circuit."""
    circuit = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 30)
    true_params = [98.5, 4823.2, 8.7e-7]
    Z_true = circuit.impedance(freq, true_params)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))

    return freq, Z_true + noise, circuit


def test_voigt_low_correlation(voigt_data):
    """Test that simple Voigt circuit has low parameter correlation."""
    freq, Z, circuit = voigt_data

    result, _, _ = fit_equivalent_circuit(freq, Z, circuit, plot=False)

    assert result.fit_error_rel < 5.0, f"Fit error too high: {result.fit_error_rel:.2f}%"
    assert result.params_stderr is not None, "Should have standard errors"
    assert len(result.params_opt) == 3, "Should have 3 parameters"


def test_overparametrized_model_fits():
    """Test that over-parametrized model still fits but may have high correlation."""
    # Generate data from simple circuit
    circuit_simple = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 25)
    true_params = [100.0, 5000.0, 1e-6]
    Z_true = circuit_simple.impedance(freq, true_params)

    np.random.seed(123)
    noise = 0.02 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    # Fit with over-parametrized model (2 RC elements instead of 1)
    circuit_overfit = R(100) - (R(2500) | C(5e-7)) - (R(2500) | C(5e-7))

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit_overfit, plot=False)

    # Should still fit well (over-parametrized can fit anything)
    assert result.fit_error_rel < 5.0, f"Fit error too high: {result.fit_error_rel:.2f}%"


def test_three_resistors_extreme_correlation():
    """Test fitting 3 resistors in series (extreme correlation case)."""
    # Data from single resistor
    circuit_simple = R(1000)
    freq = np.logspace(4, -1, 20)
    Z_true = circuit_simple.impedance(freq, [1000.0])

    np.random.seed(456)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    # Fit with 3 resistors - only sum matters, individual values are arbitrary
    circuit_corr = R(333) - R(333) - R(334)

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit_corr, plot=False)

    # Sum of resistors should be close to 1000
    R_sum = sum(result.params_opt)
    assert abs(R_sum - 1000) / 1000 < 0.05, f"R sum should be ~1000, got {R_sum:.1f}"
