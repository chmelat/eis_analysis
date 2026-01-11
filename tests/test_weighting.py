#!/usr/bin/env python3
"""Test weighting types for circuit fitting."""

import numpy as np
import pytest
from eis_analysis.fitting.circuit_elements import R, C
from eis_analysis.fitting import fit_equivalent_circuit


@pytest.fixture
def synthetic_voigt_data():
    """Generate synthetic Voigt circuit data: R_s + (R || C)."""
    freq = np.logspace(4, -1, 50)  # 10 kHz to 0.1 Hz
    omega = 2 * np.pi * freq

    # True parameters
    R_s, R1, C1 = 100.0, 5000.0, 1e-6
    Z = R_s + R1 / (1 + 1j * omega * R1 * C1)

    # Add 1% noise (reproducible)
    np.random.seed(42)
    noise = 0.01 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))

    return freq, Z + noise, (R_s, R1, C1)


@pytest.mark.parametrize("weighting", ['uniform', 'sqrt', 'proportional', 'modulus'])
def test_weighting_produces_valid_fit(synthetic_voigt_data, weighting):
    """Test that each weighting type produces a valid fit."""
    freq, Z, true_params = synthetic_voigt_data
    circuit = R(100) - (R(5000) | C(1e-6))

    result, Z_fit, fig = fit_equivalent_circuit(
        freq, Z, circuit, weighting=weighting, plot=False
    )

    assert result.fit_error_rel < 5.0, f"Fit error too high: {result.fit_error_rel:.2f}%"
    assert len(result.params_opt) == 3, "Should have 3 parameters"


@pytest.mark.parametrize("weighting", ['uniform', 'sqrt', 'proportional', 'modulus'])
def test_weighting_recovers_parameters(synthetic_voigt_data, weighting):
    """Test that fitted parameters are close to true values."""
    freq, Z, true_params = synthetic_voigt_data
    R_s_true, R1_true, C1_true = true_params
    circuit = R(100) - (R(5000) | C(1e-6))

    result, _, _ = fit_equivalent_circuit(
        freq, Z, circuit, weighting=weighting, plot=False
    )

    R_s_fit, R1_fit, C1_fit = result.params_opt

    # Allow 15% tolerance due to noise
    assert abs(R_s_fit - R_s_true) / R_s_true < 0.15, f"R_s error too high"
    assert abs(R1_fit - R1_true) / R1_true < 0.15, f"R1 error too high"
    assert abs(C1_fit - C1_true) / C1_true < 0.15, f"C1 error too high"
