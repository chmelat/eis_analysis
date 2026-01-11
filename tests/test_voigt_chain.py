#!/usr/bin/env python3
"""Test Voigt chain linear fitting and helper functions."""

import numpy as np
import pytest
from eis_analysis.fitting import fit_equivalent_circuit, fit_voigt_chain_linear
from eis_analysis.fitting.voigt_chain import (
    generate_tau_grid,
    compute_voigt_matrix,
    estimate_R_linear,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def two_voigt_data():
    """Generate synthetic data: R_s + K(R1,tau1) + K(R2,tau2)."""
    freq = np.logspace(4, 0, 40)
    omega = 2 * np.pi * freq

    R_s, R1, tau1, R2, tau2 = 100.0, 1000.0, 0.001, 2000.0, 0.01
    Z = R_s + R1 / (1 + 1j * omega * tau1) + R2 / (1 + 1j * omega * tau2)

    return freq, Z, (R_s, [R1, R2], [tau1, tau2])


@pytest.fixture
def three_voigt_data():
    """Generate synthetic data with 3 well-separated time constants."""
    freq = np.logspace(5, -2, 50)
    omega = 2 * np.pi * freq

    R_s = 50.0
    R_values = [500.0, 1500.0, 3000.0]
    tau_values = [1e-4, 1e-2, 1.0]

    Z = R_s + np.zeros(len(freq), dtype=complex)
    for R_i, tau_i in zip(R_values, tau_values):
        Z += R_i / (1 + 1j * omega * tau_i)

    return freq, Z, (R_s, R_values, tau_values)


# =============================================================================
# Tests: Internal helper functions
# =============================================================================

def test_generate_tau_grid_coverage():
    """Test tau grid covers frequency range with extension."""
    freq = np.logspace(4, -1, 30)  # 10 kHz to 0.1 Hz
    tau = generate_tau_grid(freq, n_per_decade=3, extend_decades=1.0)

    assert len(tau) > 0, "Should generate tau values"

    # Check extension - should cover below f_min
    f_from_tau_min = 1 / (2 * np.pi * tau.max())
    expected_f_min = freq.min() / 10  # 1 decade extension
    assert f_from_tau_min < expected_f_min * 1.5, "Extension should cover lower frequencies"


def test_compute_voigt_matrix_shape():
    """Test Voigt matrix has correct shape and values."""
    freq = np.array([0.1, 1.0, 10.0, 100.0])
    tau = np.array([0.01, 0.1, 1.0])
    A = compute_voigt_matrix(freq, tau, include_Rs=True)

    assert A.shape == (len(freq), len(tau) + 1), "Matrix shape mismatch"
    assert np.allclose(A[:, 0], 1.0), "First column should be 1.0 for R_s"
    assert np.all(A[:, 1:] >= 0) and np.all(A[:, 1:] <= 1), "h(omega, tau) should be in [0, 1]"


def test_estimate_R_linear_perfect_recovery(two_voigt_data):
    """Test NNLS regression recovers parameters from clean data."""
    freq, Z, (R_s_true, R_true, tau_true) = two_voigt_data

    # Estimate with exact tau values
    R_est, residual, L_value = estimate_R_linear(
        freq, Z, np.array(tau_true), include_Rs=True, include_L=False
    )

    assert np.allclose(R_est[0], R_s_true, rtol=0.01), f"R_s mismatch: {R_est[0]} vs {R_s_true}"
    assert np.allclose(R_est[1:], R_true, rtol=0.01), f"R_i mismatch: {R_est[1:]} vs {R_true}"
    assert L_value is None, "L should be None when include_L=False"


# =============================================================================
# Tests: fit_voigt_chain_linear with different options
# =============================================================================

def test_nnls_fit_produces_valid_result(two_voigt_data):
    """Test NNLS fitting (allow_negative=False)."""
    freq, Z, _ = two_voigt_data

    np.random.seed(42)
    noise = 0.01 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z + noise

    circuit, params = fit_voigt_chain_linear(
        freq, Z_noisy, n_per_decade=3, allow_negative=False
    )

    Z_fit = circuit.impedance(freq, params)
    error_rel = 100 * np.sqrt(np.mean(np.abs(Z_noisy - Z_fit)**2)) / np.mean(np.abs(Z_noisy))

    assert error_rel < 10.0, f"Fit error too high: {error_rel:.2f}%"
    assert len(params) >= 3, "Should have at least R_s + one K element"


def test_nnls_produces_nonnegative_R(two_voigt_data):
    """Test that NNLS produces only non-negative R values."""
    freq, Z, _ = two_voigt_data

    np.random.seed(42)
    noise = 0.01 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z + noise

    circuit, params = fit_voigt_chain_linear(
        freq, Z_noisy, n_per_decade=3, allow_negative=False
    )

    # R values: R_s at index 0, then R_i at odd indices (R_s, R1, tau1, R2, tau2, ...)
    # Skip last element if it's L
    param_count = len(params)
    if 'L(' in str(circuit):
        param_count -= 1

    R_values = [params[0]] + [params[i] for i in range(1, param_count, 2)]

    assert all(R >= 0 for R in R_values), f"Found negative R: {R_values}"


def test_pseudoinverse_fit_works(two_voigt_data):
    """Test pseudoinverse fitting with mu optimization."""
    freq, Z, _ = two_voigt_data

    np.random.seed(42)
    noise = 0.01 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z + noise

    circuit, params = fit_voigt_chain_linear(
        freq, Z_noisy, auto_optimize_M=True, mu_threshold=0.85
    )

    Z_fit = circuit.impedance(freq, params)
    error_rel = 100 * np.sqrt(np.mean(np.abs(Z_noisy - Z_fit)**2)) / np.mean(np.abs(Z_noisy))

    assert error_rel < 15.0, f"Fit error too high: {error_rel:.2f}%"


# =============================================================================
# Tests: Integration with nonlinear fitting
# =============================================================================

def test_voigt_chain_as_initial_guess(two_voigt_data):
    """Test using linear fit as initial guess for nonlinear optimization."""
    freq, Z, _ = two_voigt_data

    np.random.seed(42)
    noise = 0.02 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z + noise

    circuit, params = fit_voigt_chain_linear(
        freq, Z_noisy, n_per_decade=2, extend_decades=0.5, prune_threshold=0.05
    )

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, plot=False)

    assert result.fit_error_rel < 5.0, f"Fit error too high: {result.fit_error_rel:.2f}%"


def test_three_time_constants(three_voigt_data):
    """Test fitting with three well-separated time constants."""
    freq, Z, _ = three_voigt_data

    np.random.seed(123)
    noise = 0.01 * np.abs(Z) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z + noise

    circuit, params = fit_voigt_chain_linear(
        freq, Z_noisy, n_per_decade=3, extend_decades=1.0, prune_threshold=0.02
    )

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, plot=False)

    assert result.fit_error_rel < 3.0, f"Fit error too high: {result.fit_error_rel:.2f}%"
