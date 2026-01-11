#!/usr/bin/env python3
"""Test hybrid lambda selection (GCV + L-curve) for DRT analysis."""

import numpy as np
import pytest
from eis_analysis.drt.gcv import (
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
)


# =============================================================================
# Helper functions
# =============================================================================

def generate_voigt_impedance(frequencies, R_inf, voigt_elements):
    """Generate impedance for Voigt circuit."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in voigt_elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


def generate_warburg_impedance(frequencies, R_inf, R_ct, C_dl, sigma_w):
    """Generate impedance for Randles circuit with Warburg."""
    omega = 2 * np.pi * frequencies
    Z_W = sigma_w * (1 - 1j) / np.sqrt(omega)
    Z_ct = R_ct / (1 + 1j * omega * R_ct * C_dl)
    return R_inf + Z_ct + Z_W


def build_drt_matrices(frequencies, Z, R_inf, n_tau=100):
    """Build matrices A, b, L for DRT problem."""
    omega = 2 * np.pi * frequencies

    f_max = frequencies.max()
    f_min = frequencies.min()
    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    d_ln_tau = np.mean(np.diff(np.log(tau)))

    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom
    A = np.vstack([A_re, A_im])

    b = np.concatenate([Z.real - R_inf, Z.imag])

    L = np.zeros((n_tau - 2, n_tau))
    for i in range(n_tau - 2):
        L[i, i] = 1
        L[i, i + 1] = -2
        L[i, i + 2] = 1

    return A, b, L


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def frequencies():
    """Test frequencies: 100kHz to 10mHz."""
    return np.logspace(5, -2, 71)


@pytest.fixture
def voigt_data(frequencies):
    """Generate Voigt circuit data."""
    R_inf = 100
    voigt_elements = [(1000, 1e-3), (2000, 1e-1)]
    Z = generate_voigt_impedance(frequencies, R_inf, voigt_elements)
    return frequencies, Z, R_inf


@pytest.fixture
def warburg_data(frequencies):
    """Generate Warburg (Randles) circuit data."""
    R_inf = 50
    Z = generate_warburg_impedance(frequencies, R_inf, R_ct=500, C_dl=1e-5, sigma_w=100)
    return frequencies, Z, R_inf


# =============================================================================
# Tests
# =============================================================================

def test_gcv_returns_positive_lambda(voigt_data):
    """Test that GCV returns a positive lambda value."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = build_drt_matrices(frequencies, Z, R_inf)

    lambda_gcv, gcv_score = find_optimal_lambda_gcv(A, b, L)

    assert lambda_gcv > 0, "Lambda should be positive"
    assert np.isfinite(lambda_gcv), "Lambda should be finite"
    assert np.isfinite(gcv_score), "GCV score should be finite"


def test_hybrid_returns_positive_lambda(voigt_data):
    """Test that Hybrid method returns a positive lambda value."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = build_drt_matrices(frequencies, Z, R_inf)

    lambda_hybrid, hybrid_score, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"
    assert 'method_used' in diag, "Diagnostics should include method_used"


def test_gcv_and_hybrid_similar_for_clean_voigt(voigt_data):
    """Test that GCV and Hybrid give similar results for clean Voigt data."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = build_drt_matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    lambda_hybrid, _, _ = find_optimal_lambda_hybrid(A, b, L)

    # For clean data, both methods should be within 1 order of magnitude
    ratio = lambda_hybrid / lambda_gcv
    assert 0.1 < ratio < 10, f"Methods differ too much: ratio={ratio:.2f}"


def test_hybrid_works_for_noisy_data(voigt_data):
    """Test that Hybrid method works for noisy data."""
    frequencies, Z, R_inf = voigt_data

    np.random.seed(42)
    noise = 0.02 * np.abs(Z) * (np.random.randn(len(Z)) + 1j * np.random.randn(len(Z)))
    Z_noisy = Z + noise

    A, b, L = build_drt_matrices(frequencies, Z_noisy, R_inf)

    lambda_hybrid, _, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"


def test_hybrid_works_for_warburg_data(warburg_data):
    """Test that Hybrid method works for Warburg (diffusion) data."""
    frequencies, Z, R_inf = warburg_data
    A, b, L = build_drt_matrices(frequencies, Z, R_inf)

    lambda_hybrid, _, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"


def test_lambda_in_reasonable_range(voigt_data):
    """Test that lambda values are in a reasonable range."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = build_drt_matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    lambda_hybrid, _, _ = find_optimal_lambda_hybrid(A, b, L)

    # Lambda should typically be between 1e-6 and 1e2
    assert 1e-8 < lambda_gcv < 1e3, f"GCV lambda out of range: {lambda_gcv}"
    assert 1e-8 < lambda_hybrid < 1e3, f"Hybrid lambda out of range: {lambda_hybrid}"
