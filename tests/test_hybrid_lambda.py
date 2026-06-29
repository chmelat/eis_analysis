#!/usr/bin/env python3
"""Test hybrid lambda selection (GCV + L-curve) for DRT analysis."""

import numpy as np
import pytest
from eis_analysis.drt.gcv import (
    compute_gcv_score,
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
)
from eis_analysis.drt.core import _build_drt_matrices


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


def _matrices(frequencies, Z, R_inf):
    """A, b, L from the production matrix builder (no test re-implementation)."""
    m = _build_drt_matrices(frequencies, Z, R_inf)
    return m.A, m.b, m.L


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
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_gcv, gcv_score = find_optimal_lambda_gcv(A, b, L)

    assert lambda_gcv > 0, "Lambda should be positive"
    assert np.isfinite(lambda_gcv), "Lambda should be finite"
    assert np.isfinite(gcv_score), "GCV score should be finite"


def test_hybrid_returns_positive_lambda(voigt_data):
    """Test that Hybrid method returns a positive lambda value."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_hybrid, hybrid_score, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"
    assert 'method_used' in diag, "Diagnostics should include method_used"


def test_gcv_and_hybrid_similar_for_clean_voigt(voigt_data):
    """Test that GCV and Hybrid give similar results for clean Voigt data."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

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

    A, b, L = _matrices(frequencies, Z_noisy, R_inf)

    lambda_hybrid, _, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"


def test_hybrid_works_for_warburg_data(warburg_data):
    """Test that Hybrid method works for Warburg (diffusion) data."""
    frequencies, Z, R_inf = warburg_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_hybrid, _, diag = find_optimal_lambda_hybrid(A, b, L)

    assert lambda_hybrid > 0, "Lambda should be positive"
    assert np.isfinite(lambda_hybrid), "Lambda should be finite"


def test_lambda_in_reasonable_range(voigt_data):
    """Test that lambda values are in a reasonable range."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    lambda_hybrid, _, _ = find_optimal_lambda_hybrid(A, b, L)

    # Lambda should typically be between 1e-6 and 1e2
    assert 1e-8 < lambda_gcv < 1e3, f"GCV lambda out of range: {lambda_gcv}"
    assert 1e-8 < lambda_hybrid < 1e3, f"Hybrid lambda out of range: {lambda_hybrid}"


# =============================================================================
# compute_gcv_score correctness (F13)
# =============================================================================

def test_gcv_score_finite_positive(voigt_data):
    """GCV score is finite and positive across the search range."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    for lam in np.logspace(-5, 0, 12):
        score = compute_gcv_score(lam, A, b, L)
        assert np.isfinite(score), f"GCV score not finite at lambda={lam:.1e}"
        assert score > 0, f"GCV score not positive at lambda={lam:.1e}"


def test_gcv_score_has_minimum(voigt_data):
    """The lambda selected by GCV scores no worse than the range endpoints.

    Verifies find_optimal_lambda_gcv genuinely minimizes compute_gcv_score
    (the selector and the score function are consistent), not just returns
    something in range.
    """
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    score_opt = compute_gcv_score(lambda_gcv, A, b, L)
    score_lo = compute_gcv_score(1e-5, A, b, L)
    score_hi = compute_gcv_score(1.0, A, b, L)

    assert score_opt <= score_lo + 1e-12, "selected lambda worse than lower edge"
    assert score_opt <= score_hi + 1e-12, "selected lambda worse than upper edge"
