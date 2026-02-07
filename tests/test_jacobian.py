#!/usr/bin/env python3
"""Systematic tests: analytic Jacobian vs numerical central differences.

Three levels:
1. Individual elements (R, C, L, Q, W, Wo, K, G)
2. Simple compositions (Series, Parallel, nested)
3. Complex nested circuits (multiple levels of nesting)
"""

import numpy as np
import pytest
from eis_analysis.fitting import R, C, L, Q, W, Wo, K, G
from eis_analysis.fitting.jacobian import circuit_jacobian


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def numerical_jacobian(circuit, freq, params, eps=1e-7):
    """Central-difference numerical Jacobian for any circuit/element.

    Uses relative step sizing: h = eps * |param| (or eps if param is 0).

    Returns
    -------
    dZ_num : ndarray of complex, shape (n_freq, n_params)
    """
    params = list(params)
    n_params = len(params)
    n_freq = len(freq)
    dZ = np.zeros((n_freq, n_params), dtype=complex)

    for j in range(n_params):
        p_fwd = params.copy()
        p_bwd = params.copy()
        h = eps * abs(params[j]) if params[j] != 0 else eps
        p_fwd[j] += h
        p_bwd[j] -= h
        Z_fwd = circuit.impedance(freq, p_fwd)
        Z_bwd = circuit.impedance(freq, p_bwd)
        dZ[:, j] = (Z_fwd - Z_bwd) / (2 * h)

    return dZ


def assert_jacobian_close(circuit, freq, params, tol=1e-5):
    """Assert analytic Jacobian matches numerical for every column."""
    params = list(params)
    Z_anal, dZ_anal = circuit_jacobian(circuit, freq, params)
    dZ_num = numerical_jacobian(circuit, freq, params)

    for j in range(dZ_anal.shape[1]):
        col_anal = dZ_anal[:, j]
        col_num = dZ_num[:, j]
        scale = np.max(np.abs(col_num))
        if scale < 1e-30:
            # Both should be near zero
            assert np.max(np.abs(col_anal)) < 1e-20, (
                f"Column {j}: numerical ~0 but analytic is not"
            )
            continue
        rel_err = np.max(np.abs(col_anal - col_num)) / scale
        assert rel_err < tol, (
            f"Column {j}: relative error {rel_err:.2e} exceeds tolerance {tol}"
        )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def freq():
    """50 points from 0.1 Hz to 100 kHz (log-spaced)."""
    return np.logspace(-1, 5, 50)


# ---------------------------------------------------------------------------
# Level 1: Individual elements
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("element, params", [
    (R(100), [100.0]),
    (C(1e-6), [1e-6]),
    (L(1e-3), [1e-3]),
    (Q(1e-5, 0.85), [1e-5, 0.85]),
    (W(50), [50.0]),
    (Wo(100, 0.01), [100.0, 0.01]),
    (K(500, 1e-4), [500.0, 1e-4]),
    (G(100, 1e-3), [100.0, 1e-3]),
], ids=["R", "C", "L", "Q", "W", "Wo", "K", "G"])
def test_element_jacobian(freq, element, params):
    """Analytic Jacobian of each element matches numerical."""
    assert_jacobian_close(element, freq, params)


# ---------------------------------------------------------------------------
# Level 2: Simple compositions
# ---------------------------------------------------------------------------

def test_series_R_C(freq):
    """Series: R - C."""
    circuit = R(100) - C(1e-6)
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


def test_parallel_R_C(freq):
    """Parallel: R | C."""
    circuit = R(100) | C(1e-6)
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


def test_randles_base(freq):
    """R - (R | C) -- basic Randles structure."""
    circuit = R(10) - (R(100) | C(1e-6))
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


# ---------------------------------------------------------------------------
# Level 3: Complex nested circuits
# ---------------------------------------------------------------------------

def test_two_rc_series(freq):
    """R - (R | C) - (R | C) -- two RC arcs in series."""
    circuit = R(10) - (R(100) | C(1e-6)) - (R(200) | C(1e-5))
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


def test_randles_cpe_warburg(freq):
    """R - (R | Q) - W -- Randles with CPE and Warburg."""
    circuit = R(10) - (R(100) | Q(1e-5, 0.85)) - W(50)
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


def test_nested_series_in_parallel(freq):
    """R - ((R | C) - (R | C)) | C -- nested series inside parallel."""
    inner = (R(100) | C(1e-6)) - (R(200) | C(1e-5))
    circuit = R(10) - (inner | C(1e-4))
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)


def test_three_level_nesting(freq):
    """R - (R | (C - (R | C))) -- three levels of nesting."""
    innermost = R(50) | C(1e-5)
    mid = C(1e-6) - innermost
    circuit = R(10) - (R(100) | mid)
    params = circuit.get_all_params()
    assert_jacobian_close(circuit, freq, params)
