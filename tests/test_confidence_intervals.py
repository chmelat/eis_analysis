#!/usr/bin/env python3
"""Test confidence intervals computation for circuit fitting."""

import numpy as np
import pytest
from scipy.stats import t
from dataclasses import replace

from eis_analysis.fitting import R, C, fit_equivalent_circuit
from eis_analysis.fitting.circuit import FitResult
from eis_analysis.fitting.covariance import compute_confidence_interval


def testcompute_confidence_interval_margins():
    """Test that CI margins are computed correctly using t-distribution."""
    params = np.array([100.0, 5000.0, 1e-6])
    stderr = np.array([1.0, 50.0, 1e-8])
    n_data = 50
    dof = n_data - len(params)

    ci_low_95, ci_high_95 = compute_confidence_interval(params, stderr, n_data, 0.95)

    # Verify against scipy t-distribution
    t_95 = t.ppf(0.975, dof)
    expected_margin = t_95 * stderr[0]
    actual_margin = ci_high_95[0] - params[0]

    assert np.isclose(expected_margin, actual_margin, rtol=1e-10), \
        f"Margin mismatch: expected {expected_margin}, got {actual_margin}"


def test_ci_99_wider_than_95():
    """Test that 99% CI is wider than 95% CI."""
    params = np.array([100.0, 5000.0, 1e-6])
    stderr = np.array([1.0, 50.0, 1e-8])
    n_data = 50

    ci_low_95, ci_high_95 = compute_confidence_interval(params, stderr, n_data, 0.95)
    ci_low_99, ci_high_99 = compute_confidence_interval(params, stderr, n_data, 0.99)

    width_95 = ci_high_95 - ci_low_95
    width_99 = ci_high_99 - ci_low_99

    assert np.all(width_99 > width_95), "99% CI should be wider than 95% CI"


def test_fit_result_ci_properties():
    """Test FitResult.params_ci_95 and params_ci_99 properties."""
    circuit = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 30)
    true_params = [98.5, 4823.2, 8.7e-7]
    Z_true = circuit.impedance(freq, true_params)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, plot=False)

    ci_low_95, ci_high_95 = result.params_ci_95
    ci_low_99, ci_high_99 = result.params_ci_99

    # CI should be arrays of same length as params
    assert len(ci_low_95) == len(result.params_opt)
    assert len(ci_high_95) == len(result.params_opt)

    # 99% CI should be wider
    width_95 = ci_high_95 - ci_low_95
    width_99 = ci_high_99 - ci_low_99
    assert np.all(width_99 > width_95)


def test_invalid_stderr_returns_inf_ci():
    """Test that invalid (inf) stderr returns +/-inf CI."""
    circuit = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 30)
    Z = circuit.impedance(freq, [100.0, 5000.0, 1e-6])

    result, _, _ = fit_equivalent_circuit(freq, Z, circuit, plot=False)

    # Replace stderr with inf
    result_bad = replace(result, params_stderr=np.full_like(result.params_opt, np.inf))
    ci_low, ci_high = result_bad.params_ci_95

    assert np.all(np.isinf(ci_low)) and np.all(ci_low < 0), "Expected -inf for low CI"
    assert np.all(np.isinf(ci_high)) and np.all(ci_high > 0), "Expected +inf for high CI"


def test_backwards_compatibility_no_n_data():
    """Test FitResult can be created without _n_data (uses default)."""
    circuit = R(100) - (R(5000) | C(1e-6))

    result = FitResult(
        circuit=circuit,
        params_opt=np.array([100.0, 5000.0, 1e-6]),
        params_stderr=np.array([1.0, 50.0, 1e-8]),
        fit_error_rel=1.0,
        fit_error_abs=10.0,
        quality='good'
        # _n_data omitted - uses default 0
    )

    assert result._n_data == 0, "Default _n_data should be 0"

    # CI should still be computable (with dof=1, very wide)
    ci_low, ci_high = result.params_ci_95
    assert len(ci_low) == 3
    assert len(ci_high) == 3


def test_t_distribution_conservative_for_small_n():
    """Test that t-distribution is more conservative than normal for small samples."""
    n_small = 10
    n_large = 100
    n_params = 3

    t_small = t.ppf(0.975, n_small - n_params)
    t_large = t.ppf(0.975, n_large - n_params)
    z_normal = 1.96

    # Small sample should have larger t-critical
    assert t_small > z_normal, "t-critical should be > 1.96 for small samples"
    assert t_small > t_large, "t-critical should decrease with more data"
    # Large sample should approach normal
    assert abs(t_large - z_normal) < 0.1, "t-critical should approach 1.96 for large samples"
