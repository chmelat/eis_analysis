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
    dof = 47

    ci_low_95, ci_high_95 = compute_confidence_interval(params, stderr, dof, 0.95)

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
    dof = 47

    ci_low_95, ci_high_95 = compute_confidence_interval(params, stderr, dof, 0.95)
    ci_low_99, ci_high_99 = compute_confidence_interval(params, stderr, dof, 0.99)

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


def test_ci_dof_matches_covariance_residual_dof():
    """CI uses the residual dof of the variance estimate: 2*n_freq - n_free.

    Pre-fix the CI used n_freq - n_total, which is wrong for complex (split
    real/imag) residuals.
    """
    circuit = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 30)  # n_freq = 30
    true_params = [98.5, 4823.2, 8.7e-7]
    Z_true = circuit.impedance(freq, true_params)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq))
                                     + 1j * np.random.randn(len(freq)))
    result, _, _ = fit_equivalent_circuit(freq, Z_true + noise, circuit, plot=False)

    # 3 free params, 2*30 = 60 real residuals -> dof = 57.
    expected_dof = 2 * len(freq) - len(result.params_opt)
    assert result._dof == expected_dof

    ci_low, ci_high = result.params_ci_95
    # Recover the t-multiplier from a parameter with finite, positive stderr.
    # R and C are positive scale parameters, so the CI is log-space:
    # ci_high = p * exp(t * se / p) -> t = ln(ci_high / p) * p / se.
    i = int(np.argmax(np.isfinite(result.params_stderr) & (result.params_stderr > 0)))
    p, se = result.params_opt[i], result.params_stderr[i]
    t_used = np.log(ci_high[i] / p) * p / se
    assert t_used == pytest.approx(t.ppf(0.975, expected_dof), rel=1e-9)


def test_log_space_ci_positive_for_scale_params():
    """Regression: log-space CI of a positive scale parameter never crosses
    zero, even when stderr is comparable to the parameter value.

    Pre-fix the symmetric Wald interval p +/- t*se gave e.g.
    R = 14.2 +/- 7.4 -> CI [-0.43, 28.7] (negative resistance).
    """
    params = np.array([14.2])
    stderr = np.array([7.38])
    dof = 133

    ci_low, ci_high = compute_confidence_interval(
        params, stderr, dof, 0.95, log_scale=[True]
    )

    assert ci_low[0] > 0, f"Log-space CI lower bound must be positive, got {ci_low[0]}"
    assert ci_high[0] > params[0]
    # Multiplicative symmetry: ci_high / p == p / ci_low
    assert ci_high[0] / params[0] == pytest.approx(params[0] / ci_low[0], rel=1e-12)


def test_log_space_ci_converges_to_linear_for_small_stderr():
    """For se/p << 1 the log-space CI converges to the symmetric interval."""
    params = np.array([1e5])
    stderr = np.array([1e2])  # se/p = 0.1%
    dof = 100

    ci_lin_low, ci_lin_high = compute_confidence_interval(params, stderr, dof, 0.95)
    ci_log_low, ci_log_high = compute_confidence_interval(
        params, stderr, dof, 0.95, log_scale=[True]
    )

    # Agreement to first order in se/p (relative difference ~ (t*se/p)^2 / 2)
    assert ci_log_low[0] == pytest.approx(ci_lin_low[0], rel=1e-5)
    assert ci_log_high[0] == pytest.approx(ci_lin_high[0], rel=1e-5)


def test_log_space_mask_selects_per_parameter():
    """Only masked parameters get log-space CI; unmasked keep symmetric CI."""
    params = np.array([100.0, 0.6])   # scale param R, exponent n
    stderr = np.array([80.0, 0.05])
    dof = 50
    t_crit = t.ppf(0.975, dof)

    ci_low, ci_high = compute_confidence_interval(
        params, stderr, dof, 0.95, log_scale=[True, False]
    )

    # R (masked): positive, multiplicative
    assert ci_low[0] > 0
    # n (unmasked): exactly symmetric p +/- t*se
    assert ci_low[1] == pytest.approx(params[1] - t_crit * stderr[1], rel=1e-12)
    assert ci_high[1] == pytest.approx(params[1] + t_crit * stderr[1], rel=1e-12)


def test_fit_result_ci_log_scale_mask_from_bounds():
    """fit_equivalent_circuit builds the log-scale mask from bounds: R/C/Q
    get log-space CI (positive lower bound), CPE exponent n stays linear."""
    from eis_analysis.fitting import Q

    circuit = R(100) - (R(5000) | Q(1e-6, 0.8))
    freq = np.logspace(4, -1, 30)
    Z_true = circuit.impedance(freq, [100.0, 5000.0, 1e-6, 0.8])

    np.random.seed(7)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq))
                                     + 1j * np.random.randn(len(freq)))
    result, _, _ = fit_equivalent_circuit(freq, Z_true + noise, circuit, plot=False)

    # Labels: R, R, Q, n -> mask True, True, True, False
    assert result._ci_log_scale == [True, True, True, False]

    ci_low, _ = result.params_ci_95
    # Scale parameters (R, R, Q): CI strictly positive
    for i in (0, 1, 2):
        assert ci_low[i] > 0, f"Param {i}: log-space CI must be positive"


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


def test_backwards_compatibility_no_dof():
    """Test FitResult can be created without _dof (uses default)."""
    circuit = R(100) - (R(5000) | C(1e-6))

    result = FitResult(
        circuit=circuit,
        params_opt=np.array([100.0, 5000.0, 1e-6]),
        params_stderr=np.array([1.0, 50.0, 1e-8]),
        fit_error_rel=1.0,
        fit_error_abs=10.0,
        quality='good'
        # _dof omitted - uses default 0
    )

    assert result._dof == 0, "Default _dof should be 0"

    # CI should still be computable (with dof clamped to 1, very wide)
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
