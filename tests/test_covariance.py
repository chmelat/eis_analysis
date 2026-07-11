#!/usr/bin/env python3
"""Tests for covariance computation (fitting/covariance.py).

Focus: condition_number reports cond(J^T J) = cond(J)^2 of the *column-scaled*
Jacobian, and the is_well_conditioned threshold is applied to that
scale-invariant value (audit finding #3, corrected for column scaling).

Column scaling is essential for EIS: parameters span many orders of magnitude
(R ~ 1e5, Q ~ 1e-6, n ~ 0.5), so the raw Jacobian columns differ by >10 decades
and cond(J) exceeds 1e10 purely from those units -- a scaling artifact that must
NOT be reported as ill-conditioning or rank deficiency.
"""

import numpy as np
import pytest

from eis_analysis.fitting.covariance import compute_covariance_matrix


def _equal_norm_correlated_J(d):
    # Two unit-norm columns with Gram off-diagonal d, so column scaling is a
    # no-op and cond(J^T J) = (1 + d) / (1 - d) reflects genuine correlation.
    a = np.sqrt((1.0 + d) / 2.0)
    b = np.sqrt((1.0 - d) / 2.0)
    return np.column_stack([[a, b], [a, -b]])


def test_condition_number_is_cond_of_JtJ():
    # Equal-norm correlated columns, d = 0.5 -> cond(J^T J) = 1.5 / 0.5 = 3.
    J = _equal_norm_correlated_J(0.5)
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number == pytest.approx(3.0)


def test_scale_disparate_orthogonal_is_well_conditioned():
    # Orthogonal columns differing by 6 decades in units (cond of RAW J^T J is
    # 1e12) are perfectly identifiable: after column scaling they are the
    # identity, so the problem is well-conditioned and both stderrs are finite.
    # This is the EIS regression: unit disparity must not be read as ill
    # conditioning.
    J = np.diag([1.0, 1e-6])
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number == pytest.approx(1.0)
    assert result.is_well_conditioned
    assert np.all(np.isfinite(result.stderr))
    # cov = s^2 (J^T J)^{-1} = s^2 diag([1, 1e12]); the small-sensitivity
    # parameter gets a large but finite (estimable) stderr.
    assert result.stderr[1] / result.stderr[0] == pytest.approx(1e6, rel=1e-6)


def test_genuine_ill_conditioning_flagged():
    # Correlation-driven ill conditioning that column scaling cannot remove:
    # d -> 1 gives cond(J^T J) = (1 + d) / (1 - d) > 1e10.
    d = 1.0 - 2e-11
    J = _equal_norm_correlated_J(d)
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number > 1e10
    assert not result.is_well_conditioned


def test_well_conditioned_true_for_benign_jacobian():
    # Equal-norm, weakly correlated columns -> cond(J^T J) = 3, below 1e10.
    J = _equal_norm_correlated_J(0.5)
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.is_well_conditioned


def _rank_deficient_jacobian():
    # 3 columns, col2 = col0 + col1 -> rank 2, one zero singular value.
    c0 = np.array([1.0, 0.0, 0.0, 1.0])
    c1 = np.array([0.0, 1.0, 0.0, 1.0])
    return np.column_stack([c0, c1, c0 + c1])


def test_rank_deficient_returns_inf():
    J = _rank_deficient_jacobian()
    r = np.array([0.1, 0.2, 0.1, 0.05])
    result = compute_covariance_matrix(J, r, n_params=3)
    assert result.rank == 2
    assert not result.is_well_conditioned
    assert np.all(np.isinf(result.stderr))
    assert np.all(np.isinf(np.diag(result.cov)))
    assert "Rank-deficient" in result.warning_message


def test_rank_deficient_with_fixed_params():
    # Free block (columns) is rank-deficient; one parameter is fixed.
    J = _rank_deficient_jacobian()  # 3 free columns, rank 2
    r = np.array([0.1, 0.2, 0.1, 0.05])
    fixed = [True, False, False, False]  # 4 params total, 1 fixed, 3 free
    result = compute_covariance_matrix(J, r, n_params=4, fixed_params=fixed)
    assert len(result.stderr) == 4
    assert result.stderr[0] == 0.0          # fixed -> known exactly
    assert np.all(np.isinf(result.stderr[1:]))  # free, rank-deficient


def test_full_rank_stderr_finite():
    # Full rank -> covariance estimable, finite standard errors.
    J = np.diag([2.0, 0.5])
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert np.all(np.isfinite(result.stderr))


def test_integration_real_fit_has_finite_stderr():
    """Regression (v0.17.1): a good fit of a realistic circuit must yield
    finite standard errors.

    Replicates the demo scenario Rs - (R0||Q0) - (R1||Q1): the parameters span
    >10 decades (R ~ 1e5, Q ~ 1e-6, n ~ 0.5), so the raw Jacobian columns are
    scale-disparate and the pre-fix rank test on the unscaled Jacobian reported
    rank 5/7 -> all stderr inf, despite every parameter being identifiable.
    """
    from eis_analysis.fitting import R, Q, fit_equivalent_circuit

    circuit = R(10) - (R(1e5) | Q(1e-6, 0.6)) - (R(8e5) | Q(3e-5, 0.43))
    freq = np.logspace(5, -2, 70)
    true_params = [10.0, 1e5, 1e-6, 0.6, 8e5, 3e-5, 0.43]
    Z_true = circuit.impedance(freq, true_params)

    rng = np.random.RandomState(42)
    noise = 0.01 * np.abs(Z_true) * (rng.randn(len(freq)) + 1j * rng.randn(len(freq)))
    Z_noisy = Z_true + noise

    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit, plot=False)

    assert result.fit_error_rel < 5.0, \
        f"Fit itself failed (error {result.fit_error_rel:.1f}%), test is inconclusive"
    assert np.all(np.isfinite(result.params_stderr)), \
        f"Good fit reported non-finite stderr: {result.params_stderr}"
    n_params = len(result.params_opt)
    assert result.diagnostics.covariance_rank == n_params, \
        f"Full-rank problem reported rank {result.diagnostics.covariance_rank}/{n_params}"
    # Stderr must also be informative, not just finite: well below the
    # parameter magnitudes for 1% noise on 70 points.
    rel_err = result.params_stderr / np.abs(result.params_opt)
    assert np.all(rel_err < 1.0), f"Implausibly large relative stderr: {rel_err}"
