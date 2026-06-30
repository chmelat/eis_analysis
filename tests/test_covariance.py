#!/usr/bin/env python3
"""Tests for covariance computation (fitting/covariance.py).

Focus: condition_number reports cond(J^T J) = cond(J)^2 and the
is_well_conditioned threshold is applied to that value (audit finding #3).

Diagonal Jacobians are used so the singular values (hence cond(J)) are known
exactly.
"""

import numpy as np
import pytest

from eis_analysis.fitting.covariance import compute_covariance_matrix


def test_condition_number_is_cond_of_JtJ():
    # Singular values {2.0, 0.5} -> cond(J) = 4 -> cond(J^T J) = 16.
    J = np.diag([2.0, 0.5])
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number == pytest.approx(16.0)


def test_well_conditioned_threshold_on_JtJ():
    # cond(J) = 1e6 (< 1e10) but cond(J^T J) = 1e12 (> 1e10): the problem is
    # NOT well-conditioned for the covariance inverse.
    J = np.diag([1.0, 1e-6])
    r = np.array([0.1, 0.2])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number == pytest.approx(1e12, rel=1e-6)
    assert not result.is_well_conditioned


def test_well_conditioned_true_for_benign_jacobian():
    # cond(J) = 4 -> cond(J^T J) = 16, comfortably below 1e10.
    J = np.diag([2.0, 0.5])
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
