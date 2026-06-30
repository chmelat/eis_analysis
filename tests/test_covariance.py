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
