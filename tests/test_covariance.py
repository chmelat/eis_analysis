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
    # The zero row leaves J^T J as it is and gives dof = 1 > 0.
    a = np.sqrt((1.0 + d) / 2.0)
    b = np.sqrt((1.0 - d) / 2.0)
    return np.column_stack([[a, b, 0.0], [a, -b, 0.0]])


def test_condition_number_is_cond_of_JtJ():
    # Equal-norm correlated columns, d = 0.5 -> cond(J^T J) = 1.5 / 0.5 = 3.
    J = _equal_norm_correlated_J(0.5)
    r = np.array([0.1, 0.2, 0.0])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number == pytest.approx(3.0)


def test_genuine_ill_conditioning_flagged():
    # Correlation-driven ill conditioning that column scaling cannot remove:
    # d -> 1 gives cond(J^T J) = (1 + d) / (1 - d) > 1e10.
    d = 1.0 - 2e-11
    J = _equal_norm_correlated_J(d)
    r = np.array([0.1, 0.2, 0.0])
    result = compute_covariance_matrix(J, r, n_params=2)
    assert result.condition_number > 1e10
    assert not result.is_well_conditioned


def test_well_conditioned_true_for_benign_jacobian():
    # Equal-norm, weakly correlated columns -> cond(J^T J) = 3, below 1e10.
    J = _equal_norm_correlated_J(0.5)
    r = np.array([0.1, 0.2, 0.0])
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


def test_partial_rank_deficiency_per_parameter():
    """Only parameters overlapping the null space get inf; the rest get
    finite stderr from the pseudo-inverse over the identifiable subspace.

    col1 = col2 (collinear pair), col0 orthogonal to both. The null vector
    (0, 1, -1)/sqrt(2) involves params 1 and 2 only, so param 0 must keep a
    finite, exact stderr.
    """
    c0 = np.array([1.0, 0.0, 0.0, 0.0])
    c1 = np.array([0.0, 1.0, 1.0, 0.0])
    J = np.column_stack([c0, c1, c1.copy()])
    r = np.array([0.1, 0.2, 0.1, 0.05])
    result = compute_covariance_matrix(J, r, n_params=3)

    assert result.rank == 2
    # Params 1, 2 span the null space -> non-identifiable
    assert np.isinf(result.stderr[1]) and np.isinf(result.stderr[2])
    assert np.all(np.isinf(result.cov[1, :])) and np.all(np.isinf(result.cov[:, 2]))
    # Param 0 is orthogonal to the null space -> exact finite stderr:
    # s^2 = rss/dof = 0.0625/1, unit column -> stderr = sqrt(0.0625) = 0.25
    assert result.stderr[0] == pytest.approx(0.25, rel=1e-12)
    assert result.cov[0, 0] == pytest.approx(0.0625, rel=1e-12)
    # Warning names the non-identifiable parameters
    assert "[1, 2]" in result.warning_message


def test_partial_rank_deficiency_with_fixed_params():
    """Non-identifiable indices in the warning are in FULL parameter space
    when fixed parameters shift the free-parameter indexing."""
    c0 = np.array([1.0, 0.0, 0.0, 0.0])
    c1 = np.array([0.0, 1.0, 1.0, 0.0])
    J = np.column_stack([c0, c1, c1.copy()])   # 3 free columns
    r = np.array([0.1, 0.2, 0.1, 0.05])
    fixed = [True, False, False, False]         # full params: 0 fixed, 1..3 free
    result = compute_covariance_matrix(J, r, n_params=4, fixed_params=fixed)

    assert result.stderr[0] == 0.0                     # fixed -> known exactly
    assert result.stderr[1] == pytest.approx(0.25)     # identifiable free
    assert np.isinf(result.stderr[2]) and np.isinf(result.stderr[3])
    assert "[2, 3]" in result.warning_message          # full-space indices


def test_partial_rank_deficiency_identifiable_variance_matches_pinv():
    """The identifiable parameter's variance is s^2 [pinv(J^T J)]_00.

    The test above has a unit column orthogonal to the others, where S = 1
    and the column scale is 1: a mutation run (doc/MUTATION_TESTING.md, D2)
    found every slip in this branch's arithmetic passing. Here c0 is
    correlated with c1, c2 = 5 c1, and the norms differ. e0 is orthogonal
    to the null space (0, 5, -1), so its variance is the same for any
    generalized inverse of J^T J: the raw pseudo-inverse is an exact
    reference for the scaled one, here 1 / (|c0|^2 - (c0.c1)^2 / |c1|^2).
    """
    c0 = np.array([3.0, 1.0, 0.0, 0.5, 0.0])
    c1 = np.array([0.0, 2.0, 1.0, 0.0, 1.0])
    J = np.column_stack([c0, c1, 5.0 * c1])
    r = np.array([0.1, -0.2, 0.3, 0.05, -0.1])
    result = compute_covariance_matrix(J, r, n_params=3)

    s2 = r @ r / (5 - 3)
    expected = s2 * np.linalg.pinv(J.T @ J)[0, 0]
    assert expected == pytest.approx(s2 / (10.25 - 4.0 / 6.0), rel=1e-12)
    assert result.rank == 2
    assert np.isinf(result.stderr[1]) and np.isinf(result.stderr[2])
    assert result.cov[0, 0] == pytest.approx(expected, rel=1e-10)
    assert result.stderr[0] == pytest.approx(np.sqrt(expected), rel=1e-10)


def test_full_rank_stderr_finite():
    # Full rank -> covariance estimable, finite standard errors.
    J = np.array([[2.0, 0.0], [0.0, 0.5], [0.0, 0.0]])  # 3 residuals, dof 1
    r = np.array([0.1, 0.2, 0.0])
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

    result, _ = fit_equivalent_circuit(freq, Z_noisy, circuit)

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
    # Scale parameters (all except the two CPE exponents n) get log-space
    # CIs, which must never cross zero.
    ci_low, _ = result.params_ci_95
    scale_idx = [i for i, is_log in enumerate(result._ci_log_scale) if is_log]
    assert len(scale_idx) == 5  # Rs, R0, Q0, R1, Q1
    for i in scale_idx:
        assert ci_low[i] > 0, f"Param {i}: scale-parameter CI crosses zero"


@pytest.mark.parametrize("n_res, n_free", [(2, 2), (2, 3), (0, 2)])
def test_too_few_residuals_gives_inf_stderr(n_res, n_free):
    """n_res <= n_free leaves no residual variance (dof <= 0): stderr inf, not
    a finite value from a dof clamped to 1 (cold review F1, F2, F5)."""
    rng = np.random.default_rng(0)
    J = rng.normal(size=(n_res, n_free))
    fixed = [False] * n_free + [True]
    result = compute_covariance_matrix(J, rng.normal(size=n_res), n_free + 1, fixed_params=fixed)
    assert np.all(np.isinf(result.stderr[:n_free])) and result.stderr[-1] == 0
    assert result.dof == n_res - n_free
    assert not result.is_well_conditioned and result.warning_message


@pytest.mark.parametrize("where", ["jacobian", "residuals"])
@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_non_finite_input_gives_no_covariance(where, bad):
    """NaN or inf in the input: cov None and every stderr inf, fixed ones too,
    where a NaN residual gave a NaN cov with no warning (cold review F3, F4)."""
    rng = np.random.default_rng(0)
    J, r = rng.normal(size=(6, 2)), rng.normal(size=6)
    (J if where == "jacobian" else r)[2] = bad
    result = compute_covariance_matrix(J, r, 3, fixed_params=[False, True, False])
    assert result.cov is None and np.all(np.isinf(result.stderr))
    assert "Non-finite" in result.warning_message
