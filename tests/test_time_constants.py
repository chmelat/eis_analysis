"""
Tests for fitting/time_constants.py against its CONTRACT (cold review,
doc/CONTRACT_REVIEW.md, round 2).
"""

import numpy as np

from eis_analysis.fitting.circuit_elements import C, G, Q, R
from eis_analysis.fitting.result import FitResult
from eis_analysis.fitting.time_constants import arc_resistance_sums, arc_time_constants


def _result(circuit, params, cov):
    params = np.asarray(params, dtype=float)
    return FitResult(circuit=circuit, params_opt=params,
                     params_stderr=np.zeros_like(params), fit_error_rel=0.0,
                     cov=np.asarray(cov, dtype=float),
                     param_labels=circuit.get_param_labels())


def test_tau_at_correlation_minus_one_is_a_rounding_zero():
    """Equal relative errors of R and C at correlation -1 fix tau exactly:
    stderr 0, as for the sum, where tau gave None in ~40 % of cases (F2)."""
    rng = np.random.default_rng(0)
    for _ in range(200):
        R_val, C_val = 10 ** rng.uniform(0, 5), 10 ** rng.uniform(-12, -3)
        rel = 10 ** rng.uniform(-4, -1)
        s_R, s_C = rel * R_val, rel * C_val
        cov = [[s_R ** 2, -s_R * s_C], [-s_R * s_C, s_C ** 2]]
        arc, = arc_time_constants(_result(R(1) | C(1), [R_val, C_val], cov))
        assert arc.tau_stderr is not None and arc.tau_stderr < 1e-6 * arc.tau


def test_non_finite_parameter_gives_no_tau_and_no_sum():
    """n = inf read as tau = 1 s with stderr 0 (F3); G = inf as a finite
    1/G = 0 with its variance erased (F7); a huge tau with a large relative
    error gave tau_stderr inf (/code-review)."""
    arc_q = R(1) | Q(1, 0.8)
    assert arc_time_constants(_result(arc_q, [100, 1e-5, np.inf], np.eye(3) * 1e-30)) == []

    chain = R(1) - (G(1) | C(1)) - (G(1) | C(1))
    result = _result(chain, [1, np.inf, 1e-5, 0.01, 1e-4], np.diag([0, 1e-4, 0, 1e-4, 0]))
    assert arc_resistance_sums(result) == []

    arc, = arc_time_constants(_result(R(1) | C(1), [1e200, 1e100], np.diag([0, 1e218])))
    assert arc.tau_stderr is None and arc.stderr_note
