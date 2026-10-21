"""Tests for the canonical order of interchangeable series blocks.

Equal blocks in series, two (R|Q) in R-(R|Q)-(R|Q), can swap values without
changing the impedance; DE returned them at random, least_squares in the
order of its start. The fit now reports them in ascending tau
(doc/archive/BLOCK_ORDER_PLAN.md).
"""

import numpy as np

import pytest

from eis_analysis.fitting import fit_equivalent_circuit
from eis_analysis.fitting.bounds import generate_simple_bounds
from eis_analysis.fitting.circuit_elements import C, K, Q, R, Wo
from eis_analysis.fitting.block_order import canonical_block_order

FREQ = np.logspace(-2, 6, 40)


def _order(circuit, lower=None, upper=None):
    params = circuit.get_all_params()
    if lower is None:
        lower, upper = generate_simple_bounds(circuit.get_param_labels())
    perm, free_perm, _ = canonical_block_order(circuit, params, circuit.get_all_fixed_params(),
                                               lower, upper)
    Z_before = circuit.impedance(FREQ, params)
    Z_after = circuit.impedance(FREQ, list(np.asarray(params)[perm]))
    np.testing.assert_allclose(Z_after, Z_before, rtol=1e-12)
    return perm.tolist(), free_perm.tolist()


@pytest.mark.parametrize('circuit, expected', [
    # slow arc written first -> swapped
    (R(10) - (R(1e4) | Q(1e-4, 0.8)) - (R(1e3) | Q(1e-7, 0.9)), [0, 4, 5, 6, 1, 2, 3]),
    # already ascending -> unchanged
    (R(10) - (R(1e3) | Q(1e-7, 0.9)) - (R(1e4) | Q(1e-4, 0.8)), [0, 1, 2, 3, 4, 5, 6]),
    # a fixed R makes the blocks different models
    (R(10) - (R("1e4") | Q(1e-4, 0.8)) - (R(1e3) | Q(1e-7, 0.9)), [0, 1, 2, 3, 4, 5, 6]),
    # different structure
    (R(10) - (R(1e4) | Q(1e-4, 0.8)) - (R(1e3) | C(1e-7)), [0, 1, 2, 3, 4, 5]),
    # no tau for (R|Wo)
    (R(10) - (R(1e4) | Wo(10, 1)) - (R(1e3) | Wo(10, 1e-3)), [0, 1, 2, 3, 4, 5, 6]),
    # a Series inside a Series is the same chain
    ((R(1e4) | Q(1e-4, 0.8)) - (R(10) - (R(1e3) | Q(1e-7, 0.9))), [4, 5, 6, 3, 0, 1, 2]),
    # (C|R) is the same block as (R|C); the values land on their own kind
    (R(10) - (R(1e4) | C(1e-4)) - (C(1e-7) | R(1e3)), [0, 4, 3, 2, 1]),
    # K by its own tau, three blocks
    (R(10) - K(1e3, 1.0) - K(1e2, 1e-3) - K(5, 1e-5), [0, 5, 6, 3, 4, 1, 2]),
])
def test_order(circuit, expected):
    perm, free_perm = _order(circuit)
    assert perm == expected
    if not any(circuit.get_all_fixed_params()):
        assert free_perm == expected


def test_different_bounds_not_swapped():
    circuit = R(10) - (R(1e4) | C(1e-4)) - (R(1e3) | C(1e-7))
    lower, upper = generate_simple_bounds(circuit.get_param_labels())
    upper[1] = 1e6       # only the first block's R
    assert _order(circuit, lower, upper)[0] == [0, 1, 2, 3, 4]


def test_lm_fit_same_from_either_order():
    """A start with the slow arc first gives the same fit, stderr and
    covariance as one in ascending order: the Jacobian columns move with
    the parameters."""
    true = R(10) - (R(1e3) | Q(1e-7, 0.9)) - (R(5e3) | Q(1e-4, 0.8))
    Z = true.impedance(FREQ, true.get_all_params())
    rng = np.random.default_rng(0)
    Z = Z + 0.005 * np.abs(Z) * (rng.normal(size=Z.size) + 1j * rng.normal(size=Z.size))

    fast_first = R(12) - (R(8e2) | Q(2e-7, 0.85)) - (R(6e3) | Q(2e-4, 0.85))
    slow_first = R(12) - (R(6e3) | Q(2e-4, 0.85)) - (R(8e2) | Q(2e-7, 0.85))
    a, _ = fit_equivalent_circuit(FREQ, Z, fast_first)
    b, _ = fit_equivalent_circuit(FREQ, Z, slow_first)
    assert a.block_order == list(range(7))
    assert b.block_order == [0, 4, 5, 6, 1, 2, 3]

    np.testing.assert_allclose(b.params_opt, a.params_opt, rtol=1e-4)
    np.testing.assert_allclose(b.params_stderr, a.params_stderr, rtol=1e-3)
    corr = [r.cov / np.outer(r.params_stderr, r.params_stderr) for r in (a, b)]
    np.testing.assert_allclose(corr[1], corr[0], atol=1e-3)
