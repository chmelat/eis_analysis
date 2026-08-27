"""Tests for the DE search space (fitting/diffevo.py).

Regression tests (2026-08-26): DE sampled its population uniformly over bounds
spanning up to 14 decades (R: 1e-4..1e10, C: 1e-15..1e-1), so nearly every
member was drawn from the top decade. The parallel branches were then shorted,
every member predicted the series resistance alone, and the population energies
were so nearly equal that scipy's convergence test (std <= tol * |mean|) fired
after the first generation: on a high-impedance sample DE reported "converged"
after 1 iteration at 99.97% error and the whole fit came from the local
refinement. DE now searches log10 of those parameters.
"""

import pickle

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.fitting import R, Q, CC
from eis_analysis.fitting.bounds import generate_simple_bounds, log_scale_ci_mask
from eis_analysis.fitting.diagnostics import compute_weights
from eis_analysis.fitting.diffevo import (
    _DECostFunction,
    _LogSpaceCost,
    _to_linear,
    fit_circuit_diffevo,
)

FREQ = np.logspace(5, -2, 40)
# High-impedance oxide sample - the regime where the linear population collapses
TRUE = [7.0, 3e8, 3.6e-7, 0.75]  # Rs, R_ct, Q, n


def make_circuit(values=(100.0, 100.0, 1e-4, 0.8)):
    """Fresh Rs-(R||Q) circuit at the element defaults, as `R()-(R()|Q())` gives.

    The starting point matters: an initial guess that already describes the
    spectrum hides the bug, because scipy seeds it into the population as x0.
    These are the values a user gets from the CLI without explicit numbers -
    six decades below this sample's R_ct, i.e. on the degenerate plateau.
    (fit_circuit_diffevo mutates the circuit via update_params, hence a fresh
    one per call.)
    """
    rs, rct, q, n = values
    return R(rs) - (R(rct) | Q(q, n))


def true_impedance():
    return make_circuit(TRUE).impedance(FREQ, TRUE)


# =============================================================================
# Search-space transform
# =============================================================================

def test_log_mask_leaves_cpe_exponent_linear():
    """R, Q, C are searched in log space; the CPE exponent n is not."""
    lower = [1e-4, 1e-12, 0.3, 1e-15]   # R, Q, n, C
    upper = [1e10, 1e-1, 1.0, 1e-1]
    assert log_scale_ci_mask(lower, upper) == [True, True, False, True]


def test_log_mask_leaves_cole_cole_exponent_linear():
    """CC capacitances and tau are searched in log space; alpha is not."""
    lower, upper = generate_simple_bounds(CC().get_param_labels())
    assert log_scale_ci_mask(lower, upper) == [True, True, True, False]


def test_to_linear_transforms_only_masked_entries():
    mask = np.array([True, False, True])
    out = _to_linear([3.0, 0.75, -6.0], mask)
    np.testing.assert_allclose(out, [1e3, 0.75, 1e-6])


def test_to_linear_roundtrip():
    mask = np.array([True, True, False, True])
    values = np.array([1e-4, 5e7, 0.62, 2.5e-9])
    search = np.where(mask, np.log10(values), values)
    np.testing.assert_allclose(_to_linear(search, mask), values, rtol=1e-12)


def test_log_space_cost_is_picklable():
    """The wrapper exists to survive pickling for --de-workers > 1."""
    Z = true_impedance()
    mask = np.array([True, True, True, False])  # Rs, R_ct, Q log; n linear
    cost = _LogSpaceCost(_DECostFunction(make_circuit(), FREQ, Z, np.ones(len(Z))),
                         mask)
    restored = pickle.loads(pickle.dumps(cost))
    point = [np.log10(7.0), np.log10(3e8), np.log10(3.6e-7), 0.75]
    assert restored(point) == cost(point)


# =============================================================================
# End-to-end: DE must explore instead of converging on the first generation
# =============================================================================

def test_de_alone_finds_basin_on_high_impedance_data():
    """The DE stage itself must reach the data, not leave it to the refinement."""
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=300)
    plt.close('all')
    diag = result.diagnostics

    assert diag.de_iterations > 1        # used to stop early on the plateau
    assert diag.de_error < 30.0          # used to be 97.7% with linear sampling


def test_log_searched_params_reported():
    """Diagnostics name which parameters were searched logarithmically."""
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=100)
    plt.close('all')

    assert result.diagnostics.log_search_params == ['R0', 'R1', 'Q0']  # not n0


def test_de_result_x_is_in_linear_space():
    """de_result.x must come back as physical parameters, not log10 of them."""
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=100)
    plt.close('all')

    w = compute_weights(Z, 'modulus')
    cost = _DECostFunction(make_circuit(), FREQ, Z, w)
    assert cost(result.de_result.x) == pytest.approx(result.diagnostics.de_cost,
                                                     rel=1e-9)


def test_fixed_param_excluded_from_log_search():
    """The mask lives in free-parameter space; labels must still line up."""
    Z = true_impedance()
    circuit = R("7") - (R(1e5) | Q(1e-6, 0.9))
    result, _, _ = fit_circuit_diffevo(circuit, FREQ, Z, seed=42, maxiter=100)
    plt.close('all')

    assert result.diagnostics.log_search_params == ['R1', 'Q0']
    assert result.best_result.params_opt[0] == 7.0


# =============================================================================
# Warning when the global stage contributed nothing
# =============================================================================

def test_stalled_de_warns():
    """maxiter=1 leaves DE far from the data - the fit is then the local run."""
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=1)
    plt.close('all')

    assert any('Global search contributed nothing' in w
               for w in result.best_result.all_warnings)


def test_healthy_de_does_not_warn():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=300)
    plt.close('all')

    assert not any('Global search contributed nothing' in w
                   for w in result.best_result.all_warnings)
