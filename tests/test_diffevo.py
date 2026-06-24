#!/usr/bin/env python3
"""Unit and integration tests for Differential Evolution circuit fitting
(fitting/diffevo.py).

Covers:
- DE_STRATEGIES mapping and fallback
- _DECostFunction: free/fixed parameter reconstruction, scalar cost,
  picklability (its stated purpose for workers > 1)
- fit_circuit_diffevo: parameter recovery, return contract, reproducibility
  via the seed argument, strategy/jacobian selection, fixed parameters,
  refinement choice, and diagnostics population

DE is stochastic by default (seed=None). Integration tests pass an explicit
`seed` so they are deterministic; recovery uses noise-free synthetic data.
"""

import pickle

import numpy as np
import pytest

# Suppress matplotlib GUI (fit_circuit_diffevo builds a figure)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.fitting import R, C
from eis_analysis.fitting.circuit import FitResult
from eis_analysis.fitting.diffevo import (
    fit_circuit_diffevo,
    _DECostFunction,
    DE_STRATEGIES,
    DiffEvoResult,
)

FREQ = np.logspace(5, -1, 25)
TRUE = [100.0, 5000.0, 1e-6]  # Rs, R_ct, C_dl


def make_circuit(values=(120.0, 4000.0, 2e-6)):
    """Fresh Rs-(R||C) circuit (fit_circuit_diffevo mutates it via update_params)."""
    rs, rct, cdl = values
    return R(rs) - (R(rct) | C(cdl))


def true_impedance():
    """Noise-free spectrum generated from the TRUE parameters."""
    return make_circuit(TRUE).impedance(FREQ, TRUE)


# =============================================================================
# DE_STRATEGIES
# =============================================================================

def test_de_strategies_mapping():
    assert DE_STRATEGIES == {
        1: 'randtobest1bin',
        2: 'best1bin',
        3: 'rand1bin',
    }


def test_de_strategies_unknown_falls_back_to_default():
    assert DE_STRATEGIES.get(99, 'randtobest1bin') == 'randtobest1bin'


# =============================================================================
# _DECostFunction (deterministic units)
# =============================================================================

def test_costfunction_reconstruct_no_fixed_returns_free_unchanged():
    Z = true_impedance()
    cost = _DECostFunction(make_circuit(), FREQ, Z, np.ones(len(Z)))
    assert cost._reconstruct_params([1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]


def test_costfunction_reconstruct_inserts_fixed_values():
    Z = true_impedance()
    cost = _DECostFunction(
        make_circuit(), FREQ, Z, np.ones(len(Z)),
        fixed_params=[True, False, False],
        full_initial_guess=[10.0, 0.0, 0.0],
    )
    # Free vector has 2 entries; fixed index 0 takes the initial-guess value.
    assert cost._reconstruct_params([5000.0, 1e-6]) == [10.0, 5000.0, 1e-6]


def test_costfunction_zero_at_true_params():
    Z = true_impedance()
    cost = _DECostFunction(make_circuit(), FREQ, Z, np.ones(len(Z)))
    assert cost(TRUE) == pytest.approx(0.0, abs=1e-18)


def test_costfunction_positive_for_wrong_params():
    Z = true_impedance()
    cost = _DECostFunction(make_circuit(), FREQ, Z, np.ones(len(Z)))
    assert cost([200.0, 5000.0, 1e-6]) > 0.0


def test_costfunction_is_picklable():
    # The class exists to be picklable for workers > 1.
    Z = true_impedance()
    cost = _DECostFunction(make_circuit(), FREQ, Z, np.ones(len(Z)))
    restored = pickle.loads(pickle.dumps(cost))
    assert restored(TRUE) == cost(TRUE)
    assert restored([200.0, 5000.0, 1e-6]) == cost([200.0, 5000.0, 1e-6])


# =============================================================================
# fit_circuit_diffevo - recovery and contract
# =============================================================================

def test_recovers_parameters_on_noise_free_data():
    Z = true_impedance()
    result, _, fig = fit_circuit_diffevo(
        make_circuit(), FREQ, Z, seed=42, maxiter=300
    )
    try:
        recovered = np.asarray(result.best_result.params_opt, dtype=float)
        np.testing.assert_allclose(recovered, TRUE, rtol=0.05)
        assert result.final_error < 1.0  # percent
    finally:
        plt.close('all')


def test_return_contract():
    Z = true_impedance()
    result, Z_fit, fig = fit_circuit_diffevo(
        make_circuit(), FREQ, Z, seed=42, maxiter=200
    )
    try:
        assert isinstance(result, DiffEvoResult)
        assert isinstance(result.best_result, FitResult)
        assert Z_fit.shape == FREQ.shape
        assert np.all(np.isfinite(Z_fit))
        assert isinstance(fig, plt.Figure)
    finally:
        plt.close('all')


def test_seed_makes_runs_reproducible():
    Z = true_impedance()
    r1, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=7, maxiter=200)
    r2, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=7, maxiter=200)
    plt.close('all')
    np.testing.assert_array_equal(
        np.asarray(r1.best_result.params_opt, dtype=float),
        np.asarray(r2.best_result.params_opt, dtype=float),
    )


# =============================================================================
# Strategy and Jacobian selection
# =============================================================================

def test_strategy_reflected_in_diagnostics():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(
        make_circuit(), FREQ, Z, seed=42, maxiter=100, strategy=2
    )
    plt.close('all')
    assert result.strategy == 'best1bin'
    assert result.diagnostics.strategy == 'best1bin'


def test_unknown_strategy_falls_back_to_default():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(
        make_circuit(), FREQ, Z, seed=42, maxiter=100, strategy=99
    )
    plt.close('all')
    assert result.diagnostics.strategy == 'randtobest1bin'


def test_numeric_jacobian_path():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(
        make_circuit(), FREQ, Z, seed=42, maxiter=200,
        use_analytic_jacobian=False
    )
    plt.close('all')
    assert result.diagnostics.jacobian_type == 'numeric'
    assert result.final_error < 1.0  # still converges


# =============================================================================
# Fixed parameters (string parameter -> fixed)
# =============================================================================

def test_fixed_param_stays_constant():
    # Data generated with Rs = 10; Rs fixed via string "10".
    true = [10.0, 5000.0, 1e-6]
    Z = make_circuit(true).impedance(FREQ, true)
    circuit = R("10") - (R(5000) | C(1e-6))
    result, _, _ = fit_circuit_diffevo(circuit, FREQ, Z, seed=42, maxiter=300)
    plt.close('all')
    params = np.asarray(result.best_result.params_opt, dtype=float)
    assert params[0] == pytest.approx(10.0)            # fixed, unchanged
    np.testing.assert_allclose(params[1:], true[1:], rtol=0.05)  # free recovered


def test_fixed_param_indices_in_diagnostics():
    Z = make_circuit([10.0, 5000.0, 1e-6]).impedance(FREQ, [10.0, 5000.0, 1e-6])
    circuit = R("10") - (R(5000) | C(1e-6))
    result, _, _ = fit_circuit_diffevo(circuit, FREQ, Z, seed=42, maxiter=100)
    plt.close('all')
    assert result.diagnostics.n_fixed_params == 1
    assert result.diagnostics.fixed_param_indices == [0]


# =============================================================================
# Refinement choice and diagnostics population
# =============================================================================

def test_refinement_never_worse_than_de():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=200)
    plt.close('all')
    # The code keeps the better of DE / least_squares, so final <= de.
    assert result.final_error <= result.de_error + 1e-9


def test_diagnostics_populated():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=200)
    plt.close('all')
    d = result.diagnostics
    assert d.de_iterations > 0
    assert d.de_evaluations > 0
    assert d.total_evaluations > 0
    assert isinstance(d.de_converged, bool)
    assert len(d.initial_guess) == 3  # full parameter vector


def test_diffevoresult_fields():
    Z = true_impedance()
    result, _, _ = fit_circuit_diffevo(make_circuit(), FREQ, Z, seed=42, maxiter=200)
    plt.close('all')
    assert isinstance(result.improvement, float)
    assert result.n_evaluations > 0
    assert result.de_error >= 0.0
    assert result.final_error >= 0.0
