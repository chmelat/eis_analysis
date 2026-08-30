#!/usr/bin/env python3
"""Tests for AIC/BIC model selection (fitting/diagnostics.py) and the circuit
comparison built on it (cli/handlers/model_comparison.py).

The point of the feature is that a fit error cannot distinguish a real
improvement from fitting the noise: adding an element almost always lowers
the residual. The central test here is exactly that - on data generated from
a known circuit, the criterion must reject a superfluous extra branch even
though the bigger model fits better.

DE is stochastic by default, so every fit passes an explicit seed.
"""

import numpy as np
import pytest

# Suppress matplotlib GUI (the optimizers build figures)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.cli.handlers.model_comparison import (
    score_candidates,
    log_comparison,
)
from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting.circuit_elements import Q, R
from eis_analysis.fitting.diagnostics import compute_information_criteria
from eis_analysis.fitting.diffevo import fit_circuit_diffevo
from eis_analysis.fitting.multistart import fit_circuit_multistart
from eis_analysis.fitting.circuit import fit_equivalent_circuit

FREQ = np.logspace(-2, 5, 80)


def _synthesize(expression, noise=0.02, seed=0):
    """Impedance of `expression` with proportional Gaussian noise."""
    circuit = parse_circuit_expression(expression)
    Z = circuit.impedance(FREQ, circuit.get_all_params())
    if not noise:
        return Z
    rng = np.random.default_rng(seed)
    sigma = noise * np.abs(Z)
    return Z + rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)


def _fit(expression, Z):
    """Fit one candidate with DE and return its FitResult."""
    result, _, _ = fit_circuit_diffevo(
        parse_circuit_expression(expression), FREQ, Z, seed=0
    )
    plt.close('all')
    return result.best_result


# =============================================================================
# The criteria themselves
# =============================================================================

def test_criteria_match_the_closed_form():
    """AIC and BIC must equal n*ln(RSS/n) + penalty, computed by hand."""
    # Uniform weighting normalizes to w = 1, so RSS is the plain SSR and the
    # expected value can be written down without reproducing compute_weights.
    Z = np.array([3.0 + 4.0j, 6.0 + 8.0j])
    Z_fit = np.array([2.0 + 4.0j, 6.0 + 6.0j])
    # residuals: re 1, 0; im 0, 2  ->  RSS = 1 + 4 = 5
    rss, aic, bic = compute_information_criteria(Z, Z_fit, 'uniform', 3)

    n = 4  # 2 points x (Re, Im)
    assert rss == pytest.approx(5.0)
    assert aic == pytest.approx(n * np.log(5.0 / n) + 2 * 3)
    assert bic == pytest.approx(n * np.log(5.0 / n) + 3 * np.log(n))


def test_perfect_fit_stays_finite():
    """RSS = 0 is the correct limit but must not poison the delta arithmetic."""
    Z = np.array([1.0 + 1.0j, 2.0 + 2.0j])
    rss, aic, bic = compute_information_criteria(Z, Z.copy(), 'modulus', 2)

    assert rss == 0.0
    assert np.isfinite(aic) and np.isfinite(bic)


def test_weighting_changes_rss_but_not_the_ranking():
    """RSS depends on the weighting; the ordering of candidates must not.

    Criteria are only comparable within one weighting, so the guarantee that
    matters is internal consistency: whichever candidate wins under modulus
    weighting must also win under uniform weighting on this clean a case.
    """
    Z = _synthesize('R(10)-(R(200)|C(2e-5))')
    good = _fit('R(10)-(R(200)|C(2e-5))', Z)
    bad = _fit('R(10)-(C(2e-5))', Z)

    rankings = []
    for weighting in ('uniform', 'modulus'):
        scores = score_candidates(
            FREQ, Z,
            [('good', good), ('bad', bad)],
            weighting,
        )
        rankings.append(scores[0].expression)

    assert rankings == ['good', 'good']


# =============================================================================
# Model selection on real fits - the core of the feature
# =============================================================================

def test_bic_rejects_a_superfluous_branch():
    """On data from R-(R|C), BIC must not buy a second RC branch.

    The two-branch model has strictly more freedom and therefore a lower
    residual; if the comparison were made on fit error alone it would win.
    This is the test the whole feature exists for.
    """
    Z = _synthesize('R(10)-(R(200)|C(2e-5))')
    simple = _fit('R(10)-(R(200)|C(2e-5))', Z)
    complicated = _fit('R(10)-(R(150)|C(2e-5))-(R(50)|C(1e-6))', Z)

    # Precondition: the bigger model really does fit better, otherwise this
    # test would pass for the wrong reason.
    assert complicated.fit_error_rel < simple.fit_error_rel

    scores = score_candidates(
        FREQ, Z,
        [('simple', simple), ('complicated', complicated)],
        'modulus',
    )
    assert scores[0].expression == 'simple'
    assert scores[1].bic - scores[0].bic > 2.0


def test_bic_buys_a_cpe_when_the_data_has_one():
    """The penalty must not make the criterion reject every extra parameter.

    On data from a CPE circuit, the extra exponent earns its place: the
    complement of the test above.
    """
    Z = _synthesize('R(10)-(R(200)|Q(2e-5,0.85))')
    with_cpe = _fit('R(10)-(R(200)|Q(2e-5,0.85))', Z)
    with_c = _fit('R(10)-(R(200)|C(2e-5))', Z)

    scores = score_candidates(
        FREQ, Z,
        [('cpe', with_cpe), ('c', with_c)],
        'modulus',
    )
    assert scores[0].expression == 'cpe'
    assert scores[1].bic - scores[0].bic > 10.0


# =============================================================================
# n_free_params - the k that goes into the penalty
# =============================================================================

@pytest.mark.parametrize('optimizer', ['single', 'de', 'multistart'])
def test_n_free_params_counts_free_parameters(optimizer):
    """Every optimizer must report k, not leave it at the default 0."""
    Z = _synthesize('R(10)-(R(200)|Q(2e-5,0.85))', noise=0.0)
    circuit = R(10) - (R(200) | Q(2e-5, 0.85))

    if optimizer == 'single':
        result = fit_equivalent_circuit(FREQ, Z, circuit, plot=False)[0]
    elif optimizer == 'de':
        result = fit_circuit_diffevo(circuit, FREQ, Z, maxiter=30, seed=0)[0].best_result
    else:
        result = fit_circuit_multistart(circuit, FREQ, Z, n_restarts=2)[0].best_result
    plt.close('all')

    assert result.n_free_params == 4


@pytest.mark.parametrize('optimizer', ['single', 'de', 'multistart'])
def test_n_free_params_excludes_fixed_parameters(optimizer):
    """A fixed parameter costs no degrees of freedom and must not be charged."""
    Z = _synthesize('R(10)-(R(200)|Q(2e-5,0.85))', noise=0.0)
    circuit = R("10") - (R(200) | Q(2e-5, 0.85))  # string value = fixed

    if optimizer == 'single':
        result = fit_equivalent_circuit(FREQ, Z, circuit, plot=False)[0]
    elif optimizer == 'de':
        result = fit_circuit_diffevo(circuit, FREQ, Z, maxiter=30, seed=0)[0].best_result
    else:
        result = fit_circuit_multistart(circuit, FREQ, Z, n_restarts=2)[0].best_result
    plt.close('all')

    assert result.n_free_params == 3


# =============================================================================
# Ranking mechanics
# =============================================================================

def test_failed_candidate_sorts_last_and_does_not_break_the_table():
    """A candidate that failed to fit is reported, not fatal."""
    Z = _synthesize('R(10)-(R(200)|C(2e-5))')
    good = _fit('R(10)-(R(200)|C(2e-5))', Z)

    scores = score_candidates(
        FREQ, Z, [('failed', None), ('good', good)], 'modulus'
    )

    assert scores[0].expression == 'good'
    assert scores[-1].result is None
    assert np.isinf(scores[-1].bic)

    # Must not raise on the None row
    log_comparison(scores, len(FREQ), 'modulus')


def test_non_finite_score_is_reported_as_failed():
    """A NaN score must not win: NaN compares false against everything, so it
    would neither sort last nor be rejected, and would poison the table's
    reference values."""
    Z = _synthesize('R(10)-(R(200)|C(2e-5))')
    good = _fit('R(10)-(R(200)|C(2e-5))', Z)

    broken = _fit('R(10)-(R(200)|C(2e-5))', Z)
    broken.params_opt = np.full(len(broken.params_opt), np.nan)

    scores = score_candidates(
        FREQ, Z, [('broken', broken), ('good', good)], 'modulus'
    )

    assert scores[0].expression == 'good'
    assert scores[-1].result is None
    log_comparison(scores, len(FREQ), 'modulus')


def test_unpopulated_n_free_params_is_loud():
    """k = 0 would charge no complexity penalty and win every comparison.

    The Voigt-chain FitResult leaves the field at its default, so this must
    raise rather than silently return a criterion of zero complexity.
    """
    Z = np.ones(80, dtype=complex) * (10 + 5j)
    with pytest.raises(ValueError, match='n_free_params must be positive'):
        compute_information_criteria(Z, Z * 1.01, 'modulus', 0)


def test_score_keeps_the_command_line_position():
    """Figures are saved by input position while the table sorts by BIC, so
    the position has to survive the sort."""
    Z = _synthesize('R(10)-(R(200)|Q(2e-5,0.85))')
    worse = _fit('R(10)-(R(200)|C(2e-5))', Z)
    better = _fit('R(10)-(R(200)|Q(2e-5,0.85))', Z)

    scores = score_candidates(
        FREQ, Z, [('worse', worse), ('better', better)], 'modulus'
    )

    assert scores[0].expression == 'better'
    assert scores[0].index == 2      # second on the command line
    assert scores[1].index == 1


def test_all_candidates_failed_produces_no_table():
    """Nothing to rank means nothing is logged, and no exception."""
    scores = score_candidates(
        FREQ, np.ones(80, dtype=complex), [('a', None), ('b', None)], 'modulus'
    )
    log_comparison(scores, len(FREQ), 'modulus')

    assert all(s.result is None for s in scores)
