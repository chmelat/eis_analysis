"""Test the DE archive check (fitting/de_archive.py and its use in diffevo).

1. Candidate selection takes the best distinct point of each early window,
   so a basin visited early survives a late collapse of the population.
2. Distinctness is a decade in a scale parameter or 0.15 in a CPE exponent.
3. Ambiguity needs both a statistically equal cost and a different spectrum.
4. End to end: a DE run trapped in a local minimum is repaired and says so;
   a correct fit is left alone; with workers > 1 the check is skipped.
"""

import numpy as np
import pytest

from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting import fit_circuit_diffevo
from eis_analysis.fitting.de_archive import (
    ARCHIVE_WINDOWS, Refinement, ambiguous_alternatives, choose, select_archive_candidates)


def _archive(n_pop=10, n_gen=100):
    """Basin B (R ~ 1e3) explored in the first 20 generations, then everything in A (R ~ 1)."""
    rng = np.random.default_rng(0)
    n = n_pop * n_gen
    params = np.column_stack([10 ** rng.normal(0, 0.05, n), np.full(n, 0.8)])
    costs = 1.0 + rng.random(n)
    early = np.arange(n) < 20 * n_pop
    params[early, 0] = 10 ** rng.normal(3, 0.05, early.sum())
    costs[early] = 5.0 + rng.random(early.sum())          # B worse than A as DE sees it
    return costs, params


def test_early_basin_survives_the_late_collapse():
    costs, params = _archive()
    picked = select_archive_candidates(costs, params, 10, np.array([False, True]))
    assert picked[0] == int(np.argmin(costs))
    assert any(params[i, 0] > 100 for i in picked[1:])   # a point from basin B
    assert len(picked) <= ARCHIVE_WINDOWS + 1


def test_picked_points_are_pairwise_distinct():
    costs, params = _archive()
    picked = select_archive_candidates(costs, params, 10, np.array([False, True]))
    logs = np.log10(params[picked, 0])
    assert all(abs(a - b) >= 1 for i, a in enumerate(logs) for b in logs[i + 1:])


@pytest.mark.parametrize("dn, distinct", [(0.1, False), (0.2, True)])
def test_exponent_distance_is_scaled(dn, distinct):
    params = np.array([[1.0, 0.7], [1.0, 0.7 + dn]])
    costs = np.array([1.0, 2.0])
    picked = select_archive_candidates(costs, params, 1, np.array([False, True]))
    assert (len(picked) == 2) == distinct


def test_non_finite_costs_are_never_picked():
    """An overflowed impedance gives a NaN cost; argmin would return it as 'DE's best'."""
    costs, params = _archive()
    costs[5] = np.nan
    params[5, 0] = 1e9                                    # far from everything: would be picked
    picked = select_archive_candidates(costs, params, 10, np.array([False, True]))
    assert 5 not in picked and np.isfinite(costs[picked]).all()


def test_empty_archive_gives_no_candidates():
    assert select_archive_candidates(np.array([]), np.zeros((0, 2)), 10, np.array([False, False])) == []


# 50 points, 2 free parameters: s^2 = best / 98 = 0.01, unit weights
Z50 = np.full(50, 10 - 5j)
W50 = np.ones(50)


def _ref(cost, Z):
    return Refinement(cost, np.asarray(Z, dtype=complex), [cost], None)


def test_ambiguity_needs_equal_cost_and_a_distinguishable_spectrum():
    best = 0.98
    far = Z50 + 0.2                                       # separation 50*0.04/0.01 = 200
    refined = [_ref(best + 0.05, far),                    # delta-chi2 = 5: ambiguous
               _ref(best + 0.50, far),                    # delta-chi2 = 50: clearly worse
               _ref(best + 0.01, Z50)]                    # same spectrum (permutation)
    alts = ambiguous_alternatives(_ref(best, Z50), refined, W50, 2)
    assert len(alts) == 1 and alts[0][0] == pytest.approx(5.0)


def test_spectrum_inside_the_residual_level_is_not_another_model():
    """A poor fit (large s^2) reached at a slightly different valley point."""
    best = 98.0                                           # s^2 = 1: |Z| ~ 11, residual ~ 1
    nearby = Z50 * 1.02                                   # 2 % off: separation 50*0.05/1 = 2.5
    assert ambiguous_alternatives(_ref(best, Z50), [_ref(best + 0.5, nearby)], W50, 2) == []


def test_rounding_level_differences_are_not_another_model():
    """A noise-free fit: candidates refined into the same minimum differ at ~1e-15."""
    nearby = Z50 * (1 + 1e-14)
    assert ambiguous_alternatives(_ref(2e-29, Z50), [_ref(1e-29, nearby)], W50, 2) == []


def test_exact_fit_still_sees_a_genuinely_equal_model():
    """Two models fitting noise-free data exactly but predicting differently elsewhere."""
    other = Z50 * 1.05
    assert len(ambiguous_alternatives(_ref(0.0, Z50), [_ref(0.0, other)], W50, 2)) == 1


def test_archive_candidate_in_the_same_minimum_does_not_churn():
    from_de = _ref(1.0, Z50)
    sel = choose(5.0, from_de, [_ref(1.0 - 1e-9, Z50)], W50, 2)
    assert sel.chosen is from_de and not sel.local_minimum


def test_better_archive_candidate_wins_and_is_reported():
    from_de, cand = _ref(1.0, Z50), _ref(0.1, Z50 * 1.05)
    sel = choose(5.0, from_de, [cand], W50, 2)
    assert sel.chosen is cand and sel.best_refined is cand and sel.local_minimum


def test_de_point_kept_when_every_refinement_is_worse():
    sel = choose(0.5, _ref(1.0, Z50), [_ref(0.8, Z50 * 1.05)], W50, 2)
    assert sel.chosen is None and sel.best_refined.cost == 0.8


def test_archive_alone_when_the_refinement_from_de_failed():
    cand = _ref(0.1, Z50)
    sel = choose(5.0, None, [cand], W50, 2)
    assert sel.chosen is cand and not sel.local_minimum


FREQ = np.logspace(5, -3, 129)


def test_trapped_run_is_repaired_and_reported(monkeypatch):
    """randtobest1bin with seed 3 stops after ~49 generations in the W-like basin of Wo.

    The bounds are pinned: the trajectory of a seeded DE run depends on them,
    and a later change of PARAMETER_BOUNDS must not decide whether DE is trapped.
    """
    from eis_analysis.fitting import bounds
    for label, rng in {'R': (1e-4, 1e10), 'R_W': (1e-4, 1e10), 'τ_W': (1e-6, 1e4)}.items():
        monkeypatch.setitem(bounds.PARAMETER_BOUNDS, label, rng)
    truth = parse_circuit_expression("R(10)-(R(100)|Wo(100,1))")
    Z = truth.impedance(FREQ, truth.get_all_params())
    start = parse_circuit_expression("R(5)-(R(50)|Wo(50,0.3))")
    result, _ = fit_circuit_diffevo(start, FREQ, Z, seed=3, strategy=1)
    assert result.diagnostics.de_error > 1.0              # DE itself was trapped
    assert result.diagnostics.archive_used
    assert result.best_result.params_opt == pytest.approx([10, 100, 100, 1], rel=1e-6)
    warnings = result.diagnostics.warnings
    assert any("local minimum" in w for w in warnings)
    # Diagnostics describe the fit that is returned, not the discarded one
    assert result.diagnostics.refined_error == pytest.approx(result.final_error)
    assert result.improvement > 99.9
    assert not any("using DE result" in w or "contributed nothing" in w for w in warnings)


def test_correct_fit_is_left_alone():
    rng = np.random.default_rng(3)
    truth = parse_circuit_expression("R(10)-(R(200)|C(2e-5))")
    Z0 = truth.impedance(FREQ, truth.get_all_params())
    Z = Z0 + 0.01 * np.abs(Z0) * (rng.standard_normal(len(FREQ)) + 1j * rng.standard_normal(len(FREQ)))
    result, _ = fit_circuit_diffevo(parse_circuit_expression("R(1)-(R(1)|C(1e-6))"), FREQ, Z, seed=0)
    diag = result.diagnostics
    assert diag.archive_candidates > 0 and not diag.archive_used
    assert not any("Ambiguous" in w or "local minimum" in w for w in diag.warnings)


@pytest.mark.parametrize("truth, expect_warning", [
    ("R(10)-(R(200)|C(2e-5))", True),           # the model fits: ambiguity is reported
    ("R(10)-(R(200)|Q(2e-5,0.6))", False),      # an R|C cannot fit n = 0.6: Poor, not reported
])
def test_ambiguity_reported_only_for_an_acceptable_fit(monkeypatch, truth, expect_warning):
    """An alternative is injected; whether it is shown depends on the fit quality alone."""
    from eis_analysis.fitting import diffevo
    real_choose = diffevo.choose

    def choose_with_alternative(*args, **kwargs):
        sel = real_choose(*args, **kwargs)
        sel.alternatives = [(1.0, np.array([1.0, 2.0, 3.0]))]
        return sel

    monkeypatch.setattr(diffevo, "choose", choose_with_alternative)
    tc = parse_circuit_expression(truth)
    Z = tc.impedance(FREQ, tc.get_all_params())
    result, _ = fit_circuit_diffevo(parse_circuit_expression("R(1)-(R(1)|C(1e-6))"), FREQ, Z, seed=0)
    assert (result.final_error < 5.0) == expect_warning
    assert any("Ambiguous" in w for w in result.diagnostics.warnings) == expect_warning


def test_archive_check_runs_with_workers():
    """Worker processes keep their own archives; the parent records through the map."""
    truth = parse_circuit_expression("R(10)-(R(100)|Wo(100,1))")
    Z = truth.impedance(FREQ, truth.get_all_params())
    start = parse_circuit_expression("R(5)-(R(50)|Wo(50,0.3))")
    result, _ = fit_circuit_diffevo(start, FREQ, Z, seed=1, strategy=1, workers=2)
    assert result.diagnostics.archive_checked and result.diagnostics.archive_candidates > 0
    assert result.final_error < 0.01


def test_archive_check_can_be_switched_off():
    truth = parse_circuit_expression("R(10)-(R(200)|C(2e-5))")
    Z = truth.impedance(FREQ, truth.get_all_params())
    result, _ = fit_circuit_diffevo(parse_circuit_expression("R(1)-(R(1)|C(1e-6))"),
                                    FREQ, Z, seed=0, archive_check=False)
    diag = result.diagnostics
    assert not diag.archive_checked and diag.archive_candidates == 0


def test_same_alternative_reached_twice_counts_once():
    best = 0.98
    far = Z50 + 0.2
    refined = [_ref(best + 0.05, far), _ref(best + 0.06, far * (1 + 1e-6))]
    assert len(ambiguous_alternatives(_ref(best, Z50), refined, W50, 2)) == 1
