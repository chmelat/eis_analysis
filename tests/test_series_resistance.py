"""Tests for the series resistance determined only together with an arc.

An (R|Q) arc whose time constant lies far above the measured frequencies is
a constant resistance: the data fix only Rs + R_k. A fit with that arc as a
small low-frequency feature gave Rs a CI of +/- 1 % around the sum, missing
the true Rs (AutoECM benchmark, 22 of 47 covered). series_resistance.py
widens it to the flat profile.
"""

import numpy as np

from eis_analysis.fitting import fit_circuit_multistart, fit_equivalent_circuit
from eis_analysis.fitting.circuit_elements import Q, R

FREQ = np.logspace(-1, 5, 40)


def _noisy(circuit, seed):
    Z = circuit.impedance(FREQ, circuit.get_all_params())
    rng = np.random.default_rng(seed)
    return Z + 0.01 * np.abs(Z) * (rng.standard_normal(Z.size) + 1j * rng.standard_normal(Z.size))


# tau of the first arc ~1e-15 s: 200 Ohm in series with Rs = 20 Ohm
CONSTANT_ARC = R(20) - (R(200) | Q(1e-11, 0.6)) - (R(1000) | Q(1e-5, 0.85))


def _lf_start():
    return R(50) - (R(30) | Q(1e-3, 0.9)) - (R(1000) | Q(1e-5, 0.85))


def test_constant_arc_widens_rs_ci():
    # Seed 4: the fit puts Rs at 219.8 +/- 0.7, the arc as an LF feature
    result, _ = fit_equivalent_circuit(FREQ, _noisy(CONSTANT_ARC, 4), _lf_start())
    split = result.series_split
    assert split is not None and split.index == 0 and split.dchi2 < 3.84
    # The arc's R (true 200 Ohm, fitted as the small LF feature) takes the rest
    i_arc = split.arc_index
    assert i_arc is not None
    for low, high in (result.params_ci_95, result.params_ci_99):
        assert low[0] == split.lower
        assert low[0] < 20 < high[0]
        assert high[0] > 215
        assert low[i_arc] == split.arc_lower and high[i_arc] > 200
    assert any(w.startswith("R0 is determined only as R0 + ") for w in result.all_warnings)


def test_multistart_checks_the_best_start():
    result, _ = fit_circuit_multistart(FREQ, _noisy(CONSTANT_ARC, 4), _lf_start(),
                                       n_restarts=3, rng=0)
    assert result.best_result.series_split is not None
    assert result.best_result.params_ci_95[0][0] < 20


def test_arcs_in_window_keep_rs_ci():
    true = R(20) - (R(300) | Q(1e-7, 0.8)) - (R(1000) | Q(1e-5, 0.85))
    start = R(30) - (R(200) | Q(2e-7, 0.85)) - (R(800) | Q(2e-5, 0.8))
    result, _ = fit_equivalent_circuit(FREQ, _noisy(true, 0), start)
    assert result.series_split is None
    assert not any("is determined only as" in w for w in result.all_warnings)


def test_no_series_resistance():
    true = (R(300) | Q(1e-7, 0.8)) - (R(1000) | Q(1e-5, 0.85))
    result, _ = fit_equivalent_circuit(FREQ, _noisy(true, 0), true)
    assert result.series_split is None
