"""Tests for the series resistance determined only together with an arc.

An (R|Q) arc whose time constant lies far above the measured frequencies is
a constant resistance: the data fix only Rs + R_k. A fit with that arc as a
small low-frequency feature gave Rs a CI of +/- 1 % around the sum, missing
the true Rs (AutoECM benchmark, 22 of 47 covered). series_resistance.py
widens it to the flat profile.
"""

from unittest import mock

import numpy as np

import eis_analysis.fitting.circuit as circuit_module
from eis_analysis.fitting import fit_circuit_multistart, fit_equivalent_circuit
from eis_analysis.fitting.circuit_elements import C, Q, R
from eis_analysis.fitting.series_resistance import series_resistance_split

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


def test_degenerate_arc_fits_as_well_and_a_nan_candidate_is_skipped():
    """The constant arc fits as well with dchi2 at rounding level of either
    sign, which the warning read as "not the least-squares minimum" (F5);
    a first candidate with a NaN dchi2 blocked the later one (F3)."""
    f = np.logspace(5, -1, 61)
    true = R(10) - (R(50) | C(1e-5)) - (R(5) | C(1e-12))
    rng = np.random.default_rng(4)
    noise = (rng.normal(size=f.size) + 1j * rng.normal(size=f.size)) / np.sqrt(2)
    Z = true.impedance(f, true.get_all_params()) * (1 + 0.01 * noise)
    result, _ = fit_equivalent_circuit(f, Z, R(10) - (R(50) | C(1e-5)) - (R(5) | C(1e-9)))
    split = result.series_split
    assert split is not None and split.arc_resistance == 'R2' and abs(split.dchi2) < 1e-3
    warning, = (w for w in result.all_warnings if "is determined only as" in w)
    assert "fits as well" in warning and "not the least-squares minimum" not in warning

    fit_local = circuit_module._fit_local

    def nan_for_r1(*args, **kwargs):
        merged, Z_fit = fit_local(*args, **kwargs)
        # The merge of R1 starts from [R0 + R1, R2, C2]
        if kwargs['initial_guess'][1] == result.params_opt[3]:
            merged.params_opt = merged.params_opt * np.nan
        return merged, Z_fit

    with mock.patch.object(circuit_module, '_fit_local', nan_for_r1):
        split = series_resistance_split(result, f, Z)
    assert split is not None and split.arc_resistance == 'R2'
