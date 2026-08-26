"""Tests for fixed parameters (values passed as strings: R("15")).

Regression tests (2026-08-26): fixing parameters must behave identically in
every optimizer.
- A circuit with *all* parameters fixed has an empty optimization vector;
  scipy then failed with an opaque message ("bounds should be a sequence
  containing finite real valued (min, max) pairs" from DE, "index -1 is out
  of bounds" from least_squares) instead of saying what was wrong.
- The local-fit path clipped fixed values into the free-parameter bounds, so
  R("1e-9") was silently fitted as R = 1e-4 (DE kept 1e-9).
- A fixed value outside those bounds (R("-5"), CPE n = 1.5) passed silently;
  fixing bypasses the bounds, so it now warns.
"""

import numpy as np
import pytest

from eis_analysis.fitting import (
    fit_circuit_diffevo,
    fit_circuit_multistart,
    fit_equivalent_circuit,
)
from eis_analysis.fitting.circuit_elements import Q, R


def _data(n_freq=30):
    """Noise-free R_s - (R_ct || Q) spectrum."""
    freq = np.logspace(-1, 5, n_freq)
    omega = 2 * np.pi * freq
    Z = 15.0 + 1 / (1 / 5e4 + 1e-6 * (1j * omega) ** 0.88)
    return freq, Z


def _fit(optimizer, circuit, freq, Z):
    """Run one optimizer, return its FitResult."""
    if optimizer == 'single':
        return fit_equivalent_circuit(freq, Z, circuit, plot=False)[0]
    if optimizer == 'de':
        return fit_circuit_diffevo(circuit, freq, Z, maxiter=30, seed=0)[0].best_result
    return fit_circuit_multistart(circuit, freq, Z, n_restarts=2)[0].best_result


OPTIMIZERS = ['single', 'de', 'multistart']


@pytest.mark.parametrize('optimizer', OPTIMIZERS)
def test_all_params_fixed_reports_nothing_to_fit(optimizer):
    """An all-fixed circuit must fail with an explanation, not a scipy message."""
    freq, Z = _data()
    circuit = R("15") - (R("5e4") | Q("1e-6", "0.88"))

    with pytest.raises((ValueError, RuntimeError), match="nothing to fit"):
        _fit(optimizer, circuit, freq, Z)


@pytest.mark.parametrize('optimizer', OPTIMIZERS)
@pytest.mark.parametrize('fixed_value', [1e-9, 1e12])  # below / above R bounds
def test_fixed_value_outside_bounds_is_kept(optimizer, fixed_value):
    """A fixed value must reach the fit unchanged, bounds notwithstanding."""
    freq, Z = _data()
    circuit = R(str(fixed_value)) - (R(1e4) | Q(1e-6, 0.9))

    result = _fit(optimizer, circuit, freq, Z)

    assert result.params_opt[0] == fixed_value
    assert result.bound_status[0] == 'fixed'


@pytest.mark.parametrize('optimizer', OPTIMIZERS)
def test_fixed_value_outside_bounds_warns(optimizer):
    """R = -5 Ohm is nonphysical; fixing bypasses the bounds, so warn."""
    freq, Z = _data()
    circuit = R("-5") - (R(1e4) | Q(1e-6, 0.9))

    result = _fit(optimizer, circuit, freq, Z)

    assert any('Fixed parameter R0' in w and 'outside the range' in w
               for w in result.all_warnings)


@pytest.mark.parametrize('optimizer', OPTIMIZERS)
def test_fixed_value_inside_bounds_is_silent(optimizer):
    """A sane fixed value must not produce a warning - and must be honored."""
    freq, Z = _data()
    circuit = R("15") - (R(1e4) | Q(1e-6, 0.9))

    result = _fit(optimizer, circuit, freq, Z)

    assert result.params_opt[0] == 15.0
    assert not any('Fixed parameter' in w for w in result.all_warnings)
    assert result.params_stderr[0] == 0.0  # fixed -> known exactly
