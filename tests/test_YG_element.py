#!/usr/bin/env python3
"""Test YG element (Young-Göhr passive layer with an exponential conductivity profile)."""

import numpy as np
import pytest

from eis_analysis.fitting import C, R, YG
from eis_analysis.fitting.circuit_elements.composite import YG_P_MIN


@pytest.fixture
def freq():
    """Test frequencies: 1 MHz to 1 mHz, the oxide measurement window."""
    return np.logspace(6, -3, 60)


@pytest.fixture
def yg_params():
    """Zahner's own simulated example (manual 11/2023, section 2.3.9).

    C = 10 uF, tau = 100 ms, p = 0.05 - a film capacitance and a strong
    conductivity gradient, the regime the element was built for.
    """
    return 1e-5, 0.05, 0.1


def _naive_impedance(freq, C_val, p_val, tau_val):
    """The formula as printed, built straight from exp(1/p).

    An independent oracle for the overflow-free implementation, valid only
    while exp(1/p) fits in a float64 - hence p >= 0.05 at every call site.
    """
    omega = 2 * np.pi * freq
    a = 1j * omega * tau_val
    E = np.exp(1.0 / p_val)
    return p_val / (1j * omega * C_val) * np.log((1 + a * E) / (1 + a))


# --- 1. The impedance itself ---

def test_yg_matches_the_printed_formula(freq, yg_params):
    """The log-space rewrite must reproduce the closed form where both work.

    The rewrite exists only to survive small p; it must not change the value
    anywhere the naive form is valid. Measured agreement is at machine
    precision, so 1e-14 relative is a real constraint, not a loose one.
    """
    yg = YG(*yg_params)
    Z = yg.impedance(freq, list(yg_params))
    Z_naive = _naive_impedance(freq, *yg_params)

    assert np.max(np.abs(Z - Z_naive) / np.abs(Z_naive)) < 1e-14


def test_yg_stays_finite_where_the_naive_formula_overflows(freq):
    """p = 1e-3 is inside the bounds and overflows exp(1/p); Z must survive.

    This is the whole reason `_yg_log_terms` exists. The naive form returns
    inf/nan here, so a regression would be silent in a fit: least_squares
    sees nan residuals and stops.
    """
    yg = YG(1e-5, 1e-3, 0.1)
    Z = yg.impedance(freq, [1e-5, 1e-3, 0.1])
    assert np.all(np.isfinite(Z))

    with np.errstate(over='ignore', invalid='ignore'):
        Z_naive = _naive_impedance(freq, 1e-5, 1e-3, 0.1)
    assert not np.all(np.isfinite(Z_naive)), "oracle no longer overflows"


def test_yg_high_frequency_limit_is_the_capacitance(yg_params):
    """Above the capacitive corner the element is C - the fitted parameter.

    This is what lets the oxide analysis read the thickness straight off C
    with no Hsu-Mansfeld/Brug conversion in between.
    """
    C_val, p_val, tau_val = yg_params
    f_hi = np.array([1e9])          # far above 1/(2*pi*tau) = 1.59 Hz
    Z = YG(*yg_params).impedance(f_hi, list(yg_params))
    Z_ideal = C(C_val).impedance(f_hi, [C_val])

    assert np.abs(Z[0] - Z_ideal[0]) / np.abs(Z_ideal[0]) < 1e-8


def test_yg_low_frequency_limit_is_R_dc(yg_params):
    """Below the resistive corner the element is a real resistance p*tau*(E-1)/C."""
    yg = YG(*yg_params)
    f_lo = np.array([yg.dc_corner_freq * 1e-4])
    Z = yg.impedance(f_lo, list(yg_params))

    assert np.abs(Z[0].real - yg.R_dc) / yg.R_dc < 1e-3
    assert abs(Z[0].imag) < 1e-3 * yg.R_dc


def test_yg_degenerates_to_an_ideal_capacitor_as_p_vanishes(freq, yg_params):
    """p -> 0 is a plain C: zero penetration depth, no conductivity at all.

    Checked on both sides of YG_P_MIN so the explicit guard and the general
    branch agree - the guard must be the limit, not a different element.
    """
    C_val = yg_params[0]
    Z_ideal = C(C_val).impedance(freq, [C_val])

    Z_above = YG(C_val, 1e-3, 0.1).impedance(freq, [C_val, 1e-3, 0.1])
    assert np.max(np.abs(Z_above - Z_ideal) / np.abs(Z_ideal)) < 0.02

    p_below = YG_P_MIN / 2
    Z_below = YG(C_val, p_below, 0.1).impedance(freq, [C_val, p_below, 0.1])
    assert np.allclose(Z_below, Z_ideal, rtol=1e-15)


def test_yg_phase_follows_zahners_cpe_approximation(yg_params):
    """The mid-band phase must track phi = -90*(1 - q), q = 1/(ln(w*tau) + 1/p).

    Zahner's approximation, and only an approximation: measured deviation
    reaches 6 degrees, so the tolerance is loose on purpose. It still pins
    the sign, the CPE-like plateau and the parametrisation of p - a swapped
    p and tau, or a dropped 1/p, moves the phase by tens of degrees.
    """
    C_val, p_val, tau_val = yg_params
    f_mid = np.logspace(2, -1, 10)
    omega = 2 * np.pi * f_mid

    phase = np.angle(YG(*yg_params).impedance(f_mid, list(yg_params)), deg=True)
    q = 1.0 / (np.log(omega * tau_val) + 1.0 / p_val)
    phase_approx = -90.0 * (1.0 - q)

    assert np.max(np.abs(phase - phase_approx)) < 8.0
    assert np.all(phase < 0) and np.all(phase > -90.0)


# --- 2. Derived properties ---

def test_yg_corner_frequencies_bracket_the_cpe_band(yg_params):
    """The capacitive corner is 1/(2*pi*tau); the resistive one is e^(1/p) below it."""
    yg = YG(*yg_params)
    assert yg.characteristic_freq == pytest.approx(1.0 / (2 * np.pi * 0.1))
    assert yg.dc_corner_freq == pytest.approx(
        yg.characteristic_freq * np.exp(-1.0 / 0.05))
    assert yg.dc_corner_freq < yg.characteristic_freq


def test_yg_R_dc_is_infinite_at_the_degenerate_p(yg_params):
    """p below YG_P_MIN is an ideal capacitor: no DC path, R = inf, not nan.

    Reachable despite the bounds - YG(1e-5, "5e-4", 0.1) fixes p by string
    and skips them.
    """
    assert YG(1e-5, YG_P_MIN / 2, 0.1).R_dc == float('inf')
    assert np.isfinite(YG(*yg_params).R_dc)


# --- 3. Construction, repr, fixed parameters ---

def test_yg_repr_and_fixed_params():
    """String arguments fix a parameter and are echoed quoted in the repr."""
    yg = YG("1e-5", 0.05, "0.1")
    assert yg.fixed_params == [True, False, True]
    assert yg.get_param_labels() == ['C_YG', 'p_YG', 'τ_YG']
    assert 'YG(' in repr(yg) and '"1e-05"' in repr(yg)


def test_yg_composes_with_the_operators(freq, yg_params):
    """A YG in series must add its impedance to the rest of the circuit."""
    circuit = R(20) - YG(*yg_params)
    Z = circuit.impedance(freq, [20.0, *yg_params])
    Z_yg = YG(*yg_params).impedance(freq, list(yg_params))

    assert np.allclose(Z - Z_yg, 20.0)
