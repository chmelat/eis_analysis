#!/usr/bin/env python3
"""Test YG element (Young-Göhr passive layer with an exponential conductivity profile)."""

import numpy as np
import pytest

from eis_analysis.fitting import C, R, YG, fit_equivalent_circuit
from eis_analysis.fitting.bounds import (
    classify_bound_status, generate_simple_bounds, log_scale_ci_mask)
from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting.circuit_elements.composite import YG_P_MIN
from eis_analysis.fitting.jacobian import element_jacobian


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


# --- 4. Analytic Jacobian ---

@pytest.mark.parametrize("p_val", [0.01, 0.05, 0.1, 0.5])
def test_yg_jacobian_matches_central_differences(freq, p_val):
    """All three columns must match finite differences across the p range.

    The metric is max|analytic - numeric| / max|numeric| over the whole
    sweep, not a per-point ratio: where a derivative passes through zero the
    central difference is pure cancellation and a per-point ratio reports
    1e-2 for a derivative that is in fact exact to 1e-10.
    """
    C_val, tau_val = 1e-5, 0.1
    params = [C_val, p_val, tau_val]
    yg = YG(*params)

    Z, dZ = element_jacobian(yg, freq, params)
    assert np.max(np.abs(Z - yg.impedance(freq, params))) < 1e-15

    for col in range(3):
        step = 1e-5 * params[col]
        up, down = list(params), list(params)
        up[col] += step
        down[col] -= step
        numeric = (yg.impedance(freq, up) - yg.impedance(freq, down)) / (2 * step)
        error = np.max(np.abs(dZ[:, col] - numeric)) / np.max(np.abs(numeric))
        assert error < 1e-6, f"column {col} (p = {p_val}): {error:.2e}"


def test_yg_jacobian_is_finite_at_the_degenerate_p(freq):
    """Below YG_P_MIN the element is a capacitor: dZ/dp and dZ/dtau are zero.

    Not nan. least_squares treats a nan column as a failure of the whole
    fit, so the guard has to produce the limit here too, exactly as
    `impedance` does.
    """
    params = [1e-5, YG_P_MIN / 2, 0.1]
    Z, dZ = element_jacobian(YG(*params), freq, params)

    assert np.all(np.isfinite(dZ))
    assert np.allclose(dZ[:, 0], -Z / 1e-5)
    assert np.all(dZ[:, 1] == 0) and np.all(dZ[:, 2] == 0)


# --- 5. Bounds and CLI parsing ---

def test_yg_bounds_are_registered(yg_params):
    """Every label needs its own entry; a missing one falls back silently.

    generate_simple_bounds looks bounds up by the bare label string and
    returns DEFAULT_BOUNDS = (1e-15, 1e15) for anything it does not know -
    30 decades, no error. The upper-bound assertions below are what catches
    that.
    """
    lower, upper = generate_simple_bounds(YG(*yg_params).get_param_labels())

    assert upper[0] < 1e15 and upper[1] < 1e15 and upper[2] < 1e15
    # C and tau are scale parameters, p is a bounded ratio like n and alpha_CC
    assert log_scale_ci_mask(lower, upper) == [True, False, True]


def test_yg_defaults_lie_inside_their_bounds():
    """A default outside its own bounds gets clipped before the fit starts."""
    yg = YG()
    lower, upper = generate_simple_bounds(yg.get_param_labels())

    for value, lo, hi, label in zip(yg.get_all_params(), lower, upper,
                                    yg.get_param_labels()):
        assert lo <= value <= hi, f"{label} = {value:g} outside ({lo:g}, {hi:g})"


def test_yg_p_bounds_do_not_flag_a_well_determined_p():
    """p = 0.01 is Zahner's own example and must not read as "at its bound".

    classify_bound_status uses 1% of the range on a linear parameter, so the
    upper bound of 0.5 is what keeps that threshold below 0.01. Raising it to
    1.0 would make every strong conductivity gradient look constrained.
    """
    lower, upper = generate_simple_bounds(['p_YG'])

    assert classify_bound_status(0.01, lower[0], upper[0]) == ''
    assert classify_bound_status(0.05, lower[0], upper[0]) == ''
    # the genuinely degenerate ends still report
    assert classify_bound_status(1.5e-3, lower[0], upper[0]) == 'lower'
    assert classify_bound_status(0.499, lower[0], upper[0]) == 'upper'


def test_yg_parses_from_circuit_string():
    """The CLI reaches elements only through parse_circuit_expression."""
    circuit = parse_circuit_expression("L(1e-6) - R(20) - YG(1e-5, 0.05, 0.1)")

    assert circuit.get_param_labels()[-3:] == ['C_YG', 'p_YG', 'τ_YG']
    assert circuit.get_all_params()[-3:] == [1e-5, 0.05, 0.1]


# --- 6. Round-trip fit ---

def test_yg_round_trip_fit_recovers_the_layer(freq, yg_params):
    """A noisy synthetic oxide spectrum must return C and p from a bad guess.

    Also the end-to-end proof that the analytic Jacobian is wired: a missing
    branch raises RuntimeError out of least_squares rather than falling back,
    and params_significance would come back None.
    """
    C_val, p_val, tau_val = yg_params
    truth = R(20) - YG(*yg_params)
    Z_true = truth.impedance(freq, [20.0, *yg_params])

    rng = np.random.default_rng(42)
    Z = Z_true * (1 + 0.01 * rng.standard_normal(len(freq)))

    guess = R(50) - YG(3e-5, 0.15, 0.03)      # deliberately off
    result, _, _ = fit_equivalent_circuit(freq, Z, guess, weighting='modulus',
                                         plot=False)

    assert result.fit_error_rel < 2.0
    # None when any element in the circuit lacks an analytic derivative
    assert result.params_significance is not None

    _, C_fit, p_fit, tau_fit = result.params_opt
    assert C_fit == pytest.approx(C_val, rel=0.05)
    assert p_fit == pytest.approx(p_val, rel=0.10)
    assert tau_fit == pytest.approx(tau_val, rel=0.20)
