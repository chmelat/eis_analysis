#!/usr/bin/env python3
"""Test DQ element (bounded power-law DRT, truncated CPE)."""

import numpy as np
import pytest
from scipy.integrate import quad

from eis_analysis.fitting import K, Q, R, DQ, fit_equivalent_circuit
from eis_analysis.fitting.jacobian import element_jacobian
from eis_analysis.fitting.bounds import generate_simple_bounds, log_scale_ci_mask
from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.analysis.oxide import analyze_oxide_layer
from eis_analysis.fitting.circuit import FitResult


@pytest.fixture
def freq():
    """Test frequencies: 1 MHz to 1 mHz, the oxide measurement window."""
    return np.logspace(6, -3, 60)


@pytest.fixture
def dq_params():
    """Oxide-like parameters: n = 0.57, distribution inside the window.

    A is Ohm*s^-n, not a capacitance: A = (1-n)*tau_min^(1-n)/C_eff puts
    C_eff at ~1e-7 F and R_pol at ~4e7 Ohm, the scale of a ZrO2 film.
    """
    return 1.2e6, 0.57, 5e-2, 8.0


def test_dq_matches_adaptive_quadrature(dq_params):
    """Gauss-Legendre must match QUADPACK, an independent integrator.

    The point of the element is to explain ~1e-3 structure in the residuals,
    so a quadrature error of that size would be indistinguishable from the
    effect being measured. This pins the integration itself.
    """
    A, n, tau_min, U = dq_params
    dq = DQ(A, n, tau_min, U)
    s_min, s_max = np.log(tau_min), np.log(tau_min) + U

    for f in (1e-3, 1e-1, 3.07, 1e2, 1e5):
        omega = 2 * np.pi * f

        def integrand(s, part):
            z = np.exp(n * s) / (1 + 1j * omega * np.exp(s))
            return z.real if part == 're' else z.imag

        ref = A * complex(quad(integrand, s_min, s_max, args=('re',), limit=200)[0],
                          quad(integrand, s_min, s_max, args=('im',), limit=200)[0])
        Z = dq.impedance(np.array([f]), list(dq_params))[0]
        assert abs(Z - ref) / abs(ref) < 1e-9, f"quadrature off at {f} Hz"


def test_dq_wide_limit_is_ideal_cpe(freq):
    """Bounds far outside the window: DQ -> CPE with Q = sin(pi*n)/(pi*A).

    This is the defining relation to the element DQ generalises; it fails if
    the amplitude convention or the d(ln tau) measure ever changes.

    The tolerance is a *truncation* floor, not a quadrature one: the CPE is
    the limit of infinite bounds, and the mass beyond tau_max still missing
    at 1 mHz is (w*tau_max)^(n-1)*sin(pi*n)/(pi*(1-n)) ~ 2e-4. Tightening it
    would need a wider distribution than the quadrature can resolve, and a
    wrong amplitude convention is off by tens of percent, not by 1e-3.
    """
    A, n = 1e-3, 0.6
    tau_min, U = 1e-15, 60.0  # ~26 decades, window sits well inside
    Z_dq = DQ(A, n, tau_min, U).impedance(freq, [A, n, tau_min, U])

    Q_cpe = np.sin(np.pi * n) / (np.pi * A)
    Z_cpe = Q(Q_cpe, n).impedance(freq, [Q_cpe, n])

    rel_err = np.abs(Z_dq - Z_cpe) / np.abs(Z_cpe)
    assert np.max(rel_err) < 1e-3, f"wide DQ is not an ideal CPE: {np.max(rel_err)}"
    assert np.median(rel_err) < 1e-4, "wide DQ drifts from the CPE mid-window"


def test_dq_narrow_limit_is_single_rc(freq):
    """U -> 0 collapses the distribution onto one Voigt element."""
    A, n, tau_min, U = 1e-3, 0.6, 1e-2, 1e-6
    Z_dq = DQ(A, n, tau_min, U).impedance(freq, [A, n, tau_min, U])

    # All the mass sits at tau_min with weight A*tau_min^n*U
    R = A * tau_min ** n * U
    Z_rc = K(R, tau_min).impedance(freq, [R, tau_min])

    rel_err = np.max(np.abs(Z_dq - Z_rc) / np.abs(Z_rc))
    assert rel_err < 1e-6, f"narrow DQ is not a single RC: {rel_err}"


def test_dq_three_regimes(dq_params):
    """DC -> R_pol, high frequency -> C_eff, and a CPE slope in between."""
    A, n, tau_min, U = dq_params
    dq = DQ(A, n, tau_min, U)
    params = list(dq_params)

    # Well below 1/tau_max: finite, real polarisation resistance
    Z_dc = dq.impedance(np.array([1e-8]), params)[0]
    assert Z_dc.real == pytest.approx(dq.R_pol, rel=1e-4)
    assert abs(Z_dc.imag) < 1e-4 * abs(Z_dc.real)

    # Well above 1/tau_min: capacitive
    f_hi = 1e8
    Z_hi = dq.impedance(np.array([f_hi]), params)[0]
    Z_cap = 1 / (1j * 2 * np.pi * f_hi * dq.C_eff)
    assert abs(Z_hi - Z_cap) / abs(Z_cap) < 1e-3

    # Inside the bounds: |Z| falls with slope -n (fit over the central decades)
    f_mid = np.logspace(np.log10(1 / (2 * np.pi * dq.tau_max)) + 1,
                        np.log10(1 / (2 * np.pi * tau_min)) - 1, 30)
    slope = np.polyfit(np.log10(f_mid),
                       np.log10(np.abs(dq.impedance(f_mid, params))), 1)[0]
    assert slope == pytest.approx(-n, abs=0.03)


def test_dq_derived_properties(dq_params):
    """tau_max, R_pol and C_eff against their closed forms."""
    A, n, tau_min, U = dq_params
    dq = DQ(A, n, tau_min, U)

    assert dq.tau_max == pytest.approx(tau_min * np.exp(U), rel=1e-12)
    assert dq.R_pol == pytest.approx(
        A * (dq.tau_max ** n - tau_min ** n) / n, rel=1e-12)
    assert dq.C_eff == pytest.approx(
        (1 - n) / (A * (tau_min ** (n - 1) - dq.tau_max ** (n - 1))), rel=1e-12)

    # n = 1 is reachable (it is the upper bound) and must not divide by zero
    assert DQ(A, 1.0, tau_min, U).C_eff == pytest.approx(1 / (A * U), rel=1e-12)


def test_dq_fixed_params_and_repr(dq_params):
    """Strings fix parameters, as for every other element."""
    A, n, tau_min, U = dq_params
    dq = DQ(A, str(n), tau_min, U)

    assert dq.fixed_params == [False, True, False, False]
    assert dq.n == pytest.approx(n)
    assert 'DQ(' in repr(dq)


def test_dq_jacobian_matches_central_differences(freq, dq_params):
    """Analytic Jacobian against central differences, all four columns.

    The (tau_min, U) parametrisation is the subtle one: tau_min shifts both
    limits of the integral, so its derivative carries the integrand at *both*
    ends. A formula written for independent limits passes every other test in
    this file and fails only here.
    """
    dq = DQ(*dq_params)
    params = list(dq_params)

    Z_anal, dZ_anal = element_jacobian(dq, freq, params)
    assert np.max(np.abs(Z_anal - dq.impedance(freq, params))) < 1e-15

    for j in range(4):
        p_fwd, p_bwd = params.copy(), params.copy()
        h = 1e-7 * abs(params[j])
        p_fwd[j] += h
        p_bwd[j] -= h
        dZ_num = (dq.impedance(freq, p_fwd) - dq.impedance(freq, p_bwd)) / (2 * h)

        scale = np.max(np.abs(dZ_num))
        rel_err = np.max(np.abs(dZ_anal[:, j] - dZ_num)) / scale
        assert rel_err < 1e-6, f"column {j} off by {rel_err}"


def test_dq_fit_runs_without_numeric_fallback(freq, dq_params):
    """A circuit containing DQ must fit: there is no numeric-Jacobian path.

    make_jacobian_function() only returns a closure, so a missing analytic
    branch surfaces as a RuntimeError out of least_squares rather than a
    fallback. This is the test that catches that.
    """
    A, n, tau_min, U = dq_params
    truth = R(20) - DQ(A, n, tau_min, U)
    Z = truth.impedance(freq, truth.get_all_params())

    guess = R(10) - DQ(5e5, 0.5, 1e-1, 6.0)
    result, _, _ = fit_equivalent_circuit(freq, Z, guess, plot=False)

    assert result.fit_error_rel < 0.1
    R_fit, A_fit, n_fit, tau_fit, U_fit = result.params_opt
    assert n_fit == pytest.approx(n, rel=1e-3)
    assert tau_fit == pytest.approx(tau_min, rel=1e-2)
    assert U_fit == pytest.approx(U, rel=1e-2)

    # Significance is None when any element lacks an analytic derivative
    assert result.params_significance is not None


def test_dq_bounds_keep_width_linear(dq_params):
    """All four DQ labels must be in PARAMETER_BOUNDS, U on a linear scale.

    An unknown label falls back silently to (1e-15, 1e15), 30 decades, which
    would put A_DQ's initial guess in the wrong place and make "U at its
    upper bound" - the report that says the distribution runs past the
    measured window - unreadable.
    """
    labels = DQ(*dq_params).get_param_labels()
    lower, upper = generate_simple_bounds(labels)

    assert log_scale_ci_mask(lower, upper) == [True, False, True, False]
    assert (lower[3], upper[3]) == (0.1, 30.0)
    assert upper[0] < 1e15, "A_DQ fell back to DEFAULT_BOUNDS"


def test_dq_parses_from_circuit_string():
    """The CLI reaches elements only through parse_circuit_expression()."""
    circuit = parse_circuit_expression("L(1e-6) - R(20) - DQ(1.2e6, 0.57, 5e-2, 8)")

    assert circuit.get_param_labels()[-4:] == ['A_DQ', 'n_DQ', 'τ_DQ', 'U_DQ']
    assert circuit.get_all_params()[-4:] == [1.2e6, 0.57, 5e-2, 8.0]


def _fit_result(circuit):
    """Minimal FitResult carrying a circuit, as the oxide analysis wants it."""
    params = np.array(circuit.get_all_params())
    return FitResult(circuit=circuit, params_opt=params,
                     params_stderr=np.zeros_like(params), fit_error_rel=0.1)


def test_dq_feeds_the_oxide_analysis(freq, dq_params):
    """The oxide layer must read C_eff off DQ directly, with no CPE model.

    A Q in the same place goes through Hsu-Mansfeld (and Brug) to guess an
    effective capacitance; DQ already has one as a limit of the fitted
    distribution, which is the point of using it on an oxide film.
    """
    A, n, tau_min, U = dq_params
    circuit = R(20) - DQ(A, n, tau_min, U)
    Z = circuit.impedance(freq, circuit.get_all_params())

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result(circuit))

    assert oxide.element_type == 'DQ'
    dq = DQ(A, n, tau_min, U)
    assert oxide.element_params['C'] == pytest.approx(dq.C_eff, rel=1e-12)
    assert oxide.element_R == pytest.approx(dq.R_pol, rel=1e-12)
    assert oxide.element_params['tau_max'] == pytest.approx(dq.tau_max, rel=1e-12)

    # Thickness from the plate-capacitor formula on that exact capacitance
    expected_d = 8.854e-14 * 22.0 * 1.0 / dq.C_eff  # cm
    assert oxide.thickness_nm == pytest.approx(expected_d * 1e7, rel=1e-3)


def test_dq_warns_when_the_plateau_is_out_of_window(dq_params):
    """C_eff above the highest measured frequency is an extrapolation."""
    A, n, _, U = dq_params
    tau_min = 1e-9  # plateau starts at ~1.6e8 Hz, far above the window
    freq = np.logspace(5, -3, 40)
    circuit = R(20) - DQ(A, n, tau_min, U)
    Z = circuit.impedance(freq, circuit.get_all_params())

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result(circuit))

    assert any('extrapolation' in w for w in oxide.warnings), oxide.warnings
