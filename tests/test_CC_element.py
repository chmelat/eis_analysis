#!/usr/bin/env python3
"""Test CC element (Cole-Cole dielectric relaxation, permittivity plane)."""

import numpy as np
import pytest
from eis_analysis.fitting import R, C, Q, CC, fit_equivalent_circuit
from eis_analysis.fitting.bounds import generate_simple_bounds, log_scale_ci_mask
from eis_analysis.fitting.jacobian import element_jacobian


@pytest.fixture
def freq():
    """Test frequencies: 1 MHz to 0.01 Hz."""
    return np.logspace(6, -2, 60)


@pytest.fixture
def cc_element_params():
    """CC element parameters: oxide-film-like dielectric."""
    C_inf = 1e-8
    dC = 1e-7
    tau = 1e-3
    alpha = 0.25
    return C_inf, dC, tau, alpha


def test_cc_impedance_matches_analytical(freq, cc_element_params):
    """Test CC element impedance matches the Cole-Cole formula."""
    C_inf, dC, tau, alpha = cc_element_params

    cc_elem = CC(C_inf, dC, tau, alpha)
    Z_CC = cc_elem.impedance(freq, [C_inf, dC, tau, alpha])

    omega = 2 * np.pi * freq
    C_star = C_inf + dC / (1 + (1j * omega * tau) ** (1 - alpha))
    Z_expected = 1 / (1j * omega * C_star)

    max_diff = np.max(np.abs(Z_CC - Z_expected))
    assert max_diff < 1e-10, f"Impedance differs from analytical: {max_diff}"


def test_cc_equals_composite_circuit(freq, cc_element_params):
    """CC(C_inf, dC, tau, alpha) == C(C_inf) | (C(dC) - Q(dC/tau^(1-a), a)).

    The Cole-Cole element is exactly realisable from existing elements; this
    is the oracle that pins down the impedance formula. The composite is not
    a substitute for fitting -- its CPE coefficient couples three parameters
    non-linearly -- but it is an independent implementation of the same maths.
    """
    C_inf, dC, tau, alpha = cc_element_params

    cc_elem = CC(C_inf, dC, tau, alpha)
    Z_CC = cc_elem.impedance(freq, [C_inf, dC, tau, alpha])

    composite = C(C_inf) | (C(dC) - Q(dC / tau ** (1 - alpha), alpha))
    Z_composite = composite.impedance(freq, composite.get_all_params())

    rel_err = np.max(np.abs(Z_CC - Z_composite)) / np.max(np.abs(Z_CC))
    assert rel_err < 1e-12, f"CC differs from equivalent composite: {rel_err}"


def test_cc_alpha_zero_is_debye(freq, cc_element_params):
    """alpha = 0 reduces to Debye: C(C_inf) | (C(dC) - R(tau/dC))."""
    C_inf, dC, tau, _ = cc_element_params

    cc_elem = CC(C_inf, dC, tau, 0.0)
    Z_CC = cc_elem.impedance(freq, [C_inf, dC, tau, 0.0])

    debye = C(C_inf) | (C(dC) - R(tau / dC))
    Z_debye = debye.impedance(freq, debye.get_all_params())

    rel_err = np.max(np.abs(Z_CC - Z_debye)) / np.max(np.abs(Z_CC))
    assert rel_err < 1e-12, f"alpha=0 is not Debye: {rel_err}"


def test_cc_limiting_behavior(cc_element_params):
    """High frequency -> pure C_inf; low frequency -> pure C_static."""
    C_inf, dC, tau, alpha = cc_element_params
    cc_elem = CC(C_inf, dC, tau, alpha)
    params = [C_inf, dC, tau, alpha]

    # Well above 1/tau the relaxing branch is shorted out: Z -> 1/(jw*C_inf)
    f_high = np.array([1e9])
    Z_high = cc_elem.impedance(f_high, params)
    Z_c_inf = 1 / (1j * 2 * np.pi * f_high * C_inf)
    assert np.abs(Z_high[0] - Z_c_inf[0]) / np.abs(Z_c_inf[0]) < 0.01

    # Well below 1/tau the dielectric is fully relaxed: Z -> 1/(jw*C_static)
    f_low = np.array([1e-8])
    Z_low = cc_elem.impedance(f_low, params)
    Z_c_static = 1 / (1j * 2 * np.pi * f_low * cc_elem.C_static)
    assert np.abs(Z_low[0] - Z_c_static[0]) / np.abs(Z_c_static[0]) < 0.01

    # Blocking dielectric: |Z| diverges as omega -> 0
    assert np.abs(Z_low[0]) > np.abs(Z_high[0])


def test_cc_element_properties(cc_element_params):
    """Test CC computed properties."""
    C_inf, dC, tau, alpha = cc_element_params
    cc_elem = CC(C_inf, dC, tau, alpha)

    assert cc_elem.C_static == pytest.approx(C_inf + dC, rel=1e-12)
    assert cc_elem.characteristic_freq == pytest.approx(1 / (2 * np.pi * tau), rel=1e-12)


def test_cc_bounds_keep_alpha_linear(cc_element_params):
    """The broadening exponent must not be searched/CI'd in log space.

    generate_simple_bounds() falls back to (1e-15, 1e15) for unknown labels,
    which spans 30 decades and would classify alpha as a positive scale
    parameter. This test fails if the Cole-Cole labels go missing from
    PARAMETER_BOUNDS.
    """
    labels = CC(*cc_element_params).get_param_labels()
    lower, upper = generate_simple_bounds(labels)

    assert log_scale_ci_mask(lower, upper) == [True, True, True, False]
    assert (lower[3], upper[3]) == (0.0, 0.9)


def test_cc_jacobian(freq, cc_element_params):
    """Analytic Jacobian matches central differences."""
    C_inf, dC, tau, alpha = cc_element_params
    cc_elem = CC(C_inf, dC, tau, alpha)
    params = [C_inf, dC, tau, alpha]

    Z_jac, dZ_jac = element_jacobian(cc_elem, freq, params)
    assert np.max(np.abs(Z_jac - cc_elem.impedance(freq, params))) < 1e-15

    # Looser thresholds on tau/alpha: their finite-difference truncation error
    # dominates (~1e-5), the analytic formula itself is exact.
    tolerances = [1e-5, 1e-5, 1e-4, 1e-4]
    for i, tol in enumerate(tolerances):
        h = 1e-8 * max(abs(params[i]), 1e-9)
        p_plus, p_minus = list(params), list(params)
        p_plus[i] += h
        p_minus[i] -= h
        dZ_num = (cc_elem.impedance(freq, p_plus) -
                  cc_elem.impedance(freq, p_minus)) / (2 * h)
        rel_err = np.max(np.abs(dZ_jac[:, i] - dZ_num)) / np.max(np.abs(dZ_jac[:, i]))
        assert rel_err < tol, f"Jacobian column {i} off by {rel_err}"


def test_cc_element_fitting(freq):
    """Round-trip fit on synthetic Cole-Cole data with 1% noise."""
    R_s_true = 20.0
    C_inf_true, dC_true, tau_true, alpha_true = 2e-9, 5e-8, 2e-3, 0.30

    omega = 2 * np.pi * freq
    C_star = C_inf_true + dC_true / (1 + (1j * omega * tau_true) ** (1 - alpha_true))
    Z_true = R_s_true + 1 / (1j * omega * C_star)

    np.random.seed(42)
    noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) +
                                     1j * np.random.randn(len(freq)))
    Z_noisy = Z_true + noise

    # Deliberately offset initial guess
    circuit = R(10) - CC(1e-9, 1e-7, 5e-4, 0.15)
    result, _, _ = fit_equivalent_circuit(freq, Z_noisy, circuit,
                                          weighting='modulus', plot=False)

    assert result.fit_error_rel < 2.0, f"Fit error too high: {result.fit_error_rel}"

    _, C_inf_fit, dC_fit, tau_fit, alpha_fit = result.params_opt
    assert abs(C_inf_fit - C_inf_true) / C_inf_true < 0.15
    assert abs(dC_fit - dC_true) / dC_true < 0.15
    assert abs(tau_fit - tau_true) / tau_true < 0.15
    assert abs(alpha_fit - alpha_true) < 0.05


def test_cc_fixed_parameters():
    """String parameters are fixed, and scaling preserves that."""
    cc_elem = CC("1e-8", 1e-7, 1e-3, "0.2")
    assert cc_elem.fixed_params == [True, False, False, True]

    # _scale must not silently free the fixed parameters (R/Q/K do lose them)
    scaled = 2 * cc_elem
    assert scaled.fixed_params == [True, False, False, True]
    assert scaled.C_inf == pytest.approx(2e-8)
    assert scaled.dC == pytest.approx(2e-7)
    assert scaled.tau == pytest.approx(1e-3)
    assert scaled.alpha == pytest.approx(0.2)
