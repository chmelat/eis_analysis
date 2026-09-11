#!/usr/bin/env python3
"""Test DQ element (bounded power-law DRT, truncated CPE)."""

import numpy as np
import pytest
from scipy.integrate import quad

from eis_analysis.fitting import K, Q, DQ


@pytest.fixture
def freq():
    """Test frequencies: 1 MHz to 1 mHz, the oxide measurement window."""
    return np.logspace(6, -3, 60)


@pytest.fixture
def dq_params():
    """Oxide-like parameters: n = 0.57, distribution inside the window."""
    return 1e-3, 0.57, 5e-2, 8.0


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
