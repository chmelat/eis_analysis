"""
Tests for R_inf estimation (rinf_estimation/rlk_fit.py).

Regression tests for a P1 finding (2026-07-26): the purely-capacitive branch
fit a 2nd-degree polynomial (Re = f(Im)) to extrapolate R_inf, but never
checked that the top-frequency-decade selection returned enough points for
that fit. With fewer than 3 points the polynomial is underdetermined and
np.polyfit silently returns a least-norm solution that looks like a valid
fit but isn't. Reproduced on synthetic data with a known Rs = 5.0 Ohm and 2
points in the top decade: the unfixed code reported R_inf = 4.32 Ohm (-14%)
with fit_success=True and an empty warnings list. The fix falls back to the
measured highest-frequency Re(Z) and records a warning when the point count
is too low for the polynomial fit.
"""

import numpy as np

from eis_analysis.rinf_estimation import estimate_rinf_with_inductance


def _capacitive_Z(f, Rs=5.0, Rct=100.0, C=1e-5):
    """Synthetic purely-capacitive Rs-Rct||C spectrum with known Rs."""
    return Rs + Rct / (1 + 1j * 2 * np.pi * f * Rct * C)


def test_capacitive_few_points_falls_back_to_hf():
    # Sparse grid: only 2 points fall in the top frequency decade
    # (>= f_max / 10 = 1e3), so the polynomial fit is underdetermined.
    f = np.array([1e-2, 1e-1, 1, 2, 5, 10, 20, 50, 100, 200, 6e3, 1e4])
    Z = _capacitive_Z(f)

    res, _ = estimate_rinf_with_inductance(f, Z)

    assert abs(res.R_inf - 5.0) / 5.0 < 0.02, (
        f"Expected R_inf within 2% of true Rs=5.0, got {res.R_inf}"
    )
    assert 'hf_too_few_points' in res.method, res.method
    assert res.warnings, "Expected a non-empty warnings list"


def test_capacitive_enough_points_uses_polynomial():
    # Dense grid: at least 3 points land in the top frequency decade,
    # so the polynomial extrapolation path is exercised normally.
    f = np.logspace(-2, 4, 40)
    Z = _capacitive_Z(f)

    res, _ = estimate_rinf_with_inductance(f, Z)

    assert 'polynomial' in res.method or 'hf_fallback' in res.method, res.method
    assert res.R_inf > 0
    assert res.poly_coeffs is not None


def test_few_points_result_has_no_poly_coeffs():
    f = np.array([1e-2, 1e-1, 1, 2, 5, 10, 20, 50, 100, 200, 6e3, 1e4])
    Z = _capacitive_Z(f)

    res, _ = estimate_rinf_with_inductance(f, Z)

    assert res.poly_coeffs is None
    assert res.R_inf_poly is None


def test_inductive_data_unaffected():
    # Top decade is inductive (Im(Z) > 0), so the capacitive guard must not
    # affect this branch at all.
    f = np.logspace(-2, 4, 40)
    Z = 5.0 + 1j * 2 * np.pi * f * 1e-7

    res, _ = estimate_rinf_with_inductance(f, Z)

    assert res.method.startswith('rlk_linear') or res.method.startswith('zero_crossing'), (
        res.method
    )
