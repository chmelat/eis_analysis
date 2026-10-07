"""
Tests for R_inf estimation (rinf_estimation/estimate.py).

Cases follow doc/archive/AUDIT_ri_fit_2026-09-25.md, where the former three-branch
estimator (zero crossing / polynomial / R-L-K with fixed tau) was off by
+10 % to +100 % on inductive spectra (C, C2, A3) and by +27 % on a strongly
open CPE arc (D), always with fit_success=True and no warning.
"""

import numpy as np
import pytest

from eis_analysis.io import load_csv_data
from eis_analysis.rinf_estimation import estimate_rinf


def _zarc(f, R0, tau, n=1.0):
    return R0 / (1 + (1j * 2 * np.pi * f * tau) ** n)


def _spectrum(Rs, L, R0, tau, n=1.0, f_max=1e6):
    """Rs + jwL + ZARC, 10 points/decade down to 0.1 Hz."""
    f = np.logspace(np.log10(f_max), -1, int(10 * (np.log10(f_max) + 1)) + 1)
    return f, Rs + 1j * 2 * np.pi * f * L + _zarc(f, R0, tau, n)


def _f_arc(f_peak):
    return 1 / (2 * np.pi * f_peak)


# name -> (Rs, spectrum). Rs = 10 Ohm, R = 100 Ohm unless stated.
CASES = {
    'flat': (10, _spectrum(10, 0, 100, _f_arc(16))),
    'A2 separate arc': (10, _spectrum(10, 1e-6, 100, _f_arc(1.6e3))),
    'A3 L + arc 160 kHz': (10, _spectrum(10, 1e-5, 100, _f_arc(160e3))),
    'C L + RC 16 kHz': (10, _spectrum(10, 1e-5, 100, 1e-5)),
    'C2 L + ZARC 0.7': (10, _spectrum(10, 1e-5, 100, 1e-5, 0.7)),
    'B n=0.7 capacitive': (10, _spectrum(10, 0, 100, 1e-4, 0.7, f_max=1e5)),
    'D Rs=0.05, n=0.6': (0.05, _spectrum(0.05, 0, 100, 1e-3, 0.6, f_max=1e5)),
}


def _noisy(Z, seed, level=0.01):
    rng = np.random.default_rng(seed)
    return Z * (1 + level * (rng.standard_normal(len(Z))
                             + 1j * rng.standard_normal(len(Z))))


@pytest.mark.parametrize('name', CASES)
def test_noiseless_case_recovers_rs(name):
    Rs, (f, Z) = CASES[name]
    res = estimate_rinf(f, Z)
    assert res.method == 'rlq_fit', res.warnings
    assert res.R_inf == pytest.approx(Rs, rel=0.02)
    assert res.warnings == []


@pytest.mark.parametrize('name', ['C L + RC 16 kHz', 'B n=0.7 capacitive'])
def test_noisy_determinable_case_uses_fit(name):
    # 1 % noise: the audit's C and B cases stay within a few percent.
    Rs, (f, Z) = CASES[name]
    res = estimate_rinf(f, _noisy(Z, seed=1))
    assert res.method == 'rlq_fit', res.warnings
    assert res.R_inf == pytest.approx(Rs, rel=0.05)


def _undeterminable(name):
    """Arc entirely above f_max (A1), strongly open CPE arc with noise (D),
    two overlapping CPEs the model does not describe (example CSV)."""
    if name == 'A1':
        f, Z = _spectrum(10, 1e-7, 100, _f_arc(16e6))
        return 10, f, _noisy(Z, seed=0)
    if name == 'D':
        Rs, (f, Z) = CASES['D Rs=0.05, n=0.6']
        return Rs, f, _noisy(Z, seed=0)
    csv = load_csv_data('example/example_eis_data.csv')
    return 10, csv.frequencies, csv.Z


@pytest.mark.parametrize('name', ['A1', 'D', 'CSV'])
def test_undeterminable_falls_back_to_upper_bound_with_warning(name):
    Rs, f, Z = _undeterminable(name)
    res = estimate_rinf(f, Z)
    assert res.method == 'hf_bound'
    # Re(Z) at f_max: an upper bound, tighter than the 5-point HF median
    # on an open arc (D: +2483 % instead of +3295 %).
    assert res.R_inf == Z.real[np.argmax(f)] >= Rs
    assert res.fit is not None and 'not determined' in res.warnings[-1]


def test_accepts_lists_and_drops_non_finite_points():
    _, (f, Z) = CASES['flat']
    Z_bad = Z.copy()
    Z_bad[5] = np.nan
    res = estimate_rinf(list(f), list(Z_bad))
    assert res.R_inf == pytest.approx(10, rel=0.02)
    assert 'non-finite' in res.warnings[0]


def test_rejects_shape_mismatch():
    _, (f, Z) = CASES['flat']
    with pytest.raises(ValueError):
        estimate_rinf(f, Z[:-1])


def test_too_few_window_points_fall_back_without_fit():
    # Not an underdetermined fit: 3 points in the window, 5 parameters.
    _, (f, Z) = CASES['flat']
    res = estimate_rinf(f[::10], Z[::10])
    assert res.method == 'hf_bound' and res.fit is None
    assert 'need >=' in res.warnings[0]
