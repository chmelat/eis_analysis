"""
Tests for R_inf estimation (rinf_estimation/estimate.py).

Cases follow doc/AUDIT_ri_fit_2026-09-25.md, where the former three-branch
estimator (zero crossing / polynomial / R-L-K with fixed tau) was off by
+10 % to +100 % on inductive spectra (C, C2, A3) and by +27 % on a strongly
open CPE arc (D), always with fit_success=True and no warning.
"""

import numpy as np
import pytest

from eis_analysis.io.data_loading import load_csv_data
from eis_analysis.rinf_estimation import estimate_rinf


def _zarc(f, R0, tau, n=1.0):
    return R0 / (1 + (1j * 2 * np.pi * f * tau) ** n)


def _spectrum(Rs, L, R0, tau, n=1.0, f_max=1e6):
    """Rs + jwL + ZARC, 10 points/decade down to 0.1 Hz."""
    f = np.logspace(np.log10(f_max), -1, int(10 * (np.log10(f_max) + 1)) + 1)
    return f, Rs + 1j * 2 * np.pi * f * L + _zarc(f, R0, tau, n)


def _f_arc(f_peak):
    return 1 / (2 * np.pi * f_peak)


# (name, Rs, spectrum). Rs = 10 Ohm, R = 100 Ohm unless stated.
CASES = [
    ('flat', 10, _spectrum(10, 0, 100, _f_arc(16))),
    ('A2 separate arc', 10, _spectrum(10, 1e-6, 100, _f_arc(1.6e3))),
    ('A3 L + arc 160 kHz', 10, _spectrum(10, 1e-5, 100, _f_arc(160e3))),
    ('C L + RC 16 kHz', 10, _spectrum(10, 1e-5, 100, 1e-5)),
    ('C2 L + ZARC 0.7', 10, _spectrum(10, 1e-5, 100, 1e-5, 0.7)),
    ('B n=0.7 capacitive', 10, _spectrum(10, 0, 100, 1e-4, 0.7, f_max=1e5)),
    ('D Rs=0.05, n=0.6', 0.05, _spectrum(0.05, 0, 100, 1e-3, 0.6, f_max=1e5)),
]


def _noisy(Z, seed, level=0.01):
    rng = np.random.default_rng(seed)
    return Z * (1 + level * (rng.standard_normal(len(Z))
                             + 1j * rng.standard_normal(len(Z))))


def test_noiseless_cases_recover_rs():
    for name, Rs, (f, Z) in CASES:
        res = estimate_rinf(f, Z)
        assert res.method == 'rlq_fit', (name, res.warnings)
        assert res.R_inf == pytest.approx(Rs, rel=0.02), name
        assert res.warnings == [], name


def test_noisy_determinable_cases_use_fit():
    # 1 % noise: the audit's C and B cases stay within a few percent.
    for name, Rs, (f, Z) in (CASES[3], CASES[5]):
        res = estimate_rinf(f, _noisy(Z, seed=1))
        assert res.method == 'rlq_fit', (name, res.warnings)
        assert res.R_inf == pytest.approx(Rs, rel=0.05), name


def test_undeterminable_falls_back_to_median_with_warning():
    # Arc entirely above f_max (A1), strongly open CPE arc with noise (D),
    # and two overlapping CPEs the model does not describe (example CSV).
    f_a1, Z_a1 = _spectrum(10, 1e-7, 100, _f_arc(16e6))
    csv = load_csv_data('example/example_eis_data.csv')
    for name, f, Z in (('A1', f_a1, _noisy(Z_a1, seed=0)),
                       ('D', CASES[6][2][0], _noisy(CASES[6][2][1], seed=0)),
                       ('CSV', csv.frequencies, csv.Z)):
        res = estimate_rinf(f, Z)
        assert res.method == 'hf_median', name
        assert res.R_inf == res.R_inf_median, name
        assert res.fit is not None and 'not determined' in res.warnings[-1], name


def test_input_handling():
    f, Z = CASES[0][2]
    # Plain lists are accepted; non-finite points are dropped and reported.
    Z_bad = Z.copy()
    Z_bad[5] = np.nan
    res = estimate_rinf(list(f), list(Z_bad))
    assert res.R_inf == pytest.approx(10, rel=0.02)
    assert 'non-finite' in res.warnings[0]
    with pytest.raises(ValueError):
        estimate_rinf(f, Z[:-1])
    # Too few points in the window: median, not an underdetermined fit.
    res = estimate_rinf(f[::10], Z[::10])  # 3 points in the window
    assert res.method == 'hf_median' and res.fit is None
    assert 'need >=' in res.warnings[0]
