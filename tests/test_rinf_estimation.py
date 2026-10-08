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


def test_negative_upper_bound_is_clipped_to_zero():
    """At a point dominated by noise Re(Z) can be negative; R_s >= 0, so the
    bound is 0, not the negative value (stress test, constant noise:
    rc/196 gave R_inf = -1.9e4 Ohm against R_s = 2.0e4)."""
    f, Z = _spectrum(Rs=5, L=0, R0=100, tau=_f_arc(1e3))
    Z[0] = -3.0 - 1.0j    # f_max swamped by noise
    est = estimate_rinf(f[[0, 20, 30, 40, 50]], Z[[0, 20, 30, 40, 50]])
    assert est.method == 'hf_bound'
    assert est.R_inf == 0.0
    assert est.R_inf_hf == -3.0
    assert 'negative' in est.warnings[-1]
    # No fit ran to read the noise, yet the range must not collapse to (0, 0)
    lo, hi = est.R_inf_range
    assert lo == 0.0 and hi >= 5


def test_range_of_a_noisy_bound_holds_rs():
    """Constant noise at ~0.6x |Z(f_max)|: the window fit cannot determine
    R_s and the top point reads -4e3 Ohm against R_s = 2e4. R_inf_range,
    which local_exponent takes for its R_inf sensitivity, must still hold
    R_s: (0, 0) from the clipped bound made 68 of 81 n(f) points
    'determined' that were up to 0.79 off (code review)."""
    f = np.logspace(5, -1, 61)
    w = 2 * np.pi * f
    Z = 2e4 + 1 / (1e-9 * (1j * w) ** 0.85 + 1j * w * 2e-11)
    rng = np.random.default_rng(8)
    Z = Z + 1.5e4 * (rng.normal(size=61) + 1j * rng.normal(size=61))
    est = estimate_rinf(f, Z)
    assert est.method == 'hf_bound' and est.R_inf_hf < 0 and est.R_inf == 0.0
    lo, hi = est.R_inf_range
    assert lo == 0.0 and hi >= 2e4


# Spectra of the stress test (tests/stress.py) that broke unit invariance:
# (expression, f_max, f_min, points), noise-free
UNIT_CASES = {
    # open film arc in the window: R_k ran into R <= 1e10 Ohm at Z x 1000,
    # the fit failed and R_inf fell back to the HF bound, 1.2e9 instead of 1.6e3
    'oxide/10': ('R(1.5060315555665122) - (R(24856893.804588407)'
                 '|Q(4.7578522543970145e-11,0.9080075824366625))',
                 1858.4534441671751, 6.656689604430729, 33),
    # blocking tail: LM stopped on xtol, which mixes the parameters' units,
    # at a different point for Z and for exactly Z/2
    'blocking/38': ('R(4254241.56646422) - (R(232721572.02547705)'
                    '|Q(1.838929474718818e-10,0.6633267992074743))'
                    ' - Q(4.0252334682656605e-11,0.8996571087461844)',
                    107962.72875713777, 0.24267724373552552, 52),
}


@pytest.mark.parametrize('name', list(UNIT_CASES))
@pytest.mark.parametrize('k', [0.5, 1e-3, 1e3])
def test_rinf_does_not_depend_on_units(name, k):
    """Z -> k*Z gives k*R_inf by the same method.

    Regression: the window fit used the absolute PARAMETER_BOUNDS and
    least_squares' unit-mixing stop criteria; R_inf shifted on 51 % of the
    stress test's random spectra, by over 10 % on 7 % of them.
    """
    from eis_analysis.cli.utils import parse_circuit_expression

    expression, f_max, f_min, n = UNIT_CASES[name]
    f = np.logspace(np.log10(f_max), np.log10(f_min), n)
    circuit = parse_circuit_expression(expression)
    Z = circuit.impedance(f, circuit.get_all_params())

    ref, scaled = estimate_rinf(f, Z), estimate_rinf(f, k * Z)
    assert scaled.method == ref.method
    assert scaled.R_inf == pytest.approx(k * ref.R_inf, rel=1e-6)


def test_small_rinf_in_front_of_large_film():
    """Noise-free 0.05 Ohm R_s in front of a GOhm film, |Z| 3e3..2.5e5 Ohm.

    Regression (code review of the unit fix): lower bounds of R and L at
    1e-6 of the window's |Z| biased R_s to 0.064 Ohm, and the absolute gtol
    on the now dimensionless cost stopped the fit early; R_inf fell back to
    the HF bound, 243 Ohm.
    """
    from eis_analysis.cli.utils import parse_circuit_expression

    f = np.logspace(5, 1, 41)
    circuit = parse_circuit_expression("R(0.05)-(R(1e9)|Q(1e-9,0.95))")
    result = estimate_rinf(f, circuit.impedance(f, circuit.get_all_params()))
    assert result.method == 'rlq_fit'
    assert result.R_inf == pytest.approx(0.05, rel=1e-3)
