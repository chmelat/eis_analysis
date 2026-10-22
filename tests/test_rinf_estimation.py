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
    assert res.notes == []


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



def test_range_covers_the_model_error_on_an_open_end():
    """A CPE arc peaking at 10 kHz, still open at f_max = 100 kHz: the
    two-decade window fit is accepted 10 % below R_s = 10 Ohm (no noise),
    the top decade's 1.2 % below with a far smaller stderr, so it is used.

    Regression: R_inf_range was R_inf +- max(5 %, 3 stderr), 8.36..9.57 Ohm,
    which the model error escapes; the stress test found R_s outside it on
    39 of 1417 spectra. The window spread now widens it.
    """
    f = np.logspace(5, -2, 71)
    jw = 2j * np.pi * f
    Z = 10 + 1 / (1e-4 + 5.26e-7 * jw**0.6) + 1 / (1e-4 + 1.59e-5 * jw)
    est = estimate_rinf(f, Z)
    assert est.method == 'rlq_fit'
    assert est.f_window.min() >= f.max() / 10 * (1 - 1e-12)   # the top decade
    assert abs(est.R_inf - 10) < 0.02 * 10
    assert est.R_inf_window_spread > 0.05 * est.R_inf      # the full window's model error
    lo, hi = est.R_inf_range
    assert lo <= 10 <= hi
    assert hi <= est.R_inf_upper


def test_range_stays_tight_on_a_closed_end():
    """A closed arc (peak at 160 Hz, f_max 100 kHz): the narrow window gives
    R_s too, so the window spread adds nothing to the +-5 % floor, and the
    upper end is Re Z(f_max), the tighter physical bound."""
    f = np.logspace(5, -2, 71)
    est = estimate_rinf(f, 10 + 1000 / (1 + 2j * np.pi * f * 1e-3))
    assert est.method == 'rlq_fit'
    assert est.R_inf_window_spread < 1e-3 * est.R_inf
    lo, hi = est.R_inf_range
    assert (lo, hi) == pytest.approx((0.95 * est.R_inf, est.R_inf_upper), rel=1e-9)
    assert est.R_inf_upper < 1.05 * est.R_inf


def test_range_ignores_an_arc_below_the_window():
    """Two arcs two decades apart, the top end closed (code review): a
    3-decade comparison window reached into the second arc and widened the
    range to 3.5..10.2 Ohm; the narrow window stays inside the first."""
    f = np.logspace(7, -2, 91)
    jw = 2j * np.pi * f
    est = estimate_rinf(f, 10 + 100 / (1 + jw * 1e-6) + 1000 / (1 + jw * 1e-4))
    assert est.method == 'rlq_fit'
    lo, hi = est.R_inf_range
    assert lo <= 10 <= hi
    assert lo > 0.9 * est.R_inf


def test_zero_point_below_the_window_does_not_break_the_spread():
    """Z = 0 at a point 2.5 decades below f_max (code review: the 3-decade
    window then divided by min|Z| = 0 and raised ZeroDivisionError)."""
    f = np.logspace(5, -2, 71)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * 1e-3)
    Z[25] = 0
    est = estimate_rinf(f, Z)
    assert est.method == 'rlq_fit'
    assert est.R_inf == pytest.approx(10, rel=1e-6)


def test_duplicate_top_frequency_gives_the_lowest_re_in_any_order():
    """f_max measured twice, Re 1 Ohm apart, on an arc still open at the top
    (contract review round 5, F1): the HF bound came from whichever point
    came first, so R_inf moved by 1 Ohm with the order of the input. Every
    point is a bound; the lowest is the tightest."""
    f = np.logspace(5, -1, 61)
    Z = 10 + 100 / (1 + 100 * 1e-8 * (2j * np.pi * f) ** 0.7)
    f2, Z2 = np.concatenate([[f[0]], f]), np.concatenate([[Z[0] + 1], Z])
    results = [estimate_rinf(f2[p], Z2[p]) for p in
               (np.arange(len(f2)), np.arange(len(f2))[::-1],
                np.random.default_rng(3).permutation(len(f2)))]
    for est in results:
        assert est.method == 'hf_bound'
        assert est.R_inf_hf == Z[0].real
        assert est.R_inf == results[0].R_inf
        assert est.R_inf_range == pytest.approx(results[0].R_inf_range, rel=1e-9)


def test_all_inductive_spectrum_warns_that_there_is_no_bound():
    """No point with Im(Z) <= 0 (round 5, F2): Re(Z) at f_max was reported
    as the HF upper bound, though it lies below R_s = 10 here."""
    for n in (4, 8):    # fallback without a fit, and after the window fit
        f = np.logspace(6, 5.5, n)
        Z = 10 - np.linspace(1, 0.1, n) + 2j * np.pi * f * 1e-7
        est = estimate_rinf(f, Z)
        assert not est.hf_is_bound and est.method == 'hf_bound'
        assert 'No point with Im(Z) <= 0' in est.warnings[0]
        assert est.warnings[-1].endswith('no bound')
        assert not any('upper bound Re(Z)' in w for w in est.warnings)


def test_inductive_point_at_the_top_frequency_is_noted_as_skipped():
    """f_max measured twice, once inductive with the lower Re (code review
    of round 5): the bound is the capacitive point, and the warning says
    an inductive point was skipped although f_hf == f_max."""
    f = np.array([1e5, 1e5, 1e4])
    Z = np.array([9 + 0.5j, 11 - 0.1j, 12 - 1j])
    est = estimate_rinf(f, Z)
    assert est.hf_is_bound and est.R_inf_hf == 11 and est.f_hf == 1e5
    assert 'inductive top above Im(Z) = 0 skipped' in est.warnings[-1]


def test_hf_r_inf_does_not_depend_on_the_order():
    """hf_r_inf, calculate_rpol's R_inf (the DRT's default until v0.61; code
    review of round 5): with three points at f_max its two 'top' points
    depended on the input order, and with them whether the median or the
    bound was returned."""
    from eis_analysis.rinf_estimation import hf_r_inf

    f = np.concatenate([[1e5, 1e5, 1e5], np.logspace(4.9, 0, 50)])
    Z = np.concatenate([[10 - 1j, 11 - 1j, 14 - 1j], 40 - 5j + 0 * f[3:]])
    ref = hf_r_inf(f, Z)
    for seed in range(20):
        p = np.random.default_rng(seed).permutation(len(f))
        assert hf_r_inf(f[p], Z[p]) == ref


def test_fitted_rs_at_zero_falls_back(monkeypatch):
    """R_s on its lower bound 0 has no relative stderr (round 5, F3):
    stderr / R_s raised ZeroDivisionError instead of the fallback."""
    import eis_analysis.rinf_estimation.estimate as module

    fit_window = module._fit_window

    def fit_at_zero(f_win, Z_win):
        fit, Z_fit = fit_window(f_win, Z_win)
        fit.params_opt = fit.params_opt.copy()
        fit.params_opt[0] = 0.0
        fit.params_stderr = fit.params_stderr.copy()
        fit.params_stderr[0] = 0.0
        return fit, Z_fit

    f = np.logspace(5, -2, 71)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * 1e-3)
    monkeypatch.setattr(module, '_fit_window', fit_at_zero)
    est = estimate_rinf(f, Z)
    assert est.method == 'hf_bound' and est.fit is not None
    assert 'not determined' in est.warnings[-1]


def _undetermined_fit(monkeypatch):
    """Make the window fit report R_s with a 100 % stderr."""
    import eis_analysis.rinf_estimation.estimate as module

    fit_window = module._fit_window

    def fit_undetermined(f_win, Z_win):
        fit, Z_fit = fit_window(f_win, Z_win)
        fit.params_stderr = fit.params_stderr.copy()
        fit.params_stderr[0] = fit.params_opt[0]
        return fit, Z_fit

    monkeypatch.setattr(module, '_fit_window', fit_undetermined)


def test_closed_end_notes_instead_of_warning(monkeypatch):
    """|Im| 1.6 % of Re at f_max: the bound needs no fit, so an undetermined
    R_s is a note, not the 'not determined' warning; the result is the
    same fallback."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * 1e-2)
    est = estimate_rinf(f, Z)
    assert est.method == 'hf_bound' and est.warnings == []
    assert est.notes[-1].startswith('High-frequency end closed')
    assert est.R_inf == Z.real[0] and est.L == 0.0
    assert est.R_inf_range == (0.0, est.R_inf_upper)


@pytest.mark.parametrize('im_share', [0.10, 0.06])
def test_open_end_still_warns(monkeypatch, im_share):
    """|Im| above 5 % of Re at f_max: the bound may be off, so warn."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    tau = 1000 / (2 * np.pi * 1e5 * 10 * im_share)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * tau)
    est = estimate_rinf(f, Z)
    assert est.method == 'hf_bound' and est.notes == []
    assert 'not determined' in est.warnings[-1]


def test_inductive_top_still_warns(monkeypatch):
    """An L cancels the CPE's Im at f_hf (~40 kHz): |Im(f_hf)| is ~0.02 %
    of Re, yet the bound is 17 % above Rs (n = 0.6 leaves Re), so warn."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    Z = 10 + 2j * np.pi * f * 9.15e-6 + 1000 / (1 + (2j * np.pi * f * 0.07) ** 0.6)
    est = estimate_rinf(f, Z)
    i_hf = int(np.argmin(np.abs(f - est.f_hf)))
    assert est.f_hf < f.max() and abs(Z.imag[i_hf]) < 0.01 * est.R_inf_hf
    assert est.R_inf_hf > 1.1 * 10
    assert est.method == 'hf_bound' and est.notes == []
    assert 'not determined' in est.warnings[-1]


def test_inductive_point_within_noise_still_warns(monkeypatch):
    """Code review: an L-cancelled end passed when its inductive points
    stayed below 3x the rejected fit's noise. Any point with Im > 0 above
    f_hf now blocks the note, however small its Im."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * 1e-2)
    Z[0] = Z[0].real + 1e-6j
    est = estimate_rinf(f, Z)
    assert est.f_hf == f[1]
    assert est.method == 'hf_bound' and est.notes == []
    assert 'not determined' in est.warnings[-1]


def test_noisy_closed_end_still_warns(monkeypatch):
    """Code review: one noisy |Im| can fall below 5 % of Re on an end that
    is not closed. |Im| 1.6 % plus 3x a 2 % noise is above 5 %: warn."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    Z = _noisy(10 + 1000 / (1 + 2j * np.pi * f * 1e-2), seed=0, level=0.02)
    Z[0] = Z[0].real - 0.001j  # a measured |Im| far below 5 %
    est = estimate_rinf(f, Z)
    assert est.f_hf == f[0]
    assert est.method == 'hf_bound' and est.notes == []
    assert 'not determined' in est.warnings[-1]


def test_inductive_duplicate_at_f_hf_still_warns(monkeypatch):
    """Code review: a repeated f_max with Im > 0 next to the capacitive
    point hf_bound picks may be an L cancelling the arc there: warn."""
    _undetermined_fit(monkeypatch)
    f = np.logspace(5, -2, 71)
    Z = 10 + 1000 / (1 + 2j * np.pi * f * 1e-2)
    f, Z = np.append(f, f[0]), np.append(Z, Z[0].real + 2j)
    est = estimate_rinf(f, Z)
    assert est.f_hf == f[0]
    assert est.method == 'hf_bound' and est.notes == []
    assert 'not determined' in est.warnings[-1]


def _rc143(R2=3.3830971543626727, C2=0.004788848173323314):
    """Stress case rc/143 (noise-free): R_s, a series L, two (R|C) arcs."""
    f = np.logspace(np.log10(15614.98), np.log10(1.2204), 30)
    jw = 2j * np.pi * f
    Rs = 1.2389074735828003
    Z = (Rs + jw * 8.188237996267678e-08
         + 0.4770687162622361 / (1 + jw * 0.4770687162622361 * 0.0002792966695560041)
         + R2 / (1 + jw * R2 * C2))
    return f, Z, Rs

def test_top_decade_fit_wins_over_a_second_arc_in_the_window():
    """rc/143 of the stress test (noise-free): R_s 1.239, a series L and two
    arcs, the slower one reaching into the two-decade window. One (R|Q)
    took both, n 0.70, R_s 5.5 % low at 2.6 % stderr; the top decade alone
    gives R_s within 0.2 % at 0.05 %, so it is used."""
    f, Z, Rs = _rc143()
    est = estimate_rinf(f, Z)
    assert est.method == 'rlq_fit' and est.fit_window == 'top'
    assert abs(est.R_inf - Rs) < 0.002 * Rs
    assert est.f_window.min() >= f.max() / 10 * (1 - 1e-12)
    assert est.R_inf_window_spread > 0.05 * Rs
    # The spread is the discarded fit's error: the range is +-5 % without it
    lo, hi = est.R_inf_range
    assert lo == pytest.approx(0.95 * est.R_inf) and hi <= 1.05 * est.R_inf
    assert lo <= Rs <= hi


def test_top_decade_replaces_the_bound_when_two_decades_do_not_determine_r_s():
    """Code review: rc/143 with a 100 Ohm second arc spoils the two-decade
    fit (R_s -17 % at 11.7 %), which fell back to the HF bound while the
    top decade gives R_s within 0.3 %. Now that is used; the range stays the
    bound's, nothing measures its model error, and the warning stays."""
    f, Z, Rs = _rc143(R2=100.0, C2=0.002)
    est = estimate_rinf(f, Z)
    assert est.method == 'rlq_fit' and est.fit_window == 'top_only'
    assert abs(est.R_inf - Rs) < 0.005 * Rs
    assert est.R_inf_range == (0.0, max(est.R_inf_upper, est.R_inf))
    assert 'not determined' in est.warnings[-1] and 'top 1 decade gives' in est.warnings[-1]
