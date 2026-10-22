#!/usr/bin/env python3
"""
End-to-end correctness tests for DRT reconstruction (audit finding F13).

The lambda-selection tests were smoke-only (lambda positive/finite). These pin
that the full `calculate_drt` pipeline recovers *known* spectra: for an ideal
Voigt (R||C) circuit the DRT must show peaks at the true time constants
tau = R*C, and the integral of gamma must recover R_pol = sum(R_i).

Tolerances are empirically grounded (two-peak Voigt, fixed DRT_LAMBDA_DEFAULT,
measured 2026-10-10): detected tau within ~0.005 decade of truth, R_pol within
~0.4%. We assert the safer 0.15
decade / 3% to stay robust across noise realizations.
"""

import numpy as np

from eis_analysis.drt import calculate_drt
from eis_analysis.fitting.config import DRT_MIN_EFFECTIVE_BINS


def _voigt_impedance(frequencies, R_inf, elements):
    """Ideal Voigt: R_inf + sum_i R_i / (1 + j*omega*tau_i). DRT peak at tau_i."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


FREQUENCIES = np.logspace(5, -2, 71)  # 100 kHz .. 10 mHz, covers tau 1e-3..1e-1


def _peak_taus(result):
    return sorted(p['tau'] for p in (result.diagnostics.scipy_peaks or []))


def test_two_peak_recovery():
    """Two well-separated RC peaks recovered at the correct tau, R_pol exact."""
    R_inf = 100.0
    elements = [(1000.0, 1e-3), (2000.0, 1e-1)]  # (R, tau)
    Z = _voigt_impedance(FREQUENCIES, R_inf, elements)

    # Fixed DRT_LAMBDA_DEFAULT: the smoothness check below is calibrated on it
    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy', auto_lambda=False)
    taus = _peak_taus(r)

    assert len(taus) == 2, f"expected 2 peaks, got {len(taus)}: {taus}"
    for got, (R, tau_true) in zip(taus, elements):
        assert abs(np.log10(got / tau_true)) < 0.15, (
            f"peak tau {got:.2e} off from true {tau_true:.2e}"
        )

    R_pol_true = sum(R for R, _ in elements)
    assert abs(r.R_pol - R_pol_true) / R_pol_true < 0.03, (
        f"R_pol {r.R_pol:.1f} != true {R_pol_true:.1f}"
    )
    # Smooth DRT at this lambda - shape analysis is meaningful (F3 metric).
    assert r.diagnostics.n_effective_bins > DRT_MIN_EFFECTIVE_BINS


def test_scipy_peak_frequency_convention():
    """Scipy peak 'frequency' is the characteristic frequency 1/(2*pi*tau).

    Regression for audit 2026-07-02 finding 2.1: scipy peaks reported
    f = 1/tau while GMM peaks and auto_suggest use f = 1/(2*pi*tau), so the
    CLI showed frequencies differing by 2*pi for the same process.
    """
    R_inf = 50.0
    R, tau_true = 500.0, 1e-2
    Z = _voigt_impedance(FREQUENCIES, R_inf, [(R, tau_true)])

    r = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='scipy')
    peaks = r.diagnostics.scipy_peaks
    assert peaks, "expected at least one scipy peak"

    # Structural identity: frequency key is exactly 1/(2*pi*tau) per peak.
    for p in peaks:
        assert p['frequency'] == 1 / (2 * np.pi * p['tau'])

    # Physical check: detected frequency matches the true characteristic
    # frequency within the same 0.15-decade tolerance as tau recovery.
    f_char_true = 1 / (2 * np.pi * tau_true)
    assert abs(np.log10(peaks[0]['frequency'] / f_char_true)) < 0.15


def test_peak_r_estimate_invariant_to_normalize_rpol():
    """Peak R_estimates stay in Ohm when gamma is normalized by R_pol.

    Regression for audit 2026-07-02 finding 2.2: peak detection ran on the
    normalized gamma, so with normalize_rpol=True the R_estimates were
    dimensionless fractions (sum = 1) instead of resistances in Ohm.
    """
    R_inf = 100.0
    elements = [(1000.0, 1e-3), (2000.0, 1e-1)]
    Z = _voigt_impedance(FREQUENCIES, R_inf, elements)

    r_plain = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='scipy')
    r_norm = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='scipy',
                           normalize_rpol=True)

    R_plain = [p['R_estimate'] for p in r_plain.diagnostics.scipy_peaks]
    R_norm = [p['R_estimate'] for p in r_norm.diagnostics.scipy_peaks]
    assert np.allclose(R_plain, R_norm), (
        f"R_estimate changed under normalize_rpol: {R_plain} vs {R_norm}"
    )
    # Physical sanity: recovered R close to the true element resistances.
    for got, R_true in zip(sorted(R_norm), sorted(R for R, _ in elements)):
        assert abs(got - R_true) / R_true < 0.05

    # The returned gamma is still the normalized one (integral = 1) and the
    # unnormalized original is preserved alongside it.
    d_ln_tau = float(np.mean(np.diff(np.log(r_norm.tau))))
    assert abs(np.sum(r_norm.gamma) * d_ln_tau - 1.0) < 1e-6
    assert r_norm.gamma_original is not None

    # Same invariance for GMM peaks (R_i = weight_i * R_pol must be in Ohm).
    g_plain = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='gmm')
    g_norm = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='gmm',
                           normalize_rpol=True)
    assert g_plain.peaks and g_norm.peaks
    R_g_plain = [p['R_estimate'] for p in g_plain.peaks]
    R_g_norm = [p['R_estimate'] for p in g_norm.peaks]
    assert np.allclose(R_g_plain, R_g_norm)


def test_single_peak_recovery():
    """Single RC: one peak at the right tau, R_pol recovered."""
    R_inf = 50.0
    R, tau_true = 500.0, 1e-2
    Z = _voigt_impedance(FREQUENCIES, R_inf, [(R, tau_true)])

    r = calculate_drt(FREQUENCIES, Z, auto_lambda=False, peak_method='scipy')
    taus = _peak_taus(r)

    assert len(taus) == 1, f"expected 1 peak, got {len(taus)}: {taus}"
    assert abs(np.log10(taus[0] / tau_true)) < 0.15
    assert abs(r.R_pol - R) / R < 0.03, f"R_pol {r.R_pol:.1f} != true {R:.1f}"


def test_default_r_inf_is_the_window_fit():
    """Without a preset the DRT takes R_inf from estimate_rinf: an arc still
    open at f_max, where the HF median (and with it R_pol) is off."""
    f = np.logspace(5, -1, 61)
    Z = 10 + 100 / (1 + (2j * np.pi * f * 1e-4) ** 0.8)
    r = calculate_drt(f, Z)
    rinf = r.diagnostics.rinf
    assert rinf.method == 'rlq_fit' and rinf.estimate is not None
    assert r.R_inf == rinf.estimate.R_inf
    assert abs(r.R_inf - 10) < 0.05 < abs(rinf.R_inf_median - 10)


def test_default_r_inf_warning_joins_the_drt_warnings():
    """A bound instead of the fit qualifies R_pol too."""
    from eis_analysis.io import load_csv_data

    # Two overlapping CPEs: the fit cannot determine R_s (test_rinf_estimation)
    data = load_csv_data('example/example_eis_data.csv')
    r = calculate_drt(data.frequencies, data.Z)
    est = r.diagnostics.rinf.estimate
    assert est.method == 'hf_bound' and est.warnings
    assert set(est.warnings) <= set(r.diagnostics.nnls.warnings)


def test_nan_frequency_is_rejected():
    """NaN passed the f <= 0 test and estimate_rinf raised on it."""
    f = np.full(61, np.nan)
    assert not calculate_drt(f, np.ones(61, dtype=complex)).success


def test_voigt_suggestion_uses_the_drt_r_inf():
    """Code review: the suggestion took R_inf from the HF median (11.8 Ohm
    here) while the DRT's peaks stood on the window fit's 10.0, so R_pol
    from the data and the peaks' resistances disagreed."""
    from eis_analysis.fitting import analyze_voigt_elements

    f = np.logspace(5, -1, 61)
    Z = 10 + 100 / (1 + (2j * np.pi * f * 1e-4) ** 0.8)
    r = calculate_drt(f, Z)
    s = analyze_voigt_elements(r.tau, r.gamma, f, Z, r.R_inf,
                               peak_indices=[p['index'] for p in r.diagnostics.scipy_peaks])
    assert s.R_inf == r.R_inf
    assert abs(s.R_pol - 100) < 0.5


def test_default_r_inf_leaks_no_runtime_warning():
    """Code review: the window fit, now in every default DRT, put L at
    4e-35 H with stderr 2e-7, and stderr / L overflowed in the parameter
    check of fit_equivalent_circuit (seed 27, 3 % noise)."""
    import warnings

    f = np.logspace(5, -1, 61)
    Zc = 10 + 100 / (1 + 2j * np.pi * f * 1e-2)
    rng = np.random.default_rng(27)
    Z = Zc + 0.03 * np.abs(Zc) * (rng.standard_normal(61) + 1j * rng.standard_normal(61)) / np.sqrt(2)
    # Recorded, not raised: an exception inside the fit would be caught there
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        calculate_drt(f, Z)
    assert not [w for w in caught if issubclass(w.category, RuntimeWarning)]
