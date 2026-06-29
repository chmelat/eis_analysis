#!/usr/bin/env python3
"""
End-to-end correctness tests for DRT reconstruction (audit finding F13).

The lambda-selection tests were smoke-only (lambda positive/finite). These pin
that the full `calculate_drt` pipeline recovers *known* spectra: for an ideal
Voigt (R||C) circuit the DRT must show peaks at the true time constants
tau = R*C, and the integral of gamma must recover R_pol = sum(R_i).

Tolerances are empirically grounded (two-peak Voigt, lambda=0.1): detected tau
within ~0.04 decade of truth, R_pol within ~0.4%. We assert the safer 0.15
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

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')  # default lambda=0.1
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


def test_single_peak_recovery():
    """Single RC: one peak at the right tau, R_pol recovered."""
    R_inf = 50.0
    R, tau_true = 500.0, 1e-2
    Z = _voigt_impedance(FREQUENCIES, R_inf, [(R, tau_true)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')
    taus = _peak_taus(r)

    assert len(taus) == 1, f"expected 1 peak, got {len(taus)}: {taus}"
    assert abs(np.log10(taus[0] / tau_true)) < 0.15
    assert abs(r.R_pol - R) / R < 0.03, f"R_pol {r.R_pol:.1f} != true {R:.1f}"
