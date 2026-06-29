#!/usr/bin/env python3
"""
Regression tests for peak resistance estimation (audit finding F4).

Both the scipy path (`_estimate_peak_resistance`) and the GMM path
(`gmm_peak_detection`) used to integrate the *total* gamma over each peak's
window, double-counting the overlap region of adjacent peaks so that
sum(R_i) > R_pol. The fixes partition the resistance: scipy splits at the
valleys between peaks, GMM decomposes R_pol by component weight (sum w = 1).
These tests lock sum(R_i) == R_pol for overlapping peaks.
"""

import numpy as np
import pytest

from eis_analysis.drt.core import _estimate_peak_resistance, _rpol_from_gamma
from eis_analysis.drt.peaks import gmm_peak_detection


def _rpol(gamma, tau):
    """Rectangle-rule R_pol, matching the DRT kernel and production (F10)."""
    d_ln_tau = float(np.mean(np.diff(np.log(tau))))
    return _rpol_from_gamma(gamma, d_ln_tau)


def _two_gaussian_gamma(tau, centers, sigma, amps):
    """Two overlapping Gaussians in log10(tau) space."""
    log_tau = np.log10(tau)
    gamma = np.zeros_like(tau)
    for mu, a in zip(centers, amps):
        gamma += a * np.exp(-0.5 * ((log_tau - mu) / sigma) ** 2)
    return gamma


@pytest.fixture
def overlapping_peaks():
    """Two overlapping DRT peaks (1.2 decades apart, sigma 0.3)."""
    tau = np.logspace(-5, 1, 400)
    gamma = _two_gaussian_gamma(tau, centers=[-3.4, -2.2], sigma=0.3,
                                amps=[80.0, 50.0])
    return tau, gamma


def test_scipy_peak_resistance_sums_to_rpol(overlapping_peaks):
    """sum(R_i) over valley-partitioned peaks equals R_pol of spanned range."""
    tau, gamma = overlapping_peaks

    from scipy.signal import find_peaks
    peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * 0.05)
    assert len(peaks_idx) >= 2, "fixture should produce >=2 detectable peaks"

    resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)

    # Disjoint half-open segments cover the whole array, so the rectangle-rule
    # sum equals the full-range rectangle R_pol exactly (F10).
    R_pol_full = _rpol(gamma, tau)
    assert np.isclose(sum(resistances), R_pol_full, rtol=1e-9), (
        f"sum(R_i)={sum(resistances):.4f} != R_pol={R_pol_full:.4f}"
    )

    # Each peak resistance must be positive and physically below the total.
    assert all(r > 0 for r in resistances)
    assert all(r < R_pol_full for r in resistances)


def test_scipy_peak_resistance_no_double_count(overlapping_peaks):
    """The fixed sum must not exceed R_pol (the old code double-counted)."""
    tau, gamma = overlapping_peaks
    from scipy.signal import find_peaks
    peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * 0.05)

    resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)
    R_pol_full = _rpol(gamma, tau)

    assert sum(resistances) <= R_pol_full * (1 + 1e-9)


def test_scipy_peak_resistance_empty():
    """No peaks -> empty list."""
    tau = np.logspace(-5, 1, 100)
    assert _estimate_peak_resistance(tau, np.ones_like(tau), np.array([], dtype=int)) == []


def test_gmm_peak_resistance_sums_to_rpol(overlapping_peaks):
    """GMM R_estimate decomposition: sum(R_i) == R_pol exactly (sum w = 1)."""
    tau, gamma = overlapping_peaks
    peaks, gmm_model, _ = gmm_peak_detection(tau, gamma, bic_threshold=10.0)

    assert gmm_model is not None and len(peaks) >= 1
    R_pol_full = _rpol(gamma, tau)

    total = sum(p['R_estimate'] for p in peaks)
    # GMM weights sum to 1, so weight_i * R_pol sums to R_pol exactly.
    assert np.isclose(total, R_pol_full, rtol=1e-9), (
        f"sum(R_estimate)={total:.4f} != R_pol={R_pol_full:.4f}"
    )
    assert all(p['R_estimate'] > 0 for p in peaks)


def test_rpol_unified_integration(overlapping_peaks):
    """F10: scipy sum(R_i) and GMM R_pol both equal the rectangle R_pol.

    All three derive R_pol from gamma with the same rectangle rule (consistent
    with the DRT kernel), so they must agree exactly for the same gamma.
    """
    tau, gamma = overlapping_peaks
    rect = _rpol(gamma, tau)

    from scipy.signal import find_peaks
    idx, _ = find_peaks(gamma, height=np.max(gamma) * 0.05)
    scipy_sum = sum(_estimate_peak_resistance(tau, gamma, idx))

    peaks, model, _ = gmm_peak_detection(tau, gamma, bic_threshold=10.0)
    gmm_sum = sum(p['R_estimate'] for p in peaks)

    assert np.isclose(scipy_sum, rect, rtol=1e-9)
    assert np.isclose(gmm_sum, rect, rtol=1e-9)
