#!/usr/bin/env python3
"""
DRT peak significance (drt.significance): a local maximum of gamma is a peak
when refitting it as a shoulder of the nearest taller peak costs at least
DRT_PEAK_DCHI2_MIN noise variances.

The two cases the former 3 % height threshold got wrong: a resolved process
three orders smaller in R (dropped), and the lobes regularization splits one
broad process into (kept).
"""

import numpy as np

from eis_analysis.drt import calculate_drt
from eis_analysis.fitting import analyze_voigt_elements

FREQUENCIES = np.logspace(5, -2, 71)
OMEGA = 2 * np.pi * FREQUENCIES


def _zarc(R, tau, n=1.0):
    return R / (1 + (1j * OMEGA * tau) ** n)


def _noisy(Z, sigma, seed=0):
    e = np.random.default_rng(seed).standard_normal((2, len(Z)))
    return Z + sigma * np.abs(Z) * (e[0] + 1j * e[1])


def test_small_resolved_process_is_a_peak_and_an_arc():
    """100 Ohm next to 100 kOhm, 0.1 % noise: the fast arc (0.15 % of the
    gamma maximum) is kept, and the suggestion has both arcs."""
    Z = _noisy(10 + _zarc(100, 1e-5) + _zarc(1e5, 1.0), 1e-3)
    result = calculate_drt(FREQUENCIES, Z, auto_lambda=True)
    peaks = result.diagnostics.scipy_peaks
    taus = sorted(p['tau'] for p in peaks)

    assert len(taus) == 2
    assert abs(np.log10(taus[0] / 1e-5)) < 0.1
    assert abs(np.log10(taus[1] / 1.0)) < 0.1

    suggestion = analyze_voigt_elements(result.tau, result.gamma, FREQUENCIES, Z,
                                        peak_indices=[p['index'] for p in peaks])
    assert len(suggestion.elements) == 2


def test_broad_process_lobes_are_not_peaks():
    """One ZARC, n = 0.8, 1 % noise: the side lobes of its DRT are rejected
    (the 3 % threshold reported two peaks here)."""
    Z = _noisy(10 + _zarc(1000, 1e-3, 0.8), 1e-2)
    result = calculate_drt(FREQUENCIES, Z, auto_lambda=True)

    assert len(result.diagnostics.scipy_peaks) == 1
    assert result.diagnostics.rejected_peaks
