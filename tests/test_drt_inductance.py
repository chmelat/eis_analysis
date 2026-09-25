#!/usr/bin/env python3
"""
Series inductance L in the DRT model.

Without it the DRT kernel yields only Im(Z) < 0, so an inductive
high-frequency end deforms gamma: with the exact R_inf, the reconstruction
error of Rs + jwL + RC was 5-18 % (doc/DRT_RINF_L_ANALYSIS_2026-09-25.md).
The L column is unregularized and solved with gamma; 'auto' adds it only when
the top decade has a point with Im(Z) > 0.
"""

import numpy as np
import pytest

from eis_analysis.drt import calculate_drt

# 0.1 Hz .. 1 MHz, 10 points/decade
FREQUENCIES = np.logspace(6, -1, 71)
OMEGA = 2 * np.pi * FREQUENCIES
R_S = 10.0


def _zarc(R, tau, n=1.0):
    return R / (1 + (1j * OMEGA * tau) ** n)


def _noisy(Z, seed=0, level=0.01):
    rng = np.random.default_rng(seed)
    return Z * (1 + level * (rng.standard_normal(len(Z)) + 1j * rng.standard_normal(len(Z))))


@pytest.mark.parametrize('n', [1.0, 0.7], ids=['C RC', 'C2 ZARC 0.7'])
def test_auto_models_inductive_end(n):
    # Audit cases C and C2: L = 10 uH, 1 % noise.
    Z = _noisy(R_S + 1j * OMEGA * 10e-6 + _zarc(100.0, 1e-5, n))
    r = calculate_drt(FREQUENCIES, Z, r_inf_preset=R_S, auto_lambda=True)
    assert r.diagnostics.inductance_used
    assert r.L_series == pytest.approx(10e-6, rel=0.05)
    assert r.reconstruction_error < 5.0  # 17 % / 15 % without L
    assert r.R_pol == pytest.approx(100.0, rel=0.03)


# Capacitive end (ZARC n = 0.8): no point with Im(Z) > 0 in the top decade.
Z_CAPACITIVE = _noisy(R_S + _zarc(100.0, 1e-4, 0.8))


def test_auto_leaves_capacitive_end_unchanged():
    auto = calculate_drt(FREQUENCIES, Z_CAPACITIVE, r_inf_preset=R_S, auto_lambda=True)
    off = calculate_drt(FREQUENCIES, Z_CAPACITIVE, r_inf_preset=R_S, auto_lambda=True,
                        inductance=False)
    assert not auto.diagnostics.inductance_used
    assert auto.L_series == 0.0
    np.testing.assert_array_equal(auto.gamma, off.gamma)


def test_forced_inductance_on_capacitive_end_warns():
    # The L it finds on such data is model error, not a measured inductance.
    forced = calculate_drt(FREQUENCIES, Z_CAPACITIVE, r_inf_preset=R_S, auto_lambda=True,
                           inductance=True)
    assert forced.L_series > 0
    assert any('without inductive points' in w for w in forced.warnings)


def test_rejects_unknown_inductance_mode():
    with pytest.raises(ValueError):
        calculate_drt(FREQUENCIES, Z_CAPACITIVE, inductance='yes')
