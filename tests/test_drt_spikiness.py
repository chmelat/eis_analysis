#!/usr/bin/env python3
"""
Tests for DRT shape-quality detection (audit finding F3).

Auto-lambda optimizes data fit, so for low-noise data it drives lambda toward 0,
producing a sparse/spiky DRT on which peak-shape analysis (scipy find_peaks, GMM)
is meaningless. The fix is advisory: a participation-ratio metric N_eff flags a
too-sparse DRT, and lambda landing at the search-range edge is surfaced. Neither
alters gamma or the detected peaks.

These tests pin:
1. _effective_bins metric correctness on known gamma.
2. Degenerate auto-lambda DRT emits the sparse/spiky warning (low N_eff).
3. Healthy DRT (manual lambda) does not warn (high N_eff).
4. lambda-edge detection on a clean low-noise RC fitted with auto-lambda.
"""

import numpy as np

from eis_analysis.drt.core import _effective_bins, calculate_drt
from eis_analysis.io import generate_synthetic_data
from eis_analysis.fitting.config import DRT_MIN_EFFECTIVE_BINS


def test_effective_bins_metric():
    """N_eff ~1 for a single spike, ~k for k equal bins, 0 for empty."""
    spike = np.zeros(50)
    spike[10] = 5.0
    assert abs(_effective_bins(spike) - 1.0) < 1e-9

    flat = np.zeros(50)
    flat[:8] = 3.0  # 8 equal nonzero bins -> N_eff == 8
    assert abs(_effective_bins(flat) - 8.0) < 1e-9

    assert _effective_bins(np.zeros(50)) == 0.0


def _two_cpe_data(seed=0):
    np.random.seed(seed)
    return generate_synthetic_data(noise=0.01)


def test_degenerate_autolambda_warns():
    """Auto-lambda on low-noise CPE data -> sparse DRT -> warning + low N_eff."""
    f, Z = _two_cpe_data(seed=0)
    r = calculate_drt(f, Z, peak_method='scipy', auto_lambda=True)

    assert r.diagnostics.n_effective_bins < DRT_MIN_EFFECTIVE_BINS
    assert any("sparse/spiky" in w for w in r.warnings), r.warnings


def test_healthy_lambda_no_warn():
    """Manual moderate lambda -> smooth DRT -> no sparse warning, high N_eff."""
    f, Z = _two_cpe_data(seed=0)
    r = calculate_drt(f, Z, peak_method='scipy', auto_lambda=False)  # lambda=0.1

    assert r.diagnostics.n_effective_bins > DRT_MIN_EFFECTIVE_BINS
    assert not any("sparse/spiky" in w for w in r.warnings), r.warnings
    assert not any("search-range edge" in w for w in r.warnings), r.warnings


def _clean_rc(noise=0.005, seed=0):
    """Single clean R||C (sharp DRT) -> auto-lambda collapses toward 0."""
    f = np.logspace(-2, 5, 70)
    w = 2 * np.pi * f
    Z = 10 + 100.0 / (1 + 1j * w * 0.1)
    np.random.seed(seed)
    Z = Z + noise * np.abs(Z) * (np.random.randn(70) + 1j * np.random.randn(70))
    return f, Z


def test_lambda_at_edge_detected():
    """Clean low-noise RC with auto-lambda: lambda hits range edge -> warning."""
    f, Z = _clean_rc()
    r = calculate_drt(f, Z, peak_method='scipy', auto_lambda=True)

    lam = r.diagnostics.lambda_sel
    assert lam.lambda_at_edge or lam.corner_at_edge, (
        f"expected edge flag; lambda={lam.lambda_value:.2e}"
    )
    assert any("search-range edge" in w for w in r.warnings), r.warnings
