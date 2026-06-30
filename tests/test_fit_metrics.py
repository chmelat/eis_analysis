#!/usr/bin/env python3
"""Tests for compute_fit_metrics (fitting/diagnostics.py).

Focus: the reported relative error is weighting-consistent
(sum(w*|dZ|) / sum(w*|Z|)) and does not double-count 1/|Z| (audit finding #6).
"""

import numpy as np
import pytest

from eis_analysis.fitting.diagnostics import compute_fit_metrics

# Varying magnitude so weighted and unweighted metrics differ.
# |dZ| = [0.1, 0.5, 10], |Z| = [1, 10, 100], rel = [0.10, 0.05, 0.10].
Z = np.array([1.0, 10.0, 100.0], dtype=complex)
Z_FIT = np.array([1.1, 10.5, 110.0], dtype=complex)


def test_modulus_equals_mean_relative_error():
    # Modulus weighting (w = 1/|Z|) must reduce to mean(|dZ|/|Z|), not the
    # double-counted 1/|Z|^2 value (which was ~9.55%).
    fit_error_rel, _, _ = compute_fit_metrics(Z, Z_FIT, 'modulus')
    expected = np.mean(np.abs(Z - Z_FIT) / np.abs(Z)) * 100  # 8.333...
    assert fit_error_rel == pytest.approx(expected, rel=1e-9)


def test_uniform_equals_aggregate_relative_error():
    # Uniform weighting (w = 1) -> sum(|dZ|) / sum(|Z|).
    fit_error_rel, _, _ = compute_fit_metrics(Z, Z_FIT, 'uniform')
    expected = np.sum(np.abs(Z - Z_FIT)) / np.sum(np.abs(Z)) * 100  # 9.549...
    assert fit_error_rel == pytest.approx(expected, rel=1e-9)


def test_perfect_fit_zero_error():
    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z, 'modulus')
    assert fit_error_rel == pytest.approx(0.0, abs=1e-12)
    assert fit_error_abs == pytest.approx(0.0, abs=1e-12)
    assert quality == 'excellent'
