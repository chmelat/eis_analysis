#!/usr/bin/env python3
"""
Data weighting of the DRT least-squares term.

The unweighted DRT fits residuals in absolute Ohm, so the low-frequency points
with the largest |Z| dominate and a small arc next to a large one is smoothed
away. Weighting rescales each frequency's rows of A and b by compute_weights,
renormalized so that ||w*Z|| = ||Z|| (lambda keeps its meaning). 'sqrt' is the
default: it resolves such arcs and, unlike 'modulus', does not blow up under
constant (|Z|-independent) noise.
"""

import itertools

import numpy as np
import pytest

from eis_analysis.drt import calculate_drt
from eis_analysis.drt.linear_system import _build_drt_matrices


FREQUENCIES = np.logspace(6, -2, 81)  # 1 MHz .. 10 mHz
R_INF = 20.0
# Three RC arcs, the fast one 40x smaller than the slow one.
ELEMENTS = [(50.0, 1e-5), (200.0, 1e-3), (2000.0, 1e-1)]  # (R, tau)
R_POL = sum(R for R, _ in ELEMENTS)


def _voigt_impedance(frequencies, R_inf, elements):
    """Ideal Voigt: R_inf + sum_i R_i / (1 + j*omega*tau_i)."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


Z_CLEAN = _voigt_impedance(FREQUENCIES, R_INF, ELEMENTS)


def _constant_noise(Z, seed=0):
    """Noise of 0.5 % of max|Z| on every point, independent of |Z|."""
    rng = np.random.default_rng(seed)
    sigma = 0.005 * np.max(np.abs(Z))
    return Z + sigma * (rng.standard_normal(len(Z)) + 1j * rng.standard_normal(len(Z)))


def _matched_peaks(result):
    """For each true element, (log10 tau error, relative R error) of the nearest scipy peak."""
    peaks = result.diagnostics.scipy_peaks or []
    matched = []
    for R, tau in ELEMENTS:
        best = min(peaks, key=lambda p: abs(np.log10(p['tau'] / tau)), default=None)
        if best is None or abs(np.log10(best['tau'] / tau)) > 0.5:
            matched.append(None)
        else:
            matched.append((abs(np.log10(best['tau'] / tau)),
                            abs(best['R_estimate'] - R) / R))
    return matched


# =============================================================================
# Matrix construction
# =============================================================================

def test_rows_are_weighted_and_norm_matched():
    """A and b rows are w * unweighted rows; ||w*Z|| == ||Z||; A_re/A_im untouched."""
    for weighting in ['sqrt', 'modulus', 'proportional']:
        plain = _build_drt_matrices(FREQUENCIES, Z_CLEAN, R_INF, 100)
        weighted = _build_drt_matrices(FREQUENCIES, Z_CLEAN, R_INF, 100, weighting)
        w = weighted.weights
        w2 = np.concatenate([w, w])

        assert np.allclose(weighted.A, w2[:, None] * plain.A, rtol=1e-12), weighting
        assert np.allclose(weighted.b, w2 * plain.b, rtol=1e-12)
        assert np.isclose(np.linalg.norm(w * Z_CLEAN), np.linalg.norm(Z_CLEAN), rtol=1e-12)
        # Reconstruction matrices stay physical
        assert np.array_equal(weighted.A_re, plain.A_re)
        assert np.array_equal(weighted.A_im, plain.A_im)
        # Weight shape: larger weight where |Z| is smaller
        assert w[np.argmin(np.abs(Z_CLEAN))] > w[np.argmax(np.abs(Z_CLEAN))]


def test_uniform_is_the_unweighted_system():
    """'uniform' reproduces the pre-weighting system exactly (backward compatibility)."""
    m = _build_drt_matrices(FREQUENCIES, Z_CLEAN, R_INF, 100, 'uniform')
    assert np.array_equal(m.weights, np.ones(len(FREQUENCIES)))
    assert np.array_equal(m.b, np.concatenate([Z_CLEAN.real - R_INF, Z_CLEAN.imag]))
    assert np.array_equal(m.A, np.vstack([m.A_re, m.A_im]))


# =============================================================================
# End-to-end effect
# =============================================================================

def test_sqrt_resolves_small_fast_arc_that_uniform_misses():
    """At the default lambda=0.1, uniform merges the 50 Ohm arc; sqrt recovers all three.

    Measured (noise-free): uniform finds 2 peaks, sqrt finds 3 with tau errors
    <= 0.02 decade and R errors 1.6 / 3.4 / 1.2 %.
    """
    uniform = calculate_drt(FREQUENCIES, Z_CLEAN, weighting='uniform')
    assert None in _matched_peaks(uniform), "uniform unexpectedly resolved the fast arc"

    sqrt = calculate_drt(FREQUENCIES, Z_CLEAN)  # default weighting
    assert sqrt.diagnostics.weighting == 'sqrt'
    for match, (R, tau) in zip(_matched_peaks(sqrt), ELEMENTS):
        assert match is not None, f"sqrt missed the peak at tau={tau}"
        tau_err, R_err = match
        assert tau_err < 0.15, f"tau={tau}: {tau_err:.3f} decade off"
        assert R_err < 0.10, f"tau={tau}: R off by {R_err:.1%}"
    assert abs(sqrt.R_pol - R_POL) / R_POL < 0.03


def test_sqrt_is_robust_to_constant_noise():
    """Constant noise: modulus weights amplify the HF noise, auto-lambda runs to the edge.

    Measured (seed 0): modulus lambda=4.0 (at edge), R_pol off 2.6 %;
    sqrt lambda=0.12 (inside the range), R_pol off 0.9 %.
    """
    Z = _constant_noise(Z_CLEAN)

    modulus = calculate_drt(FREQUENCIES, Z, auto_lambda=True, weighting='modulus')
    assert modulus.diagnostics.lambda_sel.lambda_at_edge

    sqrt = calculate_drt(FREQUENCIES, Z, auto_lambda=True, weighting='sqrt')
    assert not sqrt.diagnostics.lambda_sel.lambda_at_edge
    assert abs(sqrt.R_pol - R_POL) / R_POL < 0.02


def test_weighted_drt_is_scale_invariant():
    """Z * 1000 gives the same lambda and gamma * 1000 (norm-matched weights)."""
    Z = _constant_noise(Z_CLEAN)
    r1 = calculate_drt(FREQUENCIES, Z, auto_lambda=True)
    r2 = calculate_drt(FREQUENCIES, Z * 1000, auto_lambda=True)

    assert r1.lambda_used == pytest.approx(r2.lambda_used, rel=1e-9)
    assert np.allclose(r2.gamma, 1000 * r1.gamma, rtol=1e-6, atol=1e-9 * np.max(r2.gamma))


def test_non_finite_impedance_fails_gracefully():
    """One bad point must yield an empty result, not an exception.

    Regression: under a non-uniform weighting one NaN made every weight NaN
    (compute_weights normalizes by the mean) and np.linalg.cond(A) raised an
    uncaught LinAlgError.
    """
    for weighting, bad in itertools.product(['uniform', 'sqrt', 'modulus'], [np.nan, np.inf]):
        Z = Z_CLEAN.copy()
        Z[5] = bad
        result = calculate_drt(FREQUENCIES, Z, weighting=weighting)
        assert not result.success, (weighting, bad)


def test_unknown_weighting_is_rejected():
    """A misspelled weighting must not silently fall back to uniform."""
    with pytest.raises(ValueError, match="modulos"):
        calculate_drt(FREQUENCIES, Z_CLEAN, weighting='modulos')
