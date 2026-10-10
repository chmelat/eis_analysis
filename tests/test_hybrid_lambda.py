#!/usr/bin/env python3
"""Test hybrid lambda selection (GCV + L-curve) for DRT analysis."""

import os

import numpy as np
import pytest
from eis_analysis.drt.gcv import (
    compute_gcv_score,
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
)
from eis_analysis.drt.linear_system import _build_drt_matrices
from eis_analysis.fitting.config import DRT_LAMBDA_DEFAULT, DRT_LAMBDA_RANGE

EXAMPLE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "example"
)


# =============================================================================
# Helper functions
# =============================================================================

def generate_voigt_impedance(frequencies, R_inf, voigt_elements):
    """Generate impedance for Voigt circuit."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in voigt_elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


def _matrices(frequencies, Z, R_inf):
    """A, b, L from the production matrix builder (no test re-implementation)."""
    m = _build_drt_matrices(frequencies, Z, R_inf)
    return m.A, m.b, m.L


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def frequencies():
    """Test frequencies: 100kHz to 10mHz."""
    return np.logspace(5, -2, 71)


@pytest.fixture
def voigt_data(frequencies):
    """Generate Voigt circuit data."""
    R_inf = 100
    voigt_elements = [(1000, 1e-3), (2000, 1e-1)]
    Z = generate_voigt_impedance(frequencies, R_inf, voigt_elements)
    return frequencies, Z, R_inf


# =============================================================================
# Tests
# =============================================================================

def test_gcv_and_hybrid_similar_for_clean_voigt(voigt_data):
    """Test that GCV and Hybrid give similar results for clean Voigt data."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    lambda_hybrid, _, _ = find_optimal_lambda_hybrid(A, b, L)

    # For clean data, both methods should be within 1 order of magnitude
    ratio = lambda_hybrid / lambda_gcv
    assert 0.1 < ratio < 10, f"Methods differ too much: ratio={ratio:.2f}"


# =============================================================================
# compute_gcv_score correctness (F13)
# =============================================================================

def test_gcv_score_finite_positive(voigt_data):
    """GCV score is finite and positive across the search range."""
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    for lam in np.logspace(np.log10(DRT_LAMBDA_RANGE[0]), np.log10(DRT_LAMBDA_RANGE[1]), 12):
        score = compute_gcv_score(lam, A, b, L)
        assert np.isfinite(score), f"GCV score not finite at lambda={lam:.1e}"
        assert score > 0, f"GCV score not positive at lambda={lam:.1e}"


def test_gcv_score_has_minimum(voigt_data):
    """The lambda selected by GCV scores no worse than the range endpoints.

    Verifies find_optimal_lambda_gcv genuinely minimizes compute_gcv_score
    (the selector and the score function are consistent), not just returns
    something in range.
    """
    frequencies, Z, R_inf = voigt_data
    A, b, L = _matrices(frequencies, Z, R_inf)

    lambda_gcv, _ = find_optimal_lambda_gcv(A, b, L)
    score_opt = compute_gcv_score(lambda_gcv, A, b, L)
    score_lo = compute_gcv_score(DRT_LAMBDA_RANGE[0], A, b, L)
    score_hi = compute_gcv_score(DRT_LAMBDA_RANGE[1], A, b, L)

    assert score_opt <= score_lo + 1e-12, "selected lambda worse than lower edge"
    assert score_opt <= score_hi + 1e-12, "selected lambda worse than upper edge"


# =============================================================================
# Decision rule: the larger of the two lambdas, disagreement reported
# =============================================================================

def _two_zarc(noise, seed=0):
    """Two ZARCs over R_inf = 10 Ohm with proportional complex noise."""
    f = np.logspace(5, -2, 71)
    w = 2 * np.pi * f
    Z = (10 + 100 / (1 + (1j * w * 1e-4) ** 0.85)
         + 300 / (1 + (1j * w * 1e-1) ** 0.9))
    rng = np.random.default_rng(seed)
    return f, Z + noise * np.abs(Z) * (rng.standard_normal(71) + 1j * rng.standard_normal(71))


@pytest.mark.parametrize("spectrum, winner", [
    ('EISPOT-M136113-4', 'gcv'),   # corner 0.24 decades below GCV
    ('two_zarc_0.1%', 'lcurve'),   # corner 1.2 decades above GCV
])
def test_larger_lambda_wins(spectrum, winner):
    """The hybrid takes whichever of GCV and the L-curve corner is larger.

    One spectrum per direction, so min() in place of max(), an always-L-curve
    rule or swapped stage labels each fail one case.
    """
    from eis_analysis.drt import calculate_drt
    from eis_analysis.io import load_data

    if spectrum == 'two_zarc_0.1%':
        f, Z = _two_zarc(0.001)
        r = calculate_drt(f, Z, auto_lambda=True, r_inf_preset=10.0, inductance=False)
    else:
        path = os.path.join(EXAMPLE_DIR, f'{spectrum}.DTA')
        if not os.path.exists(path):
            pytest.skip(f"example/{spectrum}.DTA missing")
        data = load_data(path)
        r = calculate_drt(data.frequencies, data.Z, auto_lambda=True)

    lam = r.diagnostics.lambda_sel
    assert lam.hybrid_stage == winner
    assert lam.lambda_value == (lam.lambda_gcv if winner == 'gcv' else lam.lambda_lcurve)


@pytest.mark.parametrize("corner_below_gcv", [True, False])
def test_corner_below_gcv_is_warned(monkeypatch, corner_below_gcv):
    """The search's corner_below_gcv flag reaches DRTResult.warnings.

    The search is stubbed: a corner below GCV was not reached on any measured
    spectrum.
    """
    import eis_analysis.drt.linear_system as ls
    from eis_analysis.drt import calculate_drt

    diag = {'lambda_gcv': 1e-5, 'lambda_lcurve': 1e-7, 'method_used': 'gcv',
            'corner_at_edge': False, 'corner_below_gcv': corner_below_gcv}
    monkeypatch.setattr(ls, 'find_optimal_lambda_hybrid',
                        lambda A, b, L, **kw: (1e-5, 1.0, diag))

    f, Z = _two_zarc(0.005)
    r = calculate_drt(f, Z, auto_lambda=True, r_inf_preset=10.0, inductance=False)

    assert any("L-curve corner" in w for w in r.warnings) == corner_below_gcv


@pytest.mark.parametrize("kwargs, method, lambda_value", [
    ({}, 'hybrid', None),                                   # documented default
    ({'lambda_reg': 1e-4}, 'user', 1e-4),                   # explicit lambda wins...
    ({'lambda_reg': 1e-4, 'auto_lambda': True}, 'user', 1e-4),  # ...even over auto
    ({'auto_lambda': False}, 'default', DRT_LAMBDA_DEFAULT),  # opt-out
])
def test_lambda_selection_defaults(voigt_data, kwargs, method, lambda_value):
    """calculate_drt auto-selects lambda unless the caller gives one."""
    from eis_analysis.drt import calculate_drt

    f, Z, _ = voigt_data
    lam = calculate_drt(f, Z, **kwargs).diagnostics.lambda_sel

    assert lam.method == method
    if lambda_value is not None:
        assert lam.lambda_value == lambda_value


@pytest.mark.parametrize("bad", [0.0, -1e-6, float('nan'), float('inf')])
def test_invalid_lambda_is_rejected(voigt_data, bad):
    """An explicit lambda is used as is, so it must be usable: sqrt(lambda)
    of a negative or non-finite value would poison the system silently, and
    0 drops the regularization."""
    from eis_analysis.drt import calculate_drt

    f, Z, _ = voigt_data
    with pytest.raises(ValueError, match="lambda_reg"):
        calculate_drt(f, Z, lambda_reg=bad)
