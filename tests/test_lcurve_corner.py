#!/usr/bin/env python3
"""
Tests for L-curve corner detection (audit finding F6).

`find_lcurve_corner` picks the point of maximum *signed* curvature. The sign
of the curvature depends on the orientation of the (rho, eta) curve and the
traversal direction (increasing lambda), so a wrong convention would make
argmax select a flat tail instead of the corner. These tests pin down:

- the curvature formula + sign convention (circle of known radius/orientation),
- that the corner is localized at the actual bend of a synthetic L-curve
  oriented exactly like a real one (increasing lambda -> rho up, eta down),
- that on a real DRT L-curve the orientation assumption holds and the corner
  lands in the interior,
- the edge-corner diagnostic emitted by find_optimal_lambda_hybrid.

A sign/orientation regression (e.g. argmax -> argmin) makes
test_corner_on_synthetic_L fail.
"""

import numpy as np

from eis_analysis.drt.gcv import (
    compute_lcurve_curvature,
    compute_lcurve_point,
    find_lcurve_corner,
    find_optimal_lambda_hybrid,
)
from eis_analysis.drt.core import _build_drt_matrices


# =============================================================================
# A. compute_lcurve_curvature - formula and sign convention
# =============================================================================

def test_curvature_circle_magnitude_and_sign():
    """CCW circle of radius r has constant signed curvature +1/r."""
    r = 2.0
    theta = np.linspace(0.0, np.pi / 2, 200)  # increasing -> counter-clockwise
    rho = r * np.cos(theta)
    eta = r * np.sin(theta)

    kappa = compute_lcurve_curvature(rho, eta)

    # Skip the ends where np.gradient uses one-sided differences.
    interior = kappa[5:-5]
    assert np.allclose(interior, 1.0 / r, rtol=0.05), (
        f"CCW circle curvature should be +{1.0/r:.3f}, got mean "
        f"{interior.mean():.3f}"
    )
    assert np.all(interior > 0), "CCW circle must have positive curvature"


def test_curvature_circle_cw_negative():
    """Reversing traversal (CW) flips the sign to -1/r."""
    r = 2.0
    theta = np.linspace(np.pi / 2, 0.0, 200)  # decreasing -> clockwise
    rho = r * np.cos(theta)
    eta = r * np.sin(theta)

    kappa = compute_lcurve_curvature(rho, eta)
    interior = kappa[5:-5]
    assert np.allclose(interior, -1.0 / r, rtol=0.05)
    assert np.all(interior < 0), "CW circle must have negative curvature"


def test_curvature_straight_line_zero():
    """Collinear points have ~zero curvature."""
    rho = np.linspace(0.0, 5.0, 50)
    eta = 2.0 * rho + 1.0
    kappa = compute_lcurve_curvature(rho, eta)
    assert np.allclose(kappa, 0.0, atol=1e-9)


def test_curvature_too_few_points():
    """n < 3 returns zeros of length n."""
    out = compute_lcurve_curvature(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    assert out.shape == (2,)
    assert np.all(out == 0.0)


# =============================================================================
# B. find_lcurve_corner - localization and correct extremum
# =============================================================================

def _synthetic_L(n=51, k=25, step=0.2):
    """
    Sharp L-curve oriented like a real one: increasing index (= increasing
    lambda) gives rho up / eta down, corner convex toward the origin at index k.

    Vertical branch (i < k): rho const, eta decreasing.
    Horizontal branch (i > k): eta const, rho increasing.
    """
    rho = np.zeros(n)
    eta = np.zeros(n)
    rho[k:] = np.arange(n - k) * step       # 0 on [0,k], rising after
    eta[:k + 1] = np.arange(k, -1, -1) * step  # falling to 0 at k, then 0
    lambda_values = np.logspace(-5, 1, n)
    return lambda_values, rho, eta, k


def test_corner_on_synthetic_L():
    """The corner is detected at the bend (KEY F6 regression test)."""
    lambda_values, rho, eta, k = _synthetic_L()
    _, corner_idx, _ = find_lcurve_corner(lambda_values, rho, eta)
    assert abs(corner_idx - k) <= 1, (
        f"corner_idx={corner_idx} not at the bend k={k} "
        f"(a sign/orientation flip would land elsewhere)"
    )


def test_corner_is_interior_not_endpoint():
    """The detected corner is strictly interior to the array."""
    lambda_values, rho, eta, _ = _synthetic_L()
    _, corner_idx, _ = find_lcurve_corner(lambda_values, rho, eta)
    assert 1 < corner_idx < len(rho) - 2


def test_corner_orientation_from_real_lcurve():
    """On a real DRT L-curve: rho increases / eta decreases with lambda, and
    the corner is interior."""
    f = np.logspace(5, -1, 60)
    w = 2 * np.pi * f
    R_inf, R, tau0 = 10.0, 100.0, 1e-3
    Z = R_inf + R / (1 + 1j * w * tau0)
    rng = np.random.default_rng(0)
    Z = Z + (0.01 * np.abs(Z)) * (rng.standard_normal(len(Z))
                                  + 1j * rng.standard_normal(len(Z)))

    m = _build_drt_matrices(f, Z, R_inf, n_tau=80)
    lambdas = np.logspace(-5, 0, 30)
    rho = np.empty(len(lambdas))
    eta = np.empty(len(lambdas))
    for i, lam in enumerate(lambdas):
        rho[i], eta[i], _ = compute_lcurve_point(lam, m.A, m.b, m.L)

    # Orientation assumption the corner-finder relies on.
    assert rho[-1] > rho[0], "residual norm should grow with lambda"
    assert eta[-1] < eta[0], "solution seminorm should shrink with lambda"

    _, corner_idx, _ = find_lcurve_corner(lambdas, rho, eta)
    assert 0 < corner_idx < len(lambdas) - 1


def test_corner_fallback_few_points():
    """n < 3 -> graceful middle fallback, no crash."""
    lambda_values = np.array([1e-3, 1e-1])
    rho = np.array([0.1, 0.5])
    eta = np.array([0.5, 0.1])
    lam_c, corner_idx, curv = find_lcurve_corner(lambda_values, rho, eta)
    assert 0 <= corner_idx < len(lambda_values)
    assert lam_c in lambda_values


# =============================================================================
# C. find_optimal_lambda_hybrid - edge-corner diagnostic
# =============================================================================

def _voigt_matrices(noise=0.01, n_tau=80, seed=1):
    f = np.logspace(5, -1, 60)
    w = 2 * np.pi * f
    R_inf, R, tau0 = 10.0, 100.0, 1e-3
    Z = R_inf + R / (1 + 1j * w * tau0)
    rng = np.random.default_rng(seed)
    Z = Z + (noise * np.abs(Z)) * (rng.standard_normal(len(Z))
                                   + 1j * rng.standard_normal(len(Z)))
    m = _build_drt_matrices(f, Z, R_inf, n_tau=n_tau)
    return m.A, m.b, m.L


def test_hybrid_corner_interior_clean():
    """Normal noisy Voigt: lambda_lcurve is interior, corner_at_edge False."""
    A, b, L = _voigt_matrices(noise=0.02)
    _, _, diag = find_optimal_lambda_hybrid(A, b, L)
    assert diag['corner_at_edge'] is False
    lo, hi = diag['lambda_values'][0], diag['lambda_values'][-1]
    assert lo < diag['lambda_lcurve'] < hi


def test_edge_corner_condition_detected():
    """
    When the bend sits at the boundary of the searched range, find_lcurve_corner
    returns a boundary index -- this is exactly the condition
    find_optimal_lambda_hybrid flags as corner_at_edge=True.

    Forcing the flag through the full hybrid is data-dependent (GCV centers the
    L-curve window on the corner, so it rarely reaches the boundary), so the
    detection condition is tested deterministically here; the False branch of
    the plumbing is exercised by test_hybrid_corner_interior_clean.
    """
    lambda_values, rho, eta, k = _synthetic_L(n=20, k=18)  # bend near top edge
    _, corner_idx, _ = find_lcurve_corner(lambda_values, rho, eta)

    n = len(lambda_values)
    assert corner_idx == k
    # The hybrid hook flags exactly this: corner_idx <= 1 or >= n - 2.
    assert corner_idx <= 1 or corner_idx >= n - 2
