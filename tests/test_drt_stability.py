#!/usr/bin/env python3
"""
Tests for lambda-probe peak stability diagnostics (drt/stability.py).

A peak of the main DRT solution is tracked across re-solves at lambdas around
the selected lambda*; a peak that survives every probe with small tau drift is
'stable', one missing from more than half of the probes is an 'artifact'.
"""

import numpy as np
import pytest

from eis_analysis.drt import calculate_drt, LambdaProbePoint, NNLSSolution
from eis_analysis.drt import stability as stability_mod
from eis_analysis.drt.stability import (
    _assess_peaks,
    _match_tolerances,
    MATCH_TOLERANCE_CAP,
    MATCH_TOLERANCE_FLOOR,
)


def _voigt_impedance(frequencies, R_inf, elements):
    """Ideal Voigt: R_inf + sum_i R_i / (1 + j*omega*tau_i)."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


FREQUENCIES = np.logspace(5, -2, 71)  # 100 kHz .. 10 mHz


# =============================================================================
# End-to-end: calculate_drt with lambda_probe
# =============================================================================

def test_well_separated_peaks_are_stable():
    """Two RC peaks 2 decades apart must be stable across the lambda probe."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy', lambda_probe=True)
    stability = r.diagnostics.stability

    assert stability is not None
    assert len(stability.probe_points) == 4
    assert all(p.success for p in stability.probe_points)
    assert len(stability.peak_stability) == 2
    for peak in stability.peak_stability:
        assert peak.verdict == 'stable', (
            f"peak tau={peak.tau_ref:.2e}: {peak.verdict}, "
            f"persistence {peak.persistence}/{peak.n_probes}, "
            f"drift {peak.tau_drift_decades:.2f}"
        )
        assert peak.persistence == peak.n_probes
        assert peak.tau_drift_decades < 0.15


def test_no_probe_by_default():
    """Backward compatibility: without lambda_probe there is no stability data."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3)])
    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')
    assert r.diagnostics.stability is None


def test_pure_resistor_no_peaks():
    """gamma ~ 0 (pure resistor): probe runs, empty peak stability, no crash."""
    Z = np.full_like(FREQUENCIES, 500.0, dtype=complex)
    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy', lambda_probe=True)
    stability = r.diagnostics.stability
    assert stability is not None
    assert stability.peak_stability == []


def test_probe_survives_solver_failure(monkeypatch):
    """One failing probe lambda is recorded, the others are still evaluated."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    real_solve = stability_mod._solve_nnls
    fail_lambda = {'value': None}

    def failing_solve(A, b, L, lambda_reg, n_tau, Z_arg):
        if fail_lambda['value'] is None:
            fail_lambda['value'] = lambda_reg  # fail the first probe only
        if lambda_reg == fail_lambda['value']:
            return NNLSSolution(gamma=None, success=False,
                                warnings=['forced failure'])
        return real_solve(A, b, L, lambda_reg, n_tau, Z_arg)

    monkeypatch.setattr(stability_mod, '_solve_nnls', failing_solve)

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy', lambda_probe=True)
    stability = r.diagnostics.stability

    failed = [p for p in stability.probe_points if not p.success]
    assert len(failed) == 1
    assert 'forced failure' in failed[0].error
    assert any('failed' in w for w in stability.warnings)
    # Remaining probes still classify both peaks
    assert len(stability.peak_stability) == 2
    assert all(p.n_probes == 3 for p in stability.peak_stability)


def test_gmm_reference_peaks():
    """GMM main peaks are tracked against scipy probe peaks."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    r = calculate_drt(FREQUENCIES, Z, peak_method='gmm', lambda_probe=True)
    stability = r.diagnostics.stability
    assert stability is not None
    assert len(stability.peak_stability) == len(r.peaks)
    assert all(p.verdict == 'stable' for p in stability.peak_stability)


# =============================================================================
# Unit: matching tolerance
# =============================================================================

def test_single_peak_tolerance_is_cap():
    tol = _match_tolerances(np.log10(np.array([1e-3])))
    assert tol[0] == MATCH_TOLERANCE_CAP


def test_close_peaks_tolerance_floor():
    # 0.3 decade apart -> half-gap 0.15 clamps to floor
    tol = _match_tolerances(np.array([-3.0, -2.7]))
    assert np.allclose(tol, MATCH_TOLERANCE_FLOOR)


def test_distant_peaks_tolerance_cap():
    # 3 decades apart -> half-gap 1.5 clamps to cap
    tol = _match_tolerances(np.array([-4.0, -1.0]))
    assert np.allclose(tol, MATCH_TOLERANCE_CAP)


# =============================================================================
# Unit: verdict assessment on synthetic probe points
# =============================================================================

def _probe_point(lam, peaks):
    return LambdaProbePoint(lambda_value=lam, success=True, peaks=peaks)


def test_artifact_verdict_low_persistence():
    """Peak present in only 1 of 4 probes -> artifact."""
    reference = [(1e-3, 100.0)]
    points = [
        _probe_point(1e-2, [{'tau': 1e-3, 'R_estimate': 100.0}]),
        _probe_point(3e-2, []),
        _probe_point(3e-1, []),
        _probe_point(1.0, []),
    ]
    result = _assess_peaks(reference, points)
    assert result[0].verdict == 'artifact'
    assert result[0].persistence == 1


def test_marginal_verdict_large_drift():
    """Peak present everywhere but drifting beyond the stable limit -> marginal."""
    reference = [(1e-3, 100.0)]
    # drift 0.3 decade < tolerance (0.5) so it matches, but > stable limit (0.2)
    tau_drifted = 10 ** (-3 + 0.3)
    points = [
        _probe_point(lam, [{'tau': tau_drifted, 'R_estimate': 100.0}])
        for lam in (1e-2, 3e-2, 3e-1, 1.0)
    ]
    result = _assess_peaks(reference, points)
    assert result[0].verdict == 'marginal'
    assert result[0].persistence == 4
    assert result[0].tau_drift_decades == pytest.approx(0.3, abs=0.01)


def test_drift_beyond_tolerance_not_matched():
    """A probe peak a full decade away must not match the reference."""
    reference = [(1e-3, 100.0)]
    points = [
        _probe_point(lam, [{'tau': 1e-2, 'R_estimate': 100.0}])
        for lam in (1e-2, 3e-2, 3e-1, 1.0)
    ]
    result = _assess_peaks(reference, points)
    assert result[0].persistence == 0
    assert result[0].verdict == 'artifact'


def test_probe_peak_matched_to_nearest_reference_only():
    """One probe peak cannot satisfy two reference peaks at once."""
    reference = [(1e-3, 100.0), (10 ** -2.4, 100.0)]  # 0.6 decade apart
    points = [
        _probe_point(lam, [{'tau': 1e-3, 'R_estimate': 100.0}])
        for lam in (1e-2, 3e-2, 3e-1, 1.0)
    ]
    result = _assess_peaks(reference, points)
    assert result[0].persistence == 4   # exact match
    assert result[1].persistence == 0   # 0.6 dec away, tolerance floor 0.3
