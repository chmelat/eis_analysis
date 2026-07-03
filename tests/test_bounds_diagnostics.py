"""
Tests for at-bound parameter diagnostics (fitting/bounds.py + circuit.py).

Regression tests for audit findings (2026-07-03):
- D1: FitDiagnostics.bounds_warnings/params_at_bounds and
  FitResult.bound_status must use one criterion (classify_bound_status)
  and therefore cannot disagree.
- D2: warning indices/labels are in full parameter space (fixed params
  included), consistent with param labels shown to the user.
"""

import numpy as np

from eis_analysis.fitting import fit_equivalent_circuit
from eis_analysis.fitting.bounds import build_bound_status
from eis_analysis.fitting.circuit_elements import C, Q, R


# --- Unit tests: build_bound_status ---

R_BOUNDS_LO, R_BOUNDS_HI = 1e-4, 1e10  # PARAMETER_BOUNDS['R']


def test_build_bound_status_interior():
    status = build_bound_status(
        np.array([100.0]), [R_BOUNDS_LO], [R_BOUNDS_HI], None)
    assert status == ['']


def test_build_bound_status_near_lower_log():
    # 5e-4 is 0.7 decades above 1e-4 (< 1 decade threshold for wide bounds)
    status = build_bound_status(
        np.array([5e-4]), [R_BOUNDS_LO], [R_BOUNDS_HI], None)
    assert status == ['lower']


def test_build_bound_status_near_upper_log():
    status = build_bound_status(
        np.array([2e9]), [R_BOUNDS_LO], [R_BOUNDS_HI], None)
    assert status == ['upper']


def test_build_bound_status_linear_branch():
    # n bounds (0.3, 1.0) span < 6 decades -> linear 1%-of-range criterion
    status = build_bound_status(
        np.array([0.65, 0.995]), [0.3, 0.3], [1.0, 1.0], None)
    assert status == ['', 'upper']


def test_build_bound_status_fixed_param():
    status = build_bound_status(
        np.array([100.0, 5e-4]), [R_BOUNDS_LO] * 2, [R_BOUNDS_HI] * 2,
        [True, False])
    assert status == ['fixed', 'lower']


def test_build_bound_status_no_bounds():
    status = build_bound_status(np.array([1.0, 2.0]), None, None, None)
    assert status == ['', '']


# --- Regression D1: one criterion for warnings and bound_status ---

def _at_bound_indices(result):
    return [i for i, s in enumerate(result.bound_status)
            if s in ('lower', 'upper')]


def test_bounds_warnings_consistent_with_bound_status():
    """Fitted R below 1e-3 is 'lower' per bound_status -> warning must exist.

    The pre-fix Step 4 criterion (|p-b|/|b| < 0.01) would stay silent for
    R = 5e-4 while bound_status said 'lower' (contradictory CLI output).
    """
    freq = np.logspace(4, 0, 20)
    Z = np.full_like(freq, 5e-4, dtype=complex)  # pure resistor at 0.5 mOhm

    result, _, _ = fit_equivalent_circuit(freq, Z, R(1e-3), plot=False)

    assert abs(result.params_opt[0] - 5e-4) / 5e-4 < 1e-3
    assert result.bound_status == ['lower']
    assert result.diagnostics.params_at_bounds == _at_bound_indices(result)
    assert len(result.diagnostics.bounds_warnings) == 1
    assert 'lower' in result.diagnostics.bounds_warnings[0]


def test_interior_fit_no_bounds_warnings():
    """Well-conditioned Voigt fit: no parameter near a bound, no warnings."""
    circuit = R(100.0) - (R(5000.0) | C(1e-6))
    freq = np.logspace(5, -1, 40)
    Z = circuit.impedance(freq, [100.0, 5000.0, 1e-6])

    result, _, _ = fit_equivalent_circuit(freq, Z, circuit, plot=False)

    assert result.bound_status == ['', '', '']
    assert result.diagnostics.params_at_bounds == []
    assert result.diagnostics.bounds_warnings == []


# --- Regression D2: full-space indices and labels with fixed params ---

def test_bounds_warning_full_space_index_with_fixed_param():
    """n0 driven to its upper bound behind a fixed R0.

    Full-space parameter order is [R0 (fixed), R1, Q0, n0]; the warning must
    refer to index 3 / label 'n0'. The pre-fix code reported the free-space
    index 2, which the user would read as Q0.
    """
    freq = np.logspace(5, -1, 40)
    # Data from an ideal capacitor (n = 1) -> fitted n0 ends at upper bound 1.0
    omega = 2 * np.pi * freq
    Z = 100.0 + 5000.0 / (1 + 1j * omega * 5000.0 * 1e-6)

    circuit = R("100") - (R(5000.0) | Q(1e-6, 0.95))
    result, _, _ = fit_equivalent_circuit(freq, Z, circuit, plot=False)

    assert result.bound_status[0] == 'fixed'
    assert result.bound_status[3] == 'upper'
    assert 3 in result.diagnostics.params_at_bounds
    assert result.diagnostics.params_at_bounds == _at_bound_indices(result)
    n_warnings = [w for w in result.diagnostics.bounds_warnings if 'n0' in w]
    assert len(n_warnings) == 1
    assert 'upper' in n_warnings[0]
