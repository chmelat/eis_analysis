"""
Tests for at-bound parameter diagnostics (fitting/bounds.py + circuit.py).

Regression tests for audit findings (2026-07-03):
- D1: FitDiagnostics.bounds_warnings/params_at_bounds and
  FitResult.bound_status must use one criterion (classify_bound_status)
  and therefore cannot disagree.
- D2: warning indices/labels are in full parameter space (fixed params
  included), consistent with param labels shown to the user.
- D3: initial guess silently clipped into bounds must be reported in
  FitDiagnostics.warnings.
"""

import numpy as np
import pytest

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

    result, _ = fit_equivalent_circuit(freq, Z, R(1e-3))

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

    result, _ = fit_equivalent_circuit(freq, Z, circuit)

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
    result, _ = fit_equivalent_circuit(freq, Z, circuit)

    assert result.bound_status[0] == 'fixed'
    assert result.bound_status[3] == 'upper'
    assert 3 in result.diagnostics.params_at_bounds
    assert result.diagnostics.params_at_bounds == _at_bound_indices(result)
    n_warnings = [w for w in result.diagnostics.bounds_warnings if 'n0' in w]
    assert len(n_warnings) == 1
    assert 'upper' in n_warnings[0]


# --- Regression D3: clipped initial guess is reported ---

def test_clipped_initial_guess_warns():
    """R(0) is below the lower bound 1e-4 -> clipped, warning must appear.

    Pre-fix, _prepare_optimization recorded clipped_params but nothing
    propagated it: the fit silently started from a different point than
    the user specified.
    """
    circuit = R(0.0) - (R(5000.0) | C(1e-6))
    freq = np.logspace(5, -1, 40)
    Z = circuit.impedance(freq, [100.0, 5000.0, 1e-6])

    result, _ = fit_equivalent_circuit(freq, Z, circuit)

    clip_warnings = [w for w in result.diagnostics.warnings if 'clipped' in w]
    assert len(clip_warnings) == 1
    assert 'R0' in clip_warnings[0]
    assert clip_warnings[0] in result.all_warnings


def test_in_bounds_guess_no_clip_warning():
    """Guess inside bounds: no clipping warning."""
    circuit = R(100.0) - (R(5000.0) | C(1e-6))
    freq = np.logspace(5, -1, 40)
    Z = circuit.impedance(freq, [100.0, 5000.0, 1e-6])

    result, _ = fit_equivalent_circuit(freq, Z, circuit)

    assert [w for w in result.diagnostics.warnings if 'clipped' in w] == []


# --- Units: a fit must not depend on them ---

def _blocking_spectrum():
    """Stress test case blocking/38: R - (R|Q) - Q, noise-free."""
    f = np.logspace(np.log10(107962.72875713777), np.log10(0.24267724373552552), 52)
    truth = [4254241.56646422, 232721572.02547705, 1.838929474718818e-10,
             0.6633267992074743, 4.0252334682656605e-11, 0.8996571087461844]
    circuit = R(truth[0]) - (R(truth[1]) | Q(truth[2], truth[3])) - Q(truth[4], truth[5])
    return f, circuit.impedance(f, truth)


def _start():
    return R(8.5e6) - (R(9.3e7) | Q(5.5e-10, 0.7)) - Q(2e-11, 0.85)


def test_fit_does_not_depend_on_units():
    """Exactly halved data (a binary scaling, no rounding) gives exactly
    scaled parameters when the bounds scale too.

    Regression: least_squares stops on xtol and gtol, which mix the
    parameters' units (norm(step) over R ~ 1e7 and C ~ 1e-12 alike) and
    the cost's (Ohm^2), and the x_scale floor of 1e-10 rescaled the small C
    differently: the same spectrum in other units stopped elsewhere (stress
    test case rc/24, here 79x apart). The absolute PARAMETER_BOUNDS still
    steer the trust region (scipy's trf scales by the distance to the
    bounds), so this needs scaled bounds.
    """
    f = np.logspace(np.log10(19071.51670080574), np.log10(0.0035536906006914185), 95)
    truth = [530028.4117574899, 20022831.072997387, 1.7818095897842149e-12,
             5310970.904689248, 1.4246395569165385e-10]
    start = [739208.3200348469, 34032330.40650132, 1.1750296971834518e-12,
             11024351.515477987, 1.487868908363994e-10]
    noise = np.random.default_rng(0).standard_normal((2, f.size))
    powers = np.array([1, 1, -1, 1, -1])
    lower = np.array([1e-4, 1e-4, 1e-15, 1e-4, 1e-15])
    upper = np.array([1e10, 1e10, 1e-1, 1e10, 1e-1])

    results = []
    for k in (1.0, 0.5):
        circuit = R(truth[0]) - (R(truth[1]) | C(truth[2])) - (R(truth[3]) | C(truth[4]))
        Z = circuit.impedance(f, truth)
        Z = k * (Z + 1e-3 * np.abs(Z) * (noise[0] + 1j * noise[1]))
        x0 = np.array(start) * k ** powers
        circuit = R(x0[0]) - (R(x0[1]) | C(x0[2])) - (R(x0[3]) | C(x0[4]))
        bounds = (lower * k ** powers, upper * k ** powers)
        fit, _ = fit_equivalent_circuit(f, Z, circuit, bounds=bounds)
        results.append(fit.params_opt / k ** powers)
    np.testing.assert_allclose(results[1], results[0], rtol=1e-9)


def test_bounds_override():
    """Explicit bounds replace PARAMETER_BOUNDS, also in bound_status."""
    f, Z = _blocking_spectrum()
    lower = [1.0, 1.0, 1e-13, 0.3, 1e-13, 0.3]
    upper = [5e6, 1e7, 1e-6, 1.0, 1e-6, 1.0]   # R_k capped below its truth
    fit, _ = fit_equivalent_circuit(f, Z, _start(), bounds=(lower, upper))
    assert fit.params_opt[1] <= 1e7
    assert fit.bound_status[1] == 'upper'


def test_bounds_override_rejects_empty_interval():
    """lower >= upper fails with the parameter named, not deep in scipy."""
    f, Z = _blocking_spectrum()
    lower = [1.0, 1.0, 1e-13, 0.3, 1e-13, 0.3]
    upper = [5e6, 1e7, 1e-6, 0.3, 1e-6, 1.0]   # n0: 0.3 .. 0.3
    with pytest.raises(ValueError, match='n0'):
        fit_equivalent_circuit(f, Z, _start(), bounds=(lower, upper))
