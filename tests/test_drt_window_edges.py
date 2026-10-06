#!/usr/bin/env python3
"""
Tests for the measured-window edge diagnostics of the DRT.

By default the tau grid spans exactly the measured window (tau = 1/(2*pi*f)),
so both failure modes here are about the edges of that window:

- a peak closer than DRT_PEAK_EDGE_DECADES to an end has one flank that no
  measurement constrains, and the basin integral behind its R_estimate is
  truncated by the end of the array;
- non-negative NNLS cannot represent response whose time constant lies
  outside the grid, so it piles gamma into the outermost bin instead.

With tau_extend_decades the grid runs past the slow end of the window; the
peak flag still measures against the window, and mass placed past it is
reported as extrapolated (tests at the end of the file).

Neither is reported unless these diagnostics fire, which is what the tests
lock down. Also covers the lambda probe reporting the span it actually
achieved after clipping at the probe bounds.
"""

import numpy as np
import pytest

from eis_analysis.drt import calculate_drt
from eis_analysis.drt.core import _detect_peaks
from eis_analysis.drt.significance import PeakSignificanceResult
from eis_analysis.drt.estimation import (
    _edge_pile_up,
    _flag_boundary_peaks,
    _lf_rc_ratio,
)
from eis_analysis.drt.stability import (
    probe_lambda_stability,
    MIN_PROBES_FOR_STABLE,
    PROBE_LAMBDA_MAX,
)
from eis_analysis.drt.linear_system import _build_drt_matrices
from eis_analysis.fitting.config import (
    DRT_EDGE_BIN_RPOL_FRACTION,
    DRT_EXTRAPOLATED_RPOL_FRACTION,
    DRT_LF_RC_RATIO_MIN,
    DRT_PEAK_EDGE_DECADES,
    DRT_TAU_EXTEND_STEPS,
)


def _voigt_impedance(frequencies, R_inf, elements):
    """Ideal Voigt: R_inf + sum_i R_i / (1 + j*omega*tau_i)."""
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)
    for R, tau in elements:
        Z += R / (1 + 1j * omega * tau)
    return Z


FREQUENCIES = np.logspace(5, -2, 71)  # 100 kHz .. 10 mHz
TAU_GRID = np.logspace(-6, 1, 100)    # 7 decades, ascending
TAU_WINDOW = (TAU_GRID[0], TAU_GRID[-1])


# =============================================================================
# Boundary-peak flag
# =============================================================================

# TAU_WINDOW spans log10(tau) from -6 to 1, so the distance to the nearer end is
# min(log_tau + 6, 1 - log_tau) and the threshold sits 0.7 decade in.
FLAG_CASES = [  # (log_tau, expected_distance, expected_flag)
    (-2.5, 3.5, False),                                # middle of the window
    (-5.8, 0.2, True),                                 # near the fast end
    (0.8, 0.2, True),                                  # near the slow end
    (1.3, -0.3, True),                                 # past the slow end
    (-6 + DRT_PEAK_EDGE_DECADES + 0.05,
     DRT_PEAK_EDGE_DECADES + 0.05, False),             # just inside the threshold
    (-6 + DRT_PEAK_EDGE_DECADES - 0.05,
     DRT_PEAK_EDGE_DECADES - 0.05, True),              # just outside it
]


def test_peak_is_flagged_by_distance_to_the_window_edge():
    """Both ends are treated the same, and the cut is the documented constant."""
    for log_tau, expected_distance, expected_flag in FLAG_CASES:
        peaks = [{'tau': 10 ** log_tau}]
        _flag_boundary_peaks(TAU_WINDOW, peaks, 'tau')

        assert peaks[0]['edge_distance_decades'] == pytest.approx(
            expected_distance, abs=0.01
        ), log_tau
        assert peaks[0]['boundary_sensitive'] is expected_flag, log_tau
        assert peaks[0]['outside_window'] is (expected_distance < 0), log_tau


def test_detect_peaks_annotates_both_paths():
    """Every peak dict leaving _detect_peaks carries the flag."""
    gamma = np.exp(-0.5 * ((np.log10(TAU_GRID) + 3.0) / 0.3) ** 2) * 100.0
    # A single maximum is never tested (delta_chi2 = inf), so it is significant
    peak = np.array([int(np.argmax(gamma))])
    significance = PeakSignificanceResult(
        candidates=peak, delta_chi2=np.array([np.inf]), significant=peak,
        noise_sigma=0.01, threshold=16.0)

    gmm_peaks, _, scipy_peaks = _detect_peaks(
        TAU_GRID, gamma, 'gmm', n_data=len(FREQUENCIES), tau_window=TAU_WINDOW,
        significance=significance
    )

    assert scipy_peaks
    for peak in scipy_peaks:
        assert 'boundary_sensitive' in peak
        assert 'edge_distance_decades' in peak
    if gmm_peaks:
        for peak in gmm_peaks:
            assert 'boundary_sensitive' in peak


# =============================================================================
# Edge-bin pile-up
# =============================================================================

def _gaussian_gamma(grid, centre=-2.5, sigma=0.4, amp=100.0):
    return np.exp(-0.5 * ((np.log10(grid) - centre) / sigma) ** 2) * amp


def test_smooth_distribution_has_no_pile_up():
    """Gamma rises from both edges towards its peak, so neither end is loaded."""
    gamma = _gaussian_gamma(TAU_GRID)
    d_ln_tau = float(np.mean(np.diff(np.log(TAU_GRID))))
    R_pol = float(np.sum(gamma) * d_ln_tau)

    fraction, end = _edge_pile_up(gamma, d_ln_tau, R_pol)
    assert fraction == 0.0
    assert end is None


def test_pile_up_is_detected_and_attributed_to_its_end():
    """Mass heaped against either boundary is reported, with the side named."""
    for end_index, expected_end in [(-1, 'high'), (0, 'low')]:
        gamma = _gaussian_gamma(TAU_GRID)
        d_ln_tau = float(np.mean(np.diff(np.log(TAU_GRID))))
        gamma[end_index] = float(np.sum(gamma))  # a lobe as heavy as the rest
        R_pol = float(np.sum(gamma) * d_ln_tau)

        fraction, end = _edge_pile_up(gamma, d_ln_tau, R_pol)
        assert fraction > DRT_EDGE_BIN_RPOL_FRACTION
        assert end == expected_end, end_index


def test_pile_up_measure_is_independent_of_grid_density():
    """
    Regression: a per-bin measure diluted as n_tau grew, so the same data
    warned at -n 100 and went silent at -n 400. The falling-run measure must
    read the same lobe at any grid density.
    """
    fractions = []
    for n_tau in (100, 200, 400, 800):
        grid = np.logspace(-6, 1, n_tau)
        gamma = _gaussian_gamma(grid)
        # A lobe of fixed width in decades, falling away from the slow end.
        decades_from_end = np.log10(grid[-1]) - np.log10(grid)
        gamma += 200.0 * np.exp(-decades_from_end / 0.3)
        d_ln_tau = float(np.mean(np.diff(np.log(grid))))
        R_pol = float(np.sum(gamma) * d_ln_tau)
        fractions.append(_edge_pile_up(gamma, d_ln_tau, R_pol)[0])

    assert min(fractions) > DRT_EDGE_BIN_RPOL_FRACTION
    # Relative, because absolute spread means nothing without the level. The
    # per-bin measure this replaced ran 0.068 -> 0.038 over the same 4x
    # refinement, a ratio of 1.8 that carried it across the threshold; what is
    # left here is rectangle-rule convergence, not grid dependence.
    assert max(fractions) / min(fractions) < 1.1, (
        f"pile-up fraction drifts with n_tau: {fractions}"
    )


def test_zero_rpol_is_not_a_division():
    """An empty DRT reports no pile-up rather than raising."""
    assert _edge_pile_up(np.zeros(100), 0.1, 0.0) == (0.0, None)


# =============================================================================
# End-to-end through calculate_drt
# =============================================================================

def test_clean_interior_spectrum_reports_no_edge_problems():
    """Two well-centred RC peaks: no flags, no pile-up warning."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')

    assert r.diagnostics.n_boundary_peaks == 0
    assert r.diagnostics.edge_pile_up_fraction < DRT_EDGE_BIN_RPOL_FRACTION
    assert not any('window' in w for w in r.diagnostics.nnls.warnings)
    assert not any(p.get('edge_contaminated') for p in r.diagnostics.scipy_peaks)


def test_process_slower_than_the_window_is_reported():
    """
    An RC process an order of magnitude slower than the lowest measured
    frequency cannot be placed on the grid; NNLS piles it into the last bin.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (5000.0, 50.0)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')
    diag = r.diagnostics

    assert diag.edge_pile_up_fraction > DRT_EDGE_BIN_RPOL_FRACTION
    assert diag.edge_pile_up_end == 'high'
    assert any('slow end of the tau grid' in w for w in diag.nnls.warnings)


def test_the_peak_that_swallowed_the_pile_up_is_marked():
    """
    Regression: the basin partition runs to the end of the array, so the
    outermost peak absorbs out-of-window mass however far away it sits - here
    2.8 decades from the edge, well clear of any distance test. Reporting the
    pile-up while printing its inflated R_estimate unmarked was the gap.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (5000.0, 50.0)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='scipy')
    peaks = r.diagnostics.scipy_peaks

    assert peaks, "expected at least one detected peak"
    contaminated = [p for p in peaks if p.get('edge_contaminated')]
    assert len(contaminated) == 1
    assert contaminated[0] is peaks[-1], "the slow-end peak absorbs a 'high' lobe"
    assert contaminated[0]['boundary_sensitive'] is False, (
        "this is exactly the case a distance test cannot catch"
    )
    assert contaminated[0]['R_estimate'] > 1.5 * 1000.0, "R should be visibly inflated"
    assert any('inflated' in w for w in r.diagnostics.nnls.warnings)


def test_gmm_marks_every_peak_because_rpol_is_the_divisor():
    """GMM splits R_pol by component weight, so pile-up taints all of them."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (5000.0, 50.0)])

    r = calculate_drt(FREQUENCIES, Z, peak_method='gmm')
    if not r.peaks:
        pytest.skip("GMM did not converge on this spectrum")

    assert all(p.get('edge_contaminated') for p in r.peaks)


# =============================================================================
# Lambda-probe span after clipping
# =============================================================================

def test_probe_reports_full_span_when_nothing_clips():
    """An interior lambda* gets the full two-decade sweep."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    matrices = _build_drt_matrices(FREQUENCIES, Z, 100.0, 100)

    stability = probe_lambda_stability(
        matrices, 1e-3, [(1e-3, 1000.0)], Z, 100.0
    )

    assert stability.n_clipped == 0
    assert stability.span_decades == pytest.approx(2.0, abs=0.01)
    assert not any('clipped' in w for w in stability.warnings)


def test_probe_near_the_upper_bound_reports_the_narrowed_span():
    """
    At lambda* = PROBE_LAMBDA_MAX / 2 both upward probes clip to
    PROBE_LAMBDA_MAX and collapse into one, so 'stable' covers barely over a
    decade, not two. The old code clipped silently.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    matrices = _build_drt_matrices(FREQUENCIES, Z, 100.0, 100)
    lambda_star = PROBE_LAMBDA_MAX / 2

    stability = probe_lambda_stability(
        matrices, lambda_star, [(1e-3, 1000.0)], Z, 100.0
    )

    assert stability.n_clipped == 2
    expected = np.log10(PROBE_LAMBDA_MAX / (lambda_star / 10.0))
    assert stability.span_decades == pytest.approx(expected, abs=0.01)
    assert stability.span_decades < 2.0
    assert any('clipped' in w for w in stability.warnings)


def test_a_single_surviving_probe_cannot_certify_stability():
    """
    Regression: at lambda* = 10 * PROBE_LAMBDA_MAX all four requested probes
    clip onto PROBE_LAMBDA_MAX and dedup to one, and 'persistence == n_probes' then
    awarded STABLE off a single re-solve.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    matrices = _build_drt_matrices(FREQUENCIES, Z, 100.0, 100)

    stability = probe_lambda_stability(
        matrices, 10 * PROBE_LAMBDA_MAX, [(1e-3, 1000.0)], Z, 100.0
    )

    n_successful = sum(1 for p in stability.probe_points if p.success)
    assert n_successful < MIN_PROBES_FOR_STABLE
    assert stability.peak_stability
    assert all(p.verdict != 'stable' for p in stability.peak_stability)
    assert any('marginal' in w for w in stability.warnings)


# =============================================================================
# Tau-grid extension past the slow end of the window
# =============================================================================

# RC 0.5 decade past the window (tau_max = 1/(2*pi*10 mHz) = 15.9 s).
SLOW_ELEMENTS = [(1000.0, 1e-3), (5000.0, 50.0)]
TAU_MAX = 1 / (2 * np.pi * FREQUENCIES.min())


def test_extension_keeps_the_spacing_and_the_window_points():
    """n_tau points stay on the window; the extension continues at the same step."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, SLOW_ELEMENTS)
    plain = _build_drt_matrices(FREQUENCIES, Z, 100.0, 100)
    extended = _build_drt_matrices(FREQUENCIES, Z, 100.0, 100, tau_extend_decades=1.0)

    assert np.array_equal(extended.tau[:100], plain.tau)
    assert extended.d_ln_tau == plain.d_ln_tau
    assert np.allclose(np.diff(np.log(extended.tau)), plain.d_ln_tau, rtol=1e-9)
    assert extended.tau[-1] >= TAU_MAX * 10 * (1 - 1e-9)  # at least one decade
    assert len(extended.tau) - 100 == int(np.ceil(np.log(10) / plain.d_ln_tau))
    assert extended.L.shape == (len(extended.tau) - 2, len(extended.tau))
    assert extended.tau_window == pytest.approx((plain.tau[0], plain.tau[-1]))


def test_extension_places_the_slow_process_past_the_window():
    """
    Without the extension the 50 s RC piles into the last bin (see above); with
    one decade and auto-lambda it becomes a peak past the window. Measured:
    tau off by 0.07 decade, R off by 18 %, pile-up gone.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, SLOW_ELEMENTS)
    r = calculate_drt(FREQUENCIES, Z, auto_lambda=True, tau_extend_decades=1.0)
    diag = r.diagnostics

    slow = max(diag.scipy_peaks, key=lambda p: p['tau'])
    assert abs(np.log10(slow['tau'] / 50.0)) < 0.15
    assert abs(slow['R_estimate'] - 5000.0) / 5000.0 < 0.25
    assert slow['outside_window'] and slow['boundary_sensitive']
    assert diag.edge_pile_up_fraction < DRT_EDGE_BIN_RPOL_FRACTION
    assert diag.R_pol_extrapolated_fraction > DRT_EXTRAPOLATED_RPOL_FRACTION
    assert any('past the measured window' in w for w in diag.nnls.warnings)
    assert diag.tau_extend_decades == 1.0 and diag.tau_extend_note is None


def test_clean_spectrum_extrapolates_nothing():
    """Peaks well inside the window leave the extension empty and quiet."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    r = calculate_drt(FREQUENCIES, Z, tau_extend_decades=1.0)

    assert r.diagnostics.R_pol_extrapolated_fraction < 0.01
    assert not any('past the measured window' in w for w in r.diagnostics.nnls.warnings)


def test_auto_leaves_a_clean_spectrum_alone():
    """No slow-end pile-up, no extension: identical to the default grid."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, [(1000.0, 1e-3), (2000.0, 1e-1)])
    plain = calculate_drt(FREQUENCIES, Z, auto_lambda=True)
    auto = calculate_drt(FREQUENCIES, Z, auto_lambda=True, tau_extend_decades='auto')

    assert auto.diagnostics.tau_extend_decades == 0.0
    assert 'not needed' in auto.diagnostics.tau_extend_note
    assert np.array_equal(auto.gamma, plain.gamma)


def test_auto_takes_the_smallest_step_that_resolves_the_pile_up():
    """The 50 s RC: 0.5 decade still piles up, 1.0 closes it (measured)."""
    Z = _voigt_impedance(FREQUENCIES, 100.0, SLOW_ELEMENTS)
    r = calculate_drt(FREQUENCIES, Z, auto_lambda=True, tau_extend_decades='auto')

    assert r.diagnostics.tau_extend_decades == 1.0
    assert 'resolved' in r.diagnostics.tau_extend_note
    assert r.diagnostics.edge_pile_up_fraction < DRT_EDGE_BIN_RPOL_FRACTION


def test_auto_keeps_the_plain_grid_when_no_step_resolves():
    """
    At the fixed lambda = 0.1 the slow mass runs to the end of every extended
    grid (pile-up 62 -> 80 -> 88 -> 92 %); 'auto' must not pick a step that
    only moves the heap further out.
    """
    Z = _voigt_impedance(FREQUENCIES, 100.0, SLOW_ELEMENTS)
    plain = calculate_drt(FREQUENCIES, Z)
    auto = calculate_drt(FREQUENCIES, Z, tau_extend_decades='auto')

    assert auto.diagnostics.tau_extend_decades == 0.0
    assert f'up to {DRT_TAU_EXTEND_STEPS[-1]}' in auto.diagnostics.tau_extend_note
    assert np.array_equal(auto.gamma, plain.gamma)


def test_auto_does_not_extend_a_capacitive_end():
    """
    Blocking electrode (series C): extending turned the capacitor into a huge
    extrapolated peak and pushed the real 1000 Ohm arc below detection.
    """
    omega = 2 * np.pi * FREQUENCIES
    Z = _voigt_impedance(FREQUENCIES, 10.0, [(1000.0, 1e-3)]) + 1 / (1j * omega * 1e-3)
    # Fixed lambda: against the capacitive heap the RC peak is only 3.1-3.3 %
    # of max gamma for lambda 1e-8..1e-5, right at the 3 % detection
    # threshold, so auto-lambda would decide the test. At 1e-9 it is 5.1 %.
    r = calculate_drt(FREQUENCIES, Z, lambda_reg=1e-9, tau_extend_decades='auto')

    assert r.diagnostics.tau_extend_decades == 0.0
    assert 'capacitive' in r.diagnostics.tau_extend_note
    assert any(abs(np.log10(p['tau'] / 1e-3)) < 0.15 for p in r.diagnostics.scipy_peaks)


def test_lf_ratio_matches_a_single_rc():
    """
    r = 2x/(1+x^2) for one RC, x = omega * tau taken at the geometric centre
    of the four lowest frequencies the slope is fitted over (not at f_min).
    """
    for log_x in [0.0, 0.5, 1.0]:
        tau = 10 ** log_x / (2 * np.pi * FREQUENCIES.min())
        Z = _voigt_impedance(FREQUENCIES, 10.0, [(1000.0, tau)])
        omega_centre = 2 * np.pi * np.exp(np.mean(np.log(np.sort(FREQUENCIES)[:4])))
        x = omega_centre * tau

        assert _lf_rc_ratio(FREQUENCIES, Z) == pytest.approx(2 * x / (1 + x ** 2), rel=0.05), log_x


def test_lf_ratio_of_a_capacitor_and_of_an_inductive_end():
    """Series C: r = 0, well under the guard. No capacitive response: inf."""
    omega = 2 * np.pi * FREQUENCIES
    assert _lf_rc_ratio(FREQUENCIES, 10.0 + 1 / (1j * omega * 1e-3)) < DRT_LF_RC_RATIO_MIN / 10
    assert _lf_rc_ratio(FREQUENCIES, 10.0 + 1j * omega * 1e-6) == float('inf')


def test_invalid_extension_is_rejected():
    for bad in [-1.0, float("inf"), 'always', None]:
        Z = _voigt_impedance(FREQUENCIES, 100.0, SLOW_ELEMENTS)
        with pytest.raises(ValueError, match="tau_extend_decades"):
            calculate_drt(FREQUENCIES, Z, tau_extend_decades=bad)


def test_cli_tau_extend_parses_and_rejects():
    import argparse
    from eis_analysis.cli.parser import _tau_extend
    for text, expected in [('auto', 'auto'), ('1', 1.0), ('0.5', 0.5)]:
        assert _tau_extend(text) == expected, text
    for text in ['-1', 'maybe', 'inf', 'nan']:
        with pytest.raises(argparse.ArgumentTypeError):
            _tau_extend(text)
