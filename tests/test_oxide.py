"""
Tests for oxide layer analysis (analysis/oxide.py).

Regression tests for audit findings (2026-07-02):
- O2: estimate_permittivity() must not log a bogus "Oxide thickness"
  computed from a dummy epsilon_r.
- O3: silent assumptions must be visible (candidate listing, CPE n warning,
  high-frequency fallback checks and median estimate).
"""

import logging

import numpy as np

from eis_analysis.analysis.config import BRUG_HM_DIVERGENCE_MAX, EPSILON_0
from eis_analysis.analysis.oxide import analyze_oxide_layer, estimate_permittivity
from eis_analysis.fitting.circuit import FitResult
from eis_analysis.fitting.circuit_elements import C, CC, K, L, Q, R

OXIDE_LOGGER = 'eis_analysis.analysis.oxide'

R_S = 100.0    # Series resistance [Ohm]
R_P = 5000.0   # Parallel (oxide) resistance [Ohm]
C_P = 1e-6     # Oxide capacitance [F]


def _synthetic_voigt():
    """R_S - (R_P || C_P) impedance over 10 mHz .. 100 kHz."""
    freq = np.logspace(5, -2, 50)
    omega = 2 * np.pi * freq
    Z = R_S + R_P / (1 + 1j * omega * R_P * C_P)
    return freq, Z


def _fit_result_voigt():
    """FitResult wrapper around the known circuit (oxide.py reads only .circuit)."""
    circuit = R(R_S) - (R(R_P) | C(C_P))
    params = np.array([R_S, R_P, C_P])
    return FitResult(
        circuit=circuit,
        params_opt=params,
        params_stderr=np.zeros_like(params),
        fit_error_rel=0.1,
    )


def test_analyze_oxide_layer_thickness_from_circuit():
    freq, Z = _synthetic_voigt()
    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt())

    assert oxide is not None
    assert oxide.element_type == 'C'
    assert abs(oxide.element_R - R_P) < 1e-9
    assert abs(oxide.capacitance - C_P) / C_P < 1e-9

    d_nm_expected = EPSILON_0 * 22.0 / C_P * 1e7
    assert abs(oxide.thickness_nm - d_nm_expected) / d_nm_expected < 1e-9


def test_estimate_permittivity_does_not_log_thickness(caplog):
    """Regression (audit O2): no 'Oxide thickness' line from dummy epsilon_r."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        result = estimate_permittivity(
            freq, Z, thickness_nm=19.5, fit_result=_fit_result_voigt()
        )

    assert result is not None and result.permittivity is not None
    assert 'Oxide thickness' not in caplog.text
    assert 'Permittivity' in caplog.text


def test_estimate_permittivity_fallback_does_not_log_thickness(caplog):
    """Regression (audit O2): same for the high-frequency fallback path."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        result = estimate_permittivity(freq, Z, thickness_nm=19.5)

    assert result is not None and result.permittivity is not None
    assert 'Oxide thickness' not in caplog.text


def test_permittivity_thickness_roundtrip():
    """estimate_permittivity() is the exact inverse of analyze_oxide_layer()."""
    freq, Z = _synthetic_voigt()
    fit_result = _fit_result_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)
    inverse = estimate_permittivity(
        freq, Z, thickness_nm=oxide.thickness_nm, fit_result=fit_result
    )

    assert abs(inverse.permittivity - 22.0) / 22.0 < 1e-9


def test_permittivity_result_carries_input_thickness():
    """In inverse mode the thickness is the input, not a derived value."""
    freq, Z = _synthetic_voigt()

    result = estimate_permittivity(
        freq, Z, thickness_nm=19.5, fit_result=_fit_result_voigt()
    )

    assert result is not None
    assert result.thickness_nm == 19.5
    assert result.thickness_brug_nm is None
    assert result.element_type == 'C'


def test_permittivity_brug_comparison():
    """Brug ε_r is reported alongside Hsu-Mansfeld and scales with C."""
    freq, Z = _synthetic_voigt()

    result = estimate_permittivity(
        freq, Z, thickness_nm=19.5, fit_result=_fit_result_voigt_q(0.9)
    )

    assert result is not None
    assert result.permittivity_brug is not None
    # ε_r is linear in C_specific, so both ratios must agree exactly
    eps_ratio = result.permittivity_brug / result.permittivity
    C_ratio = result.capacitance_specific_brug / result.capacitance_specific
    assert abs(eps_ratio - C_ratio) / C_ratio < 1e-9


def test_permittivity_brug_none_without_series_R():
    """No series R -> no Brug ε_r, mirroring the thickness path."""
    freq, Z = _synthetic_voigt()
    n = 0.9
    fit_result = _fit_result(R(R_P) | Q(Q_VAL, n), [R_P, Q_VAL, n])

    result = estimate_permittivity(
        freq, Z, thickness_nm=19.5, fit_result=fit_result
    )

    assert result is not None
    assert result.permittivity is not None
    assert result.permittivity_brug is None


# --- Audit O3: candidate listing and selection assumption ---

def _fit_result_two_voigts():
    """Two Voigt elements; the second (R=5000) is the dominant barrier."""
    circuit = R(R_S) - (R(1000.0) | C(1e-5)) - (R(R_P) | C(C_P))
    params = np.array([R_S, 1000.0, 1e-5, R_P, C_P])
    return FitResult(
        circuit=circuit,
        params_opt=params,
        params_stderr=np.zeros_like(params),
        fit_error_rel=0.1,
    )


def test_candidates_listed_and_assumption_noted(caplog):
    """Regression (audit O3): all candidates logged, selection assumption stated."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(
            freq, Z, epsilon_r=22.0, fit_result=_fit_result_two_voigts()
        )

    assert oxide is not None
    assert oxide.element_R == R_P  # larger R wins
    assert '[1] C: R = 1000.0' in caplog.text
    assert '[2] C: R = 5000.0' in caplog.text
    assert 'Selection assumes the largest-R element' in caplog.text


# --- Audit O3: CPE exponent warning ---

def _fit_result_voigt_q(n):
    circuit = R(R_S) - (R(R_P) | Q(2e-6, n))
    params = np.array([R_S, R_P, 2e-6, n])
    return FitResult(
        circuit=circuit,
        params_opt=params,
        params_stderr=np.zeros_like(params),
        fit_error_rel=0.1,
    )


def test_cpe_low_n_warns(caplog):
    """Regression (audit O3): n < 0.8 -> C_eff not well-defined warning."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(
            freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt_q(0.7)
        )

    assert oxide is not None
    assert 'not well-defined' in caplog.text


def test_cpe_high_n_no_warning(caplog):
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(
            freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt_q(0.9)
        )

    assert oxide is not None
    assert 'not well-defined' not in caplog.text


# --- Traversal and CPE conversion (audit 2026-07-02, priority 4) ---

K_TAU = 1e-3   # K element time constant [s] -> C = K_TAU/R_P = 2e-7 F
Q_VAL = 2e-6   # CPE coefficient used by _fit_result_voigt_q()


def _fit_result(circuit, params):
    return FitResult(
        circuit=circuit,
        params_opt=np.asarray(params),
        params_stderr=np.zeros(len(params)),
        fit_error_rel=0.1,
    )


def test_k_element_traversal():
    """K element is found and converted via C = tau/R."""
    freq, Z = _synthetic_voigt()
    fit_result = _fit_result(R(R_S) - K(R_P, K_TAU), [R_S, R_P, K_TAU])

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_type == 'K'
    assert abs(oxide.element_R - R_P) < 1e-9
    assert abs(oxide.element_tau - K_TAU) < 1e-12
    C_expected = K_TAU / R_P
    assert abs(oxide.capacitance - C_expected) / C_expected < 1e-9


def test_cpe_hsu_mansfeld_conversion():
    """C_eff of a dominant CPE follows Hsu-Mansfeld: (R*Q)^(1/n)/R."""
    freq, Z = _synthetic_voigt()
    n = 0.9

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt_q(n))

    assert oxide is not None
    assert oxide.element_type == 'Q'
    C_expected = (R_P * Q_VAL) ** (1.0 / n) / R_P
    assert abs(oxide.capacitance - C_expected) / C_expected < 1e-9
    tau_expected = R_P * C_expected
    assert abs(oxide.element_tau - tau_expected) / tau_expected < 1e-9


def test_cpe_brug_conversion():
    """Brug (2D) comparison: C = Q^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)."""
    freq, Z = _synthetic_voigt()
    n = 0.9

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt_q(n))

    assert oxide is not None
    C_brug_expected = Q_VAL ** (1.0 / n) * (1.0 / R_S + 1.0 / R_P) ** ((n - 1.0) / n)
    assert abs(oxide.capacitance_brug - C_brug_expected) / C_brug_expected < 1e-9

    d_brug_expected = EPSILON_0 * 22.0 / C_brug_expected * 1e7
    assert abs(oxide.thickness_brug_nm - d_brug_expected) / d_brug_expected < 1e-9


def test_cpe_brug_equals_hsu_mansfeld_at_n_one():
    """At n = 1 both conversions must reduce to C = Q."""
    freq, Z = _synthetic_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt_q(1.0))

    assert oxide is not None
    assert abs(oxide.capacitance - Q_VAL) / Q_VAL < 1e-9
    assert abs(oxide.capacitance_brug - Q_VAL) / Q_VAL < 1e-9


def test_cpe_brug_unavailable_without_series_R(caplog):
    """No series R in circuit -> Brug fields None, informative log."""
    freq, Z = _synthetic_voigt()
    n = 0.9
    circuit = R(R_P) | Q(Q_VAL, n)
    fit_result = _fit_result(circuit, [R_P, Q_VAL, n])

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.capacitance_brug is None
    assert oxide.thickness_brug_nm is None
    assert 'Brug (2D) estimate not available' in caplog.text


def test_voigt_c_element_has_no_brug_fields():
    """Brug conversion applies only to Q elements - ideal C gets None."""
    freq, Z = _synthetic_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_voigt())

    assert oxide is not None
    assert oxide.capacitance_brug is None
    assert oxide.capacitance_specific_brug is None
    assert oxide.thickness_brug_nm is None


def test_mixed_voigt_k_traversal():
    """Voigt and K candidates in series are both found; larger R wins."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(1000.0) | C(1e-5)) - K(R_P, K_TAU)
    fit_result = _fit_result(circuit, [R_S, 1000.0, 1e-5, R_P, K_TAU])

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_type == 'K'
    assert oxide.element_R == R_P


# --- Audit O4: traversal robustness ---

def test_k_element_zero_R_skipped(caplog):
    """Regression (audit O4): K with R=0 must not raise ZeroDivisionError."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - K(0.0, 1e-4) - (R(R_P) | C(C_P))
    fit_result = _fit_result(circuit, [R_S, 0.0, 1e-4, R_P, C_P])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_type == 'C'
    assert oxide.element_R == R_P
    assert 'non-positive R' in caplog.text


def test_multiple_R_in_parallel_warns(caplog):
    """Regression (audit O4): (R1|R2|C) warns instead of silently taking the last R."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(1000.0) | R(2000.0) | C(C_P))
    fit_result = _fit_result(circuit, [R_S, 1000.0, 2000.0, C_P])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_R == 2000.0  # last one wins (documented behavior)
    assert 'Multiple R elements' in caplog.text


def test_multiple_cap_in_parallel_warns(caplog):
    """(R|C1|C2): the larger capacitance wins, not the last one written.

    Was audit O4's "using the last one". Position in the expression is not a
    physical criterion; both capacitors are now candidates in their own right
    and are ranked by size, since they share the one resistance and the
    largest-R heuristic cannot separate them.
    """
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(R_P) | C(1e-5) | C(C_P))
    fit_result = _fit_result(circuit, [R_S, R_P, 1e-5, C_P])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert abs(oxide.capacitance - 1e-5) / 1e-5 < 1e-9   # the larger one
    assert 'share one parallel resistance' in caplog.text
    assert 'not separately identifiable' in caplog.text


# --- Audit O3: high-frequency fallback (Mode 2) ---

def test_hf_fallback_median_estimate():
    """Fallback C is the median over the top frequency decade, close to C_P."""
    freq, Z = _synthetic_voigt()
    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert oxide is not None
    assert oxide.element_type == 'estimate'
    # omega*R_P*C_P >= ~314 in the top decade, so C_i ~ C_P within ~1e-5
    assert abs(oxide.capacitance - C_P) / C_P < 1e-3


def test_hf_fallback_series_combination_warning(caplog):
    """Regression (audit O3): fallback warns about series capacitance combination."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert 'series combination' in caplog.text


def test_hf_fallback_spread_warning(caplog):
    """Regression (audit O3): warn when omega*R*C >> 1 does not hold in the decade.

    R_P*C = 5e-6 s puts the characteristic frequency (~32 kHz) inside the
    top decade, so C_i = -1/(omega*Z'') drifts by ~10x across it.
    """
    freq = np.logspace(5, -2, 50)
    omega = 2 * np.pi * freq
    Z = R_S + R_P / (1 + 1j * omega * R_P * 1e-9)

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert oxide is not None
    assert 'may not hold' in caplog.text


def test_hf_fallback_settled_estimate_no_spread_warning(caplog):
    """Series R does not invalidate C = -1/(omega*Z'') -> no spread warning."""
    freq = np.logspace(5, -2, 50)
    omega = 2 * np.pi * freq
    # Series R-C: C_i is exact at every frequency despite resistive phase
    Z = 1000.0 - 1j / (omega * 1e-6)

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert abs(oxide.capacitance - 1e-6) / 1e-6 < 1e-9
    assert 'may not hold' not in caplog.text


def test_hf_fallback_inductive_data(caplog):
    """No capacitive point in the top decade -> single-point path with warning."""
    freq = np.logspace(5, -2, 50)
    omega = 2 * np.pi * freq
    Z = 100.0 + 1j * omega * 1e-6  # inductive everywhere

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert oxide is not None  # pre-0.16.16 behavior preserved
    assert 'inductive' in caplog.text


def test_brug_suppressed_when_series_R_at_optimizer_floor(caplog):
    """R_s driven to the optimizer floor -> Brug suppressed, not silently absurd.

    Regression for the ZrO2 permittivity report of 2026-08-18: a CPE with
    n < 1 mimics a series resistance at high frequency, so R_s collapsed to
    its 1e-4 Ohm lower bound. Brug's C ~ R_s^((1-n)/n) then differed from
    Hsu-Mansfeld by 233x and yielded eps_r = 1.55 (near vacuum) with no
    warning at all.
    """
    freq, Z = _synthetic_voigt()
    n = 0.819
    R_s_floored = 1e-4
    circuit = R(R_s_floored) - (R(R_P) | Q(Q_VAL, n))
    fit_result = _fit_result(circuit, [R_s_floored, R_P, Q_VAL, n])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = estimate_permittivity(freq, Z, thickness_nm=1925.0,
                                      fit_result=fit_result)

    # Hsu-Mansfeld (3D) is unaffected - it never uses R_s
    assert oxide is not None
    assert oxide.permittivity is not None
    # Brug is suppressed rather than reported as a bogus comparison
    assert oxide.capacitance_brug is None
    assert oxide.capacitance_specific_brug is None
    assert oxide.permittivity_brug is None
    assert 'did not identify it' in caplog.text


def test_brug_divergence_warning_above_threshold(caplog):
    """Plausible R_s but huge R_ct/R_s -> value kept, divergence flagged."""
    freq, Z = _synthetic_voigt()
    n = 0.819
    R_s_small, R_ct_large = 0.05, 2.6e7  # ratio^((1-n)/n) is far above 10x
    circuit = R(R_s_small) - (R(R_ct_large) | Q(Q_VAL, n))
    fit_result = _fit_result(circuit, [R_s_small, R_ct_large, Q_VAL, n])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.capacitance_brug is not None  # above the R_s floor, still reported
    ratio = oxide.capacitance / oxide.capacitance_brug
    assert ratio > BRUG_HM_DIVERGENCE_MAX
    # The ratio is exactly (1 + R_ct/R_s)^((1-n)/n)
    expected = (1.0 + R_ct_large / R_s_small) ** ((1.0 - n) / n)
    assert abs(ratio - expected) / expected < 1e-9
    assert 'do not bracket a single C_eff' in caplog.text


def test_brug_no_divergence_warning_for_healthy_fit(caplog):
    """Well-determined R_s and moderate R_ct -> no divergence warning."""
    freq, Z = _synthetic_voigt()
    n = 0.95  # near-ideal CPE keeps the two models close

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_voigt_q(n))

    assert oxide is not None
    assert oxide.capacitance_brug is not None
    assert oxide.capacitance / oxide.capacitance_brug < BRUG_HM_DIVERGENCE_MAX
    assert 'do not bracket a single C_eff' not in caplog.text


# ---------------------------------------------------------------------------
# Cole-Cole element (added 2026-08-27): a blocking dielectric, so it carries
# no parallel resistance and the largest-R selection rule cannot apply to it.
# ---------------------------------------------------------------------------

CC_C_INF = 5e-8    # High-frequency capacitance [F]
CC_DC = 9.5e-7     # Relaxation strength [F]
CC_TAU = 2e-3      # Relaxation time [s]
CC_ALPHA = 0.25    # Broadening exponent


def _fit_result_cc(extra=None):
    """R_S - CC(...), optionally with an extra element in series."""
    cc = CC(CC_C_INF, CC_DC, CC_TAU, CC_ALPHA)
    circuit = R(R_S) - cc if extra is None else R(R_S) - cc - extra
    return _fit_result(circuit, circuit.get_all_params())


def test_cc_element_traversal_uses_static_capacitance():
    """CC is found, and its capacitance is C_s = C_inf + dC, not C_inf."""
    freq, Z = _synthetic_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_cc())

    assert oxide is not None
    assert oxide.element_type == 'CC'
    C_expected = CC_C_INF + CC_DC
    assert abs(oxide.capacitance - C_expected) / C_expected < 1e-12
    assert abs(oxide.element_tau - CC_TAU) < 1e-15


def test_cc_element_reports_no_resistance():
    """A blocking dielectric has no DC path; element_R must be None, not 0."""
    freq, Z = _synthetic_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=_fit_result_cc())

    assert oxide is not None
    assert oxide.element_R is None


def test_cc_thickness_round_trips():
    """A film built to be 20 nm at eps_r = 22 is measured back as 20 nm."""
    freq, Z = _synthetic_voigt()
    eps_r, area_cm2, d_nm = 22.0, 1.0, 20.0
    C_s = EPSILON_0 * eps_r * area_cm2 / (d_nm * 1e-7)

    circuit = R(R_S) - CC(0.05 * C_s, 0.95 * C_s, CC_TAU, CC_ALPHA)
    fit_result = _fit_result(circuit, circuit.get_all_params())

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=eps_r, area_cm2=area_cm2,
                                fit_result=fit_result)

    assert oxide is not None
    assert abs(oxide.thickness_nm - d_nm) / d_nm < 1e-9


def test_cc_wins_over_a_larger_r_voigt_element():
    """An explicit CC is the dielectric; the largest-R heuristic does not override it."""
    freq, Z = _synthetic_voigt()
    # A Voigt element with a resistance far larger than anything else present
    circuit = R(R_S) - CC(CC_C_INF, CC_DC, CC_TAU, CC_ALPHA) - (R(1e9) | C(1e-9))
    fit_result = _fit_result(circuit, circuit.get_all_params())

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_type == 'CC'
    assert abs(oxide.capacitance - (CC_C_INF + CC_DC)) / (CC_C_INF + CC_DC) < 1e-12


def test_multiple_cc_elements_warn_and_pick_largest(caplog):
    """Two dielectric relaxations: largest C_s wins, and the user is told."""
    freq, Z = _synthetic_voigt()
    small = CC(1e-9, 1e-8, 1e-4, 0.1)
    big = CC(CC_C_INF, CC_DC, CC_TAU, CC_ALPHA)
    circuit = R(R_S) - small - big
    fit_result = _fit_result(circuit, circuit.get_all_params())

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert abs(oxide.capacitance - (CC_C_INF + CC_DC)) / (CC_C_INF + CC_DC) < 1e-12
    assert any("2 Cole-Cole elements" in r.message for r in caplog.records)


def test_cc_uses_fitted_values_not_the_initial_guess():
    """Regression: the traversal must read node.params, not node.C_inf/.dC.

    update_params() rewrites params but leaves the named attributes (and the
    C_static property) at their construction-time values, so reading the
    attributes reported the initial guess and the thickness came out wrong.
    """
    freq, Z = _synthetic_voigt()
    # Initial guess an order of magnitude off, so the two are unmistakable
    guess_C_inf, guess_dC, guess_tau = 1e-9, 1e-8, 1e-5
    circuit = R(R_S) - CC(guess_C_inf, guess_dC, guess_tau, 0.05)
    fitted = [R_S, CC_C_INF, CC_DC, CC_TAU, CC_ALPHA]
    circuit.update_params(fitted)

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result(circuit, fitted))

    assert oxide is not None
    C_fitted = CC_C_INF + CC_DC
    assert abs(oxide.capacitance - C_fitted) / C_fitted < 1e-12
    assert oxide.capacitance > 10 * (guess_C_inf + guess_dC)   # not the guess
    assert abs(oxide.element_tau - CC_TAU) < 1e-15


def test_cc_inside_a_parallel_leakage_branch_is_found():
    """R_leak || CC still reports the CC - the traversal recurses into Parallel."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(1e7) | CC(CC_C_INF, CC_DC, CC_TAU, CC_ALPHA))
    fit_result = _fit_result(circuit, circuit.get_all_params())

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert oxide.element_type == 'CC'
    assert oxide.element_R is None


# ---------------------------------------------------------------------------
# Cole-Cole capacitance vs. the measured frequency window (added 2026-08-29).
#
# C_s = C_inf + dC is the omega -> 0 limit of C*(omega) and is only what the
# data determine when the relaxation is actually traced. _synthetic_voigt()
# spans 1e-2 .. 1e5 Hz; CC_TAU = 2e-3 s puts f_char at ~80 Hz, well inside it,
# so every test above keeps exercising the static path.
#
# Reported failure: a fit landed on tau = 1e4 s (the upper bound of
# PARAMETER_BOUNDS['tau_CC']), i.e. f_char = 1.6e-5 Hz, three decades below
# f_min. dC was unidentified (95% CI spanning a factor of 4) yet dominated
# C_s, giving eps_r = 3491 where the high-frequency value C_inf gives 17.6.
# ---------------------------------------------------------------------------

# Values from the reported fit
BUG_C_INF = 8.117e-08   # High-frequency capacitance [F]
BUG_DC = 1.603e-05      # Relaxation strength [F] - unidentified in that fit
BUG_TAU = 1e4           # Relaxation time [s] - at its upper fitting bound


def _fit_result_cc_tau(tau, C_inf=BUG_C_INF, dC=BUG_DC):
    """R_S - CC(C_inf, dC, tau, alpha) with tau given (may be a str = fixed)."""
    circuit = R(R_S) - CC(C_inf, dC, tau, CC_ALPHA)
    return _fit_result(circuit, circuit.get_all_params())


def test_cc_below_window_reports_c_inf_not_static():
    """tau on the upper bound: f_char is below f_min, so only C_inf is measured."""
    freq, Z = _synthetic_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result_cc_tau(BUG_TAU))

    assert oxide is not None
    assert oxide.element_type == 'CC'
    assert abs(oxide.capacitance - BUG_C_INF) / BUG_C_INF < 1e-12
    # The static value would have been ~200x larger - that is the whole bug
    assert oxide.capacitance < 0.01 * (BUG_C_INF + BUG_DC)


def test_cc_below_window_permittivity_matches_reference_model():
    """The reported case end to end: eps_r must come out ~17.6, not 3491.

    192 nm at 1 cm2 is the geometry implied by the report - it turns the two
    candidate capacitances into the two permittivities that were quoted there,
    C_s -> 3491 (what was printed) and C_inf -> 17.6 (what the reference model
    gives, 17.4).
    """
    freq, Z = _synthetic_voigt()

    oxide = estimate_permittivity(freq, Z, thickness_nm=192.0, area_cm2=1.0,
                                  fit_result=_fit_result_cc_tau(BUG_TAU))

    assert oxide is not None
    assert 17.0 < oxide.permittivity < 18.0
    # What the static capacitance would have given, for contrast
    eps_r_static = (BUG_C_INF + BUG_DC) * (192.0 * 1e-7) / (1.0 * EPSILON_0)
    assert eps_r_static > 3000.0


def test_cc_below_window_warns_about_window_and_bound(caplog):
    """Both diagnostics fire, and they are separate statements."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                            fit_result=_fit_result_cc_tau(BUG_TAU))

    text = caplog.text
    assert 'BELOW the measured window' in text
    assert 'extrapolation to DC' in text
    assert 'upper fitting bound' in text


def test_cc_above_window_keeps_static_capacitance(caplog):
    """tau on the LOWER bound: the window sits at omega*tau << 1, so C_s holds."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_cc_tau(1e-9))

    assert oxide is not None
    C_s = BUG_C_INF + BUG_DC
    assert abs(oxide.capacitance - C_s) / C_s < 1e-12
    # The window rule drives the value; the bound test only warns
    assert 'ABOVE the measured window' in caplog.text
    assert 'lower fitting bound' in caplog.text
    assert 'BELOW the measured window' not in caplog.text


def test_cc_fixed_tau_warns_about_window_but_not_about_the_bound(caplog):
    """A tau pinned by the user is a choice, not an undetermined parameter."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_cc_tau(str(BUG_TAU)))

    assert oxide is not None
    assert abs(oxide.capacitance - BUG_C_INF) / BUG_C_INF < 1e-12
    assert 'BELOW the measured window' in caplog.text
    assert 'fitting bound' not in caplog.text


def test_cc_near_window_edge_keeps_static_but_warns(caplog):
    """f_char half a decade inside f_min: C_s stands, its determination does not."""
    freq, Z = _synthetic_voigt()
    f_min = float(np.min(freq))
    tau_edge = 1.0 / (2 * np.pi * f_min * 10 ** 0.5)   # f_char = f_min * 10^0.5

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_cc_tau(tau_edge))

    assert oxide is not None
    C_s = BUG_C_INF + BUG_DC
    assert abs(oxide.capacitance - C_s) / C_s < 1e-12
    assert 'marginally determined' in caplog.text


def test_cc_inside_window_logs_no_capacitance_warning(caplog):
    """The unchanged path: a traced relaxation reports C_s with no complaint."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_cc())

    assert oxide is not None
    C_s = CC_C_INF + CC_DC
    assert abs(oxide.capacitance - C_s) / C_s < 1e-12
    assert 'measured window' not in caplog.text
    assert 'fitting bound' not in caplog.text


# ---------------------------------------------------------------------------
# Selecting the dielectric element (added 2026-08-29).
#
# The old traversal only registered a capacitance that shared a Parallel with
# a resistance - a Voigt element. Reported failure: in L - R0 - (Q|C) the C
# has no resistance beside it, so the whole analysis fell through to the
# high-frequency spectral estimate (which itself reported a spread of 1.37
# across the top decade) even though C was a fitted parameter with a 0.7 %
# confidence interval. It landed 2 % from the reference value by luck.
#
# The criterion is now physical: the dielectric element is the one whose
# admittance rises as omega^n with n near 1, whatever its type and wherever it
# sits in the expression.
# ---------------------------------------------------------------------------

REPORTED_C = 9.06e-08   # The fitted capacitance the module used to ignore [F]


def _fit_result_l_r_qc(n=0.85, C_val=REPORTED_C):
    """The reported model L - R0 - (Q|C): a capacitance with no parallel R."""
    circuit = L(1e-7) - R(R_S) - (Q(3e-6, n) | C(C_val))
    return _fit_result(circuit, circuit.get_all_params())


def test_capacitance_without_parallel_resistance_is_found(caplog):
    """Regression: a fitted C with no R beside it must not fall through."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_l_r_qc())

    assert oxide is not None
    assert oxide.element_type == 'C'
    assert abs(oxide.capacitance - REPORTED_C) / REPORTED_C < 1e-12
    assert oxide.element_R is None          # there is no parallel resistance
    assert oxide.element_tau is None        # so there is no RC time constant
    # The whole point: no fallback
    assert 'NOT FROM THE FIT' not in caplog.text
    assert 'Falling back' not in caplog.text


def test_series_capacitance_is_found():
    """The same holds for a plain series C - no Parallel node at all."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - C(REPORTED_C)

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result(circuit,
                                                       circuit.get_all_params()))

    assert oxide is not None
    assert oxide.element_type == 'C'
    assert abs(oxide.capacitance - REPORTED_C) / REPORTED_C < 1e-12


def test_ideal_c_beats_a_cpe_in_the_same_parallel(caplog):
    """C (n = 1 exactly) outranks a CPE whose capacitance needs a model."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(R_P) | Q(3e-6, 0.95) | C(REPORTED_C))

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result(
                                        circuit, circuit.get_all_params()))

    assert oxide is not None
    assert oxide.element_type == 'C'        # not the Q, despite sharing R
    assert abs(oxide.capacitance - REPORTED_C) / REPORTED_C < 1e-12


def test_low_n_cpe_is_not_a_dielectric_but_still_beats_the_fallback(caplog):
    """n = 0.57 is transport, not a dielectric - say so, but stay on the fit."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_voigt_q(0.57))

    assert oxide is not None
    assert oxide.element_type == 'Q'                    # still used
    assert 'No dielectric element in circuit' in caplog.text
    assert 'no dielectric meaning' in caplog.text
    assert 'NOT FROM THE FIT' not in caplog.text        # not the spectral guess


def test_near_ideal_cpe_is_a_dielectric(caplog):
    """n = 0.9 is a near-ideal CPE - a dielectric, no "not a dielectric" warning."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result_voigt_q(0.9))

    assert oxide is not None
    assert oxide.element_type == 'Q'
    assert 'No dielectric element in circuit' not in caplog.text


def test_cpe_without_parallel_resistance_is_not_a_candidate(caplog):
    """Hsu-Mansfeld and Brug both need R; without it a Q cannot be converted."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - Q(3e-6, 0.95)

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                    fit_result=_fit_result(
                                        circuit, circuit.get_all_params()))

    assert oxide is not None
    assert oxide.element_type == 'estimate'     # nothing convertible -> fallback
    assert 'no parallel resistance' in caplog.text
    assert 'NOT FROM THE FIT' in caplog.text


def test_hf_fallback_says_the_value_is_not_from_the_fit(caplog):
    """The fallback must be as loud as a parameter sitting on its bound."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        analyze_oxide_layer(freq, Z, epsilon_r=22.0)

    assert 'NOT FROM THE FIT' in caplog.text
    assert 'no confidence interval' in caplog.text


def test_cc_still_wins_over_a_plain_capacitance():
    """C is CC's degenerate case (dC = 0); the general model must not lose."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - CC(CC_C_INF, CC_DC, CC_TAU, CC_ALPHA) - (R(R_P) | C(1e-5))

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0,
                                fit_result=_fit_result(circuit,
                                                       circuit.get_all_params()))

    assert oxide is not None
    assert oxide.element_type == 'CC'
