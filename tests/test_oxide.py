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

from eis_analysis.analysis.config import EPSILON_0
from eis_analysis.analysis.oxide import analyze_oxide_layer, estimate_permittivity
from eis_analysis.fitting.circuit import FitResult
from eis_analysis.fitting.circuit_elements import C, K, Q, R

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
        eps_r = estimate_permittivity(
            freq, Z, thickness_nm=19.5, fit_result=_fit_result_voigt()
        )

    assert eps_r is not None
    assert 'Oxide thickness' not in caplog.text
    assert 'Permittivity' in caplog.text


def test_estimate_permittivity_fallback_does_not_log_thickness(caplog):
    """Regression (audit O2): same for the high-frequency fallback path."""
    freq, Z = _synthetic_voigt()

    with caplog.at_level(logging.INFO, logger=OXIDE_LOGGER):
        eps_r = estimate_permittivity(freq, Z, thickness_nm=19.5)

    assert eps_r is not None
    assert 'Oxide thickness' not in caplog.text


def test_permittivity_thickness_roundtrip():
    """estimate_permittivity() is the exact inverse of analyze_oxide_layer()."""
    freq, Z = _synthetic_voigt()
    fit_result = _fit_result_voigt()

    oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)
    eps_r = estimate_permittivity(
        freq, Z, thickness_nm=oxide.thickness_nm, fit_result=fit_result
    )

    assert abs(eps_r - 22.0) / 22.0 < 1e-9


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
    """Regression (audit O4): (R|C1|C2) warns instead of silently taking the last C."""
    freq, Z = _synthetic_voigt()
    circuit = R(R_S) - (R(R_P) | C(1e-5) | C(C_P))
    fit_result = _fit_result(circuit, [R_S, R_P, 1e-5, C_P])

    with caplog.at_level(logging.WARNING, logger=OXIDE_LOGGER):
        oxide = analyze_oxide_layer(freq, Z, epsilon_r=22.0, fit_result=fit_result)

    assert oxide is not None
    assert abs(oxide.capacitance - C_P) / C_P < 1e-9  # last one wins
    assert 'Multiple C/Q elements' in caplog.text


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
