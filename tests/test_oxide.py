"""
Tests for oxide layer analysis (analysis/oxide.py).

Regression tests for audit finding O2 (2026-07-02): estimate_permittivity()
must not log a bogus "Oxide thickness" computed from a dummy epsilon_r.
"""

import logging

import numpy as np

from eis_analysis.analysis.config import EPSILON_0
from eis_analysis.analysis.oxide import analyze_oxide_layer, estimate_permittivity
from eis_analysis.fitting.circuit import FitResult
from eis_analysis.fitting.circuit_elements import C, R

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
