#!/usr/bin/env python3
"""Parameter significance: S_i = max_n |dln|Z_n|/dln P_i|.

Two things are worth pinning down, and they are not the same thing:

1. that the closed-form expression is the derivative it claims to be
   (checked against central differences, like tests/test_jacobian.py), and
2. that the number lands on the scale its interpretation assumes -
   ~1 for a dominating element, below SIGNIFICANCE_NEGLIGIBLE for one that
   can be dropped. A formula can be correct and still be useless if the
   threshold in config.py does not mean what it says.
"""

import numpy as np
import pytest

from eis_analysis.fitting import R, C, Q, W
from eis_analysis.fitting.config import SIGNIFICANCE_NEGLIGIBLE
from eis_analysis.fitting.diagnostics import compute_significance


@pytest.fixture
def freq():
    return np.logspace(5, -2, 60)


def numerical_significance(circuit, freq, params, eps=1e-6):
    """S by central differences in log-log space, independent of the analytic path."""
    params = list(params)
    S = np.zeros(len(params))

    for j in range(len(params)):
        if params[j] == 0:
            continue  # dln P undefined; the analytic form returns 0 by construction
        p_fwd, p_bwd = params.copy(), params.copy()
        p_fwd[j] = params[j] * (1 + eps)
        p_bwd[j] = params[j] * (1 - eps)
        log_Z_fwd = np.log(np.abs(circuit.impedance(freq, p_fwd)))
        log_Z_bwd = np.log(np.abs(circuit.impedance(freq, p_bwd)))
        d_log_P = np.log1p(eps) - np.log1p(-eps)
        S[j] = np.max(np.abs((log_Z_fwd - log_Z_bwd) / d_log_P))

    return S


# ---------------------------------------------------------------------------
# 1. Correctness
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("circuit", [
    R(10) - (R(1000) | C(1e-6)),
    R(10) - (R(1000) | Q(1e-6, 0.8)) - W(50),
    (R(10) | C(1e-9)) - (R(500) | Q(2e-5, 0.7)) - (R(2000) | C(1e-4)),
])
def test_matches_central_differences(circuit, freq):
    """The closed form is the log-log derivative it claims to be."""
    params = circuit.get_all_params()
    S = compute_significance(circuit, freq, params)

    assert S is not None
    assert np.allclose(S, numerical_significance(circuit, freq, params), rtol=1e-4)


def test_one_value_per_parameter(freq):
    """S is aligned with the parameter vector, not with elements."""
    # Q contributes two parameters, so five elements give six values.
    circuit = R(10) - (R(1000) | Q(1e-6, 0.8))
    S = compute_significance(circuit, freq, circuit.get_all_params())

    assert len(S) == len(circuit.get_all_params()) == 4
    assert np.all(np.isfinite(S))
    assert np.all(S >= 0)


# ---------------------------------------------------------------------------
# 2. Calibration - the interpretation the threshold rests on
# ---------------------------------------------------------------------------

def test_negligible_element_falls_below_threshold(freq):
    """Zahner's own example: a resistor too small to matter scores ~0.

    Their results window shows a second resistance of 338 nOhm at S = 0.002
    next to an arc of ordinary size, and concludes it can be removed.
    """
    circuit = R(10) - (R(1000) | C(1e-6)) - R(338e-9)
    S = compute_significance(circuit, freq, circuit.get_all_params())

    R_negligible = S[-1]
    assert R_negligible < SIGNIFICANCE_NEGLIGIBLE
    # Same order of magnitude as the value the manual reports for this case.
    assert R_negligible < 1e-3


def test_dominating_resistor_approaches_one(freq):
    """A resistor that owns part of the window scores ~1, and never above it.

    For an element entering linearly, S is bounded by 1: with dZ/dR = 1 the
    expression reduces to R*cos(phi)/|Z|, the fraction of |Z| it accounts for.
    """
    # R_s dominates at high frequency, R_ct at low - each owns one end.
    circuit = R(50) - (R(5000) | C(1e-6))
    S = compute_significance(circuit, freq, circuit.get_all_params())

    R_s, R_ct = S[0], S[1]
    assert R_s == pytest.approx(1.0, abs=0.05)
    assert R_ct == pytest.approx(1.0, abs=0.05)
    assert R_s <= 1.0 + 1e-9 and R_ct <= 1.0 + 1e-9


def test_cpe_exponent_may_exceed_one(freq):
    """A non-linearly entering parameter is allowed past 1 - it is not a bug.

    d ln|Z|/d ln(alpha) = -alpha*ln(omega/omega_0) grows without bound away
    from the normalisation frequency, so the "fraction of |Z|" reading that
    bounds a resistor does not apply here.
    """
    circuit = R(10) - (R(1000) | Q(1e-6, 0.8))
    S = compute_significance(circuit, freq, circuit.get_all_params())

    assert S[3] > 1.0  # the exponent n of Q


def test_scale_invariant(freq):
    """S is dimensionless: scaling the whole circuit leaves it unchanged.

    This is what makes R [Ohm] and C [F] comparable on one axis.
    """
    circuit = R(10) - (R(1000) | C(1e-6))
    params = circuit.get_all_params()
    # Scale every impedance by 1000: R -> 1000R, C -> C/1000.
    scaled = [params[0] * 1000, params[1] * 1000, params[2] / 1000]

    S = compute_significance(circuit, freq, params)
    S_scaled = compute_significance(circuit, freq, scaled)

    assert np.allclose(S, S_scaled, rtol=1e-9)


def test_zero_parameter_scores_zero(freq):
    """A parameter that is exactly zero cannot influence the impedance."""
    circuit = R(10) - (R(1000) | C(1e-6)) - R(0.0)
    S = compute_significance(circuit, freq, circuit.get_all_params())

    assert S[-1] == 0.0
