"""Test the finite Warburg elements Ws (tanh) and Wo (coth).

Up to and including v0.52.0 the tanh form was called Wo, under the
docstring "open" - the name ZView and impedance.py give to the coth form. The tests pin which
boundary each name stands for, by its low-frequency limit:

1. Ws: constant concentration at the far boundary -> Z -> R_W (resistive).
2. Wo: zero flux at the far boundary -> Z -> R_W/3 + R_W/(jw*tau) (capacitive).
3. Both reduce to semi-infinite W(R_W/sqrt(2*tau)) at high frequency.
4. Both fit back their own parameters from noise-free data written with
   cosh/sinh, independent of the element code.
"""

import numpy as np
import pytest

from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting import W, Wo, Ws, fit_equivalent_circuit

R_W, TAU = 100.0, 1.0
FREQ = np.logspace(5, -3, 129)   # 16 points/decade, as example/diffusion_*.csv
U = np.sqrt(2j * np.pi * FREQ * TAU)


def test_Ws_low_frequency_limit_is_the_diffusion_resistance():
    Z = Ws(R_W, TAU).impedance(np.array([1e-6]), [R_W, TAU])
    assert abs(Z[0] - R_W) < 1e-3   # Im ~ -R_W*w*tau/3 = 2e-4 Ohm at 1 uHz


def test_Wo_low_frequency_limit_is_resistance_third_plus_capacitance():
    f = np.array([1e-6])
    Z = Wo(R_W, TAU).impedance(f, [R_W, TAU])
    assert Z[0].real == pytest.approx(R_W / 3, rel=1e-6)
    assert Z[0].imag == pytest.approx(-R_W / (2 * np.pi * f[0] * TAU), rel=1e-6)


@pytest.mark.parametrize("element", [Ws, Wo])
def test_high_frequency_limit_is_semi_infinite_warburg(element):
    f = np.array([1e5])
    sigma = R_W / np.sqrt(2 * TAU)
    Z = element(R_W, TAU).impedance(f, [R_W, TAU])
    assert Z[0] == pytest.approx(W(sigma).impedance(f, [sigma])[0], rel=1e-9)


@pytest.mark.parametrize("name, Z_diffusion", [
    ("Ws", R_W * np.sinh(U) / (U * np.cosh(U))),
    ("Wo", R_W * np.cosh(U) / (U * np.sinh(U))),
])
def test_fit_recovers_parameters(name, Z_diffusion):
    # cosh/sinh overflow above |u| ~ 710 (|u| ~ 790 at 1e5 Hz with tau = 1 s),
    # so the reference data stop at 1e4 Hz
    keep = FREQ <= 1e4
    Z = 10.0 + Z_diffusion[keep]
    start = parse_circuit_expression(f"R(5) - {name}(50, 0.3)")
    result, _ = fit_equivalent_circuit(FREQ[keep], Z, start)
    assert result.params_opt == pytest.approx([10.0, R_W, TAU], rel=1e-4)


def test_Wo_survives_large_argument():
    """1/(u*tanh u) must stay finite where coth's cosh/sinh would overflow."""
    Z = Wo(R_W, 1e4).impedance(np.array([1e6]), [R_W, 1e4])
    assert np.all(np.isfinite(Z))


# Bounds of sigma and R_W span what R spans: a mOhm battery and a
# high-impedance oxide both fit without coming within the one decade
# classify_bound_status flags. (C stays a decade below its own 0.1 F bound,
# a separate limit.)
@pytest.mark.parametrize("truth, start", [
    ("R(0.002)-(R(0.005)|C(0.005))-W(0.0003)", "R(0.001)-(R(0.001)|C(0.001))-W(0.01)"),
    ("R(10)-(R(1e6)|C(1e-9))-W(1e6)", "R(5)-(R(1e5)|C(1e-8))-W(1e5)"),
    ("R(0.002)-(R(0.005)|C(0.002))-Ws(0.003,20)", "R(0.001)-(R(0.001)|C(0.001))-Ws(0.01,5)"),
    ("R(10)-(R(1e6)|C(1e-9))-Ws(3e8,20)", "R(5)-(R(1e5)|C(1e-8))-Ws(1e8,5)"),
])
def test_warburg_scale_covers_batteries_and_coatings(truth, start):
    tc = parse_circuit_expression(truth)
    f = np.logspace(5, -3, 81)
    Z = tc.impedance(f, tc.get_all_params())
    result, _ = fit_equivalent_circuit(f, Z, parse_circuit_expression(start))
    assert result.params_opt == pytest.approx(tc.get_all_params(), rel=1e-4)
    assert not any(result.bound_status)
