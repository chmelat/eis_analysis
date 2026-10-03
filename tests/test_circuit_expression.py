#!/usr/bin/env python3
"""Tests for parse_circuit_expression() (cli/utils.py).

Regression tests (2026-08-27, plan 002): the element registry is duplicated
between circuit_elements/__init__.py's __all__ and the hand-written
safe_namespace dict in cli/utils.py, so adding an element requires updating
both. GE drifted -- it was fully implemented, had an analytic Jacobian, its
own test file and a row in the README table, but was missing from the
namespace, so --circuit "GE(100,1e-3)" failed with "name 'GE' is not defined"
for eight months. The parametrized table below is the guard: it must list
every element the README documents.
"""

import numpy as np
import pytest

import eis_analysis
from eis_analysis.cli.utils import parse_circuit_expression


# Exactly as the README's element table writes them.
DOCUMENTED_ELEMENTS = [
    "R(100)",
    "C(1e-6)",
    "L(1e-6)",
    "Q(1e-4, 0.8)",
    "W(50)",
    "Ws(100, 1.0)",
    "Wo(100, 1.0)",
    "K(1000, 1e-4)",
    "GE(100, 1e-3)",
    "CC(1e-8, 1e-7, 1e-3, 0.2)",
]
NAMES = ["R", "C", "L", "Q", "W", "Ws", "Wo", "K", "GE", "CC"]


def test_all_documented_elements_parse():
    """Every element the README documents is reachable from --circuit."""
    for expr in DOCUMENTED_ELEMENTS:
        circuit = parse_circuit_expression(expr)
        assert circuit is not None, expr
        assert len(circuit.get_all_params()) > 0, expr


def test_all_documented_elements_are_exported_top_level():
    """`from eis_analysis import <element>` works for each of them."""
    for name in NAMES:
        assert hasattr(eis_analysis, name), f"eis_analysis.{name} missing"
        assert name in eis_analysis.__all__, f"{name} missing from __all__"


def test_gerischer_in_series_circuit():
    """GE composes and computes: R_s + GE has DC limit R_s + sigma."""
    circuit = parse_circuit_expression("R(10) - GE(100, 1e-3)")
    assert len(circuit.get_all_params()) == 3

    freq = np.logspace(4, -1, 40)
    Z = circuit.impedance(freq, circuit.get_all_params())
    assert np.all(np.isfinite(Z))
    # Z_G(0) = sigma, so Re(Z) -> R_s + sigma = 110 Ohm
    assert abs(Z[-1].real - 110.0) / 110.0 < 0.01


def test_series_and_parallel_composition():
    """The two documented operators still build the expected parameter list."""
    circuit = parse_circuit_expression("R(100) - (R(5000) | C(1e-6))")
    assert circuit.get_all_params() == pytest.approx([100.0, 5000.0, 1e-6])


def test_fixed_parameter_syntax():
    """Quoted values are fixed during fitting (README convention)."""
    circuit = parse_circuit_expression('R("100") - C(1e-6)')
    assert circuit.get_all_fixed_params() == [True, False]


def test_invalid_expression_raises_valueerror():
    """The documented contract: unparseable input raises ValueError."""
    with pytest.raises(ValueError):
        parse_circuit_expression("NotAnElement(1)")


def test_builtins_are_not_reachable():
    """eval() runs with no builtins, so the namespace is elements only."""
    with pytest.raises(ValueError):
        parse_circuit_expression("__import__('os').getcwd()")


@pytest.mark.parametrize("expr", DOCUMENTED_ELEMENTS + [
    'G(1e-9)', 'DQ(1e-3, 0.6, 1e-6, 10)', 'YG(1e-6, 0.1, 1e-3)',
    'R(10) - (R("1000") | Q(1e-5, "0.9")) - YG("1e-6", 0.1, 1e-3)',
    'R("1234.56") - Q(1.23456789e-5, "0.87654321")',
])
def test_repr_parses_back(expr):
    """repr(circuit) is a circuit expression: same elements, values and fixes.

    The fit export stores str(circuit) for the Voigt chain; named arguments
    such as K(R=..., τ=...) used to make that unparseable.
    """
    circuit = parse_circuit_expression(expr)
    again = parse_circuit_expression(repr(circuit))
    assert again.get_all_fixed_params() == circuit.get_all_fixed_params()
    fixed = circuit.get_all_fixed_params()
    for value, original, is_fixed in zip(again.get_all_params(),
                                         circuit.get_all_params(), fixed):
        if is_fixed:
            # Never refitted, so it must come back exactly
            assert value == original
        else:
            # Only a starting point; repr rounds it to 4 digits
            assert value == pytest.approx(original, rel=1e-3)


def test_repr_of_fixed_numpy_value_parses_back():
    """After a fit the values are np.float64; numpy >= 2 reprs those as
    'np.float64(1000.0)', which must not leak into the expression."""
    class Float64Like(float):
        def __repr__(self):
            return f"np.float64({float(self)!r})"

    circuit = parse_circuit_expression('R("1000")')
    circuit.update_params([Float64Like(1000.0)])
    assert repr(circuit) == 'R("1000.0")'
