#!/usr/bin/env python3
"""Named parameter attributes on circuit elements.

Regression tests (2026-08-27): every element copied its parameters into
named attributes in __init__ (self.R = self.params[0], ...). update_params()
- which the fitter calls with the optimized values - rewrote params but left
those copies untouched, so K.R, GE.sigma, CC.C_inf and the properties derived
from them (K.capacitance, GE.characteristic_freq, CC.C_static) still returned
the initial guess after a fit. It surfaced as --analyze-oxide reporting a
thickness computed from the guess rather than the fit. The names are now
read-only views of params, so they cannot fall out of step.
"""

import numpy as np
import pytest
from eis_analysis.fitting import R, C, L, Q, W, Wo, K, GE, CC


# factory, attribute names in parameter order, fitted values to push in.
# Factories, not instances: update_params() mutates, and a shared instance
# would leak state from one element's check into the next.
ELEMENTS = {
    "R": (lambda: R(100), ["R"], [333.0]),
    "C": (lambda: C(1e-6), ["C"], [4.2e-6]),
    "L": (lambda: L(1e-6), ["L"], [7.5e-6]),
    "Q": (lambda: Q(1e-4, 0.8), ["Q", "n"], [3e-4, 0.65]),
    "W": (lambda: W(50), ["sigma"], [125.0]),
    "Wo": (lambda: Wo(100, 1.0), ["R_W", "tau_W"], [250.0, 2.5]),
    "K": (lambda: K(1000, 1e-4), ["R", "tau"], [2000.0, 5e-4]),
    "GE": (lambda: GE(100, 1e-3), ["sigma", "tau"], [200.0, 5e-3]),
    "CC": (lambda: CC(1e-8, 1e-7, 1e-3, 0.2), ["C_inf", "dC", "tau", "alpha"],
           [2e-8, 2e-7, 5e-3, 0.3]),
}


def test_named_attributes_track_update_params():
    """After update_params(), every named attribute reports the new value."""
    for label, (make, names, fitted) in ELEMENTS.items():
        element = make()
        element.update_params(fitted)

        for name, expected in zip(names, fitted):
            assert getattr(element, name) == pytest.approx(expected), f"{label}.{name}"
        assert element.get_all_params() == pytest.approx(fitted), label


def test_named_attributes_are_read_only():
    """params is the single source of truth; the names cannot be assigned."""
    for label, (make, names, fitted) in ELEMENTS.items():
        with pytest.raises(AttributeError):
            setattr(make(), names[0], fitted[0])
            pytest.fail(f"{label}.{names[0]} is assignable")


def test_repr_reflects_updated_params():
    """__repr__ reads the attributes, so it must show fitted values too."""
    for label, (make, names, fitted) in ELEMENTS.items():
        element = make()
        before = repr(element)
        element.update_params(fitted)
        assert repr(element) != before, label


def test_derived_properties_track_update_params():
    """The properties computed from the names follow as well."""
    k = K(1000, 1e-4)
    k.update_params([2000.0, 5e-4])
    assert k.capacitance == pytest.approx(5e-4 / 2000.0)
    assert k.characteristic_freq == pytest.approx(1 / (2 * np.pi * 5e-4))

    g = GE(100, 1e-3)
    g.update_params([200.0, 5e-3])
    assert g.characteristic_freq == pytest.approx(1 / (2 * np.pi * 5e-3))

    cc = CC(1e-8, 1e-7, 1e-3, 0.2)
    cc.update_params([2e-8, 2e-7, 5e-3, 0.3])
    assert cc.C_static == pytest.approx(2e-8 + 2e-7)
    assert cc.characteristic_freq == pytest.approx(1 / (2 * np.pi * 5e-3))


def test_impedance_after_update_matches_the_new_params():
    """The end-to-end consequence: a refitted element computes new impedance."""
    freq = np.logspace(4, -1, 20)
    k = K(1000, 1e-4)
    k.update_params([2000.0, 5e-4])

    Z_actual = k.impedance(freq, k.get_all_params())
    Z_expected = 2000.0 / (1 + 1j * 2 * np.pi * freq * 5e-4)
    assert np.max(np.abs(Z_actual - Z_expected)) < 1e-12
