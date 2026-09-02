"""
Basic lumped circuit elements: resistor (R), conductance (G), capacitor (C),
inductor (L).
"""
from __future__ import annotations

import numpy as np
from typing import List, Union
from numpy.typing import NDArray

from .base import CircuitElement, param_property


class R(CircuitElement):
    """
    Resistor.

    Z_R = R

    Parameters
    ----------
    R : float or str, optional
        Resistance [Ω] (default: 100)
        If passed as string, the parameter is fixed during fitting.

    Examples
    --------
    >>> r = R(100)    # 100 Ω resistor (free parameter)
    >>> r = R()       # default 100 Ω (free)
    >>> r = R("0.86") # Fixed R_inf = 0.86 Ω
    """

    R = param_property(0)

    def __init__(self, R: Union[float, str] = 100.0):
        super().__init__(R)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_val = params[0]
        return R_val * np.ones_like(freq, dtype=complex)

    def get_param_labels(self) -> List[str]:
        return ['R']

    def __repr__(self) -> str:
        if self.fixed_params[0]:
            return f'R("{self.R:.4g}")'
        return f"R({self.R:.4g})"


class C(CircuitElement):
    """
    Capacitor.

    Z_C = 1 / (jωC)

    Parameters
    ----------
    C : float or str, optional
        Capacitance [F] (default: 1e-6)
        If passed as string, the parameter is fixed during fitting.

    Examples
    --------
    >>> c = C(1e-6)    # 1 μF capacitor (free parameter)
    >>> c = C()        # default 1 μF (free)
    >>> c = C("1e-6")  # Fixed capacitance
    """

    C = param_property(0)

    def __init__(self, C: Union[float, str] = 1e-6):
        super().__init__(C)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        C_val = params[0]
        omega = 2 * np.pi * freq
        return 1 / (1j * omega * C_val)

    def get_param_labels(self) -> List[str]:
        return ['C']

    def __repr__(self) -> str:
        if self.fixed_params[0]:
            return f'C("{self.C:.4g}")'
        return f"C({self.C:.4g})"


class L(CircuitElement):
    """
    Inductor.

    Z_L = jωL

    Parameters
    ----------
    L : float or str, optional
        Inductance [H] (default: 1e-6)
        If passed as string, the parameter is fixed during fitting.

    Examples
    --------
    >>> l = L(1e-6)    # 1 μH inductor (free parameter)
    >>> l = L()        # default 1 μH (free)
    >>> l = L("1e-6")  # Fixed inductance
    """

    L = param_property(0)

    def __init__(self, L: Union[float, str] = 1e-6):
        super().__init__(L)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        L_val = params[0]
        omega = 2 * np.pi * freq
        return 1j * omega * L_val

    def get_param_labels(self) -> List[str]:
        return ['L']

    def __repr__(self) -> str:
        if self.fixed_params[0]:
            return f'L("{self.L:.4g}")'
        return f"L({self.L:.4g})"


class G(CircuitElement):
    """
    Conductance - a resistor parametrized by its admittance.

    Y_G = G,  Z_G = 1 / G

    Electrically identical to ``R(1/G)``, but the two are not
    interchangeable for *fitting*. A parallel resistance the data cannot
    determine from above drives R towards its upper bound (10 GOhm), which
    is an artificial edge of the parameter box, not a feature of the model:
    the Jacobian column dZ/dR ~ 1/R^2 collapses there, the covariance
    becomes numerically meaningless, and the reported error is spurious
    precision rather than a large uncertainty. The same physical limit is
    G -> 0, a finite point where dZ/dG = -Z^2 stays well conditioned. The
    fit then returns G = (0 +- s) S, whose one-sided interval is a genuine
    lower bound on R.

    Use G in place of a parallel R when the resistance may be too large to
    resolve (blocking coatings, intact oxide layers). For a resistance the
    data determines, R is the more readable parametrization.

    Parameters
    ----------
    G : float or str, optional
        Conductance [S] (default: 1e-6, i.e. 1 MOhm)
        If passed as string, the parameter is fixed during fitting.

    Notes
    -----
    G = 0 is inside the allowed range, so the impedance is clipped at
    G_MIN to keep 1/G finite. The clip is invisible physically (1e-30 S is
    1e30 Ohm) and exact in the parallel chain rule: the Jacobian factor
    (Z_total/Z_G)^2 * dZ_G/dG evaluates to -Z_total^2, the correct
    derivative of Z = 1/(G + Y_rest).

    Examples
    --------
    >>> g = G(1e-9)     # 1 nS = 1 GOhm (free parameter)
    >>> g = G()         # default 1 uS (free)
    >>> g = G("0")      # fixed open circuit
    >>> circuit = R(10) - (G(1e-9) | Q(1e-7, 0.9))
    """

    # ponytail: 1e-30 S (1e30 Ohm) stands in for G = 0 so that 1/G never
    # overflows to inf - complex 1/(0+0j) is inf+nanj in numpy, and the NaN
    # would poison Parallel.impedance. Ceiling: a conductance genuinely
    # below 1e-30 S is indistinguishable from zero here. Upgrade path is an
    # admittance-native Parallel (elements expose Y, no reciprocal taken).
    G_MIN = 1e-30

    G = param_property(0)

    def __init__(self, G: Union[float, str] = 1e-6):
        super().__init__(G)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        G_val = max(params[0], self.G_MIN)
        return (1 / G_val) * np.ones_like(freq, dtype=complex)

    def get_param_labels(self) -> List[str]:
        return ['G']

    @property
    def resistance(self) -> float:
        """Equivalent resistance R = 1/G [Ohm]; inf for G = 0."""
        return 1 / self.G if self.G > 0 else np.inf

    def __repr__(self) -> str:
        if self.fixed_params[0]:
            return f'G("{self.G:.4g}")'
        return f"G({self.G:.4g})"
