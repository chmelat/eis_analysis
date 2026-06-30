"""
Basic lumped circuit elements: resistor (R), capacitor (C), inductor (L).
"""
from __future__ import annotations

import numpy as np
from typing import List, Union
from numpy.typing import NDArray

from .base import CircuitElement


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

    def __init__(self, R: Union[float, str] = 100.0):
        super().__init__(R)
        self.R = self.params[0]  # Store as float

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_val = params[0]
        return R_val * np.ones_like(freq, dtype=complex)

    def _scale(self, scalar: float) -> 'R':
        return R(scalar * self.R)

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

    def __init__(self, C: Union[float, str] = 1e-6):
        super().__init__(C)
        self.C = self.params[0]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        C_val = params[0]
        omega = 2 * np.pi * freq
        return 1 / (1j * omega * C_val)

    def _scale(self, scalar: float) -> 'C':
        return C(scalar * self.C)

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

    def __init__(self, L: Union[float, str] = 1e-6):
        super().__init__(L)
        self.L = self.params[0]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        L_val = params[0]
        omega = 2 * np.pi * freq
        return 1j * omega * L_val

    def _scale(self, scalar: float) -> 'L':
        return L(scalar * self.L)

    def get_param_labels(self) -> List[str]:
        return ['L']

    def __repr__(self) -> str:
        if self.fixed_params[0]:
            return f'L("{self.L:.4g}")'
        return f"L({self.L:.4g})"
