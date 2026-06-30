"""
Distributed circuit elements: constant phase element (Q) and Warburg
diffusion (semi-infinite W, finite/open Wo).
"""
from __future__ import annotations

import numpy as np
from typing import List, Union
from numpy.typing import NDArray

from .base import CircuitElement


class Q(CircuitElement):
    """
    Constant Phase Element (CPE).

    Z_Q = 1 / (Q * (jω)^n)

    Parameters
    ----------
    Q : float, optional
        CPE coefficient [F·s^(n-1)] (default: 1e-4)
    n : float, optional
        CPE exponent [dimensionless], 0 < n < 1 (default: 0.8)

    Notes
    -----
    Special cases:
    - n = 1: ideal capacitor
    - n = 0.5: Warburg diffusion
    - 0.5 < n < 1: distributed relaxation

    Examples
    --------
    >>> q = Q(1e-4, 0.8)  # typical CPE
    >>> q = Q()           # default values
    >>> q = Q(1e-4)**0.9  # change exponent
    """

    def __init__(self, Q_val: Union[float, str] = 1e-4, n: Union[float, str] = 0.8):
        super().__init__(Q_val, n)
        self.Q = self.params[0]
        self.n = self.params[1]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        Q_val, n_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return 1 / (Q_val * (1j * omega) ** n_val)

    def _scale(self, scalar: float) -> 'Q':
        return Q(scalar * self.Q, self.n)

    def __pow__(self, exponent: float) -> 'Q':
        """Change CPE exponent: Q()**0.9"""
        return Q(self.Q, exponent)

    def get_param_labels(self) -> List[str]:
        return ['Q', 'n']

    def __repr__(self) -> str:
        Q_str = f'"{self.Q:.4g}"' if self.fixed_params[0] else f"{self.Q:.4g}"
        n_str = f'"{self.n:.4g}"' if self.fixed_params[1] else f"{self.n:.4g}"
        return f"Q({Q_str}, {n_str})"


class W(CircuitElement):
    """
    Warburg semi-infinite diffusion element.

    Z_W = σ/√ω * (1 - j) = σ(1-j)/√ω

    Parameters
    ----------
    sigma : float or str, optional
        Warburg coefficient [Ω·s^(-1/2)] (default: 50.0)
        If passed as string, the parameter is fixed during fitting.

    Examples
    --------
    >>> w = W(50)     # Warburg with σ=50 (free parameter)
    >>> w = W()       # default σ=50 (free)
    >>> w = W("50")   # Fixed σ=50
    """

    def __init__(self, sigma: Union[float, str] = 50.0):
        super().__init__(sigma)
        self.sigma = self.params[0]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        sigma_val = params[0]
        omega = 2 * np.pi * freq
        return sigma_val / np.sqrt(omega) * (1 - 1j)

    def _scale(self, scalar: float) -> 'W':
        return W(scalar * self.sigma)

    def get_param_labels(self) -> List[str]:
        return ['σ']

    def __repr__(self) -> str:
        sigma_str = f'"{self.sigma:.4g}"' if self.fixed_params[0] else f"{self.sigma:.4g}"
        return f"W(σ={sigma_str})"


class Wo(CircuitElement):
    """
    Warburg open (bounded) diffusion element.

    Z_Wo = R_W * tanh(√(jωτ_W)) / √(jωτ_W)

    Parameters
    ----------
    R_W : float or str, optional
        Warburg resistance [Ω] (default: 100.0)
        If passed as string, the parameter is fixed during fitting.
    tau_W : float or str, optional
        Diffusion time constant [s] (default: 1.0)
        If passed as string, the parameter is fixed during fitting.

    Examples
    --------
    >>> wo = Wo(100, 1.0)      # Both parameters free
    >>> wo = Wo()              # default values (both free)
    >>> wo = Wo("100", 1.0)    # R_W fixed, tau_W free
    >>> wo = Wo("100", "1.0")  # Both parameters fixed
    """

    def __init__(self, R_W: Union[float, str] = 100.0, tau_W: Union[float, str] = 1.0):
        super().__init__(R_W, tau_W)
        self.R_W = self.params[0]
        self.tau_W = self.params[1]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_W_val, tau_W_val = params[0], params[1]
        omega = 2 * np.pi * freq
        arg = np.sqrt(1j * omega * tau_W_val)
        return R_W_val * np.tanh(arg) / arg

    def _scale(self, scalar: float) -> 'Wo':
        return Wo(scalar * self.R_W, self.tau_W)

    def get_param_labels(self) -> List[str]:
        return ['R_W', 'τ_W']

    def __repr__(self) -> str:
        R_W_str = f'"{self.R_W:.4g}"' if self.fixed_params[0] else f"{self.R_W:.4g}"
        tau_W_str = f'"{self.tau_W:.4g}"' if self.fixed_params[1] else f"{self.tau_W:.4g}"
        return f"Wo(R={R_W_str}, τ={tau_W_str})"
