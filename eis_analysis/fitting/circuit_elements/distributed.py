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


class CC(CircuitElement):
    """
    Cole-Cole dielectric relaxation element (permittivity plane).

    C*(ω) = C_inf + ΔC / (1 + (jωτ)^(1-α))

    Z_CC = 1 / (jω · C*(ω))

    Models a dielectric with a *distribution* of relaxation times, as opposed
    to Q (a distribution of RC time constants seen in the impedance plane).
    The arc is depressed in the complex capacitance (equivalently permittivity)
    plane, which is where the physics of an oxide or polymer film lives.

    Parameters
    ----------
    C_inf : float or str, optional
        High-frequency limit capacitance [F] (default: 1e-8)
    dC : float or str, optional
        Relaxation strength ΔC = C_s - C_inf [F] (default: 1e-7),
        where C_s is the static (low-frequency) capacitance
    tau : float or str, optional
        Relaxation time [s] (default: 1e-3)
    alpha : float or str, optional
        Broadening exponent [dimensionless], 0 <= α < 1 (default: 0.2)

    Notes
    -----
    Special cases:
    - α = 0: Debye relaxation (a single relaxation time)
    - 0 < α < 1: the C* arc is depressed, its centre lies below the real
      axis by α·90°; larger α means a broader distribution

    Parametrised by ΔC rather than C_s so that ΔC > 0 (from the bounds)
    enforces C_s > C_inf on its own. With (C_inf, C_s) that inequality would
    be a coupling between two parameters, which box bounds cannot express.

    Capacitances, not permittivities: the element carries no geometry.
    Convert with ε_r = C·d/(ε₀·A) in the analysis layer, where the film
    thickness and electrode area are already known.

    Equivalent to a composite of existing elements (exact, verified to 2e-16):

        CC(C_inf, ΔC, τ, α) == C(C_inf) | (C(ΔC) - Q(ΔC/τ**(1-α), α))
        CC(C_inf, ΔC, τ, 0) == C(C_inf) | (C(ΔC) - R(τ/ΔC))       # Debye

    The composite form is not a substitute for fitting: its CPE coefficient
    ΔC/τ^(1-α) couples three parameters non-linearly, so a fit would report
    Q and n with confidence intervals on the wrong quantities, and ΔC would
    appear twice as two independent free parameters.

    Generalising to a second exponent gives Havriliak-Negami; see
    doc/LEVM_CIRCUITS.md (LEVM NDE = 6, 7).

    Examples
    --------
    >>> cc = CC(1e-8, 1e-7, 1e-3, 0.2)  # depressed dielectric arc
    >>> cc = CC()                       # default values
    >>> cc = CC(1e-8, 1e-7, 1e-3, 0.0)  # Debye limit
    """

    def __init__(self, C_inf: Union[float, str] = 1e-8,
                 dC: Union[float, str] = 1e-7,
                 tau: Union[float, str] = 1e-3,
                 alpha: Union[float, str] = 0.2):
        super().__init__(C_inf, dC, tau, alpha)
        self.C_inf = self.params[0]
        self.dC = self.params[1]
        self.tau = self.params[2]
        self.alpha = self.params[3]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        C_inf_val, dC_val = params[0], params[1]
        tau_val, alpha_val = params[2], params[3]
        omega = 2 * np.pi * freq
        C_star: NDArray[np.complex128] = (
            C_inf_val + dC_val / (1 + (1j * omega * tau_val) ** (1.0 - alpha_val))
        )
        return 1 / (1j * omega * C_star)

    def _scale(self, scalar: float) -> 'CC':
        """Scale both capacitances, preserving tau and alpha.

        Scaling a dielectric element means putting an identical one in
        parallel, i.e. scaling the electrode area: the capacitances are
        extensive, tau and alpha are intensive material properties.
        """
        new_C_inf = scalar * self.C_inf
        new_dC = scalar * self.dC
        # Preserve fixed status: the base class marks a parameter fixed when
        # it arrives as a str, so pass the number through repr, not through
        # quote characters (those would reach float() and raise).
        C_inf_arg = str(new_C_inf) if self.fixed_params[0] else new_C_inf
        dC_arg = str(new_dC) if self.fixed_params[1] else new_dC
        tau_arg = str(self.tau) if self.fixed_params[2] else self.tau
        alpha_arg = str(self.alpha) if self.fixed_params[3] else self.alpha
        return CC(C_inf_arg, dC_arg, tau_arg, alpha_arg)

    def get_param_labels(self) -> List[str]:
        return ['C_inf', 'ΔC', 'τ_CC', 'α_CC']

    def __repr__(self) -> str:
        C_inf_str = f'"{self.C_inf:.4g}"' if self.fixed_params[0] else f"{self.C_inf:.4g}"
        dC_str = f'"{self.dC:.4g}"' if self.fixed_params[1] else f"{self.dC:.4g}"
        tau_str = f'"{self.tau:.4g}"' if self.fixed_params[2] else f"{self.tau:.4g}"
        alpha_str = f'"{self.alpha:.4g}"' if self.fixed_params[3] else f"{self.alpha:.4g}"
        return f"CC(C_inf={C_inf_str}, ΔC={dC_str}, τ={tau_str}, α={alpha_str})"

    @property
    def C_static(self) -> float:
        """Static (low-frequency) capacitance C_s = C_inf + ΔC [F]"""
        return self.C_inf + self.dC

    @property
    def characteristic_freq(self) -> float:
        """Peak-loss frequency f = 1/(2*pi*tau) [Hz]"""
        return 1.0 / (2 * np.pi * self.tau)
