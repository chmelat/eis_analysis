"""
Distributed circuit elements: constant phase element (Q), Warburg
diffusion (semi-infinite W, finite/open Wo), Cole-Cole relaxation (CC)
and the bounded power-law distribution (DQ).
"""
from __future__ import annotations

import numpy as np
from typing import List, Tuple, Union
from numpy.typing import NDArray

from .base import CircuitElement, param_property


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
    """

    Q = param_property(0)
    n = param_property(1)

    def __init__(self, Q_val: Union[float, str] = 1e-4, n: Union[float, str] = 0.8):
        super().__init__(Q_val, n)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        Q_val, n_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return 1 / (Q_val * (1j * omega) ** n_val)

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

    sigma = param_property(0)

    def __init__(self, sigma: Union[float, str] = 50.0):
        super().__init__(sigma)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        sigma_val = params[0]
        omega = 2 * np.pi * freq
        return sigma_val / np.sqrt(omega) * (1 - 1j)

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

    R_W = param_property(0)
    tau_W = param_property(1)

    def __init__(self, R_W: Union[float, str] = 100.0, tau_W: Union[float, str] = 1.0):
        super().__init__(R_W, tau_W)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_W_val, tau_W_val = params[0], params[1]
        omega = 2 * np.pi * freq
        arg = np.sqrt(1j * omega * tau_W_val)
        return R_W_val * np.tanh(arg) / arg

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

    C_inf = param_property(0)
    dC = param_property(1)
    tau = param_property(2)
    alpha = param_property(3)

    def __init__(self, C_inf: Union[float, str] = 1e-8,
                 dC: Union[float, str] = 1e-7,
                 tau: Union[float, str] = 1e-3,
                 alpha: Union[float, str] = 0.2):
        super().__init__(C_inf, dC, tau, alpha)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        C_inf_val, dC_val = params[0], params[1]
        tau_val, alpha_val = params[2], params[3]
        omega = 2 * np.pi * freq
        C_star: NDArray[np.complex128] = (
            C_inf_val + dC_val / (1 + (1j * omega * tau_val) ** (1.0 - alpha_val))
        )
        return 1 / (1j * omega * C_star)

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


# Gauss-Legendre nodes for the DQ integral over ln tau. Measured, not guessed:
# what sets the accuracy is not the width in decades but the node spacing
# relative to the kernel's transition, which is ~1 neper wide. Worst-case
# relative error over a 1 MHz - 1 mHz window, against N=2048:
#
#     U (nepers) | N=48    N=64    N=96    N=128
#     -----------|--------------------------------
#     10         | 3e-13   3e-13   3e-13   3e-13
#     20         | 5e-7    3e-9    5e-13   5e-13
#     30         | 8e-5    3e-6    3e-9    4e-12
#     60         | 8e-3    1e-3    4e-5    1e-6
#
# N=96 holds 3e-9 across the whole allowed range of U_DQ (<= 30, see
# PARAMETER_BOUNDS), three orders below the 1e-3 residual structure this
# element exists to explain; N=64 would only reach 3e-6 there. The extra cost
# is 45 us per model call, ~5 s over a full DE run.
#
# The midpoint rule - equivalently, a chain of N discrete RC elements - is only
# O(h^2) and at a realistic N leaves a smooth 1e-4..1e-3 deviation, i.e. an
# artefact the size of the effect. Hence quadrature, not discretisation.
DQ_QUAD_NODES = 96

_GL_X, _GL_W = np.polynomial.legendre.leggauss(DQ_QUAD_NODES)


def dq_quadrature(omega: NDArray[np.float64], n: float, tau_min: float,
                  U: float) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """Gauss-Legendre nodes in ln(tau) and the DQ integrand at them.

    Shared by ``DQ.impedance`` and the analytic Jacobian so the two can never
    integrate over different nodes.

    Returns
    -------
    s : ndarray of float, shape (DQ_QUAD_NODES,)
        Nodes in s = ln(tau), spanning [ln(tau_min), ln(tau_min) + U]
    integrand : ndarray of complex, shape (n_freq, DQ_QUAD_NODES)
        e^(n*s) / (1 + j*omega*e^s)
    """
    s = np.log(tau_min) + 0.5 * U * (_GL_X + 1.0)
    integrand: NDArray[np.complex128] = (
        np.exp(n * s) / (1.0 + 1j * omega[:, None] * np.exp(s))
    )
    return s, integrand


class DQ(CircuitElement):
    """
    Bounded power-law distribution of relaxation times (truncated CPE).

    γ(τ) = A·τⁿ for τ_min <= τ <= τ_max, zero outside

    Z_DQ(ω) = ∫ γ(τ)/(1 + jωτ) d ln τ

    An ideal CPE is the same power law with no bounds, and that is exactly
    what makes it unphysical: γ(τ) = A·τⁿ is not normalisable, so the CPE has
    no distribution of relaxation times to recover (a DRT of a CPE-containing
    fit does not converge) and no DC limit, which forces a separate parallel
    conductance into the model. Truncating the distribution fixes both. The
    element has three regimes in one:

    - below 1/τ_max: a finite polarisation resistance R_pol
    - between 1/τ_max and 1/τ_min: CPE behaviour, slope -n
    - above 1/τ_min: capacitive, C_eff

    So one DQ does the work of G | Q | C - but with R_pol and C_eff *derived*
    from the distribution rather than independent, which is what removes the
    degeneracy that makes a fitted C drift with the frequency window.

    Parameters
    ----------
    A : float or str, optional
        Distribution amplitude [Ω·s^-n] (default: 1e-3). For an unbounded
        distribution A = sin(πn)/(π·Q), where Q is the CPE coefficient.
    n : float or str, optional
        Power-law exponent [dimensionless] (default: 0.6), the same exponent
        the CPE reports: inside the bounds |Z| has slope -n.
    tau_min : float or str, optional
        Lower bound of the distribution [s] (default: 1e-6)
    U : float or str, optional
        Log-width of the distribution, U = ln(τ_max/τ_min) [dimensionless]
        (default: 10.0)

    Notes
    -----
    Parametrised by (τ_min, U) rather than (τ_min, τ_max) so that U > 0 keeps
    the bounds ordered by construction - box bounds cannot express
    τ_max > τ_min, and a fitter handed the pair directly will cross them.
    τ_max is available as a derived property.

    A bound that lies outside the measured window is still identifiable, but
    the fit leans on the power law to place it: expect U to correlate with n
    (measured: -0.85). Read n together with the bound status of U, not alone.

    Concept from LEVM (Macdonald), where the same truncation appears as the
    DWC models with limits U1, U2.

    Examples
    --------
    >>> dq = DQ(1e-3, 0.6, 1e-6, 10.0)  # 10 nepers ~ 4.3 decades wide
    >>> dq = DQ()                       # default values
    >>> dq = DQ(1e-3, "0.6", 1e-6, 10)  # exponent fixed, rest free
    """

    A = param_property(0)
    n = param_property(1)
    tau_min = param_property(2)
    U = param_property(3)

    def __init__(self, A: Union[float, str] = 1e-3,
                 n: Union[float, str] = 0.6,
                 tau_min: Union[float, str] = 1e-6,
                 U: Union[float, str] = 10.0):
        super().__init__(A, n, tau_min, U)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        A_val, n_val = params[0], params[1]
        tau_min_val, U_val = params[2], params[3]
        omega = 2 * np.pi * freq
        _, integrand = dq_quadrature(omega, n_val, tau_min_val, U_val)
        # 0.5*U is the Jacobian of the map from [-1, 1] to [s_min, s_max]
        Z: NDArray[np.complex128] = A_val * 0.5 * U_val * (integrand @ _GL_W)
        return Z

    def get_param_labels(self) -> List[str]:
        return ['A_DQ', 'n_DQ', 'τ_DQ', 'U_DQ']

    def __repr__(self) -> str:
        A_str = f'"{self.A:.4g}"' if self.fixed_params[0] else f"{self.A:.4g}"
        n_str = f'"{self.n:.4g}"' if self.fixed_params[1] else f"{self.n:.4g}"
        tau_str = f'"{self.tau_min:.4g}"' if self.fixed_params[2] else f"{self.tau_min:.4g}"
        U_str = f'"{self.U:.4g}"' if self.fixed_params[3] else f"{self.U:.4g}"
        return f"DQ(A={A_str}, n={n_str}, τ_min={tau_str}, U={U_str})"

    @property
    def tau_max(self) -> float:
        """Upper bound of the distribution τ_max = τ_min·e^U [s]"""
        return self.tau_min * np.exp(self.U)

    @property
    def R_pol(self) -> float:
        """DC limit Z(0) = A·(τ_max^n - τ_min^n)/n [Ω], A·U at n = 0

        n = 0 is outside PARAMETER_BOUNDS, but bounds only constrain what is
        fitted: a parameter fixed as a string enters the fit as given (see
        validate_fixed_params), and DQ(..., n="0", ...) is a legitimate flat
        distribution. Without this branch R_pol is 0/0 = nan, which then
        fails the "has a parallel resistance" test silently.
        """
        if self.n == 0:
            return self.A * self.U
        return self.A * (self.tau_max ** self.n - self.tau_min ** self.n) / self.n

    @property
    def C_eff(self) -> float:
        """High-frequency capacitance, from Z -> 1/(jω·C_eff) above 1/τ_min [F]

        C_eff = (1-n)/(A·(τ_min^(n-1) - τ_max^(n-1)))
        """
        if self.n == 1:
            return 1.0 / (self.A * self.U)
        return (1.0 - self.n) / (
            self.A * (self.tau_min ** (self.n - 1) - self.tau_max ** (self.n - 1))
        )
