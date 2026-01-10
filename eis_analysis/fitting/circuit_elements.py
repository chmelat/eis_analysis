"""
Circuit elements for EIS analysis using operator overloading.

This module replaces the old string parser from impedance.py with a modern
operator overloading approach inspired by EISAnalysis.jl (Julia).

Usage:
    from eis_analysis.fitting.circuit_elements import R, C, Q, W

    # Build circuits using operators
    circuit = R(100) - (R(5000) | C(1e-6))  # Voigt element
    randles = R(10) - (R(100) - W(50)) | Q(1e-4, 0.8)

    # Values serve as initial guess for fitting
    params = circuit.get_all_params()  # [100, 5000, 1e-6]

Operators:
    -  : Series connection (Z_total = Z1 + Z2)
    |  : Parallel connection (1/Z_total = 1/Z1 + 1/Z2)
    *  : Scale parameter (2*R(100) = R(200))
    ** : Exponent for Q/CPE (Q()**0.9)

Author: EIS Analysis Toolkit
Date: 2025-12-14
"""
from __future__ import annotations

import numpy as np
from typing import List, Union, TYPE_CHECKING
from abc import ABC, abstractmethod
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .circuit_builder import Series, Parallel
    Circuit = Union[Series, Parallel, 'CircuitElement']


# ============================================================================
# Base Class
# ============================================================================

class CircuitElement(ABC):
    """
    Abstract base class for all circuit elements.

    Attributes
    ----------
    params : list of float
        Element parameters (serve as initial guess for fitting)
    fixed_params : list of bool
        Whether each parameter is fixed (True) or free to fit (False)

    Notes
    -----
    Parameters passed as strings are automatically treated as fixed.
    Example: R("0.86") creates a fixed resistor, R(0.86) is free to fit.
    """

    def __init__(self, *params):
        # Convert string parameters to float and mark as fixed
        self.params = []
        self.fixed_params = []
        for param in params:
            if isinstance(param, str):
                self.params.append(float(param))
                self.fixed_params.append(True)
            else:
                self.params.append(float(param))
                self.fixed_params.append(False)

    @abstractmethod
    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        """
        Calculate impedance Z(ω) with given parameters.

        Parameters
        ----------
        freq : ndarray of float
            Frequencies [Hz]
        params : list of float
            Element parameters

        Returns
        -------
        Z : ndarray of complex
            Impedance [Ω]
        """
        pass

    @abstractmethod
    def _scale(self, scalar: float) -> 'CircuitElement':
        """Return scaled version of element (for * operator)"""
        pass

    @abstractmethod
    def __repr__(self) -> str:
        """String representation"""
        pass

    # Operator overloading for circuit building
    def __sub__(self, other: Union['CircuitElement', 'Circuit']) -> 'Circuit':
        """Series connection: Z_total = Z1 + Z2"""
        from .circuit_builder import Series
        return Series([self, other])

    def __or__(self, other: Union['CircuitElement', 'Circuit']) -> 'Circuit':
        """Parallel connection: 1/Z_total = 1/Z1 + 1/Z2"""
        from .circuit_builder import Parallel
        return Parallel([self, other])

    def __mul__(self, scalar: float) -> 'CircuitElement':
        """Scale element parameter"""
        return self._scale(scalar)

    def __rmul__(self, scalar: float) -> 'CircuitElement':
        """Allow 2*R syntax"""
        return self.__mul__(scalar)

    def get_all_params(self) -> List[float]:
        """Get all parameters (for compatibility with Circuit class)"""
        return self.params

    def get_all_fixed_params(self) -> List[bool]:
        """Get fixed parameter flags (for compatibility with Circuit class)"""
        return self.fixed_params

    def update_params(self, params: List[float]) -> int:
        """
        Update element parameters with fitted values.

        Parameters
        ----------
        params : list of float
            New parameter values (must match number of element parameters)

        Returns
        -------
        n_consumed : int
            Number of parameters consumed from the list
        """
        n_params = len(self.params)
        self.params = list(params[:n_params])
        return n_params

    @abstractmethod
    def get_param_labels(self) -> List[str]:
        """
        Get parameter labels for display.

        Returns
        -------
        labels : list of str
            Parameter labels (e.g., ['R'], ['Q', 'n'], etc.)
        """
        pass


# ============================================================================
# Resistor
# ============================================================================

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


# ============================================================================
# Capacitor
# ============================================================================

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


# ============================================================================
# Constant Phase Element (Q)
# ============================================================================

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


# ============================================================================
# Inductor
# ============================================================================

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


# ============================================================================
# Warburg Element
# ============================================================================

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


# ============================================================================
# Warburg Open (finite-length bounded diffusion)
# ============================================================================

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


# ============================================================================
# K Element (Voigt with tau parametrization)
# ============================================================================

class K(CircuitElement):
    """
    Voigt element parametrized by resistance R and time constant τ.

    This is an alternative parametrization of the parallel R-C circuit (Voigt element)
    that uses the time constant τ = R×C instead of capacitance C. This parametrization
    is particularly useful because:

    1. τ directly relates to the characteristic frequency: f = 1/(2πτ)
    2. R and τ are more independent (better numerical conditioning)
    3. Consistent with DRT (Distribution of Relaxation Times) notation
    4. Used in Lin-KK test (Schönleber et al. 2014)

    Z_K = R / (1 + jωτ)

    Equivalent to (R || C) where C = τ/R.

    Parameters
    ----------
    R : float or str, optional
        Resistance [Ω] (default: 1000.0)
        If passed as string, the parameter is fixed during fitting.
    tau : float or str, optional
        Time constant [s] (default: 1e-4)
        If passed as string, the parameter is fixed during fitting.

    Notes
    -----
    Characteristic frequency: f_c = 1 / (2π τ)
    Capacitance: C = τ / R

    At ω = 1/τ (characteristic frequency):
        Z = R/2 * (1 - j)  (half-power point, phase = -45°)

    For ω << 1/τ: Z → R (capacitor blocks, full resistance)
    For ω >> 1/τ: Z → 0 (capacitor shorts)

    Examples
    --------
    >>> k = K(1000, 1e-4)      # R=1kΩ, τ=100μs → f=1.59 kHz, C=100 nF
    >>> k = K()                # default R=1kΩ, τ=100μs
    >>> k = K("1000", 1e-4)    # R fixed, τ free
    >>> k = K("1000", "1e-4")  # Both parameters fixed

    Convert from (R || C) to K:
    >>> # (R(1000) | C(1e-7)) is equivalent to K(1000, 1e-4)
    >>> # because τ = R×C = 1000 × 1e-7 = 1e-4 s

    References
    ----------
    Schönleber, M. et al. "A Method for Improving the Robustness of
    linear Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20–27 (2014)
    """

    def __init__(self, R: Union[float, str] = 1000.0, tau: Union[float, str] = 1e-4):
        super().__init__(R, tau)
        self.R = self.params[0]
        self.tau = self.params[1]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_val, tau_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return R_val / (1 + 1j * omega * tau_val)

    def _scale(self, scalar: float) -> 'K':
        """Scale R while preserving τ (changes C = τ/R)"""
        return K(scalar * self.R, self.tau)

    def get_param_labels(self) -> List[str]:
        return ['R', 'τ']

    def __repr__(self) -> str:
        R_str = f'"{self.R:.4g}"' if self.fixed_params[0] else f"{self.R:.4g}"
        tau_str = f'"{self.tau:.4g}"' if self.fixed_params[1] else f"{self.tau:.4g}"
        return f"K(R={R_str}, τ={tau_str})"

    def to_RC(self) -> 'Circuit':
        """
        Convert K element to equivalent (R || C) circuit.

        Returns
        -------
        circuit : Circuit
            Parallel R-C circuit equivalent to this K element

        Examples
        --------
        >>> k = K(1000, 1e-4)
        >>> rc = k.to_RC()  # Returns (R(1000) | C(1e-7))
        """
        from .circuit_builder import Parallel
        C_val = self.tau / self.R  # C = τ/R
        return Parallel([R(self.R), C(C_val)])

    @property
    def capacitance(self) -> float:
        """Get equivalent capacitance C = τ/R"""
        return self.tau / self.R

    @property
    def characteristic_freq(self) -> float:
        """Get characteristic frequency f = 1/(2πτ) [Hz]"""
        return 1.0 / (2 * np.pi * self.tau)


# ============================================================================
# Gerischer Element (coupled reaction-diffusion)
# ============================================================================

class G(CircuitElement):
    """
    Gerischer element for coupled reaction-diffusion processes.

    Z_G = sigma / sqrt(1 + j*omega*tau)

    Models systems where diffusion is coupled with a first-order chemical
    reaction, such as:
    - SOFC cathodes (oxygen reduction)
    - Porous electrodes with surface reactions
    - Mixed ionic-electronic conductors (MIECs)

    Parameters
    ----------
    sigma : float or str, optional
        Pre-factor [Ohm*s^(1/2)] (default: 100.0)
        Related to the diffusion-reaction resistance.
        If passed as string, the parameter is fixed during fitting.
    tau : float or str, optional
        Reaction time constant [s] (default: 1e-3)
        If passed as string, the parameter is fixed during fitting.

    Notes
    -----
    Limiting behavior:
    - omega -> 0: Z -> sigma (real resistance)
    - omega -> inf: Z -> 0 (Warburg-like decay)

    Nyquist plot shows an asymmetric arc with a characteristic
    "tail" toward low frequencies, distinct from a symmetric
    RC semicircle.

    The Gerischer element differs from Warburg:
    - Warburg: pure diffusion (Z ~ 1/sqrt(omega))
    - Gerischer: diffusion + reaction (Z ~ 1/sqrt(1 + j*omega*tau))

    Examples
    --------
    >>> g = G(100, 1e-3)       # sigma=100, tau=1ms, both free
    >>> g = G()                # default values (both free)
    >>> g = G("100", 1e-3)     # sigma fixed, tau free
    >>> g = G("100", "1e-3")   # Both parameters fixed

    Typical circuit: R_s - G (series resistance + Gerischer)
    >>> circuit = R(10) - G(100, 1e-3)

    References
    ----------
    Gerischer, H. "Wechselstrompolarisation von Elektroden mit einem
    potentialbestimmenden Schritt beim Gleichgewichtspotential I"
    Zeitschrift fur Physikalische Chemie, 198, 286-313 (1951)
    """

    def __init__(self, sigma: Union[float, str] = 100.0,
                 tau: Union[float, str] = 1e-3):
        super().__init__(sigma, tau)
        self.sigma = self.params[0]
        self.tau = self.params[1]

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        sigma_val, tau_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return sigma_val / np.sqrt(1 + 1j * omega * tau_val)

    def _scale(self, scalar: float) -> 'G':
        """Scale sigma while preserving tau"""
        new_sigma = scalar * self.sigma
        # Preserve fixed status
        sigma_arg = f'"{new_sigma}"' if self.fixed_params[0] else new_sigma
        tau_arg = f'"{self.tau}"' if self.fixed_params[1] else self.tau
        return G(sigma_arg, tau_arg)

    def get_param_labels(self) -> List[str]:
        return ['σ_G', 'τ_G']

    def __repr__(self) -> str:
        sigma_str = f'"{self.sigma:.4g}"' if self.fixed_params[0] else f"{self.sigma:.4g}"
        tau_str = f'"{self.tau:.4g}"' if self.fixed_params[1] else f"{self.tau:.4g}"
        return f"G(σ={sigma_str}, τ={tau_str})"

    @property
    def characteristic_freq(self) -> float:
        """Get characteristic frequency f = 1/(2*pi*tau) [Hz]"""
        return 1.0 / (2 * np.pi * self.tau)


__all__ = ['R', 'C', 'Q', 'L', 'W', 'Wo', 'K', 'G', 'CircuitElement']
