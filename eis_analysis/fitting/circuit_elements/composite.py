"""
Composite circuit elements: Voigt element in R-τ parametrization (K) and
the Gerischer element for coupled reaction-diffusion (G).
"""
from __future__ import annotations

import numpy as np
from typing import List, Union, TYPE_CHECKING
from numpy.typing import NDArray

from .base import CircuitElement
from .basic import R, C

if TYPE_CHECKING:
    from ..circuit_builder import Series, Parallel
    Circuit = Union[Series, Parallel, 'CircuitElement']


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
        from ..circuit_builder import Parallel
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
