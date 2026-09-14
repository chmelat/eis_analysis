"""
Composite circuit elements: Voigt element in R-τ parametrization (K), the
Gerischer element for coupled reaction-diffusion (GE) and the Young-Göhr
passive-layer element (YG).
"""
from __future__ import annotations

import numpy as np
from typing import List, Tuple, Union, TYPE_CHECKING
from numpy.typing import NDArray

from .base import CircuitElement, param_property
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

    R = param_property(0)
    tau = param_property(1)

    def __init__(self, R: Union[float, str] = 1000.0, tau: Union[float, str] = 1e-4):
        super().__init__(R, tau)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        R_val, tau_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return R_val / (1 + 1j * omega * tau_val)

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


class GE(CircuitElement):
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
        Pre-factor [Ohm] (default: 100.0)
        The DC limit of the element: Z(omega -> 0) = sigma. (The sqrt is
        dimensionless in this tau parametrization, unlike the equivalent
        Z = Z_0 / sqrt(k + j*omega) form where Z_0 is in Ohm*s^(-1/2).)
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
    >>> g = GE(100, 1e-3)       # sigma=100, tau=1ms, both free
    >>> g = GE()                # default values (both free)
    >>> g = GE("100", 1e-3)     # sigma fixed, tau free
    >>> g = GE("100", "1e-3")   # Both parameters fixed

    Typical circuit: R_s - GE (series resistance + Gerischer)
    >>> circuit = R(10) - GE(100, 1e-3)

    References
    ----------
    Gerischer, H. "Wechselstrompolarisation von Elektroden mit einem
    potentialbestimmenden Schritt beim Gleichgewichtspotential I"
    Zeitschrift fur Physikalische Chemie, 198, 286-313 (1951)
    """

    sigma = param_property(0)
    tau = param_property(1)

    def __init__(self, sigma: Union[float, str] = 100.0,
                 tau: Union[float, str] = 1e-3):
        super().__init__(sigma, tau)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        sigma_val, tau_val = params[0], params[1]
        omega = 2 * np.pi * freq
        return sigma_val / np.sqrt(1 + 1j * omega * tau_val)

    def get_param_labels(self) -> List[str]:
        return ['σ_GE', 'τ_GE']

    def __repr__(self) -> str:
        sigma_str = f'"{self.sigma:.4g}"' if self.fixed_params[0] else f"{self.sigma:.4g}"
        tau_str = f'"{self.tau:.4g}"' if self.fixed_params[1] else f"{self.tau:.4g}"
        return f"GE(σ={sigma_str}, τ={tau_str})"

    @property
    def characteristic_freq(self) -> float:
        """Get characteristic frequency f = 1/(2*pi*tau) [Hz]"""
        return 1.0 / (2 * np.pi * self.tau)


# Above 1/p = 709 the factor exp(1/p) overflows float64 (exp(709) ~ 8e307,
# exp(710) = inf). `_yg_log_terms` never forms that factor, so the impedance
# itself stays finite below this p; the constant marks where YG's remaining
# closed forms (R_dc) and the p -> 0 degeneracy take over instead.
YG_P_MIN = 1.0 / 709.0


def _yg_log_terms(
    omega: NDArray[np.float64], p: float, tau: float
) -> Tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Return ln(1 + jωτ·e^(1/p)) and ln(1 + jωτ) without overflowing.

    The formula as printed overflows for p < 1/709 -- inside YG's own
    parameter bounds -- because it builds e^(1/p) as a number. Factoring that
    exponential out of the logarithm removes it:

        ln(1 + a·e^u) = ln(e^u·(e^-u + a)) = u + ln(e^-u + a)

    e^-u underflows to zero harmlessly, and a = jωτ is purely imaginary, so
    the sum adds into an empty real part and never cancels. The remaining
    u + ln(...) does cancel once e^-u dominates, costing 4e-11 relative at
    p = 0.5 and nothing measurable below p = 0.1 -- far under the 1e-6 the
    Jacobian is held to.

    Splitting ln(A/B) into ln(A) - ln(B) is safe here rather than a branch-cut
    hazard: with ω > 0 and τ > 0 both arguments lie in the first quadrant, so
    arg(A/B) stays in (-π/2, π/2) and no cut is crossed. The split is done for
    the overflow, not for the branch.

    `impedance` and the analytic Jacobian both call this, so the two can never
    evaluate different logarithms.
    """
    a = 1j * omega * tau
    u = 1.0 / p
    return u + np.log(np.exp(-u) + a), np.log1p(a)


class YG(CircuitElement):
    """
    Young-Göhr element: a passive layer whose conductivity decays exponentially.

    Z_YG = p/(jωC) · ln[ (1 + jωτ·e^(1/p)) / (1 + jωτ) ]

    Models a dielectric film into which conductivity penetrates from one side
    and falls off exponentially with depth - oxide layers on Fe, Al, Ti and Ta,
    and organic coatings under soaking. It is the physically consistent
    substitute for a CPE wherever the origin of the dispersion is known to be
    such a profile: instead of a bare exponent it carries the two quantities
    the film actually has, a capacitance and a penetration depth.

    Parametrised as (C, p, τ) rather than by the profile itself because all
    three are separately observable. C is the high-frequency plateau, p is a
    ratio the geometry fixes, and τ is read off the capacitive corner. The
    underlying (ε, ρ, d, A) are four coupled quantities the impedance cannot
    separate, which is Zahner's reason for normalising them away.

    Parameters
    ----------
    C : float or str, optional
        Total capacity of the layer neglecting its conductivity [F]
        (default: 1e-5). This is the high-frequency limit capacitance,
        C = ε₀·ε_r·A/d.
        If passed as string, the parameter is fixed during fitting.
    p : float or str, optional
        Relative penetration depth δ/d of the conductivity within the
        dielectric [dimensionless] (default: 0.05). p << 1 is a strong
        gradient in conductivity.
        If passed as string, the parameter is fixed during fitting.
    tau : float or str, optional
        Time constant of the virtual RC element built from the slice of
        dielectric at the site of highest conductivity [s] (default: 0.1),
        τ = ε₀·ε_r·ρ(x=0).
        If passed as string, the parameter is fixed during fitting.

    Notes
    -----
    Limiting behaviour, with E = e^(1/p):

    - ω >> 1/τ:     Z -> 1/(jωC), an ideal capacitor. The corner sits at
      `characteristic_freq` = 1/(2πτ).
    - ω << 1/(τE):  Z -> R_dc = p·τ·(E - 1)/C, a real resistance. That corner,
      `dc_corner_freq`, lies E times below the capacitive one, so for any
      p <~ 0.1 it is decades outside a real measurement window and R_dc is an
      extrapolation of the model rather than a fitted arc.
    - between the two the phase is nearly constant, which is the CPE-like
      band the element exists to explain.
    - p -> 0:       Z -> 1/(jωC). The model degenerates to a plain capacitor
      and p and τ stop being identifiable, which is why p sitting at its
      lower bound is worth reporting.

    Zahner gives a closed approximation for the phase in that middle band,
    useful as a sanity check rather than as a definition:

        φ ≈ -90°·(1 - q),    q = 1 / (ln(ωτ) + 1/p)

    Capacitance, not permittivity: the element carries no geometry. Convert
    with d = ε₀·ε_r·A/C in the analysis layer, where the electrode area is
    already known; the penetration depth is then δ = p·d.

    Examples
    --------
    >>> yg = YG(1e-5, 0.05, 0.1)    # Zahner's simulated example
    >>> yg = YG()                   # default values
    >>> yg = YG("1e-5", 0.05, 0.1)  # C fixed, p and tau free

    Typical oxide circuit: L - R_s - YG
    >>> circuit = L(1e-6) - R(20) - YG(1e-5, 0.05, 0.1)

    References
    ----------
    H. Göhr, "Impedance modelling of porous electrodes", Electrochemical
    Applications 1/97.
    Zahner Analysis user manual, 11/2023, section 2.3.9 "Young-Göhr impedance".
    """

    C = param_property(0)
    p = param_property(1)
    tau = param_property(2)

    def __init__(self, C: Union[float, str] = 1e-5,
                 p: Union[float, str] = 0.05,
                 tau: Union[float, str] = 0.1):
        super().__init__(C, p, tau)

    def impedance(self, freq: NDArray[np.float64],
                  params: List[float]) -> NDArray[np.complex128]:
        C_val, p_val, tau_val = params[0], params[1], params[2]
        omega = 2 * np.pi * freq
        if p_val <= YG_P_MIN:
            # The p -> 0 limit, an ideal capacitor. Reachable despite the
            # bounds: a string-fixed parameter, YG(1e-5, "5e-4", 0.1), skips
            # them entirely. Without this the element would return nan there
            # (0 · inf) instead of its own limit.
            return 1 / (1j * omega * C_val)
        t1, t2 = _yg_log_terms(omega, p_val, tau_val)
        return p_val / (1j * omega * C_val) * (t1 - t2)

    def get_param_labels(self) -> List[str]:
        return ['C_YG', 'p_YG', 'τ_YG']

    def __repr__(self) -> str:
        C_str = f'"{self.C:.4g}"' if self.fixed_params[0] else f"{self.C:.4g}"
        p_str = f'"{self.p:.4g}"' if self.fixed_params[1] else f"{self.p:.4g}"
        tau_str = f'"{self.tau:.4g}"' if self.fixed_params[2] else f"{self.tau:.4g}"
        return f"YG(C={C_str}, p={p_str}, τ={tau_str})"

    @property
    def R_dc(self) -> float:
        """DC resistance of the layer, R = p·τ·(e^(1/p) - 1)/C [Ω].

        The omega -> 0 limit of the impedance. Returns inf for p at or below
        YG_P_MIN, where e^(1/p) overflows float64 and the element is a pure
        capacitor with no DC path - the same value the limit approaches.

        Almost always an extrapolation: it is only measured below
        `dc_corner_freq`, which sits e^(1/p) times under the capacitive
        corner. Callers that rank elements by resistance must not treat it
        as a fitted one.
        """
        if self.p <= YG_P_MIN:
            return float('inf')
        return float(self.p * self.tau * np.expm1(1.0 / self.p) / self.C)

    @property
    def characteristic_freq(self) -> float:
        """Capacitive corner f = 1/(2πτ) [Hz]; above it the element is C."""
        return 1.0 / (2 * np.pi * self.tau)

    @property
    def dc_corner_freq(self) -> float:
        """Resistive corner f = e^(-1/p)/(2πτ) [Hz]; below it the element is R_dc.

        Underflows to 0.0 for very small p, which is the correct statement:
        no finite frequency reaches the DC plateau.
        """
        return float(self.characteristic_freq * np.exp(-1.0 / self.p))
