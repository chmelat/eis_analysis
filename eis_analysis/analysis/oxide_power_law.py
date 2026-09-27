"""
Power-law model of a CPE film (Hirschorn, Orazem et al.), for oxide analysis.

A resistivity falling as a power law, rho ~ (x/d)^(-1/(1-n)), from rho_0 at
the metal to rho_d at the electrolyte, with a uniform permittivity, makes the
film a CPE between f_0 = 1/(2*pi*rho_0*eps*eps0) and f_d = 1/(2*pi*rho_d*eps*eps0),
with

    Q = (eps*eps0)^n / (g * d * rho_d^(1-n)),   g = 1 + 2.88 * (1-n)^2.375

Above f_d every layer is capacitive and the film is an ideal capacitor.
rho_d is not determined by the spectrum, so it is an input; it enters only
as rho_d^(1-n), weakly for n near 1. Unlike Hsu-Mansfeld and Brug, the model needs no resistance from the
circuit, and its equivalent capacitance eps*eps0/d depends on eps.

References: Hirschorn et al., "Determination of effective capacitance and
film thickness from constant-phase-element parameters", Electrochim. Acta 55,
6218 (2010); "Constant-Phase-Element Behavior Caused by Resistivity
Distributions in Films" I/II, J. Electrochem. Soc. (2010),
doi:10.1149/1.3499565 (II).
"""

from typing import Any, Dict, List, Optional

import numpy as np

from .config import EPSILON_0


def _power_law_g(n: float) -> float:
    """Hirschorn's interpolation for the numerical factor g(n); g(1) = 1."""
    return 1.0 + 2.88 * (1.0 - n) ** 2.375


def _power_law_thickness_cm(Q_s: float, n: float, eps_r: float, rho_d: float) -> float:
    """d = (eps*eps0)^n / (g * Q_s * rho_d^(1-n)); Q_s per cm², rho_d in Ohm cm."""
    return (eps_r * EPSILON_0) ** n / (_power_law_g(n) * Q_s * rho_d ** (1.0 - n))


def _power_law_permittivity(Q_s: float, n: float, d_cm: float, rho_d: float) -> float:
    """eps = (d * g * Q_s * rho_d^(1-n))^(1/n) / eps0, the inverse of the above."""
    return (d_cm * _power_law_g(n) * Q_s * rho_d ** (1.0 - n)) ** (1.0 / n) / EPSILON_0


def _power_law(
    element_params: Dict[str, Any],
    area_cm2: float,
    rho_d: Optional[float],
    f_max: float,
    warnings: List[str],
    *,
    eps_r: Optional[float] = None,
    d_cm: Optional[float] = None
) -> Optional[float]:
    """
    Power-law thickness [cm] (given eps_r) or permittivity (given d_cm) of the
    dominant element, or None when rho_d was not given or the element is not
    a Q. Caveats are appended to `warnings`, among them a sweep reaching
    above f_d, where the film is no longer a CPE (the lower corner f_0 needs
    rho_0, which is not an input, so it is not checked).
    """
    if rho_d is None:
        return None
    if element_params.get('type') != 'Q':
        # The high-frequency estimate has no element, hence no 'type'
        which = ("the dominant element is not one" if 'type' in element_params
                 else "the capacitance is a high-frequency estimate, not a fitted Q")
        warnings.append(f"Power-law model applies to a CPE (Q) only - {which}, "
                        "so no power-law value is reported.")
        return None
    Q_s, n = element_params['Q'] / area_cm2, element_params['n']
    if eps_r is not None:
        result = _power_law_thickness_cm(Q_s, n, eps_r, rho_d)
    else:
        assert d_cm is not None
        result = eps_r = _power_law_permittivity(Q_s, n, d_cm, rho_d)
    f_d = 1.0 / (2.0 * np.pi * rho_d * eps_r * EPSILON_0)
    if f_max > f_d:
        warnings.append(
            f"Power-law model: f_delta = 1/(2*pi*rho_delta*eps*eps0) = {f_d:.3g} Hz "
            f"lies below the highest measured frequency {f_max:.3g} Hz. Above "
            "f_delta the film is an ideal capacitor, not a CPE, so a Q fitted "
            "across it is not the power-law Q - check rho_delta.")
    return result
