#!/usr/bin/env python3
"""Vypocet relativni permitivity oxidicke vrstvy z efektivni kapacity.

Model deskoveho kondenzatoru: eps_r = (C_eff / A) * d / eps_0

Pouziti:
    python3 oxide_permittivity.py <C_eff [F]> <plocha [cm2]> <tloustka [nm]>

Priklad:
    python3 oxide_permittivity.py 9.74e-7 1.0 20
"""

import sys

from eis_analysis.analysis.config import EPSILON_0  # permitivita vakua [F/cm]


def relative_permittivity(C_eff: float, area_cm2: float, thickness_nm: float) -> float:
    """Relativni permitivita z efektivni kapacity, plochy a tloustky vrstvy."""
    C_specific = C_eff / area_cm2  # [F/cm2]
    d_cm = thickness_nm * 1e-7     # nm -> cm
    return d_cm * C_specific / EPSILON_0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    C_eff, area_cm2, thickness_nm = map(float, sys.argv[1:])
    eps_r = relative_permittivity(C_eff, area_cm2, thickness_nm)
    print(f"eps_r = {eps_r:.2f}")
