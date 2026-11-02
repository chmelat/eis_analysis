"""
Parameter range of the R_inf window fits, relative to the window's |Z|.

Shared by the R-L-(R|Q) fit (estimate.py) and the R - (R|L) - (R|Q) check
(parallel_rl.py), so that both fits run under the same bounds and their BIC
comparison measures the model, not the bounds.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting.bounds import PARAMETER_BOUNDS

# Range of the window fit, relative to the window's impedance: every term from
# min|Z|/RINF_BOUND_RANGE to max|Z|*RINF_BOUND_RANGE over the window. Above
# 1e6 times the largest |Z| an arc resistance is an open arc, and Q the same
# for a short, to better than 1e-6 - under any instrument's resolution (~1e-4
# of |Z|) - so these are the upper bounds. R and L are bounded below by their
# physical limit 0 instead: a relative floor cut off what the data still
# determine (noise-free 0.05 Ohm R_s in front of a GOhm film, |Z| 3e3..2.5e5
# Ohm: an L floor at 1e-6 |Z| biased R_s to 0.064 Ohm, and R_inf fell back
# to 243 Ohm). The lower end of the range only places the starts, so that
# none is 0. Absolute PARAMETER_BOUNDS made R_inf
# depend on the units: the open arc's R_k ran into R <= 1e10 Ohm (or Q into
# 0.1, L into 1e-4 H), and Z -> 1000 Z shifted R_inf on 51 % of the stress
# test's random spectra, by over 10 % on 7 % of them.
RINF_BOUND_RANGE = 1e6


def window_range(f: NDArray, Z: NDArray) -> Dict[str, Tuple[float, float]]:
    """Range per parameter type, each term's impedance within
    min|Z|/RANGE..max|Z|*RANGE.

    R: the value itself. L: its reactance omega*L across the window. Q: its
    impedance 1/(Q omega^n) across the window and the whole n range. All scale
    with Z, so R_inf does not depend on the units.
    """
    Z_abs = np.abs(Z)
    lo, hi = float(Z_abs.min()) / RINF_BOUND_RANGE, float(Z_abs.max()) * RINF_BOUND_RANGE
    w_lo, w_hi = 2 * np.pi * float(f.min()), 2 * np.pi * float(f.max())
    n_lo, n_hi = PARAMETER_BOUNDS['n']
    w_n = [w**n for w in (w_lo, w_hi) for n in (n_lo, n_hi)]
    return {
        'R': (lo, hi),
        'L': (lo / w_hi, hi / w_lo),
        'Q': (1 / (hi * max(w_n)), 1 / (lo * min(w_n))),
        'n': (n_lo, n_hi),
    }


def window_bounds(labels: Sequence[str], ranges: Dict[str, Tuple[float, float]]
                  ) -> Tuple[List[float], List[float]]:
    """Fit bounds for parameters `labels`: the relative range, with R and L
    bounded below by their physical limit 0 (see RINF_BOUND_RANGE)."""
    lower = [0.0 if label in ('R', 'L') else ranges[label][0] for label in labels]
    upper = [ranges[label][1] for label in labels]
    return lower, upper


def clip_to(value: float, label: str, ranges: Dict[str, Tuple[float, float]]) -> float:
    """A start value clipped into the range of its parameter type."""
    return float(np.clip(value, *ranges[label]))
