"""
R_inf estimation module for EIS analysis.

Uses R-L-K() linear fit for all data types.
Model: Z(omega) = R_s + j*omega*L + R_k/(1+j*omega*tau)

Modules:
- rlk_fit.py: Core R-L-K() linear fitting
- data_selection.py: High-frequency data selection

Public API
----------
- estimate_rinf_with_inductance: Robust R_inf estimation using R-L-K() fit
- fit_rlk_model: R-L-K() linear fitting (lower-level API)
- select_highest_decade: Select highest frequency decade
"""

from .rlk_fit import fit_rlk_model, estimate_rinf_with_inductance
from .data_selection import select_highest_decade

__all__ = [
    'estimate_rinf_with_inductance',
    'fit_rlk_model',
    'select_highest_decade',
]
