"""
R_inf estimation module for EIS analysis.

Uses R-L-K() linear fit for all data types.
Model: Z(omega) = R_s + j*omega*L + R_k/(1+j*omega*tau)

Clean design: No logging in core functions, all diagnostics returned as data.
"""

from .rlk_fit import fit_rlk_model, estimate_rinf_with_inductance, RLKFitResult
from .data_selection import select_highest_decade, DataSelectionResult

__all__ = [
    # Main functions
    'estimate_rinf_with_inductance',
    'fit_rlk_model',
    'select_highest_decade',
    # Result dataclasses
    'RLKFitResult',
    'DataSelectionResult',
]
