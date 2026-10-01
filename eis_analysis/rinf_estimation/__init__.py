"""
R_inf estimation module for EIS analysis.

R-L-(R|Q) fit over the top frequency decades, with Re(Z) at the highest
frequency with Im(Z) <= 0 (an upper bound) as the fallback when the fit does
not determine R_inf. The HF median
is the DRT's default R_inf.

Clean design: No logging in core functions, all diagnostics returned as data.
"""

from .estimate import estimate_rinf, hf_median, RinfResult

__all__ = [
    'estimate_rinf',
    'hf_median',
    'RinfResult',
]
