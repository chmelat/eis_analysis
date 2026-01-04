"""
Data validation module for EIS analysis.
"""

from .kramers_kronig import kramers_kronig_validation, lin_kk_native

__all__ = [
    'kramers_kronig_validation',
    'lin_kk_native',
]
