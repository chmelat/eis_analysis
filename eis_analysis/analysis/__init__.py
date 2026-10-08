"""
Analysis module for EIS analysis.
"""

from .oxide import analyze_oxide_layer, estimate_permittivity, OxideAnalysisResult
from .local_exponent import local_exponent, LocalExponentResult

__all__ = [
    'analyze_oxide_layer',
    'estimate_permittivity',
    'OxideAnalysisResult',
    'local_exponent',
    'LocalExponentResult',
]
