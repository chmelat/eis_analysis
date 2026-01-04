"""
Analysis module for EIS analysis.
"""

from .oxide import analyze_oxide_layer, estimate_permittivity, OxideAnalysisResult

__all__ = [
    'analyze_oxide_layer',
    'estimate_permittivity',
    'OxideAnalysisResult',
]
