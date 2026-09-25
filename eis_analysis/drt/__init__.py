"""
DRT (Distribution of Relaxation Times) analysis module.
"""

from .core import (
    calculate_drt,
    DRTResult,
    DRTDiagnostics,
    RinfEstimate,
    LambdaSelection,
    NNLSSolution,
    DRTMatrices,
    LambdaProbePoint,
    PeakStability,
    StabilityDiagnostics,
)
from .gcv import (
    compute_gcv_score,
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
    compute_lcurve_point,
    find_lcurve_corner,
)
from .peaks import gmm_peak_detection

__all__ = [
    # Main function
    'calculate_drt',
    # Result dataclasses
    'DRTResult',
    'DRTDiagnostics',
    'RinfEstimate',
    'LambdaSelection',
    'NNLSSolution',
    'DRTMatrices',
    'LambdaProbePoint',
    'PeakStability',
    'StabilityDiagnostics',
    # GCV functions
    'compute_gcv_score',
    'find_optimal_lambda_gcv',
    'find_optimal_lambda_hybrid',
    'compute_lcurve_point',
    'find_lcurve_corner',
    # Peak detection
    'gmm_peak_detection',
]
