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
)
from .gcv import (
    compute_gcv_score,
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
    compute_lcurve_point,
    find_lcurve_corner,
)
from ..rinf_estimation import estimate_rinf_with_inductance
from .peaks import gmm_peak_detection, GMM_AVAILABLE

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
    # GCV functions
    'compute_gcv_score',
    'find_optimal_lambda_gcv',
    'find_optimal_lambda_hybrid',
    'compute_lcurve_point',
    'find_lcurve_corner',
    # R_inf estimation
    'estimate_rinf_with_inductance',
    # Peak detection
    'gmm_peak_detection',
    'GMM_AVAILABLE',
]
