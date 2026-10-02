"""
Data validation module for EIS analysis.
"""

from .kramers_kronig import (
    kramers_kronig_validation,
    lin_kk_native,
    KKResult,
    LinKKResult,
    compute_pseudo_chisqr,
    estimate_noise_percent,
    find_optimal_extend_decades,
    reconstruct_impedance,
)
from .zhit import (
    zhit_validation,
    zhit_reconstruct_magnitude,
    ZHITResult,
)
from .thd import (
    thd_check,
    THDChannel,
    THDResult,
    THD_THRESHOLD,
)
from .outliers import (
    find_outliers,
    OutlierPoint,
    OutlierReport,
)

__all__ = [
    'kramers_kronig_validation',
    'lin_kk_native',
    'KKResult',
    'LinKKResult',
    'compute_pseudo_chisqr',
    'estimate_noise_percent',
    'find_optimal_extend_decades',
    'reconstruct_impedance',
    'zhit_validation',
    'zhit_reconstruct_magnitude',
    'ZHITResult',
    'thd_check',
    'THDChannel',
    'THDResult',
    'THD_THRESHOLD',
    'find_outliers',
    'OutlierPoint',
    'OutlierReport',
]
