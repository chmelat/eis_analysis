"""
Data validation module for EIS analysis.
"""

from .kramers_kronig import (
    kramers_kronig_validation,
    lin_kk_native,
    low_frequency_slope,
)
from .kk_result import KKResult, LinKKResult
from .kk_grid import find_optimal_extend_decades, reconstruct_impedance
from .noise import compute_pseudo_chisqr, estimate_noise_percent, estimate_noise_percent_median
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
    'estimate_noise_percent_median',
    'low_frequency_slope',
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
