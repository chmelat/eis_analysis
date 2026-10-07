"""
I/O module for loading and generating EIS data.
"""

from .gamry import (
    load_data,
    read_gamry_native,
    parse_dta_metadata,
    parse_ocv_curve,
    expected_points,
)
from .csv_loading import load_csv_data
from .spectrum import LoadResult
from .synthetic import generate_synthetic_data
from .export import fit_result_record, save_fit_result

__all__ = [
    'load_data',
    'load_csv_data',
    'read_gamry_native',
    'parse_dta_metadata',
    'parse_ocv_curve',
    'expected_points',
    'LoadResult',
    'generate_synthetic_data',
    'fit_result_record',
    'save_fit_result',
]
