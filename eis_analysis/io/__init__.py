"""
I/O module for loading and generating EIS data.
"""

from .data_loading import (
    load_data,
    load_csv_data,
    read_gamry_native,
    parse_dta_metadata,
    parse_ocv_curve,
    expected_points,
    LoadResult,
)
from .synthetic import generate_synthetic_data

__all__ = [
    'load_data',
    'load_csv_data',
    'read_gamry_native',
    'parse_dta_metadata',
    'parse_ocv_curve',
    'expected_points',
    'LoadResult',
    'generate_synthetic_data',
]
