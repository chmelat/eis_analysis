"""
I/O module for loading and generating EIS data.
"""

from .data_loading import (
    load_data,
    load_csv_data,
    read_gamry_native,
    parse_dta_metadata,
    parse_ocv_curve,
    log_metadata
)
from .synthetic import generate_synthetic_data

__all__ = [
    'load_data',
    'load_csv_data',
    'read_gamry_native',
    'parse_dta_metadata',
    'parse_ocv_curve',
    'log_metadata',
    'generate_synthetic_data',
]
