"""
Utility functions for EIS analysis.

This package contains shared utility functions used across multiple modules,
implementing the DRY (Don't Repeat Yourself) principle.
"""

from .impedance import calculate_rpol, sort_by_frequency
from .validation import validate_eis_data

__all__ = [
    'calculate_rpol',
    'sort_by_frequency',
    'validate_eis_data',
]
