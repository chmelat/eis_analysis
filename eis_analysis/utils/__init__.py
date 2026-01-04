"""
Utility functions for EIS analysis.

This package contains shared utility functions used across multiple modules,
implementing the DRY (Don't Repeat Yourself) principle.
"""

from .compat import np_trapz
from .impedance import calculate_rpol, sort_by_frequency

__all__ = [
    'np_trapz',
    'calculate_rpol',
    'sort_by_frequency',
]
