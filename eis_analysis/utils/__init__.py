"""
Utility functions for EIS analysis.

This package contains shared utility functions used across multiple modules,
implementing the DRY (Don't Repeat Yourself) principle.
"""

from .impedance import calculate_rpol, sort_by_frequency

__all__ = [
    'calculate_rpol',
    'sort_by_frequency',
]
