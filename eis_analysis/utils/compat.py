"""
NumPy compatibility utilities.

Handles differences between NumPy versions to ensure consistent behavior
across different installations.

NumPy Version Compatibility
----------------------------
- NumPy >= 2.0: Uses 'trapezoid' (new name, recommended)
- NumPy < 2.0: Uses 'trapz' (deprecated in 2.0 but still available)

This module provides a single import point that works with both versions,
following the Single Source of Truth principle from CLAUDE.md.

Examples
--------
>>> from eis_analysis.utils.compat import np_trapz
>>> import numpy as np
>>> x = np.linspace(0, 1, 100)
>>> y = x**2
>>> area = np_trapz(y, x)
>>> print(f"Area under curve: {area:.4f}")
Area under curve: 0.3333
"""

try:
    from numpy import trapezoid as np_trapz
except ImportError:
    # Fallback for NumPy < 2.0
    from numpy import trapz as np_trapz

__all__ = ['np_trapz']
