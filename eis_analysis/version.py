"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.17.1'
__version_info__ = (0, 17, 1)
__release_date__ = '2026-07-11'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "estimate_R_linear returns a 4-tuple (elements, residual, L_value, C_value)",
    "find_optimal_M_mu returns a 6-tuple (+ C_value)",
    "find_optimal_extend_decades returns a 6-tuple (+ C_value)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
