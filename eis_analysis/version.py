"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.63.0'
__version_info__ = (0, 63, 0)
__release_date__ = '2026-10-25'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "Python 3.13 or newer is required (was 3.11)",
    "DRT R_inf comes from the top-decade fit where it determines R_s clearly "
    "better: R_inf and R_pol change there",
    "find_optimal_lambda_gcv raises ValueError when no lambda gives a finite score",
    "validate_eis_data raises ValueError for Z = 0 (|Z| below ~1.5e-154 Ohm)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"
