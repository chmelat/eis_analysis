"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.40.0'
__version_info__ = (0, 40, 0)
__release_date__ = '2026-09-25'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    'rinf_estimation: estimate_rinf() -> RinfResult replaces estimate_rinf_with_inductance()',
    'calculate_drt: use_rl_fit removed, DRTResult.figure_rinf removed',
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
