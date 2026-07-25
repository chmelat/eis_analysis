"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.19.0'
__version_info__ = (0, 19, 0)
__release_date__ = '2026-07-25'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "plot_rl_fit_diagnostics() removed from eis_analysis and "
    "eis_analysis.visualization; the R-L diagnostic figure is returned by "
    "estimate_rinf_with_inductance(..., plot=True)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
