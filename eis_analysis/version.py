"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.21.0'
__version_info__ = (0, 21, 0)
__release_date__ = '2026-07-29'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "estimate_permittivity() returns OxideAnalysisResult instead of float; "
    "use result.permittivity for the previous return value",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
