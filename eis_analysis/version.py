"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.55.0'
__version_info__ = (0, 55, 0)
__release_date__ = '2026-10-06'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "Python 3.11 or newer is required (was 3.9)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"
