"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.10.2'
__version_info__ = (0, 10, 2)
__release_date__ = '2026-01-08'

# Breaking changes in this version
__breaking_changes__ = [
    "Weighting option 'square' renamed to 'modulus' (w=1/|Z|^2)",
    "Default weighting changed from 'sqrt' to 'proportional'",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
