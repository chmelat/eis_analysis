"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.11.0'
__version_info__ = (0, 11, 0)
__release_date__ = '2026-01-08'

# Breaking changes in this version
__breaking_changes__ = [
    "Weighting names swapped: 'modulus' now means w=1/|Z|, 'proportional' now means w=1/|Z|^2",
    "Default weighting is now 'modulus' (w=1/|Z|)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
