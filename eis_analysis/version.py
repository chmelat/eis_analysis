"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.25.0'
__version_info__ = (0, 25, 0)
__release_date__ = '2026-08-27'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "Element parameter attributes (R.R, K.tau, CC.alpha, ...) are now "
    "read-only views of .params. Assigning to one raises AttributeError "
    "instead of silently desynchronising it from the impedance.",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
