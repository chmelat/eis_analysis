"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.23.0'
__version_info__ = (0, 23, 0)
__release_date__ = '2026-08-27'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "Removed the element scaling operator (2 * R(100)) and the CPE exponent "
    "operator (Q(1e-4) ** 0.9). Write the value directly: R(200), Q(1e-4, 0.9).",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
