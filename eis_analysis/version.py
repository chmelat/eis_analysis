"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.14.0'
__version_info__ = (0, 14, 0)
__release_date__ = '2026-06-27'

# Breaking changes in this version
__breaking_changes__ = [
    "Removed DRT term-type classification (--classify-terms flag, "
    "term_classification module, classify_terms parameter). Peak width in "
    "DRT is set by regularization, not by element physics, so the "
    "Voigt/CPE classification was not a reliable measurement. Determine "
    "element type (C vs Q) from circuit fitting instead.",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
