"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.47.0'
__version_info__ = (0, 47, 0)
__release_date__ = '2026-09-29'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "DRT lambda is grid-independent: manual --lambda / lambda_reg values "
    "are on a new scale (typical 1e-9 to 1e-3)",
    "LinKKResult is an alias of KKResult (positional construction breaks); "
    "KK and Z-HIT validation return no figures (use plot_kk_validation / "
    "plot_zhit_validation)",
    "LambdaSelection.hybrid_stage is 'lcurve' or 'gcv' "
    "('lcurve_correction', 'geometric_mean' removed)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
