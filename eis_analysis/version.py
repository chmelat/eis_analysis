"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.15.0'
__version_info__ = (0, 15, 0)
__release_date__ = '2026-06-29'

# Breaking changes in this version
__breaking_changes__ = [
    "GMM peak detection rewritten to use a weighted EM fitted directly to "
    "gamma(tau) instead of integer point replication (audit finding F1). "
    "Removed the sklearn dependency and the public GMM_AVAILABLE flag; GMM is "
    "now always available (pure numpy/scipy). gmm_peak_detection returns a "
    "WeightedGMMResult instead of a sklearn GaussianMixture.",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
