"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.49.0'
__version_info__ = (0, 49, 0)
__release_date__ = '2026-09-30'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "fit_equivalent_circuit, fit_circuit_multistart and fit_circuit_diffevo "
    "return (result, Z_fit) without a figure; plot= removed "
    "(use visualization.plot_circuit_fit(frequencies, Z, result))",
    "plot_circuit_fit(frequencies, Z, result, title=None) takes a FitResult "
    "(Z_fit, circuit, figsize, Z_fit_at_data removed)",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
