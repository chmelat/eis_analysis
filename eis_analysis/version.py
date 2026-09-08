"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.32.0'
__version_info__ = (0, 32, 0)
__release_date__ = '2026-09-08'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "load_data(), load_csv_data() and read_gamry_native() return a LoadResult, "
    "not a (frequencies, Z) tuple",
    "find_optimal_M_mu() returns a MuOptimization, not a 6-tuple",
    "analyze_voigt_elements() returns a VoigtSuggestion, not a dict",
    "fit_voigt_chain_linear() returns a VoigtChainFit, not a (circuit, params) tuple",
    "log_metadata() and format_voigt_report() are gone - both were console "
    "output, and now live in the CLI layer",
    "Library modules no longer print: analyze_oxide_layer(), load_data() and "
    "the rest report through their result instead of the log",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
