"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.62.0'
__version_info__ = (0, 62, 0)
__release_date__ = '2026-10-22'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "DRT R_inf defaults to the R-L-(R|Q) window fit (estimate_rinf) instead "
    "of the HF median: R_inf and R_pol change",
    "CLI --ri-fit removed (the fit is the default); --rinf-plot draws it",
    "RinfEstimate: method 'median' gone, n_points_used and warning removed, "
    "estimate added",
    "run_rinf_estimation removed; run_drt_analysis lost R_inf_computed",
    "analyze_voigt_elements takes the DRT's R_inf as a required argument",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"
