"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.20.0'
__version_info__ = (0, 20, 0)
__release_date__ = '2026-07-25'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "estimate_rinf_with_inductance() returns (RLKFitResult, fig) instead of "
    "the 5-tuple (R_inf, L, circuit, diagnostics, fig); read the fields off "
    "the dataclass",
    "use_voigt_fit parameter removed from calculate_drt() and "
    "estimate_rinf_with_inductance(); use use_rl_fit",
    "RinfEstimate.R_ct / .C_nF / .f_characteristic removed (always None)",
    "Ignored verbose parameter removed from fit_equivalent_circuit(), "
    "fit_circuit_diffevo(), fit_circuit_multistart() and "
    "estimate_rinf_with_inductance(); the -v/--verbose CLI flag is unaffected",
    "eis_analysis.utils.compat module and the eis_analysis.utils.np_trapz "
    "re-export removed",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"

# For compatibility
VERSION = __version__
