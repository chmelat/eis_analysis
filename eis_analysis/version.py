"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.50.0'
__version_info__ = (0, 50, 0)
__release_date__ = '2026-09-30'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "fit_equivalent_circuit, fit_circuit_multistart and fit_circuit_diffevo "
    "need a circuit with get_param_labels/get_all_fixed_params "
    "(no 1e-15..1e15 / all-free fallback; AttributeError otherwise)",
    "eis_analysis.version.VERSION removed (use __version__)",
    "load_csv_data raises ValueError on a header that names only some of "
    "frequency, Z_real and Z_imag, or one of them twice",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"
