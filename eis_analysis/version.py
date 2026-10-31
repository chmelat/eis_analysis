"""
Version information for EIS Analysis Toolkit.

This is the SINGLE SOURCE OF TRUTH for version information.
All other files should import from here.
"""

__version__ = '0.64.0'
__version_info__ = (0, 64, 0)
__release_date__ = '2026-10-31'

# Breaking changes in this version
__breaking_changes__: list[str] = [
    "Gamry files recorded with drift correction (ZrealDrCor/ZimagDrCor columns) "
    "give the drift-corrected impedance; --raw-impedance / drift_corrected=False "
    "reads the measured one",
]

# Human-readable version string
def get_version_string():
    """Return formatted version string."""
    return f"v{__version__} ({__release_date__})"
