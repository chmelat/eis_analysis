"""
Configuration constants for circuit fitting.

All values are physically justified or empirically determined based on
typical EIS experiments. See AUDIT_REPORT.md for detailed rationale.

References
----------
.. [1] B. Boukamp, Solid State Ionics 20 (1986) 31-44
       "A package for impedance/admittance data analysis"
.. [2] T. Reshetenko et al., J. Power Sources 269 (2014) 344-362
       "Determination of polarization resistance from DRT"
.. [3] M. Orazem, B. Tribollet, "Electrochemical Impedance Spectroscopy" (2008)
       Wiley, ISBN: 978-0-470-04140-6
"""

# =============================================================================
# DRT Peak Detection
# =============================================================================

DRT_PEAK_HEIGHT_THRESHOLD = 0.03
"""
Minimum relative peak height for detection as a separate process.

A peak must be at least 3% of maximum DRT height to be considered
a separate relaxation process. Lower values capture weak processes
but increase risk of false detections (noise, artifacts).

Reference
---------
Typical threshold 5-15% for process separation [2].
Value 3% is a very sensitive setting for detecting very weak processes.
WARNING: May capture noise in lower quality data.
"""

DRT_PEAK_PROMINENCE_THRESHOLD = 0.015
"""
Minimum peak prominence as fraction of maximum (1.5%).

Prominence measures how much a peak stands out from its surroundings.
Prevents detection of small bumps as separate peaks.
Value 1.5% is a very sensitive setting for very subtle peaks.
WARNING: Small fluctuations in data may be detected as peaks.
"""

DRT_PEAK_MIN_SPACING_DECADES = 0.5
"""
Minimum distance between peaks in log(tau) [decades].

Two peaks closer than 0.5 decades (factor 3.16x) are difficult to separate
experimentally and likely represent a single process.
"""

GMM_PEAK_HEIGHT_FACTOR = 0.05
"""
Minimum GMM peak height as fraction of maximum (5%).

GMM (Gaussian Mixture Model) peak detection is more sensitive than scipy.
Value 5% is a very sensitive setting for very weak relaxation processes.
WARNING: High sensitivity - may detect noise as separate peaks
in data with lower SNR (signal-to-noise ratio).
"""

# =============================================================================
# Fit Quality Assessment
# =============================================================================

FIT_QUALITY_EXCELLENT_ERROR = 1.0
"""
Threshold for excellent fit [%].

Relative error <1% indicates excellent model-data agreement.
"""

FIT_QUALITY_GOOD_ERROR = 10.0
"""
Threshold for good fit [%].

Relative error 1-10% is typical for good fits in real systems.
"""

# =============================================================================
# Automatic Circuit Suggestion
# =============================================================================

MAX_VOIGT_ELEMENTS = 4
"""
Maximum number of parallel RC (Voigt) elements in auto-suggested circuit.

More than 4 elements often leads to:
- Overfitting (too many parameters)
- Loss of physical meaning
- Unstable fit
- Parameter correlation

Recommendations per [1]:
- 1-2 elements: simple system (bulk + interface)
- 3-4 elements: complex system (bulk + 2-3 interfaces)
- >4 elements: likely overfit, consider DRT analysis
"""

PEAK_INTEGRATION_TOLERANCE = 0.1
"""
Tolerance for peak integration in DRT (+/-10%).

When computing R_i from peak integral, include the region
where gamma(tau) > peak_height * 0.1.
"""

RPOL_RATIO_WARNING_THRESHOLD_LOW = 0.5
"""
Lower threshold for R_pol ratio warning (50%).

If sum of R_i from peaks is <50% of R_pol from data, some processes
are likely missing (low frequencies not measured).
"""

RPOL_RATIO_WARNING_THRESHOLD_HIGH = 2.0
"""
Upper threshold for R_pol ratio warning (200%).

If sum of R_i from peaks is >200% of R_pol from data:
- Peaks are poorly integrated
- DRT normalization is incorrect
- Background noise
"""

# =============================================================================
# Initial Guess Estimation
# =============================================================================

DEFAULT_R0_GUESS = 100
"""
Default estimate for series resistance R0 [Ohm].

Typical value for aqueous electrolytes at cm scale:
- Very good electrolyte: ~10 Ohm
- Average electrolyte: ~100 Ohm
- Poor electrolyte: ~1000 Ohm
"""

DEFAULT_Q_N_GUESS = 0.7
"""
Default estimate for Q (CPE) exponent n (dimensionless).

Typical values in electrochemical systems [3]:
- Double layer capacitance: n ~ 0.9-1.0 (nearly ideal)
- Porous electrodes: n ~ 0.7-0.9
- Pore diffusion: n ~ 0.5-0.7
- Warburg diffusion: n = 0.5

n=0.7 is a reasonable compromise for initial guess.
"""

# =============================================================================
# Grid and Plotting
# =============================================================================

PLOT_GRID_ALPHA = 0.3
"""
Grid transparency in plots (30%).

Grid should be visible but unobtrusive.
"""

# =============================================================================
# Export all constants
# =============================================================================

__all__ = [
    # DRT Peak Detection
    'DRT_PEAK_HEIGHT_THRESHOLD',
    'DRT_PEAK_PROMINENCE_THRESHOLD',
    'DRT_PEAK_MIN_SPACING_DECADES',
    'GMM_PEAK_HEIGHT_FACTOR',

    # Fit Quality Assessment
    'FIT_QUALITY_EXCELLENT_ERROR',
    'FIT_QUALITY_GOOD_ERROR',

    # Automatic Circuit Suggestion
    'MAX_VOIGT_ELEMENTS',
    'PEAK_INTEGRATION_TOLERANCE',
    'RPOL_RATIO_WARNING_THRESHOLD_LOW',
    'RPOL_RATIO_WARNING_THRESHOLD_HIGH',

    # Initial Guess Estimation
    'DEFAULT_R0_GUESS',
    'DEFAULT_Q_N_GUESS',

    # Plotting
    'PLOT_GRID_ALPHA',
]
