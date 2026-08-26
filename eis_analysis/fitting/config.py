"""
Configuration constants for circuit fitting.

All values are physically justified or empirically determined based on
typical EIS experiments.

References
----------
.. [1] B. Boukamp, Solid State Ionics 20 (1986) 31-44
       "A package for impedance/admittance data analysis"
.. [2] T. Reshetenko et al., J. Power Sources 269 (2014) 344-362
       "Determination of polarization resistance from DRT"
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

DRT_MIN_EFFECTIVE_BINS = 7.0
"""
Minimum effective number of gamma bins for meaningful peak-shape analysis.

Measured as the participation ratio N_eff = (sum gamma)^2 / sum(gamma^2),
which is ~1 for a single-bin spike and grows to tens for a smooth
distribution. Below this threshold the DRT is too sparse/spiky to analyze
peak shape reliably (typically auto-lambda collapsing toward 0 on
low-noise data; see audit finding F3).

Calibration
-----------
Healthy DRT: N_eff ~ 9-20. Degenerate (auto-lambda -> 0): N_eff ~ 4-5.5.
A threshold of 7 cleanly separates the two. Advisory only (emits a
warning; does not alter gamma or detected peaks).
"""

DRT_PEAK_PROMINENCE_THRESHOLD = 0.015
"""
Minimum peak prominence as fraction of maximum (1.5%).

Prominence measures how much a peak stands out from its surroundings.
Prevents detection of small bumps as separate peaks.
Value 1.5% is a very sensitive setting for very subtle peaks.
WARNING: Small fluctuations in data may be detected as peaks.
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
# Differential Evolution Diagnostics
# =============================================================================

DE_STALLED_ERROR_PCT = 50.0
"""
Fit error [%] above which the DE stage is considered to have found nothing.

A weighted mean relative error of 50% means the model reproduces neither the
magnitude nor the shape of the spectrum: the population never left the region
where the prediction is dominated by a single element. The reported fit then
rests entirely on the local refinement, i.e. on a single starting point, which
is exactly what the global optimizer was supposed to avoid. Paired with a
refinement at least DE_STALLED_IMPROVEMENT_FACTOR times better, this is
reported as a warning rather than passing silently.
"""

DE_STALLED_IMPROVEMENT_FACTOR = 10.0
"""
How much better the refinement must be to call the DE stage stalled.

An order of magnitude separates "DE landed near the basin and least_squares
polished it" from "DE contributed nothing and least_squares did the fitting".
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
    'DRT_MIN_EFFECTIVE_BINS',
    'DRT_PEAK_PROMINENCE_THRESHOLD',
    'GMM_PEAK_HEIGHT_FACTOR',

    # Fit Quality Assessment
    'FIT_QUALITY_EXCELLENT_ERROR',
    'FIT_QUALITY_GOOD_ERROR',
    'DE_STALLED_ERROR_PCT',
    'DE_STALLED_IMPROVEMENT_FACTOR',

    # Automatic Circuit Suggestion
    'MAX_VOIGT_ELEMENTS',
    'PEAK_INTEGRATION_TOLERANCE',
    'RPOL_RATIO_WARNING_THRESHOLD_LOW',
    'RPOL_RATIO_WARNING_THRESHOLD_HIGH',

    # Plotting
    'PLOT_GRID_ALPHA',
]
