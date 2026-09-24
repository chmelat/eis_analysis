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

GMM_N_COMPONENTS_RANGE = (1, 6)
"""
Range of GMM component counts tried during BIC model selection.

Lower bound 1: a single relaxation process is a legitimate DRT.
Upper bound 6: an EIS spectrum spanning the usual 6-7 decades cannot resolve
more than roughly one process per decade, so beyond six components BIC is
fitting noise rather than physics. Widen it only when the selected count lands
on the upper bound (the CLI warns when it does).
"""

DRT_PEAK_EDGE_DECADES = 0.7
"""
Distance from the edge of the measured window below which a peak is flagged
as boundary-sensitive [decades of tau].

The measured window is tau = 1/(2*pi*f) for f_max and f_min. A peak near either
end is only half-supported by data: the kernel 1/(1+j*omega*tau) has no
measurement on one flank, the regularization is free to shape that flank, and
without a grid extension the basin integral that yields R_estimate is
truncated by the end of the array. A peak past the slow end (possible with
`tau_extend_decades`) has a negative distance and is flagged as well. The peak
position and area are then far less certain than the numbers alone suggest.

Value 0.7: a Gaussian-like DRT peak of typical width carries most of its area
within roughly +-0.7 decade of its maximum, so a peak closer than that to the
edge has a materially truncated basin. It is a conservative automation
heuristic motivated by sampling and localization limits (Macdonald 2000,
DOI 10.1088/0266-5611/16/5/324), not a physical constant - a peak flagged here
is a peak to interpret carefully, not a peak to discard.
"""

DRT_EDGE_BIN_RPOL_FRACTION = 0.05
"""
Share of R_pol heaped against an end of the tau grid above which pile-up is
reported.

Non-negative NNLS cannot represent response whose time constant lies outside
the tau grid (series inductance, an unresolved diffusion tail, a process slower
than the lowest measured frequency - or, with `tau_extend_decades`, slower
than the extended grid). It disposes of that response by heaping gamma up
against the first or last bin, which inflates the R_estimate of the peak that
absorbs it.

Measured as the mass of the falling run leaving the boundary, not of the
outermost bin: a relaxation inside the window makes gamma rise from the edge
towards its maximum, so a descending run at the edge means the maximum lies
outside it. That makes the number independent of `n_tau` - the lobe has a
width in decades, and a finer grid only spreads the same mass over more bins.
A per-bin measure does not: on a spectrum with a 50 s process measured to
10 mHz it read 0.068 at n_tau = 100 but 0.038 at n_tau = 400, so `-n 400`
silently lost a warning that `-n 100` emitted on identical data.

Value 0.05: a clean spectrum scores exactly 0, because gamma rises from both
edges towards its peaks, so the threshold only has to separate a small
boundary lobe from a real one. 5% of R_pol is small enough to catch pile-up
early and large enough not to fire on the shoulder of a legitimate peak that
happens to sit at the edge.
"""

DRT_EXTRAPOLATED_RPOL_FRACTION = 0.05
"""
Share of R_pol placed beyond the slow end of the measured window above which
the DRT warns that part of it is extrapolated.

Only possible with `tau_extend_decades > 0`. Mass past tau_max is constrained
by the high-frequency flank of its response alone, so its position and size
are extrapolations. Value 0.05, measured (sqrt weighting, 1 % proportional
noise): spectra whose peaks all sit inside the window leak at most 2.1 % of
R_pol past it at a 1-decade extension; a genuine peak 0.3 decade inside the
slow edge puts 7 % there, and warning about that is right.
"""

DRT_TAU_EXTEND_STEPS = (0.5, 1.0, 1.5)
"""
Extensions tried by `tau_extend_decades='auto'` [decades beyond tau_max].

'auto' extends only when the unextended DRT piles up at the slow end, and then
takes the smallest step that clears the pile-up; if none does, it keeps the
unextended grid. Picking by fit quality instead would always pick the largest
step: every extra decade adds freedom and lowers the residual even where it is
wrong (a blocking electrode). Capped at 1.5: in benchmarks an RC 0.5 decade
past the window closes at 1.0, one a full decade past closes at none of
these, and a wider grid only lets low-frequency noise leak further.
"""

DRT_LF_RC_RATIO_MIN = 0.2
"""
Low-frequency ratio r = (-dZ'/d ln omega) / (-Z'') below which
`tau_extend_decades='auto'` does not extend: the end is capacitive.

For a single RC, r = 2x/(1+x^2) with x = omega*tau. r is fitted over the four
lowest frequencies, so x is taken at their geometric centre; at 10 points per
decade that is 0.15 decade above f_min, and an RC 0.5 decade past the window
reads r = 0.44, one a full decade past r = 0.15. A series capacitor gives
r = 0 and a CPE r = n*cot(n*pi/2) (0.14 at n = 0.9). Extending the grid for
those turns the capacitor into a huge extrapolated peak and, in benchmarks,
pushed the real peak below the detection threshold. Value 0.2: below it the
process lies about a decade or more past the window, which no step of
DRT_TAU_EXTEND_STEPS closes anyway. It cannot separate a CPE with n < ~0.75
from an RC just past the window - both look alike at the lowest frequencies.
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

SIGNIFICANCE_NEGLIGIBLE = 0.01
"""
Below this significance a parameter may be dropped from the model.

From the Zahner Analysis manual (11/2023), section 2.2.2: "Significance values
much less than 0.01 usually indicate that the corresponding impedance element
may be omitted." The scale it sits on is fixed by the definition, not chosen:
for an element entering the impedance linearly the significance is bounded by
1 and is roughly the largest fraction of |Z| the parameter accounts for, so
0.01 means the parameter never moves the modulus by more than a percent
anywhere in the measured window. See compute_significance in diagnostics.py.
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
    'DRT_PEAK_EDGE_DECADES',
    'DRT_EDGE_BIN_RPOL_FRACTION',
    'DRT_EXTRAPOLATED_RPOL_FRACTION',
    'DRT_TAU_EXTEND_STEPS',
    'DRT_LF_RC_RATIO_MIN',
    'DRT_MIN_EFFECTIVE_BINS',
    'DRT_PEAK_PROMINENCE_THRESHOLD',
    'GMM_PEAK_HEIGHT_FACTOR',
    'GMM_N_COMPONENTS_RANGE',

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
