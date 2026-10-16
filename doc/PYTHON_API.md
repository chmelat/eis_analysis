# Python API

EIS Analysis Toolkit can be used as a Python library for integration into custom scripts and workflows.

**Design principle:** Core modules return structured data (dataclasses) without logging. All user output is handled by CLI layer. This makes library usage clean and predictable.

## Installation

```bash
# System packages
sudo apt install python3-numpy python3-scipy python3-matplotlib
```

GMM peak detection is implemented on numpy/scipy; it needs no extra package.

## Basic Workflow

```python
from eis_analysis import (
    # I/O
    load_data,
    load_csv_data,
    parse_dta_metadata,
    generate_synthetic_data,
    # Validation
    kramers_kronig_validation,
    zhit_validation,
    # R_inf estimation
    estimate_rinf,
    RinfResult,
    # DRT
    calculate_drt,
    DRTResult,               # Dataclass with DRT results
    DRTDiagnostics,          # Dataclass with DRT diagnostics
    # Fitting
    fit_equivalent_circuit,
    fit_circuit_multistart,  # Multi-start optimization
    fit_circuit_diffevo,     # Differential evolution (global optimization)
    FitResult,               # Dataclass with fit results
    FitDiagnostics,          # Dataclass with fit diagnostics
    MultistartResult,        # Dataclass with multi-start results
    MultistartDiagnostics,   # Dataclass with multi-start diagnostics
    DiffEvoResult,           # Dataclass with DE results
    DiffEvoDiagnostics,      # Dataclass with DE diagnostics
    # Voigt analysis
    analyze_voigt_elements,
    VoigtSuggestion,        # Dataclass with the Voigt suggestion
    # Oxide analysis
    analyze_oxide_layer,
    # Circuit elements for manual building
    R, C, Q, L, W, Ws, Wo, K, CC, DQ, YG,
)
from eis_analysis.fitting import (
    compute_information_criteria,  # AIC/BIC for model selection
    analyze_residuals,       # Are the residuals noise, or missing structure?
    ResidualDiagnostics,     # Dataclass with residual shape tests
)

# 1. Load data
data = load_data('data.DTA')
frequencies, Z = data.frequencies, data.Z

# 2. Kramers-Kronig validation (Lin-KK)
kk_result = kramers_kronig_validation(frequencies, Z)
print(f"Lin-KK: M={kk_result.M}, μ={kk_result.mu:.3f}, noise~{kk_result.noise_estimate:.1f}%")

# 2b. Z-HIT validation (non-parametric, faster)
zhit_result = zhit_validation(frequencies, Z)
print(f"Z-HIT: residuals={zhit_result.mean_residual_mag:.2f}%, noise<={zhit_result.noise_estimate:.1f}%")

# 3. DRT analysis with auto-lambda and GMM
drt_result = calculate_drt(
    frequencies, Z,
    auto_lambda=True,
    peak_method='gmm'
)
tau, gamma = drt_result.tau, drt_result.gamma
peaks_gmm = drt_result.peaks
# Significant DRT maxima (scipy path); the fallback when GMM finds nothing
peak_indices = [p['index'] for p in drt_result.diagnostics.scipy_peaks or []]

# 4. Voigt element analysis
suggestion = analyze_voigt_elements(
    tau, gamma, frequencies, Z, peaks_gmm=peaks_gmm, peak_indices=peak_indices
)
for elem in suggestion.elements:
    print(f"tau={elem.tau:.2e} s, R={elem.R:.1f} Ohm, C={elem.C:.2e} F")

# 5. Manual circuit building based on Voigt analysis
# For example for 2 elements:
circuit = R(100) - (R(926) | C(9.84e-07)) - (R(5037) | C(1.06e-05))

# 6. Fitting (single-start)
result, Z_fit = fit_equivalent_circuit(frequencies, Z, circuit)

# 6b. Or multi-start optimization (more robust)
ms_result, Z_fit = fit_circuit_multistart(
    frequencies, Z, circuit,
    n_restarts=10,
    weighting='modulus'
)
print(f"Improvement: {ms_result.improvement:.1f}%")
result = ms_result.best_result   # a FitResult like the single-start one

# 7. Fit results
print(f"Parameters: {result.params_opt}")
print(f"Std errors: {result.params_stderr}")
print(f"95% CI: {result.params_ci_95}")
print(f"Fit error: {result.fit_error_rel:.2f}%")
print(f"Quality: {result.quality}")
```

## Modules and Functions

### eis_analysis.io

**Data loading:**

```python
from eis_analysis.io import (
    load_data,
    load_csv_data,
    parse_dta_metadata,
    parse_ocv_curve,
)

# load_data reads Gamry .DTA only; CSV goes through load_csv_data (load_data
# on a CSV raises "No ZCURVE section found").
# Both loaders return a LoadResult: the spectrum, the DTA metadata (None
# for CSV), and any caveat about the data. Nothing is printed - a caveat
# such as a truncated sweep is data you can act on or ignore. A LoadResult
# does not unpack: freq, Z = load_data(...) raises TypeError.
result = load_data('data.DTA')
frequencies, Z = result.frequencies, result.Z
for w in result.warnings:
    print(w)

# Gamry files recorded with the THD option also carry per-point THD
# (fraction, aligned with frequencies; None when the columns are absent,
# always None for CSV)
result.current_thd, result.voltage_thd

# CSV loading (columns freq, Z_re, Z_im; see README for the variants)
result = load_csv_data('data.csv')

# Parse metadata separately (load_data already puts it on the result)
metadata = parse_dta_metadata('data.DTA')
# metadata contains: AREA, VDC, VAC, FREQ, DATE, TIME, NOTES, PSTAT

# Parse OCV (Open Circuit Voltage) curve from DTA file
ocv_data = parse_ocv_curve('data.DTA')
# ocv_data contains: 'time', 'voltage' arrays (or None if not present)
```

**Synthetic data:**

```python
from eis_analysis.io import generate_synthetic_data

# Circuit: Rs - (R0||Q0) - (R1||Q1) where Q is CPE
frequencies, Z = generate_synthetic_data(
    Rs=10,                     # Series resistance [Ohm]
    R0=1e5,                    # First parallel resistance [Ohm]
    Q0=(1e-6, 0.6),            # First CPE: (Y0 [S*s^n], n)
    R1=8e5,                    # Second parallel resistance [Ohm]
    Q1=(3e-5, 0.43),           # Second CPE: (Y0 [S*s^n], n)
    noise=0.01                 # Noise level (0.01 = 1%)
)
```

**Fit result export:**

```python
from eis_analysis.io import fit_result_record, save_fit_result
from eis_analysis.cli.utils import parse_circuit_expression

expr = 'R(10) - (R(1000) | Q(1e-5, 0.9))'
result, Z_fit = fit_equivalent_circuit(frequencies, Z, parse_circuit_expression(expr))

# The record is a plain dict; context is stored as given (input file, ...)
record = fit_result_record(frequencies, Z, Z_fit, result, weighting='modulus',
                           context={'input': 'data.DTA'})

# Writes results_fit.csv, then results_fit.json; returns both paths.
# A record that does not serialize raises before anything is written.
json_path, csv_path = save_fit_result('results_fit', record, frequencies, Z, Z_fit)
```

The record holds the fitted circuit (`str(result.circuit)`), every parameter (`value`, `stderr`,
`ci95`, `significance`, `status`), the metrics including dof and AIC/BIC,
the covariance matrix and the warnings. Values without meaning are `null`:
the uncertainty of a fixed or at-bound parameter, dof and AIC/BIC of the
linear Voigt chain (no free-parameter count), and any inf or NaN, so the
file is strict JSON. The CSV has `freq_Hz, Z_real_Ohm, Z_imag_Ohm,
Z_fit_real_Ohm, Z_fit_imag_Ohm` and loads back with `load_csv_data` (the
fit columns are skipped).

The circuit parses back with `parse_circuit_expression`, fixed values
included exactly; free values are rounded to 4 digits there, so take them
from `parameters`.
The fitted model is that expression with the stored values:

```python
import json
record = json.load(open('results_fit.json'))
circuit = parse_circuit_expression(record['circuit'])
values = [p['value'] for p in record['parameters']]
Z_model = circuit.impedance(frequencies, values)        # the fitted curve
result2, _ = fit_equivalent_circuit(f2, Z2, circuit,    # or a start for
                                    initial_guess=values)  # the next spectrum
```

An oxide result has a record of its own (`results_oxide.json` in the CLI):

```python
from eis_analysis.io import oxide_result_record, save_oxide_result

json_path = save_oxide_result('results_oxide',
                              oxide_result_record(oxide, context={'input': 'data.DTA'}))
```

### eis_analysis.validation

**Kramers-Kronig validation (Lin-KK):**

```python
from eis_analysis.validation import kramers_kronig_validation, KKResult

result = kramers_kronig_validation(
    frequencies,
    Z,
    mu_threshold=0.85,       # mu metric threshold (default)
    max_M=50,                # Upper limit on Voigt elements
    auto_extend_decades=True,  # Try extending the tau grid past the data
    extend_decades_range=(0.0, 1.0),
    include_C=False          # Series capacitance (blocking low-f tail)
)
# The model always has a series L and allows negative R_i (Lin-KK).

# result: KKResult dataclass
result.M                   # Number of Voigt elements used
result.M_lower             # First M the mu search tried (chi^2 plateau)
result.mu                  # Lin-KK stop value (< mu_threshold on normal
                           # termination; not a quality metric - use residuals).
                           # With extend_decades > 0 the returned model's own
                           # mu is at least this value
result.Z_fit               # Reconstructed impedance
result.residuals_real      # Real part residuals (fraction)
result.residuals_imag      # Imaginary part residuals (fraction)
result.pseudo_chisqr       # Pseudo chi-squared (Boukamp 1995)
result.noise_estimate      # Estimated noise [%] (Yrjana & Bobacka 2024)
result.noise_estimate_median  # The same from the median residual: a few large
                           # residuals do not move it. Above 2x apart a
                           # warning says a few points dominate noise_estimate
result.sparse_fit          # M > N/2: the chain passes through most points,
                           # the median is meaningless (no warning, and no
                           # relative criterion in find_outliers)
result.inductance          # Fitted series inductance [H]
result.capacitance         # Fitted series capacitance [F] (if include_C=True)
result.elements            # Fitted elements [R_s, R_1, ..., R_M, L]
result.tau                 # Time constants [s]
result.warnings            # Caveats (max_M reached, extension rejected,
                           # a failed test on a capacitive low-frequency end
                           # without include_C, noise dominated by a few points)
result.capacitive_end      # That last case: residuals and noise_estimate
                           # describe the Voigt chain, not the data
result.with_series_c       # ...so the test is rerun with include_C=True:
                           # a KKResult whose residuals and noise describe
                           # the data (None if the rerun raised)
result.series_c_passes     # Its verdict: False = drift or another violation
                           # beyond the end's shape; None = no rerun

# Convenience properties
result.mean_residual_real  # Mean |res_real| [%]
result.mean_residual_imag  # Mean |res_imag| [%]
result.n_above_threshold   # Points with max(|res_real|, |res_imag|) > 5%
result.is_valid            # True if at most 5% of the points are above 5%

# Figure (the validation itself builds none)
from eis_analysis.visualization import plot_kk_validation
fig = plot_kk_validation(frequencies, Z, result, flagged_frequencies=[])
```

**Z-HIT validation (non-parametric):**

```python
from eis_analysis.validation import zhit_validation, ZHITResult

result = zhit_validation(
    frequencies,
    Z,
    quality_threshold=5.0,   # Threshold on the mean |res_mag| for is_valid [%]
)

# result: ZHITResult dataclass
result.Z_mag_reconstructed  # Reconstructed |Z| [Ohm]
result.Z_fit                # Reconstructed complex impedance
result.residuals_mag        # Magnitude residuals [%]
result.residuals_real       # Real part residuals (fraction)
result.residuals_imag       # Imaginary part residuals (fraction)
result.pseudo_chisqr        # Pseudo chi-squared
result.noise_estimate       # Upper bound noise estimate [%]

# Convenience properties
result.success              # False if the reconstruction failed (empty arrays)
result.mean_residual_real   # Mean |res_real| [%]
result.mean_residual_imag   # Mean |res_imag| [%]
result.mean_residual_mag    # Mean |res_mag| [%]
result.is_valid             # True if mean_residual_mag < 5% - a global
                            # verdict: drift at the low-frequency end raises
                            # a few points and is diluted in the mean (25 %
                            # below 0.1 Hz passed with a 4.7 % mean). Read
                            # residuals_mag, or the Lin-KK verdict, which
                            # counts points above 5 %.
result.quality_label        # 'excellent' ... 'poor', graded on the mean

# Figure (the validation itself builds none)
from eis_analysis.visualization import plot_zhit_validation
fig = plot_zhit_validation(frequencies, Z, result)
```

**Per-point check (both validations):**

```python
from eis_analysis.validation import find_outliers

report = find_outliers(frequencies, kk_result, zhit_result, max_residual=5.0)
report.points    # OutlierPoint list, worst first: frequency [Hz],
                 # residual_kk and residual_zhit [%] (None where that method
                 # did not contribute), worst [%], methods ('KK', 'Z-HIT',
                 # 'KK+Z-HIT')
report.skipped   # Methods left out: median residual above max_residual (the
                 # spectrum fails as a whole), or Lin-KK with capacitive_end
                 # and no with_series_c (with it, its residuals are used)
report.frequencies_flagged_by('kk')   # [Hz], ascending: 'kk' flagged by
                 # Lin-KK, 'both' by both, 'any' all - what the CLI's
                 # --exclude-flagged leaves out; drop them from the arrays
                 # with ~np.isin(frequencies, ...) and refit
```

A point must also exceed four times its method's median residual (README,
Data quality), so the list can hold fewer KK points than
`kk_result.n_above_threshold` counts; the verdict is `is_valid`, the list
says where. A Lin-KK point over 8 times the median Lin-KK residual is listed
under `max_residual` as well (isolated bad points on a low-noise spectrum),
so the list can also hold KK points `is_valid` does not count.

**Harmonic distortion (Gamry THD option):**

```python
from eis_analysis.validation import thd_check, THD_THRESHOLD

result = load_data('data.DTA')
thd = thd_check(result.frequencies, result.current_thd, result.voltage_thd)
# None when the file has no THD columns; otherwise THDResult:
thd.current.median, thd.current.maximum, thd.current.f_at_max   # fractions, Hz
thd.current.n_above                     # points above THD_THRESHOLD (0.01)
thd.current.n_valid                     # ... out of this many points with a value
thd.current.f_above_min, thd.current.f_above_max
thd.voltage                             # same for the voltage (or None)
thd.warnings                            # e.g. a column present but empty

# Figure: THD [%] vs frequency, both channels, threshold line
from eis_analysis.visualization import plot_thd
fig = plot_thd(result.frequencies, result.current_thd, result.voltage_thd, THD_THRESHOLD)
```

**Comparison: Lin-KK vs Z-HIT:**

| Aspect | Lin-KK | Z-HIT |
|--------|--------|-------|
| Method | Parametric (Voigt chain fitting) | Non-parametric (numerical integration) |
| Speed | Slower (iterative optimization) | Faster (single pass) |
| Output | M, mu, inductance | Quality metric |
| Noise estimate | Accurate | Upper bound |

Both methods are complementary. Z-HIT is faster and model-free, Lin-KK provides
more detailed diagnostics (M, mu, inductance).

### eis_analysis.rinf_estimation

**R_inf estimation, R-L-(R|Q) fit over the top two decades:**

```python
from eis_analysis.rinf_estimation import estimate_rinf, hf_median, hf_r_inf

est = estimate_rinf(frequencies, Z)   # ValueError on shape mismatch / no finite point

# est: RinfResult dataclass
est.R_inf          # value to use [Ohm]: fitted R_s, or the HF upper bound (clipped at 0)
est.R_inf_upper    # HF bound + 3 sigma of the noise at f_hf: R_s cannot exceed it
est.R_inf_range    # (lo, hi) meant to hold R_s, the fit's model error included
                   # (a fitted range missed R_s in 7 of 878 stress-test
                   # spectra, see the docstring); est.L: the L to subtract
                   # with R_inf
est.method         # 'rlq_fit' | 'hf_bound'
est.R_inf_fit      # fitted R_s [Ohm], also when not used (None if no fit ran)
est.R_inf_stderr   # standard error of R_s [Ohm] - identifiability flag, not a CI
est.R_inf_hf       # fallback upper bound: Re(Z) at the highest f with Im(Z) <= 0 [Ohm]
est.f_hf           # its frequency [Hz]; below f_max if the top is inductive
est.fit            # FitResult, params_opt = [R_s, L, R_k, Q, n]; None if no fit ran
est.f_window, est.Z_window  # data of the fit window
est.warnings       # why the fallback was used, dropped non-finite points

R_inf, n = hf_median(frequencies, Z)  # median of Re(Z) over the top points
R_inf, warning = hf_r_inf(frequencies, Z)  # the DRT's default: that median, or
                                      # Re(Z) at the top capacitive point when the
                                      # median lies >5 % above the two top points
                                      # (warning says why; 'hf_bound' in
                                      # diagnostics.rinf.method)

from eis_analysis.visualization import plot_rinf_fit
fig = plot_rinf_fit(est)
```

### eis_analysis.drt

**DRT analysis:**

```python
from eis_analysis.drt import calculate_drt, DRTResult, DRTDiagnostics

result = calculate_drt(
    frequencies,
    Z,
    n_tau=100,             # Number of points on tau axis
    lambda_reg=None,       # Regularization parameter (None = auto)
    auto_lambda=True,      # Automatic lambda selection via GCV
    normalize_rpol=False,  # Normalize gamma(tau) by R_pol
    peak_method='scipy',   # Peak detection: 'scipy' or 'gmm'
    r_inf_preset=None,     # Preset R_inf, e.g. estimate_rinf(...).R_inf;
                           # None = HF median
    gmm_bic_threshold=10.0, # BIC threshold for GMM (default: 10.0)
    lambda_probe=False,     # Peak stability across lambda (optional); results
                            # in result.diagnostics.stability
    weighting='sqrt',       # Data weighting: 'sqrt' (default), 'uniform',
                            # 'modulus' or 'proportional'
    tau_extend_decades=0.0, # Grid past the slow window end: decades or 'auto'
    inductance='auto'       # Series L in the model: True, False or 'auto'
                            # (only when the top decade has Im(Z) > 0)
)

# result: DRTResult dataclass
result.tau                 # Time constant axis [s]
result.gamma               # Distribution function gamma(tau) [Ohm]
result.peaks               # GMM peaks (peak_method='gmm'), list of dicts;
                           # None with 'scipy' (the default) or when GMM finds
                           # none - then the maxima are in
                           # diagnostics.scipy_peaks
result.diagnostics.n_peaks # Peaks reported, whichever list holds them
result.Z_reconstructed     # Impedance rebuilt from gamma, L_series, R_inf [Ohm]
result.R_inf               # High-frequency resistance [Ohm]
result.L_series            # Series inductance [H], 0 when not modeled
result.R_pol               # Polarization resistance [Ohm]
result.lambda_used         # Regularization parameter used
result.diagnostics         # DRTDiagnostics with detailed info

# DRTDiagnostics contains:
result.diagnostics.freq_min           # Minimum frequency [Hz]
result.diagnostics.inductance_used    # Series L column in the model
result.diagnostics.inductance_note    # Why 'auto' decided as it did
result.diagnostics.freq_max           # Maximum frequency [Hz]
result.diagnostics.n_points           # Number of data points
result.diagnostics.n_tau              # Number of tau points
result.diagnostics.condition_number   # Matrix condition number
result.diagnostics.reconstruction_error_rel  # Reconstruction error [%]
result.diagnostics.rinf               # RinfEstimate dataclass
result.diagnostics.lambda_sel         # LambdaSelection dataclass
result.diagnostics.nnls               # NNLSSolution dataclass
result.diagnostics.stability          # With lambda_probe=True, else None:
                                      # StabilityDiagnostics; .peak_stability
                                      # is a list of PeakStability (tau_ref,
                                      # R_ref, persistence, n_probes,
                                      # tau_drift_decades, R_variation_rel,
                                      # verdict 'stable'/'marginal'/'artifact')
result.diagnostics.scipy_peaks        # Significant maxima (None when the
                                      # DRT failed, result.success False),
                                      # list of dicts:
                                      # index, tau, frequency, R_estimate,
                                      # delta_chi2, edge_distance_decades,
                                      # boundary_sensitive, outside_window

# Figure (calculate_drt builds none)
from eis_analysis.visualization import plot_drt
fig = plot_drt(Z, result)
```

**GCV automatic lambda selection:**

```python
from eis_analysis.drt import find_optimal_lambda_gcv, compute_gcv_score

# A, b, L: the weighted DRT system that calculate_drt builds internally.
# Searches DRT_LAMBDA_RANGE = (1e-10, 1e-2) by default; a minimum at an
# end of it is refined one decade further, so lambda_opt may lie in
# 1e-11 - 1e-1 (calculate_drt then warns "at search-range edge"). lambda is
# dimensionless (see doc/GCV_IMPLEMENTATION.md).
lambda_opt, gcv_min = find_optimal_lambda_gcv(A, b, L)

# GCV score for a specific lambda
gcv = compute_gcv_score(lambda_opt, A, b, L)
```

**GMM peak detection:**

```python
from eis_analysis.drt import calculate_drt, gmm_peak_detection

result = calculate_drt(
    frequencies, Z,
    peak_method='gmm'
)
peaks_gmm = result.peaks     # None when GMM finds no peak

# peaks_gmm is a list of dicts, each containing:
# - 'tau_center': peak position [s]
# - 'f_center': corresponding frequency [Hz]
# - 'tau_bounds': (tau_lower, tau_upper) [s]
# - 'log_tau_std': peak width (sigma) [decades]
# - 'weight': fraction of R_pol
# - 'R_estimate': weight * R_pol [Ohm]
# - 'edge_distance_decades', 'boundary_sensitive', 'outside_window':
#   position against the measured tau window

# Or use gmm_peak_detection directly; returns (peaks, model, bic_scores).
# n_data as calculate_drt passes it (BIC penalty per measured point); the
# direct call adds no window keys (edge_distance_decades, ...).
peaks, gmm_model, bic_scores = gmm_peak_detection(
    result.tau, result.gamma, n_data=len(frequencies))

# Tune GMM sensitivity with BIC threshold
result = calculate_drt(
    frequencies, Z,
    peak_method='gmm',
    gmm_bic_threshold=5.0  # More sensitive (lower = more peaks)
)
```

### eis_analysis.fitting

**Equivalent circuit fitting (operator overloading API):**

```python
from eis_analysis.fitting import (
    # Circuit elements
    R, C, Q, L, W, Ws, Wo, K, CC, DQ, YG,
    # Main functions
    fit_equivalent_circuit,
    fit_circuit_multistart,
    fit_circuit_diffevo,
    fit_voigt_chain_linear,
    # Result types
    FitResult,
    VoigtChainFit,
    MultistartResult,
    DiffEvoResult,
)

# Build circuit using operators
# - : series (Z = Z1 + Z2)
# | : parallel (1/Z = 1/Z1 + 1/Z2)
circuit = R(100) - (R(1000) | C(1e-6)) - (R(5000) | C(10e-6))

# Fit (values in circuit serve as initial guess)
result, Z_fit = fit_equivalent_circuit(
    frequencies,
    Z,
    circuit,
    weighting='modulus',      # 'uniform', 'sqrt', 'proportional', 'modulus'
    use_analytic_jacobian=True,    # Analytic Jacobian (default, faster)
    bounds=None               # (lower, upper) per parameter; default: absolute PARAMETER_BOUNDS
)

# Weighting types:
#   'uniform':     w = 1          - equal weight for all points
#   'sqrt':        w = 1/sqrt|Z|  - compromise
#   'modulus':     w = 1/|Z|      - default, balances relative errors
#   'proportional': w = 1/|Z|^2   - emphasizes low-|Z| points (high frequencies of a typical spectrum)
# See WEIGHTING_AND_STATISTICS.md for detailed guide

# result: FitResult dataclass
# Z_fit: Fitted impedance

# Figure: Nyquist and residuals (the fit itself builds none; for multistart
# and DE pass .best_result)
from eis_analysis.visualization import plot_circuit_fit
fig = plot_circuit_fit(frequencies, Z, result)
```

**FitResult dataclass:**

```python
result.params_opt      # Optimal parameters (ndarray)
result.params_stderr   # Standard errors (ndarray; inf = not identifiable)
result.params_ci_95    # 95% confidence intervals (tuple of ndarray;
                       #   log-space (always positive) for scale parameters
                       #   R/C/Q/L/sigma/tau, symmetric for exponent n)
result.params_ci_99    # 99% confidence intervals (same semantics)
result.param_labels    # Parameter labels ['R0', 'R1', 'Q0', 'n0', ...]
result.n_free_params   # Number of freely optimized parameters (k for AIC/BIC;
                       #   fixed parameters excluded)
result.fit_error_rel   # Relative fit error [%]
result.fit_error_abs   # Absolute fit error [Ohm]
result.quality         # 'excellent', 'good', 'acceptable', 'poor'
result.condition_number # cond(J_s^T J_s) of the column-scaled (unit-norm
                        #   columns) Jacobian; scale-invariant
result.is_well_conditioned # True if cond(J_s^T J_s) < 1e10
result.cov             # Covariance matrix (ndarray or None)
result.diagnostics     # FitDiagnostics dataclass
result.all_warnings    # List of all warnings (convenience property)
```

### Time constants of the arcs (`arc_time_constants`)

```python
from eis_analysis.fitting import arc_time_constants

for arc in arc_time_constants(result):   # result of a fit (.best_result for
                                         # multistart and DE)
    arc.elements     # 'R1|C0', 'R2|Q0': the arc by its parameter labels
    arc.tau          # [s]: R*C, or (R*Q)^(1/n) for a CPE (the ZARC peak,
                     #   = R*C_eff with the Hsu-Mansfeld C_eff)
    arc.tau_stderr   # [s], delta method on the covariance (R and C are
                     #   correlated); None without a covariance or with a
                     #   parameter at a bound
    arc.stderr_note  # why tau_stderr is None ('no covariance',
                     #   'n0 near its upper bound', 'covariance not
                     #   usable'), else None
    arc.frequency    # 1/(2*pi*tau) [Hz], where -Z'' peaks
```

An arc is a parallel pair of one R (or G) and one C or Q. A C or Q beside a longer
branch is not: in the Randles `C_dl | (R_ct - Ws)` the DC resistance beside
C_dl includes R_W. A `K` has tau as a parameter.

### Are the residuals noise? (`analyze_residuals`)

The fit error says how large the residuals are, not what shape they have. A
circuit missing an element can fit to a few percent while its residuals move
smoothly across the spectrum, and only the shape distinguishes the two.

```python
from eis_analysis.fitting import analyze_residuals

Z_fit = result.circuit.impedance(freq, list(result.params_opt))
d = analyze_residuals(freq, Z, Z_fit, weighting='modulus')

d.is_systematic            # True when either part fails the runs test
d.real.lag1_autocorr       # ~0 for noise, ->1 for a systematic offset
d.real.runs_p              # p-value of the runs test against the median
d.real.slope, d.real.slope_p   # the trend, per decade, with its p-value
d.real.amplitude, d.real.power # the structure left once the trend is removed
d.real.n_eff               # AR(1) effective sample size (see below)
```

The last two pairs say *what kind* of failure it is, and they are read
together. `abs(slope) * d.window_decades` is the trend's total swing across
the window and is in the same units as `amplitude`, so the ratio says which
shape dominates: a trend means an element is missing or is of the wrong type,
structure surviving the trend means the right elements and too few of them.
`power` is the normalized Lomb-Scargle peak, and white noise stays below
`MIN_PERIODOGRAM_POWER` (0.2).

No period is reported. `residual_structure` explains why: a period only means
something if the residual repeats at a fixed spacing, and measured residuals
generally do not - their extrema spread out towards high frequency while the
amplitude decays, so a fitted period is an amplitude-weighted average of local
spacings and differs between the real and imaginary parts of the same fit.

The real and imaginary parts are tested separately, and the frequencies are
sorted internally - both statistics read consecutive points as neighbours.

`n_eff = n(1-rho)/(1+rho)` is reported but is **not** substituted for `n` in
`compute_information_criteria`; see doc/MODEL_SELECTION_AIC_BIC.md for why.

What the numbers mean and how to act on them:
[RESIDUAL_DIAGNOSTICS.md](RESIDUAL_DIAGNOSTICS.md).

### Is the fit at the noise level? (`compare_fit_to_noise`)

The residual tests ask about the shape of the misfit; this asks about its
size relative to the noise. The fit is turned into a noise level with the
estimator Lin-KK reports its noise with, then divided by the Lin-KK noise of
the same points.

```python
from eis_analysis.fitting import compare_fit_to_noise
from eis_analysis.validation import kramers_kronig_validation

kk = kramers_kronig_validation(freq, Z)
# On a capacitive low-frequency end the plain run's noise is the Voigt
# chain's misfit (tens of %), and any fit would look adequate against it;
# the rerun with the series C describes the data
ref = kk.with_series_c if kk.capacitive_end and kk.with_series_c is not None else kk
a = compare_fit_to_noise(Z, Z_fit, ref.noise_estimate)

a.fit_noise_percent   # the fit's residual as noise [%]
a.ratio               # fit / max(noise, NOISE_FLOOR_PERCENT); ~1 for a correct model
a.adequate            # False when ratio > FIT_NOISE_RATIO_WARN (2)
a.warnings            # set when inadequate or not comparable
```

Pass the measured `Z` the circuit was fitted to - what Lin-KK saw. A fit to a
Z-HIT reconstruction has no usable reference: the Lin-KK noise does not
describe the reconstruction's own error, and the CLI does not compare it. If
the fit covers only part of the validated spectrum, use the Lin-KK noise of
those points, `ref.noise_estimate_on(mask)` with `mask` selecting them from
the validated frequencies. Non-finite input gives `ratio = nan` and
`adequate = True` with a warning - a missing comparison is no evidence of a
misfit.

### Comparing candidate circuits (AIC/BIC)

A lower fit error does not mean a better model: adding an element almost
always lowers the residual. `compute_information_criteria` charges for each
free parameter, so competing circuits become comparable.

```python
from eis_analysis.fitting import compute_information_criteria

Z_fit = result.circuit.impedance(freq, list(result.params_opt))
rss, aic, bic = compute_information_criteria(
    Z, Z_fit, weighting='modulus', n_free_params=result.n_free_params
)
```

Only differences between models carry meaning (below 2: indistinguishable,
4-7: noticeable, above 10: decisive), and only between models fitted to the
**same data with the same weighting**. BIC penalizes complexity harder than
AIC and is what the CLI uses to pick a winner.

The CLI ranking built on this is available too:

```python
from eis_analysis.cli.handlers.model_comparison import (
    score_candidates, log_comparison
)

scores = score_candidates(
    freq, Z, [('R-(R|C)', result_a), ('R-(R|Q)', result_b)], 'modulus',
    Z_measured=Z, noise_percent=kk.noise_estimate   # optional: fit/KK column
)
log_comparison(scores, len(freq), 'modulus')
best = scores[0].result        # sorted by BIC; failed candidates sort last
scores[0].noise                # NoiseAdequacyResult; None without noise
```

### FitDiagnostics

```python
result.diagnostics.optimizer_status    # Optimizer exit status
result.diagnostics.optimizer_message   # Optimizer message
result.diagnostics.optimizer_success   # True if converged
result.diagnostics.n_function_evals    # Number of function evaluations
result.diagnostics.jacobian_type       # 'analytic' or 'numeric'
result.diagnostics.condition_number    # cond(J_s^T J_s), column-scaled
result.diagnostics.covariance_rank     # Numerical rank of the (scaled) Jacobian
result.diagnostics.covariance_warning  # Warning if ill-conditioned or
                                       #   rank-deficient (names the
                                       #   non-identifiable parameters)
result.diagnostics.params_at_bounds    # List of param indices at bounds
result.diagnostics.warnings            # List of general warnings
```

**Multi-start optimization (for more robust fit):**

```python
from eis_analysis.fitting import fit_circuit_multistart

ms_result, Z_fit = fit_circuit_multistart(
    frequencies,
    Z,
    circuit,
    n_restarts=10,               # Number of restarts
    scale=2.0,                   # Perturbation = 2 sigma
    weighting='modulus',
    parallel=True,               # Parallel execution
    max_workers=4,
    use_analytic_jacobian=True,  # Analytic Jacobian (default)
    rng=None                     # int seed or numpy Generator for reproducible restarts
)

# ms_result: MultistartResult dataclass
ms_result.best_result   # Best FitResult
ms_result.all_results   # All FitResult objects
ms_result.n_starts      # Number of starts
ms_result.n_successful  # Number of successful fits
ms_result.improvement   # Improvement over initial fit, on the weighted SSR [%]
ms_result.diagnostics   # MultistartDiagnostics dataclass

# MultistartDiagnostics contains:
ms_result.diagnostics.n_restarts          # Total restarts
ms_result.diagnostics.n_successful        # Successful fits
ms_result.diagnostics.scale               # Perturbation scale (sigma)
ms_result.diagnostics.perturbation_method # 'covariance', 'stderr', 'log_uniform'
ms_result.diagnostics.weighting           # Weighting type used
ms_result.diagnostics.jacobian_type       # 'analytic' or 'numeric'
ms_result.diagnostics.initial_error       # Initial fit error [%]
ms_result.diagnostics.best_error          # Best fit error [%]
ms_result.diagnostics.best_start_index    # Index of best start
ms_result.diagnostics.all_errors          # List of all errors
```

**Differential Evolution (global optimization):**

```python
from eis_analysis.fitting import fit_circuit_diffevo

de_result, Z_fit = fit_circuit_diffevo(
    frequencies,
    Z,
    circuit,
    strategy=3,                  # 1=randtobest1bin, 2=best1bin, 3=rand1bin (default)
    popsize=15,                  # Population = popsize * number of parameters
    maxiter=1000,                # Max generations
    tol=0.01,                    # Convergence tolerance
    workers=1,                   # Parallelization (-1 = all CPUs)
    weighting='modulus',
    use_analytic_jacobian=True,  # Analytic Jacobian for refinement
    seed=None,                   # int for a reproducible run
    archive_check=True           # Refine early candidates against local
                                 # minima and ambiguous models (~30 extra
                                 # least_squares runs)
)

# de_result: DiffEvoResult dataclass
de_result.best_result    # FitResult after least_squares refinement
de_result.de_result      # Raw result from differential_evolution
de_result.de_error       # Error after DE (before refinement) [%]
de_result.final_error    # Error after refinement [%]
de_result.n_evaluations  # Total number of evaluations
de_result.strategy       # Strategy used
de_result.improvement    # Improvement DE -> refinement, on the weighted SSR [%]
de_result.diagnostics    # DiffEvoDiagnostics dataclass

# DiffEvoDiagnostics contains:
de_result.diagnostics.strategy           # Strategy name
de_result.diagnostics.popsize            # Population size multiplier
de_result.diagnostics.maxiter            # Max iterations
de_result.diagnostics.tol                # Convergence tolerance
de_result.diagnostics.workers            # Number of workers
de_result.diagnostics.weighting          # Weighting type
de_result.diagnostics.jacobian_type      # 'analytic' or 'numeric'
de_result.diagnostics.de_converged       # True if DE converged
de_result.diagnostics.de_iterations      # DE iterations
de_result.diagnostics.de_evaluations     # DE function evaluations
de_result.diagnostics.de_error           # Error after DE [%]
de_result.diagnostics.refined_error      # Error after refinement [%]
de_result.diagnostics.refinement_improved # True if refinement helped
de_result.diagnostics.total_evaluations  # Total evaluations
de_result.diagnostics.archive_checked    # False with archive_check=False
de_result.diagnostics.archive_candidates # Early-generation candidates refined
de_result.diagnostics.archive_used       # True if one of them is the result
# Repairs and ambiguous models are reported in de_result.diagnostics.warnings
```

**Voigt chain linear fitting:**

```python
from eis_analysis.fitting import fit_voigt_chain_linear

chain = fit_voigt_chain_linear(
    frequencies,
    Z,
    n_per_decade=3,           # Time constants per decade
    extend_decades=0.0,       # Extend tau range by N decades
    include_L=True,           # Include series inductance
    fit_type='complex',       # 'complex', 'real', or 'imag'
    prune_threshold=0.01,     # Threshold for removing small R_i
    allow_negative=False,     # Allow negative R_i (Lin-KK style)
    auto_optimize_M=False,    # Auto-optimize M elements using mu metric
    mu_threshold=0.85,        # mu threshold for auto_optimize_M
    max_M=50,                 # Maximum M elements for auto_optimize_M
    weighting='modulus'     # Weighting type
)

# chain is a VoigtChainFit:
#   .circuit          Circuit object with fitted parameters
#   .initial_params   list of fitted parameter values
#   .diagnostics      VoigtChainDiagnostics - the four steps behind it
#                     (tau grid, regression, pruning, final chain), which
#                     is what the CLI prints its section from
circuit, params = chain.circuit, chain.initial_params
```

**DE vs Multi-start:**

| Aspect | Multi-start | Differential Evolution |
|--------|-------------|------------------------|
| Initial guess | Requires approximate | Uses as x0 seed |
| Speed | Faster | Slower |
| Robustness | Medium | High |
| Global minimum | ~85-95% | ~99% |

**When to use:**
- **Multi-start:** You have an approximate initial guess, circuit has 2-5 parameters
- **DE:** Uncertain initial guess, complex circuit, many parameters (>5)

**Analytic Jacobian:**

All fit functions support analytic Jacobian (`use_analytic_jacobian=True`, default), which is faster and more accurate than numerical approximation.

Every built-in element supplies one: R, C, L, G, Q, W, Ws, Wo, Wa, Wat, K, GE, CC, DQ, YG.

A custom element that does not is **not** silently downgraded - the fit raises:

```
RuntimeError: Circuit fitting failed: Analytic Jacobian not implemented for ZZ.
Use numeric Jacobian (use_analytic_jacobian=False).
```

`fit_equivalent_circuit` does guard `make_jacobian_function` with a
`NotImplementedError` fallback, but that call only builds a closure; the
missing derivative is not reached until `least_squares` evaluates it, well
past the guard. So the switch has to be made by hand:

```python
result, Z_fit = fit_equivalent_circuit(
    frequencies, Z, circuit,
    use_analytic_jacobian=False  # required for an element without a branch
)
```

Note that `compute_significance` has no such switch: it returns `None` when
any element in the circuit lacks an analytic derivative, so the `S=` column
disappears from the fit output. Adding a branch in `jacobian.py` is the only
way to get it back.

**Fixed parameters:**

```python
# Fix a parameter by quoting its value
circuit = R(100) - (R("1000") | C(1e-6))

# R("1000") will not change during fit; its stderr is 0 and its CI
# [1000, 1000], and it does not count in result.n_free_params
result, Z_fit = fit_equivalent_circuit(frequencies, Z, circuit)
```

**Available circuit elements:**

| Element | Parameters | Impedance |
|---------|-----------|-----------|
| R(value) | R [Ohm] | Z = R |
| G(value) | G [S] | Y = G, Z = 1/G (conductance; G = 0 is a valid fit result) |
| C(value) | C [F] | Z = 1/(jomega*C) |
| L(value) | L [H] | Z = j*omega*L |
| Q(Q, n) | Q, n | Z = 1/(Q*(j*omega)^n) |
| W(sigma) | sigma [Ohm*s^(-1/2)] | Z = sigma*(1-j)/sqrt(omega) |
| Ws(R, tau) | R [Ohm], tau [s] | Z = R*tanh(u)/u, u = sqrt(j*omega*tau) (finite-length, transmissive boundary: constant concentration) |
| Wo(R, tau) | R [Ohm], tau [s] | Z = R*coth(u)/u (finite-space, reflective boundary: zero flux) |
| Wa(R, tau, gamma) | R [Ohm], tau [s], gamma | Z = R*coth(x)/x, x = (j*omega*tau)^(gamma/2) (anomalous finite-space; gamma = 1 is Wo) |
| Wat(R, tau, gamma) | R [Ohm], tau [s], gamma | Z = R*tanh(x)/x, x = (j*omega*tau)^(gamma/2) (anomalous finite-length; gamma = 1 is Ws) |
| K(R, tau) | R [Ohm], tau [s] | Z = R/(1+j*omega*tau) (Voigt) |
| GE(sigma, tau) | sigma [Ohm], tau [s] | Z = sigma/sqrt(1+j*omega*tau) (Gerischer) |
| CC(C_inf, dC, tau, alpha) | C_inf, dC [F], tau [s], alpha | Z = 1/(j*omega*(C_inf + dC/(1+(j*omega*tau)^(1-alpha)))) (Cole-Cole) |
| DQ(A, n, tau_min, U) | A [Ohm*s^-n], n, tau_min [s], U | truncated CPE: power-law DRT bounded to tau_min..tau_min*e^U |
| YG(C, p, tau) | C [F], p, tau [s] | Z = p/(j*omega*C)*ln[(1+j*omega*tau*e^(1/p))/(1+j*omega*tau)] (Young-Göhr) |

**Voigt element analysis from DRT:**

```python
from eis_analysis.fitting import analyze_voigt_elements

suggestion = analyze_voigt_elements(
    tau,                  # From calculate_drt()
    gamma,                # From calculate_drt()
    frequencies,
    Z,
    peaks_gmm=None,       # GMM peaks from calculate_drt() (drt_result.peaks)
    peak_indices=None     # [p['index'] for p in drt_result.diagnostics.scipy_peaks]
)
# One of peaks_gmm / peak_indices is required (ValueError otherwise); GMM
# peaks win when non-empty. The scipy peaks are the local maxima of gamma that
# pass calculate_drt's significance test (Delta chi^2 >= DRT_PEAK_DCHI2_MIN
# against a shoulder of the nearest taller peak), so the suggestion is built
# from exactly the peaks the DRT reports. R_i is each peak's valley-to-valley
# area: on the scipy path the same as its R_estimate (a GMM component's
# R_estimate is R_pol * weight instead).

# suggestion is a VoigtSuggestion:
#   .elements       list of VoigtElement (.id, .tau, .freq, .R, .C, .warnings)
#   .quality        'good', 'acceptable', 'uncertain' or 'poor' (always
#                   'poor' when sum(R_i) and R_pol differ by more than 2x)
#   .total_R        sum of the element resistances [Ohm]
#   .R_pol, .R_inf  from the data [Ohm]
#   .ratio          total_R / R_pol
#   .method         'gmm' or 'scipy'
#   .n_peaks_raw    peaks detected
#   .n_peaks_valid  peaks the suggestion is built from
#   .warnings       caveats about the analysis (quality follows their count)
#   .excluded_peaks why individual peaks were dropped

# The CLI's report block is cli/handlers/drt.py; a library caller reads
# the fields directly
for elem in suggestion.elements:
    print(f"Element {elem.id}: tau={elem.tau:.2e} s, R={elem.R:.1f} Ohm, C={elem.C:.2e} F")
```

**Note:**
The `analyze_voigt_elements()` function only reports estimated parameters from DRT - it does not generate a circuit automatically. Users build the circuit manually using operator overloading (conservative approach with full control).

### eis_analysis.analysis

**Oxide layer analysis (thickness estimation):**

```python
from eis_analysis.analysis import analyze_oxide_layer

result = analyze_oxide_layer(
    frequencies,
    Z,
    epsilon_r=22.0,          # Relative permittivity (ZrO2 ~ 20-25)
    area_cm2=1.0,            # Electrode area [cm^2]
    fit_result=fit_result,   # FitResult from fit_equivalent_circuit()
    rho_delta_ohm_cm=500.0   # Optional: enables the power-law model for a Q
)

# result is OxideAnalysisResult dataclass:
result.capacitance           # Effective capacitance [F]
result.capacitance_specific  # Specific capacitance [F/cm^2]
result.thickness_nm          # Oxide thickness [nm]
result.element_type          # 'C', 'K', 'Q', 'CC', 'DQ', 'YG', or 'estimate'
result.element_R             # Resistance of dominant element [Ohm]
result.element_tau           # Time constant [s]
result.element_label         # Its first parameter in the fit's param_labels
                             # ('C1', 'Q0'); None without a fit
result.thickness_pl_nm       # Power-law (Hirschorn-Orazem) comparison, or None
```

**Inverse direction (permittivity from known thickness):**

```python
from eis_analysis.analysis import estimate_permittivity

inverse = estimate_permittivity(
    frequencies,
    Z,
    thickness_nm=25.0,       # Known thickness [nm] (SEM/TEM)
    area_cm2=1.0,
    fit_result=fit_result
)

# Also an OxideAnalysisResult, with input and output swapped:
inverse.permittivity         # Estimated relative permittivity
inverse.permittivity_brug    # Brug (2D) comparison, or None
inverse.permittivity_pl      # Power-law comparison (needs rho_delta_ohm_cm), or None
inverse.thickness_nm         # The thickness that was passed in
```

**Formulas:**
```
d = epsilon_0 * epsilon_r / C_specific     (analyze_oxide_layer)
epsilon_r = d * C_specific / epsilon_0     (estimate_permittivity)
```

**Local CPE exponent n(f) (no circuit fit needed):**

```python
from eis_analysis.analysis import local_exponent
from eis_analysis.rinf_estimation import estimate_rinf

est = estimate_rinf(frequencies, Z)
# est.L: the window fit's L when it determined R_inf, else 0;
# est.R_inf_range: R_inf +- max(5 %, 3 stderr, 3x est.R_inf_window_spread) within
#                   0..R_inf_upper, or (0, R_inf_upper) when undetermined
res = local_exponent(frequencies, Z, R_inf=est.R_inf, L=est.L, R_inf_range=est.R_inf_range)

# res is a LocalExponentResult (arrays sorted by ascending frequency):
res.frequencies, res.n       # n(f) = d ln Re(Y) / d ln(omega), NaN where undefined
res.n_uncertainty            # regression stderr and R_inf sensitivity (R_inf_range, default +-5 %)
res.valid                    # points with uncertainty <= 0.02
res.n_min, res.f_n_min       # over the valid points
res.n_max, res.f_n_max
res.span                     # n_max - n_min; > 0.1 adds a warning
res.warnings

from eis_analysis.visualization import plot_local_exponent
fig = plot_local_exponent(res)
```

Y = 1/(Z - R_inf - j omega L). One CPE gives a flat n(f); see the README
section "Local CPE exponent n(f)" for how to read the rest.

### eis_analysis.visualization

**Basic visualization (Nyquist + Bode):**

```python
from eis_analysis.visualization import visualize_data

fig = visualize_data(
    frequencies,
    Z,
    title='EIS spectrum'
)
# Plots: Nyquist diagram, Bode amplitude, Bode phase
```

**OCV visualization:**

```python
from eis_analysis.visualization import visualize_ocv
from eis_analysis.io import parse_ocv_curve

ocv_data = parse_ocv_curve('data.DTA')
if ocv_data is not None:
    fig = visualize_ocv(ocv_data, title='OCV Curve')
```

**DRT:**

```python
from eis_analysis.visualization import plot_drt

result = calculate_drt(frequencies, Z, peak_method='gmm')
fig = plot_drt(Z, result)
# Plots: DRT spectrum, Nyquist reconstruction check;
# with GMM peaks also the per-peak deconvolution and the BIC curve
```

**Circuit fit:**

```python
from eis_analysis.visualization import plot_circuit_fit

result, Z_fit = fit_equivalent_circuit(frequencies, Z, circuit)
fig = plot_circuit_fit(frequencies, Z, result)
# Plots: Nyquist with the fitted curve, residuals vs frequency.
# Multistart / DE: pass ms_result.best_result / de_result.best_result
```

## Usage Examples

### Example 1: Basic DRT Analysis

```python
from eis_analysis import load_data, calculate_drt
from eis_analysis.visualization import plot_drt
import matplotlib.pyplot as plt

# Load data
data = load_data('my_data.DTA')
frequencies, Z = data.frequencies, data.Z

# DRT with automatic lambda
result = calculate_drt(
    frequencies, Z,
    auto_lambda=True
)

# Access results
print(f"R_inf = {result.R_inf:.2f} Ohm")
print(f"R_pol = {result.R_pol:.2f} Ohm")
print(f"lambda = {result.lambda_used:.2e}")
print(f"Found {result.diagnostics.n_peaks} peaks")

# Show figure
plot_drt(Z, result)
plt.show()
```

### Example 2: Complete Workflow with GMM

```python
from eis_analysis import (
    load_data,
    kramers_kronig_validation,
    zhit_validation,
    calculate_drt,
    analyze_voigt_elements,
    fit_equivalent_circuit,
    R, C  # Circuit elements
)

# Load and validate
data = load_data('data.DTA')
frequencies, Z = data.frequencies, data.Z

# Kramers-Kronig validation (Lin-KK)
kk_result = kramers_kronig_validation(frequencies, Z)
print(f"Lin-KK: M={kk_result.M}, μ={kk_result.mu:.3f}")
print(f"  Noise estimate: {kk_result.noise_estimate:.1f}%")
print(f"  Valid: {kk_result.is_valid}")

# Z-HIT validation (faster, non-parametric)
zhit_result = zhit_validation(frequencies, Z)
print(f"Z-HIT: residuals={zhit_result.mean_residual_mag:.2f}%")
print(f"  Noise upper bound: {zhit_result.noise_estimate:.1f}%")

# DRT with GMM
drt_result = calculate_drt(
    frequencies, Z,
    auto_lambda=True,
    peak_method='gmm'
)

# Voigt analysis
suggestion = analyze_voigt_elements(
    drt_result.tau, drt_result.gamma, frequencies, Z,
    peaks_gmm=drt_result.peaks,
    peak_indices=[p['index'] for p in drt_result.diagnostics.scipy_peaks or []]
)

# Build circuit manually based on analysis
circuit = R(drt_result.R_inf)
for elem in suggestion.elements:
    print(f"tau={elem.tau:.2e} s, R={elem.R:.1f} Ohm, C={elem.C:.2e} F")
    circuit = circuit - (R(elem.R) | C(elem.C))

# Fit
result, Z_fit = fit_equivalent_circuit(frequencies, Z, circuit)

print("\nFitted parameters:")
print(result)  # FitResult object with params_opt, params_stderr, params_ci_95
```

### Example 3: Batch Processing

```python
from pathlib import Path
from eis_analysis import load_data, calculate_drt
from eis_analysis.visualization import plot_drt
import matplotlib.pyplot as plt

data_dir = Path('data')
output_dir = Path('results')
output_dir.mkdir(exist_ok=True)

for dta_file in data_dir.glob('*.DTA'):
    print(f"Processing {dta_file.name}...")

    # Analysis
    data = load_data(dta_file)
    frequencies, Z = data.frequencies, data.Z
    result = calculate_drt(
        frequencies, Z,
        auto_lambda=True,
        peak_method='gmm'
    )

    # Save figure (plot_drt raises on a failed DRT; skip it, keep the loop)
    if result.success:
        fig = plot_drt(Z, result)
        fig.savefig(output_dir / f"{dta_file.stem}_drt.png", dpi=300)
        plt.close(fig)

    # Print summary
    print(f"  R_inf={result.R_inf:.1f} Ohm, {result.diagnostics.n_peaks} peaks")
```

### Example 4: Oxide Layer Analysis

```python
from eis_analysis import load_data, fit_equivalent_circuit, R, C
from eis_analysis.analysis import analyze_oxide_layer

# Load data
data = load_data('oxide_sample.DTA')
frequencies, Z = data.frequencies, data.Z

# Build circuit using operator overloading
circuit = R(100) - (R(500) | C(1e-6)) - (R(2000) | C(5e-6))

# Fit
fit_result, Z_fit = fit_equivalent_circuit(frequencies, Z, circuit)

# Oxide layer analysis (estimate thickness from permittivity)
oxide = analyze_oxide_layer(
    frequencies, Z,
    epsilon_r=22.0,
    area_cm2=1.0,
    fit_result=fit_result
)

print(f"Dominant element: {oxide.element_type}")
print(f"Resistance: {oxide.element_R:.1f} Ohm")
print(f"Capacitance: {oxide.capacitance*1e6:.2f} uF")
print(f"Specific cap.: {oxide.capacitance_specific*1e6:.2f} uF/cm^2")
print(f"Oxide thickness: {oxide.thickness_nm:.1f} nm")
print(f"Time constant: {oxide.element_tau*1e3:.2f} ms")
```

### Example 5: Differential Evolution for Complex Circuits

```python
from eis_analysis import (
    load_data,
    fit_circuit_diffevo,
    fit_circuit_multistart,
    R, C, Q
)

# Load data
data = load_data('complex_sample.DTA')
frequencies, Z = data.frequencies, data.Z

# Complex circuit with Q elements (many parameters)
circuit = R(100) - (R(500) | Q(1e-6, 0.9)) - (R(2000) | Q(1e-5, 0.85))

# Differential evolution - better for global minimum
de_result, Z_fit = fit_circuit_diffevo(
    frequencies,
    Z,
    circuit,
    strategy=3,                  # rand1bin (default)
    popsize=15,                  # Population size multiplier
    maxiter=1000,                # Max generations
    workers=-1,                  # Use all CPUs
    weighting='modulus',
    use_analytic_jacobian=True   # Analytic Jacobian for refinement
)

print(f"DE error: {de_result.de_error:.3f}%")
print(f"Refined error: {de_result.final_error:.3f}%")
print(f"Improvement: {de_result.improvement:.1f}%")
print(f"Evaluations: {de_result.n_evaluations}")

# Access FitResult
result = de_result.best_result
print(f"\nOptimal parameters: {result.params_opt}")
print(f"Standard errors: {result.params_stderr}")
print(f"Quality: {result.quality}")

# Comparison with multi-start
ms_result, Z_fit_ms = fit_circuit_multistart(
    frequencies, Z, circuit,
    n_restarts=20,
    weighting='modulus'
)

print(f"\nMulti-start error: {ms_result.best_result.fit_error_rel:.3f}%")
print(f"DE error: {de_result.final_error:.3f}%")

# Strategy selection:
# strategy=1: randtobest1bin - faster, can get trapped in a local minimum
# strategy=2: best1bin - faster, may miss global
# strategy=3: rand1bin - thorough, slower
```

### Example 6: Voigt Chain Linear Fitting

```python
from eis_analysis import load_data, fit_voigt_chain_linear

# Load data
data = load_data('data.DTA')
frequencies, Z = data.frequencies, data.Z

# Voigt chain linear fit (no nonlinear optimization)
chain = fit_voigt_chain_linear(
    frequencies, Z,
    n_per_decade=3,
    auto_optimize_M=True,    # Auto-select number of elements
    mu_threshold=0.85
)

print(f"Circuit: {chain.circuit}")
print(f"Parameters: {chain.initial_params}")
print(f"M chosen: {chain.diagnostics.mu_optimization.M}")

# Compute fitted impedance
Z_fit = chain.circuit.impedance(frequencies, chain.initial_params)
```

## API Reference

For detailed function documentation and parameters, see docstrings in code:

```python
import eis_analysis
help(eis_analysis.calculate_drt)
help(eis_analysis.fit_equivalent_circuit)
help(eis_analysis.estimate_rinf)
```

Or use an IDE with docstring support (VS Code, PyCharm, etc.).

## Links

+ [README.md](../README.md) - Main documentation (CLI)
+ [CHANGELOG.md](../CHANGELOG.md) - Change history
+ [WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md) - Weighting types and statistical metrics guide
