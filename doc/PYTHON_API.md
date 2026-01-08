# Python API

**Current version:** v0.10.2

EIS Analysis Toolkit can be used as a Python library for integration into custom scripts and workflows.

## Installation

```bash
# System packages
sudo apt install python3-numpy python3-scipy python3-matplotlib

# Python packages (for GMM peak detection)
pip install scikit-learn --break-system-packages
```

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
    estimate_rinf_with_inductance,
    # DRT
    calculate_drt,
    # Fitting
    fit_equivalent_circuit,
    fit_circuit_multistart,  # Multi-start optimization
    fit_circuit_diffevo,     # Differential evolution (global optimization)
    FitResult,               # Dataclass with fit results
    MultistartResult,        # Dataclass with multi-start results
    DiffEvoResult,           # Dataclass with DE results
    # Voigt analysis
    analyze_voigt_elements,
    format_voigt_report,
    # Oxide analysis
    analyze_oxide_layer,
    # Circuit elements for manual building
    R, C, Q, L, W, Wo, K,
)

# 1. Load data
frequencies, Z = load_data('data.DTA')

# 2. Kramers-Kronig validation (Lin-KK)
kk_result = kramers_kronig_validation(frequencies, Z)
print(f"Lin-KK: M={kk_result.M}, μ={kk_result.mu:.3f}, noise~{kk_result.noise_estimate:.1f}%")

# 2b. Z-HIT validation (non-parametric, faster)
zhit_result = zhit_validation(frequencies, Z)
print(f"Z-HIT: residuals={zhit_result.mean_residual_mag:.2f}%, noise<={zhit_result.noise_estimate:.1f}%")

# 3. DRT analysis with auto-lambda and GMM
tau, gamma, fig_drt, peaks_gmm, fig_rinf = calculate_drt(
    frequencies, Z,
    auto_lambda=True,
    peak_method='gmm'
)

# 4. Voigt element analysis
voigt_info = analyze_voigt_elements(
    tau, gamma, frequencies, Z, peaks_gmm=peaks_gmm
)
print(format_voigt_report(voigt_info))  # Box-style report

# 5. Manual circuit building based on Voigt analysis
# For example for 2 elements:
circuit = R(100) - (R(926) | C(9.84e-07)) - (R(5037) | C(1.06e-05))

# 6. Fitting (single-start)
result, Z_fit, fig = fit_equivalent_circuit(frequencies, Z, circuit)

# 6b. Or multi-start optimization (more robust)
ms_result, Z_fit, fig = fit_circuit_multistart(
    circuit, frequencies, Z,
    n_restarts=10,
    weighting='proportional'
)
print(f"Improvement: {ms_result.improvement:.1f}%")

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
    log_metadata,
)

# Automatic format detection (.DTA or .csv)
frequencies, Z = load_data('data.DTA')

# Explicit CSV loading
frequencies, Z = load_csv_data('data.csv')

# Parse metadata from Gamry .DTA file
metadata = parse_dta_metadata('data.DTA')
# metadata contains: AREA, VDC, VAC, FREQ, DATE, TIME, NOTES, PSTAT
log_metadata(metadata)  # Log metadata to console

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

### eis_analysis.validation

**Kramers-Kronig validation (Lin-KK):**

```python
from eis_analysis.validation import kramers_kronig_validation, KKResult

result = kramers_kronig_validation(
    frequencies,
    Z,
    mu_threshold=0.85,       # mu metric threshold (default)
    include_L=True,          # Include series inductance
    allow_negative=True      # Allow negative R_i (Lin-KK style)
)

# result: KKResult dataclass
result.M                   # Number of Voigt elements used
result.mu                  # Quality metric (>0.85 = good)
result.Z_fit               # Reconstructed impedance
result.residuals_real      # Real part residuals (fraction)
result.residuals_imag      # Imaginary part residuals (fraction)
result.pseudo_chisqr       # Pseudo chi-squared (Boukamp 1995)
result.noise_estimate      # Estimated noise [%] (Yrjana & Bobacka 2024)
result.inductance          # Fitted inductance [H] (if include_L=True)
result.figure              # matplotlib Figure

# Convenience properties
result.mean_residual_real  # Mean |res_real| [%]
result.mean_residual_imag  # Mean |res_imag| [%]
result.is_valid            # True if residuals < 5%
```

**Z-HIT validation (non-parametric):**

```python
from eis_analysis.validation import zhit_validation, ZHITResult

result = zhit_validation(
    frequencies,
    Z,
    ref_freq=None,           # Reference frequency [Hz] (default: geometric mean)
    quality_threshold=5.0,   # Threshold for quality metric [%]
    optimize_offset=False    # Use weighted offset optimization
)

# result: ZHITResult dataclass
result.Z_mag_reconstructed  # Reconstructed |Z| [Ohm]
result.Z_fit                # Reconstructed complex impedance
result.residuals_mag        # Magnitude residuals [%]
result.residuals_real       # Real part residuals (fraction)
result.residuals_imag       # Imaginary part residuals (fraction)
result.pseudo_chisqr        # Pseudo chi-squared
result.noise_estimate       # Upper bound noise estimate [%]
result.quality              # Quality metric (0-1)
result.ref_freq             # Reference frequency used [Hz]
result.figure               # matplotlib Figure

# Convenience properties
result.mean_residual_real   # Mean |res_real| [%]
result.mean_residual_imag   # Mean |res_imag| [%]
result.mean_residual_mag    # Mean |res_mag| [%]
result.is_valid             # True if mean_residual_mag < 5%
```

**Comparison: Lin-KK vs Z-HIT:**

| Aspect | Lin-KK | Z-HIT |
|--------|--------|-------|
| Method | Parametric (Voigt chain fitting) | Non-parametric (numerical integration) |
| Speed | Slower (iterative optimization) | Faster (single pass) |
| Output | M, mu, inductance | Quality metric, ref_freq |
| Noise estimate | Accurate | Upper bound |

Both methods are complementary. Z-HIT is faster and model-free, Lin-KK provides
more detailed diagnostics (M, mu, inductance).

### eis_analysis.rinf_estimation

**R_inf estimation with inductance compensation:**

```python
from eis_analysis.rinf_estimation import estimate_rinf_with_inductance

R_inf, L, circuit, diagnostics, fig = estimate_rinf_with_inductance(
    frequencies,
    Z,
    max_L_nH=1000.0,     # Maximum reasonable inductance [nH]
    verbose=True,        # Log fitting progress
    plot=True            # Create diagnostic plot
)

# R_inf: High-frequency resistance [Ohm]
# L: Inductance [H]
# circuit: None (legacy placeholder)
# diagnostics: dict with fitting details
# fig: matplotlib Figure with R-L-K fit diagnostics

# diagnostics contains:
#   'n_points_used': number of HF points used
#   'R_squared': coefficient of determination
#   'method': 'zero_crossing', 'polynomial', 'rlk_linear', etc.
#   'R_inf_poly': polynomial extrapolation result (for capacitive data)
#   'R_inf_hf': Re(Z) at highest frequency
```

### eis_analysis.drt

**DRT analysis:**

```python
from eis_analysis.drt import calculate_drt

tau, gamma, fig_drt, peaks_info, fig_rinf = calculate_drt(
    frequencies,
    Z,
    n_tau=100,             # Number of points on tau axis
    lambda_reg=None,       # Regularization parameter (None = auto)
    auto_lambda=True,      # Automatic lambda selection via GCV
    normalize_rpol=False,  # Normalize gamma(tau) by R_pol
    use_rl_fit=False,      # R+L fit for R_inf (deprecated)
    use_voigt_fit=False,   # Voigt fit for R_inf
    peak_method='scipy',   # Peak detection: 'scipy' or 'gmm'
    r_inf_preset=None,     # Preset R_inf value (optional)
    gmm_bic_threshold=10.0 # BIC threshold for GMM (default: 10.0)
)

# tau: Time constant axis [s]
# gamma: Distribution function gamma(tau) [Ohm]
# fig_drt: matplotlib Figure with DRT spectrum
# peaks_info: List of dicts with peak information
# fig_rinf: matplotlib Figure with R_inf fit (if use_voigt_fit=True)
```

**GCV automatic lambda selection:**

```python
from eis_analysis.drt import find_optimal_lambda_gcv, compute_gcv_score

# Find optimal lambda using GCV
lambda_opt, gcv_min, lambdas, gcv_scores = find_optimal_lambda_gcv(
    frequencies, Z, tau,
    lambda_range=(1e-6, 1e2),
    n_lambdas=50
)

# Compute GCV score for a specific lambda
gcv = compute_gcv_score(frequencies, Z, tau, lambda_value)
```

**GMM peak detection:**

```python
from eis_analysis.drt import calculate_drt, gmm_peak_detection, GMM_AVAILABLE

# Check if GMM is available (requires scikit-learn)
if GMM_AVAILABLE:
    tau, gamma, fig, peaks_gmm, _ = calculate_drt(
        frequencies, Z,
        peak_method='gmm'
    )

    # peaks_gmm is a list of dicts, each containing:
    # - 'tau': peak position [s]
    # - 'freq': corresponding frequency [Hz]
    # - 'gamma': peak height [Ohm]
    # - 'weight': relative weight
    # - 'std': peak width (sigma)

    # Or use gmm_peak_detection directly
    peaks = gmm_peak_detection(tau, gamma, n_components_max=5)

# Tune GMM sensitivity with BIC threshold
tau, gamma, fig, peaks, _ = calculate_drt(
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
    R, C, Q, L, W, Wo, K,
    # Main functions
    fit_equivalent_circuit,
    fit_circuit_multistart,
    fit_circuit_diffevo,
    fit_voigt_chain_linear,
    # Result types
    FitResult,
    MultistartResult,
    DiffEvoResult,
)

# Build circuit using operators
# - : series (Z = Z1 + Z2)
# | : parallel (1/Z = 1/Z1 + 1/Z2)
circuit = R(100) - (R(1000) | C(1e-6)) - (R(5000) | C(10e-6))

# Fit (values in circuit serve as initial guess)
result, Z_fit, fig = fit_equivalent_circuit(
    frequencies,
    Z,
    circuit,
    weighting='proportional',      # 'uniform', 'sqrt', 'proportional', 'modulus'
    use_analytic_jacobian=True     # Analytic Jacobian (default, faster)
)

# Weighting types:
#   'uniform':     w = 1          - equal weight for all points
#   'sqrt':        w = 1/sqrt|Z|  - compromise
#   'proportional': w = 1/|Z|     - default, balances relative errors
#   'modulus':     w = 1/|Z|^2    - emphasizes low-impedance region
# See WEIGHTING_AND_STATISTICS.md for detailed guide

# result: FitResult dataclass
# Z_fit: Fitted impedance
# fig: matplotlib Figure with Nyquist and residuals
```

**FitResult dataclass:**

```python
result.params_opt      # Optimal parameters (ndarray)
result.params_stderr   # Standard errors of parameters (ndarray)
result.params_ci_95    # 95% confidence intervals (tuple of ndarray)
result.params_ci_99    # 99% confidence intervals (tuple of ndarray)
result.fit_error_rel   # Relative fit error [%]
result.fit_error_abs   # Absolute fit error [Ohm]
result.quality         # 'excellent', 'good', 'acceptable', 'poor'
result.condition_number # Condition number of Jacobian
result.is_well_conditioned # True if cond < 1e10
result.cov             # Covariance matrix (ndarray or None)
```

**Multi-start optimization (for more robust fit):**

```python
from eis_analysis.fitting import fit_circuit_multistart

ms_result, Z_fit, fig = fit_circuit_multistart(
    circuit,
    frequencies,
    Z,
    n_restarts=10,               # Number of restarts
    scale=2.0,                   # Perturbation = 2 sigma
    weighting='proportional',
    parallel=True,               # Parallel execution
    max_workers=4,
    use_analytic_jacobian=True   # Analytic Jacobian (default)
)

# ms_result: MultistartResult dataclass
ms_result.best_result   # Best FitResult
ms_result.all_results   # All FitResult objects
ms_result.n_starts      # Number of starts
ms_result.n_successful  # Number of successful fits
ms_result.improvement   # Improvement over initial fit [%]
```

**Differential Evolution (global optimization):**

```python
from eis_analysis.fitting import fit_circuit_diffevo

de_result, Z_fit, fig = fit_circuit_diffevo(
    circuit,
    frequencies,
    Z,
    strategy=1,                  # 1=randtobest1bin, 2=best1bin, 3=rand1bin
    popsize=15,                  # Population = popsize * number of parameters
    maxiter=1000,                # Max generations
    tol=0.01,                    # Convergence tolerance
    workers=1,                   # Parallelization (-1 = all CPUs)
    weighting='proportional',
    use_analytic_jacobian=True   # Analytic Jacobian for refinement
)

# de_result: DiffEvoResult dataclass
de_result.best_result    # FitResult after least_squares refinement
de_result.de_result      # Raw result from differential_evolution
de_result.de_error       # Error after DE (before refinement) [%]
de_result.final_error    # Error after refinement [%]
de_result.n_evaluations  # Total number of evaluations
de_result.strategy       # Strategy used
de_result.improvement    # Improvement DE -> refinement [%]
```

**Voigt chain linear fitting:**

```python
from eis_analysis.fitting import fit_voigt_chain_linear

circuit, params = fit_voigt_chain_linear(
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
    weighting='proportional'     # Weighting type
)

# circuit: Circuit object with fitted parameters
# params: list of fitted parameter values
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

```python
# Supported elements for analytic Jacobian:
# R, C, L, Q, W, Wo, K

# For unsupported elements, the system automatically switches to numerical:
result, Z_fit, fig = fit_equivalent_circuit(
    frequencies, Z, circuit,
    use_analytic_jacobian=False  # Force numerical Jacobian
)
```

**Fixed parameters:**

```python
# Fix a parameter using fixed=True
circuit = R(100) - (R(1000, fixed=True) | C(1e-6))

# R(1000) will not change during fit
result, Z_fit, fig = fit_equivalent_circuit(frequencies, Z, circuit)
```

**Available circuit elements:**

| Element | Parameters | Impedance |
|---------|-----------|-----------|
| R(value) | R [Ohm] | Z = R |
| C(value) | C [F] | Z = 1/(jomega*C) |
| L(value) | L [H] | Z = j*omega*L |
| Q(Q, n) | Q, n | Z = 1/(Q*(j*omega)^n) |
| W(sigma) | sigma [Ohm*s^(-1/2)] | Z = sigma*(1-j)/sqrt(omega) |
| Wo(R, tau) | R [Ohm], tau [s] | Warburg bounded |
| K(R, tau) | R [Ohm], tau [s] | Z = R/(1+j*omega*tau) (Voigt) |

**Voigt element analysis from DRT:**

```python
from eis_analysis.fitting import analyze_voigt_elements, format_voigt_report

voigt_info = analyze_voigt_elements(
    tau,                  # From calculate_drt()
    gamma,                # From calculate_drt()
    frequencies,
    Z,
    peaks_gmm=None,       # GMM peaks from calculate_drt() (optional)
    classify_terms=False  # Classify term types (requires GMM)
)

# voigt_info is a dict with keys:
#   'elements': list of dict, each containing:
#       - 'id': int (1-based index)
#       - 'tau': float [s]
#       - 'freq': float [Hz]
#       - 'R': float [Ohm]
#       - 'C': float [F]
#       - 'warnings': list of str
#   'quality': str ('good', 'acceptable', 'uncertain', 'poor')
#   'total_R': float (sum R_i) [Ohm]
#   'R_pol': float (from data) [Ohm]
#   'R_inf': float (from data) [Ohm]
#   'ratio': float (total_R / R_pol)
#   'warnings': list of str (global warnings)
#   'method': str ('gmm' or 'scipy')

# Formatted report (as in CLI)
report = format_voigt_report(voigt_info)
print(report)
# Displays box-style table with tau, f, R, C for each element

# Or process programmatically
for elem in voigt_info['elements']:
    print(f"Element {elem['id']}: tau={elem['tau']:.2e} s, R={elem['R']:.1f} Ohm, C={elem['C']:.2e} F")
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
    fit_result=fit_result    # FitResult from fit_equivalent_circuit()
)

# result is OxideAnalysisResult dataclass:
result.capacitance           # Effective capacitance [F]
result.capacitance_specific  # Specific capacitance [F/cm^2]
result.thickness_nm          # Oxide thickness [nm]
result.element_type          # 'C', 'K', 'Q', or 'estimate'
result.element_R             # Resistance of dominant element [Ohm]
result.element_tau           # Time constant [s]
```

**Formulas:**
```
d = epsilon_0 * epsilon_r / C_specific     (analyze_oxide_layer)
```

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

**R+L fit diagnostics:**

```python
from eis_analysis.visualization import plot_rl_fit_diagnostics

fig = plot_rl_fit_diagnostics(
    frequencies,
    Z,
    R_inf,
    L,
    cutoff_freq
)
# Plots high-frequency R+L fit diagnostics
```

## Usage Examples

### Example 1: Basic DRT Analysis

```python
from eis_analysis import load_data, calculate_drt
import matplotlib.pyplot as plt

# Load data
frequencies, Z = load_data('my_data.DTA')

# DRT with automatic lambda
tau, gamma, fig, peaks, _ = calculate_drt(
    frequencies, Z,
    auto_lambda=True
)

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
    format_voigt_report,
    fit_equivalent_circuit,
    R, C  # Circuit elements
)

# Load and validate
frequencies, Z = load_data('data.DTA')

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
tau, gamma, fig_drt, peaks_gmm, _ = calculate_drt(
    frequencies, Z,
    auto_lambda=True,
    peak_method='gmm'
)

# Voigt analysis
voigt_info = analyze_voigt_elements(
    tau, gamma, frequencies, Z, peaks_gmm=peaks_gmm
)

# Display report
print(format_voigt_report(voigt_info))

# Build circuit manually based on analysis
# For example for 2 detected elements:
circuit = R(voigt_info['R_inf'])
for elem in voigt_info['elements']:
    circuit = circuit - (R(elem['R']) | C(elem['C']))

# Fit
result, Z_fit, fig_fit = fit_equivalent_circuit(frequencies, Z, circuit)

print("\nFitted parameters:")
print(result)  # FitResult object with params_opt, params_stderr, params_ci_95
```

### Example 3: Batch Processing

```python
from pathlib import Path
from eis_analysis import load_data, calculate_drt
import matplotlib.pyplot as plt

data_dir = Path('data')
output_dir = Path('results')
output_dir.mkdir(exist_ok=True)

for dta_file in data_dir.glob('*.DTA'):
    print(f"Processing {dta_file.name}...")

    # Analysis
    frequencies, Z = load_data(dta_file)
    tau, gamma, fig, peaks, _ = calculate_drt(
        frequencies, Z,
        auto_lambda=True,
        peak_method='gmm'
    )

    # Save
    fig.savefig(output_dir / f"{dta_file.stem}_drt.png", dpi=300)
    plt.close(fig)
```

### Example 4: Oxide Layer Analysis

```python
from eis_analysis import load_data, fit_equivalent_circuit, R, C
from eis_analysis.analysis import analyze_oxide_layer

# Load data
frequencies, Z = load_data('oxide_sample.DTA')

# Build circuit using operator overloading
circuit = R(100) - (R(500) | C(1e-6)) - (R(2000) | C(5e-6))

# Fit
fit_result, Z_fit, fig = fit_equivalent_circuit(frequencies, Z, circuit)

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
frequencies, Z = load_data('complex_sample.DTA')

# Complex circuit with Q elements (many parameters)
circuit = R(100) - (R(500) | Q(1e-6, 0.9)) - (R(2000) | Q(1e-5, 0.85))

# Differential evolution - better for global minimum
de_result, Z_fit, fig = fit_circuit_diffevo(
    circuit,
    frequencies,
    Z,
    strategy=1,                  # randtobest1bin (balanced)
    popsize=15,                  # Population size multiplier
    maxiter=1000,                # Max generations
    workers=-1,                  # Use all CPUs
    weighting='proportional',
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
ms_result, Z_fit_ms, fig_ms = fit_circuit_multistart(
    circuit, frequencies, Z,
    n_restarts=20,
    weighting='proportional'
)

print(f"\nMulti-start error: {ms_result.best_result.fit_error_rel:.3f}%")
print(f"DE error: {de_result.final_error:.3f}%")

# Strategy selection:
# strategy=1: randtobest1bin - balanced (default)
# strategy=2: best1bin - faster, may miss global
# strategy=3: rand1bin - thorough, slower
```

### Example 6: Voigt Chain Linear Fitting

```python
from eis_analysis import load_data, fit_voigt_chain_linear

# Load data
frequencies, Z = load_data('data.DTA')

# Voigt chain linear fit (no nonlinear optimization)
circuit, params = fit_voigt_chain_linear(
    frequencies, Z,
    n_per_decade=3,
    auto_optimize_M=True,    # Auto-select number of elements
    mu_threshold=0.85
)

print(f"Circuit: {circuit}")
print(f"Parameters: {params}")

# Compute fitted impedance
Z_fit = circuit.impedance(frequencies, params)
```

## API Reference

For detailed function documentation and parameters, see docstrings in code:

```python
import eis_analysis
help(eis_analysis.calculate_drt)
help(eis_analysis.fit_equivalent_circuit)
help(eis_analysis.estimate_rinf_with_inductance)
```

Or use an IDE with docstring support (VS Code, PyCharm, etc.).

## Links

+ [README.md](../README.md) - Main documentation (CLI)
+ [CHANGELOG.md](../CHANGELOG.md) - Change history
+ [WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md) - Weighting types and statistical metrics guide
