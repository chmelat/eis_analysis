# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

---

## Version 0.13.3 (2026-01-11)

### Code Quality Improvements

- **Refactored numerical algorithms to use scipy/numpy library functions:**
  - `zhit.py`: Replaced manual cumulative integration loops with `scipy.integrate.cumulative_trapezoid()` (10 lines -> 4 lines)
  - `gcv.py`: Replaced manual curvature calculation with `np.gradient()` and vectorized operations (17 lines -> 8 lines)
  - `core.py`: Replaced loop-based regularization matrix construction with `np.fill_diagonal()`
  - `solvers.py`: Replaced Voigt matrix loop with NumPy broadcasting

- **Added comprehensive CLI integration tests:**
  - New `tests/test_cli_integration.py` with 16 end-to-end tests
  - Tests cover: synthetic data, KK validation, DRT analysis, circuit fitting, Voigt chain, file I/O, error handling
  - Improves test coverage for complete analysis workflows

### Files Modified

- `validation/zhit.py` - Refactored cumulative integration
- `drt/gcv.py` - Refactored curvature calculation
- `drt/core.py` - Refactored regularization matrix construction
- `fitting/voigt_chain/solvers.py` - Vectorized Voigt matrix computation
- `tests/test_cli_integration.py` - New integration test suite

### Documentation

- `doc/CODE_DESIGN_ANALYSIS.md` - Code design review report
- `doc/SCIPY_NUMPY_OPTIMIZATION_REPORT.md` - Library optimization analysis

---

## Version 0.13.2 (2026-01-10)

### New Features

- **Gerischer element (G) for reaction-diffusion processes:**
  - Impedance: `Z = sigma / sqrt(1 + j*omega*tau)`
  - Models coupled diffusion with first-order chemical reaction
  - Applications: SOFC cathodes, porous electrodes, MIECs
  - Full support: operator overloading, analytic Jacobian, bounds

### Usage

```python
from eis_analysis.fitting import R, G, fit_equivalent_circuit

# Gerischer element with sigma=100, tau=1ms
g = G(100, 1e-3)

# Circuit: series resistance + Gerischer
circuit = R(10) - G(100, 1e-3)

# Fixed sigma (won't be fitted)
circuit = R(10) - G("100", 1e-3)

# Fit to data
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

### CLI

```bash
# SOFC cathode model
eis data.DTA --circuit "R(10) - G(100, 1e-3)"

# With fixed Gerischer pre-factor
eis data.DTA --circuit 'R(10) - G("100", 1e-3)'
```

### Files Modified

- `fitting/circuit_elements.py` - Added `G` class
- `fitting/jacobian.py` - Added analytic Jacobian for G
- `fitting/bounds.py` - Added bounds for sigma_G, tau_G
- `fitting/__init__.py` - Export G element
- `tests/test_G_element.py` - Comprehensive tests

---

## Version 0.13.1 (2026-01-10)

### Improvements

- **Unified error handling pattern across modules:**
  - `KKResult` dataclass extended with `error` field and `success` property
  - `kramers_kronig_validation()` returns `KKResult(error=...)` instead of `None`
  - `MultistartDiagnostics` extended with `failed_errors` field for tracking failed fits
  - Silent exception swallowing in `multistart.py` replaced with DEBUG logging

### Usage

```python
# New pattern for KK validation
result = kramers_kronig_validation(freq, Z)
if not result.success:
    print(f"Validation failed: {result.error}")

# Access multistart failure diagnostics
if ms_result.diagnostics.failed_errors:
    print(f"Failed fits: {ms_result.diagnostics.failed_errors}")
```

---

## Version 0.13.0 (2026-01-10)

### Breaking Changes

- **`calculate_drt()` now returns `DRTResult` dataclass** instead of 5-element tuple
  - Before: `tau, gamma, fig, peaks, fig_rinf = calculate_drt(...)`
  - After: `result = calculate_drt(...)`
  - Access via: `result.tau`, `result.gamma`, `result.figure`, `result.peaks`, `result.figure_rinf`
  - Additional fields: `result.R_inf`, `result.R_pol`, `result.lambda_reg`, `result.diagnostics`

### Refactoring

- **Major refactoring: Verbose logging moved from core to CLI**
  - Core modules now return structured data (dataclasses) instead of logging
  - All user output handled by CLI layer
  - Clean library usage without side effects

- **New diagnostics dataclasses:**
  - `DRTDiagnostics` - DRT analysis details (condition number, lambda selection, NNLS solution)
  - `RinfEstimate` - R_inf estimation details (method, R_squared, inductance)
  - `LambdaSelection` - Lambda selection details (GCV, L-curve, hybrid)
  - `FitDiagnostics` - Circuit fitting details (optimizer status, covariance info)
  - `MultistartDiagnostics` - Multi-start optimization details
  - `DiffEvoDiagnostics` - Differential evolution details

- **Modules refactored:**
  - `drt/core.py` - Removed 97 logger calls, added structured diagnostics
  - `rinf_estimation/rlk_fit.py` - Removed 48 logger calls
  - `fitting/circuit.py` - Removed 26 logger calls
  - `fitting/diffevo.py` - Removed 33 logger calls
  - `fitting/multistart.py` - Removed 17 logger calls
  - `validation/kramers_kronig.py` - Removed ~20 logger calls

### Improvements

- **CLI output preserved:** All diagnostic output from CLI remains unchanged
  - DE/multistart optimization progress now logged from CLI handlers
  - Fit results with parameters, stderr, and 95% CI displayed correctly
- **New exports from `eis_analysis`:**
  - `DRTResult`, `DRTDiagnostics`
  - `FitDiagnostics`, `MultistartDiagnostics`, `DiffEvoDiagnostics`

### Documentation

- Updated `doc/PYTHON_API.md` with new dataclass API

---

## Version 0.12.1 (2026-01-09)

### Improvements

- **Unified DRTResult dataclass:** Removed duplicate `CLIDRTResult` from CLI module
  - Single source of truth: `eis_analysis.drt.DRTResult`
  - Re-exported from `eis_analysis.cli` for convenience
  - Consistent field names: `figure`, `peaks`, `figure_rinf`
  - Added diagnostic fields: `R_inf`, `R_pol`, `lambda_used`, `reconstruction_error`
  - Added `success` property and `as_tuple()` method for backward compatibility

---

## Version 0.12.0 (2026-01-09)

### Refactoring

- **Major CLI refactoring:** Monolithic `eis.py` (1107 lines) split into modular structure
  - Main `eis.py` reduced to 146 lines (87% reduction)
  - New `eis_analysis/cli/` subpackage with 6 focused modules:
    - `logging.py` - Custom log formatters and setup (131 lines)
    - `parser.py` - Argument parsing with logical grouping (260 lines)
    - `data_handling.py` - Data loading and filtering (181 lines)
    - `handlers.py` - Analysis workflow handlers (602 lines)
    - `utils.py` - Helper functions and dataclasses (174 lines)
    - `__init__.py` - Public API exports (57 lines)

### Improvements

- CLI `--help` now displays arguments in logical groups:
  - Input/Output
  - DRT Analysis
  - Kramers-Kronig Validation
  - Z-HIT Validation
  - Circuit Fitting
  - Voigt Chain Fitting
  - Oxide Layer Analysis
  - Visualization

### For Developers

- CLI components can now be imported from `eis_analysis.cli`:
  ```python
  from eis_analysis.cli import (
      parse_arguments,
      setup_logging,
      run_drt_analysis,
      EISAnalysisError,
  )
  ```
- Each CLI module can be tested independently
- Better separation of concerns improves maintainability

---

## Version 0.11.3 (2026-01-09)

### Improvements

- `lin_kk_native()` now returns `LinKKResult` dataclass instead of 10-element tuple
  - Named access: `result.mu` instead of `result[1]`
  - Properties: `mean_residual_real`, `mean_residual_imag`, `is_valid`
  - Backward incompatible but improves API consistency

---

## Version 0.11.2 (2026-01-09)

### Improvements

- Simplified Z-HIT phase derivative calculation
  - Replaced 9-line manual central differences with `np.gradient()`
  - Cleaner code, same numerical behavior
  - Correctly handles non-equidistant frequency data

---

## Version 0.11.1 (2026-01-09)

### Fixes

- Fixed `--auto-extend` optimization searching unnecessary negative range
  - `extend_decades` works in TIME DOMAIN (extends tau_max toward higher values / lower frequencies)
  - Negative values have no effect (cannot extend toward higher frequencies)
  - Search range changed from `(-max, +max)` to `(0, +max)`
- Fixed optimizer returning boundary value (-4.0) when chi-squared is flat
  - Now prefers values closer to 0 when multiple values have same chi-squared
- Fixed `extend_decades` parameter documentation to clarify time domain behavior

### Documentation

- Added `doc/VALIDATION_METHOD_COMPARISON.md` - comprehensive comparison of Lin-KK vs Z-HIT validation methods
  - Explains when Lin-KK fails (capacitive/inductive data with phase near +/-90 degrees)
  - Documents structural limitations of Voigt chain model
  - Provides recommendations for choosing validation method
- Updated `doc/WEIGHTING_AND_STATISTICS.md`:
  - Fixed incorrect "modulus weighting" comment in chi-squared definition (should be "proportional")
  - Added theoretical derivation for noise estimation formula (constant 5000 = 100^2/2)
- Updated `doc/CODE_ANALYSIS_REPORT.md` to version 0.11.0

---

## Version 0.11.0 (2026-01-08)

### Breaking Changes

- **Weighting names swapped to match common EIS terminology:**
  - `modulus` now means w = 1/|Z| (previously called `proportional`)
  - `proportional` now means w = 1/|Z|^2 (previously called `modulus`)
- Default weighting is now `modulus` (w = 1/|Z|, Lin-KK standard)
- **WARNING:** Existing scripts using `--weighting proportional` or `--weighting modulus` will get DIFFERENT behavior!

### Migration Guide

If you were using:
- `--weighting proportional` (old: 1/|Z|) -> use `--weighting modulus` (new: 1/|Z|)
- `--weighting modulus` (old: 1/|Z|^2) -> use `--weighting proportional` (new: 1/|Z|^2)

The default behavior (1/|Z| weighting) is unchanged, only the name changed from `proportional` to `modulus`.

---

## Version 0.10.2 (2026-01-08)

### Fixes

- Fixed smooth curve plotting for `--voigt-chain` option (was displaying as broken line instead of smooth curve)
- Fixed smooth curve plotting for KK validation Nyquist plot
- Both fixes use 300 interpolated frequency points for smooth fitted curve display

### Documentation

- Added `doc/WEIGHTING_AND_STATISTICS.md` - comprehensive guide to weighting types and statistical metrics
- Covers: weighting types (uniform, sqrt, proportional, modulus), pseudo chi-squared, noise estimation, residuals, fit error, standard errors, confidence intervals
- Includes practical recommendations and troubleshooting guide

---

## Version 0.10.1 (2026-01-08)

### Fixes

- Changed default weighting from `modulus` to `proportional` for circuit fitting
- Testing showed `modulus` (1/|Z|^2) leads to poor parameter estimation for large resistances
- `proportional` (1/|Z|) provides more consistent results across multiple runs

---

## Version 0.10.0 (2026-01-08)

### New Features

- Added Z-HIT (Z-Hilbert Impedance Transform) validation as default data quality test
- Z-HIT provides non-parametric K-K validation using numerical integration (faster than Lin-KK)
- Z-HIT runs by default alongside Lin-KK; disable with `--no-zhit`
- New functions: `zhit_validation()`, `zhit_reconstruct_magnitude()`
- Implementation inspired by pyimpspec library, uses direct integration instead of FFT-based Hilbert transform
- Added pseudo chi-squared metric to Lin-KK validation (Boukamp 1995)
- Added noise estimation from pseudo chi-squared (Yrjana & Bobacka 2024)
- Z-HIT noise estimate shown as upper bound (includes integration approximation error)
- New helper functions: `compute_pseudo_chisqr()`, `estimate_noise_percent()`
- Added automatic extend_decades optimization for KK validation (`--auto-extend` flag)
- Added `--extend-decades-max` CLI parameter to control search range for `--auto-extend` (default: 1.0)
- New functions: `find_optimal_extend_decades()`, `reconstruct_impedance()`
- Added weighted offset optimization for Z-HIT validation (`--zhit-optimize-offset` flag)
- New Z-HIT parameters: `optimize_offset`, `offset_center`, `offset_width`
- New weighting option `modulus` (w=1/|Z|^2) for circuit fitting

### Breaking Changes

- Weighting option `square` renamed to `modulus` (w=1/|Z|^2)
- Default weighting for circuit fitting changed from `sqrt` to `proportional`
- `kramers_kronig_validation()` now returns `KKResult` dataclass instead of tuple
- KKResult provides named access to all validation results: M, mu, Z_fit, residuals, pseudo_chisqr, noise_estimate, inductance, figure
- KKResult includes convenience properties: `mean_residual_real`, `mean_residual_imag`, `is_valid`
- `zhit_validation()` now returns `ZHITResult` dataclass instead of tuple
- ZHITResult provides: Z_mag_reconstructed, Z_fit, residuals_mag, residuals_real, residuals_imag, pseudo_chisqr, noise_estimate, quality, ref_freq, figure
- ZHITResult includes convenience properties: `mean_residual_real`, `mean_residual_imag`, `mean_residual_mag`, `is_valid`
- Z-HIT now computes complex residuals (Re/Im) for better diagnostics alongside magnitude residuals

### Documentation

- Added virtual environment (venv) installation guide with Linux/macOS and Windows examples
- Installation section now recommends venv as the primary approach for isolated environments
- Updated Z-HIT implementation specification (`doc/ZHIT_IMPLEMENTATION_SPEC.md`)
  - Documents design decision: numerical integration vs FFT-based Hilbert transform
  - Includes comparison with pyimpspec implementation
  - Explains noise estimation caveat (upper bound)
- Added weighting explanation to Lin-KK documentation (`doc/LinKK_analysis.md`)

---

## Version 0.9.4 (2026-01-05)

### Documentation

- Added comprehensive time constant analysis report (`time_constant_analysis_report.md`)
- Report compares DRT+GMM vs circuit fitting methods for parameter identification
- Demonstrates effect of regularization parameter λ on DRT spectrum and peak detection
- Provides practical recommendations for choosing analysis methods

### Notes

- Retroactive note: `--gmm-bic-threshold` parameter was added in v0.9.3 but not documented in changelog

---

## Version 0.9.3 (2026-01-05)

### New Features

- Added `--gmm-bic-threshold` CLI parameter for tuning GMM peak detection sensitivity
- Allows users to control BIC threshold for adding components (default: 10.0)
- Lower values detect more peaks, higher values are more conservative

### Fixes

- Fixed missing K element export from `eis_analysis/__init__.py`
- Added K to module imports and `__all__` list

---

## Version 0.9.2 (2026-01-05)

### Fixes

- Fixed missing K element in circuit parser - K element was implemented but not available in CLI parser
- Added K to circuit element imports in `eis.py`
- Added K to `parse_circuit_expression()` safe namespace

### Documentation

- Translated CIRCUIT_PARSER.md to English (was in Czech)
- Updated CIRCUIT_PARSER.md to include K element documentation
- Condensed VERSION_MANAGEMENT.md (from 155 to 56 lines, -64%)
- Removed unnecessary note about impedance package from PYTHON_API.md
- Updated pyproject.toml documentation in VERSION_MANAGEMENT.md

---

## Version 0.9.1 (2026-01-04)

### Fixes

- Fixed `pip install -e .` - CLI command `eis` now works correctly
- Added `py-modules = ["eis"]` to pyproject.toml

### Documentation

- Updated README: all examples use `eis` command instead of `eis.py`
- Added installation instructions with `eis` command usage
- Fixed GitHub repository URL
- Updated CLI help examples

---

## Version 0.9.0 (2025-12-30)

### Initial Release

First version with new versioning scheme.

**Main features:**

- DRT analysis with Tikhonov regularization
- Automatic lambda selection (GCV)
- GMM peak detection in DRT spectrum
- Kramers-Kronig validation
- Circuit fitting with operator overloading syntax
- Voigt chain initial guess (linear regression)
- Oxide layer analysis
- Gamry DTA and CSV format support

---
