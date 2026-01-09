# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

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
