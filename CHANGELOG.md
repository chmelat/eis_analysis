# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

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
