# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

---

## Version 0.15.0 (2026-06-29)

### Fixed (DRT math audit F1)

- **GMM peak detection now uses a weighted EM** fitted directly to gamma(tau)
  (`drt/peaks.py`) — addresses `DRT_MATH_AUDIT_2026-06-27` finding F1. The old
  path treated the DRT *density* gamma(tau) as a *sample* by replicating each
  tau bin an integer number of times proportional to gamma
  (`np.repeat(X, max(1, round(gamma/mean)))`). This was methodically incorrect:
  - the BIC scale floated with the artificial replicated count `N`, so
    `bic_threshold` was not comparable across datasets and the log-likelihood
    was inflated by the number of copies;
  - the `max(1, …)` floor gave every bin >=1 sample, injecting a uniform
    background that pulled component means toward the grid center and inflated
    peak widths (sigma);
  - quantization collapsed all gamma < mean/2 to a single level.
  - New helper `_weighted_gaussian_mixture_1d` runs a weighted EM where gamma is
    the per-point weight (no replication, no floor, no quantization). For model
    selection, BIC uses `N = number of frequency measurements` (threaded as
    `n_data` from `calculate_drt`), so `bic_threshold` is comparable across
    datasets. New `tests/test_gmm_weighted.py` (4 tests): mixture recovery,
    BIC scale invariance under gamma scaling, no background-floor width
    inflation, and the n_data penalty.

### Changed / Breaking

- **Removed the optional `scikit-learn` dependency.** sklearn's
  `GaussianMixture` does not support sample weights, so the weighted EM is
  implemented in pure numpy/scipy. GMM peak detection is now always available.
- **Removed the public `GMM_AVAILABLE` flag** (from `eis_analysis` and
  `eis_analysis.drt`). It was always-true gating for the sklearn import, which
  no longer exists.
- `gmm_peak_detection` now returns a `WeightedGMMResult` (with `means_`,
  `covariances_`, `weights_`, `n_components` mirroring the sklearn interface)
  instead of a sklearn `GaussianMixture`, and accepts an optional `n_data`
  argument for the BIC penalty.

## Version 0.14.0 (2026-06-27)

### Tests & Hardening

- **Validated and hardened L-curve corner detection** (`drt/gcv.py`) —
  addresses `DRT_MATH_AUDIT_2026-06-27` finding F6. `find_lcurve_corner`
  picks the point of maximum *signed* curvature, whose sign depends on the
  curve orientation; this was previously untested, so a sign/orientation
  regression would silently pick a wrong lambda.
  - New `tests/test_lcurve_corner.py` (10 tests): pins the curvature formula
    and sign convention (circle of known radius/orientation gives +-1/r),
    localizes the corner on a synthetic L-curve oriented like a real one, and
    checks the orientation assumption (rho up / eta down with lambda) on a
    real DRT L-curve built via the production `_build_drt_matrices`. A
    sign flip (argmax->argmin) makes `test_corner_on_synthetic_L` fail
    (verified by mutation).
  - Confirmed the existing `argmax(positive curvature)` is correct (the
    L-curve corner is a CCW turn -> positive curvature); documented the sign
    convention in the `find_lcurve_corner` docstring.
  - `find_optimal_lambda_hybrid` now sets `diagnostics['corner_at_edge']` and
    logs a warning when the corner lands at the boundary of the searched
    L-range (the true corner likely lies outside the window; widen
    `lcurve_decades`).

### Breaking Changes

- **Removed DRT term-type classification** (CPE/C determination from DRT
  peaks) — addresses `DRT_MATH_AUDIT_2026-06-27` finding F2. An ideal R||C
  (Voigt) element is a Dirac delta in the DRT; the observed peak width is set
  by the regularization parameter lambda, not by the element physics. Width
  thresholds therefore conflated "more regularization" with "CPE", making the
  classification regularization-dependent rather than a measurement. Element
  type (C vs Q) should be determined from circuit fitting instead.
  - Removed CLI flag `--classify-terms` (and the auto-enable-GMM behavior it
    triggered).
  - Removed the `classify_terms` parameter from
    `fitting.analyze_voigt_elements` and the `'classification'` key from each
    returned element dict; `format_voigt_report` no longer prints the
    classification table and always recommends R||C elements.
  - Deleted the `eis_analysis.drt.term_classification` module
    (`classify_peak_type`, `classify_all_peaks`).

## Version 0.13.17 (2026-06-27)

### Code Quality

- **Fixed all mypy errors in the `drt/` package** (`core.py`, `gcv.py`,
  `term_classification.py`) — addresses `AUDIT_2026-06-23` priority 2. No
  behavior change; annotations, `float`/`int` casts, and one `assert`
  documenting an existing invariant.
  - `core.py`: `gamma` is `Optional` but `nnls_result.success == True`
    guarantees a valid array (see `_solve_nnls`); added `assert gamma is not
    None` to narrow the type, clearing the `union-attr`/`operator`/`arg-type`
    errors around normalization, peak detection, and visualization.
  - `gcv.py`: cast `np.linalg.norm` results to `float` (resolves `np.log10`
    overload), cast `np.argmin`/`np.argmax` indices to `int` (resolves
    `__getitem__` overloads), renamed list→ndarray reassignments
    (`gcv_scores`→`gcv_scores_arr`, `rho`→`rho_arr`, …), and annotated the
    hybrid `diagnostics` dict as `Dict[str, Any]`.
  - `term_classification.py`: annotated `type_counts: Dict[str, int]`.

### Bug Fixes

- **Fixed peak-resistance double-counting for overlapping DRT peaks**
  (`drt/peaks.py`, `drt/core.py`) — addresses `DRT_MATH_AUDIT_2026-06-27`
  finding F4. Both peak paths integrated the *total* gamma over each peak's
  window, so the overlap region of adjacent peaks was counted multiple times
  and `sum(R_i) > R_pol`.
  - GMM path (`gmm_peak_detection`): `R_estimate` is now the unit
    decomposition `weight_i * R_pol` (GMM weights sum to 1), guaranteeing
    `sum(R_i) == R_pol` exactly.
  - scipy path (`_estimate_peak_resistance`): the tau axis is partitioned at
    the valleys (gamma minima) between consecutive peaks; trapz additivity
    over the shared boundary node makes `sum(R_i)` equal the spanned-range
    R_pol with no overlap. Removed the now-unused `tolerance` parameter.
  - New regression tests in `tests/test_peak_resistance.py` (4) lock the
    `sum(R_i) == R_pol` property for overlapping peaks on both paths.

- **Narrowed `except Exception` in the DRT core** (`drt/core.py`, `drt/peaks.py`,
  `drt/gcv.py`) — addresses `AUDIT_2026-06-23` finding 2.5 / priority 1. The
  broad handlers swallowed any error (including programmer mistakes such as
  `KeyError`/typos) and silently returned a fallback, hiding real bugs behind
  data-shaped failures.
  - `drt/core.py` `_select_lambda`: both nested handlers around auto-lambda
    selection narrowed to `(np.linalg.LinAlgError, ValueError)`; added
    `logger.debug(..., exc_info=True)` in each branch so a genuine failure
    leaves a traceback instead of vanishing into the `lambda=0.1` fallback.
  - `drt/peaks.py`: GMM fit handler narrowed to
    `(ValueError, np.linalg.LinAlgError)` (sklearn `GaussianMixture.fit`).
  - `drt/gcv.py`: both NNLS handlers (`compute_gcv_score`,
    `compute_lcurve_point`) narrowed to `(RuntimeError, ValueError)`
    (scipy `nnls`).
  - Behavior on expected numerical failures is unchanged; existing
    warning/debug logging is preserved. Verified by the full test suite
    (149 passed) and ruff.

## Version 0.13.16 (2026-06-25)

### Improvements

- **Frequency filtering (`--f-min`/`--f-max`) now logs under its own section
  header** (`cli/data_handling.py`). Previously the filter output appeared
  immediately after the Z-HIT validation block with no header, making it look
  as if the crop applied only to validation. It now prints a
  `"Frequency filtering (analysis range)"` header (matching the other section
  blocks) and reports the resulting analysis frequency range, so it is clear
  the crop applies to all downstream analysis stages (visualization, R_inf,
  DRT, circuit fit) and not to the KK/Z-HIT validation, which run on the full
  spectrum.

## Version 0.13.15 (2026-06-24)

### Refactoring

- **Split `cli/handlers.py` (867 lines) into a `cli/handlers/` package**
  — addresses the audit file-size debt (`AUDIT_2026-06-23` section 4 /
  priority 4); `handlers.py` was the largest file in the project, 73 % over
  the 500-line limit.
  - Submodules by pipeline stage: `validation.py` (KK + Z-HIT), `rinf.py`,
    `drt.py` (DRT + Voigt-from-DRT), `fitting.py` (circuit fitting +
    optimizer diagnostics logging), `oxide.py`. Largest submodule is 418
    lines; all are under the limit.
  - Public API is unchanged: `handlers/__init__.py` re-exports the seven
    `run_*` functions, so `from eis_analysis.cli.handlers import run_*`
    (and `from eis_analysis.cli import run_*`) keep working.
  - Pure code move — no behavior change. Verified by the full test suite
    (149 passed), CLI smoke + circuit-fit/DRT runs, ruff, and mypy
    (unchanged at 63, none new in the package).

## Version 0.13.14 (2026-06-24)

### Improvements

- **`fit_circuit_diffevo` accepts a `seed` argument** (`fitting/diffevo.py`)
  - Threads an optional `seed` into `differential_evolution(seed=...)`.
    Default `None` keeps the previous non-deterministic behavior; passing an
    int makes runs reproducible (used by the new tests).

### Tests

- **Add unit + integration tests for Differential Evolution fitting**
  (`tests/test_diffevo.py`, 18 tests) — addresses the open audit item
  (`AUDIT_2026-06-23` priority 5): the project's most robust optimizer
  previously had no direct tests, only indirect CLI coverage.
  - `DE_STRATEGIES` mapping and unknown-strategy fallback.
  - `_DECostFunction`: free/fixed parameter reconstruction, scalar cost
    (zero at the true parameters), and picklability (its stated purpose for
    `workers > 1`).
  - `fit_circuit_diffevo`: exact recovery on noise-free data, return
    contract, seed reproducibility, strategy/jacobian selection, fixed
    parameters held constant, refinement never worse than DE, and
    diagnostics population. Verified with a mutation check on the
    fixed-parameter reconstruction.

## Version 0.13.13 (2026-06-24)

### Tests

- **Add unit tests for Kramers-Kronig validation** (`tests/test_kramers_kronig.py`,
  24 tests) — addresses the open audit item (`AUDIT_2026-06-23` priority 5):
  the validation core previously had no direct tests, only indirect CLI
  coverage.
  - Pure helpers with exact/analytic checks: `compute_pseudo_chisqr`
    (perfect fit = 0, known weighted value), `estimate_noise_percent`
    (exact formula, monotonicity), `reconstruct_impedance` (single-Voigt vs
    analytic, inductance term, `include_L=False`, low-frequency limit).
  - `lin_kk_native` and `kramers_kronig_validation` on KK-compliant Voigt
    spectra: validity, residual/figure shapes, mu in (0, 1].
  - `KKResult` dataclass contract (empty defaults, validity threshold, error).
  - `find_optimal_extend_decades` range/output.
  - Regression for the v0.13.10 auto-extend fix: on the example spectrum the
    bounded tau grid yields a spurious imaginary residual (>10%) that the
    default extension removes (<5%). Verified with a mutation check.

## Version 0.13.12 (2026-06-24)

### Improvements

- **Show the detected OCV value in the terminal output** (`cli/data_handling.py`)
  - The `OCV data` section previously printed only the point count and
    duration. It now also reports the stabilized open-circuit voltage
    (final filtered `Vf`), the mean over the record, and the drift
    (`|Vf[-1] - Vf[0]|`), e.g.
    `OCV = -446.8 mV (mean -451.7 mV, drift 0.61 mV)`. These values were
    already shown in the OCV plot; this surfaces them without opening the PDF.

## Version 0.13.11 (2026-06-24)

### Bug Fixes

- **Fix misleading DRT "ill-conditioned matrix" warning** (`drt/core.py`)
  - The diagnostic reported the condition number of the *bare* DRT kernel
    `A`, which is intrinsically ill-conditioned for any DRT problem (its
    singular values decay exponentially — precisely why Tikhonov
    regularization is applied). On sparse data this produced an alarming
    `Matrix A is ill-conditioned (1.6e15)` even though the solve was
    numerically sound.
  - DRT never solves `A x = b`; it solves the regularized system
    `[A; sqrt(lambda)*L] x = [b; 0]`, whose condition number on the same
    data is ~450. `DRTDiagnostics.condition_number` (and the CLI warning)
    now report that regularized value, so the false alarm disappears while a
    genuinely ill-conditioned *solve* still triggers the warning.
  - No change to DRT results (gamma, R_pol, peaks, reconstruction error).
  - Note: the separate high reconstruction error / R_pol overestimate seen
    when using `--f-min` to truncate the spectrum is a data-range
    limitation (missing low-frequency arm), not a matrix-conditioning issue.

## Version 0.13.10 (2026-06-24)

### Bug Fixes

- **Fix spurious KK imaginary-part residuals on capacitive/inductive tails**
  (`validation/kramers_kronig.py`, `cli/parser.py`)
  - The Lin-KK test fits a Voigt chain to the **real part** (`fit_type='real'`)
    over a time-constant grid bounded by the measured frequency range
    (`extend_decades=0`). On data with a strong low-frequency capacitive tail
    (or high-frequency inductive tail), that grid cannot reproduce the
    imaginary part near the frequency edges, producing a large, smooth,
    *systematic* imaginary residual (e.g. 22% on `EISPOT-test1.DTA`) even
    though the data is fully KK-compliant — as independently confirmed by
    Z-HIT (~0.4%). The artifact is the edge effect Schönleber et al. (2014)
    address via tau-range extension.
  - Fix: `--auto-extend` (extend_decades optimization) is now **on by
    default**; the library default `kramers_kronig_validation(...,
    auto_extend_decades=True)` matches. With the extended grid the imaginary
    residual collapses to match Z-HIT (0.41% on the same file). Use
    `--no-auto-extend` to restore the previous behavior.

- **Validate the full measured spectrum, not the filtered subset** (`eis.py`)
  - KK and Z-HIT validation now run on the full loaded data, *before* the
    `--f-min`/`--f-max` frequency filter is applied. The filter applies only
    to the analysis stages (R_inf, DRT, circuit fit). Previously, filtering
    first (e.g. `--f-min 1000`) truncated the spectrum that KK then validated;
    because KK is an integral relation over all frequencies, removing the
    low-frequency arm produced a spurious imaginary-part residual
    (18.6% on the filtered `EISPOT-test1.DTA` vs 0.41% on the full spectrum).

### Improvements

- **Cleaner KK terminal output** (`cli/handlers.py`)
  - The Kramers-Kronig section now prints a header and a summary block
    (M, mu, extend_decades, mean residuals, pseudo chi^2, noise estimate,
    data-quality verdict), mirroring the Z-HIT validation section.

## Version 0.13.9 (2026-06-23)

### Bug Fixes

- **Fix wrong `best_start_index` in parallel multi-start** (`fitting/multistart.py`)
  - In parallel mode (`parallel=True`), the restart that produced the best
    fit was identified via `all_errors.index(best_error)`. But `all_errors`
    is filled in *completion* order (`as_completed`) and interleaved with
    `None` for failed fits, so its positions do not correspond to restart
    numbers. The reported `best_start_index` (and the CLI line
    `Best start: #N`) was therefore wrong whenever a non-initial restart won
    and did not happen to complete first.
  - The fit result itself was never affected — only the diagnostic label.
  - Fix: track each successful fit's `start_idx` in a `result_indices` list
    aligned with `all_results`, and read the winning index directly. Robust
    to completion-order shuffling and to failed restarts.
  - Also annotated the local `all_errors` as `List[Optional[float]]`,
    clearing the related mypy `arg-type`/`append` errors.
  - Regression tests: `tests/test_multistart_best_index.py` (parallel
    completion-order, sequential, initial-fit-wins).

### Tests

- **Add IO module test suite** (`tests/test_io_data_loading.py`, 27 tests)
  - `io/data_loading.py` (650 lines) previously had zero direct unit tests
    despite tight coupling to user file formats (audit finding 2.4). New
    suite covers `read_gamry_native`/`load_data`, `load_csv_data`,
    `parse_dta_metadata`, and `parse_ocv_curve`.
  - Synthetic in-memory `.DTA`/CSV fixtures exercise edge cases
    (European decimals, `EXPERIMENTABORTED` truncation, malformed/non-finite
    rows, `MIN_DATA_POINTS` validation, CSV delimiter/column auto-detect and
    positional fallback, metadata defaults). Smoke tests anchor the parsers
    to the real export `example/EISPOT-test1.DTA` (now tracked).
  - Detection power verified via mutation (disabling European-decimal
    conversion fails the suite).

### Improvements

- **Harden IO module** (`io/data_loading.py`) — audit finding 2.4
  - Narrowed all three `except Exception` to `except OSError` (file open in
    `read_gamry_native`, `parse_ocv_curve`, `parse_dta_metadata`). Missing/
    unreadable files are still handled gracefully (`FileNotFoundError`/
    `PermissionError` are `OSError` subclasses), but a genuine parsing bug in
    `parse_dta_metadata` now surfaces instead of being silently turned into
    partial metadata. Per-field `(ValueError, IndexError)` guards inside the
    loops are unchanged (they correctly skip malformed data rows).
  - Annotated the `metadata` dict as `Dict[str, Any]`, clearing all 12 mypy
    errors in the module (project total 75 → 63).
  - Two regression tests lock the narrowed-except behavior
    (`test_metadata_missing_file_returns_defaults`,
    `test_ocv_missing_file_returns_none`).

- **Fix implicit-Optional annotations on `fixed_params`/`full_initial_guess`**
  (`fitting/diffevo.py`, `fitting/jacobian.py`, `fitting/covariance.py`)
  - These parameters default to `None` and `None` is a reachable value
    (a circuit without `get_all_fixed_params` leaves `fixed_params=None`),
    but the annotations declared them as non-optional `list`/`List[bool]`.
    All consumers already guard `is None`, so there was no runtime bug — but
    the false annotation disabled mypy's ability to catch a real None misuse
    and invited a future refactor to drop the guard. Annotations now read
    `Optional[...]`. No behavior change. Clears the related mypy errors.

---

## Version 0.13.8 (2026-04-26)

### Improvements

- **Suppress meaningless CIs for parameters at bounds** (`fitting/bounds.py`,
  `fitting/circuit.py`, `fitting/diffevo.py`, `cli/handlers.py`)
  - When a fit pushes a parameter to its lower or upper bound, the
    Jacobian-based confidence interval is invalid (the local quadratic
    approximation assumes an interior optimum). Previously the CLI still
    printed the numerical CI, sometimes spanning negative values for a
    resistance, masking that the parameter was unidentified.
  - New helper `classify_bound_status(value, lower, upper)` in
    `fitting/bounds.py` (same threshold as `check_bounds_proximity`:
    1 decade on log scale, 1% of range on linear scale).
  - `FitResult` carries a new `bound_status: Optional[List[str]]` field
    (per-parameter: '', 'lower', 'upper', or 'fixed').
  - `_log_fit_result` replaces the CI line with
    `[at <lower|upper> bound — CI not meaningful]` for bound-hit parameters
    and `[fixed]` for fixed parameters.

---

## Version 0.13.7 (2026-04-26)

### Bug Fixes

- **Fixed misleading "Initial guess" log line in DE optimizer**
  (`cli/handlers.py`, `fitting/diffevo.py`)
  - `_log_diffevo_diagnostics` read `result.circuit.get_all_params()` after
    `fit_circuit_diffevo` had already called `circuit.update_params(...)`,
    so the line printed the **final fit values** disguised as the initial
    guess (the two were always identical, regardless of how far DE moved).
  - `DiffEvoDiagnostics` now carries an `initial_guess` field captured
    before DE runs; the handler reads from there.
  - For a default `R()-(R()|Q())-(R()|Q())` circuit the line now correctly
    shows `[100, 100, 1e-4, 0.8, 100, 1e-4, 0.8]` (constructor defaults).

---

## Version 0.13.6 (2026-04-26)

### Bug Fixes (Z-HIT validation)

Based on external code review of `validation/zhit.py` (see
`doc/ZHIT_AUDIT_2026-04-26.md`):

- **Output arrays now match user input order** (`validation/zhit.py`)
  - `zhit_validation` previously returned arrays sorted by ascending frequency
    regardless of input order, causing element-wise mismatch when users plotted
    residuals against their own (often descending) frequency arrays.
  - Inverse permutation now applied to `Z_mag_reconstructed`, `Z_fit`,
    `residuals_mag`, `residuals_real`, `residuals_imag` before returning.
- **Phase unwrap added** (`validation/zhit.py`)
  - `np.unwrap` on `np.arctan2` output prevents 2*pi jumps at the [-pi, pi]
    boundary from spiking the `np.gradient` derivative used in the second-order
    correction. Relevant for inductive systems and noisy data near the wrap point.

### Improvements (Z-HIT validation)

- **Stratified quality labels** (`validation/zhit.py`)
  - New `ZHITResult.quality_label` property: excellent (<0.5%), good (<1.0%),
    acceptable (<2.5%), marginal (<5.0%), poor (>=5.0%).
  - Replaces the binary "good / may contain artifacts" log message that called
    5% residuals "good" — misleading given KK-clean data sits at 0.1-0.5%.
- **`ZHITResult.is_valid` honors `quality_threshold`** (`validation/zhit.py`)
  - Previously hardcoded `< 5.0` regardless of the user-supplied
    `quality_threshold`. New `quality_threshold` field stores the value.
- **Data-driven default for offset window center** (`validation/zhit.py`)
  - `_calculate_offset_weighted` / `zhit_reconstruct_magnitude` /
    `zhit_validation` now accept `offset_center=None` (new default), in which
    case the Gaussian window is centered at `median(log10(frequencies))`.
  - Old fixed default of `1.5` (~31.6 Hz) fell outside the spectrum for
    low-frequency scans (e.g. mHz corrosion measurements), producing
    near-zero weights.

### Documentation

- Added `doc/ZHIT_AUDIT_2026-04-26.md` with full review and remediation plan.

---

## Version 0.13.5 (2026-04-24)

### Bug Fixes

- **Fixed swapped arguments in `log_separator` calls** (`cli/handlers.py`)
  - 7 call sites used `log_separator("=", 50)` instead of `log_separator(50, "=")`.
  - Latent bug: output was correct thanks to Python's commutative `str * int`,
    but would break on any future refactor using `length` as an integer.
- **Fixed incorrect `Optional[KKResult]` return type** in `validation/kramers_kronig.py`
  - Function always returns a `KKResult` (with `error` set on failure), never `None`.
  - Callers accessing `result.success` directly were type-unsafe per the old annotation.
- **Added `None`-narrowing in DRT visualization** (`drt/core.py`)
  - `peaks_result` access inside `if use_gmm:` guarded with explicit assert.

### Infrastructure

- **Added mypy job to CI** (`.github/workflows/ci.yml`)
  - Runs on Python 3.12, `continue-on-error: true` (reports without blocking).
- **Added `[tool.mypy]` section to `pyproject.toml`** with `ignore_missing_imports = true`
  - Silences 43 noisy stub-missing errors from matplotlib/scipy.

### Cleanup

- Removed stale artifacts from repository root (all previously untracked):
  `out*.txt`, `test_G_element.png`, `zry-3d_*.pdf`, `mod2c.sh~`,
  `doc/CODE_REVIEW_REPORT.md~`.

### Documentation

- Added project audit (`doc/AUDIT_2026-04-24.md`) with priorities and remediation plan.

---

## Version 0.13.4 (2026-02-07)

### Code Quality

- **Fixed 7 lint errors** (unused imports and variables) left over from previous refactoring
  - Removed unused imports: `FitDiagnostics`, `MultistartDiagnostics`, `DiffEvoDiagnostics`, `check_bounds_proximity`, `check_parameter_diagnostics`
  - Removed unused variables: `freq_warnings` in `drt/core.py`, `n_params` in `fitting/diffevo.py`

### Infrastructure

- **Added GitHub Actions CI pipeline** (`.github/workflows/ci.yml`)
  - Lint job: ruff check on every push/PR
  - Test job: pytest on Python 3.9 and 3.12

### Documentation

- Added comprehensive code review report (`doc/CODE_REVIEW_2026.md`)

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
