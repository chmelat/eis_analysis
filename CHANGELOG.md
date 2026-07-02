# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

---

## Version 0.16.16 (2026-07-02)

### Fixed (oxide audit 2026-07-02, finding O3)

- **Silent assumptions in oxide analysis are now visible in the log**
  (`analysis/oxide.py`, thresholds documented in `analysis/config.py`):
  - All candidate capacitive elements (type, R, C/Q, tau) are listed, and
    the selection assumption ("largest R = compact oxide barrier") is
    stated explicitly so the choice can be verified against a
    charge-transfer interpretation.
  - A warning is logged when the dominant element is a CPE with n < 0.8
    (`CPE_N_RELIABLE_MIN`), where a single effective capacitance is not
    well-defined and the thickness estimate may be unreliable.
- **High-frequency fallback (no fitted circuit) is more robust**: the
  capacitance is now the median of `C_i = -1/(omega*Z'')` over the
  capacitive points in the top frequency decade (`HF_ESTIMATE_DECADE_FACTOR`)
  instead of the single highest-frequency point. **Numeric change:** Mode 2
  results can differ slightly from <= 0.16.15 (fitted-circuit results are
  unchanged). New warnings: the fallback always notes that multilayer
  (series) systems yield the series combination of layer capacitances, and
  it warns when the per-point estimates spread by more than max/min = 1.2
  (`HF_C_SPREAD_MAX_RATIO`) across the decade — i.e. `omega*R*C >> 1` does
  not hold and the estimate is unreliable. (A phase-angle check was
  considered and rejected: the phase at the highest frequency is dominated
  by the series resistance and flags data where the estimate is in fact
  exact.) If no capacitive point exists in the top decade, the previous
  single-point behavior is preserved.
  - Regression tests added to `tests/test_oxide.py` (candidate listing,
    n-warning on/off, median accuracy, spread warning on/off, series note,
    inductive-data edge case).

---

## Version 0.16.15 (2026-07-02)

### Fixed (oxide audit 2026-07-02, finding O2)

- **`estimate_permittivity` no longer logs a bogus oxide thickness**
  (`analysis/oxide.py`). It previously reused `analyze_oxide_layer` with a
  dummy `epsilon_r=1.0` to extract the capacitance, so the log contained a
  full analysis block including a meaningless "Oxide thickness" line (and a
  misleading "Oxide layer analysis" header). The shared element-selection
  and capacitance-extraction logic now lives in a private helper
  `_extract_capacitance()`; each public function logs only the quantity it
  actually derives (thickness resp. permittivity). Computed values are
  unchanged for both functions; `analyze_oxide_layer` log output is
  unchanged.
  - Regression tests in new `tests/test_oxide.py` (no thickness line in the
    permittivity log for both the fitted-circuit and the high-frequency
    fallback path, plus thickness/permittivity roundtrip consistency).

---

## Version 0.16.14 (2026-07-02)

### Fixed (oxide audit 2026-07-02, finding O1)

- **CPE conversion formula correctly attributed to Hsu-Mansfeld, not Brug**
  (`analysis/oxide.py`, `doc/OXIDE_ANALYSIS_GUIDE.md`). The formula
  `C_eff = (R*Q)^(1/n) / R` used for R||Q elements is the Hsu & Mansfeld
  (2001) conversion via `tau = (R*Q)^(1/n)`, corresponding to a normal (3D,
  through-layer) distribution of time constants — the appropriate model for
  oxide layers. It was mislabeled as the Brug (1984) formula, which is a
  different expression (`C = Q^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)`) for a
  surface (2D) distribution and involves the series resistance; for n ~ 0.8
  the two differ by tens of percent. Documentation-only fix: computed values
  are unchanged; docstrings, comments, log messages, and the guide were
  corrected, and the 3D-model assumption is now documented explicitly.

---

## Version 0.16.13 (2026-07-02)

### Changed (audit 2026-07-02, cleanup findings 2.4-2.6)

- **Removed dead code from `fitting/diagnostics.py`** (finding 2.4):
  `check_parameter_diagnostics` and `log_fit_results` were never called
  anywhere in the package or tests; the latter was a stale duplicate of
  `_log_fit_result` in `cli/handlers/fitting.py` (which additionally handles
  bound status and CI suppression). No public API impact (neither was
  exported from `fitting/__init__.py`).
- **Fit quality log line shows the correct threshold per tier**
  (`cli/handlers/fitting.py`, finding 2.5): the line printed "(<10.0%)" for
  every tier including Excellent; it now shows the bound matching the tiers
  in `compute_fit_metrics` (<1.0% / <10.0% / <20.0% / >=20.0%), driven by the
  config constants.
- **Covariance warnings no longer overwrite each other**
  (`fitting/covariance.py`, finding 2.5): the "Negative variance" message is
  appended to an earlier ill-conditioned warning instead of replacing it.
  Also removed a dead `diag_cov = np.abs(diag_cov)` assignment (stderr
  already applies abs).
- **Unknown weighting string now logs a warning** (`fitting/diagnostics.py`,
  finding 2.5): `compute_weights` fell back to uniform weights silently;
  relevant for Python API callers (the CLI restricts choices in the parser).
- **Documented the BIC gap behavior in GMM early stopping** (`drt/peaks.py`,
  finding 2.5): when a middle n fails, improvement is compared across the gap
  between valid models with the same threshold (a deliberate simplification).
- **mypy: `drt/` package is clean again** (`drt/plotting.py`, finding 2.6):
  replaced `plt.cm.tab10` (attr-defined error, regression since the v0.13.17
  cleanup; the code moved from `core.py` to `plotting.py` in the v0.16.8
  split) with the equivalent `plt.get_cmap('tab10')`.

---

## Version 0.16.12 (2026-07-02)

### Fixed (audit 2026-07-02, finding 2.3)

- **A failed least_squares refinement in diffevo no longer masquerades as a
  successful one** (`fitting/diffevo.py`). When refinement raised an
  exception, `ls_result` was aliased to the DE result, so the equal costs made
  the selection report `refinement_improved=True`, attributed the DE
  message/status to the optimizer, and double-counted DE evaluations in
  `total_evaluations`/`n_function_evals` (via the aliased `nfev`). The failure
  path is now gated by an explicit `refinement_ran` flag:
  `refinement_improved=False`, `optimizer_message='DE only (refinement
  failed)'`, `optimizer_status=-1`, and evaluation counts include only the DE
  run. Fitted parameters, covariance, and the existing "Refinement failed"
  warning are unchanged. Also removed a redundant duplicate
  `compute_fit_metrics` evaluation in the selection branches.
  - Regression test `test_refinement_failure_reported_honestly` in
    `tests/test_diffevo.py`.

---

## Version 0.16.11 (2026-07-02)

### Fixed (audit 2026-07-02, finding 2.2)

- **`--normalize-rpol` no longer breaks peak resistance estimates**
  (`drt/core.py`, `drt/plotting.py`, `cli/handlers/drt.py`). Peak detection
  ran on the normalized gamma (gamma/R_pol), so `R_estimate` of both scipy
  and GMM peaks came out as dimensionless fractions (sum = 1) yet was printed
  as "R ~ ... Ohm", and the Voigt element analysis derived C = tau/R from
  these wrong values. Peak detection, reconstruction, and shape metrics now
  always use the unnormalized (physical) gamma; only the returned/plotted
  gamma is normalized. The GMM deconvolution panel likewise plots the
  unnormalized gamma, matching its component scaling in Ohm. Behavior without
  `--normalize-rpol` is unchanged.
  - Regression test `test_peak_r_estimate_invariant_to_normalize_rpol` in
    `tests/test_drt_recovery.py`.

---

## Version 0.16.10 (2026-07-02)

### Fixed (audit 2026-07-02, finding 2.1)

- **Scipy DRT peaks now report the characteristic frequency 1/(2·π·τ)**
  (`drt/core.py`). The `'frequency'` field of `scipy_peaks` was computed as
  `1/τ`, while GMM peaks (`drt/peaks.py`) and `fitting/auto_suggest.py` use the
  standard RC characteristic frequency `f = 1/(2πτ)` — the CLI printed both as
  "f = ... Hz", so the same relaxation process showed frequencies differing by
  a factor of 2π (~6.28×) depending on the peak method. Reported frequencies
  from the default scipy method drop by 2π; τ and R estimates are unchanged.
  Also fixed the same convention in `doc/DRT_INTUITION.md`.
  - Regression test `test_scipy_peak_frequency_convention` in
    `tests/test_drt_recovery.py`.

---

## Version 0.16.9 (2026-06-30)

### Changed (refactor, AUDIT_2026-06-23 section 4 / priority 4)

- **Split `fitting/circuit_elements.py` (647 lines) into a package.** It was the
  last `fitting/` file over the 500-line limit and the last big monolit named by
  the audit. Converted to `fitting/circuit_elements/`, mirroring the existing
  `cli/handlers/` precedent, grouped by element category:
  - `base.py` (134) — `CircuitElement` abstract base + operator overloading
  - `basic.py` (132) — `R`, `C`, `L` (lumped ideal elements)
  - `distributed.py` (152) — `Q` (CPE), `W`, `Wo` (Warburg diffusion)
  - `composite.py` (209) — `K` (Voigt R-τ), `G` (Gerischer)
  - `__init__.py` (37) — re-exports all 9 public symbols + `__all__`

  Pure structural refactor — **no behavior, signature, or numerics changes.**
  All consumers import via `eis_analysis.fitting.circuit_elements`, which the
  package `__init__.py` re-exports unchanged, so no consumer or test was touched.
  Verified: 193/193 tests pass, ruff clean, whole-project mypy count unchanged
  (35 — the one pre-existing `R.impedance` error just moved file), CLI smoke OK.

  This completes audit priority 4, whose two named targets — `cli/handlers.py`
  (v0.13.15) and `circuit_elements.py` (this release) — are now both split. Six
  files remain over the limit (down from 9), none flagged as priority.

---

## Version 0.16.8 (2026-06-30)

### Changed (refactor, AUDIT_2026-06-23 section 4 / priority 4)

- **Split `drt/core.py` (891 lines) into focused submodules.** It was the
  project's largest monolith and had grown since the audit (810 → 891), making
  it the worst file-size offender. The DRT pipeline is now spread across cohesive
  modules under `eis_analysis/drt/`, each well under the 500-line limit:
  - `core.py` (322) — orchestrator: `calculate_drt` + `_detect_peaks`
  - `results.py` (155) — the 6 result dataclasses (`DRTResult`,
    `DRTDiagnostics`, `RinfEstimate`, `LambdaSelection`, `NNLSSolution`,
    `DRTMatrices`)
  - `linear_system.py` (223) — frequency validation, matrix build, lambda
    selection, NNLS solve
  - `estimation.py` (137) — R_inf, R_pol, per-peak resistance, effective-bins
  - `plotting.py` (118) — `_create_visualization`

  Pure structural refactor — **no behavior, signature, or numerics changes.**
  All symbols (public dataclasses + `calculate_drt`, plus the private helpers
  imported directly by tests) remain importable from `eis_analysis.drt.core`
  via re-export, so `drt/__init__.py` and all tests are unchanged. Verified:
  193/193 tests pass, `drt/` stays mypy-clean, ruff clean, CLI smoke OK.

---

## Version 0.16.7 (2026-06-30)

### Fixed (mypy, AUDIT_2026-06-23 section 3 / priority 2)

Fixed the genuine latent-bug mypy errors outside `drt/` flagged by the audit,
as distinct from the numpy dtype-variance noise it classifies as cosmetic.
mypy: 55 → 35 errors. Tests 193/193, ruff clean.

- **`any` (builtin) used as a type annotation** (`fitting/multistart.py`,
  `fitting/diffevo.py`). The matplotlib-figure return slot and the
  `DiffEvoResult.de_result` field were annotated with the builtin `any`
  instead of `typing.Any`. Now `Any`.
- **Implicit Optional `max_iter: int = None`** in `robust_nnls`
  (`fitting/voigt_chain/solvers.py`) — PEP 484 violation, same pattern as
  audit finding 2.3. Now `Optional[int]` (the `None` default was already
  handled at runtime).
- **`auto_suggest` diagnostics dict masked type checking**
  (`fitting/auto_suggest.py`). `diagnostics` was inferred as
  `dict[str, object]`, so `.append()`, `len()` and iteration over its values
  failed type-checking across the whole function (10 errors). Annotated
  `Dict[str, Any]`; also fixed a `list` → `ndarray` reassignment of `peaks`.
- **`None` propagation in Voigt-chain build** (`fitting/voigt_chain/fitting.py`).
  `circuit` was inferred as `Optional[R]`, letting `None` flow into
  `circuit - L(...)` and `circuit.get_all_params()` — a real crash risk for a
  circuit with no series resistance and no K elements. Annotated
  `Optional[Circuit]` with an explicit `None` guard (raises `ValueError`). The
  return type was also corrected from `Series` to the truthful `Circuit`, since
  a single bare element can be returned.

---

## Version 0.16.6 (2026-06-30)

### Fixed (diffevo math audit #6)

- **Reported relative error no longer double-counts 1/|Z|**
  (`fitting/diagnostics.py`). `compute_fit_metrics` computed
  `Σ wᵢ·(|ΔZᵢ|/|Zᵢ|) / Σ wᵢ`; with modulus weighting (`wᵢ = 1/|Zᵢ|`) the 1/|Z|
  appeared twice, giving an effective 1/|Z|² emphasis on the residual. The
  metric is now weighting-consistent: `Σ wᵢ|ΔZᵢ| / Σ wᵢ|Zᵢ| · 100`, applying the
  weight once. For the default modulus weighting this equals the mean relative
  error `mean(|ΔZ|/|Z|)`. Reporting only (since audit #1 the DE/refinement
  selection uses the weighted SSR, not this metric); quality thresholds
  unchanged.
  - Regression tests in `tests/test_fit_metrics.py`
    (`test_modulus_equals_mean_relative_error`,
    `test_uniform_equals_aggregate_relative_error`,
    `test_perfect_fit_zero_error`).

This completes the diffevo math audit (`doc/AUDIT_diffevo.md`, findings #1-#6).

---

## Version 0.16.5 (2026-06-30)

### Fixed (diffevo math audit #5)

- **Confidence intervals now use the same degrees of freedom as the variance
  estimate** (`fitting/covariance.py`, `fitting/circuit.py`). The residual
  variance `s²` is computed with `dof = 2·n_freq − n_free_params` (complex data
  split into 2·n_freq real residuals), but `compute_confidence_interval` derived
  its t-quantile from `dof = n_freq − n_total_params` — a different value (factor
  of two on observations, plus counting fixed parameters). The t-multiplier
  therefore did not match the `s²` it scaled. The residual `dof` is now the
  single source of truth: `CovarianceResult.dof` (= `n_residuals − n_free_params`)
  is propagated to `FitResult` and used directly by the CI.
  - `compute_confidence_interval` signature changed: `n_data` → `dof`.
  - `FitResult._n_data` renamed to `FitResult._dof` (residual degrees of
    freedom). Internal field; the public CI properties are unchanged.
  - Regression test `tests/test_confidence_intervals.py::
    test_ci_dof_matches_covariance_residual_dof` (fails on the pre-fix code).

---

## Version 0.16.4 (2026-06-30)

### Fixed (diffevo math audit #4)

- **Rank-deficient covariance is now reported as infinite, not regularized**
  (`fitting/covariance.py`). For singular values below the rank threshold the
  code substituted `S_inv_sq = 1/threshold²` — an arbitrary, `S[0]`-dependent
  finite value that inflated the variance of unidentifiable directions by an
  undefined amount. When `JᵀJ` is singular the covariance does not exist, so the
  covariance and standard errors of the affected (free) parameters are now set
  to `inf` (the `scipy.optimize.curve_fit` convention); fixed parameters keep
  zero variance. Full-rank fits are unchanged (the regularization branch never
  triggered for them). Downstream already handles `inf` stderr (multistart falls
  back to log-uniform perturbation, confidence intervals become `±inf`).
  - Regression tests in `tests/test_covariance.py`:
    `test_rank_deficient_returns_inf`, `test_rank_deficient_with_fixed_params`,
    `test_full_rank_stderr_finite`.

---

## Version 0.16.3 (2026-06-30)

### Fixed (diffevo math audit #3)

- **`condition_number` now reports cond(JᵀJ), not cond(J)**
  (`fitting/covariance.py`). The covariance is `s²·(JᵀJ)⁻¹`, so the reliability
  of that inverse is governed by `cond(JᵀJ) = cond(J)²`. The code computed
  `cond(J) = S_max/S_min` but the docstrings (and the `is_well_conditioned <
  1e10` test) were written for `JᵀJ`, so the threshold effectively allowed
  `cond(JᵀJ)` up to ~1e20 — a numerically unreliable covariance could be
  flagged "well-conditioned". `condition_number` is now `(S_max/S_min)²` and the
  `1e10` threshold applies to `cond(JᵀJ)` as documented.
  - Side effect: borderline fits (cond(J) between 1e5 and 1e10) are now
    correctly flagged not well-conditioned; `multistart` falls back from
    covariance-based to stderr-based perturbation for those (fallback already
    existed).
  - Docstrings/docs aligned (`covariance.py`, `fitting/circuit.py`,
    `doc/PYTHON_API.md`, `doc/WEIGHTING_AND_STATISTICS.md`,
    `doc/MULTISTART_OPTIMIZATION.md`).
  - Regression tests in `tests/test_covariance.py`:
    `test_condition_number_is_cond_of_JtJ`,
    `test_well_conditioned_threshold_on_JtJ`,
    `test_well_conditioned_true_for_benign_jacobian`.

---

## Version 0.16.2 (2026-06-30)

### Fixed (diffevo math audit #1, #2)

- **DE-vs-refinement selection now uses the optimized objective**
  (`fitting/diffevo.py`). Differential evolution and the `least_squares`
  refinement both minimize the weighted sum of squared residuals
  (S = sum w^2 |dZ|^2), but the choice between their results — and the reported
  `improvement` — was made on the weighted *mean relative error* (a different,
  L1-style metric). Because the two metrics can disagree, a genuinely better
  refined fit could be discarded with a spurious "Refinement worsened fit".
  Selection and `improvement` now compare S directly (via the existing DE cost
  function). The `de_error` / `refined_error` percentages are unchanged and
  remain for display. New diagnostic fields `de_cost` / `refined_cost` expose
  the objective values.
- **Covariance is now evaluated at the returned parameters**
  (`fitting/diffevo.py`). When the DE result was kept, the covariance combined
  residuals at the DE point with the Jacobian from the (different) `least_squares`
  point, giving an invalid covariance/standard errors. The Jacobian is now
  recomputed at the chosen point (analytic when available). This also removes a
  latent `UnboundLocalError` (`cov_result`) that could crash the refinement
  failure fallback.
  - Regression tests in `tests/test_diffevo.py`:
    `test_diagnostics_expose_objective_costs`,
    `test_selection_picks_lower_objective`,
    `test_improvement_is_objective_based`,
    `test_refinement_failure_falls_back_with_valid_covariance`,
    `test_covariance_computed_at_returned_point`.

---

## Version 0.16.1 (2026-06-29)

### Changed (circuit suggestion: log edge-excluded peaks)

- **Edge-excluded DRT peaks are now logged with a reason** (`fitting/auto_suggest.py`).
  The auto circuit suggestion drops peaks sitting within 5% of the tau-grid
  edges (likely truncation artifacts), so the peak count can shrink silently
  (e.g. 2 detected -> 1 Voigt element). Each excluded peak now emits a warning
  with its tau, frequency, and reason, e.g.
  `Peak at tau = 1.59e+01 s (f = 1.00e-02 Hz) excluded: near right edge (low f)`.
  Logging only; the filtering itself is unchanged.
  - Regression test
    `tests/test_cli_integration.py::test_voigt_edge_peak_excluded_and_logged`.

### Fixed (CLI DRT peak listing)

- **"Found N peaks" now matches the peaks actually listed** (`cli/handlers/drt.py`).
  With GMM peak detection, `n_peaks` counts the merged GMM components, but the
  handler listed `scipy_peaks` (the raw `find_peaks` maxima kept for
  diagnostics), which GMM may merge via BIC. The header therefore reported e.g.
  2 peaks while 3 were printed. The handler now lists `result.peaks` (GMM
  components) when the method is `gmm`, and `scipy_peaks` otherwise, so the
  count and the listing always agree. Display-only; DRT/gamma and downstream
  fitting were unaffected.
  - Regression test `tests/test_cli_integration.py::test_drt_peak_count_matches_listing`
    asserts the listed peak count equals the reported count for both methods.

### Tests (DRT math audit F13)

- **Correctness tests for lambda selection and DRT reconstruction** —
  addresses `DRT_MATH_AUDIT_2026-06-27` finding F13. The lambda-selection tests
  were smoke-only (lambda positive/finite/in-range) and re-implemented the DRT
  matrix construction locally, so drift between test and production matrices
  would have been invisible.
  - `tests/test_hybrid_lambda.py` now imports the production
    `_build_drt_matrices` instead of a local copy (removed). Adds
    `compute_gcv_score` tests (finite/positive across the range; the
    GCV-selected lambda scores no worse than the range endpoints, i.e. the
    selector genuinely minimizes the score function).
  - New `tests/test_drt_recovery.py`: end-to-end `calculate_drt` recovery of
    known Voigt spectra — detected peaks land at the true time constants
    (within 0.15 decade) and the gamma integral recovers R_pol (within 3%),
    for both two-peak and single-peak circuits.

### Fixed (DRT math audit F10)

- **Unified R_pol integration with the DRT kernel** (`drt/core.py`,
  `drt/peaks.py`) — addresses finding F10. `R_pol_from_gamma` and the per-peak
  resistances used the trapezoid rule, while the DRT kernel integrates with the
  rectangle rule (`A_re -> d_ln_tau` as omega -> 0), so the reported R_pol did
  not exactly match the model's own DC limit. All R_pol computations now use the
  rectangle rule via the new `_rpol_from_gamma` helper:
  - total `R_pol_from_gamma` = `sum(gamma) * d_ln_tau`;
  - `_estimate_peak_resistance` partitions the tau axis into disjoint half-open
    valley segments (rectangle sum per segment), keeping `sum(R_i) = R_pol`
    exactly (simpler than the previous shared-node trapz);
  - GMM `R_pol` (for the `weight_i * R_pol` decomposition) likewise.
  - The difference vs trapz is negligible (~0.1%, half the endpoint values),
    but R_pol is now consistent with the reconstructed model. Peak-resistance
    tests updated to the rectangle reference; new `test_rpol_unified_integration`
    pins that the scipy, GMM and rectangle R_pol all agree.

Note: `fitting/auto_suggest.py` still uses trapz over threshold windows (a
separate path with its own known double-counting, out of scope here).

## Version 0.16.0 (2026-06-29)

### Added (DRT math audit F3 / F7)

- **DRT shape-quality detection** (`drt/core.py`) — addresses
  `DRT_MATH_AUDIT_2026-06-27` finding F3 (variant a+c: detect and warn).
  Lambda selection optimizes data fit, so for low-noise data auto-lambda
  drives lambda toward 0, producing a sparse/spiky DRT on which peak-shape
  analysis (scipy `find_peaks`, GMM) is meaningless. Auto-lambda is the CLI
  default (`python3 eis.py` without `--lambda`), so this affected normal use.
  - New participation-ratio metric `N_eff = (sum gamma)^2 / sum(gamma^2)`
    (`_effective_bins`), exposed as `DRTDiagnostics.n_effective_bins` and shown
    in the CLI. ~1 for a single spike, tens for a smooth distribution.
  - Warns when `N_eff < DRT_MIN_EFFECTIVE_BINS` (=7, calibrated: healthy DRT
    ~9-20, degenerate ~4-5.5) — "DRT is sparse/spiky ... consider a higher
    lambda".
  - **Lambda-edge detection (F7)**: the `corner_at_edge` flag already computed
    in `gcv.py` was being discarded by `_select_lambda`; it is now surfaced on
    `LambdaSelection` together with a new `lambda_at_edge` flag (set when the
    selected lambda, or the GCV guess, hits a search-range bound). Emits an
    "auto-lambda at search-range edge" warning.
  - Advisory only: gamma and the detected peaks are unchanged.
  - New `tests/test_drt_spikiness.py` (4 tests): metric correctness, degenerate
    auto-lambda warns, healthy lambda does not, lambda-edge detection.

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
