# Changelog - EIS Analysis Toolkit

Complete change history for all project versions.

---

## Version 0.28.1 (2026-09-01)

### Changed

- **Every block of CLI output now sits under the section it belongs to.** Two
  blocks did not. The per-point suspicious-point table was printed bare, right
  after the Z-HIT block, so it read as part of Z-HIT validation even though it
  draws on both validations; it now opens its own `Per-point residual check`
  section. And the hybrid lambda search (GCV + L-curve) narrated its progress
  from inside the solver, which runs before any of the DRT section headers
  exist - so it appeared as a loose block between the Nyquist figure and
  `R_inf estimation`. It is now reported inside `DRT Analysis`, on one line
  carrying both stages:

  ```
  Lambda: Hybrid GCV + L-curve
    lambda = 2.20e-04  (GCV 3.79e-04 -> L-curve corner 2.20e-04, ratio 0.58)
  ```

  The search itself is still visible with `-v`. Nothing about the selected
  lambda changed, only where and how it is reported.

  GMM peak detection (`--peak-method gmm`) had the same problem for the same
  reason and got the same treatment: its BIC model search no longer narrates
  itself from inside the solver, and `Peak detection in DRT spectrum` now
  reports which model BIC settled on, warns when that count sits at the upper
  bound of the searched range or when a candidate model failed to fit, and says
  so when early stopping picked a simpler model than the raw BIC minimum. The per-peak listing gained the two numbers only
  GMM provides - component width in decades, and mixture weight.

- **`--ri-fit` no longer reports R_inf twice.** The flag runs the estimation and
  prints it, then hands the value to the DRT stage, which printed the whole
  `R_inf estimation` section again - relabelled `Method: Preset value`, with the
  same number and the same median comparison. The DRT stage now reports an
  estimate only when it made one itself; a value handed in from outside is
  stated where it is used, with the comparison that places it against the data:

  ```
  Using R_inf = 1.216 Ohm (preset; HF median = 1.265 Ohm, -3.9%)
  ```

- **`generate_synthetic_data()` no longer prints from inside the library.** It
  announced the circuit it was building at INFO level, so any library caller got
  four lines of console output along with its data. The demo's CLI path prints
  the same block from `SYNTHETIC_DATA_PARAMS`, which is where those parameters
  are chosen. Console output is unchanged.

- **Z-HIT validation no longer prints from inside the library.** Its section
  header and result summary came from `zhit_validation()` itself, while the
  Kramers-Kronig equivalents came from the CLI - so the two validations drew
  their section boundaries in different places, and calling `zhit_validation()`
  as a library function wrote to the console. Both now print from the CLI
  handler, and the computation returns a result and nothing else. Console
  output is unchanged. `ZHITResult` gained a `success` property, matching
  `KKResult`, so a failed reconstruction is told apart from a valid one without
  inspecting array lengths.

- **`Lambda: GCV (automatic)` no longer hides a hybrid run.** The method was
  reported as `gcv` whenever the L-curve corner agreed with the GCV guess, so
  the common case looked like plain GCV. A successful hybrid search now always
  reports as such, and `GCV (L-curve correction failed)` means what it says:
  the fallback taken after the hybrid search raised. The two disagreement cases
  (L-curve correcting upward, or the geometric-mean compromise) are reported as
  warnings rather than being buried in the solver's own log.

---

## Version 0.28.0 (2026-09-01)

### Added

- **Suspicious points are now listed individually.** Both validations already
  computed a residual for every measured point, but the CLI only reported the
  averages - so a bad spectrum was visible while the bad *point* was not, and
  finding it meant reading the residual plot by eye. After KK and Z-HIT run,
  points whose relative deviation `|Z - Z_fit| / |Z|` exceeds `--max-residual`
  (default 5 %) in either method are listed with their frequency, both
  residuals and which method flagged them, and marked in the residual plot of
  the method that flagged them.

  The two methods contribute different things: Lin-KK fits both impedance
  components and is sensitive to errors in phase, Z-HIT reconstructs the
  magnitude from the phase and is sensitive to errors in magnitude and to
  drift. An injected spike is caught by both; taking the worse of the two
  catches what either one sees.

  Nothing is removed or down-weighted. This is deliberate: in EIS most
  low-frequency deviations are sample drift, a real property of a
  non-stationary measurement rather than a defect, and discarding those points
  would hide a bad experiment instead of reporting it. The output says so
  when flagged points reach the low-frequency end.

  Two guards keep the list honest. A point must also exceed four times the
  median residual *of the method that flagged it*: a residual is the sum of
  the data error and the method's own reconstruction error, and Z-HIT's is
  much the larger, since it rebuilds the magnitude from a numerically
  differentiated phase. Measured over 25 noise realizations of the synthetic
  spectrum, the factor cuts false alarms on defect-free data from 8.2 % of
  points to 1.9 % while still catching 23 of 25 injected 8 % spikes; on
  measured spectra it changes nothing, because there the absolute threshold
  stays binding. Second, a method whose *median* residual already exceeds the
  threshold is skipped entirely: over half its points would be flagged, so the
  spectrum fails as a whole, which the `Data quality:` line already says, and
  enumerating ninety points would only dress a global failure up as a list of
  local ones.

### Changed

- `run_kk_validation` and `run_zhit_validation` return their full result
  object (`KKResult` / `ZHITResult`) instead of just the figure, which is what
  makes the per-point residuals reachable from outside. The figure is still
  there as `result.figure`.

---

## Version 0.27.1 (2026-08-31)

### Changed

- **The Gamry parser now reads its column positions from the ZCURVE header row**
  instead of assuming `Freq`, `Zreal` and `Zimag` sit in columns 3 to 5. Gamry
  is consistent for ZCURVE, so nothing changes for the files we see; the point
  is that a file with a different layout now fails loudly rather than shifting
  every value one column across and yielding plausible nonsense. A header that
  names none of the three falls back to the old positional order and says so.

  The row-length check follows from the same place. It was a fixed `>= 5`;
  it is now derived from the rightmost column actually read, so a layout that
  puts `Zimag` further right no longer accepts rows too short to hold it.

### Fixed

- **An `EXPERIMENTABORTED` marker ahead of `ZCURVE` is now named for what it
  is.** A run aborted while the cell was still settling writes the ZCURVE
  table header with no rows beneath it. The parser reported that as
  `No valid EIS data found ... Check if file contains ZCURVE section`, which
  points at the one thing that is not wrong - the section is there. It now
  raises `Experiment ... was aborted before the impedance sweep started`.

- **Dead half of a guard in `parse_ocv_curve()`.** The section-end test read
  `not line[0].isdigit() and not line.startswith('\t')`, but the line had
  already been stripped, so the second clause could never be true. Behaviour is
  unchanged - the loop is bounded by the header's point count either way.

---

## Version 0.27.0 (2026-08-31)

### Added

- **Truncated sweeps are now flagged when the file is loaded.** The DTA header
  states the requested sweep - `FREQINIT`, `FREQFINAL`, `PTSPERDEC` - so the
  intended point count follows from it. `load_data()` compares that against
  what the file actually contains and warns when the run came up short:

  ```
  Sweep may be truncated: 92 points, header implies 96
    Lowest measured frequency 7.94e-04 Hz, header FREQFINAL 3.00e-04 Hz
  ```

  A sweep that stops above `FREQFINAL` loses exactly the low-frequency points
  that carry the slowest processes, which is where DRT resolution and any
  diffusion branch of a circuit fit are decided. Nothing downstream can detect
  the loss, because the remaining data is perfectly valid - it just answers a
  narrower question than the one the measurement was set up to ask. Three of
  the five reference exports in the working tree turn out to be short.

  The three header fields were already parsed; only the comparison is new.

- **`expected_points(metadata)`** in `io/data_loading.py`, exported from
  `eis_analysis.io`. Returns the point count the header implies, or `None`
  when the sweep parameters are missing or inconsistent (a non-EIS file, a
  truncated header, an inverted frequency range).

  Only a shortfall is reported. Gamry stops at the first point at or below
  `FREQFINAL`, so a complete sweep may overshoot the computed count by one,
  and treating that as an anomaly would warn on healthy files.

### Fixed

- **Diacritics in DTA metadata are no longer silently deleted.**
  `parse_dta_metadata()` and `parse_ocv_curve()` opened the file as UTF-8 with
  `errors='ignore'`, so on a file written by a Czech Windows machine - cp1250,
  which is what these instruments produce - every accented byte was dropped
  without an error or a warning:

  ```
  TITLE "Mereni vzorku"          ->  'Men vzorku'
  NOTES "elektroda c. 2"         ->  'elektroda . 2'
  ```

  Sample identification and operator notes were the two fields worst affected,
  and nothing downstream could tell that characters were missing.

  All three DTA parsers now share `_read_dta_lines()`, which decodes as UTF-8,
  then cp1250, then latin-1 - the last never fails, because it maps every
  byte - and folds the result to its closest ASCII equivalent via NFKD
  normalization. `read_gamry_native()` was already reading ISO-8859-1 and so
  never lost data; it moves to the shared reader for consistency, and its own
  UTF-8 fallback goes away because ISO-8859-1 decodes any byte sequence and
  could never raise to reach it.

  Folding is deliberate rather than a limitation. It makes the result
  identical whichever of the encodings was correct, so a wrong guess degrades
  to unaccented text instead of mojibake, and plain ASCII survives the PDF
  export, whose font does not cover the full Latin-2 range. Numeric fields
  cannot contain a non-ASCII character, so measurements are untouched.

  ISO-8859-2 - the Unix Czech encoding, not what Windows writes - is not in
  the list. cp1250 decodes it without raising, so most of its text survives,
  but the few characters where the two code pages disagree fold away.

---

## Version 0.26.0 (2026-08-30)

### Added

- **`--circuit` can be repeated to compare candidate models by AIC/BIC.**
  Every candidate is fitted on the same data with the same weighting and the
  results are ranked in one table with the free parameter count, fit error,
  dAIC, dBIC and the condition number.

  Adding an element to a circuit almost always lowers the residual, so a
  comparison of fit errors systematically favours the most complex candidate.
  AIC and BIC charge for each free parameter and answer the question the fit
  error cannot: does the data support the extra element? Only differences
  carry meaning - below 2 the models are indistinguishable, above 10 the
  difference is decisive - and only within one run, since a change of
  weighting or frequency range invalidates the comparison.

  BIC selects the reported winner, which is also what `--analyze-oxide`
  receives; the selected circuit is named in the log. A candidate that fails
  to fit - including one whose expression does not parse - is reported in the
  table and skipped rather than ending the run. Figures are saved per
  candidate (`_fit_1`, `_fit_2`, ...) and the table's `-c` column gives the
  command-line position, so a reordered ranking still maps to its files.

  The condition number is in the table on purpose: an overparametrized
  circuit can win on AIC while its parameters are not jointly identifiable.
  A winner flagged `!` (cond > 1e10) is a warning, not a recommendation.

  **A single `--circuit` behaves exactly as before** - one fit, no table, and
  the figure keeps its original filename.

- **`compute_information_criteria`** in `fitting/diagnostics.py`, exported
  from `eis_analysis.fitting`. Returns the weighted RSS, AIC and BIC for a
  fitted model. It rejects `n_free_params <= 0`: a zero would charge no
  complexity penalty at all and hand that model every comparison it enters,
  which is precisely what an unpopulated field would have caused silently.

- **`doc/MODEL_SELECTION_AIC_BIC.md`** - what the criteria measure, the
  arithmetic worked through on a real case, and what they cannot tell you.

- **`FitResult.n_free_params`** - the number of freely optimized parameters
  (fixed parameters excluded), which is the `k` the criteria charge for.
  Previously only the private `_dof` carried this, and only indirectly.

### Notes

- The Voigt chain (`--voigt-chain`) deliberately does not enter the ranking.
  Its `K(R, tau)` elements sit on a fixed tau grid that is not optimized, so
  a naive parameter count overstates its complexity by M (about 35 BIC points
  at M = 7, enough to invert the ranking), and the fair count would still
  only be a lower bound because both the grid and `prune_threshold` are
  data-informed. `--voigt-chain` keeps its own separate branch unchanged.

---

## Version 0.25.5 (2026-08-30)

### Added

- **Regression benchmark on the ZScope reference circuits**
  (`tests/test_zscope_benchmark.py`, 20 tests). Parameter recovery is checked
  on the four circuits of the published ZScope v2.2.0 benchmark - Randles,
  Randles+Warburg, CPE Randles and two time constants - over their frequency
  grid (80 points, 1e-2 to 1e5 Hz).

  Three assertions per circuit: exact recovery of the generating parameters
  from noise-free data (< 0.1%), parameter recovery within a calibrated band
  at 2% and 5% proportional noise, and convergence to the noise floor. The
  ground truth was verified against the ZScope CSVs to max|dZ|/|Z| ~ 5e-7;
  the data is regenerated in the test rather than vendored, because the
  license of their CSVs is unclear (their repository ships only the CC BY 3.0
  license of its website template).

  The Two_TimeConstants comparison orders the RC branches by tau before
  comparing. A series chain of Voigt elements is invariant under permutation
  of its branches, and the fitter does return them swapped - it did so on the
  ZScope 2% data at an otherwise correct RMSE of 2.0%.

  See `doc/ZSCOPE_COMPARISON.md` for the full comparison against ZScope.

---

## Version 0.25.4 (2026-08-30)

### Fixed

- **A confidence interval no longer overflows on a parameter the data did not
  constrain** (`fitting/covariance.py`). The log-space CI is
  `(p/f, p*f)` with `f = exp(t * se / p)`, and `exp` overflows once the
  relative uncertainty passes ~709. A scale parameter driven onto its bound
  reaches that easily - the plain `eis` demo run does it, with `t*se/p = 4.7e5`
  for an R0 sitting on its 1e-4 lower bound - so every such fit printed
  `RuntimeWarning: overflow encountered in exp` to stderr.

  The reported interval was never wrong: `f = inf` gives `(p/inf, p*inf) =
  (0, inf)`, which is the exact limit. The limit is now taken directly instead
  of being produced by the overflow, so the values are unchanged and the
  warning is gone. The parameter that triggered it is one whose CI the CLI
  does not even print - it reports `[at lower bound - CI not meaningful]` -
  so the warning was noise about a discarded number.

  Regression test in `tests/test_confidence_intervals.py` turns RuntimeWarning
  into an error for the overflowing case.

### Documentation

- `README.md`: `--save` documented the wrong output filename
  (`results_nyquist.png`; the plot is written as `results_nyquist_bode`), and
  the CSV section listed `frequency`/`Z_real`/`Z_imag` as if those exact names
  were required - the loader accepts a family of aliases and falls back to
  positional columns.

---

## Version 0.25.3 (2026-08-29)

### Fixed

- **A fitted capacitance with no resistance beside it is no longer ignored**
  (`analysis/oxide.py`). The traversal only registered a `C` or `Q` that
  shared a `Parallel` with an `R` - a Voigt element. In `L - R0 - (Q|C)` the
  `C` has no resistance next to it, so `--analyze-oxide` reported
  `No Voigt (R||C), K, or R||Q elements found in circuit` and fell through to
  the high-frequency spectral estimate - which itself warned that its
  per-point values spread by a factor of 1.37 across the top decade - even
  though `C` was a fitted parameter with a 0.7 % confidence interval. It
  happened to land within 2 % of the reference value; with other data or
  another circuit the error could have been large and invisible.

  `_find_parallel_rc_elements()` is now `_find_capacitive_elements()` and
  collects every `C`, `Q`, `K` and `CC` in the tree, with or without a
  parallel resistance. A resistance, where present, is still recorded - it is
  what the Hsu-Mansfeld/Brug conversion, `tau = R*C` and the largest-R
  heuristic need.

### Changed

- **The dielectric element is now chosen by physics, not by element type or
  position in the expression.** Selection was conflating two questions -
  which element carries the dielectric response, and which value to take from
  it. They are now separate steps.

  Step 1: the dielectric element is the one whose admittance rises as
  `omega^n` with `n` near 1. `C` and `K` are `n = 1` by construction, a `CC`
  is `n = 1` in both of its limits, and a CPE qualifies only when its exponent
  is near-ideal (reusing `CPE_N_RELIABLE_MIN = 0.8`); a CPE at `n ~ 0.6` is
  transport or a distribution of resistivity, not a dielectric. Among the
  qualifying elements the order is `CC` (a plain `C` is its degenerate case
  `dC = 0`, so the general model wins), then `C`/`K` whose capacitance is a
  fitted parameter, then near-ideal `Q` whose capacitance still needs a model.
  Within a tier the largest-R heuristic is unchanged.

  This replaces the `Multiple C/Q elements in one parallel combination - using
  the last one` behaviour: position in the circuit expression was never a
  physical criterion. Elements that share one parallel resistance are now
  ranked by capacitance, with a warning that their individual values are not
  separately identifiable from the spectrum.

  A circuit whose only capacitive elements are low-`n` CPEs still uses them -
  they are fitted parameters, unlike the spectral fallback - but says plainly
  that the result has no dielectric meaning.

- **The high-frequency fallback says it is not from the fit.** It now leads
  with a warning of the same weight as a parameter sitting on its bound: the
  capacitance carries no confidence interval and any thickness or permittivity
  from it is an order-of-magnitude figure, however close it happens to land.
  A `Q` with no parallel resistance is reported as not convertible (both
  Hsu-Mansfeld and Brug need `R`) rather than silently skipped.

---

## Version 0.25.2 (2026-08-29)

### Fixed

- **`--analyze-oxide` no longer reports a Cole-Cole capacitance the data never
  measured** (`analysis/oxide.py`, `analysis/config.py`). For a dominant `CC`
  element the static capacitance `C_s = C_inf + dC` was reported
  unconditionally. `C_s` is the exact `omega -> 0` limit of `C*(omega)`, but
  only the *measured* one when the relaxation lies inside the frequency
  window. A fit that landed on `tau = 1e4 s` - the upper bound of
  `PARAMETER_BOUNDS['tau_CC']` - put `f_char = 1/(2*pi*tau) = 1.6e-5 Hz` three
  decades below `f_min`, so every measured point sat at `omega*tau >> 1`, where
  `C*(omega) -> C_inf`. `dC` was consequently unidentified (95% CI spanning a
  factor of 4) yet dominated the reported `C_s = 1.611e-05 F`, giving
  `eps_r = 3491` where the value the data support, `C_inf = 8.117e-08 F`,
  gives 17.6 - matching the reference model's 17.4.

  The reported capacitance now follows where `f_char` sits relative to
  `[f_min, f_max]`:
  - inside the window -> `C_s` (the relaxation is traced);
  - below `f_min` -> `C_inf`, with a warning that `dC` is an extrapolation to
    DC through a region with no data;
  - above `f_max` -> `C_s`, with a warning that only the *sum* is identified,
    not the `C_inf`/`dC` split;
  - within one decade of either edge -> `C_s`, flagged as only marginally
    determined (`CC_WINDOW_EDGE_MARGIN_DECADES`).

  The dominant-CC selection (`max` over the static capacitances) uses the same
  window-aware value, so the choice and the reported number cannot disagree.

- **A Cole-Cole `tau` sitting on a fitting bound is now said out loud.** The
  window test above happens to cover the reported case, but a parameter on its
  bound always means the data did not determine it, so the check is explicit
  rather than inferred from `f_char`. It reuses `classify_bound_status()` and
  `PARAMETER_BOUNDS['tau_CC']` from `fitting/bounds.py` - the project-wide
  definition of "at a bound" - and is skipped for a `tau` the user fixed
  (`CC(..., tau="1e4")`), which is a choice rather than a failed fit.

  Regression tests in `tests/test_oxide.py` cover the reported case end to end
  (`eps_r` ~ 17.6 instead of 3491), the lower-bound case where the window rule
  and the bound test point in opposite directions (`C_s` stands, the bound
  still warns), the fixed-`tau` case, the near-edge case, and the in-window
  path, which is unchanged.

---

## Version 0.25.1 (2026-08-27)

### Fixed

- **The Gerischer element `G` is finally reachable from `--circuit`**
  (`cli/utils.py`). It has been fully implemented since v0.13.2 - impedance,
  analytic Jacobian, `tests/test_G_element.py`, a row in the README element
  table and its own section in `doc/CIRCUIT_PARSER.md` - but was missing from
  the parser's `safe_namespace`, so `--circuit "G(100,1e-3)"` failed with
  `name 'G' is not defined`. A user following the README hit a wall on a
  feature the codebase already supported, for eight months.

- **`from eis_analysis import G` works.** The same v0.13.2 commit also
  omitted `G` from the top-level `__init__.py` import and `__all__`, so the
  element was only reachable through `eis_analysis.fitting`.

  The root cause of both is a duplicated registry: the element list exists in
  `circuit_elements/__init__.py`'s `__all__`, again in `fitting/__init__.py`,
  again in the top-level `__init__.py`, and again as the hand-written
  `safe_namespace` dict. Adding an element means updating all four, and `G`
  was missed in two of them. `tests/test_circuit_expression.py` is the guard
  in the meantime: it parses every element the README documents and checks
  each is exported top-level, so the next omission fails a test rather than
  reaching a user.

- `doc/CIRCUIT_PARSER.md`'s `safe_namespace` snippet listed `G` (which did
  not exist there) and omitted `CC` (which did). It now matches the code.

---

## Version 0.25.0 (2026-08-27)

### Fixed

- **Element parameter attributes no longer go stale after a fit**
  (`fitting/circuit_elements/`). Every element copied its parameters into
  named attributes in `__init__` (`self.R = self.params[0]`, ...).
  `update_params()` - which the fitter calls with the optimized values -
  rewrote `params` but left those copies untouched, so after a fit `K.R`,
  `G.sigma`, `CC.C_inf` and the properties derived from them
  (`K.capacitance`, `G.characteristic_freq`, `CC.C_static`) still returned
  the initial guess. Nothing warned; the numbers simply came from before
  the fit.

  It was found the hard way in 0.24.0: `--analyze-oxide` reported a film
  thickness of 17.7 nm against a true 25.0 because the Cole-Cole branch read
  `node.C_static`. Every other branch in `analysis/oxide.py` reads
  `node.params[i]`, which is precisely why they were unaffected - the
  convention existed because the attributes could not be trusted.

  The names are now read-only views of `params`, built by the new
  `param_property()` helper in `circuit_elements/base.py`, so they cannot
  fall out of step. `params` is the single source of truth. This also
  removes 16 lines of copying from the nine `__init__` methods, and a new
  element gets the behaviour by declaring `R = param_property(0)` rather
  than by remembering to keep a copy in sync.

### Changed

- **BREAKING**: assigning to a parameter attribute (`element.R = 200`) now
  raises `AttributeError`. It previously succeeded but changed only the
  attribute, leaving `params` - and therefore the computed impedance -
  untouched, so `repr(element)` and `element.impedance()` disagreed. Set
  values through the constructor or `update_params()`.

---

## Version 0.24.0 (2026-08-27)

### Added

- **`--analyze-oxide` now recognises the Cole-Cole element** (`analysis/oxide.py`).
  Fitting an oxide film with `CC` and then asking for its thickness used to
  report `No Voigt (R||C), K, or R||Q elements found` and fall back to the
  high-frequency estimate - the one element written for dielectric films was
  the one the oxide analysis could not see.

  The capacitance used is the static one, `C_s = C_inf + dC`. Unlike the Q
  case this is exact rather than effective: `C_s` is the `omega -> 0` limit of
  `C*(omega)` for any `alpha`, so no Hsu-Mansfeld / Brug modelling choice
  arises and no exponent-threshold warning is needed. It is also the
  capacitance that pairs with a *static* permittivity such as `eps_r = 22` for
  ZrO2; `C_inf` would give the high-frequency value.

- **A `CC` in the circuit takes precedence over the largest-R heuristic.**
  That heuristic exists to tell a compact oxide apart from a charge-transfer
  process when both appear as R||C arcs; a Cole-Cole element models the
  dielectric relaxation explicitly, so there is nothing to disambiguate. With
  more than one `CC` the largest static capacitance is used and a warning is
  logged - assigning several relaxations to layers is the operator's call.

  A `CC` is a blocking dielectric with no DC path, so it reports
  `element_R = None` (the log says `n/a`) and its time constant is the
  relaxation time, not `R*C`. `element_type` gains the value `'CC'`.

### Notes

- `_find_parallel_rc_elements()` reads `node.params`, never an element's
  named attributes. `update_params()` rewrites `params` but leaves attributes
  such as `K.R`, `G.sigma` and `CC.C_inf` - and the properties derived from
  them, `K.capacitance`, `G.characteristic_freq`, `CC.C_static` - at their
  construction-time values, so those still hold the initial guess after a fit.
  `tests/test_oxide.py::test_cc_uses_fitted_values_not_the_initial_guess`
  pins this.

---

## Version 0.23.0 (2026-08-27)

### Removed

- **The element scaling operator `*` and the CPE exponent operator `**`.**
  `2 * R(100)` was another way to write `R(200)`, and `Q(1e-4) ** 0.9`
  another way to write `Q(1e-4, 0.9)`. Neither enabled anything; both were
  a second spelling of a value you can type directly.

  Their cost was not zero. `_scale()` was an `@abstractmethod`, so all nine
  elements had to implement it just to keep the operator working - roughly
  25 lines of ceremony, plus one more obligation on every element added
  later. Nothing inside the package ever constructed a circuit by scaling.

  They also carried a bug. Scaling discarded the fixed-parameter flags:
  `2 * R("100")` returned a *free* 200 Ohm resistor, so a parameter pinned
  on purpose - an `R_inf` from a separate estimate, an oxide capacitance
  from a known thickness - would silently be refitted, with the degrees of
  freedom (and therefore the confidence intervals) computed as if it were
  free. `Q("1e-4", 0.8) ** 0.9` lost the coefficient's fixed status the same
  way, though keeping it is the entire point of that operator. `G` was worse
  still: `G._scale` wrapped values in literal quote *characters*, which
  reached `float()` in `CircuitElement.__init__`, so `2 * G("100", 1e-3)`
  raised `ValueError` outright.

  Deleting the operators removes the defect rather than fixing it in nine
  places. Building circuits with `-` (series) and `|` (parallel) is
  unaffected.

### Migration

```python
2 * R(100)          ->  R(200)
0.5 * C(1e-6)       ->  C(5e-7)
Q(1e-4) ** 0.9      ->  Q(1e-4, 0.9)
```

---

## Version 0.22.0 (2026-08-27)

### Added

- **Cole-Cole element `CC(C_inf, dC, tau, alpha)`** for dielectric relaxation
  (`fitting/circuit_elements/distributed.py`). Its impedance is
  `Z = 1/(j*omega*C*)` with `C* = C_inf + dC/(1 + (j*omega*tau)^(1-alpha))`,
  so it models a *distribution* of relaxation times as a depressed arc in the
  complex capacitance (equivalently permittivity) plane. This is a different
  quantity from what `Q` describes: `Q` depresses the arc in the impedance
  plane, which is the right picture for a rough electrode but not for the
  dielectric dispersion of an oxide or polymer film. `alpha = 0` recovers
  Debye relaxation. The element carries no geometry - its parameters are
  capacitances, so the conversion to a relative permittivity stays in the
  analysis layer where the film thickness and electrode area are known.
  Available from `--circuit`, from `eis_analysis` and `eis_analysis.fitting`,
  and it has an analytic Jacobian for all four parameters.

  Parametrised by the relaxation strength `dC = C_s - C_inf` rather than by
  the static capacitance `C_s`, because `dC > 0` from the bounds then enforces
  `C_s > C_inf` on its own; box bounds cannot express a coupling between two
  parameters.

  The element is exactly equivalent to
  `C(C_inf) | (C(dC) - Q(dC/tau**(1-alpha), alpha))`, which the test suite
  uses as an independent oracle. That composite is not a substitute for
  fitting it directly: its CPE coefficient couples three parameters
  non-linearly, so a fit would report `Q` and `n` with confidence intervals
  on the wrong quantities, and `dC` would appear twice as two independent
  free parameters.

- **Bounds for the four Cole-Cole parameters** (`fitting/bounds.py`).
  `alpha` is bounded `(0.0, 0.9)`; the lower bound is exactly `0.0` so that
  `log_scale_ci_mask()` classifies it as a linear parameter, like the CPE
  exponent `n`. Without an entry in `PARAMETER_BOUNDS` a new label silently
  falls back to `(1e-15, 1e15)`, and DE would then have searched `log10(alpha)`
  across 30 decades for an exponent that belongs in `[0, 1)`.

---

## Version 0.21.6 (2026-08-27)

### Changed

- **The DE `Improvement` line now reports the optimized objective**
  (`cli/handlers/fitting.py`). It recomputed its own value from the two
  weighted relative errors, while `DiffEvoResult.improvement` - computed
  from the weighted SSR, which is what the DE/refinement selection actually
  uses - sat unused. The two can move in opposite directions: on a
  high-impedance oxide sample the SSR fell from `1.247761e9` to `1.247604e9`
  (the refinement genuinely improved the fit) while the relative error rose
  from 2.827% to 2.830%, so the log reported `Improvement: -0.1%` for a fit
  that had got better. The line now reads `Improvement (SSR)` and uses the
  prepared field. **Only this displayed number changes** - the selection,
  the reported errors and the fit results are untouched.

### Refactored

- **`estimate_R_linear` no longer carries its own copy of `compute_weights`**
  (`fitting/voigt_chain/fitting.py`). The same four weighting branches and
  the same mean-1 normalization existed twice; the only divergent branch
  (a silent fallback to modulus) was unreachable, since the function already
  validates `weighting` against those four values and raises otherwise.
  Weights are bit-for-bit identical on all four weighting types.

- **`fit_circuit_diffevo` computes the fit metrics twice instead of three
  times** (`fitting/diffevo.py`). The chosen point is always either the DE
  or the refined one, both already evaluated, so the third pass over the
  spectrum recomputed a known result.

## Version 0.21.5 (2026-08-26)

### Changed

- **Differential Evolution now searches `log10` of the scale parameters**
  (`fitting/diffevo.py`). DE sampled its population uniformly over bounds
  spanning 6-14 decades (`R: 1e-4..1e10`, `Q: 1e-12..1e-1`,
  `C: 1e-15..1e-1`), so nearly every member was drawn from the top decade:
  the parallel branches were shorted, every member predicted the series
  resistance alone, and the population energies were so nearly equal that
  scipy's convergence test (`std <= tol * |mean|`) could fire after the
  first generation. On a high-impedance oxide sample DE reported
  "converged" after **1 iteration at 99.970% error** and the whole fit came
  from the local refinement - which is what left parameters unidentifiable
  (`R_ct = 2.99e8 +/- 2.58e8`, rank-deficient Jacobians with `+/- inf` CIs).
  Parameters whose bounds span more than `LOG_SCALE_BOUND_RATIO` (6 decades)
  are now searched logarithmically; the CPE exponent `n` (0.3-1.0) stays
  linear. Measured on that sample: 99.970% -> **19.109%** after the DE stage
  alone (and 66.1% -> 19.25% on the unscaled spectrum). The refinement and
  the analytic Jacobian are unchanged - they still work in linear space, and
  `DiffEvoResult.de_result.x` is still reported in physical units.
  **DE fit results shift** for circuits with decade-spanning parameters.

### Added

- **A warning when the global search contributed nothing**
  (`fitting/diffevo.py`, `config.py`). If DE ends above
  `DE_STALLED_ERROR_PCT` (50%) and the local refinement is at least
  `DE_STALLED_IMPROVEMENT_FACTOR` (10x) better, the fit rests on that single
  local run, not on the global search - now said out loud, with the DE
  iteration count and both errors.

- **The DE diagnostics name the log-searched parameters**
  (`cli/handlers/fitting.py`): `Search space: log10 for R1, Q0, C0`.

## Version 0.21.4 (2026-08-26)

### Fixed

- **A circuit with every parameter fixed crashed with an opaque scipy
  message** (`fitting/bounds.py`, `circuit.py`, `diffevo.py`). The
  free-parameter vector is then empty: DE reported "bounds should be a
  sequence containing finite real valued (min, max) pairs", the local fit
  "index -1 is out of bounds for axis 0 with size 0". All optimizers now
  share `validate_fixed_params()` and say what is wrong ("All 4 circuit
  parameters are fixed - there is nothing to fit").

- **The local fit silently clipped fixed values into the free-parameter
  bounds** (`fitting/circuit.py`). `R("1e-9")` was fitted as R = 1e-4 and
  `R("1e12")` as 1e10 with `--optimizer single` / `multistart`, while DE
  kept the value as given. Fixed parameters are never optimized, so the
  bounds do not apply to them; they are no longer clipped, and all three
  optimizers now return the value that was fixed.

- **A fixed value outside the physically reasonable range passed silently**
  (`fitting/bounds.py`). Fixing bypasses `PARAMETER_BOUNDS`, so `R("-5")`
  (negative resistance) or a CPE exponent `n = 1.5` entered the fit without
  comment. Such values are still honored - they are the caller's choice -
  but now produce a warning naming the parameter and the range.

- **Multi-start returned a best fit whose circuit held the wrong parameters**
  (`fitting/multistart.py`). Every restart fitted the same `Circuit` object and
  `fit_equivalent_circuit()` writes its parameters into it at the end, so the
  circuit - and therefore `best_result.circuit` - ended up holding the
  parameters of whichever restart ran last. Consumers that read parameters
  from the circuit tree rather than from `params_opt` (oxide capacitance, and
  thus `--thickness` permittivity and thickness estimates) silently used a
  non-optimal fit; on a two-Voigt test case eps_r came out 64.8 instead of
  58.8. In parallel mode (`parallel=True`, library API only - the CLI runs
  sequentially) the same sharing was also a data race: several threads wrote
  the circuit's parameters while `impedance()` read them. Each fit now runs on
  its own copy of the circuit, so every result's circuit matches its own
  parameters, and the circuit passed in by the caller is synced to the best
  fit at the end, as the single-fit and DE paths already did. Fits with
  `--optimizer single` / `de` were never affected.

## Version 0.21.3 (2026-08-18)

### Fixed

- **The Brug (2D) capacitance was reported even when the fit had not
  identified the series resistance** (`analysis/oxide.py`). A CPE with
  n < 1 mimics a series resistance at high frequency, so R_s is often
  unidentifiable and the optimizer drives it to its 0.1 mOhm lower bound.
  Brug scales as `C ~ R_s^((1-n)/n)`, so a floored R_s produced an
  arbitrarily small capacitance: on a ZrO2 sample it differed from
  Hsu-Mansfeld by 233x and yielded eps_r = 1.55 (near vacuum) with no
  warning. Below `BRUG_RS_MIN_OHM` (10 mOhm) the Brug estimate is now
  suppressed and a warning names the cause and the check (compare the
  fitted R_s against Re(Z) at the highest measured frequency).
  Hsu-Mansfeld, which never uses R_s, is unaffected.

- **A large Hsu-Mansfeld/Brug spread was presented as a comparison**
  (`analysis/oxide.py`). The ratio is exactly `(1 + R_ct/R_s)^((1-n)/n)`;
  past `BRUG_HM_DIVERGENCE_MAX` (10x) the pair no longer brackets a single
  C_eff and the spread only reflects how well R_s is known. The two values
  are still printed, now with a warning that any derived thickness or
  permittivity is an order-of-magnitude estimate.

## Version 0.21.2 (2026-07-29)

Documentation only; no behavior change.

### Fixed

- **`_extract_capacitance()` docstring understated its return value**
  (`analysis/oxide.py`). Both return paths have included `C_eff_brug` and
  `C_specific_brug` since 0.20.0, but the documented key list never
  mentioned them. Since 0.21.0 both callers read those keys, so a reader
  following the docstring would conclude the Brug values came from
  somewhere else. The list now names all eight keys and states when the
  optional ones are None.

## Version 0.21.1 (2026-07-29)

### Fixed

- **An explicit `--area 1.0` was silently overridden by DTA metadata**
  (`cli/handlers/oxide.py`). The handler used `args.area == 1.0` as a
  proxy for "the flag was not given", which argparse cannot distinguish
  from an explicitly passed 1.0. On a DTA file carrying its own `AREA`
  header, 1.0 was therefore the one value a user could not force.

  Area scales the result directly - `C_specific = C_eff / area` and
  `d = epsilon_0 * epsilon_r * area / C_eff` - so a file declaring
  `AREA = 0.5` silently halved a thickness the user had asked to be
  computed for 1 cm². The log did say `Using area from DTA metadata`,
  but that reads as routine information rather than as "your value was
  discarded".

  `--area` now defaults to `None`, the same sentinel-free pattern applied
  to `--epsilon-r` in 0.21.0: an explicit value always wins, metadata
  fills in only when the flag is absent, and `DEFAULT_AREA_CM2` (1.0, new
  documented constant in `analysis/config.py`) applies when there is
  neither. `--help` now states the real default.

  - Regression tests in `tests/test_cli_integration.py`: explicit
    `--area 1.0` beats metadata 0.5, and an omitted flag still picks the
    metadata value. Both were verified to fail against the old handler.

## Version 0.21.0 (2026-07-29)

### Added

- **Inverse oxide analysis from the CLI: `--thickness`.** With a thickness
  known independently (SEM/TEM cross-section, ellipsometry), the analysis
  now runs backwards - the thickness is the input and the relative
  permittivity the estimated quantity:

  ```bash
  eis data.DTA --circuit "R(100) - (R(5000) | C(1e-6))" \
      --analyze-oxide --thickness 25
  ```

  This is a cross-validation step: an estimated epsilon_r near the
  literature value for the expected oxide (ZrO2 ~ 20-25) indicates the
  chosen equivalent circuit is physically consistent, while a wild value
  points at the circuit or the element selection rather than the material.

  The computation itself is not new - `estimate_permittivity()` has
  implemented `epsilon_r = d * C_specific / epsilon_0` since 0.16.x and
  shares its element-selection core with `analyze_oxide_layer()`. It was
  simply unreachable from the CLI, available only through the Python API
  or the ad-hoc `oxide_permittivity.py` script.

- **Brug (2D) comparison permittivity.** The inverse direction now reports
  a second epsilon_r derived from the Brug effective capacitance, matching
  what the thickness direction has done since 0.20.0. Available for a
  dominant Q element when the circuit contains a series resistance;
  exposed as `OxideAnalysisResult.permittivity_brug` and logged as
  `ε_r (Brug)`. The spread between the two CPE-distribution models
  brackets the model uncertainty, and it is as relevant in one direction
  as the other. Printed with `.3g` rather than a fixed-point format:
  when n is far from 1 the two models diverge by orders of magnitude and
  `.1f` would render a small-but-nonzero value as a misleading `0.0`.

### Changed

- **BREAKING: `estimate_permittivity()` returns `OxideAnalysisResult`**
  instead of a bare `float` (`analysis/oxide.py`). Returning a scalar left
  nowhere to put the Brug comparison value, and made the function
  asymmetric with `analyze_oxide_layer()`, which has always returned the
  full dataclass. The result object carries the same capacitance and
  element fields as the forward direction, with the roles of input and
  output swapped: `thickness_nm` holds the value that was passed in, and
  `permittivity` holds the estimate. `thickness_brug_nm` stays `None`
  here, since in this direction the thickness is not derived.

  Migration: `eps_r = estimate_permittivity(...)` becomes
  `result = estimate_permittivity(...)`, then `result.permittivity`.

- **`--epsilon-r` now defaults to `None`**, resolved to
  `DEFAULT_EPSILON_R` (22.0, new documented constant in
  `analysis/config.py`) by the oxide handler. The previous `default=22.0`
  made an explicit `--epsilon-r 22` indistinguishable from the flag being
  absent - the same sentinel weakness `--area` still has. Distinguishing
  them is what allows the handler to warn, rather than silently ignore,
  when `--epsilon-r` is combined with `--thickness`, where the
  permittivity is the estimated quantity and any supplied value is
  meaningless. In the spirit of audit finding O3 (2026-07-02): no silent
  assumptions.

  User-visible default is unchanged; `--help` still reports 22 for ZrO2.

  - Tests added to `tests/test_oxide.py` (Brug permittivity present and
    scaling linearly with C, absent without a series R, input thickness
    carried through), `tests/test_cli_parser.py` (None defaults, and that
    conflicting flags parse rather than erroring) and
    `tests/test_cli_integration.py` (inverse mode reports permittivity and
    not thickness; conflicting `--epsilon-r` warns).

## Version 0.20.5 (2026-07-26)

Build configuration only; no code changes.

### Fixed

- **CI `typecheck` job aborting with a fatal config error instead of type
  checking.** `mypy` 2.x rejects `python_version = "3.9"` outright
  (`Python 3.9 is not supported (must be 3.10 or higher)`) and exits 2 without
  checking anything - so the job had not actually been type checking since
  `mypy` 2.0 reached CI, it had merely been failing early. The `continue-on-error`
  flag on that job hid the difference: a red mark either way.

  The `dev` extra now pins `mypy<2`. The project still declares
  `requires-python = ">=3.9"` and CI runs the test matrix on 3.9, so holding
  the checker back is the option consistent with that promise; raising
  `python_version` to `"3.10"` would have silenced the message while quietly
  dropping 3.9 from what gets verified. Verified with `mypy` 1.20.2: exit
  code 1 and the same 34 pre-existing errors, with no fatal config error.

  Those 34 remain non-blocking and are mostly numpy dtype variance the checker
  cannot see through. Two are worth a note, both in `fitting/circuit.py`
  (lines 336 and 441): an `Optional[List[float]]` used without being narrowed.
  Neither is reachable today - `lower_bounds` and `upper_bounds` are always
  assigned together in `_prepare_optimization`, and `build_bound_status()`
  returns no `'lower'`/`'upper'` status when bounds are absent - but both rely
  on invariants that nothing enforces.

## Version 0.20.4 (2026-07-26)

Build configuration only; no code changes.

### Fixed

- **CI `lint` job red since 2026-07-25 without any code change.** The project
  had no `[tool.ruff]` section in `pyproject.toml`, so it inherited whatever
  the newest `ruff` release happened to consider its default rule set - and CI
  installs `ruff` unpinned (`pip install ruff`). `ruff 0.16.0`, released
  2026-07-23, widened that default substantially, and the first push afterwards
  (`f0f55a5`, v0.20.0) turned the lint job red. Every commit since inherited the
  failure, including both release commits. The timeline confirms the cause: CI
  was green on `99b09f6` (2026-07-19, ruff 0.15.x) and has failed on every run
  after ruff 0.16.0 shipped.

  Reproduced locally with ruff 0.16.0: **578 findings**, none of them
  correctness issues. The bulk is typing modernization and import hygiene -
  `UP006` non-PEP585 annotations (193), `FA100` missing
  `from __future__ import annotations` (158), `I001` unsorted imports (57),
  `UP035` deprecated `typing` imports (56), `UP037` quoted annotations (35),
  `RUF022` unsorted `__all__` (28), `UP007` non-PEP604 unions (18).

  The rule set is now pinned explicitly to `["E4", "E7", "E9", "F"]` - which is
  what the project has in fact been linted against all along (ruff's historical
  default: pycodestyle errors plus pyflakes). No source file changed. Adopting
  the wider set remains possible, but it is a deliberate decision rather than
  something a dependency upgrade imposes: `UP007` wants `X | None`, which on the
  supported Python 3.9 floor requires `from __future__ import annotations` in
  every module - precisely what `FA100` is asking for.

### Known issues

- The CI `typecheck` job also reports a failure (exit code 2). It runs with
  `continue-on-error: true`, so it does not fail the workflow. Unchanged by
  this release; diagnosed and fixed in 0.20.5.

## Version 0.20.3 (2026-07-26)

### Fixed

- **Silently wrong `R_inf` when the top frequency decade has fewer than 3
  points.** The purely capacitive branch of `--ri-fit` extrapolates `Re(Z)`
  to `Im(Z) = 0` by fitting a 2nd-degree polynomial `Re = a*Im^2 + b*Im + c`
  to the points in the highest frequency decade. A 2nd-degree polynomial has
  3 coefficients, so with only 1 or 2 points the system is underdetermined
  and `numpy.polyfit` silently returns a least-norm solution that looks like
  a valid fit but is not one. Reproduced on synthetic data with a known
  answer: true `R_inf = 5.0` Ohm, 2 points in the top decade, reported
  `R_inf = 4.32` Ohm (a 14% error), with `fit_success=True`, an empty
  `warnings` list, and a `method` string claiming the polynomial fit
  succeeded. `R_inf` feeds directly into DRT analysis, so the wrong value
  propagated into downstream results without any indication something was
  off. The fit now checks the point count before calling `polyfit`: with
  fewer than 3 points it falls back to the measured highest-frequency
  `Re(Z)` (the same fallback already used when the polynomial overshoots),
  labels the result `capacitive_hf_too_few_points_*`, and appends a warning
  naming the point count. `R_inf_poly` and `poly_coeffs` are `None` in this
  case instead of holding values from a fit that was never actually
  performed.

## Version 0.20.2 (2026-07-25)

### Removed

- **`md2pdf.sh`** — local documentation build helper, not part of the
  toolkit and not distributed by `pyproject.toml`. Nothing in the repository
  referenced it: no README instructions, no CI step, no other script.
  Recoverable with `git show 8aade2b:md2pdf.sh`.

## Version 0.20.1 (2026-07-25)

Documentation only; no code changes.

### Fixed

- **Incorrect `--ri-fit` description in README.** It claimed the flag was
  "suitable for data with inductive loop" and used "R+L+K model fit", but the
  estimator has three branches and the R-L-K least-squares fit is only one of
  them: a zero crossing of Im(Z) inside the top decade is interpolated
  directly, and purely capacitive data are extrapolated to Im=0 by a
  2nd-degree polynomial in the Nyquist plane. The flag handles both inductive
  and capacitive high-frequency ends.

- **Undocumented default R_inf method.** The median of Re(Z) over the (up to
  5) highest-frequency points is what nearly every run uses, and it appeared
  in neither README nor any `doc/*.md`.

### Added

- **`doc/RINF_ESTIMATION.md`** — the R_inf estimation module had no
  documentation, unlike every other algorithm in the toolkit. Covers both
  methods, the decade selection window, all three `--ri-fit` branches with
  their guards, the warnings, the failure paths, and why `--ri-fit` makes the
  DRT report `Method: Preset value`.

## Version 0.20.0 (2026-07-25)

Completes the ponytail audit of 2026-07-25 (`doc/PONYTAIL_AUDIT_2026-07-25.md`)
with findings #2, #4-#7 — the group that changes public signatures. Backward
compatibility was deliberately not preserved.

### Changed

- **BREAKING: `estimate_rinf_with_inductance()` returns `(RLKFitResult, fig)`**
  instead of the 5-tuple `(R_inf, L, circuit, diagnostics, fig)`. The `L` slot
  was bound and dropped by both callers, `circuit` was documented as "Always
  None (legacy placeholder)", and `diagnostics` was a lossy re-projection of
  `RLKFitResult` into a dict (plus an `r_squared`/`R_squared` alias pair).
  Read the fields off the dataclass instead: `fit.R_inf`, `fit.L_nH`,
  `fit.R_squared`, `fit.n_points_used`, ... (audit finding #5).

  `RLKFitResult` now defaults its non-core fields and gains a
  `RLKFitResult.failed()` constructor, so the two failure paths (fewer than
  3 points, exception during the fit) return the *same type* rather than a
  differently shaped dict. `fit_success` distinguishes them; `R_inf` still
  falls back to the median of Re(Z) and stays usable.

- **Fixed as a side effect:** `run_rinf_estimation()` read
  `diagnostics['n_points_used']` unguarded, but the fit-exception fallback
  never set that key. The resulting `KeyError` was swallowed by a broad
  `except Exception` that reported "R_inf estimation failed" even though a
  usable R_inf had been computed. With a single return type the failure mode
  is gone; such a run now reports the R_inf and the warning.

### Removed

- **BREAKING: `use_voigt_fit` parameter** from `calculate_drt()`,
  `_estimate_r_inf()` and `estimate_rinf_with_inductance()`. It was ignored at
  the bottom of the chain and its only remaining effect was selecting the
  cosmetic label `method='voigt_fit'`. Since `use_rl_fit=True` was never
  passed by anything, the two flags together only chose between two strings.
  The CLI now passes `use_rl_fit=args.ri_fit`; `'voigt_fit'` is gone from
  `RinfEstimate.method` and from the CLI label table (finding #2).

- **BREAKING: `RinfEstimate.R_ct`, `.C_nF`, `.f_characteristic`** — always
  `None`, because the producer read diagnostics keys that never existed. The
  two print blocks guarded by them were unreachable, as was the
  `'voigt' in method` branch in the `--ri-fit` handler (finding #4).

- **BREAKING: ignored `verbose` parameter** from `fit_equivalent_circuit()`,
  `fit_circuit_diffevo()`, `fit_circuit_multistart()` and
  `estimate_rinf_with_inductance()`, plus the four call sites that still
  passed it. Logging level is the actual control. The `-v/--verbose` CLI flag
  is unaffected (finding #7).

- **BREAKING: `eis_analysis.utils.compat` module and the `np_trapz`
  re-export** from `eis_analysis.utils`. 32 lines (26 of them docstring)
  wrapping a 4-line import shim with a single consumer. The shim itself is
  still required — `numpy.trapezoid` only exists in NumPy >= 2.0 while the
  project floor is `numpy>=1.20` — so it now lives inline in its one consumer,
  `fitting/auto_suggest.py` (finding #6).

## Version 0.19.0 (2026-07-25)

### Removed

- **BREAKING: `plot_rl_fit_diagnostics()` removed from the public API.**
  The function and its module `eis_analysis/visualization/diagnostics.py`
  (224 lines) had zero call sites — it was superseded by `_plot_rlk_fit()`
  in `eis_analysis/rinf_estimation/rlk_fit.py`, which is what `--ri-fit`
  actually renders. Removed from `eis_analysis.__all__` and from
  `eis_analysis.visualization.__all__`; the stale `PYTHON_API.md` example
  (which documented a signature the function never had) is gone too.
  CLI behaviour is unaffected. Code importing this symbol must switch to
  the `--ri-fit` diagnostic figure returned by
  `estimate_rinf_with_inductance(..., plot=True)`
  (ponytail audit 2026-07-25 finding #1, see
  `doc/PONYTAIL_AUDIT_2026-07-25.md`).

- **`OnePerLineHelpFormatter`** from `eis_analysis/cli/parser.py` — 37 lines
  that never executed, because the parser sets `usage=` explicitly and the
  formatter's `_format_usage()` therefore always took its early-return
  branch. Replaced by `argparse.RawDescriptionHelpFormatter` (its base
  class); `--help` output is byte-identical (finding #3).

- **Unused parameter `refit_positive`** from `fit_voigt_chain_linear()`,
  whose own docstring described it as "Unused parameter for backward
  compatibility". No caller passed it (finding #9).

- **Dead constants** `DRT_TOLERANCE`, `GAMMA_MIN_REASONABLE` from
  `eis_analysis/drt/core.py` — both were annotated "Unused (pre-existing)"
  and were in no `__all__` (finding #10).

- **Unused parameters in private helpers** — `weighting` from
  `_prepare_optimization()` and `R_inf` from `_create_visualization()`;
  neither was read in the function body (finding #11).

- **Dead configuration constants** `DRT_PEAK_MIN_SPACING_DECADES`,
  `DEFAULT_R0_GUESS`, `DEFAULT_Q_N_GUESS` from
  `eis_analysis/fitting/config.py` (ponytail audit finding #1,
  see `doc/PONYTAIL_AUDIT.md`). None of them was referenced anywhere
  in the codebase.

- **Duplicate data files in repository root** `example_eis_data.csv` and
  `real_gamry_example.DTA` — identical copies remain in `example/`
  (ponytail audit finding #3).

- **`requirements.txt`** — duplicated the dependency list from
  `pyproject.toml`; install with `pip install .` instead
  (ponytail audit finding #4).

### Changed

- **CSV delimiter auto-detection simplified** — `_detect_delimiter()` is now
  a single `max(',', '\t', ';', key=header_line.count)`. Detection of comma,
  semicolon and tab headers is unchanged; a header containing none of the
  three now yields `','` (the documented default) instead of `'\t'`
  (ponytail audit 2026-07-25 finding #8).

- **CLI logging formatters consolidated** — the four per-level formatter
  classes in `eis_analysis/cli/logging.py` replaced by a single
  `PrefixFormatter` with a level-to-prefix table; output is unchanged
  (ponytail audit finding #2).

- **Docstring cleanup** — trimmed the oversized `sort_by_frequency`
  docstring and dropped a stale `AUDIT_REPORT.md` reference from
  `eis_analysis/fitting/config.py` (ponytail audit findings #5, #6).

### Added

- **Brug (2D) comparison estimate in oxide analysis.** For a dominant
  R||Q element, `analyze_oxide_layer()` now also evaluates the Brug
  (1984) effective capacitance `C = Q^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)`
  (surface/2D distribution of time constants) and the corresponding
  thickness, alongside the primary Hsu-Mansfeld (normal/3D) values.
  Rs is taken from the series path of the fitted circuit; without a
  series R the new `OxideAnalysisResult` fields (`capacitance_brug`,
  `capacitance_specific_brug`, `thickness_brug_nm`) stay `None`.
  The spread between the two estimates brackets the model uncertainty
  of the CPE-to-capacitance conversion. Theory discussion (2D vs 3D
  distributions, sensitivity to R, power-law model, Zr oxide context)
  added to `doc/OXIDE_ANALYSIS_GUIDE.md`.

- **Standalone script `oxide_permittivity.py`.** Computes oxide layer
  relative permittivity from effective capacitance, sample area, and known
  layer thickness (parallel plate capacitor model, using `EPSILON_0` from
  `eis_analysis.analysis.config`).

- **Lambda-probe peak stability diagnostics (`--lambda-probe`).** Opt-in
  DRT diagnostic inspired by (and extending) the lambda-sensitivity probe of
  Auto DRT Analyzer: the regularized NNLS problem is re-solved at
  `lambda* x 10^(+-0.5)` and `lambda* x 10^(+-1)` around the selected
  lambda*, and each peak of the main solution is tracked across the probe
  solutions by proximity in log10(tau). Each peak gets a persistence count,
  a position drift in decades, a relative R variation, and a verdict:
  `stable` (present in every probe, drift < 0.2 dec, R variation < 25%),
  `artifact` (missing from more than half of the probes), or `marginal`.
  Probe curves are drawn as thin overlays in the DRT figure. New module
  `drt/stability.py`; results in `DRTDiagnostics.stability`
  (`StabilityDiagnostics`, `PeakStability`, `LambdaProbePoint`); Python API
  via `calculate_drt(..., lambda_probe=True)`. Default behavior is
  unchanged (probe runs only when requested).

---

## Version 0.18.0 (2026-07-11)

### Changed

- **Log-space confidence intervals for positive scale parameters.** The
  symmetric Wald interval `p +/- t*se` produces physically meaningless
  negative bounds for positive scale parameters when the standard error is
  comparable to the value (e.g. `R = 14.2 +/- 7.4 -> CI [-0.43, 28.7]`,
  a resistance CI containing negative values). Scale parameters (R, C, Q, L,
  sigma, tau, ...) now get the delta-method CI on ln(p):
  `CI = (p/f, p*f)` with `f = exp(t*se/p)` -- always positive,
  multiplicatively symmetric, converging to the linear interval for
  `se/p -> 0`. Which parameters are scale-like is decided by the same
  criterion `classify_bound_status` already uses (positive lower bound and
  bounds spanning more than 6 decades), so the CPE exponent `n` (0.3-1.0)
  keeps the symmetric linear-space CI.
  - API: `compute_confidence_interval` gains an optional `log_scale` mask
    (default None = previous symmetric behavior); new helper
    `bounds.log_scale_ci_mask`; new constant `bounds.LOG_SCALE_BOUND_RATIO`.

- **Per-parameter inf on rank-deficient covariance.** Previously a
  rank-deficient Jacobian made every free parameter report `stderr = inf`,
  losing the information which parameter is the problem. Now only parameters
  whose direction overlaps the null space of J^T J (more than
  `NULLSPACE_OVERLAP_TOL = 1e-3` of the direction's norm) are reported as
  non-identifiable (`inf`); the rest get their variance from the
  Moore-Penrose pseudo-inverse restricted to the identifiable subspace. The
  covariance warning names the non-identifiable parameters. One
  non-identifiable parameter no longer destroys the confidence intervals of
  the others (`params_ci_95/99` are now per-parameter).

---

## Version 0.17.1 (2026-07-11)

### Fixed

- **Regression: parameter confidence intervals reported as `+/- inf` on
  essentially every fit (including the built-in demo).** The rank and
  condition-number checks in `compute_covariance_matrix` ran on the *raw*
  Jacobian. EIS parameters span many orders of magnitude (R ~ 1e5 Ohm,
  Q ~ 1e-6, n ~ 0.5), so the Jacobian columns differ by >10 decades and
  `cond(J)` exceeds the 1e10 threshold purely from those units — a scaling
  artifact, not genuine non-identifiability. Since audit #4 (0.16.x) reports
  rank-deficient covariance as infinite, this surfaced as `+/- inf` standard
  errors on well-determined parameters. The Jacobian is now column-scaled to
  unit norm (van der Sluis preconditioning) before the rank and conditioning
  analysis, making both scale-invariant; the covariance is inverted in the
  scaled space and rescaled via `(J^T J)^{-1} = D^{-1} (J_s^T J_s)^{-1} D^{-1}`.
  Genuinely rank-deficient Jacobians (linearly dependent columns) still report
  `inf`, preserving the audit #4 behavior. The reported `condition_number` is
  now the scale-invariant `cond(J_s^T J_s)`.

---

## Version 0.17.0 (2026-07-09)

### Added

- **Optional series capacitance in the Lin-KK model (`--kk-series-c`).**
  Blocking (capacitive) low-frequency behavior — e.g. two-electrode cells
  where the counter electrode contributes a series interfacial capacitance —
  is KK-compliant but cannot be represented by the Voigt chain: a series C
  has zero real part and diverging Z'', so the real-part Lin-KK fit produces
  imaginary residuals that grow toward low frequencies while the real fit
  stays good (observed on sample M136119: mean |res_imag| 8.7 % without the
  term, 0.13 % with it). The new term `1/(j*omega*C)` is fitted from the
  imaginary residual exactly like the series inductance L (Schonleber
  Lin-KK `add_cap`). Off by default so existing KK results stay comparable.
  - CLI: `--kk-series-c` flag; fitted C is printed in the KK summary, and a
    hint suggests the flag when the residual pattern matches (imaginary
    residuals dominate, real fit good).
  - API: `include_C` parameter on `kramers_kronig_validation`,
    `lin_kk_native`, `find_optimal_M_mu` and `estimate_R_linear`;
    new `capacitance` field on `KKResult`/`LinKKResult`.

### Breaking (API only, CLI unaffected)

- `estimate_R_linear` returns a 4-tuple `(elements, residual, L_value,
  C_value)` instead of a 3-tuple. `C_value` is never part of the
  `elements` array (unlike L).
- `find_optimal_M_mu` and `find_optimal_extend_decades` return a 6-tuple
  (appended `C_value`) instead of a 5-tuple.

## Version 0.16.23 (2026-07-08)

### Fixed (CLI parser audit 2026-07-08, finding P1)

- **`--multistart N` is no longer silently ignored.** With the default
  `--optimizer de`, passing `--multistart N` had no effect — the fit ran
  Differential Evolution and N was discarded without a word. Now
  `--multistart N` implies `--optimizer multistart`; combining it with
  an explicit `--optimizer de/single` is a clean argparse error instead
  of a silent surprise.
- **"0 = disabled" semantics removed.** The help text claimed
  `--multistart 0` disables multistart, but with `--optimizer multistart`
  the handler mapped 0 (and any negative N) to a hidden default of 16
  restarts. Non-positive N is now rejected by the parser; the default of
  16 is documented in the help text and README.
- Regression tests in `tests/test_cli_parser.py` (new file).

## Version 0.16.22 (2026-07-08)

### Fixed (fitting diagnostics audit 2026-07-03, finding D4)

- **Voigt chain linear fit no longer claims exact parameters.** The CLI
  handler (`--voigt-chain` without nonlinear refinement) filled
  `FitResult.params_stderr` with zeros, which reads as "known exactly"
  and made `params_ci_95` return CIs of +/- 0. The linear fit provides
  no uncertainty estimate, so stderr is now `inf` ("unknown" — the
  convention used by `covariance.py` and `diffevo.py` for missing
  covariance) and the CIs are honestly (-inf, inf). CLI output is
  unaffected (this path never printed per-parameter stderr); only API
  consumers of the returned `FitResult` see the change.
- Regression test in `tests/test_cli_integration.py`.

## Version 0.16.21 (2026-07-08)

### Fixed (fitting diagnostics audit 2026-07-03, finding D3)

- **Initial guess clipped into bounds is no longer silent.**
  `_prepare_optimization` clips out-of-bounds initial values (e.g. `R(0)`
  starts from the lower bound 1e-4) and recorded the indices in
  `OptimizationSetup.clipped_params`, but nothing propagated them — the
  fit silently started from a different point than the user specified.
  Each clipped parameter now produces a warning in
  `FitDiagnostics.warnings` (visible via `FitResult.all_warnings` and the
  CLI) plus a `logger.warning`, with the parameter label, original value
  and clipped value. With an explicit `initial_guess` override nothing is
  reported: the merge discards the clipped values, and an out-of-bounds
  override already fails loudly in `least_squares` ("x0 is infeasible").
- Regression tests in `tests/test_bounds_diagnostics.py`.

## Version 0.16.20 (2026-07-03)

### Fixed (Kramers-Kronig audit 2026-07-03, finding K2)

- **mu is now documented and presented as the Lin-KK stop value, not
  a data-quality metric.** By construction of the Lin-KK iteration the
  returned mu is below mu_threshold on every normal termination (it is
  the value at the first M where mu dropped below the threshold), yet
  the docs claimed the opposite: `PYTHON_API.md` said "Quality metric
  (>0.85 = good)" — making every successful KK run read as bad — and
  README said values below the threshold "indicate problematic data".
  Corrected in: `KKResult`/`LinKKResult` docstrings, `lin_kk_native` and
  `kramers_kronig_validation` Returns sections, README `--mu-threshold`
  description, PYTHON_API, and CLI `--mu-threshold`/`--voigt-mu-threshold`
  help texts. Data quality is judged by the residuals (`is_valid`).
- **Inverted threshold description fixed** (`find_optimal_M_mu`): "Lower
  values -> more conservative (fewer elements)" was backwards — the
  iteration runs while mu > threshold, so a lower threshold stops later
  and yields *more* elements.
- CLI log now gives the context inline:
  `KK: M=25, mu=0.7324 (Lin-KK stop, threshold 0.85), ...`; the residuals
  plot title labels the value as `stop mu=`.
- No numeric changes; documentation and presentation only.
  - Regression test in `tests/test_kramers_kronig.py` (CLI log carries the
    stop-value context; normal termination ends below the threshold).

## Version 0.16.19 (2026-07-03)

### Fixed (Kramers-Kronig audit 2026-07-03, finding K1)

- **R_s recovery in the imaginary-only fit was numerically broken**
  (`fitting/voigt_chain/fitting.py`, `estimate_R_linear` with
  `fit_type='imag'`). The real-part prediction used for the Boukamp R_s
  recovery was built from the *weighted* design matrix and "un-weighted"
  by multiplying with |Z| — valid only for unnormalized modulus weights,
  while the weights are normalized to mean 1 (and other weightings never
  matched at all). Present since the initial commit. Measured on
  KK-consistent synthetic data (true R_s = 100 Ohm, M=10):

  | weighting | before    | after fix |
  |-----------|-----------|-----------|
  | modulus   | -4513     | 98.7      |
  | uniform   | -36596    | 70.2*     |

  (*remaining deviation is tau-grid discretization at M=10; at M=15 all
  weightings recover R_s within 0.3 %.) The prediction is now computed
  from the unweighted design matrix (dividing the weights back out).
  **Not affected:** the default paths — KK validation (`fit_type='real'`)
  and the Voigt-chain CLI default (`complex`) — plus R_k, L and the mu
  metric. **Numeric change (correction):** results of
  `--voigt-chain --voigt-fit-type imag` and direct API calls with
  `fit_type='imag'`.
  - Regression tests in `tests/test_voigt_chain.py`: R_s recovery for all
    four weightings, the low-impedance (mOhm) case where the old bias
    looked plausible, and an end-to-end `fit_voigt_chain_linear`
    imag-fit error check.

## Version 0.16.18 (2026-07-03)

### Fixed (fitting diagnostics audit 2026-07-03, findings D1, D2, D5)

- **One criterion for "parameter at a bound"** (`fitting/circuit.py`,
  `fitting/bounds.py`): `FitDiagnostics.bounds_warnings`/`params_at_bounds`
  used to check "within 1 % of the bound value" while
  `FitResult.bound_status` (which suppresses CIs in the CLI) used the
  `classify_bound_status` criterion (1 decade on log scale for wide
  bounds), so the CLI could tag a parameter "[at lower bound — CI not
  meaningful]" without any corresponding warning. Both channels are now
  derived from a single `bound_status` vector built by the new shared
  helper `build_bound_status()` (also replaces the duplicated block in
  `fitting/diffevo.py`, whose diagnostics now carry the bounds fields
  too). **Behavioral change:** bounds warnings fire on the looser,
  CI-consistent criterion; fit numerics are unchanged.
- **Warning indices are now in full parameter space with labels**:
  messages like "Parameter 0 at lower bound" used free-space indices
  (fixed parameters excluded), so with fixed parameters the number
  pointed at the wrong parameter. Now:
  "Parameter n0 = 1.000e+00 near upper bound 1.0e+00".
  `FitDiagnostics.params_at_bounds` holds full-space indices.
- **Removed dead `check_bounds_proximity()`** (`fitting/bounds.py`) — a
  third, never-called implementation of the same check. **Breaking** only
  for direct API imports of this function; use `build_bound_status()` /
  `classify_bound_status()` instead.
- Tests added in `tests/test_bounds_diagnostics.py`: `build_bound_status`
  unit coverage (log/linear branches, fixed, no bounds), warning/status
  consistency regression, full-space index regression with a fixed
  parameter.

---

## Version 0.16.17 (2026-07-03)

### Fixed (oxide audit 2026-07-02, finding O4)

- **K element with R = 0 no longer crashes the circuit traversal**
  (`analysis/oxide.py`): `C = tau/R` raised `ZeroDivisionError` before the
  dominant-element filter could drop the element. Fitting was protected by
  parameter bounds, but a direct API call with such a circuit crashed. The
  element is now skipped with a warning.
- **Ambiguous parallel combinations are no longer silent**: circuits like
  `(R1 | R2 | C)` or `(R | C1 | C2)` used to silently take the last R and
  the last C/Q; a warning is now logged (behavior unchanged — the last
  element still wins).
- **Dead code removed**: the Z''-maximum and 1 kHz fallback branches of
  `_estimate_cpe_capacitance()` were unreachable (the only caller always
  passes R > 0 after the dominant-element filter), so the helper now
  implements only the Hsu-Mansfeld conversion with a required R. Unused
  `traverse()` parameters dropped; `Q`/`R` parameter names no longer shadow
  the circuit-element classes; `OxideAnalysisResult.element_params` is a
  copy instead of an internal mutable dict.
- **Numeric results are unchanged** for all previously working inputs.
  - Tests added to `tests/test_oxide.py`: K element traversal, Hsu-Mansfeld
    conversion accuracy, mixed Voigt+K dominance (audit priority 4), plus
    regressions for the K R=0 crash and the multiple-R/C warnings.

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
