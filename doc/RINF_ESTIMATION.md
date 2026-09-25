# R_inf Estimation

How the toolkit determines the high-frequency (ohmic) resistance R_inf — the
real-axis intercept that the DRT kernel subtracts before solving for
gamma(tau), and the R_s of an equivalent circuit.

Module: `eis_analysis/rinf_estimation/`

---

## Two methods

| Method | When | Data used |
|--------|------|-----------|
| **HF median** | default | up to 5 highest-frequency points |
| **R-L-(R\|Q) fit** | `--ri-fit` | the top two frequency decades |

The fit returns `RinfResult` (`eis_analysis/rinf_estimation/estimate.py`); the
DRT records the value it used in `RinfEstimate` (`eis_analysis/drt/results.py`),
whose `method` is `'preset'` or `'median'`.

### Why R_inf matters

The DRT model is

```
Z(omega) = R_inf + integral gamma(tau) / (1 + j*omega*tau) d ln(tau)
```

R_inf is subtracted from the data before the regularized NNLS solve, so an
error in R_inf does not cancel out — it is absorbed by gamma(tau), typically
as a spurious peak at the short-tau end or as a global offset in the
reconstruction. Overestimating R_inf can also drive Im(Z) of the residual
positive, which the non-negative solver cannot represent at all.

---

## Default: HF median

Used whenever `--ri-fit` is not given.

```
n_avg = min(HF_MEDIAN_MAX_POINTS, max(1, N // 10))     # HF_MEDIAN_MAX_POINTS = 5
R_inf = median(Re(Z) over the n_avg highest frequencies)
```

`N` is the number of data points. The median (rather than the mean) makes the
estimate insensitive to a single noisy top-frequency point, which is where
instrument noise is usually worst.

**Assumption:** the impedance has already flattened onto the real axis at the
top of the measured range, i.e. Im(Z) -> 0. This holds for many datasets and
costs nothing to compute.

**When it fails:** if the spectrum still has a significant imaginary part at
f_max — an inductive tail from the cabling, or a charge-transfer arc that is
not yet closed — the median of Re(Z) is biased, because Re(Z) has not reached
its limit. That is the case `--ri-fit` addresses.

Implementation: `hf_median()` in `eis_analysis/rinf_estimation/estimate.py`,
shared by the DRT and `--ri-fit`.

---

## `--ri-fit`: R-L-(R|Q) fit over the top two decades

Model:

```
Z(omega) = R_s + j*omega*L + R_k / (1 + R_k*Q*(j*omega)^n)
```

R_s is the estimate of R_inf. `j*omega*L` absorbs the cable/lead inductance and
the (R|Q) element the high-frequency end of the first arc, including a
depressed (CPE) one, so neither distorts the extrapolation to omega -> infinity.
One nonlinear fit (`fit_equivalent_circuit`, modulus weighting, the standard
parameter bounds, so L >= 1 pH and 0.3 <= n <= 1) covers inductive,
capacitive and mixed high-frequency ends alike.

### Window

All points with `f >= f_max / 10**RINF_FIT_DECADES`, `RINF_FIT_DECADES = 2`.
One decade resolves an arc just above f_max slightly better but scatters about
ten times more under 1 % noise; see `AUDIT_ri_fit_2026-09-25.md`. The window
must hold at least `RINF_FIT_MIN_POINTS = 5` points (5 free parameters,
10 real residuals), otherwise the HF median is used.

### Start values

| Parameter | Start |
|-----------|-------|
| R_s | lowest Re(Z) in the window (every term of the model has Re >= 0) |
| R_k | spread of Re(Z) in the window |
| tau = (R_k*Q)^(1/n) | 1/omega at the -Im(Z) maximum, where an arc peaks |
| n | 0.8 |
| L | Im(Z)/omega at f_max if the top point is inductive, else 1 nH |

A single fit from these start values is used. Multistart was measured worse:
under noise it finds degenerate minima with a lower residual (R_s at its lower
bound, the (R|Q) acting as a resistor) and takes 20x longer.

### Identifiability

The fitted R_s is used only when its standard error (from the fit's
covariance) satisfies

```
stderr(R_s) / R_s <= RINF_REL_STDERR_MAX = 0.05
```

Otherwise `R_inf` is the HF median, and a warning gives the fitted value, its
stderr and the reason. Measured on the audit's synthetic set with 1 % noise,
determinable cases have 0.2-3.6 %, non-determinable ones 14 % and more; 5 %
leaves a ~2x margin on both sides.

The stderr is a **flag, not an error bar.** When the model does not describe
the data (e.g. two overlapping CPE arcs) it understates the real error: on
`example/example_eis_data.csv` it reports 23 % while the fitted R_s is off by
+153 %.

### What `--ri-fit` cannot do

- **An arc lying entirely above f_max** looks, within the window, exactly like
  a flat high-frequency end. No method can tell the two apart from the data;
  the fit is flagged and the median used.
- **A strongly open CPE arc** (phase of tens of degrees at f_max, low n) leaves
  R_s undetermined under realistic noise; again flagged, median used, and the
  median itself is then far off. Measure to higher frequencies.
- On clean flat or purely inductive ends the fit is sometimes flagged too
  (the (R|Q) term has nothing to fit); the median it falls back to is exact
  there.

### Failure paths

`estimate_rinf()` raises `ValueError` only for arrays of different shape or
with no finite point. Non-finite points are dropped with a warning. A window
with too few points, a fit that fails (`RuntimeError` from the fitter) or an
undetermined R_s all give `method='hf_median'` with the reason in `warnings`;
`fit` is `None` in the first two cases.

---

## Interaction with the DRT

When `--ri-fit` runs, the CLI computes R_inf **first**, prints the fitted value
and the HF median, and hands the chosen one to `calculate_drt()` as
`r_inf_preset`:

```
R-L-(R|Q) fit, 506-3.99e+04 Hz (20 points): R_inf = 1.082 +- 0.011 Ohm (1 %)
  L = 326 nH, fit error 1.1 %
HF median (5 points): R_inf = 1.446 Ohm
Using R_inf = 1.082 Ohm (fit)
```

The DRT stage then skips its own `R_inf estimation` section and states the
value where it uses it, with a comparison against the HF median. The diagnostic
figure is saved as `<prefix>_ri_fit.png`.

---

## Python API

```python
from eis_analysis.rinf_estimation import estimate_rinf
from eis_analysis.visualization import plot_rinf_fit

est = estimate_rinf(frequencies, Z)

est.R_inf          # value to use [Ohm]
est.method         # 'rlq_fit' | 'hf_median'
est.R_inf_fit      # fitted R_s [Ohm], also when not used (None if no fit)
est.R_inf_stderr   # its standard error [Ohm]
est.R_inf_median   # HF median [Ohm]
est.fit            # FitResult of R-L-(R|Q): params_opt = [R_s, L, R_k, Q, n]
est.f_window, est.Z_window  # data of the fit window
est.warnings

fig = plot_rinf_fit(est)

from eis_analysis.drt import calculate_drt
result = calculate_drt(frequencies, Z, r_inf_preset=est.R_inf)
```

Full reference: [PYTHON_API.md](PYTHON_API.md)

---

## Choosing a method

Use the default HF median when the Nyquist plot already runs into the real
axis at the top of your frequency range — it is the common case and adds no
assumptions.

Reach for `--ri-fit` when the high-frequency end is still curving: an
inductive tail (Im(Z) > 0 at f_max, typical above ~100 kHz with ordinary
cabling), or a capacitive arc that has not closed. The CLI prints both numbers.
A large spread says the top of the spectrum is not flat and that the median
would have been biased; a warning that R_inf is not determined says the data
cannot settle it either way.
