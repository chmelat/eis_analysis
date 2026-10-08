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

Implementation: `hf_median()` in `eis_analysis/rinf_estimation/estimate.py`.
`--ri-fit` does not fall back to it (see "Fallback" below).

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
ten times more under 1 % noise; see `archive/AUDIT_ri_fit_2026-09-25.md`. The window
must hold at least `RINF_FIT_MIN_POINTS = 5` points (5 free parameters,
10 real residuals), otherwise the fallback below is used.

### Start values

| Parameter | Start |
|-----------|-------|
| R_s | lowest Re(Z) in the window (every term of the model has Re >= 0) |
| R_k | spread of Re(Z) in the window |
| tau = (R_k*Q)^(1/n) | 1/omega at the -Im(Z) maximum, where an arc peaks |
| n | 0.8 |
| L | Im(Z)/omega at f_max if the top point is inductive, else a reactance of 0.1 % of \|Z\| at f_max |

### Bounds

Upper bounds relative to the window's impedance, with
`RINF_BOUND_RANGE = 1e6`: R_s and R_k up to 1e6 x max|Z|, L such that its
reactance stays below that across the window, and Q such that 1/(Q omega^n)
stays within min|Z|/1e6 .. 1e6 x max|Z|, over the window and the whole n
range 0.3-1 (Q is the only parameter with a relative lower bound). Above
those an arc is open (or shorted) to better than 1e-6, under any
instrument's resolution. R and L are bounded below by 0, their physical
limit: a relative floor cut off what noise-free data still determine (an L
floor at 1e-6 |Z| biased a 0.05 Ohm R_s in front of a GOhm film to
0.064 Ohm). Start values are clipped into min|Z|/1e6 .. max|Z| x 1e6 (in
impedance), so none starts at 0. The absolute
`PARAMETER_BOUNDS` of the circuit fit made R_inf depend on the units: an open
arc's R_k ran into R <= 1e10 Ohm, and Z -> 1000 Z shifted R_inf on 51 % of
the stress test's random spectra (`tests/stress.py`).

A single fit from these start values is used. Multistart was measured worse:
under noise it finds degenerate minima with a lower residual (R_s at its lower
bound, the (R|Q) acting as a resistor) and takes 20x longer.

### Identifiability

The fitted R_s is used only when its standard error (from the fit's
covariance) satisfies

```
stderr(R_s) / R_s <= RINF_REL_STDERR_MAX = 0.05
```

Otherwise `R_inf` is the fallback below, and a warning gives the fitted value,
its stderr and the reason. Measured on the audit's synthetic set with 1 % noise,
determinable cases have 0.2-3.6 %, non-determinable ones 14 % and more; 5 %
leaves a ~2x margin on both sides.

The stderr is a **flag, not an error bar.** When the model does not describe
the data (e.g. two overlapping CPE arcs) it understates the real error: on
`example/example_eis_data.csv` it reports 23 % while the fitted R_s is off by
+153 %.

### Fallback: HF upper bound

When the fit does not determine R_inf, `R_inf = Re(Z)` at the highest
frequency with Im(Z) <= 0 (`method='hf_bound'`, value in `R_inf_hf`, its
frequency in `f_hf`). Every passive element adds Re >= 0 to R_s, so any
point is an upper bound of R_inf. A negative Re(Z) there (noise exceeding
the signal, or a lead artifact) bounds nothing, as R_s >= 0: `R_inf` is
then 0 and `R_inf_hf` keeps the measured value. On a capacitive top this is f_max, the
tightest one the data give: on an open arc Re(Z) still falls towards R_s as
the frequency rises. The 5-point HF median used before 0.41 reaches back
into the arc and overestimates more. Measured with 1 % noise:

| Case | HF median | Re(Z) at f_max |
|------|-----------|----------------|
| open CPE arc (D, Rs = 0.05 Ohm) | +3295 % | +2483 % |
| arc above f_max (A1) | +997 % | +996 % |
| `real_gamry_example.DTA` | 1402 Ohm | 826 Ohm |
| flat / pure R-L / Warburg end, fit flagged | ~0 % | -0.6 to -0.7 % |

The last row is the price: one noisy point instead of a median of five. On
these ends the fit is flagged only occasionally (4-10 of 20 noise seeds).
`min(median, Re(Z) at f_max)` gave the same numbers as Re(Z) at f_max alone.

**Inductive top.** Above the Im(Z) = 0 crossing the bound skips to the first
capacitive point. There a series L contributes nothing to Re(Z), but lead
artifacts (mutual inductance, stray capacitance of the cabling) can pull
Re(Z) below R_s, even below zero, which no passive model allows. On a
redoxED flow-cell sweep (`cell_EIS_1`, cycle 1, |Z| ~ 0.5 Ohm) Re(Z) falls
from 0.168 Ohm at the 89 kHz crossing to 0.003 Ohm at 446 kHz and -0.064 Ohm
at 500 kHz. Re(Z) at f_max handed the DRT 0.003 Ohm, and 32 % of R_pol piled
up at the fast end of the tau grid; the crossing gives 0.168 Ohm, R_pol
0.39 instead of 0.55 Ohm, and KK flags every point above it. This is the HFR
of the fuel-cell and flow-cell literature.

The price is a looser bound on a clean inductive spectrum the fit cannot
settle: open CPE arc D with L = 10 uH, 1 % noise, gives 2.22 Ohm at 40 kHz
instead of 1.29 Ohm at f_max (R_s = 0.05 Ohm, both off by thousands of
percent). With the arc far below the crossing the two differ by ~1e-4 Ohm.
Telling an artifact from a genuine arc above the crossing would need the
points the loader already drops for Re(Z) < 0. A spectrum inductive down to
its lowest frequency falls back to f_max.

### What `--ri-fit` cannot do

- **An arc lying entirely above f_max** looks, within the window, exactly like
  a flat high-frequency end. No method can tell the two apart from the data;
  the fit is flagged and the upper bound used.
- **A strongly open CPE arc** (phase of tens of degrees at f_max, low n) leaves
  R_s undetermined under realistic noise; again flagged, and the upper bound
  is then far off. Measure to higher frequencies.
- On clean flat or purely inductive ends the fit is sometimes flagged too
  (the (R|Q) term has nothing to fit); the upper bound it falls back to is
  within noise there.

### Failure paths

`estimate_rinf()` raises `ValueError` only for arrays of different shape or
with no finite point. Non-finite points are dropped with a warning. A window
with too few points, a fit that fails (`RuntimeError` from the fitter) or an
undetermined R_s all give `method='hf_bound'` with the reason in `warnings`;
`fit` is `None` in the first two cases.

---

## Interaction with the DRT

When `--ri-fit` runs, the CLI computes R_inf **first**, prints the fitted value
and the upper bound, and hands the chosen one to `calculate_drt()` as
`r_inf_preset`:

```
R-L-(R|Q) fit, 506-3.99e+04 Hz (20 points): R_inf = 1.082 +- 0.011 Ohm (1 %)
  L = 326 nH, fit error 1.1 %
HF upper bound: Re(Z) = 1.355 Ohm at 3.99e+04 Hz
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
est.method         # 'rlq_fit' | 'hf_bound'
est.R_inf_fit      # fitted R_s [Ohm], also when not used (None if no fit)
est.R_inf_stderr   # its standard error [Ohm]
est.R_inf_hf       # fallback upper bound: Re(Z) at the highest f with Im(Z) <= 0 [Ohm]
est.f_hf           # its frequency [Hz]; below f_max if the top is inductive
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
cabling), or a capacitive arc that has not closed. The DRT section compares
the value used with the HF median; a large spread says the top of the spectrum
is not flat and that the median would have been biased. A warning that R_inf
is not determined says the data cannot settle it: the value handed on is then
only an upper bound.
