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
| **R-L-K estimate** | `--ri-fit` | the highest frequency decade |

Both report their result through `RinfEstimate` (`eis_analysis/drt/results.py`),
whose `method` field is one of `'preset'`, `'median'`, `'rl_fit'`.

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
n_avg = min(5, max(1, N // 10))
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

Implementation: `_estimate_r_inf()` in `eis_analysis/drt/estimation.py`.

---

## `--ri-fit`: R-L-K estimate over the highest decade

Model:

```
Z(omega) = R_s + j*omega*L + R_k / (1 + j*omega*tau)
```

R_s is the estimate of R_inf, `j*omega*L` absorbs the cable/lead inductance
and the K element absorbs the beginning of the first arc, so the extrapolation
to omega -> infinity is not distorted by either.

### Data selection

`select_highest_decade()` (`data_selection.py`) takes all points with
`f >= f_max/10`. If that window happens to be empty, it falls back to the whole
dataset and records a warning. The selected points are classified by the sign
of Im(Z) into inductive (`Im > 0`) and capacitive (`Im < 0`) counts, which
decides the branch below.

### Three branches

The name "R-L-K fit" describes only the third branch. The estimator picks
whichever of the three is best supported by the data.

**1. Zero crossing** — `method='zero_crossing_*'`

Taken when the decade contains both inductive and capacitive points, i.e.
Im(Z) changes sign inside the measured range. This is the strongest case: the
real-axis intercept is *interpolated*, not extrapolated. The first sign change
is bracketed and Re(Z) is linearly interpolated at the crossing frequency:

```
t     = -Im_1 / (Im_2 - Im_1)
f_0   = f_1 + t * (f_2 - f_1)
R_inf = Re_1 + t * (Re_2 - Re_1)
```

`f_zero_crossing` is reported. No fit is performed, so `L` and `R_k` are 0.

**2. Capacitive extrapolation** — `method='capacitive_*'`

Taken when every point in the decade is capacitive (`n_inductive == 0`), so
there is no crossing to interpolate. A 2nd-degree polynomial is fitted in the
Nyquist plane, `Re = a*Im^2 + b*Im + c`, and extrapolated to `Im = 0`; the
constant term `c` is the estimate. Two guards apply:

| Condition | Action | Suffix |
|-----------|--------|--------|
| `c <= 0` | non-physical, clamp to 1 Ohm | `near_zero_fallback` |
| `c > Re(Z) at f_max` | extrapolation overshoots, use the top-frequency point | `hf_fallback` |
| otherwise | accept `c` | `polynomial` |

Both `R_inf_poly` (raw extrapolation) and `R_inf_hf` (Re(Z) at f_max) are
reported, so the two can be compared.

**3. R-L-K linear fit** — `method='rlk_linear_*'`

The remaining case, in practice purely inductive data. tau is estimated from
the data and clamped to the decade's own range, then R_s, R_k and L are solved
by linear least squares (`estimate_R_linear` from the Voigt-chain solver,
modulus weighting, non-negative R). This is the only branch that produces a
non-zero `L`.

### Warnings

- `L_nH > max_L_nH` (default 1000 nH) — implausibly large inductance
- `L < 0` — non-physical negative inductance
- above 500 nH, `_estimate_r_inf()` adds its own note when the estimate feeds
  the DRT

An estimate carrying warnings is still returned; nothing is silently rejected.

### Failure paths

`estimate_rinf_with_inductance()` never raises. Fewer than 3 points, or an
exception inside the fit, produce an `RLKFitResult` with `fit_success=False`,
`method='fallback_insufficient_data'` or `'fallback_fit_error'`, and `R_inf`
falling back to `median(Re(Z))` over all points. Callers distinguish the cases
via `fit_success`, not by a different return shape.

---

## Interaction with the DRT

When `--ri-fit` runs, the CLI computes R_inf **first** and hands the number to
`calculate_drt()` as `r_inf_preset`. The DRT therefore reports

```
Method: Preset value
```

rather than `R-L fit`, because from its point of view the value was supplied
from outside. The R-L-K diagnostic figure is saved separately as
`<prefix>_ri_fit.png`. The `'rl_fit'` label appears only when `calculate_drt()`
is called directly with `use_rl_fit=True` from the Python API.

---

## Python API

```python
from eis_analysis.rinf_estimation import estimate_rinf_with_inductance

fit, fig = estimate_rinf_with_inductance(frequencies, Z, plot=True)

fit.R_inf          # estimate [Ohm]
fit.method         # which branch ran
fit.behavior       # 'purely_capacitive' | 'purely_inductive' | 'mixed_with_crossing'
fit.R_squared      # quality over the fitted window
fit.n_points_used  # points in the highest decade
fit.L_nH           # inductance [nH], 0 outside the R-L-K branch
fit.warnings
fit.fit_success    # False -> R_inf came from the median fallback
```

Via the DRT entry point:

```python
from eis_analysis.drt import calculate_drt

result = calculate_drt(frequencies, Z, use_rl_fit=True)   # R-L-K estimate
result = calculate_drt(frequencies, Z, r_inf_preset=12.5) # supply it yourself
```

Full reference: [PYTHON_API.md](PYTHON_API.md)

---

## Choosing a method

Use the default HF median when the Nyquist plot already runs into the real
axis at the top of your frequency range — it is the common case and adds no
assumptions.

Reach for `--ri-fit` when the high-frequency end is still curving: an
inductive tail (Im(Z) > 0 at f_max, typical above ~100 kHz with ordinary
cabling), or a capacitive arc that has not closed. Compare the two numbers —
the CLI prints the median alongside the fitted value with their difference in
percent. A large spread is itself the diagnostic: it says the top of the
spectrum is not flat, and that the median would have been biased.
