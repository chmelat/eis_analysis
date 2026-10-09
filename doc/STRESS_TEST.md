# Stress test: results and known limits

`tests/stress.py` generates random, physically faithful spectra (six circuit
families, random frequency grids and noise) and checks invariants that must
hold whatever the truth is. Design: `doc/STRESS_TEST_PLAN.md`.

    python3 tests/stress.py                          # full run, ~40 min on 4 processes
    python3 tests/stress.py --check                  # full run against tests/stress_baseline.json
    python3 tests/stress.py --update-baseline        # record a triaged full run as the baseline
    python3 -m pytest tests/ -m stress               # smoke: first 3 cases per family, ~1 min
    python3 tests/stress.py --family oxide --index 37 -v   # replay one case

A case is reproducible from `family/index` alone (separate random streams per
purpose, one BLAS thread per process), so every number below names its seed.

## Invariants implemented

| Invariant | Meaning |
|---|---|
| A | No exception, no NaN/Inf, no empty result |
| B-, B+ | Z -> k*Z with k = 0.001, 1000 scales every result by its power of k |
| Babs-, Babs+ | Same for the fit with the library's absolute bounds; a difference is the known limit, counted as "meze" |
| Crev, Cmix | Reversed and shuffled points give the same result |
| K | A second run is bit-identical |

Consistency invariants (step 3, `tests/stress_consistency.py`) compare a
result with the truth the case was generated from. A threshold check reports
measured / allowed, with the allowance factor x sigma + floor x |Z_true|
(the case's noise plus the method's own error on exact data):

| Invariant | Meaning |
|---|---|
| D | Lin-KK and Z-HIT reproduce the (KK-consistent) spectrum within noise + floor at every point |
| Eneg | DRT: gamma >= 0 |
| Erec | DRT reconstruction, RMS over the points (rc family, closed high-frequency end) |
| Edc | DRT: R_inf + R_pol between Re Z(f_min) and the true DC resistance (closed low-frequency end) |
| Epeak | DRT: each true tau within 0.15 decade of a peak (rc, noise-free, arcs >= 5 % and >= 1 decade apart) |
| F1 | A noise-free fit started at the truth stays there |
| F2 | Multistart from truth x U(0.3, 3) ends no worse than the truth; a worse end is a local minimum, a rate ("lokmin") |
| F3 | Rate: truth inside the reported 95 % CI (proportional noise, well-conditioned fits) |
| F4 | A parameter at a bound is reported |
| G | DE (CLI defaults) from truth x 10^U(-2, 2), exponents uniform in their bounds, ends no worse than the truth (cost x (1 + 1e-3) + floor); every 5th case; a worse end is a failure |
| Irange | `R_inf_range` holds the true Rs (cases without L) |
| Ihf | R_inf <= Re Z_true(f_max) + 1 % + 3 sigma (passivity) |
| Iclosed | R_inf within 5 % + 3 sigma of Rs on a closed high-frequency end (phase > -5 deg) |
| Ifit | Rate: the window fit determines R_s on a closed high-frequency end |
| M | Rate: the n(f) map through the CLI path within 2x its uncertainty of the noise-free map |

Analyses: `calculate_drt` (auto lambda), `kramers_kronig_validation`,
`zhit_validation`, `estimate_rinf`, and one LM `fit_equivalent_circuit` of the
true circuit from truth x U(0.3, 3), and `local_exponent` (the n(f) map,
with the true Rs and L subtracted, so the map is tested and not
`estimate_rinf`; n compared only where the library marks it `valid`). For B the fit's bounds scale with Z
(PARAMETER_BOUNDS x k^power), so the optimizer alone is under test.

Tolerances: relative 1e-6, plus absolute floors where the value is rounding
noise (residuals 1e-9 |Z|; a fitted parameter 1e-3 of its stderr, on
noise-free data the change that moves the fit by 1e-9 |Z|).

## Full run 2026-10-07 (v0.57.0)

1250 cases (250 per family), 1105 s on 4 processes, 3.5 s per case
(max 72.9 s, oxide). 146 failed checks out of 40000.

- A and K: no failure in any analysis.
- DRT and Z-HIT: no failure of any invariant.
- R_inf: 2 failed checks, see the limits below.
- Lin-KK: 114 failed checks, all on noise-free spectra (known limit 1).
- Fit: 30 failed checks, almost all in ill-posed cases (known limit 3);
  Babs "meze" in 624 (k = 0.001) and 776 (k = 1000) of 1250 fits (known limit 2).

## Extension 2026-10-08: `anomalous` family and the n(f) map

Family `anomalous` (id 6): half Randles with Wa or Wat, half the ZrO2
layer model Rs-(G|Wa|Q|C) or Rs-(Wat|Q|C), drawn by its crossover
frequencies (`doc/STRESS_TEST_PLAN.md`). 1500 cases (250 per family),
1455 s on 4 processes, 3.9 s per case (anomalous 3.6 s, max 8.4 s).

- The five existing families: the same 146 failed checks as on 2026-10-07.
- n(f): 9000 checks (A, B, C, K on all six families), none failed.
- anomalous, all analyses: A and K no failure; DRT, Z-HIT, R_inf no
  failure; Lin-KK 22 failed checks in 11 of 68 noise-free cases (known
  limit 1); fit 14 failed checks in 5 cases (known limit 3). The 125
  ZrO2 cases, 7-8 parameters each, have none.

## Full run 2026-10-09: consistency invariants (step 3)

1500 cases (250 per family), 1459 s on 4 processes, 3.9 s per case (max
80.8 s), 70646 checks.

- A and K: no failure in any analysis.
- B and C: 186 failed checks, all in known limits 1, 3 and 4 (Lin-KK on
  noise-free data 138 in 52 cases, fit 46, R_inf 2 on rc/36). Four more
  than on 2026-10-08: constant noise is now capped at SNR >= 10 at the
  weakest point, which redrew those cases.

### Calibration

The first run used provisional thresholds (1171 failed consistency checks,
724 of them Z-HIT D). The thresholds were then set from this run's
per-point residuals with ~2x margin: the floor at 2x the method's largest
error on exact data, the noise factor at 2x what the noisy cases need on
top of that floor.

| Check | Exact data (largest, 99 %) | Noisy: factor needed (99 %, largest) | Threshold |
|---|---|---|---|
| Lin-KK D | 6.3e-3 \|Z\|, 3.4e-3 \|Z\| | 4.4, 9.0 (of the largest relative noise) | 10 x max(sigma/\|Z\|) x \|Z_i\| + 1.3e-2 \|Z_i\| |
| Z-HIT D | 0.15 \|Z\| (rc/83), 3.2e-2 \|Z\| | 9.3 sigma, 24 sigma (rc/41) | 25 sigma + 6e-2 \|Z_i\| |
| DRT Erec | RMS 7.0e-3 \|Z\| | largest ratio 0.60 at the threshold | RMS within 2 sqrt(2) sigma + 1.5e-2 \|Z_i\| |

Lin-KK's allowance uses the spectrum's largest relative noise, not the
point's own (known limit 6); for proportional noise the two are the same.
Edc (largest 0.28 of allowed), Ihf (0.91), Iclosed, Epeak and F are
unchanged.

### Results with the calibrated thresholds

Verification run 2026-10-09 (1461 s): the same cases, bit-identical
B/C/K results, and exactly the failures predicted from the first run's
residuals:

- D: Lin-KK 0 of 1500 failed (provisional: 191); Z-HIT 1 of 1500, rc/83
  (provisional: 724; known limit 5).
- DRT: Eneg 0/1500, Erec 0/173 (provisional: 29), Edc 0/841, Epeak 1/32
  (rc/28, known limit 8).
- Fit: F1 0/353, F4 0/230. F2 ended in a local minimum in 35 of 1500
  fits (2.3 %: anomalous 12, cpe 10, diffusion 8, blocking 2, oxide 2,
  rc 1), none a failure. F3 coverage of the nominal 95 % CI: 0.897, 0.899,
  0.904 at 0.1, 1, 3 % noise (known limit 9).
- R_inf: Ihf 0/1417, Iclosed 1/737 (cpe/85), Irange 39/1417 (both known
  limit 7). Ifit: the window fit determines R_s on a closed end in 0.989,
  0.968, 0.802, 0.635 of the cases at 0, 0.1, 1, 3 % noise.
- n(f), M: 0.975, 0.945, 0.909 of the points within 2x their uncertainty
  at 0.1, 1, 3 % noise; largest |dn| 0.075, 0.084, 0.13.

The Ifit rate decides the `--ri-fit` default: as the default it would
report "R_inf not determined" on 20 % (1 % noise) to 36 % (3 % noise) of
spectra whose high-frequency end is closed, so `--ri-fit` stays opt-in
until that false alarm is silenced.

## Invariant G and the baseline, 2026-10-09

G fits with DE from a start up to two decades off the truth, as a rough
guess for an unfamiliar sample can be (`de_start`), on every 5th case: 300
fits, ~11 s each, which lengthens a full run to ~40 min (2465 s). DE takes
the start as one member of its population, so the first version, started
like F2 at truth x U(0.3, 3), tested little more than its polish; it too
passed 150 of 150.

- G: 300 of 300 passed. On noisy spectra DE ends at 0.95-0.998x the
  truth's cost (it fits part of the noise), on noise-free ones at a
  relative cost of ~1e-13.
- All other checks: the same 228 failures as the verification run, none
  new, none gone.

`tests/stress_baseline.json` holds these 228 failures and the aggregate
rates (F2 local minima, F3, Ifit, M per noise level, with the number of
cases). `--check` fails on a failure not in it, and on a full run on a rate
that moves the wrong way by more than 3 binomial standard errors over the
cases (the points of one case move together, so they do not count
separately). The smoke test (`tests/test_stress_smoke.py`, marker `stress`)
runs the first 3 cases of every family through `--check` in a subprocess,
which keeps one BLAS thread as invariant K needs; 66 s. Check of the check:
R_inf x 1.1 injected into `estimate_rinf` gave 20 new failures (Ihf,
Iclosed, Irange) and exit code 1.

## Bugs found and fixed

The first full run (before the fixes) had 1190 failed checks.

- **R_inf depended on the units** (641 of 1250 cases under B, by more than
  10 % on 7 %; oxide/10: 1.2e9 instead of 1.6e3 Ohm at k = 1000). The window
  fit used the absolute PARAMETER_BOUNDS, which an open arc's R_k ran into.
  Fixed with upper bounds relative to the window's |Z| (`RINF_BOUND_RANGE`)
  and lower bounds 0 for R and L, see `doc/RINF_ESTIMATION.md`.
- **LM fits stopped at a different point in other units**, even for an exact
  binary scaling Z/2 (rc/24: 79x apart). scipy's xtol and gtol mix the
  parameters' and the cost's units, and an x_scale floor of 1e-10 rescaled
  small capacitances differently. Fixed by optimizing x/|x0| on a
  dimensionless residual with gtol disabled (`least_squares_normalized`,
  `compute_residual_weights`); with scaled bounds 300 of 300 fits are now
  identical under Z/2, before 42 were not. The code review of this fix
  found the absolute gtol stopping fits early on a parameter small relative
  to |Z| (a noise-free 0.05 Ohm R_s in front of 3e3 Ohm gave R_inf = 243
  Ohm); gtol is now off.

Regression tests: `tests/test_rinf_estimation.py::test_rinf_does_not_depend_on_units`,
`tests/test_bounds_diagnostics.py::test_fit_does_not_depend_on_units`.

Step 3 (2026-10-08), on points where the noise exceeds the signal; details
in CHANGELOG.md (Unreleased):

- **Z-HIT unwrapped the phase**: noise-swamped points got +-2 pi offsets,
  which the phase integral carried on (diffusion/104: |Z| up to 1e35).
- **`estimate_rinf` returned a negative R_inf** (rc/196: -1.9e4 Ohm) and
  its bound ignored the noise; now clipped at 0, with `R_inf_upper` and
  `R_inf_range`.
- **The n(f) map marked points determined when R_inf was only a bound**
  (oxide/119: 0.2 off at a claimed 0.02); it now takes `R_inf_range`.

## Known limits

### 1. Lin-KK on noise-free data: rounding decides

41 of the 285 noise-free cases fail B or C in Lin-KK; no case with noise
does. On exact data the residuals are rounding noise (~1e-9 |Z|) and the
least-squares system with up to 50 elements is ill-conditioned, so a
rescaling or reordering of the points moves:

- residuals by up to ~4e-8 |Z| (absolute),
- individual element resistances by up to 1e-3 relative (cpe/161: 9.1e-4),
- M by 1-3 in 3 cases, with mu following: rc/51 (M 40 -> 39..42, mu 0.54
  -> 0.56), rc/206 (48 -> 46, mu 0.71 -> 0.83), rc/123 (50 -> 47 at the
  max_M cap, mu 0.60 -> 0.82).

`anomalous` (2026-10-08): 11 of 68 noise-free cases, elements up to
1.8e-4 relative (anomalous/241), residuals up to 4e-9 |Z| (anomalous/79),
M unchanged.

None of this is visible on measured data, where noise is at least 1e-5 |Z|.

### 2. Absolute PARAMETER_BOUNDS: a fit depends on where Z sits within them

The default bounds (R 1e-4..1e10 Ohm, C 1e-15..0.1 F, ...) are physical
limits in base units: a spectrum in other units is converted before the
fit, so Z x 1000 is a different system, not the same one in other units.
scipy's trust-region method scales its step by the distance to the bounds,
so they steer the fit even when no parameter comes near them, and Z x k
with k = 0.001 / 1000 gives a visibly different fit in 716 / 893 of 1500
cases (Babs "meze", 2026-10-09). This is intended (2026-10-07). What must
not happen, the algorithm itself depending on the magnitude of Z (a mOhm
battery against a GOhm oxide), is invariant B with the bounds scaled
along, and there the fit is exact. `fit_equivalent_circuit(...,
bounds=...)` takes other bounds. Bounds derived from the spectrum were
rejected (2026-10-07; one spectrum can hold an R_s near zero and an oxide
R near infinity, and a bound from |Z| could cut off a physical optimum),
and so were bounds scaled by the electrode area (2026-10-08,
`doc/AREA_SCALED_BOUNDS_PLAN.md`; the area is often not recorded).

### 3. Ill-posed fits are path-dependent

30 failed fit checks in 11 cases. All but two are tight (neighbouring
arcs < 1 decade apart) or weak (smallest arc < 5 % of the summed arc R):
parameters with stderr far above their value, some at a bound, move along
a flat valley at an unchanged fit error (cpe/27, cpe/173, cpe/193,
cpe/222, oxide/213, oxide/223, rc/202), or a parameter differs at 1-5e-6
relative, the LM stopping precision (cpe/213, oxide/177). The two others:

- diffusion/53: a failed fit (82 % error), only its stderr moves (8 %);
- diffusion/203: one parameter at 1.0e-6 relative, the stopping precision.

`anomalous` (2026-10-08): 14 failed checks in 5 Randles cases with Wa at
1-3 % noise (anomalous/10, 91, 101, 114, 229). The Wa tail is lost in the
noise: R_W and tau_W have relative stderr 10 to 1e4, the fit ends far
from the truth at an unchanged fit error, and stderr (up to 20 %) or a
parameter (up to 9 %, anomalous/101) moves along the valley. The family
has no arc classes (n/a in the class table).

By class (B-, B+, Crev, Cmix together): tight 11, medium 14, loose 2, weak
20, normal 7 (the two axes overlap; see the runner's class table).

### 4. R_inf window fit: stopping precision

2 failed checks: rc/36 (0.1 % constant noise) shifts by 7.5e-5 under
k = 1000 and reversal. The window fit's R_s has a stderr of up to 5 % by
construction (`RINF_REL_STDERR_MAX`); the shift is far inside it.

### 5. Z-HIT: fast phase changes and noise gain

On exact data Z-HIT is within 3.2 % of |Z| in 99 % of the cases. The
second-order term assumes the phase varies slowly in ln omega; where it
swings within a few points it does not hold: rc/83 (R_s 0.0115 Ohm in
front of a 3.8 Ohm arc, with L: the phase turns from capacitive to
inductive) 15 %, rc/64 (35 points, with L) 5.4 %.

The phase derivative is taken unsmoothed, which multiplies the phase noise
by ~gamma / d(ln omega) (gamma = -pi/6). Beyond the floor a point is off
by up to 11.7 sigma with proportional noise and 24 sigma with constant
noise (rc/41): at 3 % noise, 10-30 % of |Z| at single points. At noise of
1 % and more, Z-HIT's pointwise residuals say little; a smoothed phase
before the derivative (Ehm et al.) would lower this, not implemented.

### 6. Lin-KK under constant noise

Lin-KK fits relative residuals (weights 1/|Z|), so the noisiest point
relative to its |Z| decides how far the mu criterion lets M grow. With
constant noise at SNR 10 on the smallest |Z|, the residuals at large |Z|
stay at 1-8 % of |Z|, up to 5e5 sigma of the local noise (oxide/123:
5.2 % at 1.35e6 x min|Z|, M = 17). Lin-KK reports a consistent spectrum
as consistent, but with its accuracy set by the worst relative noise.

### 7. R_inf_range misses R_s on open high-frequency ends

39 of 1417 cases without L: the window fit (`rlq_fit`) is accepted
(stderr <= 5 %), but R_inf is off by 6-34 %, 36 times below R_s and 3
times above (cpe/101 +12 %, cpe/220 +18 %, oxide/156 +29 %). All on
open ends (phase at f_max -4 to -82 deg) at low noise: 0 % (16 cases),
0.1 % (17), 1 % (6), none at 3 %, where the stderr and so the range are
wider. The R-L-(R|Q) window model differs from the arc, CPE or Warburg
that continues above f_max; that model error is not in the stderr, and
the range's +-5 % floor does not cover it. One of them (cpe/85, phase
-4.1 deg) also fails Iclosed. The n(f) map then counts points near f_max
as determined. A bound of |R_inf - R_s| by |Im Z(f_max)| does not hold
(exceeded in 179 of 878 window fits). Open.

### 8. DRT: R_inf on open high-frequency arcs, small arcs

The DRT's own R_inf is the median of the high-frequency Re Z, which
overestimates R_s where the arc is open at f_max (rc/149: 3469 Ohm, the
window fit 1381, R_s 1383); gamma >= 0 cannot make up for it, so Erec is
checked only on a closed high-frequency end. Epeak: rc/28 misses an arc
of 3.1 kOhm (6 % of the arcs, but 0.17 % of |Z| beside R_s = 1.8 MOhm)
by 1.95 decades: regularization smooths it out.

### 9. Fit statistics: local minima and CI coverage

A multistart from truth x U(0.3, 3) ends in a local minimum in 2.3 % of
the fits, most in `anomalous`, `cpe` and `diffusion`; all identifiability
classes have some (tight 5, medium 5, loose 3, weak 8, normal 5, no
class 22). The linearized 95 % CI holds the truth in 0.90 of the
parameters at every noise level, below the nominal 0.95.
