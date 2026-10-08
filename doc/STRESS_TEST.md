# Stress test: results and known limits

`tests/stress.py` generates random, physically faithful spectra (six circuit
families, random frequency grids and noise) and checks invariants that must
hold whatever the truth is. Design: `doc/STRESS_TEST_PLAN.md`.

    python3 tests/stress.py                          # full run, ~20 min on 4 processes
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

### 2. Absolute PARAMETER_BOUNDS make a circuit fit unit-dependent

The default bounds (R 1e-4..1e10 Ohm, C 1e-15..0.1 F, ...) are absolute.
scipy's trust-region method scales its step by the distance to the bounds,
so they steer the fit even when no parameter comes near them. The same
spectrum in other units (k = 0.001 / 1000) gives a visibly different fit
in 628 / 778 of 1250 cases (Babs "meze"). With bounds scaled along with the
data the fit is exact (invariant B). `fit_equivalent_circuit(...,
bounds=...)` takes data-derived bounds, but the default stays absolute:
bounds derived from the spectrum were rejected (2026-10-07; one spectrum can
hold an R_s near zero and an oxide R near infinity, and a bound from |Z|
could cut off a physical optimum), and so were bounds scaled by the
electrode area (2026-10-08, `doc/AREA_SCALED_BOUNDS_PLAN.md`; the area is
often not recorded). This limit is accepted.

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
