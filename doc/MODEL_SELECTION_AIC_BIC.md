# Choosing between circuits - an intuitive introduction

This document explains **why a better fit does not mean a better circuit**, and
what AIC and BIC do about it. The goal is that after reading it you can look at
the comparison table the toolkit prints and know what it is telling you.

CLI usage is in [README.md](../README.md#comparing-candidate-circuits). The
implementation is `compute_information_criteria` in
`eis_analysis/fitting/diagnostics.py`, the ranking in
`eis_analysis/cli/handlers/model_comparison.py`.

---

## 1. The problem: more parameters always fit better

Classic workflow: you propose a circuit, fit it, look at the error. Try a
second circuit, look at its error. Pick the lower one.

That last step is where it breaks. Add an element and the error almost always
drops - it has to. The bigger model contains the smaller one as a special
case, so the optimizer can always fall back on the old solution, and it only
moves away if that lowers the residual.

```
  fit error
     |
     |  *
     |      *
     |          *   *   *   *   *     <- always going down,
     |                                   never tells you where to stop
     +--------------------------------> number of parameters k
        3   4   5   6   7   8
```

Here is that bias with real numbers. The data was generated from
`R(10)-(R(200)|C(2e-5))` with 2% noise, so we know the true answer:

```
Circuit                                  k    err%     RSS
R(10)-(R(200)|C(2e-5))                   3    2.379   42.885   <- the truth
R(10)-(R(150)|C(2e-5))-(R(50)|C(1e-6))   5    2.333   41.201   <- fits better
```

**The wrong circuit won.** Not by accident - this is what always happens. The
extra branch spent its freedom on the noise, and noise is real data that can
genuinely be fitted.

So "which circuit fits best?" is the wrong question. The right one is:

> **Is the improvement bigger than what an extra parameter buys by chance?**

---

## 2. Central idea: make parameters cost something

Put the fit and the complexity into one currency and add them up:

```
  score  =  (how badly it fits)  +  (what the parameters cost)
            \________________/      \___________________/
             falls as k grows        rises as k grows
```

Lower score wins. Now there is a minimum somewhere in the middle:

```
  score
     |   \                             /
     |    \                          /        penalty
     |     \                       /
     |      \___               ___/
     |          \___       ___/
     |              \_____/
     |                 ^
     +-----------------|-------------------> k
      too simple:      |      too complex:
      real structure   |      fitting the noise
      is missing       |
                   the model the data supports
```

**That minimum is the whole idea.** Everything below is just what to write in
the two brackets.

---

## 3. The two formulas

For least squares with independent Gaussian residuals:

```
  n   = 2 * (number of frequencies)     Re and Im are separate residuals
  k   = number of free parameters       fixed ones do not count
  RSS = weighted sum of squares         what the optimizer minimizes

  AIC = n * ln(RSS/n) + 2k
  BIC = n * ln(RSS/n) + k * ln(n)
        \___________/   \_______/
          the misfit     the price
```

The first term is the Gaussian log-likelihood times -2, with the noise
variance estimated from the data, and an additive constant (`n + n*ln(2*pi)`)
thrown away. Throwing it away is fine because it is identical for every
candidate - **and it is exactly why a single AIC value means nothing.** Only
differences between models can be read.

The two criteria differ only in the price tag:

```
  cost of one extra parameter, n = 160 (80 frequencies)

    AIC:  2.00   |==|
    BIC:  5.08   |=====|         ln(160) = 5.08

  BIC charges about 2.5x more -> BIC is the stricter one
```

Note `n` counts residuals, not points: each frequency gives a real and an
imaginary residual, so 80 frequencies make n = 160.

---

## 4. Following the arithmetic through

Take the two circuits from section 1. The second RC branch adds 2 parameters:

```
  misfit term, k=3:   160 * ln(42.885/160)  =  -210.663
  misfit term, k=5:   160 * ln(41.201/160)  =  -217.073
                                               ---------
  what those 2 parameters bought:                 6.410
```

Now put it against the price:

```
  AIC charges  2 * 2         =   4.000     bought 6.410 > paid 4.000
                                           -> AIC buys the second branch

  BIC charges  2 * ln(160)   =  10.150     bought 6.410 < paid 10.150
                                           -> BIC refuses it
```

Both are right about their own question. The residual really did improve, and
AIC's threshold really was cleared. BIC just demands more evidence before it
believes another process exists - and here that is what recovers the circuit
the data actually came from.

This is what the toolkit prints for that case:

```
  rank  -c  Circuit                                  k    err%     dAIC     dBIC      cond
     1   1  R(10)-(R(200)|C(2e-5))                   3    2.38      2.4      0.0   1.2e+00
     2   3  R(10)-(R(150)|C(2e-5))-(R(50)|C(1e-6))   5    2.33      0.0      3.7   1.3e+08
     3   2  R(10)-(R(200)|Q(2e-5,0.9))               4    2.38      4.4      5.1   6.0e+01
```

Each criterion column is referenced to its own best value, so the winner shows
0.0 and the rest are distances from it. **`dAIC` and `dBIC` point at different
rows here** - that is not a bug in the table, it is an honest report of a
genuinely ambiguous case.

---

## 5. AIC vs BIC - what each one is actually asking

```
   AIC                                BIC
   ---                                ---
   "Which model will predict          "Which model most probably
    new data best?"                    generated this data?"

   Does NOT assume the truth is       Assumes the truth IS among
   among your candidates - just       your candidates, and converges
   picks the best approximation       on it as data grows

   Charges 2 per parameter            Charges ln(n) per parameter

   Tends to keep slightly too         Tends to be strict; can miss a
   many parameters                    weak but real process
```

**The toolkit selects by BIC.** For EIS that is the better default: candidates
are physically motivated circuits, so assuming the true structure is among
them is reasonable, and inventing a relaxation process that is not there is a
more expensive mistake than a slightly worse prediction.

When they disagree, the disagreement *is* the result. It says the extra
element sits right at the edge of what your data can carry - so settle it with
physics, a repeat measurement, or a wider frequency range. Not by picking
whichever criterion agrees with you.

---

## 6. Reading the table

| Column | What it means |
|--------|---------------|
| `rank` | Position by BIC; rank 1 is the selected circuit |
| `-c` | Position on the command line; figures are saved under this number |
| `k` | Free parameters - what the penalty is charged on |
| `err%` | Relative fit error; **the column you must not select on** |
| `dAIC` | Distance from the best AIC. Predictive view |
| `dBIC` | Distance from the best BIC. Selects the winner |
| `cond` | Condition number of the fit; how identifiable the parameters are |

The difference scale, for both criteria:

```
  dBIC  < 2     the data cannot tell these models apart
        2 - 4   weak preference
        4 - 7   noticeable
        > 10    decisive
```

A difference below 2 does **not** mean the two circuits are equally good
physics. It means this measurement cannot separate them - a statement about
your data, not about electrochemistry. The remedy is a wider frequency range
or less noise, not a different criterion.

### The `cond` column

It is `cond(J^T J)` of the **column-scaled** Jacobian: every column is
normalised to unit norm before the SVD, so parameters spanning many decades
(`R ~ 1e5 Ohm`, `Q ~ 1e-7`) cannot inflate it. The number is therefore
unit-free - a large value is genuine parameter correlation, never an artifact
of the units. Full definition in
[WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md#35-parameter-standard-error-se).

It is a *squared* quantity, and `sqrt(cond)` is the one with the physical
meaning: the factor by which the relative uncertainty of the worst-determined
parameter direction exceeds the best-determined one.

```
  cond(J^T J)     sqrt      reading
  -----------------------------------------------------------------
  1 - 10          1 - 3     directions near-orthogonal, CI trustworthy
  1e2 - 1e4      10 - 100   noticeable correlation, CI still usable
  1e6 - 1e8      1e3 - 1e4  strongly coupled, edge of identifiability
  > 1e10         > 1e5      `!` in the table, covariance unreliable
```

So `1e8` reads worse than it is: the worst direction is 1e4 times less certain
than the best, not 1e8. Compare candidates against each other rather than
against the threshold - three orders between two rows of one table is a signal
long before either row is flagged.

---

## 7. When an extra parameter does earn its place

The penalty is a threshold, not a bias toward small circuits. Same comparison,
but on data that really does contain a CPE
(`R(10)-(R(200)|Q(2e-5,0.85))`, 2% noise):

```
Circuit                        k       RSS      err%       BIC
R(10)-(R(200)|Q(2e-5,0.85))    4     77.58     2.377    -95.52
R(10)-(R(200)|C(2e-5))         3   1300.15     9.541   +350.43
```

`dBIC = 446`. The exponent costs 5.08 and buys 446, so it is kept without
argument. A criterion that rejected this would be useless.

The opposite extreme is just as telling. Fit a CPE to the data from section 1,
which contains no CPE, and the optimizer drives the exponent to exactly
`n0 = 1.0` - the CPE collapses into a plain capacitor and the RSS comes out
*identical* to the 3-parameter fit (42.885 both times):

```
  bought:  0.000
  AIC paid: 2.00   ->  dAIC = 2.0     the parameter is rejected
  BIC paid: 5.08   ->  dBIC = 5.1     purely on its price
```

A parameter that buys nothing costs exactly its penalty. The mechanism working
as designed.

---

## 8. What to watch out for (limits)

- **The criteria do not know physics.** They rank how efficiently the numbers
  are explained. A circuit with no physical meaning can win. The criterion
  narrows the field; it does not do the electrochemistry.

- **They do not check identifiability.** This is the big one. A circuit can be
  overparametrized so badly that its parameters are not separately
  determinable - two merged time constants, or a `Q` whose magnitude and
  exponent trade off perfectly - and still score well. Look back at the table
  in section 4: the two-branch circuit that won on AIC carries `cond = 1.3e+08`
  against `1.2e+00`, eight orders worse. That is the circuit telling you its
  parameters are guesswork.

  Note it is *not* flagged: the `!` marker fires only above `1e10`, where the
  covariance is numerically unusable. So read the column, do not wait for the
  marker.

- **`cond` cannot see a parameter that does nothing.** The column scaling that
  makes it unit-free also normalises away how *strongly* each parameter acts,
  so it reports only whether parameter directions are mutually **redundant**. A
  parameter pinned at a bound - `R1 = 1e10` standing in for "this branch is not
  really there" - has a near-zero gradient, gets scaled back up to unit norm,
  and leaves `cond` looking perfectly healthy:

  ```
    rank  -c  Circuit                    k    err%     dAIC     dBIC      cond
       3   2  L()-R()-(R()|Q()|C())      6    4.90      2.8      6.0   9.8e+00
                                              R1 = 1.00e+10  [at upper bound]
  ```

  That row is not lying, it is answering a different question. Such a model is
  still charged for a parameter that buys nothing, so it trails the same
  circuit without `R1` by roughly the price of one parameter - here `2.8`
  against AIC's `2.00` and `6.0` against BIC's `ln(174) = 5.16`. The `0.8`
  excess is the same in both columns because it comes from the shared misfit
  term: the larger model's optimizer landed on a marginally worse RSS. Read
  `cond` together with the at-the-bound warnings the fit prints: the two
  diagnostics are complementary and neither replaces the other.

- **They assume your weighting describes the noise.** The residuals are treated
  as independent and Gaussian after weighting. `modulus` (the default) means
  noise proportional to `|Z|`, which is usually right for a potentiostat. If
  the real noise is very different, these criteria degrade along with every
  other statistic built on the same assumption.

  The independence half of that assumption is now measured rather than
  assumed. Every fit prints a `Residuals:` line with the lag-1 autocorrelation
  and a runs-test p-value, and warns when the residuals are not random. When
  it warns, `n` overstates how much independent information the spectrum
  carries, and the AIC/BIC arithmetic above is built on sand - a circuit
  missing an element is being compared against another circuit missing an
  element. Fix the model first; the ranking of a shortlist of inadequate
  circuits is not worth reading.

  The obvious repair - substituting an effective sample size
  `n_eff = n(1-rho)/(1+rho)` for `n` - is deliberately *not* applied, though
  `n_eff` is reported. It is computed from the residuals, so it differs per
  candidate and their BICs would no longer sit on a common scale; it does
  nothing when the model is adequate (`rho ~ 0`), which is when model
  selection actually happens; and where it does bite it degenerates. A
  measured `rho1` of `0.987` gives `n_eff = 1.0`, hence `ln(n_eff) = 0` and no
  complexity penalty at all. Correlated residuals in EIS are model error, not
  AR(1) noise, and the formula treats them as noise.

- **They only compare what you propose.** Nothing is invented. A shortlist of
  three bad circuits produces a confident winner that is still a bad circuit.

- **Comparisons are valid only within one run.** Same data, same weighting.
  `n * ln(RSS/n)` is not comparable across weightings because RSS is not, and
  nothing will warn you. This is why the feature is one invocation with
  repeated `--circuit` rather than something you assemble from separate runs.

---

## Summary in one thought

> A lower residual is not evidence of a better circuit, because an extra
> parameter always lowers it. AIC and BIC price that parameter, so the question
> becomes **"did the extra element buy more than it cost?"** - and only the
> differences between candidates ever mean anything.

```
   Candidate circuits (you propose them - physics, DRT, experience)
        |
        |  fit each: same data, same weighting
        v
   RSS and k for every candidate
        |
        |  n*ln(RSS/n) + penalty      AIC pays 2k, BIC pays k*ln(n)
        v
   dAIC / dBIC     <2 indistinguishable ... >10 decisive
        |
        |  + cond: are the winner's parameters even identifiable?
        v
   the circuit your data actually supports
```

Practical order of work: let DRT tell you **how many processes there are**,
write down 2-4 physically defensible circuits, fit them in one run, read
`dBIC`, then check `cond` **and the at-the-bound warnings** on the winner
before believing its parameters.

```bash
eis data.DTA \
  --circuit "R(100) - (R(5000) | C(1e-6))" \
  --circuit "R(100) - (R(5000) | Q(1e-6, 0.9))" \
  --circuit "R(100) - (R(5000) | Q(1e-6, 0.9)) - (R(100) | C(1e-5))"
```

---

## References

- Akaike, H. (1974). *A new look at the statistical model identification.*
  IEEE Trans. Automat. Contr. 19, 716-723. DOI: 10.1109/TAC.1974.1100705
- Schwarz, G. (1978). *Estimating the dimension of a model.*
  Ann. Statist. 6, 461-464. DOI: 10.1214/aos/1176344136
- Burnham, K. P., Anderson, D. R. (2002). *Model Selection and Multimodel
  Inference*, 2nd ed., Springer. Source of the difference scale in section 6
  and of the AIC/BIC comparison in section 5.

---

## Related documentation

- [DRT_INTUITION.md](DRT_INTUITION.md) - how many processes are there? Decide
  this before writing the shortlist
- [GMM_BAYESIAN_INTUITION.md](GMM_BAYESIAN_INTUITION.md) - the same BIC idea
  applied to counting DRT peaks
- [WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md) - what RSS, the
  weights, the fit error and the `cond` column actually are
- [CIRCUIT_PARSER.md](CIRCUIT_PARSER.md) - circuit expression syntax
- [DIFFERENTIAL_EVOLUTION.md](DIFFERENTIAL_EVOLUTION.md) - the optimizer that
  produces each candidate's fit
