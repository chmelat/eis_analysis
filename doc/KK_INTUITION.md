# Kramers-Kronig validation - an intuitive introduction

This document explains **what the KK test actually checks, and why it can
check it without knowing anything about your sample**. The goal is that after
reading it you can look at the KK section the toolkit prints, and at the
residual plot next to it, and know whether to trust the spectrum.

Implementation details (algorithm steps, signatures, defaults) are in
[LinKK_analysis.md](LinKK_analysis.md). All numbers below come from this
toolkit: synthetic spectra from `tests/test_kramers_kronig.py` and the measured
spectra in `example/`.

---

## 1. The problem: a bad spectrum looks like a good one

A low-frequency EIS sweep can take an hour. If the sample changes during that
hour - the oxide grows, the electrolyte warms up, a bubble detaches - the last
points were measured on a different sample than the first ones. The spectrum
still comes out as a smooth, plausible curve.

Nothing downstream will notice:

```
  drifting sample  -->  smooth-looking spectrum  -->  circuit fit   -->  R, C, CPE
                                                 -->  DRT           -->  peaks
                                                          ^
                                          both happily produce numbers
```

A circuit fit will find parameters and DRT will find peaks, because both assume
the data are valid. So you need a check that comes **before** any model of the
sample, and that does not depend on one.

---

## 2. Central idea: two views of one response

The real and the imaginary part of an impedance are not two independent
measurements. For any system that is

- **linear** - the response scales with the excitation,
- **causal** - it does not respond before it is excited,
- **stable** - it does not change while you measure it,

they are two views of the same underlying response, tied together by the
Kramers-Kronig relations. Know one of them over all frequencies and the other
one follows.

```
         measured Z'  ----(KK: what Z'' must be)---->  predicted Z''
                                                            |
         measured Z''  <------------ compare --------------+

         they agree     -->  linear, causal, stable
         they disagree  -->  at least one condition was broken
```

That is the whole test. It needs no circuit and no physics of your sample: only
the three conditions, which every valid EIS measurement has to meet anyway.
Drift breaks stability, too large an amplitude breaks linearity, and a bad
cable or instrument artifact usually shows up as a disagreement too.

---

## 3. The trick: fit a model that cannot break the rules

The KK relations are integrals from zero to infinite frequency. You measured
from, say, 10 mHz to 100 kHz. Integrating over a finite window directly gives
the largest errors exactly at the ends of the spectrum, which is where
problems like drift show up.

Lin-KK (Boukamp 1995, Schoenleber et al. 2014) avoids the integral altogether.
It fits the data with a model that **obeys KK by construction**: a chain of RC
(Voigt) elements.

```
  Z_fit = R_s  +  R_1/(1+jw*tau_1)  +  R_2/(1+jw*tau_2)  + ... +  jw*L

  every term is linear, causal and stable  -->  so is their sum
```

If this model can describe your data, your data are KK-consistent. If it
cannot, no model of that kind can, and the data themselves are at fault.

Two details make it practical:

- **The time constants are fixed** on a logarithmic grid spanning the measured
  frequencies. Only the resistances are fitted, so the fit is linear: no start
  values, no local minima, one exact answer.
- **The chain is a ruler, not a claim.** Its elements are not the processes in
  your sample. It is just a flexible shape that can follow any KK-consistent
  spectrum and nothing else.

The toolkit fits **only the real part**, then predicts the imaginary part from
the fitted chain. The imaginary residuals are then the KK comparison from
section 2, done element by element.

On a stationary synthetic ZARC with 0.2% noise, that prediction lands on the
measured imaginary part to 0.20% on average - the noise level.

---

## 4. How many elements? Too few and too many both lie

The chain needs enough elements to follow the spectrum, but not so many that it
starts chasing noise. Here is the same stationary ZARC fitted with a fixed
number of elements M:

```
   M    mu      largest residual    points above 5 %
   3    0.99        42.7 %              44 of 71      too coarse: false alarm
   8    0.87        15.7 %              26 of 71      on perfect data
  12    0.93         3.2 %               0
  16    0.98         0.8 %               0            <- noise level
  30    0.40         2.7 %               0            too fine: negative R_k,
  60    0.00        32.1 %               8 of 71      the fit falls apart
```

- **Too few** elements cannot follow the spectrum, so perfectly good data fail.
- **Too many** elements make the fit unstable: neighbouring elements take huge
  opposite-sign resistances that cancel on the data and blow up between them.

The warning sign of "too many" is **negative resistances**. That is what the mu
metric measures:

```
  mu = 1 - (sum of |negative R_k|) / (sum of positive R_k)

  mu = 1      no negative R_k
  mu < 0.85   negative R_k carry real weight: stop adding elements
```

The toolkit adds elements one by one and stops at the first M where mu drops
below the threshold (0.85). On this ZARC it stops at M = 24, with a largest
residual of 1.4%.

**Where to start counting matters too.** mu does not fall smoothly with M. On a
very coarse grid a relaxation that sits between two grid points gets fitted with
alternating signs, and mu dips below the threshold far too early. On an exact
two-RC spectrum the original algorithm stopped at M = 7 with 3.8% residuals. So
the search starts only where the fit quality (pseudo chi-squared) has levelled
off, and on that spectrum the residuals come out at 0.15%.

Keep in mind that **mu says nothing about data quality**. It is only the stop
rule for M. The residuals judge the data.

---

## 5. Reading the residuals: pattern before size

A residual is measured minus fitted, divided by |Z|, so 1% means 1% of the
impedance at that frequency.

```
  residual [%]                               residual [%]
    |                                          |
  1 |  .   .       .     .                   10|  *
    | . .  .. .  .   . .  . .                  |    *
  0 |--.-.---.--.---.--.---.--.--              5|       *  *
    |  .   .  .    .   .    .                  |             *   *
 -1 |    .        .       .                   0|-------------------*--*--*--*--
    +----------------------------> f           +----------------------------> f
     low                      high               low                      high

     noise: scattered around zero,             violation: smooth, one sign,
     at the noise level everywhere             concentrated in one region
```

**The pattern is the signal.** A smooth, same-sign structure is a KK violation
even when it is small. Scatter around zero is noise even when it is large.

Where the structure sits hints at the cause:

- **Low frequencies** - drift and non-stationarity. The low-frequency points
  take longest to measure, so that is where the sample has had time to
  change. A synthetic ZARC whose resistance grows 10% during the sweep gives
  imaginary residuals of 9.1% at 10 mHz, 5.1% at 100 mHz, 2.4% at 1 Hz and
  0.0% at 10 Hz.
- **A hump in the middle** - something frequency-specific: nonlinearity around
  one time constant, or a change in the middle of the sweep.
  `example/real_gamry_example.DTA` has a hump of this kind: one sign, from
  0.03 to 4 Hz, peaking at 20% near 0.25 Hz, while the real part stays
  within 3%.
- **High frequencies** - cables, the reference electrode, instrument limits.

---

## 6. The verdict: count the points, do not average them

The toolkit turns the residuals into one yes/no answer, `is_valid`. It could
average them, but a violation usually covers only part of the spectrum, and an
average dilutes it:

```
  drifting ZARC (10 %)        mean |res_imag| = 1.08 %    6 of 71 points above 5 %
  real_gamry_example.DTA      mean |res_imag| = 3.82 %   19 of 72 points above 5 %
```

Both means are well under 5%, and both spectra are broken. So the rule counts
points instead:

> A spectrum passes when **at most 5% of its points** have a residual (the
> larger of real and imaginary) **above 5%**.

The 5% line is the one drawn in the residual plot, so what you see there is
what is judged. The allowance of a few points is for the ends of the spectrum,
where Lin-KK residuals grow even on good data. With 72 points, 3 may be above
the line.

Both numbers are empirical and deliberately loose. They mean "clearly broken",
not "as good as the noise allows". Good measured spectra here stay far below
them: `example/EISPOT-test1.DTA` peaks at 1.1%.

The CLI adds a label graded on the mean residual (excellent, good, acceptable,
marginal), but the verdict overrides it. A spectrum that fails is always
"poor", and a spectrum that passes is at worst "marginal".

---

## 7. When good data fail: the ends of the spectrum

The chain's time constants cover only the measured frequencies. A process that
continues beyond the lowest frequency looks, inside the window, like something
the chain cannot build: the rising start of a much slower relaxation, or pure
capacitance. The KK relations hold, but the test cannot show it with the chain
it has.

Here is a KK-perfect ZARC whose relaxation sits at 0.16 Hz, measured only down
to 1 Hz:

```
                                   largest residual   points above 5 %   verdict
  plain Lin-KK                          52.0 %             23 of 51        fail
  + extended tau grid                   25.0 %             20 of 51        fail
  + series capacitance                   3.2 %              0              pass
  + both                                 0.5 %              0              pass
```

Two tools deal with this:

- **Extending the tau grid** (`extend_decades`) lets the chain place elements
  slower than the lowest measured frequency. The CLI searches 0 to 1 decade
  automatically. An extension is accepted only if it does not push mu below
  the stop value, so it cannot be used to sneak in overfitting.
- **A series capacitance** (`--kk-series-c`) gives the chain the one
  KK-consistent term it cannot build from RC elements: a pure capacitor has no
  real part. This is the physics of blocking systems such as two-electrode
  cells, passive films and oxide layers. It also covers Warburg diffusion: an
  exact semi-infinite or open Warburg fails the test without it.

  It is off by default because it also absorbs drift. On a ZARC whose R grows
  during the sweep (0.2 % noise, 10 seeds), a 10 % drift fails on every seed
  without the series C (largest residual 7.9 %) and on none with it (2.2 %);
  a 20 % drift still passes with it (4.0 %). Use it when the low-frequency end
  is capacitive (-Z'' keeps growing as the frequency drops), not when the
  phase returns toward zero there.

The measured ZrO2-on-Zr spectrum `example/EISPOT-M136113-4.DTA` (two-electrode
cell) shows the same thing on real data. Without a series C it passes, but its
imaginary residuals climb steadily toward the lowest frequencies, reaching 4.3%
at 1.6 mHz. With `--kk-series-c` the climb disappears, the largest residual
drops to 0.37%, and the fitted C is 35 uF. The climb was the missing capacitor,
not a violation.

When a spectrum fails with imaginary residuals dominating while the real part
fits well, the CLI suggests `--kk-series-c` itself. At the high-frequency end the series
inductance L, which is always in the model, plays the same role for cable
inductance.

---

## 8. What to watch out for (limits)

- **Passing does not mean clean.** The thresholds catch clear violations only.
  A ZARC drifting 5% during the sweep passes, with 0 or 1 points above 5%.
  Its residuals still show the low-frequency climb: the largest residual is
  about 4.5%, against 0.6-1.4% for the stationary spectrum. **Look at the
  residual plot, not just the verdict.**

- **The noise estimate is an upper bound.** It is computed from the total
  misfit, so a violation inflates it: `real_gamry_example.DTA` reports 4.8%,
  and that number is the violation, not the potentiostat. On good data it
  approaches the real noise (0.34% on `EISPOT-test1.DTA`).

- **Passing says nothing about your circuit.** The test validates the data,
  not a model of them. A KK-valid spectrum can still be fitted with a wrong
  circuit. For choosing between circuits, see
  [MODEL_SELECTION_AIC_BIC.md](MODEL_SELECTION_AIC_BIC.md).

- **Not every violation of linearity is visible.** The test compares the two
  parts of the response at the excitation frequency. A mildly nonlinear system
  can still produce parts that agree. If nonlinearity is a concern, measure at
  two amplitudes and compare the spectra.

- **mu is not a quality score.** It is always near the threshold on normal
  termination, whatever the data. Do not read anything into its value.

- **The ends of the spectrum are the weak spot.** Good data can fail there
  (section 7) and edge points are tolerated by the verdict (section 6). A
  problem confined to the last two or three points is easy to miss or to
  over-read, so check whether it goes away with `--kk-series-c` before
  blaming the sample.

---

## Summary in one thought

> The real and imaginary parts of a valid impedance are two views of one
> response. Lin-KK fits the real part with a chain that obeys the rules by
> construction, predicts the imaginary part from it, and asks **whether the
> measured imaginary part agrees**. Where it does not, the sample or the
> measurement broke linearity, causality or stability.

```
   Measured Z(omega)
        |
        |  fit the real part with RC elements on a fixed tau grid
        |  (add elements until negative R_k appear: mu < 0.85)
        v
   KK-consistent chain
        |
        |  predict the imaginary part, subtract
        v
   residuals:  pattern -> what and where,  count above 5 % -> verdict
        |
        |  fails only at the ends? try --kk-series-c before blaming the sample
        v
   data you can hand to DRT and circuit fitting
```

```bash
eis data.DTA --no-drt --no-fit --no-zhit                  # KK validation only
eis data.DTA --no-drt --no-fit --no-zhit --kk-series-c    # blocking system
```

---

## References

- Boukamp, B. A. (1995). *A Linear Kronig-Kramers Transform Test for
  Immittance Data Validation.* J. Electrochem. Soc. 142, 1885-1894.
- Schoenleber, M., Klotz, D., Ivers-Tiffee, E. (2014). *A Method for Improving
  the Robustness of linear Kramers-Kronig Validity Tests.* Electrochim. Acta
  131, 20-27. Source of the mu metric and the series capacitance.
- Yrjana, V., Bobacka, J. (2024). *Implementing Kramers-Kronig validity
  testing using pyimpspec.* Electrochim. Acta 504, 144951. Source of the noise
  estimate.

---

## Related documentation

- [LinKK_analysis.md](LinKK_analysis.md) - the implementation: algorithm,
  signatures, defaults, CLI output
- [RESIDUAL_DIAGNOSTICS.md](RESIDUAL_DIAGNOSTICS.md) - are the residuals of a
  circuit fit noise? The same "pattern before size" idea, one step later
- [DRT_INTUITION.md](DRT_INTUITION.md) - what to do with a spectrum once it
  has passed
- [MODEL_SELECTION_AIC_BIC.md](MODEL_SELECTION_AIC_BIC.md) - choosing between
  circuits for a valid spectrum
- [WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md) - the weights and
  the pseudo chi-squared
