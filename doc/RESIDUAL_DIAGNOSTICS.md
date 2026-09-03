# Are the residuals noise? - an intuitive introduction

This document explains **why a small fit error does not mean a correct
circuit**, and what the `Residuals:` line the toolkit prints does about it. The
goal is that after reading it you can look at that line and know whether your
model is finished or merely close.

CLI usage is in [README.md](../README.md). The implementation
is `analyze_residuals` in `eis_analysis/fitting/residual_diagnostics.py`, the
printing in `eis_analysis/cli/handlers/fitting.py`.

---

## 1. The problem: the fit error measures size, not correctness

The fit error is one number: how far the model curve sits from the data on
average. It says nothing about *where* it sits far, and that is the half that
carries the diagnosis.

Two fits with an identical error can be entirely different situations:

```
  residual                          residual
     |   *     *   *                   |        * * *
     |     *  *  *    *  *             |     *         *
   0 |--*----*------*---*---*--      0 |--*-------------*------
     |  *   *    *      *              |*                 *
     |    *    *     *                 |                    * *
     +------------------------->       +------------------------>
              log f                             log f

     scattered: nothing left           smooth: a whole process
     but measurement noise             the model does not have
```

Both can average to 3%. The left one is finished. The right one is a circuit
missing an element, and the missing element is sitting there in plain sight,
drawn in the residual.

Here is that with real numbers. Data generated from two RC branches with
closely spaced time constants, fitted with a single branch:

```
  Fit error: 3.74% (rel)
  Quality: Good (<10.0%)
```

**That model is wrong.** Not marginally - it is missing an entire relaxation.
Adding the second branch takes the error to 1.09%. Nothing in the error line
told you to try.

So "how large is the error?" is an incomplete question. The other half is:

> **Does what is left over still have a shape?**

---

## 2. Central idea: a correct model leaves behind something with no memory

The reasoning is one step long.

Whatever the model does not reproduce ends up in the residual. There are only
two things it can be:

```
  measurement noise      each point drawn independently
                         -> knowing r(f_i) tells you nothing about r(f_i+1)

  missing physics        a relaxation, a diffusion tail, a wrong element
                         -> a smooth function of frequency
                         -> knowing r(f_i) tells you a lot about r(f_i+1)
```

Noise has no memory. Physics does. So the test is not "how big is the
residual" but **"does the residual remember where it has been?"**

```
   independent                       correlated
   -----------                       ----------
   +  -  +  +  -  +  -  -  +  -      +  +  +  +  +  -  -  -  -  -
   sign flips constantly             sign persists in blocks
   neighbours unrelated              neighbours nearly equal
   nothing left to model             a process left to model
```

That is the whole idea. The three statistics below are three ways of asking
that one question, and a fourth number saying what shape the memory has.

---

## 3. The three statistics

All of them are computed on the **weighted** residuals - the ones the
optimizer actually minimized - and separately for the real and imaginary
parts. The two parts are never concatenated: the join between them is not a
frequency step, and they usually carry different structure anyway.

### 3.1 `rho1` - lag-1 autocorrelation

The correlation between each residual and its neighbour:

```
  rho1 = sum( r_i * r_i+1 ) / sum( r_i^2 )      on the centred series

  rho1 ~  0     neighbours unrelated          noise
  rho1 -> 1     neighbours nearly equal       a systematic offset that
                                              drifts slowly with frequency
  rho1 -> -1    neighbours opposite           point-to-point alternation
```

It is the most direct reading of "memory" and needs no threshold to be
useful: `+0.99` and `-0.05` are not close calls. In practice a correct model
lands between `-0.2` and `+0.2`, a model missing a relaxation lands above
`+0.85`.

### 3.2 `runs p` - the Wald-Wolfowitz runs test

`rho1` is a single average and can be fooled. The runs test asks a different
question with a proper p-value attached: how many times does the residual
change sides?

Count maximal blocks on one side of the median. Independent signs give a known
expected number of blocks, and the observed count is compared with it:

```
  n = 80 residuals, balanced

  expected runs if independent :  41.0
  measured on white noise      :  44     z = +0.68   p = 0.50     ok
  measured on a bad fit        :   4     z = -8.33   p = 8e-17    not noise
```

Four runs across eighty points means the residual crossed the middle three
times in seven decades. No sequence of independent draws does that.

Two details that matter:

- **The median, not zero.** A weighted fit does not force its residuals to
  zero mean. Take a bump sitting on an offset large enough that every residual
  comes out positive: a zero-referenced test then has 80 points on one side and
  none on the other, so it cannot count a single crossing and is undefined. The
  median splits the bump from its own baseline and flags it at
  `runs = 18, p = 2.3e-07`. On residuals that are already centred the two
  references agree closely; they part company exactly on that case.
- **The threshold is 0.01, not the customary 0.05.** With 80-160 points the
  test is sensitive enough that 5% flags spectra whose residuals look visibly
  random. This diagnostic is worth acting on only when it is not borderline.

`rho1` and `runs p` together answer **"is this systematic?"** - and they assume
nothing whatever about what shape the systematic part has. That makes them the
detectors of record.

### 3.3 `trend` and `structure` - what shape the memory has

Once the answer is "yes, systematic", the useful follow-up is *what kind*,
because the two kinds ask for different repairs:

```
   a trend                            structure left over
   -------                            -------------------
   the residual marches               the residual wanders up and
   steadily across the window         down within the window

   an element is missing, or          the elements are the right kind,
   is of the wrong type               there are too few of them

   try Q instead of C,                add another parallel branch
   add a Warburg
```

Both are measured, always, and neither is allowed to hide the other:

```
  1. fit a straight line through the residuals    -> the TREND
     slope [per decade], and its p-value

  2. subtract that line

  3. measure what is left with a periodogram      -> the STRUCTURE
     amplitude [same units], and its power
```

Step 2 is the load-bearing one. Without it the drift and the wobble compete
for the same measurement, the drift wins, and any genuine second shape
disappears underneath it - which is exactly how a spectrum with both used to
be reported as a drift alone, prescribing half the repair.

The two are then directly comparable, because they are in the same units:

```
  span      = |slope| * window_decades      the trend's total swing
  amplitude =                               the structure's size

  span >> amplitude    ->  mostly a drift
  span << amplitude    ->  mostly leftover structure
  span ~  amplitude    ->  both are real; do both repairs
```

`span` is printed next to the slope so you never have to do that
multiplication yourself.

---

## 4. Reading the line

Every `--circuit` and `--voigt-chain` fit ends with:

```
  Residuals: rho1 = +0.99 / +0.99 (Re/Im), runs p = 8.3e-17 / 5.4e-16
!   Residuals are not random (Re/Im, 7.0 decade window):
!     trend:     -2.2 / -2.3 per decade, span 15 / 16 (p = 7.6e-04 / 1.6e-03)
!     structure: amplitude 14 / 17 (power 0.73 / 0.89, noise < 0.2)
!     A trend means an element is missing or of the wrong type. Structure left after
!     it means the right elements, too few of them - the residual plot shows where.
```

The first line always prints. The indented block appears **only when the runs
test rejects independence** on at least one part, so silence there is a pass.

| Field | Reading |
|-------|---------|
| `rho1` | Lag-1 autocorrelation. Near 0 is noise; above ~0.85 is a missing process |
| `runs p` | Runs test against the median. Above 0.01 is a pass |
| `trend` | Slope per decade, and `span`, its total swing across the window |
| `(p = ...)` | Significance of the slope. Small means the drift is real |
| `structure` | `amplitude` of what is left once the line is removed |
| `(power ...)` | Normalized periodogram peak in [0, 1]. Below 0.2 is noise |

Everything is printed `Re / Im`, and `n/a` appears where a statistic is not
defined.

**A `Good` quality alongside that warning is not a contradiction.** The fit is
close and the model is still wrong. The two lines answer different questions,
and the residual line is the one about correctness.

### Where the 0.2 comes from

The power threshold is calibrated, not chosen. Over 2000 white-noise series of
80 points the peak power has median 0.10, p95 0.16 and max 0.30, so 0.2 lets
about 0.9% through on its own; every systematic case tested landed at
0.33-0.98. Nothing is suppressed by it - the power is always printed - it is
the level you compare against before acting on the structure line.

---

## 5. What the shapes look like

Six residual series of 80 points over 7 decades, built so the answer is known,
run through the diagnostics:

```
shape              rho1     runs p   runs    span   ampl   power
--------------------------------------------------------------------
white noise       -0.21    5.0e-01     44    0.31   0.37    0.06
ripple, 1.5 dec   +0.86    6.8e-11     12    0.69   3.04    0.93
single bump       +0.75    2.3e-07     18    0.15   1.46    0.73
monotone drift    +0.91    3.4e-15      6    6.77   0.39    0.20
drift + ripple    +0.96    1.1e-13      8    4.49   0.84    0.80
chirp, decaying   +0.93    3.4e-15      6    1.52   1.16    0.73
```

Reading down the table:

- **White noise** is the only row that passes. Note `power = 0.06` and
  `amplitude = 0.37`: the periodogram always finds *a* peak, and it is the
  power that tells you the peak means nothing.
- **A ripple** puts everything on the structure line - `span` stays small
  because a symmetric wave has no net slope.
- **A single bump** does the same, for the same reason. The diagnostics cannot
  separate "one hump" from "several waves", and after the periodogram stopped
  reporting a period they no longer pretend to. The residual plot separates
  them at a glance.
- **A monotone drift** is the mirror image: `span = 6.77` against
  `amplitude = 0.39`, a ratio of 17.
- **Drift plus ripple** carries both, and both are reported. This row is the
  reason the output is two measurements instead of one verdict.
- **A chirp** - a wave whose spacing spreads out and whose amplitude decays -
  is what real fits actually leave behind. It is detected exactly like the
  rest. See section 7.

---

## 6. Two worked repairs

### A. Structure dominant -> add a branch

Data from two RC branches with close time constants, fitted with one branch:

```
  Fit error: 3.74% (rel)
  Residuals: rho1 = +0.90 / +0.87 (Re/Im), runs p = 4.9e-09 / 6.8e-06
!   Residuals are not random (Re/Im, 7.0 decade window):
!     trend:     -0.1 / -0.0081 per decade, span 0.73 / 0.056 (p = 1.0e-01 / 8.9e-01)
!     structure: amplitude 0.92 / 0.88 (power 0.33 / 0.34, noise < 0.2)
```

Read it: the trend's p-values are 0.10 and 0.89 - no drift worth the name. The
structure carries `0.92 / 0.88` at power `0.33 / 0.34`, above the 0.2 floor.
Right elements, too few of them. Add the second branch:

```
  Fit error: 1.09% (rel)
  Residuals: rho1 = -0.05 / -0.14 (Re/Im), runs p = 1.0e+00 / 6.5e-01
```

The warning is gone, `rho1` has collapsed to noise, and the error fell by a
factor of three.

### B. Trend dominant -> an element is missing

Data with a Warburg diffusion tail, fitted without one:

```
  Fit error: 11.35% (rel)
  Residuals: rho1 = +0.92 / +0.93 (Re/Im), runs p = 3.0e-12 / 6.8e-11
!   Residuals are not random (Re/Im, 7.0 decade window):
!     trend:     -0.66 / +1.3 per decade, span 4.6 / 8.9 (p = 1.1e-04 / 1.0e-11)
!     structure: amplitude 3.5 / 3.9 (power 0.75 / 0.89, noise < 0.2)
```

Now the trend is decisive - `p = 1e-04` and `1e-11` - and on the imaginary
part its span (8.9) is more than double the structure's amplitude (3.9). An
element is missing or of the wrong type. Add the Warburg:

```
  Fit error: 1.09% (rel)
  Residuals: rho1 = -0.07 / -0.12 (Re/Im), runs p = 6.5e-01 / 1.0e+00
```

Note that in B the structure line is significant too (power 0.75 / 0.89). That
is not a contradiction and not noise: a missing Warburg deforms the residual in
more than one way. The line to act on first is the one with the larger number,
and in this case one repair cleared both.

---

## 7. Why no period is printed

Earlier versions printed the period of the strongest wave. That number is gone,
and it is worth knowing why, because the intuition it encouraged was wrong.

A period only means something if the residual repeats at a **fixed spacing** in
`log f`. Measured residuals generally do not. On a reference oxide fit the
extrema were read off directly:

```
  part   extrema [Hz]                          spacings [decades]
  ---------------------------------------------------------------
  Re     2.5e-3, 1e-1, 3, 1e3                  1.6 / 1.5 / 2.5
  Im     2.5e-3, 2e-2, 6e-1, 1e2, 2e5          0.9 / 1.5 / 2.2 / 3.3
```

The spacing **grows** towards high frequency - in the imaginary part roughly
geometrically, ratio ~1.5 - while the amplitude decays, from +-6-9% on the
first oscillation to +2.5% on the third. There is no period there to report.

Fitting a single period to that returns an amplitude-weighted average of the
local spacings, which is why the two parts of the same fit came back with 4.0
and 4.7 decades. That difference measured where the amplitude happened to sit,
not anything about the structure - and a reader who took it at face value would
conclude the process spacing differs between Re and Im, which it does not.

`amplitude` and `power` survive the loss of periodicity because neither claims
the shape repeats. Normalized power is the share of the variance the best
single sinusoid explains **at any scale**; the chirp row in section 5 lands at
0.73, nowhere near the noise floor.

One consequence to keep: a hypothesis that requires a constant ratio
`tau_k+1 / tau_k` - a hierarchical ladder of relaxations - requires a constant
period in `log f`. The measured spacings above do not support one.

---

## 8. What to watch out for (limits)

- **A constant offset is invisible here, by design.** Shifting every residual
  by the same amount leaves the crossings where they were, so the runs test
  cannot see it, and `rho1` centres the series first, so neither can it. That
  is correct: a constant offset has no shape, and its size is exactly what
  `fit_error_rel` already reports. These three statistics exist for the
  opposite case, where the error is small but structured.

- **The structure line is the least sensitive of the three.** The share of
  variance one sinusoid can explain falls as the spacing spreads out. On a
  strongly chirped residual the power sags towards the floor while `rho1` and
  the runs test are unaffected. When they disagree, believe the runs test.

- **`power` is not usable on a sparse spectrum.** The 0.2 threshold was
  calibrated at 80 points. At 10 points over 6 decades white noise exceeds it
  in essentially every trial (p95 = 0.82), and at 20 points nearly so. On a
  short spectrum read `rho1` and `runs p`, and ignore the power.

- **`slope_p` is optimistic.** It comes from a least-squares regression that
  assumes independent points - which is the very thing the runs test is there
  to check. When the residuals are correlated it is too small. Use it to rank
  the real part against the imaginary part, not as an exact probability.

- **The diagnostic says the model is wrong, not which element is missing.** It
  points at a class of repair. Which frequency range is unaccounted for is a
  question for the residual panel of the fit plot, and how many processes there
  are is a question for [DRT_INTUITION.md](DRT_INTUITION.md).

- **It cannot tell model error from bad data.** A sample drifting during the
  measurement leaves smooth, correlated residuals that look exactly like a
  missing element, and no amount of extra circuit elements will fix it. Run the
  Kramers-Kronig validation first: if the data itself fails consistency, the
  residual diagnostics of any fit to it are unreadable.

- **The parts are tested separately, and one can fail alone.** That is
  information, not a glitch - a model can be right about the resistive
  behaviour and wrong about the capacitive one. Both parts are always printed
  so the failing one cannot be mistaken for the whole picture.

- **`n_eff` is reported but never used.** The effective sample size
  `n_eff = n(1-rho)/(1+rho)` puts a number on how much independent information
  the spectrum really carries, and it is deliberately **not** substituted into
  AIC/BIC. See
  [MODEL_SELECTION_AIC_BIC.md](MODEL_SELECTION_AIC_BIC.md#8-what-to-watch-out-for-limits)
  for the three reasons.

---

## Summary in one thought

> A small fit error says the model is close. Only the shape of what is left
> says whether it is right. If the residual still remembers where it has been,
> there is physics in it that belongs in the circuit - and **the fit error will
> never tell you that**, because a wrong model can be close everywhere.

```
   Fit converges
        |
        |  weight the residuals as the optimizer did, split Re / Im
        v
   rho1  and  runs p            "is anything systematic?"
        |
        |  runs p > 0.01 -> done, the residual is noise
        v
   runs p < 0.01: split the shape in two
        |
        +---> trend      slope, span, p        an element is missing
        |                                      or is of the wrong type
        |
        +---> structure  amplitude, power      the right elements,
                                               too few of them
        |
        |  compare span with amplitude - same units - and read the
        |  residual plot to see which frequencies are unaccounted for
        v
   repair the circuit, refit, and check this line again
```

Practical order of work: pass Kramers-Kronig first, fit, read `runs p` before
`err%`, repair whichever of `trend` and `structure` is larger, refit, and only
believe an AIC/BIC comparison once this line is quiet for every candidate.

---

## References

- Wald, A., Wolfowitz, J. (1940). *On a test whether two samples are from the
  same population.* Ann. Math. Statist. 11, 147-162.
  DOI: 10.1214/aoms/1177731909 - the runs test of section 3.2
- Lomb, N. R. (1976). *Least-squares frequency analysis of unequally spaced
  data.* Astrophys. Space Sci. 39, 447-462. DOI: 10.1007/BF00648343
- Scargle, J. D. (1982). *Studies in astronomical time series analysis. II.*
  Astrophys. J. 263, 835-853. DOI: 10.1086/160554 - the periodogram used for
  the structure line, chosen because EIS frequencies are unevenly spaced in
  `log f`
- Boukamp, B. A. (1995). *A linear Kronig-Kramers transform test for immittance
  data validation.* J. Electrochem. Soc. 142, 1885-1894. DOI: 10.1149/1.2044210
  - the source of the practice of reading EIS residuals for structure rather
  than size
- Orazem, M. E., Tribollet, B. (2008). *Electrochemical Impedance
  Spectroscopy.* Wiley.

---

## Related documentation

- [MODEL_SELECTION_AIC_BIC.md](MODEL_SELECTION_AIC_BIC.md) - why this line has
  to be quiet before a candidate ranking means anything
- [WEIGHTING_AND_STATISTICS.md](WEIGHTING_AND_STATISTICS.md) - the weights
  these residuals are computed with, and the rest of the quality workflow
- [DRT_INTUITION.md](DRT_INTUITION.md) - how many processes are in the
  spectrum, which is the question a failing structure line raises
- [PYTHON_API.md](PYTHON_API.md#are-the-residuals-noise-analyze_residuals) -
  `analyze_residuals` and the `SeriesDiagnostics` fields from Python
