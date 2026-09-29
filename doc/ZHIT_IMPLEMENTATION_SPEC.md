# Z-HIT Validation - Implementation Specification

Status: implemented since v0.10.0; this document describes v0.44.1.

Code: `eis_analysis/validation/zhit.py`. API reference: docstrings and
[PYTHON_API.md](PYTHON_API.md). Known open issues: [ZHIT_REVIEW.md](ZHIT_REVIEW.md).

## Overview

Z-HIT (Z-Hilbert Impedance Transform) is a non-parametric check of
Kramers-Kronig compliance. It reconstructs |Z| from the measured phase and
compares the reconstruction with the measured |Z|. It runs by default
alongside Lin-KK.

| Aspect | Lin-KK | Z-HIT |
|--------|--------|-------|
| Method | Parametric (Voigt chain fitting) | Non-parametric (numerical integration) |
| Speed | Slower (iterative choice of M) | Single pass |
| Model dependency | Requires model selection (M elements) | Model-free |
| What it compares | Full complex residual | Magnitude only; the phase is taken as given |

## Mathematical Background

In x = ln(omega), the magnitude of a minimum-phase impedance follows from its
phase through the Kramers-Kronig relation. Ehm et al. (2001) replace the
non-local integral by a local expansion. Written out, the series is (derivation
in [ZHIT_REVIEW.md](ZHIT_REVIEW.md), 1.3.1):

```
ln|Z(x)| = C + (2/pi) * integral[phi dx]
             - (pi/6)       * phi'(x)
             - (pi^3/360)   * phi'''(x)
             - ...
```

The implementation keeps the first two terms (gamma = -pi/6):

```
ln|Z(x)| = C + (2/pi) * integral[phi dx] + gamma * d(phi)/dx
```

### Integration constant

The phase determines ln|Z| only up to C. C is set to
`median(ln|Z_exp| - reconstruction)` over the whole spectrum (since v0.44.0).
Matching a single reference point, as older versions did, carried that point's
noise and the local approximation error into every point of the
reconstruction. The median is robust to outliers and to drift in fewer than
half of the points; it fails when most of the spectrum is bad.

### Accuracy floor

The series is asymptotic, not convergent, so the two-term formula has an
error even on exact, noise-free, K-K compliant data. It peaks where the phase
bends most, i.e. around a relaxation, not at the edges of the frequency range:

| Data (noise-free) | mean | max |
|---|---|---|
| R+RC (ideal, Debye) | 0.7 % | 3.2 % |
| R+CPE, n = 0.8 | 0.3 % | 1.2 % |
| R+CPE, n = 0.6 | 0.1 % | 0.4 % |

Discretization adds almost nothing at 10 points/decade. Adding the phi'''
term would halve the floor on clean data but amplifies noise (third
derivative), making real data worse. Details and measurements:
[ZHIT_REVIEW.md](ZHIT_REVIEW.md), 1.3.

### Why not an FFT Hilbert transform

The exact relation is a multiplication by `-i * coth(pi*k/2)` in the Fourier
domain of ln(omega). Evaluating it (or `scipy.signal.hilbert`) needs the phase
outside the measured range, i.e. padding or extrapolation, and a finite EIS
window makes the result edge-sensitive. The local expansion needs neither.

## Algorithm (`zhit_validation`)

1. Sort by ascending frequency; outputs are returned in the caller's order.
2. `phi = np.unwrap(arctan2(Z.imag, Z.real))` removes 2*pi jumps that would
   spike the derivative.
3. `zhit_reconstruct_magnitude(frequencies, phi, ln|Z|)`:
   cumulative trapezoid of `(2/pi) * phi` over ln(omega), plus
   `-pi/6 * np.gradient(phi, ln omega)`, plus the median offset.
4. `Z_fit = |Z_recon| * exp(j * phi)`: only the magnitude is reconstructed,
   the measured phase is kept.
5. Residuals: `residuals_mag = (|Z| - |Z_recon|) / |Z| * 100` [%];
   `residuals_real/imag = (Z - Z_fit) / |Z|` [fraction].
6. Pseudo chi-squared, noise estimate, quality metric, figure.

Because `Z_fit` carries the measured phase, `Z - Z_fit = (|Z| - |Z_recon|) *
exp(j*phi)`: the real and imaginary residuals are the magnitude residual
projected by cos(phi) and sin(phi), not an independent check
([ZHIT_REVIEW.md](ZHIT_REVIEW.md), 2.1).

## Result

`ZHITResult` fields and properties are documented in its docstring and in
PYTHON_API.md. The decisions built on them:

- `success`: False when the reconstruction raised; all arrays are then empty.
- `is_valid`: `mean_residual_mag < quality_threshold` (default 5 %).
- `quality`: `max(0, 1 - mean_residual_mag / quality_threshold)`.
- `quality_label`: shared with Lin-KK (`_quality_label`):

| Mean \|res_mag\| | Label |
|---|---|
| < 0.5 % | excellent |
| < 1.0 % | good |
| < 2.5 % | acceptable |
| < 5.0 % | marginal (check for drift/nonlinearity) |
| >= 5.0 % | poor |

These thresholds were set for Lin-KK. Given the accuracy floor above, an
ideal RC scores "good" at best under Z-HIT.

## Noise Estimation

The noise estimate uses the Lin-KK formula (Yrjana & Bobacka 2024):

```python
noise_estimate = sqrt(chi2_ps * 5000 / n_points)
```

It is reported as an upper bound: Z-HIT residuals contain the accuracy floor
as well as noise, and, since the real/imag residuals are projections of the
magnitude residual, chi2_ps measures the magnitude deviation only.

## CLI Integration

```
--no-zhit                    Skip Z-HIT validation
--fit-on {original,zhit,all} Fit the circuit (zhit) or every stage (all)
                             against the Z-HIT reconstruction
```

`--fit-on zhit|all` cannot be combined with `--no-zhit`.

Output (`cli/handlers/validation.py`, `run_zhit_validation`):

```
============================================================
Z-HIT validation
============================================================
Z-HIT: second order, offset = median over the spectrum
  Mean |res_real|: 1.20%
  Mean |res_imag|: 0.88%
  Pseudo chi^2: 2.42e-02
  Estimated noise (upper bound): 1.31%
Data quality: acceptable (mean |res_mag|=1.55%, threshold=5.0%)
```

The data-quality line is a warning when `is_valid` is False.

The figure has two panels: measured vs reconstructed |Z| (log-log), and the
real/imag residuals in % with fixed +-5 % guide lines. It is saved as `zhit`
with `--save`.

The per-point outlier report (`validation/outliers.py`, `find_outliers`)
reads `|residuals_mag|` from this result.

## Second use: reconstruction as a data correction (v0.34.0)

`ZHITResult.Z_fit` is not only the reference curve the residuals are measured
against - `--fit-on` makes it the data the circuit is fitted to. See README,
"Z-HIT as a correction", and [ZAHNER_ANALYSIS_REVIEW.md](ZAHNER_ANALYSIS_REVIEW.md)
section 7.

Two constraints follow from the transform and are enforced by where the CLI
calls it: the reconstruction is computed on the full spectrum (a range
truncated by `--f-min`/`--f-max` is a different reconstruction, so
`apply_zhit_reconstruction()` runs before `filter_by_frequency()`), and it
trusts the phase to replace the magnitude, which is what limits it to drift.

Two costs follow from the accuracy floor and the derivative:

- Around sharp relaxations |Z| is shifted by up to ~3 % even on clean data.
  Under `--fit-on all`, R_inf and the DRT read that error wherever a
  relaxation sits; the run warns about it.
- The gamma * phi' term differentiates phase noise. On 1 % noise the
  reconstruction does not make fitted resistances clearly better and can make
  them up to ~3x worse (`tests/test_zhit_fit_on.py`). Use it for drift, not
  for scatter.

## Comparison with pyimpspec

| Feature | eis_analysis | pyimpspec |
|---------|--------------|-----------|
| Smoothing | None | lowess, savgol, modsinc, whithend |
| Interpolation | None | akima, makima, cubic, pchip splines |
| Integration | Trapezoidal on raw data | Spline integration |
| Offset | Median over the spectrum | lmfit minimization over a window |
| Auto-optimization | No | Yes (tests all combinations) |

The missing smoothing is why phase noise reaches the reconstruction
unfiltered (open point 2 of [ZHIT_AUDIT_2026-04-26.md](archive/ZHIT_AUDIT_2026-04-26.md)).
A code-level comparison is in [ZHIT_comparison_report.md](ZHIT_comparison_report.md);
the two have not been benchmarked against each other on the same data.

## References

1. Ehm, W., Kaus, R., Schiller, C.A., Strunz, W. (2001). "The evaluation of
   electrochemical impedance spectra using a modified logarithmic Hilbert
   transform." *Journal of Electroanalytical Chemistry* 499, 216-225.

2. Schiller, C.A., Richter, F., Gulzow, E., Wagner, N. (2001). "Relaxation
   impedance as a model for the deactivation mechanism of fuel cells due to
   carbon monoxide poisoning." *Physical Chemistry Chemical Physics* 3, 374-378.

3. Yrjana, V. and Bobacka, J. (2024). "Implementing Kramers-Kronig validity
   testing using pyimpspec." *Electrochim. Acta* 504, 144951.

4. pyimpspec library: https://github.com/vyrjana/pyimpspec
