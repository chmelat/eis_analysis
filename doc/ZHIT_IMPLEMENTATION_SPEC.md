# Z-HIT Validation - Implementation Specification

Status: IMPLEMENTED (v0.10.0+)

## Overview

Z-HIT (Z-Hilbert Impedance Transform) is a non-parametric method for validating Kramers-Kronig compliance of EIS data. It complements the existing Lin-KK method by providing a faster, model-free validation.

## Motivation

| Aspect | Lin-KK | Z-HIT |
|--------|--------|-------|
| Method | Parametric (Voigt chain fitting) | Non-parametric (numerical integration) |
| Speed | Slower (iterative optimization) | Faster (single pass) |
| Model dependency | Requires model selection (M elements) | Model-free |
| Sensitivity | Better for fitting quality | Better for phase-magnitude consistency |

## Implementation Notes

### Design Decision: Numerical Integration vs FFT

The Z-HIT method can be implemented using two approaches:

1. **FFT-based Hilbert transform** (`scipy.signal.hilbert`)
   - Fast but sensitive to edge effects
   - Requires padding to mitigate boundary artifacts

2. **Direct numerical integration in log-omega space** (used here)
   - Ehm et al. (2001) showed that in log-omega space:
     `H[phi] ~ integral of phi * d(ln omega) + correction term`
   - More robust for finite frequency ranges typical in EIS
   - No padding required, simpler implementation

We chose approach (2) because EIS data has limited frequency range where
edge effects from FFT-based Hilbert transform can distort results.

This design was informed by studying the pyimpspec library implementation,
which uses a more sophisticated pipeline (smoothing + spline interpolation +
integration). Our implementation is simpler but achieves comparable results
for typical EIS data.

### Mathematical Background

The Z-HIT formula (first order):

```
ln|Z(omega)| = ln|Z(omega_ref)| + (2/pi) * integral[phi * d(ln omega)]
```

Second order correction (Ehm et al. 2001):

```
ln|Z(omega)| = first_order + gamma * d(phi)/d(ln omega)
```

where gamma = -pi/6 (derived from Taylor expansion of the Hilbert transform
kernel in log-omega space, see Ehm et al. 2001 eq. 15).

## File Structure

```
eis_analysis/
  validation/
    __init__.py          # Exports
    kramers_kronig.py    # Lin-KK implementation
    zhit.py              # Z-HIT implementation
```

## Public API

### `zhit_validation()`

```python
def zhit_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    ref_freq: Optional[float] = None,
    quality_threshold: float = 5.0,
    optimize_offset: bool = False,
    offset_center: float = 1.5,
    offset_width: float = 3.0
) -> ZHITResult:
    """
    Perform Z-HIT validation on EIS data.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    ref_freq : float, optional
        Reference frequency for magnitude anchor [Hz].
        Default: geometric mean of frequency range.
    quality_threshold : float, optional
        Reference threshold for quality metric calculation [%].
        Default: 5.0 (5% mean residual = quality 0)
    optimize_offset : bool, optional
        Use weighted least-squares offset optimization instead of fixed
        reference point. Default: False
    offset_center : float, optional
        Center of weight window on log10(f) scale (default: 1.5 -> ~31.6 Hz)
    offset_width : float, optional
        Width of weight window in decades (default: 3.0)

    Returns
    -------
    ZHITResult
        Dataclass with validation results and visualization
    """
```

### `zhit_reconstruct_magnitude()`

Low-level function for magnitude reconstruction from phase:

```python
def zhit_reconstruct_magnitude(
    frequencies: NDArray[np.float64],
    phi: NDArray[np.float64],
    ln_Z_ref: float,
    ref_idx: int,
    use_second_order: bool = True,
    optimize_offset: bool = False,
    ln_Z_exp: Optional[NDArray[np.float64]] = None,
    offset_center: float = 1.5,
    offset_width: float = 3.0
) -> NDArray[np.float64]:
    """
    Reconstruct impedance magnitude from phase using Z-HIT method.

    Uses numerical integration in log-frequency space according to
    the modified logarithmic Hilbert transform (Ehm et al. 2001).
    """
```

### `ZHITResult` Dataclass

```python
@dataclass
class ZHITResult:
    Z_mag_reconstructed: NDArray[np.float64]  # Reconstructed |Z| [Ohm]
    Z_fit: NDArray[np.complex128]             # Reconstructed Z [Ohm]
    residuals_mag: NDArray[np.float64]        # Magnitude residuals [%]
    residuals_real: NDArray[np.float64]       # Real residuals (fraction)
    residuals_imag: NDArray[np.float64]       # Imag residuals (fraction)
    pseudo_chisqr: float                       # Pseudo chi-squared
    noise_estimate: float                      # Upper bound noise estimate [%]
    quality: float                             # Quality metric (0-1)
    ref_freq: float                            # Reference frequency [Hz]
    figure: Optional[plt.Figure]               # Visualization

    # Properties
    mean_residual_real: float   # Mean |res_real| [%]
    mean_residual_imag: float   # Mean |res_imag| [%]
    mean_residual_mag: float    # Mean |res_mag| [%]
    is_valid: bool              # True if mean_residual_mag < 5%
```

## CLI Integration

### Arguments

```
--no-zhit           Disable Z-HIT validation
--zhit-optimize-offset  Use weighted offset optimization
```

Z-HIT runs by default alongside Lin-KK validation.

### Output Format

```
============================================================
Z-HIT validation
============================================================
Z-HIT: ref_freq=2.81e+01 Hz
  Mean |res_real|: 1.23%
  Mean |res_imag|: 0.89%
  Pseudo chi^2: 2.70e-02
  Estimated noise (upper bound): 1.39%
Data quality is good (residuals < 5%)
```

### Graph Title

```
Z-HIT residuals (χ²=2.70e-02, noise≤1.4%)
```

Note: Noise estimate is labeled as upper bound because Z-HIT residuals
include both noise and integration approximation errors.

## Noise Estimation

The noise estimate uses the same pseudo chi-squared formula as Lin-KK
(Yrjana & Bobacka 2024):

```python
noise_estimate = sqrt(chi2_ps * 5000 / n_points)
```

**Important caveat:** For Z-HIT, this is an upper bound because:
- Lin-KK fits a model, so residuals are mostly noise
- Z-HIT reconstructs from phase, so residuals = noise + integration error

## Comparison with pyimpspec

Our implementation is simpler than pyimpspec's Z-HIT:

| Feature | eis_analysis | pyimpspec |
|---------|--------------|-----------|
| Smoothing | None | Multiple options (lowess, savgol, modsinc, whithend) |
| Interpolation | None | Multiple splines (akima, makima, cubic, pchip) |
| Integration | Trapezoidal on raw data | Spline integration |
| Offset fitting | Weighted least-squares | lmfit minimization |
| Auto-optimization | No | Yes (tests all combinations) |
| Parallel processing | No | Yes |

For typical EIS data with reasonable signal-to-noise ratio, both approaches
yield similar results. pyimpspec's approach may be more robust for noisy data.

## Quality Interpretation

| Mean Residual | Interpretation |
|---------------|----------------|
| < 2%          | Excellent K-K compliance |
| < 5%          | Good K-K compliance |
| >= 5%         | Data may contain artifacts |

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
   (inspiration for implementation approach)
