# Validation Method Comparison: Lin-KK vs Z-HIT

**Project**: eis_analysis
**Date**: 2026-01-09

---

## 1. Introduction

EIS Analysis Toolkit offers two data validation methods:

- **Lin-KK** (Linear Kramers-Kronig) - model-based method
- **Z-HIT** (Z-Hilbert Impedance Transform) - model-free method

Both methods verify whether data satisfy causality and linearity conditions, but use different approaches. This report explains the differences and provides recommendations on when to use each method.

---

## 2. Theoretical Background

### 2.1 Lin-KK

Lin-KK approximates the impedance spectrum with a Voigt chain:

```
Z(w) = R_s + jw*L + SUM_{k=1}^{M} R_k / (1 + jw*tau_k)
```

Where:
- R_s: series resistance
- L: series inductance
- R_k: resistance of k-th RC element
- tau_k: time constant (fixed, logarithmically distributed)
- M: number of elements (automatically optimized)

**Principle:** If data satisfy Kramers-Kronig relations, the Voigt chain can approximate them well. Residuals indicate data quality.

### 2.2 Z-HIT

Z-HIT reconstructs impedance modulus from phase using Hilbert transform:

```
ln|Z(w)| = ln|Z(w_ref)| + (2/pi) * integral[phi(x) d(ln x)] - (pi/6) * d(phi)/d(ln w)
```

Where:
- phi: phase angle
- w_ref: reference frequency
- last term: second-order correction (Ehm 2001)

**Principle:** For a causal system, phase and modulus are related by Hilbert transform. If reconstructed modulus matches measured, data are valid.

---

## 3. Key Differences

| Aspect | Lin-KK | Z-HIT |
|--------|--------|-------|
| Approach | Model-based (Voigt chain) | Model-free (Hilbert transform) |
| Fits | Impedance directly | Reconstructs modulus from phase |
| Assumes | RC/RQ behavior | Only causality |
| Number of parameters | M (optimized) | None |
| Noise sensitivity | Medium | Higher (phase derivative) |

---

## 4. Voigt Chain Limitations

### 4.1 RC Element Behavior

Each RC element in the Voigt chain has impedance:

```
Z_RC = R / (1 + jw*R*C)
```

**Frequency dependence:**
- At low frequencies (w -> 0): Z_RC -> R (purely real)
- At high frequencies (w -> inf): Z_RC -> 0

**Phase dependence:**
- At low frequencies: phase -> 0°
- At high frequencies: phase -> -90°

### 4.2 Problem with Capacitive Data

The Voigt chain at **low frequencies** always becomes **resistive** (phase approaching 0°).

For data with **nearly pure capacitive behavior** (phase close to -90° across entire range), the Voigt chain cannot properly approximate the impedance because:

1. Model requires transition to resistive behavior at low frequencies
2. Data show capacitive behavior even at low frequencies
3. Large residuals arise in imaginary part

### 4.3 Problem with Time Constants

The Voigt chain uses time constants tau_k logarithmically distributed in range:

```
tau_min = 1 / (2*pi*f_max)
tau_max = 1 / (2*pi*f_min)
```

For capacitive data (e.g., oxide layer), the actual system time constant is:

```
tau = R * C
```

If R is very small (low oxide layer resistance), the time constant lies **far outside** the measured frequency range. The Voigt chain then has no element that could represent this behavior.

**Example:**
- Measured range: 1 Hz - 400 kHz
- tau range: 4e-7 s to 0.16 s
- Actual system tau: << 4e-7 s (outside range)

### 4.4 Problem with Fitting Real Part

Lin-KK by default fits the **real part** of impedance and evaluates KK compatibility from imaginary part residuals.

For capacitive data:
- |Z''| >> |Z'| (imaginary part dominates)
- Typically |Z''|/|Z'| > 10

**Problems:**
1. We fit a **small signal** (Z') with relatively higher noise
2. **Large signal** (Z'') is not directly constrained by the fit
3. Small errors in Z' fit manifest as large relative errors in Z''
4. Imaginary part residuals are normalized by |Z|, which is nearly equal to |Z''|

**Illustration:**
```
Data:  Z' = 1000 Ohm,  Z'' = -50000 Ohm,  |Z| = 50010 Ohm

Fit has 5% error in Z':
  Z'_fit = 1050 Ohm  (error 50 Ohm)

But Z''_fit can be completely off:
  Z''_fit = -40000 Ohm  (error 10000 Ohm!)

Imaginary part residual:
  res_imag = (Z'' - Z''_fit) / |Z| = 10000 / 50010 = 20%
```

---

## 5. Practical Demonstration

### 5.1 Test Data

Data from oxide layer measurement (ZRY-3d-1.dta):
- Frequency range: 1 Hz - 400 kHz
- Phase range: -87.7° to -79.6° (close to -90° everywhere)
- |Z''|/|Z| ~ 99.7% (imaginary part dominates)

### 5.2 Validation Results

| Metric | Lin-KK | Z-HIT |
|--------|--------|-------|
| Mean \|res_real\| | 0.07% | 0.05% |
| Mean \|res_imag\| | **55.59%** | 0.67% |
| Pseudo chi^2 | 18.5 | 0.004 |
| Noise estimate | 40.68% | 0.58% |

### 5.3 Interpretation

**Lin-KK:**
- Real part fits well (0.07%)
- Imaginary part has huge error (55.59%)
- Model cannot represent the data
- **False positive** indication of data problems

**Z-HIT:**
- Both components have low residuals (< 1%)
- Data are valid
- **Correct** result

### 5.4 Evidence of Structural Model Limitation

The imaginary residuals are **not randomly scattered** but form a smooth systematic curve. This proves the issue is model inadequacy, not noise:

| Frequency range | Z''_fit / Z''_data | Fit phase | Data phase |
|-----------------|-------------------|-----------|------------|
| Low (1-10 Hz) | 0.23 (23%) | -41° | -80° |
| Mid (100-1000 Hz) | 0.44 (44%) | ~-60° | -85° |
| High (>100 kHz) | 0.85 (85%) | -85° | -88° |

**Key observation:** At low frequencies, the fit phase is **-41°** while data phase is **-80°**. The Voigt chain cannot maintain capacitive behavior (phase near -90°) at low frequencies - it inevitably becomes resistive (phase toward 0°).

The smooth residual curve from -76% (low freq) to -15% (high freq) is a signature of **structural model mismatch**, not measurement noise or fitting artifacts

---

## 6. Recommendations

### 6.1 When to Use Each Method

| Data Type | Recommended Method | Reason |
|-----------|-------------------|--------|
| Typical EIS data (RC/RQ) | Lin-KK or Z-HIT | Both work |
| Capacitive data (phase ~ -90°) | **Z-HIT** | Lin-KK fails |
| Inductive data (phase ~ +90°) | **Z-HIT** | Lin-KK fails |
| Noisy data | **Lin-KK** | Z-HIT sensitive to noise |
| Quick validation | Lin-KK | Faster computation |

### 6.2 Practical Rule

**Check the phase range of data:**

```
Phase range > 30°  ->  Both Lin-KK and Z-HIT work
Phase range < 30° and phase close to +/-90°  ->  Prefer Z-HIT
```

### 6.3 CLI Usage

```bash
# Standard validation (Lin-KK + Z-HIT)
eis data.dta

# Z-HIT only (for capacitive data)
eis data.dta --zhit-only

# Lin-KK only
eis data.dta --no-zhit
```

---

## 7. Summary

### 7.1 Why Lin-KK Fails for Capacitive Data

**Primary cause: Structural model limitation**

The Voigt chain becomes resistive at low frequencies (phase -> 0°), but capacitive data maintain phase close to -90° across entire range. This is proven by:
- Smooth systematic residual curve (not random scatter)
- Fit phase reaching -41° at low frequencies while data stays at -80°
- Z''_fit systematically underestimating Z''_data (only 23% at low freq)

**Contributing factors:**

1. **Time constants outside range** - Actual system tau may lie far outside the range of tau_k used in Voigt chain.

2. **Fitting small component** - Real part Z' is fitted, which for capacitive data is smaller than Z''. However, the smooth residual pattern indicates this is secondary to the structural limitation.

### 7.2 Method Comparison

| Method | Advantages | Disadvantages |
|--------|------------|---------------|
| **Lin-KK** | Robust to noise, good for typical EIS data | Fails for purely capacitive/inductive data |
| **Z-HIT** | Model-free, works for any causal system | More sensitive to noise, requires smooth phase data |

**Recommended workflow:**
1. Run both methods (default behavior)
2. If results differ significantly, check phase range of data
3. For data with phase close to +/-90° across entire range, prefer Z-HIT

---

## References

1. Boukamp, B.A. (1995). A Linear Kronig-Kramers Transform Test for Immittance Data Validation. J. Electrochem. Soc., 142, 1885-1894.

2. Schonleber, M., Klotz, D., and Ivers-Tiffee, E. (2014). A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests. Electrochim. Acta, 131, 20-27.

3. Ehm, W., Kaus, R., Schiller, C.A., Strunz, W. (2001). The evaluation of electrochemical impedance spectra using a modified logarithmic Hilbert transform. Journal of Electroanalytical Chemistry, 499, 216-225.
