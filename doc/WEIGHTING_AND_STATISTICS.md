# Weighting and Statistical Metrics in EIS Analysis

This document explains the weighting types used in EIS data fitting and the statistical metrics used to evaluate data and fit quality.

---

## 1. Why Weighting Matters

Electrochemical Impedance Spectroscopy (EIS) measures impedance over a wide frequency range, where |Z| values can vary by several orders of magnitude (typically from a few Ohms at high frequencies to hundreds of kOhms at low frequencies).

**Problem with unweighted fitting:**
Without weighting, the optimization algorithm would minimize absolute deviations, leading to:
- Dominance of high-impedance points (low frequencies)
- Neglect of low-impedance points (high frequencies)
- Poor fit in the high-frequency region

**Solution:**
Weighting normalizes the contribution of each point so that all frequency regions have comparable influence on the final fit.

---

## 2. Weighting Types

### 2.1 Uniform

```
w = 1
```

**Characteristics:**
- All points have equal weight
- Absolute deviations are minimized
- High-impedance points dominate the fit

**When to use:**
- Data with constant absolute measurement error
- Diagnostics - comparison with weighted results
- Suspected issues in low-frequency region

**CLI:** `--weighting uniform`

---

### 2.2 Square Root

```
w = 1 / sqrt(|Z|)
```

**Characteristics:**
- Compromise between uniform and proportional
- Slightly favors high-impedance points
- Smaller penalty for deviations at low |Z|

**When to use:**
- When proportional gives too much weight to high-frequency points
- Data where error grows slower than |Z|

**CLI:** `--weighting sqrt`

---

### 2.3 Modulus - DEFAULT

```
w = 1 / |Z|
```

**Characteristics:**
- Normalizes to relative deviations
- Each point contributes proportionally to its relative error
- Standard for Lin-KK validation (Schonleber 2014)
- Balances influence of points across entire frequency range

**When to use:**
- Most standard EIS measurements
- When constant relative measurement error is expected
- Recommended as default choice

**CLI:** `--weighting modulus` (default)

---

### 2.4 Proportional

```
w = 1 / |Z|^2
```

**Characteristics:**
- Strongly favors low-impedance points
- Used in Boukamp's pseudo chi-squared
- May underweight low-frequency region

**When to use:**
- High-frequency region analysis (inductance, R_inf)
- When low-frequency data contains artifacts
- Historical compatibility with some tools

**CLI:** `--weighting proportional`

---

### 2.5 Comparison and Recommendations

| Weighting | Formula | Preference | Typical Use |
|-----------|---------|------------|-------------|
| uniform | w = 1 | High \|Z\| | Diagnostics, constant abs. error |
| sqrt | w = 1/sqrt(\|Z\|) | Slightly high \|Z\| | Compromise |
| **modulus** | w = 1/\|Z\| | **Balanced** | **Default, most data** |
| proportional | w = 1/\|Z\|^2 | Low \|Z\| | High-frequency analysis |

**Recommended workflow:**
1. Start with `modulus` (default)
2. If fit is poor at high frequencies, try `proportional`
3. If fit is poor at low frequencies, try `uniform` or `sqrt`
4. Compare results and choose best compromise

---

## 3. Statistical Metrics

### 3.1 Pseudo Chi-Squared (chi^2_ps)

**Definition (Boukamp 1995):**

```
chi^2_ps = SUM[ w_i * ((Z'_exp - Z'_fit)^2 + (Z''_exp - Z''_fit)^2) ]

where:
  w_i = 1 / (Z'_i^2 + Z''_i^2)  ... proportional weighting (1/|Z|^2)
  Z' = real part of impedance
  Z'' = imaginary part of impedance
```

**Interpretation:**
- Dimensionless fit quality indicator
- Lower value = better fit
- Used in KK validation
- Basis for noise estimation

**Typical values:**
- chi^2_ps < 0.01: Very good fit
- chi^2_ps ~ 0.01-0.1: Good fit
- chi^2_ps > 0.1: Fit may be problematic

**Reference:** Boukamp, B.A. (1995). J. Electrochem. Soc., 142, 1885-1894.

---

### 3.2 Estimated Noise

**Definition (Yrjana & Bobacka 2024):**

```
noise [%] = sqrt(chi^2_ps * 5000 / N)

where:
  chi^2_ps = pseudo chi-squared
  N = number of data points
  5000 = 100^2 / 2 (theoretical constant, see derivation below)
```

**Derivation:**

Assuming residuals are caused by noise with standard deviation sigma (as fraction of |Z|):
```
Z'_exp - Z'_fit ~ sigma * |Z|
Z''_exp - Z''_fit ~ sigma * |Z|
```

Substituting into chi^2_ps with w_i = 1/|Z|^2:
```
chi^2_ps = SUM[ 1/|Z|^2 * (sigma^2*|Z|^2 + sigma^2*|Z|^2) ]
         = SUM[ 2 * sigma^2 ]
         = 2 * N * sigma^2
```

Solving for sigma and converting to percentage:
```
sigma = sqrt(chi^2_ps / (2*N))
noise [%] = 100 * sigma = sqrt(chi^2_ps * 100^2 / (2*N))
          = sqrt(chi^2_ps * 5000 / N)
```

The constant 5000 = 100^2 / 2 comes from:
- 100^2: conversion from fraction to percentage
- 2: factor for real + imaginary components

**Interpretation:**
- Estimates noise standard deviation as percentage of |Z|
- Assumes deviations from KK fit are caused by random noise

**Typical values:**
- noise < 1%: High quality data
- noise 1-3%: Good quality data
- noise 3-5%: Acceptable quality
- noise > 5%: Data may contain artifacts

**Reference:** Yrjana, V. and Bobacka, J. (2024). Electrochim. Acta, 504, 144951.

---

### 3.3 Residuals (Mean |res|)

**Definition:**

```
res_real = (Z'_exp - Z'_fit) / |Z_exp|
res_imag = (Z''_exp - Z''_fit) / |Z_exp|

Mean |res_real| [%] = 100 * mean(|res_real|)
Mean |res_imag| [%] = 100 * mean(|res_imag|)
```

**Interpretation:**
- Average relative deviation of real/imaginary parts
- Normalized by |Z| (not by individual components)
- Used in KK validation to assess data quality

**Quality thresholds (KK validation):**
- Mean |res| < 5%: Data satisfies KK relations (valid)
- Mean |res| 5-10%: Minor deviations, data acceptable
- Mean |res| > 10%: Significant deviations, data may be invalid

**Note:** In KK validation, real and imaginary parts are expected to be related by Kramers-Kronig relations. High residuals may indicate:
- System non-stationarity during measurement
- Nonlinear response
- Experimental artifacts (noise, drift)

---

### 3.4 Fit Error

#### Relative Error [%]

```
rel_error_i = |Z_i - Z_fit,i| / |Z_i|
fit_error_rel = [SUM(w_i * rel_error_i) / SUM(w_i)] * 100

where w_i = weights according to selected weighting type
```

**Interpretation:**
- Weighted average of relative deviations
- Depends on selected weighting
- Primary indicator of circuit fit quality

#### Absolute Error [Ohm]

```
fit_error_abs = mean(|Z - Z_fit|)
```

**Interpretation:**
- Average absolute deviation in Ohms
- Independent of weighting
- Useful for comparing fits across different datasets

#### Quality Assessment

| Relative Error | Assessment | Meaning |
|----------------|------------|---------|
| < 1% | Excellent | Outstanding fit |
| 1-10% | Good | Good fit |
| 10-20% | Acceptable | Acceptable fit |
| > 20% | Poor | Poor fit, consider different model |

---

### 3.5 Parameter Standard Error (SE)

**Definition:**

```
SE(theta_i) = sqrt(diag_i(Cov(theta)))

where:
  Cov(theta) = s^2 * (J^T * J)^(-1)
  s^2 = RSS / dof ... residual variance estimate
  RSS = SUM(r_i^2) ... residual sum of squares
  dof = N - p ... degrees of freedom (data points - parameters)
  J = Jacobian matrix (matrix of derivatives)
```

**Interpretation:**
- Estimate of fitted parameter uncertainty
- Based on linear approximation around optimum
- Assumes normally distributed residuals

**Important:**
- SE is valid only if fit is "good" (random residuals)
- High SE may indicate:
  - Parameter not well determined by data
  - Model is overparameterized
  - Correlation between parameters

**Condition number:**
```
cond = S_max / S_min  (via SVD decomposition of Jacobian)
```
- cond < 10^10: Well-conditioned problem, SE are reliable
- cond > 10^10: Ill-conditioned problem, SE may be unreliable

---

### 3.6 Confidence Interval (95% CI)

**Definition:**

```
CI_low = theta - t_(alpha/2, dof) * SE(theta)
CI_high = theta + t_(alpha/2, dof) * SE(theta)

where:
  theta = fitted parameter
  SE(theta) = standard error
  t_(alpha/2, dof) = critical value of t-distribution
  alpha = 1 - confidence level (for 95% CI: alpha = 0.05)
  dof = N - p ... degrees of freedom
```

**Why t-distribution (not normal)?**
- For small datasets, t-distribution is more conservative
- Accounts for uncertainty in variance estimate
- For large N, t approaches normal distribution (t -> 1.96 for 95%)

**Interpretation of 95% CI:**
- With 95% probability, true parameter value lies within this interval
- Narrow interval = well-determined parameter
- Wide interval = high uncertainty

**Example output:**
```
R1 = 1.03e+05 +/- 1.82e+01  [95% CI: 1.03e+05, 1.03e+05]
```
- Value: 103 kOhm
- Standard error: 18.2 Ohm (0.018%)
- 95% CI: very narrow, parameter is well determined

**Warning:**
- CI assumes correct model
- If model is wrong, CI has no statistical meaning
- Always check residuals and fit error

---

## 4. Summary Table of Metrics

| Metric | Range | Unit | Good Values | Use |
|--------|-------|------|-------------|-----|
| chi^2_ps | [0, inf) | - | < 0.1 | KK validation |
| Noise est. | [0, inf) | % | < 3% | Data quality |
| Mean \|res\| | [0, 100+] | % | < 5% | KK validation |
| Fit error rel. | [0, 100+] | % | < 10% | Fit quality |
| Fit error abs. | [0, inf) | Ohm | contextual | Fit quality |
| SE | [0, inf) | [param] | < 10% of param | Parameter uncertainty |
| 95% CI | - | [param] | narrow | Parameter uncertainty |
| Cond. number | [1, inf) | - | < 10^10 | Fit stability |

---

## 5. Practical Recommendations

### Quality Assessment Workflow

1. **KK validation**
   - Check Mean |res| < 5%
   - Check estimated noise < 3%
   - If not satisfied, data may contain artifacts

2. **Circuit fitting**
   - Use `--weighting modulus` (default)
   - Check fit error < 10%
   - Check residuals (should not show systematic trends)

3. **Parameters**
   - Check 95% CI - are they reasonable?
   - Check condition number < 10^10
   - If SE is high, parameter is not well determined

### Troubleshooting

| Problem | Possible Cause | Solution |
|---------|----------------|----------|
| High KK residuals | Non-stationarity, nonlinearity | Shorten measurement, reduce amplitude |
| High fit error | Wrong circuit model | Try different circuit |
| High parameter SE | Overparameterization | Simplify circuit |
| High cond. number | Parameter correlation | Fix some parameters |
| Systematic residuals | Missing circuit element | Add element (CPE, Warburg) |

---

## 6. References

1. **Boukamp, B.A.** (1995). "A Linear Kronig-Kramers Transform Test for Immittance Data Validation." *J. Electrochem. Soc.*, 142, 1885-1894.

2. **Schonleber, M. et al.** (2014). "A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests." *Electrochimica Acta*, 131, 20-27.

3. **Yrjana, V. and Bobacka, J.** (2024). "Implementing Kramers-Kronig validity testing using pyimpspec." *Electrochim. Acta*, 504, 144951.

4. **Orazem, M.E. and Tribollet, B.** (2008). "Electrochemical Impedance Spectroscopy." Wiley.

5. **Bard, Y.** (1974). "Nonlinear Parameter Estimation." Academic Press.
