# Oxide Layer Analysis Guide

## Overview

Two complementary functions for oxide layer analysis:

| Function | Input | Output | Use case |
|----------|-------|--------|----------|
| `analyze_oxide_layer()` | epsilon_r | thickness | Known permittivity, estimate thickness |
| `estimate_permittivity()` | thickness | epsilon_r | Known thickness (SEM/TEM), estimate permittivity |

**Key features:**
- Finds dominant Voigt (R||C), K, or R||Q element (largest R = main barrier)
- Uses parallel plate capacitor model
- For Q: uses Hsu-Mansfeld formula for effective capacitance (primary,
  3D model); the Brug (2D) estimate and thickness are reported alongside
  for comparison when the circuit contains a series resistance
- Returns structured `OxideAnalysisResult` dataclass

---

## Quick Start

### Python API

```python
from eis_analysis import load_data, fit_equivalent_circuit, R, C
from eis_analysis.analysis import analyze_oxide_layer, estimate_permittivity

# Load and fit
freq, Z = load_data('sample.DTA')
circuit = R(100) - (R(5000) | C(1e-6))
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Option 1: Known permittivity -> estimate thickness
oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
print(f"Thickness: {oxide.thickness_nm:.1f} nm")

# Option 2: Known thickness (from SEM) -> estimate permittivity
eps_r = estimate_permittivity(freq, Z, thickness_nm=20.0, fit_result=result)
print(f"Permittivity: {eps_r:.1f}")
```

### CLI

```bash
# Basic oxide analysis (requires fitted circuit)
python3 eis.py data.DTA --circuit 'R(100) - (R(5000) | C(1e-6))' --analyze-oxide

# With custom permittivity and area
python3 eis.py data.DTA --circuit 'R(100) - (R(5000) | C(1e-6))' \
    --analyze-oxide --epsilon-r 23 --area 0.5

# Using Voigt chain for automatic circuit
python3 eis.py data.DTA --voigt-chain --analyze-oxide
```

**Note:** For oxide analysis with custom data, you must specify `--circuit` or `--voigt-chain`.
Without these flags, fit is skipped and oxide analysis uses less accurate fallback method.

---

## API

### analyze_oxide_layer()

```python
def analyze_oxide_layer(
    frequencies: NDArray,
    Z: NDArray,
    epsilon_r: float = 22.0,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[OxideAnalysisResult]
```

**Parameters:**
- `frequencies` - Measurement frequencies [Hz]
- `Z` - Complex impedance [Ohm]
- `epsilon_r` - Relative permittivity (default: 22 for ZrO2)
- `area_cm2` - Electrode area [cm^2] (default: 1.0)
- `fit_result` - Result from `fit_equivalent_circuit()` (recommended)

**Returns:** `OxideAnalysisResult` or `None` if analysis fails.

### estimate_permittivity()

Inverse function - estimates permittivity from known thickness.

```python
def estimate_permittivity(
    frequencies: NDArray,
    Z: NDArray,
    thickness_nm: float,
    area_cm2: float = 1.0,
    fit_result: Optional[FitResult] = None
) -> Optional[float]
```

**Parameters:**
- `frequencies` - Measurement frequencies [Hz]
- `Z` - Complex impedance [Ohm]
- `thickness_nm` - Known oxide thickness [nm] (e.g., from SEM/TEM)
- `area_cm2` - Electrode area [cm^2] (default: 1.0)
- `fit_result` - Result from `fit_equivalent_circuit()` (recommended)

**Returns:** Estimated relative permittivity epsilon_r, or `None` if failed.

**Example:**
```python
from eis_analysis.analysis import estimate_permittivity

# Thickness known from SEM measurement
eps_r = estimate_permittivity(freq, Z, thickness_nm=20.0, fit_result=result)
print(f"Permittivity: {eps_r:.1f}")  # e.g., 22.0
```

### OxideAnalysisResult

```python
@dataclass
class OxideAnalysisResult:
    capacitance: float          # Effective capacitance [F]
    capacitance_specific: float # Specific capacitance [F/cm^2]
    thickness_nm: float         # Oxide thickness [nm]
    element_type: str           # 'C', 'K', 'Q', or 'estimate'
    element_R: Optional[float]  # Associated resistance [Ohm]
    element_tau: Optional[float] # Time constant [s]
    element_params: Dict        # All element parameters
    # Brug (2D) comparison values - only for a dominant Q element
    # when the circuit contains a series resistance; otherwise None
    capacitance_brug: Optional[float]          # Brug C_eff [F]
    capacitance_specific_brug: Optional[float] # Brug C_eff/area [F/cm^2]
    thickness_brug_nm: Optional[float]         # Thickness from Brug C [nm]
```

---

## CLI Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--analyze-oxide` | - | Enable oxide layer analysis |
| `--epsilon-r` | 22.0 | Relative permittivity of oxide |
| `--area` | 1.0 | Electrode area [cm^2] |

**Note:** Area is automatically loaded from DTA file metadata whenever
`--area` is left at its default (1.0). An explicit `--area 1.0` cannot be
distinguished from the default and is therefore overridden by metadata;
any other explicit value takes precedence.

---

## How It Works

### 1. Element Detection

The function traverses the circuit tree and finds:
- `Parallel(R, C)` combinations (Voigt elements)
- `Parallel(R, Q)` combinations
- `K` elements (Voigt with tau parametrization)

Malformed or ambiguous circuits are handled with a warning:
- a K element with R <= 0 is skipped (C = tau/R is undefined there);
- when one parallel combination contains multiple R or multiple C/Q
  elements (e.g. `(R1 | R2 | C)`), the last one is used.

### 2. Dominant Element Selection

Selects element with **largest resistance R**.

**Physical justification:**
- Compact oxide = excellent insulator = dominant resistance barrier
- Other processes (double layer, pores) have smaller R
- Largest R typically corresponds to compact oxide layer

All candidate elements (type, R, C/Q, tau) are listed in the log together
with the stated selection assumption, so the choice can be verified.
**Caution:** a charge-transfer process can also have the largest R -
always check that the selected element represents the oxide.

### 3. Capacitance Extraction

**For C or K elements:**
```
C_eff = C  (direct)
C_eff = tau / R  (for K element)
```

**For Q elements (Hsu-Mansfeld formula):**
```
C_eff = (R * Q)^(1/n) / R
```

This is more accurate than simple `C_eff = Q` approximation.

When the circuit also contains a series resistance Rs, the Brug (2D)
estimate `C = Q^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)` and its thickness are
computed as well and returned in `capacitance_brug` /
`thickness_brug_nm` (see "2D vs 3D Distribution of Time Constants"
below). Hsu-Mansfeld remains the primary value.

### 4. Thickness Calculation

Parallel plate capacitor model:
```
d = epsilon_0 * epsilon_r / C_specific

where:
  epsilon_0 = 8.854e-14 F/cm
  C_specific = C / area [F/cm^2]
  d in cm, converted to nm (* 1e7)
```

### 5. Fallback Without a Fitted Circuit (Mode 2)

Without `fit_result`, capacitance is estimated as the **median** of
`C_i = -1 / (omega * Z'')` over the capacitive points in the top frequency
decade (`f >= f_max / 10`). Warnings are logged:

- always: for multilayer (series) systems the estimate yields the *series
  combination* of layer capacitances;
- when the per-point estimates spread by more than max/min = 1.2 across
  the decade - the assumption `omega*R*C >> 1` then does not hold within
  the decade and the estimate is unreliable.

If no capacitive point exists in the top decade, the single
highest-frequency point is used (pre-0.16.16 behavior).

---

## Supported Circuit Elements

| Element | Detection | Capacitance |
|---------|-----------|-------------|
| `R(x) \| C(y)` | Parallel R-C | C directly |
| `R(x) \| Q(Q, n)` | Parallel R-Q | Hsu-Mansfeld: (R*Q)^(1/n)/R |
| `K(R, tau)` | K element | tau/R |

**Note:** Series R elements are ignored (they don't form RC time constants).

---

## Examples

### Example 1: Simple Voigt Circuit

```python
from eis_analysis import load_data, fit_equivalent_circuit, R, C
from eis_analysis.analysis import analyze_oxide_layer

freq, Z = load_data('sample.DTA')

# Single Voigt element
circuit = R(100) - (R(5000) | C(1e-6))
result, Z_fit, _ = fit_equivalent_circuit(freq, Z, circuit)

oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
print(f"Type: {oxide.element_type}")      # 'C'
print(f"R: {oxide.element_R:.0f} Ohm")    # ~5000
print(f"C: {oxide.capacitance:.2e} F")    # ~1e-6
print(f"d: {oxide.thickness_nm:.1f} nm")  # ~19.5
```

### Example 2: Multiple Voigt Elements

```python
# Two Voigt elements - selects one with largest R
circuit = R(100) - (R(1000) | C(1e-6)) - (R(5000) | C(1e-5))
result, Z_fit, _ = fit_equivalent_circuit(freq, Z, circuit)

oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
# Selects R=5000 element (dominant barrier)
print(f"Selected R: {oxide.element_R:.0f} Ohm")  # 5000
```

### Example 3: Q Element

```python
from eis_analysis import Q

# Voigt with Q instead of ideal capacitor
circuit = R(100) - (R(5000) | Q(1e-6, 0.9))
result, Z_fit, _ = fit_equivalent_circuit(freq, Z, circuit)

oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
print(f"Type: {oxide.element_type}")           # 'Q'
print(f"C_eff (Hsu-Mansfeld): {oxide.capacitance:.2e} F")
print(f"C_eff (Brug):         {oxide.capacitance_brug:.2e} F")   # 2D comparison
print(f"d (H-M):  {oxide.thickness_nm:.1f} nm")
print(f"d (Brug): {oxide.thickness_brug_nm:.1f} nm")
```

### Example 4: K Element (tau parametrization)

```python
from eis_analysis import K

# Voigt with tau parametrization
circuit = R(100) - K(5000, 5e-3)  # R=5000, tau=5ms -> C=1uF
result, Z_fit, _ = fit_equivalent_circuit(freq, Z, circuit)

oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
print(f"Type: {oxide.element_type}")  # 'K'
print(f"tau: {oxide.element_tau:.2e} s")
```

### Example 5: Fallback (no fit_result)

```python
# Without fitted circuit - uses high-frequency estimate
# (median over the top frequency decade, see "How It Works" section 5)
oxide = analyze_oxide_layer(freq, Z, epsilon_r=22)
print(f"Type: {oxide.element_type}")  # 'estimate'
# Less accurate, use only for quick checks
```

### Example 6: Inverse - Permittivity from Known Thickness

```python
from eis_analysis.analysis import estimate_permittivity

# Thickness measured by SEM/TEM/ellipsometry
thickness_from_sem = 25.0  # nm

eps_r = estimate_permittivity(
    freq, Z,
    thickness_nm=thickness_from_sem,
    area_cm2=1.0,
    fit_result=result
)
print(f"Estimated permittivity: {eps_r:.1f}")
# Compare with literature: ZrO2 ~ 20-25
```

### Example 7: Cross-validation

```python
# Use both functions to cross-validate results

# Method 1: Assume epsilon_r, get thickness
oxide = analyze_oxide_layer(freq, Z, epsilon_r=22, fit_result=result)
print(f"Assuming eps_r=22: d = {oxide.thickness_nm:.1f} nm")

# Method 2: Use thickness from SEM, get epsilon_r
eps_r = estimate_permittivity(freq, Z, thickness_nm=30.0, fit_result=result)
print(f"From SEM d=30nm: eps_r = {eps_r:.1f}")

# If eps_r is reasonable (15-30 for ZrO2), the model is consistent
```

---

## Physical Background

### Oxide Layer as Capacitor

Compact oxide layer behaves as parallel plate capacitor:

```
C = epsilon_0 * epsilon_r * A / d
```

Rearranged for thickness (used by `analyze_oxide_layer`):
```
d = epsilon_0 * epsilon_r / C_specific
```

Rearranged for permittivity (used by `estimate_permittivity`):
```
epsilon_r = d * C_specific / epsilon_0
```

where `C_specific = C / A` is specific capacitance [F/cm^2].

### Typical Values for ZrO2

| Property | Value | Notes |
|----------|-------|-------|
| epsilon_r | 20-25 | Depends on phase, microstructure |
| C_specific | 1-100 uF/cm^2 | Depends on thickness |
| d | 1-100 nm | Typical passive films |

### Relationship C vs d

| C_specific [uF/cm^2] | d [nm] (epsilon_r=22) |
|---------------------|----------------------|
| 100 | 1.9 |
| 10 | 19.5 |
| 1 | 195 |

**Important:** Larger C = thinner oxide, Smaller C = thicker oxide.

---

## Q Effective Capacitance

### Why Q?

Real oxide layers often show non-ideal capacitive behavior due to:
- Surface roughness
- Inhomogeneous thickness
- Distributed time constants

Q impedance: `Z_Q = 1 / (Q * (j*omega)^n)`

### Hsu-Mansfeld Formula

For parallel R-Q circuit, effective capacitance:

```
C_eff = (R * Q)^(1/n) / R
```

derived from the characteristic time constant `tau = (R * Q)^(1/n)`.

This formula accounts for:
- Q coefficient Q
- Q exponent n (0 < n < 1)
- Parallel resistance R

### 2D vs 3D Distribution of Time Constants

CPE behavior means the sample has a *distribution* of time constants
rather than a single RC. The correct conversion from (Q, n) to an
effective capacitance depends on *where* that distribution physically
lives (classification by Hirschorn et al., 2010):

- **Surface (2D, lateral) distribution** - properties vary *along* the
  electrode surface (grain orientation, roughness, local coverage,
  lateral thickness variation). The appropriate conversion is the
  **Brug (1984)** formula, which involves the series (electrolyte)
  resistance Rs:

  ```
  C_eff = Q^(1/n) * (1/Rs + 1/Rct)^((n-1)/n)
  ```

- **Normal (3D, through-layer) distribution** - properties vary *across*
  the film thickness (resistivity/stoichiometry gradient from the
  metal/oxide interface to the outer surface). The appropriate
  conversion is the **Hsu-Mansfeld (2001)** formula, which involves the
  parallel (film) resistance R:

  ```
  C_eff = (R * Q)^(1/n) / R
  ```

These are not two approximations of the same quantity but two different
physical assumptions - using the wrong one gives a systematically
shifted value, not just a less precise one. The difference grows as n
decreases: for n close to 1 all conversions converge to C = Q, while
around n = 0.8 they can differ by tens of percent to a factor of a few,
which propagates directly into the thickness estimate.

**Which model for thin oxide layers?** For compact barrier oxides
(passive films, anodic oxides, the barrier layer on Zr alloys), the
dominant source of dispersion is usually the steep resistivity gradient
across the film - e.g. substoichiometric, more conductive oxide near
the metal interface vs nearly stoichiometric oxide outside. That is a
normal (3D) distribution, so **Hsu-Mansfeld is the physically better
default for oxide layers and is the primary value in this toolkit.**
Real films always contain some lateral (2D) component as well
(undulated metal/oxide interface, local thickness variation, cracks and
porosity after transition), which is why the Brug value is reported
alongside for comparison.

**Practical caveats:**

1. **Sensitivity to R.** An error in R propagates into C_eff with
   exponent (1-n)/n. For a protective oxide the film resistance is
   large (often 1e6-1e9 Ohm cm^2) and the low-frequency plateau may lie
   below the measured range, making R poorly conditioned. At n = 0.9 a
   factor-of-10 error in R gives ~30% error in C_eff; at n = 0.8
   already ~80%. Brug uses the well-determined Rs instead - a
   *practical* (not physical) advantage.

2. **Interpretation of the spread.** The difference between the
   Hsu-Mansfeld and Brug values brackets the model uncertainty of
   C_eff. If both agree within your tolerance, the choice of model does
   not matter for your data. If they diverge strongly, treat the
   absolute thickness with caution and prefer comparisons within a
   sample series (same model, same conditions), where the systematic
   model bias cancels.

3. **Neither formula is the last word.** For a realistic power-law
   resistivity profile across the film, Hirschorn et al. (2010) derived
   a more rigorous "power-law model"
   (`C_eff = g * Q * (rho_delta * eps * eps_0)^(1-n)`); Hsu-Mansfeld
   with the film resistance can misestimate thickness when the
   low-frequency R is not purely the dielectric film response (e.g.
   mixed with a faradaic process). The power-law model is not
   implemented in this toolkit.

4. **Cross-validate absolute values.** For Zr alloys, check EIS
   thickness at least once against weight gain (~15 mg/dm^2 per um of
   ZrO2) or metallography/SEM. Always report n: as n approaches 0.8 the
   2D/3D difference stops being academic and the thickness becomes soft.

A warning is logged when the dominant element has n < 0.8: the
distribution of time constants is then too broad for a single effective
capacitance to be well-defined, and the thickness estimate may be
unreliable.

**References:** Hsu & Mansfeld, Corrosion 57, 747 (2001);
Brug et al., J. Electroanal. Chem. 176, 275 (1984);
Hirschorn et al., Electrochim. Acta 55, 6218 (2010).

**Implementation notes:** The parallel resistance R needed for the
Hsu-Mansfeld conversion is always available: the dominant element is
selected only among candidates with R > 0, so the formula is used
unconditionally. (Fallback conversions for R-less CPEs existed up to
v0.16.16 but were unreachable and have been removed.) The Rs for the
Brug comparison is taken as the sum of R elements on the series path of
the fitted circuit (outside any parallel combination); if there is
none, the Brug estimate is skipped and the result fields stay None.

---

## Limitations

1. **Parallel plate model assumption**
   - Assumes homogeneous dielectric
   - Real oxides may have layered structure, porosity

2. **Permittivity uncertainty**
   - epsilon_r varies with oxide phase, microstructure
   - Porous oxide has lower effective epsilon_r
   - Uncertainty in epsilon_r directly affects thickness estimate

3. **Single dominant element**
   - Only analyzes element with largest R
   - For complex multilayer systems, manual analysis may be needed

4. **Fit quality dependency**
   - Poor fit = wrong capacitance = wrong thickness
   - Always check fit quality (residuals < 5%)

---

## Troubleshooting

**"No Voigt (R||C), K, or R||Q elements found in circuit"**
- Check circuit structure - needs Parallel(R, C), Parallel(R, Q), or K
- Series elements alone don't work

**Unrealistic thickness (< 1 nm or > 1000 nm)**
- Check epsilon_r value
- Check electrode area
- Verify fit quality
- Consider if selected element really represents oxide

**Different result than expected**
- Function selects largest R - verify this is correct element
- For multi-layer systems, dominant R may not be the layer of interest

---

## References

- Hsu & Mansfeld (2001): Q effective capacitance formula (3D distribution)
- Brug et al. (1984): alternative Q conversion for surface (2D) distribution
- Hirschorn et al., Electrochim. Acta 55, 6218 (2010): 2D/3D classification,
  power-law model for normal resistivity distributions
- Orazem & Tribollet (2008): "Electrochemical Impedance Spectroscopy"
- Bojinov et al.: "EIS of passive films on metals"

---

*Last updated: 2026-07-19*
