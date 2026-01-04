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
- For Q: uses Brug formula for effective capacitance
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
python eis.py data.DTA --circuit 'R(100) - (R(5000) | C(1e-6))' --analyze-oxide

# With custom permittivity and area
python eis.py data.DTA --circuit 'R(100) - (R(5000) | C(1e-6))' \
    --analyze-oxide --epsilon-r 23 --area 0.5

# Using Voigt chain for automatic circuit
python eis.py data.DTA --voigt-chain --analyze-oxide
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
```

---

## CLI Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--analyze-oxide` | - | Enable oxide layer analysis |
| `--epsilon-r` | 22.0 | Relative permittivity of oxide |
| `--area` | 1.0 | Electrode area [cm^2] |

**Note:** Area can be automatically loaded from DTA file metadata if available.

---

## How It Works

### 1. Element Detection

The function traverses the circuit tree and finds:
- `Parallel(R, C)` combinations (Voigt elements)
- `Parallel(R, Q)` combinations
- `K` elements (Voigt with tau parametrization)

### 2. Dominant Element Selection

Selects element with **largest resistance R**.

**Physical justification:**
- Compact oxide = excellent insulator = dominant resistance barrier
- Other processes (double layer, pores) have smaller R
- Largest R typically corresponds to compact oxide layer

### 3. Capacitance Extraction

**For C or K elements:**
```
C_eff = C  (direct)
C_eff = tau / R  (for K element)
```

**For Q elements (Brug formula):**
```
C_eff = (R * Q)^(1/n) / R
```

This is more accurate than simple `C_eff = Q` approximation.

### 4. Thickness Calculation

Parallel plate capacitor model:
```
d = epsilon_0 * epsilon_r / C_specific

where:
  epsilon_0 = 8.854e-14 F/cm
  C_specific = C / area [F/cm^2]
  d in cm, converted to nm (* 1e7)
```

---

## Supported Circuit Elements

| Element | Detection | Capacitance |
|---------|-----------|-------------|
| `R(x) \| C(y)` | Parallel R-C | C directly |
| `R(x) \| Q(Q, n)` | Parallel R-Q | Brug: (R*Q)^(1/n)/R |
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
print(f"C_eff (Brug): {oxide.capacitance:.2e} F")
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

### Brug Formula

For parallel R-Q circuit, effective capacitance:

```
C_eff = (R * Q)^(1/n) / R
```

This formula accounts for:
- Q coefficient Q
- Q exponent n (0 < n < 1)
- Parallel resistance R

**Reference:** Brug et al., J. Electroanal. Chem. 176, 275 (1984)

### Fallback Methods

If R is unknown:
1. From Z'' maximum: `C_eff = Q * omega_max^(n-1)`
2. At 1 kHz reference: `C_eff = Q * (2*pi*1000)^(n-1)`

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

**"No Voigt elements found"**
- Check circuit structure - needs Parallel(R, C) or Parallel(R, Q)
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

- Brug et al. (1984): Q effective capacitance formula
- Orazem & Tribollet (2008): "Electrochemical Impedance Spectroscopy"
- Bojinov et al.: "EIS of passive films on metals"

---

*Last updated: 2025-12-26*
