# K Element - Voigt element with tau parametrization

## Overview

The **K(R, tau)** element is an alternative parametrization of the Voigt
element (R||C) that uses the time constant tau instead of the capacitance C.

**Impedance:**
```
Z_K(omega) = R / (1 + j*omega*tau)
```

**Equivalence:**
```
K(R, tau) == (R || C)   where   C = tau/R
```

---

## Why K(R, tau) rather than (R||C)

### 1. Physical readability

```python
# Classic Voigt: what is the characteristic frequency?
circuit = R(1000) | C(1e-7)  # tau = R*C = 1e-4 s -> f = ? (compute it first)

# K element: the frequency is right there
circuit = K(1000, 1e-4)      # tau = 1e-4 s -> f = 1/(2*pi*tau) = 1.59 kHz
```

**What tau means:**
- tau sets the characteristic frequency: **f_c = 1/(2*pi*tau)**
- At f = f_c the phase is -45 deg and |Z| = R/sqrt(2) ~ 0.707 R. This is the
  half-*power* point: the real and imaginary parts are both R/2, so |Z| is
  0.707 R, not R/2.
- tau is the relaxation time of the process

### 2. Better parameter separation when fitting

```
With (R||C), R and the characteristic frequency are coupled:
  R = 1000 Ohm, C = 1e-7 F  ->  tau = 1e-4 s, f_c = 1.6 kHz
  R = 2000 Ohm, C = 1e-7 F  ->  tau = 2e-4 s, f_c = 800 Hz   (frequency moved!)

With K(R, tau) they are not:
  K(1000, 1e-4)  ->  f_c = 1.6 kHz
  K(2000, 1e-4)  ->  f_c = 1.6 kHz   (frequency unchanged)
```

**Each parameter has one job:**
- **R** sets the amplitude (the height of the Nyquist semicircle)
- **tau** sets the frequency (its position on the frequency axis)

### 3. Consistency with DRT analysis

DRT (Distribution of Relaxation Times) uses tau as its primary variable:

```
gamma(tau) = distribution of time constants
```

**The K element is the natural bridge between DRT and circuit fitting:**
- A DRT peak at tau_i becomes a K(R_i, tau_i) element
- Integrating the peak gives R_i = INT gamma(tau) d(ln tau)
- The result is directly usable as a circuit

### 4. Lin-KK compatibility

The Lin-KK test (Schonleber et al. 2014) uses exactly this parametrization:

```
Z_LinKK = R_0 + SUM K(R_k, tau_k)
```

so results from the Lin-KK test can be carried over directly.

---

## Usage

### Basic syntax

```python
from eis_analysis.fitting import K, R

# Create a K element
k = K(1000, 1e-4)  # R = 1 kOhm, tau = 100 us

# Defaults
k = K()  # R = 1 kOhm, tau = 100 us

# Fixed parameters (passed as strings)
k = K("1000", 1e-4)    # R fixed, tau free
k = K("1000", "1e-4")  # both fixed
```

### Properties

```python
k = K(1000, 1e-4)

# Equivalent capacitance
C = k.capacitance           # 1e-7 F (tau/R)

# Characteristic frequency
f_c = k.characteristic_freq # 1591.5 Hz (1/(2*pi*tau))

# Conversion to (R||C)
rc_circuit = k.to_RC()      # (R(1000) | C(1e-07))
```

### Computing impedance

```python
import numpy as np

freq = np.logspace(4, -1, 50)  # 10 kHz -> 0.1 Hz

# Directly
Z = k.impedance(freq, [1000, 1e-4])

# Or through a circuit object
circuit = R(100) - K(1000, 1e-4)
params = circuit.get_all_params()
Z = circuit.impedance(freq, params)
```

---

## Examples

### Example 1: Single Voigt element

```python
from eis_analysis.fitting import R, K, fit_equivalent_circuit
import numpy as np

# Circuit: R_s - K(R, tau)
circuit = R(100) - K(1000, 1e-4)

# Fit the data
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Fitted parameters
R_s, R_fit, tau = result.params_opt
C = tau / R_fit  # equivalent capacitance

print(f"R_s = {R_s:.2f} Ohm")
print(f"R   = {R_fit:.2f} Ohm")
print(f"tau = {tau:.3e} s")
print(f"f_c = {1/(2*np.pi*tau):.1f} Hz")
print(f"C   = {C:.3e} F")
```

### Example 2: Voigt chain (K elements in series)

```python
# Circuit: R_s - K_1 - K_2 - K_3
# Each K element models one relaxation process
circuit = (R(100) -
           K(500, 1e-5) -   # fast process   (~15.9 kHz)
           K(1000, 1e-4) -  # medium process (~1.6 kHz)
           K(2000, 1e-3))   # slow process   (~159 Hz)

# Fit (the default modulus weighting is right for most data;
# see WEIGHTING_AND_STATISTICS.md before changing it)
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Process analysis
params = result.params_opt
R_s = params[0]
for i in range(3):
    R_i = params[1 + 2*i]
    tau_i = params[2 + 2*i]
    f_i = 1/(2*np.pi*tau_i)
    print(f"Process {i+1}: R={R_i:.1f} Ohm, tau={tau_i:.3e} s, f={f_i:.1f} Hz")
```

### Example 3: Building a circuit from DRT results

```python
# Suppose the DRT analysis found three peaks
drt_peaks = [
    {'tau': 1e-5, 'R': 500},   # peak 1
    {'tau': 1e-4, 'R': 1000},  # peak 2
    {'tau': 1e-3, 'R': 2000},  # peak 3
]

circuit = R(100)  # R_s
for peak in drt_peaks:
    circuit = circuit - K(peak['R'], peak['tau'])

print(circuit)
# R(100) - K(R=500, τ=1e-05) - K(R=1000, τ=0.0001) - K(R=2000, τ=0.001)

# Use it as the initial guess for fitting
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

### Example 4: Randles circuit with a K element

```python
from eis_analysis.fitting import Q, W

# R_s - (K || Q) - W
# K models charge transfer, Q the non-ideal double layer, W diffusion
circuit = R(10) - (K(100, 1e-3) | Q(1e-4, 0.85)) - W(50)

result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

### Example 5: Fixing a parameter

```python
# Scenario: tau is known from DRT, only R should be fitted

circuit = R(100) - K(1000, "1e-4")  # tau fixed (string!)

result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Only R_s and R are optimized, but params_opt is always the FULL parameter
# vector - the fixed tau is still in it, and its standard error is 0:
R_s_fit, R_fit, tau_fixed = result.params_opt   # tau_fixed == 1e-4
```

---

## Relation to the (R||C) element

### Mathematical equivalence

```
K(R, tau):   Z = R / (1 + j*omega*tau)
R||C:        Z = R / (1 + j*omega*R*C)
```

**Equivalence condition:**
```
tau = R * C
```

### Conversion

**From K to (R||C):**
```python
k = K(1000, 1e-4)
C = k.capacitance  # C = tau/R = 1e-7 F
rc = k.to_RC()     # (R(1000) | C(1e-07))
```

**From (R||C) to K:**
```python
rc = R(1000) | C(1e-7)
tau = 1000 * 1e-7  # tau = R*C = 1e-4 s
k = K(1000, 1e-4)
```

---

## When to use K vs (R||C)

### Use K(R, tau) when:

- **You have DRT results** - tau comes straight from the peaks
- **You know the characteristic frequency** - tau = 1/(2*pi*f)
- **You want amplitude and frequency separated** - better conditioning
- **You work with the Lin-KK test** - same parametrization as the literature
- **You need the physical interpretation** - tau is the relaxation time

### Use (R||C) when:

- **You have an electrochemical model** - C is the physical quantity (double
  layer, oxide capacitance)
- **The capacitance is known** - e.g. from electrode geometry
- **You need compatibility with older work** - the classic parametrization
- **You are teaching** - R||C is the standard notation

---

## Worked conversions

### Between parameters

```
Given K(R = 1000 Ohm, tau = 1e-4 s):

  C       = tau / R           = 1e-7 F = 100 nF
  f_c     = 1 / (2*pi*tau)    = 1591.5 Hz
  omega_c = 1 / tau           = 10000 rad/s

At f = f_c:
  Z       = R / (1 + j) = R/2 * (1 - j) = 500 - 500j Ohm
  |Z|     = R / sqrt(2)       = 707.1 Ohm
  phase   = -45 deg
```

### Time constants for given frequencies

```python
# A K element with a characteristic frequency of 1 kHz
f_c = 1000  # Hz
tau = 1 / (2 * np.pi * f_c)  # 1.59e-4 s

# Different R, all with f_c = 1 kHz
K(100, tau)    # C = 1.59 uF
K(1000, tau)   # C = 159 nF
K(10000, tau)  # C = 15.9 nF
```

### Voigt chain covering decades

```python
# Four frequency decades (10 Hz - 10 kHz), one K element per decade
f_values = [10, 100, 1000, 10000]  # Hz
tau_values = [1/(2*np.pi*f) for f in f_values]

circuit = R(100)
for tau in tau_values:
    circuit = circuit - K(1000, tau)  # same R for all

# tau_values ~ [0.016, 0.0016, 0.00016, 0.000016] s
```

---

## Implementation notes

### Numerical stability

The K(R, tau) parametrization is generally better conditioned than (R||C):

```
(R||C) with a large R:
  R = 1e6 Ohm, C = 1e-10 F  ->  tau = 1e-4 s
  The gradient with respect to C is tiny, and C itself is 10 orders of
  magnitude away from R - the Jacobian columns differ by that much.

K(R, tau):
  K(1e6, 1e-4)
  The gradient with respect to tau is far better scaled.
```

(The fitter column-scales the Jacobian before computing the covariance, so
the unit disparity does not corrupt the reported standard errors either -
see WEIGHTING_AND_STATISTICS.md, section 3.5.)

### Default bounds

Bounds cannot be given in the circuit syntax; the fitter uses defaults keyed
on the parameter label (`fitting/bounds.py`):

```
R    : 1e-4  - 1e10 Ohm   (0.1 mOhm - 10 TOhm)
tau  : 1e-9  - 1e4  s     (1 ns - 10000 s, ~13 decades of frequency)
```

The tau bounds cover f_c from ~16 uHz to ~160 MHz, comfortably beyond any
real measurement. Note that the equivalent bounds on a capacitance would have
to be derived from the ratio tau/R, which spans a far wider (and largely
meaningless) range than the C bounds used for a real (R||C) - one more reason
the K parametrization is easier to bound.

---

## References

1. Schonleber, M. et al. "A Method for Improving the Robustness of linear
   Kramers-Kronig Validity Tests." *Electrochimica Acta* 131, 20-27 (2014).
   doi: 10.1016/j.electacta.2014.01.034

2. Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance Data
   Validation." *J. Electrochem. Soc.* 142 (6), 1885-1894 (1995).

3. Ciucci, F. "Modeling electrochemical impedance spectroscopy." *Current
   Opinion in Electrochemistry* 13, 132-139 (2019).
