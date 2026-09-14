# Equivalent Circuit Parser

This document describes the implementation of the equivalent circuit parser
in the EIS Analysis Toolkit.

## Architecture

The parser uses an **operator overloading** approach inspired by the
EISAnalysis.jl library (Julia). Instead of a text parser, circuits are defined
directly in Python using operators.

### Key Modules

```
eis_analysis/fitting/
  circuit_elements/      # Element definitions (package)
    base.py              #   CircuitElement, operators, fixed parameters
    basic.py             #   R, C, L
    distributed.py       #   Q, W, Wo, CC
    composite.py         #   K, G
  circuit_builder.py     # Combinators (Series, Parallel)
  circuit.py             # Fitting functions
  bounds.py              # Default bounds per parameter label
```

### Classes

**CircuitElement** (abstract base class):
- Represents a single circuit element
- Defines operators `-` (series) and `|` (parallel)
- Stores parameters and their fixed/free status

**Series** and **Parallel** (composite circuits):
- Represent series/parallel connections
- Recursively contain other elements or composites
- Implement `impedance()`, `get_all_params()`, `update_params()`


## Parsing with eval()

Circuit strings are parsed by `parse_circuit_expression()` in
`eis_analysis/cli/utils.py`:

```python
def parse_circuit_expression(expr: str):
    safe_namespace = {
        'R': R, 'C': C, 'Q': Q, 'L': L, 'W': W,
        'Wo': Wo, 'K': K, 'G': G, 'CC': CC, 'DQ': DQ, 'YG': YG
    }
    circuit = eval(expr, {"__builtins__": {}}, safe_namespace)
    return circuit
```

Using `eval()` with a restricted namespace ensures:
- Access only to defined elements
- No built-in functions (security)
- Full Python syntax for expressions


## Operator Precedence

**Critical information:** Because Python's `eval()` is used, **Python's operator precedence**
applies, not the intuitive precedence for electrical circuits.

### Python Precedence (highest to lowest)

| Precedence | Operator | Meaning in circuits |
|------------|----------|-------------------|
| 12         | `-`      | **Series connection** |
| 8          | `|`      | **Parallel connection** |

**Operator `-` has HIGHER precedence than `|`!**

Only `-` and `|` are overloaded (`__sub__` and `__or__` in
`circuit_elements/base.py`). `+` is *not* a series connection - `R(100) +
C(1e-6)` raises `unsupported operand type(s) for +`.

### Consequences

Without parentheses, expressions are parsed as follows:

```python
# User writes:
R(1) - R(2)|C(3) - R(4)|C(5)

# Python groups every '-' first (higher precedence), then chains the '|'
# left to right (| is left-associative):
((R(1) - R(2)) | (C(3) - R(4))) | C(5)

# Resulting structure - the top-level node is Parallel, not Series:
Parallel
  |- Series: R(1) - R(2)
  |- Series: C(3) - R(4)
  '- C(5)
```

This is NOT the intended circuit `R - (R||C) - (R||C)`! Note that R(1) does
not stay in series with the rest at all - it ends up inside the first
parallel branch.

### Correct Notation

Always use explicit parentheses around parallel combinations:

```python
# Correct:
R(1) - (R(2)|C(3)) - (R(4)|C(5))

# Resulting structure:
# R(1) - [R(2)||C(3)] - [R(4)||C(5)]
```

### Recommendations

1. **Always** use parentheses around `(R|C)` combinations
2. Don't forget parentheses for complex parallel combinations
3. When in doubt, check the structure using `print(circuit)`


## Supported Elements

### R - Resistor

```python
Z_R = R
```

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| R         | Ohm  | 100     | Resistance  |

```python
R(100)      # 100 Ohm, free parameter
R("100")    # 100 Ohm, fixed parameter
R()         # default 100 Ohm
```

### G - Conductance

```python
Y_G = G
Z_G = 1 / G
```

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| G         | S    | 1e-6    | Conductance |

```python
G(1e-9)     # 1 nS = 1 GOhm, free parameter
G("0")      # fixed open circuit
G()         # default 1 uS = 1 MOhm
```

The same resistor as `R(1/G)`, parametrized by its admittance. Use it for a
parallel resistance the measured window may not resolve.

`R` is bounded at 10 GOhm, an edge of the parameter box rather than a feature
of the model. A resistance the data cannot determine from above is driven onto
it, `dZ/dR ~ 1/R^2` collapses, and the reported error becomes spurious
precision instead of a large uncertainty. The same limit written as `G` is
`G = 0`, which lies inside the allowed range and where `dZ/dG = -Z^2` stays
well conditioned:

```
G0 = 2.07e-10 +- 1.63e-10 S     ->  R > 2.4 GOhm
```

The interval is symmetric and may include zero. That is the result, not a
failure: it says the data place a lower bound on the resistance and nothing
more. Fit quality is unaffected - `R` and `G` give the same chi^2.

Applications: blocking coatings, intact oxide layers, any barrier whose
resistance may exceed what the frequency range can resolve. Prefer `R` when
the resistance is well determined, since it reads more directly.

`--analyze-oxide` treats `G` as the parallel resistance `R = 1/G`; `G = 0` is
reported as no DC path.

### C - Capacitor

```python
Z_C = 1 / (j * omega * C)
```

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| C         | F    | 1e-6    | Capacitance |

```python
C(1e-6)     # 1 uF
C("1e-6")   # fixed
```

### L - Inductor

```python
Z_L = j * omega * L
```

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| L         | H    | 1e-6    | Inductance  |

```python
L(1e-6)     # 1 uH
```

### Q - Constant Phase Element (CPE)

```python
Z_Q = 1 / (Q * (j * omega)^n)
```

| Parameter | Unit         | Default | Description  |
|-----------|--------------|---------|--------------|
| Q         | F * s^(n-1)  | 1e-4    | Q coefficient |
| n         | -            | 0.8     | Q exponent; fitted within 0.3 - 1.0 |

Special cases:
- n = 1: ideal capacitor
- n = 0.5: Warburg diffusion
- n = 0: ideal resistor

```python
Q(1e-4, 0.8)      # typical CPE
Q(1e-4, "0.9")    # Q free, n fixed
Q("1e-4", "0.9")  # both fixed
```

### W - Warburg (semi-infinite)

```python
Z_W = sigma / sqrt(omega) * (1 - j)
```

| Parameter | Unit         | Default | Description        |
|-----------|--------------|---------|---------------------|
| sigma     | Ohm*s^(-1/2) | 50      | Warburg coefficient |

```python
W(50)       # sigma = 50
```

### Wo - Warburg Open (bounded)

```python
Z_Wo = R_W * tanh(sqrt(j*omega*tau)) / sqrt(j*omega*tau)
```

| Parameter | Unit | Default | Description              |
|-----------|------|---------|--------------------------|
| R_W       | Ohm  | 100     | Warburg resistance       |
| tau_W     | s    | 1.0     | Diffusion time constant  |

```python
Wo(100, 1.0)    # R_W=100, tau=1s
```

### K - Voigt element (R||C with tau parametrization)

```python
Z_K = R / (1 + j*omega*tau)
```

Equivalent to (R || C) where C = tau/R.

| Parameter | Unit | Default | Description     |
|-----------|------|---------|-----------------|
| R         | Ohm  | 1000    | Resistance      |
| tau       | s    | 1e-4    | Time constant   |

Advantages of tau parametrization:
- tau directly determines the characteristic frequency: f = 1/(2*pi*tau)
- R and tau are more independent (better numerical conditioning)
- Consistent with DRT notation
- Used in Lin-KK test

```python
K(1000, 1e-4)   # R=1k, tau=100us, f=1.59kHz, C=100nF
```

### GE - Gerischer Element (reaction-diffusion)

```python
Z_G = sigma / sqrt(1 + j*omega*tau)
```

Models coupled diffusion with first-order chemical reaction.

| Parameter | Unit         | Default | Description              |
|-----------|--------------|---------|--------------------------|
| sigma     | Ohm          | 100     | Pre-factor (DC limit)    |
| tau       | s            | 1e-3    | Reaction time constant   |

Applications:
- SOFC cathodes (oxygen reduction reaction)
- Porous electrodes with surface reactions
- Mixed ionic-electronic conductors (MIECs)

Key differences from other elements:
- Unlike Warburg: has finite DC resistance (Z(0) = sigma)
- Unlike Voigt (K): asymmetric arc in Nyquist plot
- Characteristic frequency: f = 1/(2*pi*tau)

```python
GE(100, 1e-3)        # sigma=100, tau=1ms
GE("100", 1e-3)      # sigma fixed, tau free
GE("100", "1e-3")    # both fixed
```

> Renamed from `G` in v0.29.0. `G` is now the conductance element.

### CC - Cole-Cole Element (dielectric relaxation)

```python
C_star = C_inf + dC / (1 + (j*omega*tau)^(1-alpha))
Z_CC   = 1 / (j*omega*C_star)
```

Models a dielectric with a *distribution* of relaxation times. Unlike Q,
which describes a depressed arc in the impedance plane, CC describes a
depressed arc in the complex capacitance (equivalently permittivity) plane -
which is where the physics of an oxide or polymer film lives.

| Parameter | Unit | Default | Description                                 |
|-----------|------|---------|---------------------------------------------|
| C_inf     | F    | 1e-8    | High-frequency limit capacitance            |
| dC        | F    | 1e-7    | Relaxation strength, dC = C_s - C_inf       |
| tau       | s    | 1e-3    | Relaxation time                             |
| alpha     | -    | 0.2     | Broadening exponent, 0 <= alpha < 1         |

Applications:
- Oxide and passive films (dielectric dispersion)
- Polymers, glasses, biological tissue
- Any system where the permittivity, not the impedance, is the quantity of
  interest

Special cases:
- alpha = 0: Debye relaxation (a single relaxation time)
- larger alpha: broader distribution; the C* arc centre sits alpha*90deg
  below the real axis

Key properties:
- Static capacitance: C_s = C_inf + dC (available as `.C_static`)
- Peak-loss frequency: f = 1/(2*pi*tau) (available as `.characteristic_freq`)
- Blocking: |Z| diverges as omega -> 0, like a capacitor

The element is parametrised by dC rather than C_s so that dC > 0 (from the
bounds) enforces C_s > C_inf on its own; box bounds cannot express a coupling
between two parameters.

It carries no geometry - the parameters are capacitances, not permittivities.
Convert with eps_r = C*d/(eps_0*A) in the analysis layer: `--analyze-oxide`
recognises a CC and uses its static capacitance C_inf + dC, so
`--circuit "R(20) - CC(...)" --analyze-oxide --epsilon-r 22 --area 1.0`
reports the film thickness directly. See
[OXIDE_ANALYSIS_GUIDE.md](OXIDE_ANALYSIS_GUIDE.md).

Exactly equivalent to a composite of existing elements:

```python
CC(C_inf, dC, tau, alpha) == C(C_inf) | (C(dC) - Q(dC/tau**(1-alpha), alpha))
CC(C_inf, dC, tau, 0)     == C(C_inf) | (C(dC) - R(tau/dC))   # Debye
```

The composite is not a substitute for fitting: its CPE coefficient
`dC/tau^(1-alpha)` couples three parameters non-linearly, so a fit would
report Q and n with confidence intervals on the wrong quantities, and dC
would appear twice as two independent free parameters.

Generalising to a second exponent gives Havriliak-Negami; see
[LEVM_CIRCUITS.md](LEVM_CIRCUITS.md) (LEVM NDE = 6, 7).

```python
CC(1e-8, 1e-7, 1e-3, 0.2)         # depressed dielectric arc
CC(1e-8, 1e-7, 1e-3, 0.0)         # Debye limit
CC("1e-8", 1e-7, 1e-3, 0.2)       # C_inf fixed, rest free
R(10) - CC(1e-8, 1e-7, 1e-3, 0.2) # with series electrolyte resistance
```


### DQ - Truncated CPE (bounded power-law distribution)

```python
gamma(tau) = A*tau^n  for tau_min <= tau <= tau_max,  0 outside
Z_DQ(w)    = Int gamma(tau)/(1 + j*w*tau) d ln(tau)
```

An ideal CPE is the same power law with no bounds, which is what makes it
unphysical: gamma(tau) = A*tau^n is not normalisable, so a fit containing a
CPE has no DRT to recover and no DC limit, forcing a separate conductance
into the model. Truncating the distribution fixes both and gives three
regimes in one element:

| Frequency range | Behaviour |
|-----------------|-----------|
| below 1/tau_max | finite polarisation resistance R_pol |
| 1/tau_max ... 1/tau_min | CPE, slope -n |
| above 1/tau_min | capacitive, C_eff |

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| A | Ohm*s^-n | 1e-3 | Amplitude of the distribution |
| n | - | 0.6 | Power-law exponent, the CPE's n |
| tau_min | s | 1e-6 | Fast end of the distribution |
| U | - | 10.0 | Log-width, U = ln(tau_max/tau_min) |

Derived (properties, not fitted): `tau_max`, `R_pol`, `C_eff`.

    R_pol = A*(tau_max^n - tau_min^n)/n
    C_eff = (1-n)/(A*(tau_min^(n-1) - tau_max^(n-1)))

So one DQ does the work of `G | Q | C` - but with R_pol and C_eff derived
from the distribution rather than independent, which is what removes the
degeneracy that makes a fitted C drift with the frequency window.

The width is parametrised as U rather than tau_max so the bounds cannot
cross: box bounds cannot express tau_max > tau_min, and a fitter handed the
pair directly will swap them. U is capped at 30 (13 decades) by the
quadrature - see `DQ_QUAD_NODES` in `circuit_elements/distributed.py`.

Relation to the ideal CPE (the wide-bounds limit):

    A = sin(pi*n)/(pi*Q)

A bound outside the measured window is still identifiable, but only through
the power law: U correlates with n at -0.85, so read n beside the bound
status of U, never alone. `U_DQ [at upper bound]` is the informative result
that the slow end of the distribution lies past the measurement.

`--analyze-oxide` recognises a DQ and reads C_eff off it directly - no
Hsu-Mansfeld or Brug conversion, unlike a Q - and warns when the capacitive
plateau starts above the highest measured frequency, where C_eff is an
extrapolation.

```python
DQ(1.2e6, 0.57, 5e-2, 8)        # oxide film, distribution inside the window
DQ(1.2e6, "0.57", 5e-2, 8)      # exponent fixed, rest free
R(20) - DQ(1.2e6, 0.57, 5e-2, 8)
```

Concept from LEVM (Macdonald), where the truncation appears as the DWC
models with limits U1, U2; see [LEVM_CIRCUITS.md](LEVM_CIRCUITS.md).

### YG - Young-Göhr passive layer

```python
Z_YG(w) = p/(j*w*C) * ln[(1 + j*w*tau*e^(1/p)) / (1 + j*w*tau)]
```

A dielectric film whose conductivity penetrates from one side and decays
exponentially with depth. Where a CPE reports only an exponent, YG reports
the two quantities the film actually has, so it is the substitute to reach
for once the origin of the dispersion is known - oxide layers on Fe, Al, Ti
and Ta, and organic coatings under soaking.

| Parameter | Unit | Default | Description |
|-----------|------|---------|-------------|
| `C` | F | 1e-5 | total capacity of the layer, C = eps_0*eps_r*A/d; the high-frequency limit |
| `p` | - | 0.05 | relative penetration depth delta/d of the conductivity; p << 1 is a strong gradient |
| `tau` | s | 0.1 | time constant at the site of highest conductivity, tau = eps_0*eps_r*rho(0) |

Three regimes, as with DQ:

| Frequency range | Behaviour |
|-----------------|-----------|
| below e^(-1/p)/(2*pi*tau) | real resistance R_dc = p*tau*(e^(1/p) - 1)/C |
| e^(-1/p)/(2*pi*tau) ... 1/(2*pi*tau) | CPE-like, phase nearly constant |
| above 1/(2*pi*tau) | capacitive, C |

The middle band is the one the element exists to explain, and Zahner gives a
closed approximation for its phase, useful as a sanity check:

```python
phi = -90 deg * (1 - q),    q = 1 / (ln(w*tau) + 1/p)
```

The two corners are e^(1/p) apart, so for any realistic p the resistive one
lies decades below any sweep: **R_dc is a model extrapolation, not a measured
resistance**, and the oxide analysis reports it as such rather than ranking
elements by it.

p -> 0 degenerates to a plain capacitor, which is why p at its lower bound
(1e-3) means the parameter is unidentifiable rather than merely small. The
approach to that limit is continuous: the implementation never forms e^(1/p)
as a number - it overflows float64 below p = 1/709, which is *inside* the
bounds - so every legal p is evaluated on its merits, and the ideal-capacitor
short circuit only takes over below p = 1e-300, where 1/p stops being
representable at all.

```python
YG(1e-5, 0.05, 0.1)             # Zahner's simulated example
YG("1e-5", 0.05, 0.1)           # capacitance fixed, rest free
L(1e-6) - R(20) - YG(1e-5, 0.05, 0.1)
```

From H. Göhr; see Zahner Analysis manual 11/2023 section 2.3.9 and
[ZAHNER_ANALYSIS_REVIEW.md](ZAHNER_ANALYSIS_REVIEW.md) section 2.


## Fixed Parameters

Parameters passed as **strings** are automatically fixed during fitting:

```python
# Free parameters (will be fitted):
R(100)
Q(1e-4, 0.8)

# Fixed parameters (will not be fitted):
R("0.86")           # R_inf from previous measurement
Q("1e-4", 0.8)      # Q fixed, n free
Q("1e-4", "0.9")    # both fixed
```

Usage:
- Fixing R_inf from high-frequency measurement
- Fixing geometric parameters
- Sequential fitting (first part, then whole)


## Operators

### Series (`-`)

```python
Z_total = Z1 + Z2 + ... + Zn
```

```python
R(100) - C(1e-6)                    # R and C in series
R(10) - (R(100)|C(1e-6))            # R in series with Voigt
```

### Parallel (`|`)

```python
1/Z_total = 1/Z1 + 1/Z2 + ... + 1/Zn
```

```python
R(1000) | C(1e-6)                   # Voigt element
R(100) | Q(1e-4, 0.8)               # R||Q
```


## Circuit Examples

### Simple Voigt Element

```python
R(100) - (R(5000) | C(1e-6))
```

Structure: R_s - (R_p || C_p)

### Randles Circuit

```python
R(10) - ((R(100) - W(50)) | Q(1e-4, 0.8))
```

Structure: R_s - ((R_ct - W) || Q_dl)

### Two Voigt Elements (oxide layer)

```python
R("1.69") - (R(4000) | Q(5e-8, 0.98)) - (R(1200) | Q(3e-8, 0.96))
```

Structure: R_s(fix) - (R1 || Q1) - (R2 || Q2)

### K Element (alternative to Voigt)

```python
R(1) - K(1000, 1e-4) - K(500, 1e-3)
```

### Gerischer Element (SOFC cathode)

```python
R(10) - GE(100, 1e-3)
```

Structure: R_s - GE (series resistance + Gerischer reaction-diffusion)

### Mixed Voigt + Gerischer

```python
R(1) - K(500, 1e-4) - GE(100, 1e-3)
```

Structure: R_s - (R||C) - G (electrolyte + charge transfer + reaction-diffusion)


## Internal Representation

After parsing, a tree structure is created:

```python
circuit = R(100) - (R(5000) | C(1e-6))
print(type(circuit))  # <class 'Series'>
print(circuit)        # R(100) - (R(5000) | C(1e-6))
```

Internal structure:
```
Series
  |- R(100)
  '- Parallel
       |- R(5000)
       '- C(1e-6)
```

### Methods

```python
# All parameters (for fitting)
circuit.get_all_params()        # [100, 5000, 1e-6]

# Which are fixed
circuit.get_all_fixed_params()  # [False, False, False]

# Parameter names
circuit.get_param_labels()      # ['R', 'R', 'C']

# Calculate impedance
Z = circuit.impedance(frequencies, params)

# Update after fitting
circuit.update_params(fitted_params)
```

**Note:** the labels are the symbols the fit output prints, not the argument
names used in this document's parameter tables. They are Greek where the
symbol is: `W` -> `σ`, `Wo` -> `R_W`, `τ_W`, `Q` -> `Q`, `n`, `K` -> `R`, `τ`,
`G` -> `σ_G`, `τ_G`, `CC` -> `C_inf`, `ΔC`, `τ_CC`, `α_CC`. `bounds.py` keys
its default bounds on these labels.


## Limitations of Current Implementation

### 1. Operator Precedence

As described above, `-` has higher precedence than `|`. Solution: use parentheses.

### 2. No Bounds in Syntax

Bounds for parameters cannot be specified in the circuit string.
Default bounds are used based on element type.

### 3. No Structure Validation

Parser does not validate if circuit is physically meaningful.
For example, `C(1e-6) | C(1e-6)` is syntactically correct but meaningless.

### 4. Order Dependency

Parameter order in `get_all_params()` depends on element order
in the expression. Changing structure changes parameter indices.


## Security

**Only ever pass trusted circuit strings to the parser.** A circuit
expression is executed, not parsed, so it must be treated like any other
piece of Python source the user hands over.

The parser calls `eval()` with a restricted namespace:

```python
eval(expr, {"__builtins__": {}}, safe_namespace)
```

- `__builtins__: {}` - built-in names are not bound
- `safe_namespace` - only circuit elements are bound

**This is not a sandbox.** Emptying `__builtins__` removes the convenient
names (`open`, `__import__`, `eval`), but not the objects they lead to: a
literal such as `()` still exposes the whole class hierarchy through
`__class__`, and any class found that way carries a `__globals__` with a live
reference to the real built-ins. Arbitrary code execution, file access and
imports are therefore all reachable from a deliberately written expression.

The restriction is worth keeping - it turns a typo into an error instead of a
surprise - but it defends against accidents, not against a hostile string. Do
not expose the parser to untrusted input (a web form, a shared job queue)
without running it in a real sandbox: a separate process with dropped
privileges, or a proper grammar-based parser (see Possible Future
Extensions), which would remove `eval()` altogether.


## Possible Future Extensions

### Custom Parser

Implementing a recursive descent parser would enable:
- Correct operator precedence (| > -)
- Structure validation
- Better error messages
- Bounds in syntax: `R(100, bounds=(10, 1000))`

Grammar:
```
expr   := term ('-' term)*
term   := factor ('|' factor)*
factor := element | '(' expr ')'
element := R(...) | C(...) | Q(...) | ...
```

Estimated complexity: ~100 lines of code.

### Alternative Operators

Changing operators to those with opposite precedence in Python:

| New | Old | Precedence |
|-----|-----|------------|
| `*` | `|` | 13 (higher) |
| `+` | `-` | 12 (lower)  |

Breaking change for users.


## References

- Schonleber, M. et al. "A Method for Improving the Robustness of linear
  Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
- Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance
  Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
