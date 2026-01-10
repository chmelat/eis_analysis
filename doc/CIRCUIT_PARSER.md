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
├── circuit_elements.py   # Basic elements (R, C, Q, L, W, Wo, K, G)
├── circuit_builder.py    # Combinators (Series, Parallel)
└── circuit.py            # Fitting functions
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

Circuit strings are parsed in the `eis.py` file using the `parse_circuit_expression()` function:

```python
def parse_circuit_expression(expr: str):
    safe_namespace = {
        'R': R, 'C': C, 'Q': Q, 'L': L,
        'W': W, 'Wo': Wo, 'K': K, 'G': G
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
| 13         | `*`, `/` | Parameter scaling |
| 12         | `-`, `+` | **Series connection** |
| 8          | `|`      | **Parallel connection** |

**Operator `-` has HIGHER precedence than `|`!**

### Consequences

Without parentheses, expressions are parsed as follows:

```python
# User writes:
R(1) - R(2)|C(3) - R(4)|C(5)

# Python interprets (due to precedence - > |):
R(1) - (R(2) | (C(3) - (R(4) | C(5))))

# Resulting structure:
#   R(1)
#     │
#   Series
#     │
#   R(2) ══╗
#          ║ Parallel
#   C(3)   ║
#     │    ║
#   Series ║
#     │    ║
#   R(4) ══╬══╗
#          ║  ║ Parallel
#   C(5) ══╝══╝
```

This is NOT the intended circuit `R - (R||C) - (R||C)`!

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
| n         | -            | 0.8     | Q exponent (0-1) |

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

### G - Gerischer Element (reaction-diffusion)

```python
Z_G = sigma / sqrt(1 + j*omega*tau)
```

Models coupled diffusion with first-order chemical reaction.

| Parameter | Unit         | Default | Description              |
|-----------|--------------|---------|--------------------------|
| sigma     | Ohm*s^(1/2)  | 100     | Pre-factor (DC limit)    |
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
G(100, 1e-3)        # sigma=100, tau=1ms
G("100", 1e-3)      # sigma fixed, tau free
G("100", "1e-3")    # both fixed
```


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

### Scaling (`*`)

```python
2 * R(100)      # = R(200)
0.5 * C(1e-6)   # = C(5e-7)
```

### Exponent (`**`) - Q only

```python
Q(1e-4) ** 0.9    # changes n to 0.9
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
R(10) - G(100, 1e-3)
```

Structure: R_s - G (series resistance + Gerischer reaction-diffusion)

### Mixed Voigt + Gerischer

```python
R(1) - K(500, 1e-4) - G(100, 1e-3)
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
├── R(100)
└── Parallel
    ├── R(5000)
    └── C(1e-6)
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

Parser uses `eval()` with restricted namespace:

```python
eval(expr, {"__builtins__": {}}, safe_namespace)
```

- `__builtins__: {}` - no built-in functions
- `safe_namespace` - only circuit elements

This prevents:
- Execution of arbitrary code
- Access to file system
- Module imports

Nevertheless, it is recommended to use only trusted inputs.


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
