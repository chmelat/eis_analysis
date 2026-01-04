"""
Equivalent circuit fitting module using operator overloading.

This module completely replaces the old impedance.py string parser with
a modern operator overloading approach inspired by EISAnalysis.jl (Julia).

Refactored architecture (v2.2):
- circuit_elements.py: Circuit elements (R, C, Q, W, etc.)
- circuit_builder.py: Series and Parallel combinators
- circuit.py: Main fitting function using new parser
- auto_suggest.py: Automatic circuit suggestion from DRT
- config.py: Configuration constants with documentation

Public API
----------
Circuit Elements:
- R: Resistor
- C: Capacitor
- Q: Constant Phase Element (CPE)
- L: Inductor
- W: Warburg semi-infinite diffusion
- Wo: Warburg bounded diffusion
- K: Voigt element with tau parametrization (R, τ)

Main Functions:
- fit_equivalent_circuit: Fit circuit to impedance data
- suggest_circuit_from_drt: Suggest circuit from DRT analysis

Usage Example
-------------
```python
from eis_analysis.fitting import R, C, Q, fit_equivalent_circuit

# Build circuit using operators
circuit = R(100) - (R(5000) | C(1e-6))  # Voigt element

# Fit to data
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
print(result.params_opt)  # [98.5, 4823.2, 8.7e-7]
```

Operators:
- `-` : Series connection (Z_total = Z1 + Z2)
- `|` : Parallel connection (1/Z_total = 1/Z1 + 1/Z2)
- `*` : Scale parameter (2*R(100) = R(200))
- `**`: Exponent for Q/CPE (Q()**0.9)

Changes from v2.1:
------------------
- REMOVED: impedance.py dependency (no more string parsing!)
- NEW: Operator overloading approach for circuit building
- REMOVED: parse_user_initial_guess (no longer needed)
- REMOVED: estimate_initial_guess (values now in circuit object)
- REMOVED: generate_bounds (simplified for now)
- REMOVED: check_fit_quality (simplified for now)

The new approach is:
1. Simpler (no string parsing bugs)
2. More powerful (arbitrary nesting works)
3. Type-safe (Python type hints)
4. Elegant (looks like math notation)
"""

# Import circuit elements
from .circuit_elements import R, C, Q, L, W, Wo, K, CircuitElement

# Import circuit builders
from .circuit_builder import Series, Parallel, Circuit

# Import main fitting function
from .circuit import fit_equivalent_circuit, FitResult

# Import covariance for advanced users
from .covariance import CovarianceResult

# Import Voigt element analysis from DRT
from .auto_suggest import analyze_voigt_elements, format_voigt_report

# Import Voigt chain linear fitting
from .voigt_chain import (
    fit_voigt_chain_linear,
    generate_tau_grid,
    compute_voigt_matrix,
    estimate_R_linear,
)

# Import multi-start optimization
from .multistart import (
    fit_circuit_multistart,
    MultistartResult,
    perturb_from_covariance,
    perturb_from_stderr,
)

# Import differential evolution optimization
from .diffevo import (
    fit_circuit_diffevo,
    DiffEvoResult,
    DE_STRATEGIES,
)

__all__ = [
    # Circuit elements
    'R',
    'C',
    'Q',
    'L',
    'W',
    'Wo',
    'K',
    'CircuitElement',

    # Circuit builders
    'Series',
    'Parallel',
    'Circuit',

    # Main API
    'fit_equivalent_circuit',
    'FitResult',
    'CovarianceResult',
    'analyze_voigt_elements',
    'format_voigt_report',

    # Voigt chain linear fitting
    'fit_voigt_chain_linear',
    'generate_tau_grid',
    'compute_voigt_matrix',
    'estimate_R_linear',

    # Multi-start optimization
    'fit_circuit_multistart',
    'MultistartResult',
    'perturb_from_covariance',
    'perturb_from_stderr',

    # Differential evolution optimization
    'fit_circuit_diffevo',
    'DiffEvoResult',
    'DE_STRATEGIES',
]
