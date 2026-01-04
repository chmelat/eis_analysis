"""
Voigt chain fitting module - automatic initial guess using linear regression.

This module provides an alternative method for fitting EIS data with a chain
of Voigt elements using K(R, tau) parametrization. Instead of random initial
guesses, it uses a physics-informed approach:

1. Generate logarithmically-spaced time constants tau (2-3 per decade)
2. Fix tau values and solve for R_i using linear least squares on Z'(omega)
3. Build circuit using K(R_i, tau_i) elements
4. Use these as initial guess for final nonlinear optimization

Key advantages:
- Linear problem -> stable, fast convergence
- Uses real part Z' (more stable than Z'')
- Logarithmic tau distribution covers spectrum naturally
- Physics-informed: models actual relaxation processes

Theory
------
For a Voigt element (R || C), the impedance is:
    Z_i(omega) = R_i / (1 + j*omega*tau_i)    where tau_i = R_i * C_i

For a chain of N Voigt elements:
    Z(omega) = R_s + sum Z_i(omega)

The real part is:
    Z'(omega) = R_s + sum R_i * tau_i^2 * omega^2 / (1 + tau_i^2 * omega^2)

This is LINEAR in [R_s, R_1, R_2, ..., R_N] for fixed tau_i!

We can write:
    Z' = A @ R    where A[k, i] = tau_i^2 * omega_k^2 / (1 + tau_i^2 * omega_k^2)

This is solved using non-negative least squares (NNLS) since R_i >= 0.

Submodule structure
-------------------
- validation.py: Input data validation
- tau_grid.py: Tau grid generation
- solvers.py: NNLS solver and design matrix computation
- fitting.py: Main fitting functions (estimate_R_linear, fit_voigt_chain_linear)
- mu_optimization.py: Mu metric for overfit detection (Lin-KK style)

Author: EIS Analysis Toolkit
"""

# Public API - tau grid generation
from .tau_grid import (
    generate_tau_grid,
    generate_tau_grid_fixed_M,
)

# Public API - solvers
from .solvers import (
    compute_voigt_matrix,
)

# Public API - main fitting functions
from .fitting import (
    estimate_R_linear,
    fit_voigt_chain_linear,
    VoigtChain,
)

# Public API - mu optimization
from .mu_optimization import (
    calc_mu,
    find_optimal_M_mu,
)

__all__ = [
    # Tau grid generation
    'generate_tau_grid',
    'generate_tau_grid_fixed_M',

    # Solvers
    'compute_voigt_matrix',

    # Main fitting
    'estimate_R_linear',
    'fit_voigt_chain_linear',
    'VoigtChain',

    # Mu optimization
    'calc_mu',
    'find_optimal_M_mu',
]
