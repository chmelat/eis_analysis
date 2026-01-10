"""
Parameter bounds generation and validation for circuit fitting.

Provides physically reasonable bounds for electrochemical circuit parameters
and validation of fitted parameters against these bounds.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import List, Tuple, Optional
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# Physically reasonable ranges for electrochemical systems
PARAMETER_BOUNDS = {
    # Resistance: 0.1 mOhm - 10 TOhm
    'R': (1e-4, 1e10),

    # Capacitance: 1 fF - 100 mF
    'C': (1e-15, 1e-1),

    # Inductance: 1 pH - 100 uH (parasitic in EIS)
    'L': (1e-12, 1e-4),

    # Q (CPE) coefficient: similar to C
    'Q': (1e-12, 1e-1),

    # Q exponent: 0.3 (strongly inhomogeneous) to 1.0 (ideal C)
    'n': (0.3, 1.0),

    # Warburg coefficient: 0.01 - 100000 Ohm*s^(-1/2)
    'σ': (1e-2, 1e5),

    # Time constant: 1 ns - 10000 s
    'τ': (1e-9, 1e4),

    # Warburg bounded - resistance
    'R_W': (1e-2, 1e8),

    # Warburg bounded - diffusion time
    'τ_W': (1e-6, 1e4),

    # Gerischer element - pre-factor (similar to Warburg)
    'σ_G': (1e-2, 1e8),

    # Gerischer element - reaction time constant
    'τ_G': (1e-9, 1e4),
}

DEFAULT_BOUNDS = (1e-15, 1e15)


def generate_simple_bounds(param_labels: List[str]) -> Tuple[List[float], List[float]]:
    """
    Generate physically reasonable bounds for parameters based on their labels.

    Bounds are narrowed to realistic ranges for electrochemical systems,
    helping the optimizer avoid non-physical values.

    Ranges:
    - R (resistance): 0.1 mOhm - 10 GOhm (batteries to coatings)
    - C (capacitance): 1 fF - 100 mF (parasitic to large electrode)
    - L (inductance): 1 pH - 100 uH (parasitic, cabling)
    - Q (CPE coefficient): 1 pF*s^(n-1) - 100 mF*s^(n-1)
    - n (Q/CPE exponent): 0.3 - 1.0
    - sigma (Warburg): 0.01 - 100000 Ohm*s^(-1/2)
    - tau (time constant): 1 ns - 10000 s (covers mHz-GHz)
    - R_W, tau_W (Warburg bounded): similar to R, tau
    - sigma_G, tau_G (Gerischer): similar to sigma, tau

    Parameters
    ----------
    param_labels : list of str
        Parameter labels (e.g., ['R', 'R', 'C', 'Q', 'n'])

    Returns
    -------
    lower_bounds : list of float
        Lower bounds
    upper_bounds : list of float
        Upper bounds

    Notes
    -----
    Narrowed bounds (vs 1e-12 to 1e12):
    - Prevent convergence to non-physical values
    - Improve numerical stability (smaller condition number)
    - Help with poor initial guess

    If a legitimate parameter hits a bound, the optimizer warns.
    Bounds can be extended manually if needed.
    """
    lower_bounds = []
    upper_bounds = []

    for label in param_labels:
        lb, ub = PARAMETER_BOUNDS.get(label, DEFAULT_BOUNDS)
        lower_bounds.append(lb)
        upper_bounds.append(ub)

    return lower_bounds, upper_bounds


def check_bounds_proximity(
    params_opt: NDArray[np.float64],
    params_opt_free: NDArray[np.float64],
    lower_bounds: List[float],
    upper_bounds: List[float],
    fixed_params: Optional[List[bool]]
) -> None:
    """
    Check and warn if parameters are near bounds.

    Parameters
    ----------
    params_opt : ndarray
        All optimized parameters (including fixed)
    params_opt_free : ndarray
        Only free (optimized) parameters
    lower_bounds : list of float
        Lower bounds for free parameters
    upper_bounds : list of float
        Upper bounds for free parameters
    fixed_params : list of bool or None
        Which parameters are fixed
    """
    if lower_bounds is None or upper_bounds is None:
        return

    near_lower, near_upper = [], []

    if fixed_params is not None and any(fixed_params):
        free_to_full_idx = [i for i, is_fixed in enumerate(fixed_params) if not is_fixed]
    else:
        free_to_full_idx = list(range(len(params_opt)))

    for free_idx in range(len(params_opt_free)):
        full_idx = free_to_full_idx[free_idx]
        param_val = params_opt_free[free_idx]
        lb, ub = lower_bounds[free_idx], upper_bounds[free_idx]

        if lb > 0 and ub / lb > 1e6:
            # Logarithmic check for wide bounds
            log_param = np.log10(param_val) if param_val > 0 else -np.inf
            log_lower, log_upper = np.log10(lb), np.log10(ub)
            if log_param - log_lower < 1.0:
                near_lower.append((full_idx, param_val, lb))
            elif log_upper - log_param < 1.0:
                near_upper.append((full_idx, param_val, ub))
        else:
            # Linear check for narrow bounds
            bound_range = ub - lb
            if np.isfinite(bound_range):
                if param_val - lb < 0.01 * bound_range:
                    near_lower.append((full_idx, param_val, lb))
                elif ub - param_val < 0.01 * bound_range:
                    near_upper.append((full_idx, param_val, ub))

    if near_lower or near_upper:
        logger.warning("=" * 50)
        logger.warning("WARNING: Some parameters are at or near bounds!")
        for idx, val, bound in near_lower:
            logger.warning(f"  Parameter {idx}: {val:.3e} near LOWER bound {bound:.3e}")
        for idx, val, bound in near_upper:
            logger.warning(f"  Parameter {idx}: {val:.3e} near UPPER bound {bound:.3e}")
        logger.warning("  Recommendation: Try simpler circuit or extend bounds")
        logger.warning("=" * 50)


__all__ = [
    'generate_simple_bounds',
    'check_bounds_proximity',
    'PARAMETER_BOUNDS',
]
