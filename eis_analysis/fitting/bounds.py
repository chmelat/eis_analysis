"""
Parameter bounds generation and validation for circuit fitting.

Provides physically reasonable bounds for electrochemical circuit parameters
and validation of fitted parameters against these bounds.

Author: EIS Analysis Toolkit
"""

import numpy as np
from typing import List, Tuple, Optional
from numpy.typing import NDArray

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

# A parameter whose bounds span more than this ratio is treated as a
# positive scale parameter (log-scale semantics): bound proximity is
# measured in decades and confidence intervals are computed in log space.
# 1e6 = 6 decades separates scale parameters (R, C, Q, tau, ...; >=9 decades
# each) from linear-range parameters like the CPE exponent n (0.3-1.0).
LOG_SCALE_BOUND_RATIO = 1e6


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


def classify_bound_status(
    value: float,
    lower: float,
    upper: float
) -> str:
    """Return 'lower', 'upper', or '' depending on whether `value` sits at a bound.

    Threshold: 1 decade on log scale (when bounds span >6 decades) or 1% of
    the range on linear scale.
    """
    if not (np.isfinite(lower) and np.isfinite(upper)):
        return ''
    if lower > 0 and upper / lower > LOG_SCALE_BOUND_RATIO:
        if value <= 0:
            return 'lower'
        log_val = np.log10(value)
        log_lo, log_hi = np.log10(lower), np.log10(upper)
        if log_val - log_lo < 1.0:
            return 'lower'
        if log_hi - log_val < 1.0:
            return 'upper'
    else:
        rng = upper - lower
        if not np.isfinite(rng) or rng <= 0:
            return ''
        if value - lower < 0.01 * rng:
            return 'lower'
        if upper - value < 0.01 * rng:
            return 'upper'
    return ''


def log_scale_ci_mask(
    lower_bounds: Optional[List[float]],
    upper_bounds: Optional[List[float]]
) -> Optional[List[bool]]:
    """
    Determine which parameters should get log-space confidence intervals.

    A parameter is a positive scale parameter -- its uncertainty is
    multiplicative, so the CI is computed in log space (see
    `compute_confidence_interval`) -- when its lower bound is strictly
    positive and its bounds span more than LOG_SCALE_BOUND_RATIO. This is
    the same criterion `classify_bound_status` uses for its log-scale
    branch, so CI semantics and bound-proximity semantics agree: R, C, Q,
    L, sigma, tau get log-space CIs; the CPE exponent n (0.3-1.0) keeps
    the symmetric linear-space CI.

    Parameters
    ----------
    lower_bounds, upper_bounds : list of float or None
        Bounds in full parameter space; None if no bound info available

    Returns
    -------
    mask : list of bool or None
        True for parameters with log-space CI; None if bounds unavailable
        (all CIs stay linear/symmetric)
    """
    if lower_bounds is None or upper_bounds is None:
        return None
    return [
        lb > 0 and ub / lb > LOG_SCALE_BOUND_RATIO
        for lb, ub in zip(lower_bounds, upper_bounds)
    ]


def build_bound_status(
    params_opt: NDArray[np.float64],
    lower_bounds: Optional[List[float]],
    upper_bounds: Optional[List[float]],
    fixed_params: Optional[List[bool]]
) -> List[str]:
    """
    Build per-parameter bound status for the full parameter vector.

    Single source of truth for the "parameter at a bound" diagnostic:
    both FitResult.bound_status and FitDiagnostics.bounds_warnings are
    derived from this vector, so they cannot disagree.

    Parameters
    ----------
    params_opt : ndarray
        All optimized parameters (including fixed)
    lower_bounds, upper_bounds : list of float or None
        Bounds in full parameter space; if None, no bound info is
        available and all statuses are ''
    fixed_params : list of bool or None
        Which parameters are fixed (True = fixed)

    Returns
    -------
    bound_status : list of str
        For each parameter: '' (interior), 'lower', 'upper'
        (per `classify_bound_status`), or 'fixed'
    """
    bound_status = []
    for i, value in enumerate(params_opt):
        if fixed_params is not None and i < len(fixed_params) and fixed_params[i]:
            bound_status.append('fixed')
        elif lower_bounds is not None and upper_bounds is not None:
            bound_status.append(classify_bound_status(
                float(value), lower_bounds[i], upper_bounds[i]
            ))
        else:
            bound_status.append('')
    return bound_status


__all__ = [
    'generate_simple_bounds',
    'build_bound_status',
    'classify_bound_status',
    'log_scale_ci_mask',
    'PARAMETER_BOUNDS',
    'LOG_SCALE_BOUND_RATIO',
]
