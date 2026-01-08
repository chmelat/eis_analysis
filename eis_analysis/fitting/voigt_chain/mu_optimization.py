"""
Mu metric optimization for Voigt chain fitting (Lin-KK style).

This module provides the mu metric calculation and automatic optimization
of the number of Voigt elements.
"""

import numpy as np
import logging
from typing import Tuple, Optional
from numpy.typing import NDArray

from .validation import validate_eis_data
from .tau_grid import generate_tau_grid_fixed_M
from .fitting import estimate_R_linear

logger = logging.getLogger(__name__)


def calc_mu(R_i: NDArray[np.float64]) -> float:
    """
    Calculate mu metric for overfit detection (Lin-KK style).

    mu = 1 - (sum of negative R) / (sum of positive R)

    Parameters
    ----------
    R_i : ndarray of float
        Resistance values (without R_s)

    Returns
    -------
    mu : float
        Mu metric
        - mu -> 1.0: all R_i positive (good fit, no overfit)
        - mu -> 0.0: large negative mass (overfit)
        - mu < 0.0: dominance of negative R_i (very bad)

    Notes
    -----
    In Lin-KK test, used as stop condition for iterative determination
    of optimal number of RC elements. Typical threshold: mu < 0.85.

    If all R_i >= 0 (from NNLS), then mu = 1.0 (ideal).

    Examples
    --------
    >>> R_i = np.array([100, 500, 1000])
    >>> calc_mu(R_i)
    1.0  # All positive

    >>> R_i = np.array([100, -50, 1000])
    >>> calc_mu(R_i)
    0.9545  # Small negative mass
    """
    neg_sum = np.sum(np.abs(R_i[R_i < 0]))
    pos_sum = np.sum(np.abs(R_i[R_i >= 0]))

    if pos_sum == 0:
        # All negative (very bad)
        logger.warning("calc_mu: all R_i are negative!")
        return -1.0

    mu = 1.0 - neg_sum / pos_sum
    return mu


def find_optimal_M_mu(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    extend_decades: float = 0.0,
    include_Rs: bool = True,
    include_L: bool = True,
    fit_type: str = 'complex',
    allow_negative: bool = True,
    weighting: str = 'modulus'
) -> Tuple[int, float, NDArray[np.float64], NDArray[np.float64], Optional[float]]:
    """
    Find optimal number of Voigt elements using mu metric (Lin-KK style).

    Iteratively increases M from 3 until mu < mu_threshold or M >= max_M.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85, as in Lin-KK)
        Lower values -> more conservative (fewer elements)
    max_M : int, optional
        Maximum number of elements to try (default: 50)
    extend_decades : float, optional
        Extend tau range toward lower frequencies (default: 0.0 for Lin-KK)
    include_Rs : bool, optional
        Include series resistance R_s (default: True)
    include_L : bool, optional
        Include series inductance L (default: True for Lin-KK)
    fit_type : str, optional
        Fit type: 'real', 'imag', or 'complex' (default: 'complex')
    allow_negative : bool, optional
        Allow negative R_i values (default: True for Lin-KK compatibility)
        Note: mu metric is designed for pseudoinverse (allow_negative=True).
        With NNLS (allow_negative=False), all R_i >= 0, so mu ~ 1 always.
    weighting : str, optional
        Point weighting scheme (default: 'modulus' = Lin-KK standard)

    Returns
    -------
    M_optimal : int
        Optimal number of Voigt elements
    mu_final : float
        Final mu value at M_optimal
    tau_optimal : ndarray of float
        Optimal tau grid
    elements : ndarray of float
        Optimal element values [R_s, R_1, ..., R_M, L] or subset
    L_value : float or None
        Estimated inductance [H] if include_L=True

    Notes
    -----
    Algorithm from Schonleber et al. (2014):
    1. Start with M = 3
    2. Generate M time constants logarithmically
    3. Fit R_i using pseudoinverse (allow_negative=True)
    4. Calculate mu metric
    5. If mu > threshold, increase M by 1 and repeat
    6. If mu <= threshold, STOP (optimal M found)

    References
    ----------
    Schonleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    # Validate inputs
    validate_eis_data(frequencies, Z, context="find_optimal_M_mu")

    weighting_labels = {
        'uniform': 'uniform (w=1)',
        'sqrt': 'sqrt (w=1/sqrt|Z|)',
        'modulus': 'modulus (w=1/|Z|)',
        'proportional': 'proportional (w=1/|Z|^2)'
    }
    logger.info(f"  mu threshold: {mu_threshold}")
    logger.info(f"  Max M: {max_M}")
    logger.info(f"  Fit type: {fit_type}")
    logger.info(f"  Weighting: {weighting_labels.get(weighting, weighting)}")
    logger.info(f"  Include L: {include_L}")
    logger.info(f"  Allow negative R_i: {allow_negative}")
    if not allow_negative:
        logger.warning("  NOTE: mu metric is designed for allow_negative=True (Lin-KK)")
        logger.warning("  With NNLS, all R_i >= 0, so mu ~ 1 always")
    logger.info("")

    M = 2  # Start with M=3 (Lin-KK standard)
    mu = 1.0
    iteration = 0
    L_value = None
    tau = None
    elements = None
    R_i = np.array([])

    while mu > mu_threshold and M < max_M:
        M += 1
        iteration += 1

        # Generate tau grid for this M
        tau = generate_tau_grid_fixed_M(frequencies, M, extend_decades)

        # Fit using specified method
        elements, residual, L_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=include_Rs,
            include_L=include_L,
            fit_type=fit_type,
            allow_negative=allow_negative,
            weighting=weighting
        )

        # Extract R_i for mu calculation (exclude R_s and L)
        R_start = 1 if include_Rs else 0
        R_end = -1 if include_L else len(elements)
        R_i = elements[R_start:R_end]

        # Calculate mu
        mu = calc_mu(R_i)

        # Log every 5 iterations or when converged
        if M % 5 == 0 or mu <= mu_threshold or M == 3:
            n_negative = np.sum(R_i < 0)
            logger.info(f"  Iter {iteration:2d}: M={M:2d}, mu={mu:.4f}, "
                       f"residual={residual:.3e}, negative R_i={n_negative}/{len(R_i)}")

    # Final result
    logger.info("")
    if mu <= mu_threshold:
        logger.info(f"Optimal M found: M = {M}")
        logger.info(f"  mu = {mu:.4f} <= {mu_threshold}")
    else:
        logger.warning(f"Reached max_M = {max_M}")
        logger.warning(f"  mu = {mu:.4f} > {mu_threshold}")
        logger.warning("  Model may still be overfit!")

    # Count negative R_i
    n_negative = np.sum(R_i < 0)
    if n_negative > 0:
        logger.info(f"  Negative R_i: {n_negative}/{len(R_i)} "
                   f"({n_negative/len(R_i)*100:.1f}%)")

    logger.info("="*60)

    return M, mu, tau, elements, L_value


__all__ = ['calc_mu', 'find_optimal_M_mu']
