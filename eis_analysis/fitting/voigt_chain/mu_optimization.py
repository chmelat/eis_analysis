"""
Mu metric optimization for Voigt chain fitting (Lin-KK style).

This module provides the mu metric calculation and automatic optimization
of the number of Voigt elements.
"""

import numpy as np
import logging
from dataclasses import dataclass, field
from typing import List, Optional
from numpy.typing import NDArray

from .validation import validate_eis_data
from .tau_grid import generate_tau_grid_fixed_M
from .fitting import estimate_R_linear

logger = logging.getLogger(__name__)


@dataclass
class MuIteration:
    """One step of the M search."""
    iteration: int
    M: int
    mu: float
    residual: float
    n_negative: int   # R_i below zero at this M
    n_R: int          # R_i in total at this M


@dataclass
class MuOptimization:
    """
    Result of the mu-metric search for the number of Voigt elements.

    The search is a library computation, so it reports rather than prints:
    Kramers-Kronig calls it too, and its progress has no business appearing
    inside the KK section of the CLI output.

    Attributes
    ----------
    M : int
        Number of Voigt elements the search stopped at
    mu : float
        Mu value at that M (below mu_threshold on normal termination,
        above it only when max_M was reached)
    tau : ndarray of float
        Time-constant grid for M
    elements : ndarray of float
        Element values [R_s, R_1, ..., R_M, L] or subset
    L_value : float or None
        Estimated inductance [H] if include_L=True
    C_value : float or None
        Estimated series capacitance [F] if include_C=True
        (not part of the elements array)
    reached_max_M : bool
        True when the search ran out of M rather than converging
    n_negative, n_R : int
        Negative and total R_i at the stopping M
    iterations : list of MuIteration
        Sampled progress of the search, for callers that report it
    warnings : list of str
        Caveats about the result
    """
    M: int
    mu: float
    tau: NDArray[np.float64]
    elements: NDArray[np.float64]
    L_value: Optional[float] = None
    C_value: Optional[float] = None
    reached_max_M: bool = False
    n_negative: int = 0
    n_R: int = 0
    iterations: List[MuIteration] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


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
    include_C: bool = False,
    fit_type: str = 'complex',
    allow_negative: bool = True,
    weighting: str = 'modulus'
) -> MuOptimization:
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
        Lower values -> the iteration stops later -> more elements
        (higher overfit tolerance); higher values -> fewer elements
    max_M : int, optional
        Maximum number of elements to try (default: 50)
    extend_decades : float, optional
        Extend tau range toward lower frequencies (default: 0.0 for Lin-KK)
    include_Rs : bool, optional
        Include series resistance R_s (default: True)
    include_L : bool, optional
        Include series inductance L (default: True for Lin-KK)
    include_C : bool, optional
        Include series capacitance C for blocking low-frequency behavior
        (default: False; Schonleber Lin-KK 'add_cap')
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
    MuOptimization
        The stopping M and its fit, plus the progress of the search and any
        caveat about it. Nothing is logged: the caller decides what, if
        anything, to report.

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

    warnings: List[str] = []
    if not allow_negative:
        warnings.append("mu metric is designed for allow_negative=True (Lin-KK)")
        warnings.append("With NNLS, all R_i >= 0, so mu ~ 1 always")

    iterations: List[MuIteration] = []

    M = 2  # Start with M=3 (Lin-KK standard)
    mu = 1.0
    iteration = 0
    L_value = None
    C_value = None
    tau = None
    elements = None
    R_i = np.array([])

    while mu > mu_threshold and M < max_M:
        M += 1
        iteration += 1

        # Generate tau grid for this M
        tau = generate_tau_grid_fixed_M(frequencies, M, extend_decades)

        # Fit using specified method
        elements, residual, L_value, C_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=include_Rs,
            include_L=include_L,
            include_C=include_C,
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

        # Sample the search every 5 steps, plus the first and the last one
        if M % 5 == 0 or mu <= mu_threshold or M == 3:
            iterations.append(MuIteration(
                iteration=iteration, M=M, mu=mu, residual=residual,
                n_negative=int(np.sum(R_i < 0)), n_R=len(R_i)))

    reached_max_M = mu > mu_threshold
    if reached_max_M:
        warnings.append(f"Reached max_M = {max_M}, mu = {mu:.4f} > {mu_threshold}; "
                        f"model may still be overfit")

    return MuOptimization(
        M=M, mu=mu, tau=tau, elements=elements,
        L_value=L_value, C_value=C_value,
        reached_max_M=reached_max_M,
        n_negative=int(np.sum(R_i < 0)), n_R=len(R_i),
        iterations=iterations, warnings=warnings)


__all__ = ['calc_mu', 'find_optimal_M_mu', 'MuOptimization', 'MuIteration']
