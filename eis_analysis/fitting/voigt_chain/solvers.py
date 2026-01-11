"""
Linear algebra solvers for Voigt chain fitting.

This module provides the NNLS solver and design matrix computation.
"""

import numpy as np
import logging
from typing import Tuple
from numpy.typing import NDArray
from scipy.optimize import nnls, lsq_linear

logger = logging.getLogger(__name__)


def robust_nnls(
    A: NDArray[np.float64],
    b: NDArray[np.float64],
    max_iter: int = None
) -> Tuple[NDArray[np.float64], float]:
    """
    Robust non-negative least squares solver with multiple fallback strategies.

    Solves: min ||A @ x - b||^2 subject to x >= 0

    Strategy:
    1. Try scipy.optimize.nnls (Lawson-Hanson algorithm)
    2. If fails, try scipy.optimize.lsq_linear with bounds (proper bounded LS)
    3. If both fail, use pseudoinverse + clipping (last resort, with warning)

    Parameters
    ----------
    A : ndarray of float, shape (m, n)
        Design matrix
    b : ndarray of float, shape (m,)
        Target vector
    max_iter : int, optional
        Maximum iterations for nnls (default: 3 * n)

    Returns
    -------
    x : ndarray of float, shape (n,)
        Non-negative solution
    residual_norm : float
        Norm of residual ||A @ x - b||

    Notes
    -----
    Unlike simple clipping (np.maximum(x, 0)), this function properly solves
    the constrained optimization problem. Clipping can give suboptimal or
    even incorrect solutions because it doesn't account for the constraint
    during optimization.

    Example of why clipping is wrong:
        A = [[1, 1], [1, -1]]
        b = [1, 2]
        Unconstrained solution: x = [1.5, -0.5]
        Clipped: x = [1.5, 0] -> residual = 1.12
        True NNLS: x = [1, 0] -> residual = 1.0 (better!)
    """
    n_vars = A.shape[1]

    if max_iter is None:
        max_iter = 3 * n_vars

    # Strategy 1: scipy.optimize.nnls (Lawson-Hanson)
    try:
        x, residual_norm = nnls(A, b, maxiter=max_iter)
        return x, residual_norm
    except (RuntimeError, ValueError) as e:
        logger.debug(f"nnls failed: {e}, trying lsq_linear")

    # Strategy 2: scipy.optimize.lsq_linear with bounds
    # This is a proper bounded least squares solver
    try:
        result = lsq_linear(
            A, b,
            bounds=(0, np.inf),  # x >= 0
            method='bvls',  # Bounded-variable least squares
            max_iter=max_iter * 10,  # More iterations for complex problems
            tol=1e-10
        )

        if result.success or result.status == 1:  # 1 = max_iter reached but OK
            x = result.x
            residual_norm = np.linalg.norm(A @ x - b)
            return x, residual_norm
        else:
            logger.debug(f"lsq_linear failed: {result.message}")
    except Exception as e:
        logger.debug(f"lsq_linear exception: {e}")

    # Strategy 3: Last resort - pseudoinverse + clipping (with warning)
    logger.warning(
        "Both nnls and lsq_linear failed. Using pseudoinverse with clipping. "
        "Results may be suboptimal. Consider reducing number of parameters."
    )

    x_unconstrained = np.linalg.lstsq(A, b, rcond=None)[0]
    x = np.maximum(x_unconstrained, 0)

    # Re-solve for optimal scaling of non-zero components
    # This improves the clipped solution
    nonzero_mask = x > 0
    if np.any(nonzero_mask):
        A_reduced = A[:, nonzero_mask]
        x_reduced, _ = np.linalg.lstsq(A_reduced, b, rcond=None)[:2]
        x_reduced = np.maximum(x_reduced, 0)  # Ensure still non-negative
        x[nonzero_mask] = x_reduced

    residual_norm = np.linalg.norm(A @ x - b)
    return x, residual_norm


def compute_voigt_matrix(
    frequencies: NDArray[np.float64],
    tau: NDArray[np.float64],
    include_Rs: bool = True
) -> NDArray[np.float64]:
    """
    Compute the design matrix A for linear least squares problem.

    For fixed tau values, the real part of Voigt chain impedance is:
        Z'(omega) = R_s + sum R_i * h_i(omega, tau_i)

    where h_i(omega, tau_i) = tau_i^2 * omega^2 / (1 + tau_i^2 * omega^2)

    This is LINEAR in [R_s, R_1, R_2, ..., R_N], so we can write:
        Z' = A @ [R_s, R_1, ..., R_N]

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz] (M points)
    tau : ndarray of float
        Time constants [s] (N points)
    include_Rs : bool, optional
        If True, include column for R_s (series resistance) (default: True)

    Returns
    -------
    A : ndarray of float
        Design matrix, shape (M, N) or (M, N+1) if include_Rs=True
        A[k, 0] = 1 (for R_s) if include_Rs
        A[k, i] = 1 / (1 + (omega_k * tau_i)^2)

    Notes
    -----
    The function h_i(omega, tau_i) is the contribution of the i-th Voigt element
    to the real part of impedance, normalized by R_i.

    For omega -> 0: h_i -> 1 (capacitor blocks DC, full R_i contributes)
    For omega -> inf: h_i -> 0 (capacitor shorts, R_i is bypassed)
    At omega = 1/tau_i: h_i = 0.5 (half-power point)

    Examples
    --------
    >>> freq = np.array([1.0, 10.0, 100.0])
    >>> tau = np.array([0.01, 0.1, 1.0])
    >>> A = compute_voigt_matrix(freq, tau, include_Rs=True)
    >>> A.shape
    (3, 4)  # 3 frequencies, 3 tau + 1 R_s
    """
    omega = 2 * np.pi * frequencies  # Angular frequency [rad/s]
    n_freq = len(frequencies)
    n_tau = len(tau)

    # Determine matrix size
    n_cols = n_tau + 1 if include_Rs else n_tau

    # Initialize matrix
    A = np.zeros((n_freq, n_cols))

    # Column index offset
    col_offset = 1 if include_Rs else 0

    # First column: R_s contribution (constant offset)
    if include_Rs:
        A[:, 0] = 1.0

    # Remaining columns: Voigt elements (vectorized using broadcasting)
    # For Voigt element: Z_i = R_i / (1 + j*omega*tau_i)
    # Real part: Z'_i = R_i / (1 + (omega*tau_i)^2)
    # So: h_i(omega, tau_i) = 1 / (1 + (omega*tau_i)^2)
    # Broadcasting: tau (n_tau,) x omega (n_freq,) -> (n_tau, n_freq)
    tau_omega_sq = (tau[:, np.newaxis] * omega[np.newaxis, :]) ** 2
    A[:, col_offset:col_offset + n_tau] = (1.0 / (1 + tau_omega_sq)).T

    return A


__all__ = ['robust_nnls', 'compute_voigt_matrix']
