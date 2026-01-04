"""
Covariance matrix computation for circuit fitting.

Provides robust SVD-based covariance estimation and confidence interval
computation for fitted parameters.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import List, Union, Tuple
from numpy.typing import NDArray
from dataclasses import dataclass
from scipy.stats import t

logger = logging.getLogger(__name__)


@dataclass
class CovarianceResult:
    """
    Result of covariance matrix computation.

    Attributes
    ----------
    cov : ndarray or None
        Covariance matrix (None if computation failed)
    stderr : ndarray
        Standard errors of parameters (inf if computation failed)
    condition_number : float
        Condition number of J^T @ J (high = ill-conditioned)
    rank : int
        Numerical rank of Jacobian
    is_well_conditioned : bool
        True if condition number < 1e10
    warning_message : str or None
        Warning message if any issues detected
    """
    cov: Union[NDArray[np.float64], None]
    stderr: NDArray[np.float64]
    condition_number: float
    rank: int
    is_well_conditioned: bool
    warning_message: Union[str, None]


def compute_covariance_matrix(
    jacobian: NDArray[np.float64],
    residuals: NDArray[np.float64],
    n_params: int,
    fixed_params: List[bool] = None,
    rcond: float = 1e-10
) -> CovarianceResult:
    """
    Compute covariance matrix from Jacobian using robust SVD-based method.

    For weighted least squares, the Jacobian should already include weights
    (i.e., J_weighted[i,j] = w_i * dZ_i/dp_j). This is automatically the case
    when weights are applied to residuals and least_squares computes the Jacobian.

    Mathematical background
    -----------------------
    For weighted least squares with residuals r_i = w_i * (Z_i - Z_model_i):

        J_w = W @ J_original   (weighted Jacobian)
        RSS = r^T @ r          (sum of squared weighted residuals)
        s^2 = RSS / (n - p)    (residual variance estimate)
        cov(theta) = s^2 * (J_w^T @ J_w)^{-1}

    This function uses SVD for numerical stability:

        J_w = U @ S @ V^T
        J_w^T @ J_w = V @ S^2 @ V^T
        (J_w^T @ J_w)^{-1} = V @ S^{-2} @ V^T

    Small singular values (< rcond * max(S)) are regularized to prevent
    numerical instability.

    Parameters
    ----------
    jacobian : ndarray of float, shape (n_residuals, n_free_params)
        Jacobian matrix (already weighted if using weighted LS)
    residuals : ndarray of float, shape (n_residuals,)
        Residual vector (already weighted if using weighted LS)
    n_params : int
        Total number of parameters (including fixed)
    fixed_params : list of bool, optional
        Which parameters are fixed (True = fixed, False = free)
        If None, all parameters are assumed free
    rcond : float, optional
        Cutoff for small singular values (default: 1e-10)

    Returns
    -------
    result : CovarianceResult
        Covariance computation result with diagnostics

    Notes
    -----
    The covariance matrix is computed in the space of free parameters,
    then expanded to include fixed parameters (with zero variance/covariance).

    For ill-conditioned problems (condition_number > 1e10), the covariance
    estimate may be unreliable. This often indicates:
    - Over-parametrized model
    - Correlated parameters
    - Insufficient data in some frequency range

    References
    ----------
    Bard, Y. "Nonlinear Parameter Estimation" Academic Press, 1974.
    """
    n_residuals = len(residuals)
    n_free_params = jacobian.shape[1]

    # Degrees of freedom
    dof = max(n_residuals - n_free_params, 1)

    # Residual variance estimate
    rss = residuals @ residuals
    residual_variance = rss / dof

    warning_message = None

    try:
        # SVD decomposition: J = U @ S @ V^T
        U, S, Vt = np.linalg.svd(jacobian, full_matrices=False)

        # Condition number
        condition_number = S[0] / S[-1] if S[-1] > 0 else np.inf

        # Numerical rank (singular values above threshold)
        threshold = rcond * S[0]
        rank = np.sum(S > threshold)

        # Check conditioning
        is_well_conditioned = condition_number < 1e10

        if not is_well_conditioned:
            warning_message = (
                f"Ill-conditioned Jacobian (cond={condition_number:.2e}). "
                f"Covariance estimates may be unreliable."
            )

        if rank < n_free_params:
            warning_message = (
                f"Rank-deficient Jacobian (rank={rank}/{n_free_params}). "
                f"Some parameters are not identifiable from data."
            )

        # Compute (J^T @ J)^{-1} using SVD
        S_inv_sq = np.zeros_like(S)
        for i, s in enumerate(S):
            if s > threshold:
                S_inv_sq[i] = 1.0 / (s * s)
            else:
                S_inv_sq[i] = 1.0 / (threshold * threshold)

        V = Vt.T
        JtJ_inv_free = V @ np.diag(S_inv_sq) @ Vt

        # Covariance matrix for free parameters
        cov_free = residual_variance * JtJ_inv_free

        # Check for negative diagonal
        diag_cov = np.diag(cov_free)
        if np.any(diag_cov < 0):
            warning_message = "Negative variance detected. Using absolute values."
            diag_cov = np.abs(diag_cov)

        # Expand to full parameter space (including fixed parameters)
        if fixed_params is not None and any(fixed_params):
            n_total = len(fixed_params)
            cov = np.zeros((n_total, n_total))

            free_indices = [i for i, is_fixed in enumerate(fixed_params) if not is_fixed]

            for i_free, i_full in enumerate(free_indices):
                for j_free, j_full in enumerate(free_indices):
                    cov[i_full, j_full] = cov_free[i_free, j_free]

            stderr = np.zeros(n_total)
            for i_free, i_full in enumerate(free_indices):
                stderr[i_full] = np.sqrt(np.abs(cov_free[i_free, i_free]))
        else:
            cov = cov_free
            stderr = np.sqrt(np.abs(np.diag(cov)))

        return CovarianceResult(
            cov=cov,
            stderr=stderr,
            condition_number=condition_number,
            rank=rank,
            is_well_conditioned=is_well_conditioned,
            warning_message=warning_message
        )

    except np.linalg.LinAlgError as e:
        logger.error(f"SVD computation failed: {e}")
        n_total = n_params
        if fixed_params is not None:
            n_total = len(fixed_params)

        return CovarianceResult(
            cov=None,
            stderr=np.full(n_total, np.inf),
            condition_number=np.inf,
            rank=0,
            is_well_conditioned=False,
            warning_message=f"SVD failed: {e}"
        )


def compute_confidence_interval(
    params_opt: NDArray[np.float64],
    params_stderr: NDArray[np.float64],
    n_data: int,
    confidence_level: float = 0.95
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Compute confidence intervals for parameters using t-distribution.

    Parameters
    ----------
    params_opt : ndarray
        Optimal parameters
    params_stderr : ndarray
        Standard errors of parameters (from covariance matrix)
    n_data : int
        Number of data points (for degrees of freedom)
    confidence_level : float, optional
        Confidence level (0.95 for 95% CI, 0.99 for 99% CI), default 0.95

    Returns
    -------
    ci_low : ndarray
        Lower bounds of confidence intervals
    ci_high : ndarray
        Upper bounds of confidence intervals

    Notes
    -----
    Uses t-distribution with (n_data - n_params) degrees of freedom.
    For large datasets (n > 30), t approximates normal distribution (1.96 for 95%).
    For small datasets, t-distribution is more conservative (wider CI).
    """
    n_params = len(params_opt)
    dof = max(n_data - n_params, 1)

    alpha = 1 - confidence_level
    t_critical = t.ppf(1 - alpha/2, dof)

    margin = t_critical * params_stderr
    ci_low = params_opt - margin
    ci_high = params_opt + margin

    return ci_low, ci_high


__all__ = [
    'CovarianceResult',
    'compute_covariance_matrix',
    'compute_confidence_interval',
]
