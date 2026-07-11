"""
Covariance matrix computation for circuit fitting.

Provides robust SVD-based covariance estimation and confidence interval
computation for fitted parameters.

Author: EIS Analysis Toolkit
"""

import numpy as np
import logging
from typing import List, Optional, Union, Tuple
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
        Condition number of the column-scaled normal matrix
        J_s^T @ J_s = cond(J_s)^2, where J_s has unit-norm columns
        (high = ill-conditioned). Scale-invariant: unit disparity between
        parameters (R ~ 1e5 vs Q ~ 1e-6) does not inflate it; only genuine
        parameter correlation does. Governs reliability of the (J^T J)^{-1}
        used for the covariance.
    rank : int
        Numerical rank of the column-scaled Jacobian (scale-invariant)
    is_well_conditioned : bool
        True if cond(J_s^T J_s) < 1e10
    warning_message : str or None
        Warning message if any issues detected
    dof : int
        Residual degrees of freedom of the variance estimate
        (n_residuals - n_free_params). Used for the t-distribution CI so it
        matches the dof used for s^2.
    """
    cov: Union[NDArray[np.float64], None]
    stderr: NDArray[np.float64]
    condition_number: float
    rank: int
    is_well_conditioned: bool
    warning_message: Union[str, None]
    dof: int = 1


def compute_covariance_matrix(
    jacobian: NDArray[np.float64],
    residuals: NDArray[np.float64],
    n_params: int,
    fixed_params: Optional[List[bool]] = None,
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

    If the Jacobian is rank-deficient (any singular value < rcond * max(S)),
    J^T J is singular and the covariance does not exist. In that case the
    covariance and standard errors of the affected (free) parameters are
    reported as infinite -- the scipy.optimize.curve_fit convention -- rather
    than substituting an arbitrary regularized value.

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

    For ill-conditioned problems (cond(J^T J) > 1e10), the covariance
    estimate may be unreliable; for rank-deficient problems it is reported as
    infinite (not estimable). This often indicates:
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
        # Column-scale the Jacobian before analysing rank and conditioning.
        # EIS parameters span many orders of magnitude (R ~ 1e5 Ohm, Q ~ 1e-6,
        # n ~ 0.5), so the raw Jacobian columns differ by >10 decades and
        # cond(J) exceeds 1e10 purely from those units -- a scaling artifact,
        # not genuine non-identifiability. Scaling each column to unit norm
        # (van der Sluis preconditioning) makes the rank and condition number
        # scale-invariant, so only true rank deficiency (linearly dependent
        # columns) is flagged. The covariance is computed in the scaled space
        # and rescaled: with D = diag(col_norms), J = J_scaled @ D and
        # (J^T J)^{-1} = D^{-1} (J_scaled^T J_scaled)^{-1} D^{-1}.
        col_norms = np.linalg.norm(jacobian, axis=0)
        # A zero-norm column means the parameter has no effect on the residual
        # (genuinely non-identifiable). Use unit scale to avoid div-by-zero;
        # the resulting zero singular value drops the rank -> inf below.
        col_scale = np.where(col_norms > 0, col_norms, 1.0)
        J_scaled = jacobian / col_scale

        # SVD of the scaled Jacobian: J_scaled = U @ S @ V^T
        U, S, Vt = np.linalg.svd(J_scaled, full_matrices=False)

        # Condition number of the scaled normal matrix J_scaled^T J_scaled =
        # cond(J_scaled)^2. The covariance is s^2 (J^T J)^{-1}; its reliability
        # in terms of relative parameter uncertainties is governed by this
        # scale-invariant condition number, not the unit-dependent raw cond(J).
        cond_J = S[0] / S[-1] if S[-1] > 0 else np.inf
        condition_number = cond_J ** 2

        # Numerical rank (scaled singular values above threshold)
        threshold = rcond * S[0]
        rank = np.sum(S > threshold)

        # Check conditioning
        is_well_conditioned = condition_number < 1e10

        # Rank-deficient: J^T J is singular, so the covariance does not exist.
        # Report it as infinite (the scipy.optimize.curve_fit convention) rather
        # than substituting an arbitrary regularized value: the honest statement
        # is that the uncertainty is not estimable for the affected parameters.
        # Fixed parameters keep zero variance (they are known exactly).
        if rank < n_free_params:
            warning_message = (
                f"Rank-deficient Jacobian (rank={rank}/{n_free_params}). "
                f"Some parameters are not identifiable; covariance not estimable."
            )
            if fixed_params is not None and any(fixed_params):
                n_total = len(fixed_params)
                free_indices = [i for i, f in enumerate(fixed_params) if not f]
            else:
                n_total = n_free_params
                free_indices = list(range(n_free_params))

            cov = np.zeros((n_total, n_total))
            for i_full in free_indices:
                for j_full in free_indices:
                    cov[i_full, j_full] = np.inf
            stderr = np.zeros(n_total)
            stderr[free_indices] = np.inf

            return CovarianceResult(
                cov=cov,
                stderr=stderr,
                condition_number=condition_number,
                rank=int(rank),
                is_well_conditioned=False,
                warning_message=warning_message,
                dof=dof,
            )

        if not is_well_conditioned:
            warning_message = (
                f"Ill-conditioned normal matrix J^T J (cond={condition_number:.2e}). "
                f"Covariance estimates may be unreliable."
            )

        # Full rank here: every scaled singular value is above the threshold.
        # Invert the scaled normal matrix, (J_scaled^T J_scaled)^{-1} =
        # V @ S^{-2} @ V^T, then rescale by D^{-1} on both sides to recover
        # (J^T J)^{-1} in the original parameter units.
        S_inv_sq = 1.0 / (S * S)
        V = Vt.T
        JsTJs_inv = V @ np.diag(S_inv_sq) @ Vt
        inv_scale = 1.0 / col_scale
        JtJ_inv_free = (inv_scale[:, None] * JsTJs_inv) * inv_scale[None, :]

        # Covariance matrix for free parameters
        cov_free = residual_variance * JtJ_inv_free

        # Check for negative diagonal (numerical artifact; the stderr below
        # already uses abs). Append rather than assign so an earlier
        # ill-conditioned warning is not silently overwritten.
        if np.any(np.diag(cov_free) < 0):
            neg_msg = "Negative variance detected. Using absolute values."
            warning_message = (f"{warning_message} {neg_msg}"
                               if warning_message else neg_msg)

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
            warning_message=warning_message,
            dof=dof,
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
            warning_message=f"SVD failed: {e}",
            dof=dof,
        )


def compute_confidence_interval(
    params_opt: NDArray[np.float64],
    params_stderr: NDArray[np.float64],
    dof: int,
    confidence_level: float = 0.95,
    log_scale: Optional[List[bool]] = None
) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Compute confidence intervals for parameters using t-distribution.

    Parameters
    ----------
    params_opt : ndarray
        Optimal parameters
    params_stderr : ndarray
        Standard errors of parameters (from covariance matrix)
    dof : int
        Residual degrees of freedom of the variance estimate
        (n_residuals - n_free_params). Must match the dof used for s^2 in
        compute_covariance_matrix; available as CovarianceResult.dof.
    confidence_level : float, optional
        Confidence level (0.95 for 95% CI, 0.99 for 99% CI), default 0.95
    log_scale : list of bool, optional
        Per-parameter mask: True for positive scale parameters (R, C, Q,
        tau, ...) whose CI is computed in log space (see Notes). None
        (default) keeps the symmetric linear-space CI for all parameters.
        Use `bounds.log_scale_ci_mask` to build the mask from bounds.

    Returns
    -------
    ci_low : ndarray
        Lower bounds of confidence intervals
    ci_high : ndarray
        Upper bounds of confidence intervals

    Notes
    -----
    Uses the t-distribution with `dof` degrees of freedom. For large dof
    (> 30) t approximates the normal distribution (1.96 for 95%); for small
    dof it is more conservative (wider CI).

    Linear (default): ci = p +/- t * se. For a positive scale parameter with
    se comparable to p this produces a physically meaningless negative lower
    bound (e.g. a resistance CI containing negative values).

    Log-space (masked parameters): delta method on ln(p) gives
    se_ln = se / p, so ci = (p / f, p * f) with f = exp(t * se / p).
    Always positive and multiplicatively symmetric; for se/p -> 0 it
    converges to the linear interval. Applied only where the mask is True,
    the parameter is positive, and se is finite.
    """
    dof = max(int(dof), 1)

    alpha = 1 - confidence_level
    t_critical = t.ppf(1 - alpha/2, dof)

    params_opt = np.asarray(params_opt, dtype=float)
    params_stderr = np.asarray(params_stderr, dtype=float)

    margin = t_critical * params_stderr
    ci_low = params_opt - margin
    ci_high = params_opt + margin

    if log_scale is not None:
        mask = (np.asarray(log_scale, dtype=bool)
                & (params_opt > 0)
                & np.isfinite(params_stderr))
        if np.any(mask):
            p = params_opt[mask]
            factor = np.exp(t_critical * params_stderr[mask] / p)
            ci_low[mask] = p / factor
            ci_high[mask] = p * factor

    return ci_low, ci_high


__all__ = [
    'CovarianceResult',
    'compute_covariance_matrix',
    'compute_confidence_interval',
]
