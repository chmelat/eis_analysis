"""
Regularized linear system for DRT analysis.

Assemble and solve the Tikhonov-regularized DRT problem: validate the frequency
grid, build the system matrices A, b and regularization operator L, select the
regularization parameter lambda, and solve the non-negative least squares system.
"""

import numpy as np
import logging
from typing import Optional, List
from numpy.typing import NDArray
from scipy.optimize import nnls

from .results import DRTMatrices, LambdaSelection, NNLSSolution
from .gcv import find_optimal_lambda_gcv, find_optimal_lambda_hybrid

logger = logging.getLogger(__name__)

# Constants
MIN_FREQUENCY_RANGE = 10
GAMMA_MAX_REASONABLE = 1e10


def _validate_frequencies(frequencies: NDArray) -> List[str]:
    """
    Validate frequency array for DRT analysis.

    Returns list of warnings. Raises ValueError for critical errors.
    """
    warnings = []

    if np.any(frequencies <= 0):
        raise ValueError("All frequencies must be positive for DRT analysis")

    f_max = frequencies.max()
    f_min = frequencies.min()

    if f_max <= 0 or f_min <= 0:
        raise ValueError(f"Invalid frequency range: min={f_min}, max={f_max}")

    if f_max <= f_min:
        raise ValueError(f"Max frequency ({f_max:.2e} Hz) must be > min ({f_min:.2e} Hz)")

    freq_range = f_max / f_min
    if freq_range < MIN_FREQUENCY_RANGE:
        warnings.append(f"Small frequency range: {freq_range:.1f}x (recommended >{MIN_FREQUENCY_RANGE}x)")

    return warnings


def _build_drt_matrices(frequencies: NDArray, Z: NDArray,
                        R_inf: float, n_tau: int = 100) -> DRTMatrices:
    """
    Build DRT system matrices A, b, and regularization matrix L.
    """
    f_max = frequencies.max()
    f_min = frequencies.min()

    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    omega = 2 * np.pi * frequencies

    d_ln_tau_array = np.diff(np.log(tau))
    d_ln_tau = float(np.mean(d_ln_tau_array))

    # Vectorized matrix construction
    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom

    A = np.vstack([A_re, A_im])
    cond_A = float(np.linalg.cond(A))

    b = np.concatenate([Z.real - R_inf, Z.imag])

    # Regularization matrix (2nd derivative) - tridiagonal [1, -2, 1]
    # Shape: (n_tau - 2, n_tau) for second derivative operator
    L = np.zeros((n_tau - 2, n_tau))
    np.fill_diagonal(L, 1)           # Main diagonal at offset 0
    np.fill_diagonal(L[:, 1:], -2)   # Diagonal at offset 1
    np.fill_diagonal(L[:, 2:], 1)    # Diagonal at offset 2

    return DRTMatrices(
        A=A, A_re=A_re, A_im=A_im, b=b, L=L,
        tau=tau, d_ln_tau=d_ln_tau, condition_number=cond_A
    )


def _select_lambda(A: NDArray, b: NDArray, L: NDArray,
                   lambda_reg: Optional[float] = None,
                   auto_lambda: bool = False) -> LambdaSelection:
    """
    Select regularization parameter lambda.
    """
    # GCV search bounds. Edge detection (F3/F7): lambda landing at a bound -
    # or the GCV guess pinning there even when L-curve corrected it - signals
    # the optimizer wants more extreme regularization than the range allows.
    lambda_range = (1e-5, 1.0)

    def _at_bound(lam: Optional[float]) -> bool:
        return lam is not None and (lam <= lambda_range[0] or lam >= lambda_range[1])

    if auto_lambda:
        try:
            lambda_opt, gcv_score, diag = find_optimal_lambda_hybrid(
                A, b, L,
                lambda_range=lambda_range,
                n_search=20,
                lcurve_decades=1.5
            )
            lambda_gcv = diag.get('lambda_gcv')
            corner_at_edge = diag.get('corner_at_edge', False)
            at_edge = _at_bound(lambda_opt) or _at_bound(lambda_gcv) or corner_at_edge
            method = 'hybrid' if diag['method_used'] == 'lcurve_correction' else 'gcv'
            return LambdaSelection(
                lambda_value=lambda_opt,
                method=method,
                # lambda_gcv only meaningful (and displayed) for L-curve correction
                lambda_gcv=lambda_gcv if method == 'hybrid' else None,
                gcv_score=gcv_score,
                corner_at_edge=corner_at_edge,
                lambda_at_edge=at_edge
            )
        except (np.linalg.LinAlgError, ValueError):
            logger.debug("Hybrid lambda selection failed, falling back to GCV",
                         exc_info=True)
            try:
                lambda_opt, gcv_score = find_optimal_lambda_gcv(
                    A, b, L, lambda_range=lambda_range, n_search=20
                )
                return LambdaSelection(
                    lambda_value=lambda_opt,
                    method='gcv',
                    gcv_score=gcv_score,
                    lambda_at_edge=_at_bound(lambda_opt)
                )
            except (np.linalg.LinAlgError, ValueError):
                logger.debug("GCV lambda selection failed, using fallback "
                             "lambda=0.1", exc_info=True)
                return LambdaSelection(lambda_value=0.1, method='fallback')

    if lambda_reg is None:
        return LambdaSelection(lambda_value=0.1, method='default')

    return LambdaSelection(lambda_value=lambda_reg, method='user')


def _solve_nnls(A: NDArray, b: NDArray, L: NDArray,
                lambda_reg: float, n_tau: int,
                Z: NDArray) -> NNLSSolution:
    """
    Solve regularized NNLS problem.
    """
    warnings = []

    # Check for inductive data
    n_inductive = int(np.sum(Z.imag > 0))
    inductive_fraction = n_inductive / len(Z) if len(Z) > 0 else 0.0
    max_inductive = float(np.max(Z.imag)) if n_inductive > 0 else 0.0

    if inductive_fraction > 0.1 or max_inductive > 50:
        warnings.append(
            f"{n_inductive} points with inductive component ({inductive_fraction*100:.1f}%), "
            f"max Z'' = {max_inductive:.2f} Ohm"
        )

    # Build regularized system
    A_reg = np.vstack([A, np.sqrt(lambda_reg) * L])
    b_reg = np.concatenate([b, np.zeros(n_tau - 2)])

    # Condition number of the system actually solved. The bare kernel A is
    # intrinsically ill-conditioned for any DRT problem (that is why Tikhonov
    # regularization is applied); the regularized system is what determines the
    # numerical health of the solve.
    cond_reg = float(np.linalg.cond(A_reg))

    # Solve NNLS
    try:
        gamma, residual = nnls(A_reg, b_reg)
    except Exception as e:
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=[f"NNLS solver error: {e}"]
        )

    # Validate solution
    if np.any(~np.isfinite(gamma)):
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=["NNLS returned NaN or Inf values"]
        )

    gamma_max = float(np.max(gamma))
    gamma_nonzero = gamma[gamma > 0]
    gamma_min_nonzero = float(np.min(gamma_nonzero)) if len(gamma_nonzero) > 0 else None

    if gamma_max > GAMMA_MAX_REASONABLE:
        warnings.append(f"DRT contains very large values (max: {gamma_max:.2e} Ohm)")

    if np.sum(gamma) < 1e-10:
        warnings.append("DRT is nearly zero - data may have no relaxation structure")

    return NNLSSolution(
        gamma=gamma,
        success=True,
        n_inductive_points=n_inductive,
        inductive_fraction=inductive_fraction,
        max_inductive_imag=max_inductive,
        gamma_max=gamma_max,
        gamma_min_nonzero=gamma_min_nonzero,
        condition_number=cond_reg,
        warnings=warnings
    )
