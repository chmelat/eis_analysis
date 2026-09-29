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
from .gcv import find_optimal_lambda_gcv, find_optimal_lambda_hybrid, DRT_NNLS_MAXITER_FACTOR
from ..fitting.diagnostics import compute_weights
from ..fitting.config import DRT_LAMBDA_DEFAULT, DRT_LAMBDA_RANGE

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
                        R_inf: float, n_tau: int = 100,
                        weighting: str = 'uniform',
                        tau_extend_decades: float = 0.0,
                        inductance: bool = False) -> DRTMatrices:
    """
    Build DRT system matrices A, b, and regularization matrix L.

    ``weighting`` scales each frequency's rows of A and b by a weight from
    ``compute_weights``, rescaled so that ||w*Z|| = ||Z||: the residual stays
    commensurate with ||L gamma||, so lambda means the same under every
    weighting. A_re and A_im stay unweighted: they map gamma back to impedance.

    ``n_tau`` points span the measured window; ``tau_extend_decades`` appends
    points at the same log spacing beyond its slow end, so a process slower
    than the lowest frequency gets a place on the grid instead of piling up
    in the last bin. The fast end is not extended: R_inf is subtracted
    beforehand and an RC with tau << 1/omega_max is indistinguishable from it.

    ``inductance`` appends an unregularized column for a series j*omega*L
    (zero column in the regularization matrix), so lambda selection and NNLS
    see it as one more non-negative unknown. Without it the model can only
    produce Im(Z) < 0, and an inductive high-frequency end deforms gamma.
    """
    f_max = frequencies.max()
    f_min = frequencies.min()

    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    omega = 2 * np.pi * frequencies

    d_ln_tau_array = np.diff(np.log(tau))
    d_ln_tau = float(np.mean(d_ln_tau_array))

    n_ext = int(np.ceil(tau_extend_decades * np.log(10) / d_ln_tau))
    if n_ext > 0:
        tau = np.concatenate([tau, tau_max * np.exp(d_ln_tau * np.arange(1, n_ext + 1))])
    n_grid = len(tau)

    # Vectorized matrix construction
    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom

    weights = compute_weights(Z, weighting)
    weighted_norm = float(np.linalg.norm(weights * Z))
    if weighted_norm > 0:
        weights = weights * (float(np.linalg.norm(Z)) / weighted_norm)
    A = np.vstack([weights[:, None] * A_re, weights[:, None] * A_im])

    b = np.concatenate([weights * (Z.real - R_inf), weights * Z.imag])

    # Regularization matrix: second derivative in ln(tau), scaled so that
    #   ||A gamma - b||^2 + lambda ||L gamma||^2
    #     = n * [ (1/n) ||A gamma - b||^2 + lambda * integral (gamma'')^2 d ln tau ]
    # with n residuals. The [1, -2, 1] difference / d^2 is gamma'', sqrt(d)
    # turns the sum of squares into the rectangle-rule integral, and sqrt(n)
    # makes the misfit a mean; a bare [1, -2, 1] left lambda proportional to
    # d^3 and to n. See DRT_LAMBDA_RANGE.
    L = np.zeros((n_grid - 2, n_grid))
    np.fill_diagonal(L, 1)           # Main diagonal at offset 0
    np.fill_diagonal(L[:, 1:], -2)   # Diagonal at offset 1
    np.fill_diagonal(L[:, 2:], 1)    # Diagonal at offset 2
    L *= np.sqrt(len(b)) / d_ln_tau ** 1.5

    L_series_scale = None
    if inductance:
        L_series_scale = 1.0 / float(omega.max())
        column = np.concatenate([np.zeros_like(omega), weights * omega * L_series_scale])
        A = np.hstack([A, column[:, None]])
        L = np.hstack([L, np.zeros((L.shape[0], 1))])
    cond_A = float(np.linalg.cond(A))

    return DRTMatrices(
        A=A, A_re=A_re, A_im=A_im, b=b, L=L,
        tau=tau, d_ln_tau=d_ln_tau, condition_number=cond_A,
        weights=weights, tau_window=(float(tau_min), float(tau_max)),
        omega=omega, L_series_scale=L_series_scale
    )


def _reconstruct(matrices: DRTMatrices, gamma: NDArray, L_series: float,
                 R_inf: float) -> NDArray:
    """Model impedance R_inf + j*omega*L + integral of gamma, unweighted [Ohm]."""
    return (R_inf + (matrices.A_re + 1j * matrices.A_im) @ gamma
            + 1j * matrices.omega * L_series)


def _select_lambda(A: NDArray, b: NDArray, L: NDArray,
                   lambda_reg: Optional[float] = None,
                   auto_lambda: bool = False) -> LambdaSelection:
    """
    Select regularization parameter lambda.
    """
    # Edge detection (F3/F7): lambda landing at a bound - or the GCV guess
    # pinning there even when L-curve corrected it - signals the optimizer
    # wants more extreme regularization than DRT_LAMBDA_RANGE allows.
    def _at_bound(lam: Optional[float]) -> bool:
        return lam is not None and not DRT_LAMBDA_RANGE[0] < lam < DRT_LAMBDA_RANGE[1]

    if auto_lambda:
        try:
            lambda_opt, gcv_score, diag = find_optimal_lambda_hybrid(
                A, b, L,
                n_search=20,
                lcurve_decades=1.5
            )
            lambda_gcv = diag.get('lambda_gcv')
            corner_at_edge = diag.get('corner_at_edge', False)
            at_edge = _at_bound(lambda_opt) or _at_bound(lambda_gcv) or corner_at_edge
            # Always 'hybrid' when this path succeeds - which of the two stages
            # won is hybrid_stage, not a different method. 'gcv' below now means
            # only the fallback after the hybrid search failed.
            return LambdaSelection(
                lambda_value=lambda_opt,
                method='hybrid',
                lambda_gcv=lambda_gcv,
                lambda_lcurve=diag.get('lambda_lcurve'),
                hybrid_stage=diag.get('method_used'),
                gcv_score=gcv_score,
                corner_at_edge=corner_at_edge,
                lambda_at_edge=at_edge
            )
        except (np.linalg.LinAlgError, ValueError):
            logger.debug("Hybrid lambda selection failed, falling back to GCV",
                         exc_info=True)
            try:
                lambda_opt, gcv_score = find_optimal_lambda_gcv(
                    A, b, L, n_search=20
                )
                return LambdaSelection(
                    lambda_value=lambda_opt,
                    method='gcv',
                    gcv_score=gcv_score,
                    lambda_at_edge=_at_bound(lambda_opt)
                )
            except (np.linalg.LinAlgError, ValueError):
                logger.debug("GCV lambda selection failed, using fallback "
                             f"lambda={DRT_LAMBDA_DEFAULT}", exc_info=True)
                return LambdaSelection(lambda_value=DRT_LAMBDA_DEFAULT, method='fallback')

    if lambda_reg is None:
        return LambdaSelection(lambda_value=DRT_LAMBDA_DEFAULT, method='default')

    return LambdaSelection(lambda_value=lambda_reg, method='user')


def _solve_nnls(matrices: DRTMatrices, lambda_reg: float,
                Z: NDArray) -> NNLSSolution:
    """
    Solve regularized NNLS problem; split off L_series if the model has it.
    """
    warnings = []

    # Check for inductive data
    n_inductive = int(np.sum(Z.imag > 0))
    inductive_fraction = n_inductive / len(Z) if len(Z) > 0 else 0.0
    max_inductive = float(np.max(Z.imag)) if n_inductive > 0 else 0.0

    # With an L column the model represents Im(Z) > 0, so this is no longer
    # a data/model mismatch worth a warning.
    if matrices.L_series_scale is None and (inductive_fraction > 0.1 or max_inductive > 50):
        warnings.append(
            f"{n_inductive} points with inductive component ({inductive_fraction*100:.1f}%), "
            f"max Z'' = {max_inductive:.2f} Ohm"
        )

    # Build regularized system
    A_reg = np.vstack([matrices.A, np.sqrt(lambda_reg) * matrices.L])
    b_reg = np.concatenate([matrices.b, np.zeros(matrices.L.shape[0])])

    # Condition number of the system actually solved. The bare kernel A is
    # intrinsically ill-conditioned for any DRT problem (that is why Tikhonov
    # regularization is applied); the regularized system is what determines the
    # numerical health of the solve.
    cond_reg = float(np.linalg.cond(A_reg))

    # Solve NNLS
    try:
        x, _ = nnls(A_reg, b_reg, maxiter=DRT_NNLS_MAXITER_FACTOR * A_reg.shape[1])
    except Exception as e:
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=[f"NNLS solver error: {e}"]
        )

    # Validate solution
    if np.any(~np.isfinite(x)):
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=["NNLS returned NaN or Inf values"]
        )

    n_grid = len(matrices.tau)
    gamma = x[:n_grid]
    L_series = (float(x[n_grid]) * matrices.L_series_scale
                if matrices.L_series_scale is not None else 0.0)

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
        L_series=L_series,
        n_inductive_points=n_inductive,
        inductive_fraction=inductive_fraction,
        max_inductive_imag=max_inductive,
        gamma_max=gamma_max,
        gamma_min_nonzero=gamma_min_nonzero,
        condition_number=cond_reg,
        warnings=warnings
    )
