"""
Tau-grid extension past the slow end of the measured window.

Solves the DRT system on a fixed extension, or lets 'auto' pick the smallest
extension that clears a slow-end pile-up (see DRT_TAU_EXTEND_STEPS).
"""

from typing import Optional, Tuple, Union

from numpy.typing import NDArray

from .estimation import _edge_pile_up_fractions, _lf_rc_ratio, _rpol_from_gamma
from .linear_system import _build_drt_matrices, _select_lambda, _solve_nnls
from .results import DRTMatrices, LambdaSelection, NNLSSolution
from ..fitting.config import (DRT_EDGE_BIN_RPOL_FRACTION, DRT_LF_RC_RATIO_MIN,
                             DRT_TAU_EXTEND_STEPS)


def _solve_on_grid(frequencies: NDArray, Z: NDArray, R_inf: float, n_tau: int,
                   weighting: str, tau_extend: float, inductance: bool,
                   lambda_reg: Optional[float], auto_lambda: bool
                   ) -> Tuple[DRTMatrices, LambdaSelection, NNLSSolution]:
    """Build the system on a grid extended by ``tau_extend`` decades, pick lambda, solve.

    Lambda is selected on the system actually solved, L column included: with
    lambda taken from a run without it, the unregularized unknowns absorb what
    the penalty pushes out of gamma (doc/DRT_RINF_L_ANALYSIS_2026-09-25.md).
    """
    matrices = _build_drt_matrices(frequencies, Z, R_inf, n_tau, weighting, tau_extend,
                                   inductance)
    lambda_sel = _select_lambda(matrices.A, matrices.b, matrices.L, lambda_reg, auto_lambda)
    nnls_result = _solve_nnls(matrices, lambda_sel.lambda_value, Z)
    return matrices, lambda_sel, nnls_result


def _slow_end_piled_up(nnls_result: NNLSSolution, d_ln_tau: float) -> bool:
    """True when the solution heaps more than the pile-up threshold at the slow grid end.

    Measured on the slow end alone: _edge_pile_up reports only the heavier
    end, so a worse fast-end heap (inductance, HF tail) would hide this one.
    """
    if not nnls_result.success or nnls_result.gamma is None:
        return False
    gamma = nnls_result.gamma
    _, slow = _edge_pile_up_fractions(gamma, d_ln_tau, _rpol_from_gamma(gamma, d_ln_tau))
    return slow > DRT_EDGE_BIN_RPOL_FRACTION


def _solve_with_extension(frequencies: NDArray, Z: NDArray, R_inf: float, n_tau: int,
                          weighting: str, tau_extend_decades: Union[float, str],
                          inductance: bool, lambda_reg: Optional[float], auto_lambda: bool
                          ) -> Tuple[DRTMatrices, LambdaSelection, NNLSSolution,
                                     float, Optional[str]]:
    """
    Solve on a fixed extension, or let 'auto' choose one.

    'auto' extends only when the unextended solution piles up at the slow end
    and the low-frequency end is not capacitive (DRT_LF_RC_RATIO_MIN), then
    takes the smallest of DRT_TAU_EXTEND_STEPS that clears the pile-up. If
    none does, the unextended solution stands: a wider grid that still piles
    up only moves the heap further out.

    Returns (matrices, lambda_sel, nnls_result, extension applied, note).
    """
    if tau_extend_decades != 'auto':
        return (*_solve_on_grid(frequencies, Z, R_inf, n_tau, weighting,
                                float(tau_extend_decades), inductance, lambda_reg, auto_lambda),
                float(tau_extend_decades), None)

    base = _solve_on_grid(frequencies, Z, R_inf, n_tau, weighting, 0.0, inductance,
                          lambda_reg, auto_lambda)
    if not _slow_end_piled_up(base[2], base[0].d_ln_tau):
        return (*base, 0.0, "not needed, no slow-end pile-up")

    ratio = _lf_rc_ratio(frequencies, Z)
    if ratio < DRT_LF_RC_RATIO_MIN:
        return (*base, 0.0, f"not extended, capacitive low-frequency end "
                            f"(r = {ratio:.2f} < {DRT_LF_RC_RATIO_MIN})")

    for step in DRT_TAU_EXTEND_STEPS:
        trial = _solve_on_grid(frequencies, Z, R_inf, n_tau, weighting, step, inductance,
                               lambda_reg, auto_lambda)
        if trial[2].success and not _slow_end_piled_up(trial[2], trial[0].d_ln_tau):
            return (*trial, step, "resolved the slow-end pile-up")
    return (*base, 0.0, f"not extended, no extension up to "
                        f"{DRT_TAU_EXTEND_STEPS[-1]} decades resolves the "
                        f"slow-end pile-up")
