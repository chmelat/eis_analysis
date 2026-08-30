"""
Ranking of competing circuit models by information criteria.

Adding an element to a circuit almost always lowers the residual, so a
comparison of fit errors systematically favours the most complex candidate.
AIC and BIC charge for each free parameter and answer the question the fit
error cannot: does the data support the extra element?

Kept separate from fitting.py so the ranking is testable without the CLI.
"""

import logging
from typing import List, NamedTuple, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ...fitting import FitResult
from ...fitting.diagnostics import compute_information_criteria

logger = logging.getLogger(__name__)


class ModelScore(NamedTuple):
    """One candidate circuit and its information criteria."""

    index: int                    # position on the command line, 1-based
    expression: str
    result: Optional[FitResult]   # None when the fit failed or scored non-finite
    aic: float
    bic: float


def score_candidates(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    candidates: List[Tuple[str, Optional[FitResult]]],
    weighting: str
) -> List[ModelScore]:
    """
    Score fitted candidates and sort them best-first by BIC.

    Parameters
    ----------
    frequencies : ndarray of float
        Measurement frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ohm]
    candidates : list of (str, FitResult or None)
        Circuit expression and its fit, in command-line order; None marks a
        candidate that failed
    weighting : str
        Weighting used for every fit - the comparison is only valid when it
        is the same for all candidates

    Returns
    -------
    scores : list of ModelScore
        Sorted by BIC ascending. Failed candidates carry inf and sort last;
        `index` keeps the original command-line position.
    """
    scores = []
    for index, (expression, result) in enumerate(candidates, start=1):
        if result is None:
            scores.append(ModelScore(index, expression, None, np.inf, np.inf))
            continue

        Z_fit = result.circuit.impedance(frequencies, list(result.params_opt))
        _, aic, bic = compute_information_criteria(
            Z, Z_fit, weighting, result.n_free_params
        )

        # A fit whose predicted impedance contains NaN scores NaN, and NaN
        # compares false against everything - it would neither sort last nor
        # be rejected, and would poison the reference values in the table.
        # Such a model cannot be scored at all, so report it as failed.
        if not (np.isfinite(aic) and np.isfinite(bic)):
            logger.warning(f"Candidate {index} ({expression}) produced a "
                           "non-finite score and cannot be ranked")
            scores.append(ModelScore(index, expression, None, np.inf, np.inf))
            continue

        scores.append(ModelScore(index, expression, result, aic, bic))

    return sorted(scores, key=lambda s: s.bic)


def log_comparison(scores: List[ModelScore], n_points: int, weighting: str) -> None:
    """
    Log the ranked comparison table.

    Parameters
    ----------
    scores : list of ModelScore
        Output of score_candidates (already sorted)
    n_points : int
        Number of measurement frequencies; the residual count is twice this
    weighting : str
        Weighting used, reported so the table states its own validity domain
    """
    n_residuals = 2 * n_points
    fitted = [s for s in scores if s.result is not None]
    if not fitted:
        return

    # Only differences are meaningful, so reference everything to the best
    best_aic = min(s.aic for s in fitted)
    best_bic = min(s.bic for s in fitted)

    width = min(max(max(len(s.expression) for s in scores), 24), 44)

    logger.info("=" * 60)
    logger.info(f"Circuit comparison (n = {n_residuals} residuals, "
                f"weighting = {weighting})")
    logger.info("=" * 60)
    # 'rank' orders by BIC, '-c' is the position on the command line - the
    # figures are saved under that number, not under the rank.
    logger.info(f"  {'rank':>4}  {'-c':>2}  {'Circuit':<{width}}  {'k':>2}  "
                f"{'err%':>6}  {'dAIC':>7}  {'dBIC':>7}  {'cond':>8}")

    rank = 0
    for score in scores:
        expression = score.expression
        if len(expression) > width:
            expression = expression[:width - 3] + '...'

        if score.result is None:
            logger.info(f"  {'':>4}  {score.index:>2}  {expression:<{width}}  "
                        f"{'-':>2}  {'-':>6}  {'-':>7}  {'-':>7}  {'-':>8}  failed")
            continue

        rank += 1
        flag = '' if score.result.is_well_conditioned else '  !'
        logger.info(
            f"  {rank:>4}  {score.index:>2}  {expression:<{width}}  "
            f"{score.result.n_free_params:>2}  "
            f"{score.result.fit_error_rel:>6.2f}  "
            f"{score.aic - best_aic:>7.1f}  {score.bic - best_bic:>7.1f}  "
            f"{score.result.condition_number:>8.1e}{flag}"
        )

    logger.info("")
    logger.info("  dAIC/dBIC < 2: indistinguishable, 4-7: noticeable, "
                "> 10: decisive")
    if any(s.result is not None and not s.result.is_well_conditioned
           for s in scores):
        logger.info("  ! cond > 1e10 - covariance unreliable, model likely "
                    "overparametrized")
    logger.info("")
    logger.info(f"Selected by BIC: {scores[0].expression}  "
                f"(candidate {scores[0].index})")
    logger.info("=" * 60)
