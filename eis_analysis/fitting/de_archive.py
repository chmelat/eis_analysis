"""
Second look at a differential evolution run through its evaluation archive.

DE keeps only its best point. When the population collapses into a wrong
basin, the right one was usually visited in the early generations, while the
population was still spread out - and then forgotten. The cost function
records every evaluation; afterwards the best distinct point of each early
window of generations is refined with least_squares, and the best refined
point wins. A refined candidate with a different spectrum but a statistically
equal fit marks an ambiguous model.

Calibration (5 circuits with 6-10 parameters, noise-free and 1 %, 20 seeds:
200 fits, 40 of which DE got wrong; doc/DIFFERENTIAL_EVOLUTION.md):

    candidates from              repaired   flagged   false flags
    lowest cost overall            1/40       9/40       0/160
    4 windows (10/20/40/70 %)      7/40      21/40       0/160
    15 windows, first half        24/40      33/40       0/160
    30 windows, first half        39/40      40/40       0/160
    50 windows, first half        39/40      40/40       0/160

Choosing the lowest-cost points, measured by parameter distance, by spectrum
or with interchangeable blocks sorted, all gave 1/40: they all come from the
late, collapsed population. When a point was evaluated matters, not how
distinctness is measured.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from .config import FIT_QUALITY_ACCEPTABLE_ERROR

# Windows over the early part of the run; 30 is where the repair rate above
# saturates, more only doubles the refinement time (~1 s at 30 windows against
# 10-45 s of DE for 6-10 parameters).
ARCHIVE_WINDOWS = 30

# Fraction of the run's own generations the windows cover. The late half is
# the collapsed population: its points differ only in poorly determined
# parameters and refine back into the basin DE already found.
ARCHIVE_SPAN = 0.5

# Two candidates are distinct when some parameter differs by at least a
# decade (log10 for scale parameters) - or, for a linearly searched one (CPE
# exponent n, n_DQ, alpha_CC, p_YG), by 0.15, about a fifth of the CPE
# exponent's 0.3-1.0 range.
DISTINCT_DECADES = 1.0
EXPONENT_PER_DECADE = 0.15

# A refined candidate must beat the chosen fit by this relative margin to
# replace it, so a candidate refining into the same minimum (equal to ~1e-12)
# does not churn the result.
ADOPT_MARGIN = 1e-6

# An alternative is ambiguous when (a) its cost exceeds the best by less than
# this many residual variances, delta-chi2 = (S_alt - S) / s^2 < 10 with
# s^2 = S / (N - p) - about a 1 % chance for a 1-parameter difference; the
# benchmark did not depend on it between 4 and 25 - and (b) its spectrum
# differs from the best one by MORE than that, sum w^2 |Z_alt - Z|^2 > 10 s^2:
# another model the data could tell apart, not the same minimum reached at a
# slightly different point of a flat valley. (b) scales with the residuals: a
# fixed 1 % in |Z| flagged a two-R|Q fit with 12 % error whose "alternative"
# lay inside its own confidence intervals (EISPOT-M136113-4). Permuted series
# blocks give an identical spectrum and never pass (b).
AMBIGUITY_DCHI2 = 10.0

# The repair message says "local minimum" only when the adopted fit is another
# model, its spectrum different by more than 1 % of |Z| somewhere.
AMBIGUITY_SPECTRUM = 0.01

# Relative precision below which costs are rounding, not data: 1e-6 of |Z|
# per point. The best instrument sigma we have seen is ~1e-5 of |Z| (Zahner
# IM7, doc/ZAHNER_LINK_COMPARISON.md), float rounding is ~1e-15. Without it a
# noise-free fit (cost ~1e-29) made candidates that refined into the same
# minimum look "better" or "distinguishable" at the level of rounding
# (26 false ambiguity flags in 160 correct fits).
PRECISION_FLOOR = 1e-6


def cost_floor(Z: NDArray[np.complexfloating], weights: NDArray[np.float64]) -> float:
    """Weighted SSR of a fit that is off by PRECISION_FLOOR of |Z| at every point."""
    return float(PRECISION_FLOOR ** 2 * np.sum(weights ** 2 * np.abs(Z) ** 2))


def _coords(params: NDArray[np.float64], linear_mask: NDArray[np.bool_]) -> NDArray[np.float64]:
    """Parameters in units where 1 = one decade (or 0.15 of a linear parameter)."""
    p = np.asarray(params, dtype=float)
    log_p = np.log10(np.maximum(np.abs(p), np.finfo(float).tiny))
    return np.where(linear_mask, p / EXPONENT_PER_DECADE, log_p) / DISTINCT_DECADES


def select_archive_candidates(
    costs: NDArray[np.float64],
    params: NDArray[np.float64],
    n_pop: int,
    linear_mask: NDArray[np.bool_],
) -> List[int]:
    """
    Indices of the archive points to refine.

    The first is the lowest-cost point overall (what DE returned); then, for
    each of ARCHIVE_WINDOWS equal windows over the first ARCHIVE_SPAN of the
    generations, the lowest-cost point at least one decade from every point
    picked so far. A window narrower than one generation is widened to one;
    overlaps are harmless, the distance test drops repeats. Points with a
    non-finite cost (an impedance that overflowed) are never picked.

    Parameters
    ----------
    costs : ndarray, shape (n_eval,)
        Cost of every evaluation, in evaluation order
    params : ndarray, shape (n_eval, n_free)
        The evaluated (free, physical) parameter vectors
    n_pop : int
        Evaluations per generation (DE population size)
    linear_mask : ndarray of bool, shape (n_free,)
        True for linearly searched parameters (exponents), False for scale
        parameters (DE's log-search mask, inverted)

    Returns
    -------
    list of int
        At most ARCHIVE_WINDOWS + 1 indices into the archive
    """
    costs = np.asarray(costs, dtype=float)
    finite = np.isfinite(costs)
    if not finite.any():
        return []
    coords = _coords(params, linear_mask)
    picked = [int(np.flatnonzero(finite)[np.argmin(costs[finite])])]
    generation = np.arange(len(costs)) // max(n_pop, 1)
    n_gen = generation[-1] + 1
    edges = np.linspace(0.0, ARCHIVE_SPAN * n_gen, ARCHIVE_WINDOWS + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        window = np.flatnonzero(finite & (generation >= lo) & (generation < max(hi, lo + 1)))
        if len(window) == 0:
            continue
        window = window[np.argsort(costs[window], kind="stable")]
        # Distance of every window point to its nearest picked point, at once
        gap = np.abs(coords[window][:, None, :] - coords[picked][None, :, :]).max(axis=2).min(axis=1)
        far = np.flatnonzero(gap >= 1.0)
        if len(far):
            picked.append(int(window[far[0]]))
    return picked


@dataclass
class Refinement:
    """A least_squares run from one start: its weighted SSR, spectrum and result."""
    cost: float
    Z: NDArray[np.complexfloating]
    params: List[float]          # full parameter vector (fixed ones included)
    result: Any                  # scipy OptimizeResult


@dataclass
class Selection:
    """Outcome of choosing among DE's point and every refinement."""
    chosen: Optional[Refinement]          # None: DE's own point is kept
    best_refined: Optional[Refinement]    # lowest-cost refinement (chosen or not)
    alternatives: List[Tuple[float, Sequence[float]]]   # (delta-chi2, params), see ambiguous_alternatives
    local_minimum: bool                   # an archive candidate replaced another model


def choose(
    de_cost: float,
    from_de: Optional[Refinement],
    from_archive: Sequence[Refinement],
    weights: NDArray[np.float64],
    n_free: int,
) -> Selection:
    """
    One decision among DE's point, the refinement started there and the
    refined archive candidates.

    An archive candidate replaces the refinement from DE only when it is
    better by more than ADOPT_MARGIN and more than rounding (cost_floor), so
    one that refined into the same minimum does not churn the result. The
    best refinement is then kept unless it is worse than DE's own point.
    Ambiguous alternatives are the other refinements, compared with the chosen
    fit (see ambiguous_alternatives).
    """
    best = from_de
    local_minimum = False
    if from_archive:
        cand = min(from_archive, key=lambda r: r.cost)
        if best is None:
            best = cand
        elif best.cost - cand.cost > max(ADOPT_MARGIN * best.cost, cost_floor(best.Z, weights)):
            local_minimum = bool(np.max(np.abs(cand.Z - best.Z) / np.abs(best.Z)) > AMBIGUITY_SPECTRUM)
            best = cand
    chosen = best if best is not None and best.cost <= de_cost else None
    alternatives = []
    if chosen is not None:
        others = [r for r in ([from_de] if from_de else []) + list(from_archive) if r is not chosen]
        alternatives = ambiguous_alternatives(chosen, others, weights, n_free)
    return Selection(chosen, best, alternatives, local_minimum and chosen is not None)


def ambiguous_alternatives(
    best: Refinement,
    refined: Sequence[Refinement],
    weights: NDArray[np.float64],
    n_free: int,
) -> List[Tuple[float, Sequence[float]]]:
    """
    Refined candidates that fit as well as the best one but are another model.

    Parameters
    ----------
    best : Refinement
        The chosen fit
    refined : sequence of Refinement
        The other refined candidates
    weights : ndarray of float
        Per-point weights of the fit (cost = sum w^2 |Z - Z_data|^2)
    n_free : int
        Number of free parameters

    Returns
    -------
    list of (delta_chi2, params)
        delta_chi2 in units of the residual variance, smallest first; several
        candidates that refined into the same alternative count once
    """
    s2 = max(best.cost, cost_floor(best.Z, weights)) / max(2 * len(best.Z) - n_free, 1)

    def separated(a, b):
        return np.sum(weights**2 * np.abs(a.Z - b.Z) ** 2) / s2 > AMBIGUITY_DCHI2

    kept: List[Refinement] = []
    for r in sorted(refined, key=lambda r: r.cost):
        if (r.cost - best.cost) / s2 < AMBIGUITY_DCHI2 and separated(r, best) and all(separated(r, k) for k in kept):
            kept.append(r)
    return [(float((r.cost - best.cost) / s2), r.params) for r in kept]


def selection_warnings(sel: Selection, fit_error_rel: float, labels: Sequence[str]) -> List[str]:
    """
    What the archive check has to say about the chosen fit.

    A repair from the archive, when it is another model; an ambiguous model,
    when the fit is of acceptable quality - for a model that does not fit,
    the residuals are not noise and the delta-chi2 test turns lenient (a
    two-R|Q fit of EISPOT-M136113-4 at 12 % error was flagged with an
    alternative on the edge of its own confidence intervals).
    """
    out = []
    if sel.local_minimum:
        out.append(
            "DE ended in a local minimum; a candidate from its early generations "
            f"refined to a weighted SSR of {sel.chosen.cost:.3g} and is used instead")  # type: ignore[union-attr]
    if sel.alternatives and fit_error_rel < FIT_QUALITY_ACCEPTABLE_ERROR:
        dchi2, alt = sel.alternatives[0]
        out.append(
            f"Ambiguous model: {len(sel.alternatives)} other parameter set(s) with a different "
            f"spectrum fit as well (best delta-chi2 = {dchi2:.1f}): "
            + ", ".join(f"{lab}={v:.4g}" for lab, v in zip(labels, alt)))
    return out
