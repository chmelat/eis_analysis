"""
Mu metric optimization for Voigt chain fitting (Lin-KK style).

This module provides the mu metric calculation and automatic optimization
of the number of Voigt elements.

CONTRACT

Metric. For the resistances R_1..R_M of a Voigt chain fitted with
negative values allowed, N = sum |R_k| over R_k < 0 (the negative mass)
and P = sum R_k over R_k >= 0. mu = 1 - N / P (Schoenleber 2014): 1 with
no negative R_k, falling as the negative mass grows, below 0 when it
exceeds the positive one. More negative mass never raises mu: mu is
non-increasing in N at fixed P and non-decreasing in P at fixed N, so the
all-negative case (P = 0, N > 0) is the lowest value, -inf. With no
resistance at all (N = P = 0) there is no overfit: 1. Lin-KK stops adding
elements once negative R_k appear in quantity, the sign of fitting noise.

Search. find_optimal_M_mu fits M = min_M, min_M + 1, ... on
generate_tau_grid_fixed_M(frequencies, M, extend_decades) with
estimate_R_linear(allow_negative as passed) and stops at the first M
whose mu <= mu_threshold, or at max_M. At least one M is fitted (min_M),
whatever the threshold.

Callers.
- validation.kramers_kronig.lin_kk_native: min_M from
  kk_grid._chi2_lower_M (3 <= min_M <= max_M), max_M <= N - 2 (- [L] -
  [C] for 'imag'), extend_decades 0, include_Rs and allow_negative True,
  mu_threshold from the user (CLI --mu-threshold, default 0.85; only NaN
  is rejected here), fit_type and weighting as passed ('real' and
  'modulus' through kramers_kronig_validation). It reads M, mu, tau,
  elements, L_value, C_value and warnings; mu becomes min_mu of
  kk_grid.find_optimal_extend_decades, KKResult.mu and the CLI's stop
  line. f and Z are validated before (validate_eis_data).
- fitting.voigt_chain.fitting.fit_voigt_chain_linear (CLI --voigt-chain
  --voigt-auto-M): min_M 3, max_M from --voigt-max-M (default 50, not
  capped by the number of points), mu_threshold from
  --voigt-mu-threshold, extend_decades from --voigt-extend-decades,
  include_C False, allow_negative True, fit_type default 'complex'. It
  keeps only tau and refits it with NNLS; the CLI prints iterations,
  reached_max_M, warnings, n_negative / n_R and the stop M and mu.
- kk_grid.find_optimal_extend_decades calls calc_mu(R_1..R_M) per
  extension candidate and compares it with min_mu.
Public API (eis_analysis.fitting.voigt_chain: calc_mu, find_optimal_M_mu,
MuOptimization): any caller, nothing guaranteed beyond the rules below.

calc_mu(R_i) -> float.
  The metric above over a 1-D float64 array (R_s and L not included;
  the fit returns float64: float32 sums overflow sooner, an int64 abs
  wraps at its minimum, and nothing is promised for them). A
  NaN or infinite R_k gives NaN: an unfitted or overflowed element is no
  resistance to count, and a mu that left it out or let it swamp the
  sums would pass for a judgement. An empty array has no resistance: 1.
  Invariant under permuting R_i and scaling it by c > 0 (while the R_k
  stay normal floats: scaled to ~1e-320 they underflow). No
  RuntimeWarning while the sums N and P stay finite (|R_k| summing below
  ~1e308; the fit of |Z| within 1e-150..1e150 keeps them far below).

find_optimal_M_mu(frequencies, Z, mu_threshold=0.85, max_M=50,
extend_decades=0, include_Rs, include_L, include_C, fit_type, allow_negative,
weighting, min_M=3) -> MuOptimization.
  M: the stop M, min_M <= M <= max_M. mu: calc_mu of its R_1..R_M (after
  R_s where include_Rs, before L where include_L; C is never in
  elements). tau, elements, L_value, C_value: the fit at M, as
  estimate_R_linear returned it, so that the caller can rebuild Z_fit
  from them. reached_max_M: True exactly when no fitted M had mu <=
  mu_threshold, i.e. the search ended at max_M for want of M (with it,
  a warning saying M is the cap, not a mu choice); False when the stop
  M meets the threshold. n_negative, n_R: R_k < 0 and all R_k at M.
  iterations: one MuIteration per sampled step - min_M, every M
  divisible by 5, and the stop M always (the CLI prints them as the
  progress of the search); iteration counts from 1 at min_M, residual is
  estimate_R_linear's. warnings: also two lines when allow_negative is
  False (mu is then 1 by construction, so the search runs to max_M).
  mu_threshold >= 1 stops at min_M where mu there is a number (mu <= 1
  always); inf as well;
  -inf stops only at an M whose R_k are all negative (mu -inf), else at
  max_M. ValueError for NaN mu_threshold, for
  min_M < 1 or min_M > max_M, and where estimate_R_linear raises
  (validate_eis_data: NaN, inf, f <= 0, Z = 0, unequal lengths; unknown
  fit_type or weighting); LinAlgError where its SVD does not converge.
  A mu that is not a number (NaN) is not a stop: it does not meet the
  threshold.
  Beyond a determined fit (more unknowns than equations: 1 + M + [L] +
  [C] above N for 'real', above 2N for 'complex') the R_k are not
  unique and mu means nothing; the caller keeps max_M below that
  (lin_kk_native does, fit_voigt_chain_linear does not).
  The caller passes integer min_M and max_M, a float mu_threshold,
  extend_decades finite and >= 0, and f with at least two distinct
  values; nothing is promised otherwise: max_M 7.5 ends at M = 8 (in
  the iterations, with the warning naming 7.5), a float min_M and a
  mu_threshold of None are a TypeError,
  max_M = inf with a threshold never met does not end; a negative or NaN
  extend_decades fits the grid of 0, inf is the fit's ValueError after
  a RuntimeWarning; equal f collapse the tau grid to one value (mu 1,
  the search runs to max_M).

What must hold:
  - The inputs are not changed. The result is deterministic.
  - On what the caller guarantees, raises only ValueError and
    LinAlgError, only where stated; no RuntimeWarning on validated data.
  - Units: Z scaled by c > 0 leaves M, mu, reached_max_M, n_negative and
    the iterations' mu unchanged (to rounding) and scales the R in
    elements and L_value by c, C_value by 1/c, while |Z| stays within
    1e-150..1e150 (validate_eis_data); frequencies scaled by c leave M
    and mu unchanged and scale tau, L_value and C_value by 1/c for
    fit_type 'real' and 'complex' (for 'imag' the conditioning of the fit
    decides, see validation/kk_grid.py).
  - Permuting the points (f and Z together) changes no result (to
    rounding).
  - mu_threshold: a lower threshold never stops at a smaller M (the
    search is the first M at or below it).
"""

import numpy as np
import logging
from dataclasses import dataclass, field
from typing import List, Optional
from numpy.typing import NDArray

from ...utils.validation import validate_eis_data
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
        Mu value at that M (at or below mu_threshold on normal
        termination, above it or NaN only when max_M was reached)
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
        - -inf: no positive R_i, only negative ones; 1.0 with no
          resistance at all (empty or all zero); NaN where an R_i is NaN
          or infinite

    Notes
    -----
    In Lin-KK test, used as stop condition for iterative determination
    of optimal number of RC elements. Typical threshold: stop at
    mu <= 0.85.

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
    # An unfitted or overflowed element is no resistance to count: NaN fails
    # both tests below and would drop out of both sums, inf would swamp
    # them ([inf, 1, -0.5] gave 1.0, [-inf, 1, 2] -inf)
    if not np.isfinite(R_i).all():
        return float('nan')

    neg_sum = float(np.sum(np.abs(R_i[R_i < 0])))
    pos_sum = float(np.sum(np.abs(R_i[R_i >= 0])))

    if pos_sum == 0:
        # The limit of 1 - N/P: -1 here ranked an all-negative fit above
        # one with a little positive mass (-99 for [1, -100]). -inf still
        # stops the M search or rejects an extend_decades candidate. No
        # resistance at all is no overfit.
        logger.debug("calc_mu: no positive R_i")
        return 1.0 if neg_sum == 0 else -np.inf

    mu = 1.0 - neg_sum / pos_sum
    return mu


def find_optimal_M_mu(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    extend_decades: float = 0.0,
    include_Rs: bool = True,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'complex',
    allow_negative: bool = True,
    weighting: str = 'modulus',
    min_M: int = 3
) -> MuOptimization:
    """
    Find optimal number of Voigt elements using mu metric (Lin-KK style).

    Iteratively increases M from min_M until mu <= mu_threshold or M >= max_M.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85, as in Lin-KK)
        Lower values -> the iteration stops later -> more elements
        (higher overfit tolerance); higher values -> fewer elements.
        mu never exceeds 1, so any value >= 1 stops at min_M. NaN raises.
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
    min_M : int, optional
        First M to try (default: 3, Lin-KK standard). Lin-KK validation
        raises it to where pseudo chi^2 levels off, so that a mu dip of an
        under-resolved tau grid does not stop the search.

    Returns
    -------
    MuOptimization
        The stopping M and its fit, plus the progress of the search and any
        caveat about it. Nothing is logged: the caller decides what, if
        anything, to report.

    Notes
    -----
    Algorithm from Schonleber et al. (2014):
    1. Start with M = min_M (3 in the original)
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
    if not 1 <= min_M <= max_M:
        raise ValueError(f"find_optimal_M_mu: need 1 <= min_M <= max_M, "
                         f"got min_M={min_M}, max_M={max_M}")
    if np.isnan(mu_threshold):
        raise ValueError("find_optimal_M_mu: mu_threshold is NaN")

    warnings: List[str] = []
    if not allow_negative:
        warnings.append("mu metric is designed for allow_negative=True (Lin-KK)")
        warnings.append("With NNLS, all R_i >= 0, so mu ~ 1 always")

    iterations: List[MuIteration] = []

    M = min_M - 1  # the loop increments before fitting
    # Not 1.0: mu never exceeds 1, so a threshold >= 1 skipped the loop and
    # returned no fit at all. Every threshold >= 1 now stops at min_M; the
    # M < min_M clause covers mu_threshold = inf, which mu = inf does not.
    mu = np.inf
    iteration = 0
    L_value = None
    C_value = None
    tau = None
    elements = None
    R_i = np.array([])

    # `not mu <= threshold`: a NaN mu is no stop (mu > threshold is False)
    while M < min_M or (not mu <= mu_threshold and M < max_M):
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
        # (the last is either the stop or max_M)
        if M % 5 == 0 or mu <= mu_threshold or M == min_M or M >= max_M:
            iterations.append(MuIteration(
                iteration=iteration, M=M, mu=mu, residual=residual,
                n_negative=int(np.sum(R_i < 0)), n_R=len(R_i)))

    assert tau is not None and elements is not None  # the loop ran at least once

    reached_max_M = not mu <= mu_threshold
    if reached_max_M:
        # A NaN mu ends here too, and NaN > threshold is False
        mu_text = f"{mu:.4f} > {mu_threshold}" if mu > mu_threshold else "nan (not a number)"
        warnings.append(f"Reached max_M = {max_M}, mu = {mu_text}; "
                        f"mu criterion never triggered, M is the cap, not a mu choice")

    return MuOptimization(
        M=M, mu=mu, tau=tau, elements=elements,
        L_value=L_value, C_value=C_value,
        reached_max_M=reached_max_M,
        n_negative=int(np.sum(R_i < 0)), n_R=len(R_i),
        iterations=iterations, warnings=warnings)


__all__ = ['calc_mu', 'find_optimal_M_mu', 'MuOptimization', 'MuIteration']
