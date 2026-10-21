"""
Statistical significance of DRT peaks against the noise of the data.

A local maximum of gamma counts as a separate process only if the data need
it: the DRT is refitted under the null hypothesis that the maximum is merely
a shoulder of its taller neighbour, and the loss of fit is compared with the
noise. A height threshold relative to the tallest peak cannot do this - it
drops a resolved process three orders smaller in R and keeps the lobes into
which regularization splits one broad process.

The null hypothesis is a shoulder, not gamma = 0 around the candidate: a
lobe of a broad process (ZARC, n = 0.8) carries real mass, so removing it
worsens the fit (Delta chi^2 46-1746 in testing) and it would pass. A shoulder
keeps the mass and forbids only the separate maximum: gamma must not
increase from the neighbour's apex through the candidate and the valley
beyond it.

CONTRACT

Callers. calculate_drt (core.py) once per run: the run's matrices and
lambda, gamma = the run's NNLS solution in Ohm (before any R_pol
normalization), the Z, R_inf and L_series the system was built and solved
with, noise_sigma None. probe_lambda_stability (stability.py) per probe:
the same matrices and Z, the probe's lambda and NNLS solution, and
noise_sigma = the main run's noise_sigma. peaks.py reads the result to
report the scipy peaks, core.py lists the rejected maxima with their
delta_chi2. From DRTMatrices it reads A, b, L, tau, A_re, A_im, omega and
L_series_scale; the caller guarantees gamma and L_series minimize
||A x - b||^2 + lambda ||L x||^2 over x >= 0 (up to NNLS convergence).

assess_peak_significance(matrices, lambda_reg, gamma, Z, R_inf, L_series,
noise_sigma, threshold) -> PeakSignificanceResult.
  chi2_full: pseudo chi^2 (Boukamp, weight 1/|Z|^2 per point, real and
  imaginary parts) of R_inf + j omega L_series + A gamma against Z.
  sigma: noise_sigma when given, else sqrt(chi2_full / 2N), N = len(Z) (the
  relative residual per component); floored at machine epsilon.
  Candidates: the strict local maxima of gamma as scipy.signal.find_peaks
  finds them (a flat top counts once, at its middle; the first and last
  grid point are never candidates), ascending.
  Candidates are tested lowest first (ties in index order), so that the
  smaller maxima are decided before the taller ones.
  The taller neighbours of j: the nearest strictly taller candidate on
  each side. With none on either side, j is not tested: delta inf. They
  are tested in the order of their col (the minimum of gamma between j
  and it, both included), the higher first; the second only when the
  first gives delta >= threshold. delta is the smallest of those tested:
  a separate process must not pass as a shoulder of either.
  The boundary beyond j on each side: the nearest candidate that is at
  least as tall as j or not rejected (delta not < threshold), else the
  grid end. Smaller maxima rejected on their own lie inside the valley.
  The monotone run for a neighbour q: q on the left - from q's apex to the
  lowest point of gamma from j up to the right boundary (excluded), gamma
  non-increasing; q on the right - from the lowest point from the left
  boundary to j to q's apex, gamma non-decreasing. Among equal minima (a
  stretch of NNLS zeros) the one farthest from j: the whole valley belongs
  to the run. Smaller candidates in the run are forced into it too.
  Refit: the same regularized NNLS (same A, b, L, lambda; L_series
  re-solved, R_inf fixed) with gamma monotone on the run and >= 0
  everywhere. delta = (chi2_shoulder - chi2_full) / sigma^2, both pseudo
  chi^2 against Z. A refit that raises RuntimeError or ValueError is left
  out; delta is NaN when no refit of j succeeded. delta may be negative:
  the refit minimizes the regularized objective, not chi^2.
  significant: the candidates with delta not < threshold (inf and NaN
  kept: a maximum that cannot be tested is not rejected), ascending.
  noise_sigma and threshold in the result are the sigma and threshold used.

What must hold:
  - The inputs are not changed.
  - candidates, delta_chi2 align element by element; significant is a
    subset of candidates.
  - A candidate with no strictly taller candidate on either side (the
    tallest, and any as tall) is always significant.
  - delta does not depend on the unit of Z: Z, R_inf, L_series and gamma
    scaled by one constant (with matrices built from the scaled Z, as
    _build_drt_matrices does) give the same delta and sigma.
  - A gamma with no interior maximum (constant, monotone, all zero) gives
    empty candidates and significant, never an exception.
  - Only RuntimeError and ValueError of the refit are caught; nothing else
    raises for a gamma that the caller's guarantee covers.
  - threshold and sigma decide which smaller maxima are rejected, and with
    them the boundaries and which neighbours are tested; so delta_chi2
    depends on them beyond the factor 1 / sigma^2, and significant is not
    monotone in threshold (rejecting a small maximum can lengthen the run
    of its taller neighbour and let it pass).
  - The test is symmetric under reversing gamma, except where find_peaks
    places the candidate of a flat top of even width.
  - A finite noise_sigma > 0 and a finite threshold are the caller's
    guarantee (NaN in either keeps every candidate).
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
from numpy.typing import NDArray
from scipy.optimize import nnls
from scipy.signal import find_peaks

from .results import DRTMatrices
from .linear_system import _reconstruct
from .gcv import DRT_NNLS_MAXITER_FACTOR
from ..fitting.config import DRT_PEAK_DCHI2_MIN
from ..validation.noise import compute_pseudo_chisqr


@dataclass
class PeakSignificanceResult:
    """
    Outcome of the shoulder test for every local maximum of gamma.

    Attributes
    ----------
    candidates : ndarray of int
        Indices of all local maxima of gamma, ascending.
    delta_chi2 : ndarray of float
        Per candidate: increase of pseudo chi^2 under the shoulder hypothesis,
        in units of noise_sigma^2. inf for a peak with no taller peak on
        either side (never tested), NaN when the constrained refit failed.
    significant : ndarray of int
        The candidates not below threshold (inf and NaN included: a peak
        that cannot be tested is not rejected), ascending.
    noise_sigma : float
        Relative noise per component used to scale delta_chi2.
    threshold : float
        Delta chi^2 needed for significance.
    """
    candidates: NDArray[np.int_]
    delta_chi2: NDArray[np.float64]
    significant: NDArray[np.int_]
    noise_sigma: float
    threshold: float


def _shoulder_chi2(matrices: DRTMatrices, lambda_reg: float, Z: NDArray,
                   R_inf: float, start: int, stop: int, falling: bool) -> float:
    """
    Pseudo chi^2 of the DRT refitted with gamma monotone on [start, stop].

    falling=True: non-increasing (the taller peak is at start), so
    gamma_k = sum_{i>=k} e_i; otherwise non-decreasing, gamma_k = sum_{i<=k} e_i.
    With e_i >= 0 both are non-negative, so the substitution keeps the
    problem an NNLS: the segment's columns of A and of the regularization
    matrix are multiplied by U, the rest (other bins, L column) stay.
    """
    n = stop - start + 1
    U = np.triu(np.ones((n, n))) if falling else np.tril(np.ones((n, n)))

    def substitute(M: NDArray) -> NDArray:
        return np.hstack([M[:, :start], M[:, start:stop + 1] @ U, M[:, stop + 1:]])

    A_reg = np.vstack([substitute(matrices.A), np.sqrt(lambda_reg) * substitute(matrices.L)])
    b_reg = np.concatenate([matrices.b, np.zeros(matrices.L.shape[0])])
    x, _ = nnls(A_reg, b_reg, maxiter=DRT_NNLS_MAXITER_FACTOR * A_reg.shape[1])  # may raise
    x = np.concatenate([x[:start], U @ x[start:stop + 1], x[stop + 1:]])

    n_tau = len(matrices.tau)
    L_series = (float(x[n_tau]) * matrices.L_series_scale
                if matrices.L_series_scale is not None else 0.0)
    return compute_pseudo_chisqr(Z, _reconstruct(matrices, x[:n_tau], L_series, R_inf))


def assess_peak_significance(matrices: DRTMatrices, lambda_reg: float,
                             gamma: NDArray, Z: NDArray, R_inf: float,
                             L_series: float, noise_sigma: Optional[float] = None,
                             threshold: float = DRT_PEAK_DCHI2_MIN
                             ) -> PeakSignificanceResult:
    """
    Test every local maximum of gamma against the shoulder hypothesis.

    Parameters
    ----------
    matrices : DRTMatrices
        System of the run that produced gamma (weighted A, regularization L).
    lambda_reg : float
        Regularization parameter of that run.
    gamma : ndarray
        Physical (not R_pol-normalized) DRT [Ohm].
    Z : ndarray of complex
        Measured impedance [Ohm].
    R_inf, L_series : float
        High-frequency resistance [Ohm] and series inductance [H] of the run.
    noise_sigma : float, optional
        Relative noise per component to scale Delta chi^2 by. Default: the
        DRT's own residual sqrt(chi^2 / 2N). Not the Lin-KK noise estimate:
        derived from pseudo chi^2, it absorbs genuine misfit (drift reached
        29 % on a measured spectrum) and is computed on the spectrum before
        any frequency trimming. The lambda probe passes the main run's value,
        so that every probe solution is judged against the same noise.
    threshold : float
        Delta chi^2 needed for significance (default DRT_PEAK_DCHI2_MIN).

    Returns
    -------
    PeakSignificanceResult
    """
    chi2_full = compute_pseudo_chisqr(Z, _reconstruct(matrices, gamma, L_series, R_inf))
    sigma = (float(noise_sigma) if noise_sigma is not None
             else float(np.sqrt(chi2_full / (2 * len(Z)))))
    # Exact reconstruction would divide by zero; a DRT never gets there
    # (>= 0.6 % on noise-free synthetics), the floor only keeps it finite.
    sigma = max(sigma, float(np.finfo(float).eps))

    candidates, _ = find_peaks(gamma)
    delta = np.full(len(candidates), np.inf)
    heights = gamma[candidates]
    # Lowest first: the run of a candidate absorbs the smaller maxima beyond
    # it that are not significant themselves, decided by then
    for j in np.argsort(heights, kind='stable'):
        # The peaks it could be a shoulder of: the nearest strictly taller one
        # on each side. Smaller peaks in between are forced into the same
        # monotone run.
        sides = [k for k in (next((k for k in range(j - 1, -1, -1) if heights[k] > heights[j]), None),
                             next((k for k in range(j + 1, len(candidates)) if heights[k] > heights[j]), None))
                 if k is not None]
        if not sides:
            continue
        a = int(candidates[j])
        # Beyond the candidate the run reaches up to the next maximum that is
        # taller or significant (the grid end without one): a smaller maximum
        # rejected on its own is part of the valley, not a boundary
        def kept(k: int) -> bool:
            return bool(heights[k] >= heights[j] or not delta[k] < threshold)
        left = next((int(candidates[k]) for k in range(j - 1, -1, -1) if kept(k)), 0)
        right = next((int(candidates[k]) for k in range(j + 1, len(candidates)) if kept(k)),
                     len(gamma))
        # A separate process must not pass as a shoulder of either taller
        # neighbour: the one behind the higher col first (as for topographic
        # prominence, usually the cheaper merge), the other only when the
        # first does not already reject it; equal cols (both valleys at the
        # NNLS zero, common in a sparse DRT) need not be detected exactly.
        sides.sort(key=lambda k: -float(np.min(
            gamma[min(a, candidates[k]):max(a, candidates[k]) + 1])))
        tested = []
        for q in sides:
            # The run ends at the lowest point before the boundary, the one
            # farthest out among equal minima: a stretch of NNLS zeros belongs
            # to the run, or the refit moves the maximum just past its end.
            # A rise after that point (a pile-up from beyond the window at the
            # grid end) is not part of the shoulder.
            if q < j:  # taller peak on the left: falling through the valley
                beyond = gamma[a:right]
                stop = a + len(beyond) - 1 - int(np.argmin(beyond[::-1]))
                start, falling = int(candidates[q]), True
            else:
                start = left + int(np.argmin(gamma[left:a + 1]))
                stop, falling = int(candidates[q]), False
            try:
                chi2_shoulder = _shoulder_chi2(matrices, lambda_reg, Z, R_inf,
                                               start, stop, falling)
            except (RuntimeError, ValueError):
                # scipy's nnls raises at the iteration limit; leave untested
                continue
            tested.append((chi2_shoulder - chi2_full) / sigma ** 2)
            if tested[-1] < threshold:
                break
        delta[j] = min(tested) if tested else np.nan

    return PeakSignificanceResult(
        candidates=candidates, delta_chi2=delta,
        significant=candidates[~(delta < threshold)],
        noise_sigma=sigma, threshold=threshold)
