"""
Statistical significance of DRT peaks against the noise of the data.

A local maximum of gamma counts as a separate process only if the data need
it: the DRT is refitted under the null hypothesis that the maximum is merely
a shoulder of its taller neighbour, and the loss of fit is compared with the
noise. A height threshold relative to the tallest peak cannot do this - it
drops a resolved process three orders smaller in R and keeps the lobes into
which regularization splits one broad process.

The null hypothesis is a shoulder, not gamma = 0 over the candidate's basin:
a lobe of a broad process (ZARC, n = 0.8) carries real mass, so removing it
worsens the fit (Delta chi^2 46-1746 in testing) and it would pass. A shoulder
keeps the mass and forbids only the separate maximum: gamma must not
increase from the neighbour's apex through the candidate's basin.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional
from numpy.typing import NDArray
from scipy.optimize import nnls
from scipy.signal import find_peaks

from .results import DRTMatrices
from .estimation import _peak_basins
from .linear_system import _reconstruct
from .gcv import DRT_NNLS_MAXITER_FACTOR
from ..fitting.config import DRT_PEAK_DCHI2_MIN
from ..validation.kramers_kronig import compute_pseudo_chisqr


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
    basins = _peak_basins(gamma, candidates)
    for j in range(len(candidates)):
        # The peak it would be a shoulder of: the nearest strictly taller one
        # on each side, of which the one behind the higher valley (the col a
        # merge has to cross), as for topographic prominence. Smaller peaks
        # in between are forced into the same monotone run.
        sides = [k for k in (next((k for k in range(j - 1, -1, -1) if heights[k] > heights[j]), None),
                             next((k for k in range(j + 1, len(candidates)) if heights[k] > heights[j]), None))
                 if k is not None]
        if not sides:
            continue
        a = int(candidates[j])
        q = max(sides, key=lambda k: float(np.min(
            gamma[min(a, candidates[k]):max(a, candidates[k]) + 1])))
        # The run ends at the basin's lowest point beyond the candidate (the
        # one farthest out among equal minima, e.g. a stretch of NNLS zeros),
        # not at the basin edge: for the outermost candidate the basin runs
        # to the end of the grid, and a rise there (a pile-up from beyond the
        # window) is not part of the shoulder.
        if q < j:  # taller peak on the left: falling through the basin
            beyond = gamma[a:basins[j][1]]
            stop = a + len(beyond) - 1 - int(np.argmin(beyond[::-1]))
            start, falling = int(candidates[q]), True
        else:
            start = basins[j][0] + int(np.argmin(gamma[basins[j][0]:a + 1]))
            stop, falling = int(candidates[q]), False
        try:
            chi2_shoulder = _shoulder_chi2(matrices, lambda_reg, Z, R_inf,
                                           start, stop, falling)
        except (RuntimeError, ValueError):
            # scipy's nnls raises at the iteration limit; leave untested
            delta[j] = np.nan
            continue
        delta[j] = (chi2_shoulder - chi2_full) / sigma ** 2

    return PeakSignificanceResult(
        candidates=candidates, delta_chi2=delta,
        significant=candidates[~(delta < threshold)],
        noise_sigma=sigma, threshold=threshold)
