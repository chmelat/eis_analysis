"""
The Lin-KK model on its tau grid: impedance of a fitted Voigt chain, the
choice of the grid extension, and where the mu search starts.

Used by lin_kk_native (kramers_kronig.py); reconstruct_impedance also by the
KK plot.
"""

from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .noise import compute_pseudo_chisqr


# The mu search starts at the first M whose log10(pseudo chi^2) is within
# CHI2_PLATEAU_DECADES of the lowest value over the next CHI2_PLATEAU_WINDOW
# M. Below that point the tau grid is too coarse and mu dips from
# discretization, not from fitting noise (doc/KRAMERS_KRONIG_REVIEW.md 1.1).
# 0.3 = a factor of 2 in chi^2, RMS residuals within sqrt(2) of what more
# elements reach; 1.0 was measured to stop too early on 1 % noise.
# The window is local on purpose: on drifting (non-KK) data chi^2 keeps
# creeping down as negative R_k absorb the drift, and the minimum over all M
# pushed the start to M ~ 45-49 in 2 of 6 synthetic drift spectra, hiding a
# 10 % drift. It must still span the chi^2 bumps of a coarse grid (a few M):
# 5 stopped early on 1 % noise, 8 matched the all-M minimum on every
# KK-compliant spectrum tested.
# The condition must hold for CHI2_PLATEAU_RUN consecutive M: on a coarse grid
# a single M can land in a chi^2 dip that the next 8 M do not undercut. An
# exact 3-RC spectrum (its RC in 2 of 10 decades) dips at M = 12 and drops
# 6 more decades after M = 19; with 0.1 % noise a single-M plateau stopped at
# M = 21 and estimated 1 % noise. 2 still stopped early on one 1 % noise seed,
# 3 and 4 gave identical results; real example spectra and drift spectra do
# not change, except that the noise estimate on drift-free data gets closer
# to the true 0.2 %.
CHI2_PLATEAU_DECADES = 0.3
CHI2_PLATEAU_WINDOW = 8
CHI2_PLATEAU_RUN = 3

# find_optimal_extend_decades grid: 11 points are 0.1-decade steps over the
# default 0-1 decade range. Candidates whose chi^2 is within EXTEND_CHI2_TIE
# (relative) of the lowest are a tie, and the smallest extension wins: it
# departs least from the standard, unextended Lin-KK grid. 0.1 % is a
# tie-breaker for numerically equal chi^2 (an extension that adds nothing),
# not a tuned value.
EXTEND_N_EVALUATIONS = 11
EXTEND_CHI2_TIE = 1e-3


def reconstruct_impedance(
    frequencies: NDArray[np.float64],
    elements: NDArray[np.float64],
    tau: NDArray[np.float64],
    L_value: Optional[float],
    include_L: bool = True,
    C_value: Optional[float] = None
) -> NDArray[np.complexfloating]:
    """
    Reconstruct impedance from fitted Voigt elements.

    Parameters
    ----------
    frequencies : array
        Frequencies [Hz]
    elements : array
        Fitted elements [R_s, R_1, ..., R_M, (L)]
    tau : array
        Time constants [s]
    L_value : float or None
        Inductance value [H]
    include_L : bool
        Whether inductance is included in elements array
    C_value : float or None, optional
        Series capacitance [F]; adds 1/(j*omega*C). Unlike L, it is never
        part of the elements array (default: None = no series C term).

    Returns
    -------
    Z_fit : array
        Reconstructed complex impedance

    Raises
    ------
    ValueError
        If elements does not hold one R per tau (plus R_s and, with
        include_L, L). zip() would silently drop the surplus, e.g. the last
        R_k when include_L=True is passed for an array without L.
    """
    omega = 2 * np.pi * frequencies
    R_s = elements[0]
    R_i_end = -1 if include_L else len(elements)
    R_i = elements[1:R_i_end]
    if len(R_i) != len(tau):
        raise ValueError(f"reconstruct_impedance: {len(R_i)} resistances for {len(tau)} "
                         f"time constants (include_L={include_L}, {len(elements)} elements)")

    Z_fit = np.full_like(frequencies, R_s, dtype=np.complex128)
    for r, t in zip(R_i, tau):
        Z_fit += r / (1 + 1j * omega * t)

    if include_L and L_value is not None:
        Z_fit += 1j * omega * L_value

    if C_value is not None:
        Z_fit += 1.0 / (1j * omega * C_value)

    return Z_fit


def find_optimal_extend_decades(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    M: int,
    search_range: Tuple[float, float] = (0.0, 1.0),
    n_evaluations: int = EXTEND_N_EVALUATIONS,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'real',
    weighting: str = 'modulus',
    min_mu: float = -np.inf
) -> Optional[Tuple[float, float, NDArray[np.float64], NDArray[np.float64], Optional[float], Optional[float]]]:
    """
    Find optimal extend_decades that minimizes pseudo chi-squared.

    Uses grid search over the specified range.

    Parameters
    ----------
    frequencies : array
        Measured frequencies [Hz]
    Z : array
        Measured impedance
    M : int
        Number of Voigt elements
    search_range : tuple
        Min and max extend_decades to search
    n_evaluations : int
        Number of grid points
    include_L : bool
        Include series inductance
    include_C : bool
        Include series capacitance (blocking low-frequency behavior)
    fit_type : str
        Fit type ('real', 'imag', 'complex')
    weighting : str
        Weighting scheme
    min_mu : float
        Reject candidates whose own mu falls below this value. Lin-KK passes
        its stop mu: at fixed M, a wider tau grid can fit with oscillating
        negative R_i, the overfit the mu criterion exists to prevent.
        -inf (default) keeps every candidate.

    Returns
    -------
    None
        Only with min_mu, when every candidate falls below it.
    optimal_extend_decades : float
        Value that minimizes chi^2
    min_chi2 : float
        Minimum chi^2 achieved
    tau : array
        Time constants for optimal extend_decades
    elements : array
        Fitted elements for optimal extend_decades
    L_value : float or None
        Inductance for optimal extend_decades
    C_value : float or None
        Series capacitance for optimal extend_decades
    """
    from ..fitting.voigt_chain import generate_tau_grid_fixed_M, estimate_R_linear, calc_mu

    candidates = np.linspace(search_range[0], search_range[1], n_evaluations)
    results = []

    for ext_dec in candidates:
        tau = generate_tau_grid_fixed_M(frequencies, M, extend_decades=ext_dec)
        elements, residual, L_value, C_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=True, include_L=include_L, include_C=include_C,
            fit_type=fit_type, allow_negative=True,
            weighting=weighting
        )

        # R_1..R_M sit right after R_s whether or not L follows them
        if calc_mu(elements[1:1 + len(tau)]) < min_mu:
            continue

        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)
        chi2 = compute_pseudo_chisqr(Z, Z_fit)
        results.append((ext_dec, chi2, tau, elements, L_value, C_value))

    if not results:
        return None

    # Find minimum chi^2
    min_chi2 = min(r[1] for r in results)
    tolerance = EXTEND_CHI2_TIE * min_chi2
    near_optimal = [r for r in results if r[1] <= min_chi2 + tolerance]
    best = min(near_optimal, key=lambda x: abs(x[0]))
    return best[0], best[1], best[2], best[3], best[4], best[5]


def _chi2_lower_M(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    max_M: int,
    include_L: bool,
    include_C: bool,
    fit_type: str,
    weighting: str
) -> int:
    """First M from which CHI2_PLATEAU_RUN consecutive M have pseudo chi^2
    within CHI2_PLATEAU_DECADES of the lowest over their next
    CHI2_PLATEAU_WINDOW M (unextended tau grid, as the mu search)."""
    from ..fitting.voigt_chain import generate_tau_grid_fixed_M, estimate_R_linear

    Ms = np.arange(3, max_M + 1)
    log_chi2 = np.empty(len(Ms))
    for i, M in enumerate(Ms):
        tau = generate_tau_grid_fixed_M(frequencies, int(M))
        elements, _, L_value, C_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=True, include_L=include_L, include_C=include_C,
            fit_type=fit_type, allow_negative=True, weighting=weighting
        )
        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)
        # Exact synthetic data can fit to chi^2 = 0
        log_chi2[i] = np.log10(np.maximum(compute_pseudo_chisqr(Z, Z_fit), np.finfo(float).tiny))

    # A non-finite chi^2 makes the comparisons False: with no plateau run,
    # start where the original Lin-KK does
    level = [log_chi2[i] - log_chi2[i:i + CHI2_PLATEAU_WINDOW + 1].min() <= CHI2_PLATEAU_DECADES
             for i in range(len(Ms))]
    return int(next((M for i, M in enumerate(Ms) if all(level[i:i + CHI2_PLATEAU_RUN])),
                    Ms[0]))
