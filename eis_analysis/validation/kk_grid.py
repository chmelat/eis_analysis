"""
The Lin-KK model on its tau grid: impedance of a fitted Voigt chain, the
choice of the grid extension, and where the mu search starts.

Used by lin_kk_native (kramers_kronig.py); reconstruct_impedance also by the
KK plot.

CONTRACT

Model. A Voigt chain on a tau grid: Z(f) = R_s + sum_k R_k / (1 + j omega
tau_k) + j omega L + 1 / (j omega C), omega = 2 pi f; the L and C terms
only where the model has them. Its linear fit is
fitting.voigt_chain.estimate_R_linear(include_Rs=True, allow_negative=True)
on the tau grid of generate_tau_grid_fixed_M(frequencies, M,
extend_decades): M log-spaced tau from 1/(2 pi f_max) to 1/(2 pi f_min /
10^extend_decades), a negative extension treated as 0. The fit returns
elements [R_s, R_1..R_M, (L)], L_value (None without L) and C_value (None
without C, and where the fitted 1/C is 0); C is never in elements. mu is
fitting.voigt_chain.calc_mu(R_1..R_M), 1 - (sum of negative R_k) / (sum of
positive R_k); chi^2 is validation.noise.compute_pseudo_chisqr(Z, Z_fit) of
the fit's reconstruction.

Callers.
- validation.kramers_kronig.lin_kk_native: _chi2_lower_M(f, Z, max_M,
  include_L, include_C, fit_type, weighting) with 3 <= max_M <= N - 2
  (it raises before otherwise); the result is min_M of the mu search
  (find_optimal_M_mu, unextended grid, the same fit settings). With
  auto_extend_decades, find_optimal_extend_decades(f, Z, M, search_range,
  EXTEND_N_EVALUATIONS, ..., min_mu=mu) with M and mu of the mu search's
  stop; on None it keeps the unextended grid and warns. Then
  reconstruct_impedance of the chosen fit -> KKResult.Z_fit. f and Z are
  1-D of equal length, as the user passed them (not sorted); the fit
  rejects NaN, inf, f <= 0 and Z = 0 (validate_eis_data, ValueError).
  kramers_kronig_validation (the CLI path) calls it with
  auto_extend_decades=True, search_range (0, 1) by default or the user's,
  fit_type 'real', weighting 'modulus', L included, and turns ValueError
  and LinAlgError into KKResult.error.
- visualization.plots (KK plot): reconstruct_impedance(300 log-spaced f
  over the measured range, result.elements, result.tau,
  result.inductance, include_L=result.inductance is not None,
  C_value=result.capacitance).
Public API (eis_analysis.validation: find_optimal_extend_decades,
reconstruct_impedance): any caller, nothing guaranteed beyond the rules
below. The caller passes ndarrays.

reconstruct_impedance(frequencies, elements, tau, L_value, include_L=True,
C_value=None) -> complex128 array, the shape of frequencies.
  The model above with R_s = elements[0], R_k = elements[k]. With
  include_L the last entry of elements is the L slot: it is not an R, and
  the L term is j omega L_value (none where L_value is None); the value in
  the slot is not used. The C term where C_value is not None. ValueError
  when elements does not hold exactly 1 + len(tau) + [include_L] entries,
  an empty elements included. The caller passes f > 0 and C_value finite
  and nonzero (the fit gives None, not 0); nothing is promised otherwise.
  Linear in (elements, L_value, 1/C_value).

find_optimal_extend_decades(frequencies, Z, M, search_range=(0, 1),
n_evaluations=EXTEND_N_EVALUATIONS, include_L, include_C, fit_type,
weighting, min_mu=-inf) -> None or (extend_decades, chi2, tau, elements,
L_value, C_value).
  Candidates: n_evaluations values evenly spaced over search_range, both
  ends included (n_evaluations = 1: search_range[0] only). Each is
  fitted on its own grid; a candidate is kept when its mu >= min_mu and
  its chi^2 is finite (a fit that gave no number does not compete; mu is
  NaN only with a NaN or infinite R_k, whose chi^2 is not finite
  either; -inf, all R_k negative, falls below any finite min_mu; a NaN
  min_mu keeps no candidate).
  Among the kept ones, those with chi^2 <= chi2_min (1 +
  EXTEND_CHI2_TIE) tie, chi2_min the lowest, and the smallest extension
  of the tie wins. The tuple is that candidate's: its
  extension, its own chi^2 (within the tie of chi2_min, not chi2_min
  itself), tau, elements, L_value, C_value, so that reconstruct_impedance
  of them gives that chi^2. None when no candidate is kept; with min_mu =
  -inf only through non-finite fits. At lin_kk_native's call with
  search_range[0] = 0 the extension 0 is the mu search's stop fit itself,
  its mu equals min_mu and it is kept where its chi^2 is finite: None
  there needs search_range[0] > 0 or a stop fit without a finite chi^2
  (|Z| beyond the known limit below).
  ValueError where the fit raises (invalid data, weighting or fit_type),
  LinAlgError where its SVD does not converge. The caller passes 0 <=
  search_range[0] <= search_range[1] finite, n_evaluations >= 1, M >= 1,
  f with at least two distinct values and M not above their number;
  nothing is promised otherwise (a negative extension fits the grid of 0
  but is returned as given; equal f collapse the tau grid to one value
  and M > N fits, both return numbers). f <= 0 is the fit's ValueError,
  but the tau grid is built from f first: a RuntimeWarning (divide by
  zero, log10) comes before it.

_chi2_lower_M(frequencies, Z, max_M, include_L, include_C, fit_type,
weighting) -> int.
  For M = 3..max_M, l(M) = log10 chi^2 of the fit on the unextended grid,
  chi^2 floored at the smallest positive float (exact data fits to 0).
  M levels when l(M) - min l(M..M + CHI2_PLATEAU_WINDOW) <=
  CHI2_PLATEAU_DECADES, the window cut at max_M (no fit beyond it). The
  result is the first M from which CHI2_PLATEAU_RUN consecutive M level,
  the run cut at max_M as well: max_M with a finite chi^2 always levels
  and forms a run on its own, so the result lies in [3, max_M]. Where
  chi^2 falls up to max_M (spectra of up to ~20 points, max_M = N - 2)
  the start is near max_M; a full run would fall back to 3 there, where
  the coarse grid stops the mu search early (15 KK-compliant points at
  0.1 % noise: start 13, noise 0.8 %, valid, against start 3, 12.8 %,
  not valid). A non-finite chi^2 (NaN, inf) never levels, and a NaN in
  M's window keeps M from levelling: a failed fit is no plateau. Only
  where no run is left (non-finite chi^2 at max_M), 3, where the
  original Lin-KK starts. Errors as for find_optimal_extend_decades. The
  caller passes max_M >= 3.

What must hold:
  - The inputs are not changed. The results are deterministic.
  - Raises only ValueError and LinAlgError (the wrapper catches these),
    and only where stated; no RuntimeWarning on what the callers pass
    (f > 0: the loaders drop other points).
  - Units: Z scaled by c > 0 leaves the chosen extension, M, chi^2 and mu
    unchanged (to rounding) and scales elements and L_value by c,
    C_value by 1/c, for every weighting (the fit normalizes the weights),
    while |Z| stays within 1e-150..1e150 (the known limit,
    validate_eis_data). Above it chi^2 is not to be trusted: once a
    component of Z exceeds ~1.3e154 (the square root of the largest
    float) its square overflows in chi^2, the weight 1/|Z|^2 turns 0 and
    chi^2 drops (noise underestimated, down to 0 % with the KK test
    passing) or turns NaN where a residual overflows too; which one
    depends on the data, not on a fixed scale. Frequencies scaled
    by c leave them unchanged and scale tau, L_value and C_value by 1/c
    for fit_type 'real' and 'complex' over measured ranges (with C,
    'real' held to c = 1e3, not at 1e6); for 'imag' this rests on the
    conditioning of the fit and is not promised (with C, f x 1e3 moved
    the extension 0.0 -> 0.2 and chi^2 by a factor 2.4).
  - Permuting the points (f and Z together) changes no result (to
    rounding); Z_fit is permuted with them.
  - "To rounding" holds for the numbers, not for the choices made on
    them: where chi^2 is at rounding level (exact data the model
    represents, a pure R or R-L: 1e-31 to 1e-28), the tie between
    extensions and the plateau over M are decided by rounding, and the
    extension and M change with the point order and the scale of Z. The
    noise estimate stays at ~1e-14 %.
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
    # Counted before indexing: elements[0] of an empty array was an
    # IndexError, and with include_L elements[1:-1] of one entry matched tau = []
    if len(elements) != 1 + len(tau) + int(include_L):
        raise ValueError(f"reconstruct_impedance: {len(elements)} elements for {len(tau)} "
                         f"time constants (include_L={include_L})")
    omega = 2 * np.pi * frequencies
    R_s = elements[0]
    R_i = elements[1:1 + len(tau)]

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
        -inf (default) keeps every candidate with a finite chi^2.

    Returns
    -------
    None
        When every candidate falls below min_mu or gives no finite chi^2.
    optimal_extend_decades : float
        Smallest extension whose chi^2 is within EXTEND_CHI2_TIE of the
        lowest
    min_chi2 : float
        chi^2 of that extension (not necessarily the lowest)
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
        # `not >=`: a NaN stop mu keeps no candidate (`< NaN` kept all)
        if not calc_mu(elements[1:1 + len(tau)]) >= min_mu:
            continue

        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)
        chi2 = compute_pseudo_chisqr(Z, Z_fit)
        # A fit that gave no number does not compete: a NaN made min() and
        # the tie below depend on the candidate order
        if not np.isfinite(chi2):
            continue
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

    # A non-finite chi^2 falsifies the comparisons it enters (as itself or as
    # the window minimum); where no plateau run is left, start where the
    # original Lin-KK does
    level = [log_chi2[i] - log_chi2[i:i + CHI2_PLATEAU_WINDOW + 1].min() <= CHI2_PLATEAU_DECADES
             for i in range(len(Ms))]
    return int(next((M for i, M in enumerate(Ms) if all(level[i:i + CHI2_PLATEAU_RUN])),
                    Ms[0]))
