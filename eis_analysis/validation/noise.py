"""
Pseudo chi-squared of a model against measured impedance, and the noise
level it implies.

Not specific to Kramers-Kronig: Lin-KK, Z-HIT, the circuit fit
(fitting.diagnostics) and the DRT peak test (drt.significance) all turn
their residuals into a noise estimate with the same two functions, so that
the numbers compare.

CONTRACT

Noise model. Measured Z = Z_true + e, the real and imaginary part of e
independent N(0, (sigma |Z|)^2), sigma the same relative level for every
point and both components. The residual of a point, normalized by its
|Z_exp|, then has components ~ N(0, sigma^2) each. Every estimate here
returns sigma in percent; a model misfit or a KK violation adds to the
residuals, so each is an upper bound on the measurement noise in
expectation. A single estimate scatters around sigma: with the true model
as Z_fit it lies below sigma in about half of the spectra, and a fit with
many parameters absorbs part of the noise.

Callers. All pass ndarrays.
- validation.kramers_kronig.lin_kk_native: compute_pseudo_chisqr(Z, Z_fit)
  over the measured points and the Voigt chain's reconstruction, then
  estimate_noise_percent(chi2, len(Z)) -> KKResult.noise_estimate. Its
  linear fit rejects Z = 0, NaN and inf before (validate_eis_data,
  ValueError).
- validation.zhit: the same with Z-HIT's Z_fit (reconstructed |Z|, measured
  phase) -> ZHITResult.noise_estimate. Z = 0, NaN and inf are rejected
  before (validate_eis_data; empty result).
- validation.kk_grid: compute_pseudo_chisqr per Lin-KK candidate (extension
  choice, chi^2 plateau over M); compares the values and takes log10 with a
  floor at the smallest positive float. Z = 0 never reaches it (rejected
  by the linear fit); a chi^2 inf everywhere sends it to the smallest M.
- validation.kk_result.KKResult: noise_estimate_on(mask) passes
  sum(r_re^2 + r_im^2) over the selected residuals and their count to
  estimate_noise_percent, and promises it equals noise_estimate with every
  point selected. It holds for |Z| within 1e-150..1e150 (below, |Z|^2 here
  underflows, |Z| there not): both divide by |Z|, and Lin-KK rejects Z = 0.
  noise_estimate_median passes the stored residuals.
  validation.outliers compares noise_estimate_median with
  NOISE_FLOOR_PERCENT; lin_kk_native warns when noise_estimate exceeds
  NOISE_DOMINATED_RATIO x the median (a few points dominate).
- drt.significance: compute_pseudo_chisqr(Z, DRT reconstruction) for the
  full and the constrained solution; Delta chi^2 is scaled by default by
  sigma = sqrt(chi2 / 2N), i.e. estimate_noise_percent / 100, floored at
  machine eps (a noise_sigma passed in replaces it: the lambda probe
  passes the main run's). calculate_drt rejects Z = 0 before
  (validate_eis_data; empty result).
- fitting.diagnostics.compare_fit_to_noise: estimate_noise_percent of the
  fit's pseudo chi^2 against the Lin-KK noise; non-empty, equal lengths
  checked there, a non-finite result reported, not compared.
Public API (eis_analysis.validation): any caller, nothing guaranteed beyond
the rules below. The caller passes ndarrays (a list fails in
compute_pseudo_chisqr) of equal shape, 1-D for the median; the shapes are
not checked, other ones broadcast to a meaningless number.

compute_pseudo_chisqr(Z_exp, Z_fit) -> float.
  sum over points of |Z_exp - Z_fit|^2 / |Z_exp|^2 (Boukamp 1995), >= 0.
  0 for identical arrays without a zero and for empty ones (NaN with one). Unchanged (to
  rounding) when both are scaled by the same nonzero complex factor,
  conjugated, or permuted together, as long as |Z| and |Z_exp - Z_fit|
  stay within 1e-150..1e150 (squares are formed and turn subnormal below;
  impedances span about 1e-6..1e13 Ohm). A point with Z_exp = 0 has no
  relative residual: the result is inf, NaN where Z_fit = 0 there too (inf
  weight x 0), with a RuntimeWarning; rejecting Z = 0 is the caller's job
  (Lin-KK, Z-HIT and the DRT do). A NaN
  or inf in either array gives NaN or inf, never a finite number that
  leaves the point out.

estimate_noise_percent(chi2_ps, n_points) -> float.
  100 sqrt(chi2_ps / (2 n_points)): with the noise model, E[chi2_ps] =
  2 N sigma^2. inf for chi2_ps = inf, NaN for NaN. ValueError when
  n_points <= 0. The caller passes a positive integer n_points and
  chi2_ps >= 0 (or inf, NaN); nothing is promised otherwise (n_points =
  inf gives 0, a negative chi2_ps NaN with a RuntimeWarning).

estimate_noise_percent_median(residuals_real, residuals_imag) -> float.
  100 median(hypot(r_re, r_im)) / sqrt(2 ln 2), the residuals normalized by
  |Z| as fractions: |r| is Rayleigh with median sigma sqrt(2 ln 2), so under
  the noise model it estimates the same sigma as estimate_noise_percent.
  Robust, not immune: k large residuals among N points move it by about
  k/2 order statistics of the others, and it stays bounded by their values
  while k < N/2 (5 of 81 points off by ~2.5 % at 0.1 % noise: 0.20 %
  against 0.70 % from the mean square). nan for empty residuals and where
  a residual is NaN (an unfitted point is not ignored), except a point
  whose other component is inf: hypot(inf, nan) = inf counts as a large
  residual. Invariant under permuting the points (together) and under a
  sign change of either component.
"""

import numpy as np
from numpy.typing import NDArray


def compute_pseudo_chisqr(
    Z_exp: NDArray[np.complexfloating],
    Z_fit: NDArray[np.complexfloating]
) -> float:
    """
    Compute pseudo chi-squared (Boukamp 1995).

    Parameters
    ----------
    Z_exp : array
        Experimental impedance
    Z_fit : array
        Fitted impedance

    Returns
    -------
    float
        Pseudo chi-squared value

    References
    ----------
    Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance
    Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
    """
    weight = 1.0 / (Z_exp.real**2 + Z_exp.imag**2)
    return float(np.sum(weight * (
        (Z_exp.real - Z_fit.real)**2 +
        (Z_exp.imag - Z_fit.imag)**2
    )))


def estimate_noise_percent(chi2_ps: float, n_points: int) -> float:
    """
    Estimate noise standard deviation from pseudo chi-squared.

    Based on Yrjana & Bobacka (2024). With the same relative noise sigma in
    the real and the imaginary part, chi2_ps ~ 2 * N * sigma^2, so
    sigma [%] = 100 * sqrt(chi2_ps / (2N)) = sqrt(chi2_ps * 5000 / N).

    It is an upper bound: chi2_ps also contains the misfit of the model and
    any KK violation, not only the measurement noise.

    Parameters
    ----------
    chi2_ps : float
        Pseudo chi-squared value
    n_points : int
        Number of data points (must be positive)

    Returns
    -------
    float
        Estimated noise in percent

    Raises
    ------
    ValueError
        If n_points is not positive

    References
    ----------
    Yrjana, V. and Bobacka, J. "Implementing Kramers-Kronig validity testing
    using pyimpspec." Electrochim. Acta 504, 144951 (2024)
    """
    if n_points <= 0:
        raise ValueError(f"estimate_noise_percent: n_points must be positive, got {n_points}")
    return float(np.sqrt(chi2_ps * 5000 / n_points))


# |r| = hypot(re, im) of two independent N(0, sigma) components is Rayleigh
# distributed, with median sigma * sqrt(2 ln 2)
RAYLEIGH_MEDIAN = float(np.sqrt(2 * np.log(2)))


def estimate_noise_percent_median(
    residuals_real: NDArray[np.float64],
    residuals_imag: NDArray[np.float64]
) -> float:
    """
    Noise standard deviation [%] from the median residual.

    The same sigma per component as estimate_noise_percent, but from the
    median of |r| instead of the mean square, so a few large residuals do
    not move it. On a spectrum with 5 of 81 points off by ~2.5 % at 0.1 %
    noise the mean square gives 0.70 %, the median 0.20 %.

    Parameters
    ----------
    residuals_real, residuals_imag : ndarray
        Residuals normalized by |Z| (fractions)

    Returns
    -------
    float
        Estimated noise in percent (nan for empty residuals)
    """
    if len(residuals_real) == 0:
        return float('nan')
    return float(100 * np.median(np.hypot(residuals_real, residuals_imag)) / RAYLEIGH_MEDIAN)
