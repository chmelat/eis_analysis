"""
Pseudo chi-squared of a model against measured impedance, and the noise
level it implies.

Not specific to Kramers-Kronig: Lin-KK, Z-HIT, the circuit fit
(fitting.diagnostics) and the DRT peak test (drt.significance) all turn
their residuals into a noise estimate with the same two functions, so that
the numbers compare.
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
