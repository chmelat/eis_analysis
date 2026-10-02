"""
Z-HIT (Z-Hilbert Impedance Transform) validation for EIS data quality assessment.

Provides non-parametric K-K validation using numerical integration:
1. zhit_reconstruct_magnitude() - Core magnitude reconstruction from phase
2. zhit_validation() - High-level wrapper (plot it with
   visualization.plot_zhit_validation)
3. ZHITResult - Dataclass with validation results

Implementation notes
--------------------
The Z-HIT method reconstructs |Z| from phase using the Kramers-Kronig relation:

    ln|Z(omega)| = C + (2/pi) * H[phi(omega)]

where H is the Hilbert transform. Two computational approaches exist:

1. FFT-based Hilbert transform (scipy.signal.hilbert)
   - Fast but sensitive to edge effects
   - Requires padding to mitigate boundary artifacts

2. Direct numerical integration in log-omega space (used here)
   - Ehm et al. (2001) showed that in log-omega space:
     H[phi] ~ integral of phi * d(ln omega) + correction term
   - More robust for finite frequency ranges typical in EIS
   - No padding required, simpler implementation

We use approach (2) because EIS data has limited frequency range where
edge effects from FFT-based Hilbert transform can distort results.

References
----------
Ehm, W. et al. (2001) "The evaluation of electrochemical impedance spectra
using a modified logarithmic Hilbert transform."
Journal of Electroanalytical Chemistry 499, 216-225
"""

import numpy as np
import logging
from dataclasses import dataclass
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid

from .kramers_kronig import compute_pseudo_chisqr, estimate_noise_percent

logger = logging.getLogger(__name__)

# The second-order term takes d(phi)/d(ln omega) between points at least 5 %
# apart in frequency. A smaller step - a frequency measured twice (an up/down
# sweep, several sweeps in one file) or a very dense sweep - turns phase noise
# into a spike: 1e-3 rad over a 1e-4 step is 10 rad. Grids up to ~47
# points/decade step more than 5 % and keep np.gradient's result exactly;
# the relaxations the term corrects for span about a decade, so a 10 %
# baseline costs nothing in truncation error. A step just over 5 % is used as
# is; it adds the noise of a regular 47 points/decade sweep (measured: +1.1 to
# +2.1 percentage points on the max residual at 0.1 % noise), not more.
MIN_DERIVATIVE_STEP = np.log(1.05)


def _phase_derivative(phi: NDArray[np.float64], ln_omega: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    d(phi)/d(ln omega) with np.gradient's formulas over neighbours at least
    MIN_DERIVATIVE_STEP away: second-order three-point inside, one-sided at
    the ends. Identical to np.gradient where every step exceeds it. The
    spectrum must span more than MIN_DERIVATIVE_STEP.
    """
    n = len(ln_omega)
    lower = np.searchsorted(ln_omega, ln_omega - MIN_DERIVATIVE_STEP, side='right') - 1
    upper = np.searchsorted(ln_omega, ln_omega + MIN_DERIVATIVE_STEP, side='left')
    has_lower, has_upper = lower >= 0, upper < n
    lower, upper = np.where(has_lower, lower, 0), np.where(has_upper, upper, n - 1)
    dx1, dx2 = ln_omega - ln_omega[lower], ln_omega[upper] - ln_omega

    with np.errstate(divide='ignore', invalid='ignore'):
        central = (-dx2 / (dx1 * (dx1 + dx2)) * phi[lower]
                   + (dx2 - dx1) / (dx1 * dx2) * phi
                   + dx1 / (dx2 * (dx1 + dx2)) * phi[upper])
        forward = (phi[upper] - phi) / dx2
        backward = (phi - phi[lower]) / dx1
    return np.where(has_lower & has_upper, central, np.where(has_upper, forward, backward))


def _quality_label(mean_abs_residual_mag: float) -> str:
    """Stratified label for KK/Z-HIT magnitude residuals (in percent)."""
    if mean_abs_residual_mag < 0.5:
        return "excellent"
    if mean_abs_residual_mag < 1.0:
        return "good"
    if mean_abs_residual_mag < 2.5:
        return "acceptable"
    if mean_abs_residual_mag < 5.0:
        return "marginal (check for drift/nonlinearity)"
    return "poor"


@dataclass
class ZHITResult:
    """Result of Z-HIT validation.

    Attributes
    ----------
    Z_mag_reconstructed : NDArray[np.float64]
        Reconstructed impedance magnitude [Ohm]
    Z_fit : NDArray[np.complex128]
        Reconstructed complex impedance (|Z_recon| * exp(j*phi)) [Ohm]
    residuals_mag : NDArray[np.float64]
        Magnitude residuals [%]
    residuals_real : NDArray[np.float64]
        Real part residuals (fraction, normalized by |Z|)
    residuals_imag : NDArray[np.float64]
        Imaginary part residuals (fraction, normalized by |Z|)
    pseudo_chisqr : float
        Pseudo chi-squared (Boukamp 1995)
    noise_estimate : float
        Estimated noise [%] (Yrjana & Bobacka 2024)
    quality : float
        Quality metric (0-1 scale, based on magnitude residuals)
    quality_threshold : float
        Pass/fail threshold for `is_valid` and the `quality` metric [%].
    """
    Z_mag_reconstructed: NDArray[np.float64]
    Z_fit: NDArray[np.complex128]
    residuals_mag: NDArray[np.float64]
    residuals_real: NDArray[np.float64]
    residuals_imag: NDArray[np.float64]
    pseudo_chisqr: float
    noise_estimate: float
    quality: float
    quality_threshold: float = 5.0

    @property
    def success(self) -> bool:
        """Check if Z-HIT validation completed successfully."""
        return self.Z_mag_reconstructed.size > 0

    @property
    def mean_residual_real(self) -> float:
        """Mean absolute real residual [%]."""
        return float(np.mean(np.abs(self.residuals_real)) * 100)

    @property
    def mean_residual_imag(self) -> float:
        """Mean absolute imaginary residual [%]."""
        return float(np.mean(np.abs(self.residuals_imag)) * 100)

    @property
    def mean_residual_mag(self) -> float:
        """Mean absolute magnitude residual [%]."""
        return float(np.mean(np.abs(self.residuals_mag)))

    @property
    def is_valid(self) -> bool:
        """Check if data passes validation (mean |residual_mag| < quality_threshold)."""
        return self.mean_residual_mag < self.quality_threshold

    @property
    def quality_label(self) -> str:
        """Stratified label for the magnitude residuals.

        Absolute scale calibrated for KK/Z-HIT residuals: clean reference cells
        typically sit in the 0.1-0.5% band, while >=5% indicates likely drift,
        nonlinearity, or other K-K violations.
        """
        return _quality_label(self.mean_residual_mag)


def zhit_reconstruct_magnitude(
    frequencies: NDArray[np.float64],
    phi: NDArray[np.float64],
    ln_Z_exp: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Reconstruct impedance magnitude from phase using Z-HIT method.

    Uses numerical integration in log-frequency space according to
    the modified logarithmic Hilbert transform (Ehm et al. 2001).

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz], sorted ascending; repeated or very close
        frequencies are allowed (see MIN_DERIVATIVE_STEP)
    phi : ndarray of float
        Phase angles [rad], sorted by ascending frequency
    ln_Z_exp : ndarray of float
        Experimental ln|Z|, sorted by ascending frequency. Only fixes the
        integration constant (see Notes).

    Returns
    -------
    ln_Z_reconstructed : ndarray of float
        Reconstructed ln|Z| values

    Raises
    ------
    ValueError
        If the frequencies span less than MIN_DERIVATIVE_STEP (5 %)

    Notes
    -----
    The Z-HIT formula (first order):

        ln|Z(omega)| = C + (2/pi) * integral[phi * d(ln omega)]

    Second order correction (Ehm et al. 2001):

        ln|Z(omega)| = first_order + gamma * d(phi)/d(ln omega),  gamma = -pi/6

    The phase determines ln|Z| only up to the constant C. It is set to
    median(ln_Z_exp - reconstruction) over the whole spectrum rather than by
    matching a single reference point: a single point carries its own noise and
    the local approximation error of the second-order term into every point of
    the reconstruction (1.7% mean residual on a noise-free R+RC when the point
    sits at the relaxation). The median is robust to outliers and to drift in
    fewer than half of the points; it fails when most of the spectrum is bad.
    """
    ln_omega = np.log(2 * np.pi * frequencies)
    if len(ln_omega) < 2 or ln_omega[-1] - ln_omega[0] < MIN_DERIVATIVE_STEP:
        raise ValueError(f"Z-HIT needs frequencies spanning at least "
                         f"{np.expm1(MIN_DERIVATIVE_STEP):.0%}, got {len(ln_omega)} point(s)")

    # First order: cumulative integration of phase
    ln_Z_reconstructed = cumulative_trapezoid((2.0 / np.pi) * phi, ln_omega, initial=0)

    # Second order: np.gradient's formulas over a minimum step, so repeated
    # or very close frequencies cannot spike it
    d_phi_d_ln_omega = _phase_derivative(phi, ln_omega)

    # Second-order correction coefficient gamma = -pi/6
    # Derived from Taylor expansion of the Hilbert transform kernel in log-omega space.
    # See Ehm et al. (2001) eq. 15 and Schiller et al. (2001) for derivation.
    gamma = -np.pi / 6.0  # ≈ -0.524
    ln_Z_reconstructed += gamma * d_phi_d_ln_omega

    return ln_Z_reconstructed + np.median(ln_Z_exp - ln_Z_reconstructed)


def zhit_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    quality_threshold: float = 5.0
) -> ZHITResult:
    """
    Perform Z-HIT (Z-Hilbert Impedance Transform) validation on EIS data.

    Z-HIT uses numerical integration to validate Kramers-Kronig compliance
    without model fitting (non-parametric method).

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    quality_threshold : float, optional
        Reference threshold for quality metric calculation [%].
        Default: 5.0 (5% mean residual = quality 0)

    Returns
    -------
    ZHITResult
        Dataclass containing:
        - Z_mag_reconstructed: Reconstructed impedance magnitude [Ohm]
        - Z_fit: Reconstructed complex impedance [Ohm]
        - residuals_mag: Magnitude residuals [%]
        - residuals_real: Real part residuals (fraction)
        - residuals_imag: Imaginary part residuals (fraction)
        - pseudo_chisqr: Pseudo chi-squared (Boukamp 1995)
        - noise_estimate: Estimated noise [%]
        - quality: Quality metric (0-1 scale)

    Notes
    -----
    The Z-HIT method reconstructs |Z| from phase using:

        ln|Z(omega)| = C + (2/pi) * integral[phi * d(ln omega)]
                       + gamma * d(phi)/d(ln omega)

    with C fixed by the median over the spectrum, see
    `zhit_reconstruct_magnitude`.

    Advantages over Lin-KK:
    - No model fitting required (truly non-parametric)
    - Faster computation
    - Different sensitivity to certain non-compliance types

    References
    ----------
    Ehm, W. et al. (2001) "The evaluation of electrochemical impedance spectra
    using a modified logarithmic Hilbert transform."
    Journal of Electroanalytical Chemistry 499, 216-225
    """
    # Sort by ascending frequency for the integration; remember inverse
    # permutation so output arrays match the user's original order.
    sort_idx = np.argsort(frequencies)
    inv_idx = np.argsort(sort_idx)
    frequencies = frequencies[sort_idx]
    Z = Z[sort_idx]

    # Extract magnitude and phase. np.unwrap removes 2*pi jumps from arctan2
    # at the [-pi, pi] boundary (relevant for inductive systems or noisy data
    # near the wrap point), which would otherwise spike the phase derivative.
    Z_mag = np.abs(Z)
    phi = np.unwrap(np.arctan2(Z.imag, Z.real))

    # Perform Z-HIT magnitude reconstruction using numerical integration
    try:
        ln_Z_reconstructed = zhit_reconstruct_magnitude(
            frequencies, phi, np.log(Z_mag)
        )
        Z_mag_reconstructed = np.exp(ln_Z_reconstructed)
    except Exception as e:
        logger.error(f"Z-HIT computation failed: {e}", exc_info=True)
        # Return empty result on failure
        empty = np.array([])
        return ZHITResult(
            Z_mag_reconstructed=empty,
            Z_fit=np.array([], dtype=np.complex128),
            residuals_mag=empty,
            residuals_real=empty,
            residuals_imag=empty,
            pseudo_chisqr=0.0,
            noise_estimate=0.0,
            quality=0.0,
            quality_threshold=quality_threshold
        )

    # Reconstruct complex impedance: Z_fit = |Z_recon| * exp(j*phi)
    # Phase is preserved from original data, only magnitude is reconstructed
    Z_fit = Z_mag_reconstructed * np.exp(1j * phi)

    # Calculate magnitude residuals (in %)
    residuals_mag = (Z_mag - Z_mag_reconstructed) / Z_mag * 100.0

    # Calculate complex residuals (normalized by |Z|, as fraction)
    residuals_real = (Z.real - Z_fit.real) / Z_mag
    residuals_imag = (Z.imag - Z_fit.imag) / Z_mag

    # Calculate pseudo chi-squared and noise estimate
    pseudo_chisqr = compute_pseudo_chisqr(Z, Z_fit)
    noise_estimate = estimate_noise_percent(pseudo_chisqr, len(Z))

    # Calculate quality metric (based on magnitude residuals)
    mean_abs_residual_mag = np.mean(np.abs(residuals_mag))
    quality = max(0.0, 1.0 - mean_abs_residual_mag / quality_threshold)

    # Restore the user's original frequency ordering on output arrays so they
    # pair element-wise with the input `frequencies` / `Z`.
    return ZHITResult(
        Z_mag_reconstructed=Z_mag_reconstructed[inv_idx],
        Z_fit=Z_fit[inv_idx],
        residuals_mag=residuals_mag[inv_idx],
        residuals_real=residuals_real[inv_idx],
        residuals_imag=residuals_imag[inv_idx],
        pseudo_chisqr=pseudo_chisqr,
        noise_estimate=noise_estimate,
        quality=quality,
        quality_threshold=quality_threshold
    )
