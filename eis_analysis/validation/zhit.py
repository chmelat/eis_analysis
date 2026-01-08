"""
Z-HIT (Z-Hilbert Impedance Transform) validation for EIS data quality assessment.

Provides non-parametric K-K validation using numerical integration:
1. zhit_reconstruct_magnitude() - Core magnitude reconstruction from phase
2. zhit_validation() - High-level wrapper with visualization
3. ZHITResult - Dataclass with validation results

Implementation notes
--------------------
The Z-HIT method reconstructs |Z| from phase using the Kramers-Kronig relation:

    ln|Z(omega)| = ln|Z(omega_ref)| + (2/pi) * H[phi(omega)]

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
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass
from typing import Optional
from numpy.typing import NDArray

from .kramers_kronig import compute_pseudo_chisqr, estimate_noise_percent

logger = logging.getLogger(__name__)


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
    ref_freq : float
        Reference frequency used [Hz]
    figure : Optional[plt.Figure]
        Visualization figure
    """
    Z_mag_reconstructed: NDArray[np.float64]
    Z_fit: NDArray[np.complex128]
    residuals_mag: NDArray[np.float64]
    residuals_real: NDArray[np.float64]
    residuals_imag: NDArray[np.float64]
    pseudo_chisqr: float
    noise_estimate: float
    quality: float
    ref_freq: float
    figure: Optional[plt.Figure] = None

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
        """Check if data passes validation (magnitude residuals < 5%)."""
        return self.mean_residual_mag < 5.0


def _calculate_offset_weighted(
    ln_Z_fit: NDArray[np.float64],
    ln_Z_exp: NDArray[np.float64],
    frequencies: NDArray[np.float64],
    center: float = 1.5,
    width: float = 3.0
) -> float:
    """
    Calculate optimal offset using weighted least-squares.

    Uses a Gaussian window in log-frequency space to weight the fit,
    focusing on a specific frequency range (default: 1-1000 Hz).

    Parameters
    ----------
    ln_Z_fit : array
        Reconstructed ln|Z| (before offset adjustment)
    ln_Z_exp : array
        Experimental ln|Z|
    frequencies : array
        Frequencies [Hz]
    center : float
        Center of weight window on log10(f) scale (default: 1.5 → ~31.6 Hz)
    width : float
        Width of weight window in decades (default: 3.0)

    Returns
    -------
    offset : float
        Optimal offset to add to ln_Z_fit
    """
    log_f = np.log10(frequencies)

    # Gaussian window: ~95% of weight within 'width' decades
    sigma = width / 4
    weights = np.exp(-0.5 * ((log_f - center) / sigma)**2)

    # Clip weights to [0, 1]
    weights = np.clip(weights, 0.0, 1.0)

    # Analytical weighted least-squares solution
    diff = ln_Z_exp - ln_Z_fit
    offset = float(np.sum(weights * diff) / np.sum(weights))

    return offset


def zhit_reconstruct_magnitude(
    frequencies: NDArray[np.float64],
    phi: NDArray[np.float64],
    ln_Z_ref: float,
    ref_idx: int,
    use_second_order: bool = True,
    optimize_offset: bool = False,
    ln_Z_exp: Optional[NDArray[np.float64]] = None,
    offset_center: float = 1.5,
    offset_width: float = 3.0
) -> NDArray[np.float64]:
    """
    Reconstruct impedance magnitude from phase using Z-HIT method.

    Uses numerical integration in log-frequency space according to
    the modified logarithmic Hilbert transform (Ehm et al. 2001).

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz], sorted ascending
    phi : ndarray of float
        Phase angles [rad], sorted by ascending frequency
    ln_Z_ref : float
        Natural logarithm of impedance magnitude at reference frequency
    ref_idx : int
        Index of reference frequency in the data
    use_second_order : bool, optional
        Use second-order correction term (derivative of phase). Default: True
    optimize_offset : bool, optional
        Use weighted least-squares offset optimization instead of fixed reference
        point. Default: False
    ln_Z_exp : array, optional
        Experimental ln|Z| values, required if optimize_offset=True
    offset_center : float, optional
        Center of weight window on log10(f) scale (default: 1.5 → ~31.6 Hz)
    offset_width : float, optional
        Width of weight window in decades (default: 3.0)

    Returns
    -------
    ln_Z_reconstructed : ndarray of float
        Reconstructed ln|Z| values

    Notes
    -----
    The Z-HIT formula (first order):

        ln|Z(omega)| = ln|Z(omega_ref)| + (2/pi) * integral[phi * d(ln omega)]

    Second order correction (Ehm et al. 2001):

        ln|Z(omega)| = first_order - gamma * d(phi)/d(ln omega)

    where gamma is a weighting factor (typically ~0.2-0.5).
    """
    n = len(phi)
    ln_omega = np.log(2 * np.pi * frequencies)

    # First order: cumulative integration of phase
    # Use cumulative trapezoidal integration
    ln_Z_first_order = np.zeros(n)
    ln_Z_first_order[ref_idx] = ln_Z_ref

    # Integrate forward from reference point (ref_idx to end)
    for i in range(ref_idx + 1, n):
        d_ln_omega = ln_omega[i] - ln_omega[i - 1]
        phi_avg = 0.5 * (phi[i] + phi[i - 1])
        ln_Z_first_order[i] = ln_Z_first_order[i - 1] + (2.0 / np.pi) * phi_avg * d_ln_omega

    # Integrate backward from reference point (ref_idx to start)
    for i in range(ref_idx - 1, -1, -1):
        d_ln_omega = ln_omega[i + 1] - ln_omega[i]
        phi_avg = 0.5 * (phi[i] + phi[i + 1])
        ln_Z_first_order[i] = ln_Z_first_order[i + 1] - (2.0 / np.pi) * phi_avg * d_ln_omega

    if not use_second_order:
        return ln_Z_first_order

    # Second order correction: subtract weighted derivative of phase
    # d(phi)/d(ln omega) computed via central differences
    d_phi_d_ln_omega = np.zeros(n)

    # Central differences for interior points
    for i in range(1, n - 1):
        d_phi_d_ln_omega[i] = (phi[i + 1] - phi[i - 1]) / (ln_omega[i + 1] - ln_omega[i - 1])

    # Forward/backward differences at boundaries
    d_phi_d_ln_omega[0] = (phi[1] - phi[0]) / (ln_omega[1] - ln_omega[0])
    d_phi_d_ln_omega[-1] = (phi[-1] - phi[-2]) / (ln_omega[-1] - ln_omega[-2])

    # Second-order correction coefficient gamma = -pi/6
    # Derived from Taylor expansion of the Hilbert transform kernel in log-omega space.
    # See Ehm et al. (2001) eq. 15 and Schiller et al. (2001) for derivation.
    # Formula: ln|Z| = (2/pi) * integral + gamma * d(phi)/d(ln omega)
    gamma = -np.pi / 6.0  # ≈ -0.524

    ln_Z_reconstructed = ln_Z_first_order + gamma * d_phi_d_ln_omega

    # Calculate offset
    if optimize_offset and ln_Z_exp is not None:
        # Weighted least-squares optimization
        offset = _calculate_offset_weighted(
            ln_Z_reconstructed, ln_Z_exp, frequencies,
            center=offset_center, width=offset_width
        )
    else:
        # Original behavior: fixed offset at reference point
        offset = ln_Z_ref - ln_Z_reconstructed[ref_idx]

    ln_Z_reconstructed += offset

    return ln_Z_reconstructed


def zhit_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    ref_freq: Optional[float] = None,
    quality_threshold: float = 5.0,
    optimize_offset: bool = False,
    offset_center: float = 1.5,
    offset_width: float = 3.0
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
    ref_freq : float, optional
        Reference frequency for magnitude anchor [Hz].
        Default: geometric mean of frequency range.
    quality_threshold : float, optional
        Reference threshold for quality metric calculation [%].
        Default: 5.0 (5% mean residual = quality 0)
    optimize_offset : bool, optional
        Use weighted least-squares offset optimization instead of fixed reference
        point. Default: False
    offset_center : float, optional
        Center of weight window on log10(f) scale (default: 1.5 → ~31.6 Hz)
    offset_width : float, optional
        Width of weight window in decades (default: 3.0)

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
        - ref_freq: Reference frequency used [Hz]
        - figure: Visualization figure

    Notes
    -----
    The Z-HIT method reconstructs |Z| from phase using:

        ln|Z(omega)| = ln|Z(omega_ref)| + (2/pi) * H[phi(omega)]

    where H is the Hilbert transform operator.

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
    logger.info("=" * 60)
    logger.info("Z-HIT validation")
    logger.info("=" * 60)

    # Ensure data is sorted by ascending frequency
    sort_idx = np.argsort(frequencies)
    frequencies = frequencies[sort_idx]
    Z = Z[sort_idx]

    # Extract magnitude and phase
    Z_mag = np.abs(Z)
    phi = np.arctan2(Z.imag, Z.real)  # Phase in radians

    # Select reference frequency (geometric mean by default)
    if ref_freq is None:
        ref_freq = np.sqrt(frequencies[0] * frequencies[-1])

    # Find index closest to reference frequency
    ref_idx = np.argmin(np.abs(frequencies - ref_freq))
    actual_ref_freq = frequencies[ref_idx]
    ln_Z_ref = np.log(Z_mag[ref_idx])

    # Perform Z-HIT magnitude reconstruction using numerical integration
    try:
        ln_Z_reconstructed = zhit_reconstruct_magnitude(
            frequencies, phi, ln_Z_ref, ref_idx,
            optimize_offset=optimize_offset,
            ln_Z_exp=np.log(Z_mag),
            offset_center=offset_center,
            offset_width=offset_width
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
            ref_freq=actual_ref_freq,
            figure=None
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

    # Logging (format consistent with KK validation)
    logger.info(f"Z-HIT: ref_freq={actual_ref_freq:.2e} Hz")
    logger.info(f"  Mean |res_real|: {np.mean(np.abs(residuals_real))*100:.2f}%")
    logger.info(f"  Mean |res_imag|: {np.mean(np.abs(residuals_imag))*100:.2f}%")
    logger.info(f"  Pseudo chi^2: {pseudo_chisqr:.2e}")
    logger.info(f"  Estimated noise (upper bound): {noise_estimate:.2f}%")

    # Quality assessment
    if mean_abs_residual_mag < 5.0:
        logger.info("Data quality is good (residuals < 5%)")
    else:
        logger.warning("Data may contain artifacts (residuals >= 5%)")

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Left panel: Magnitude comparison (log-log)
    ax1 = axes[0]
    ax1.loglog(frequencies, Z_mag, 'o', label='Measured', markersize=4)
    ax1.loglog(frequencies, Z_mag_reconstructed, '-', label='Z-HIT reconstruction',
               linewidth=2, color='red')
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("|Z| [Ohm]")
    ax1.set_title("Z-HIT validation")
    ax1.legend()
    ax1.grid(True, alpha=0.3, which='both')

    # Mark reference frequency
    ax1.axvline(x=actual_ref_freq, color='gray', linestyle='--', alpha=0.5)

    # Right panel: Complex residuals
    ax2 = axes[1]
    ax2.semilogx(frequencies, residuals_real * 100, 'o', label='Real', markersize=4, alpha=0.7)
    ax2.semilogx(frequencies, residuals_imag * 100, 's', label='Imag', markersize=4, alpha=0.7)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(y=5, color='r', linestyle=':', alpha=0.5)
    ax2.axhline(y=-5, color='r', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Residuals [%]")
    ax2.set_title(f"Z-HIT residuals (χ²={pseudo_chisqr:.2e}, noise≤{noise_estimate:.1f}%)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return ZHITResult(
        Z_mag_reconstructed=Z_mag_reconstructed,
        Z_fit=Z_fit,
        residuals_mag=residuals_mag,
        residuals_real=residuals_real,
        residuals_imag=residuals_imag,
        pseudo_chisqr=pseudo_chisqr,
        noise_estimate=noise_estimate,
        quality=quality,
        ref_freq=actual_ref_freq,
        figure=fig
    )
