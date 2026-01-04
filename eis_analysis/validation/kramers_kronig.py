"""
Kramers-Kronig validation for EIS data quality assessment.

Provides two implementations:
1. lin_kk_native() - Native implementation using Voigt chain (no external dependencies)
2. kramers_kronig_validation() - High-level wrapper with visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from typing import Tuple, Optional
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'proportional'
) -> Tuple[int, float, NDArray[np.complex128], NDArray[np.float64], NDArray[np.float64], Optional[float]]:
    """
    Native Lin-KK implementation using Voigt chain fitting.

    Implements the linear Kramers-Kronig test from Schönleber et al. (2014)
    without external dependencies (replaces impedance.py linKK).

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ω]
    mu_threshold : float, optional
        Threshold for μ metric (default: 0.85)
        Lower values → more conservative (fewer elements)
    max_M : int, optional
        Maximum number of Voigt elements to try (default: 50)
    include_L : bool, optional
        Include series inductance L (default: True)
    fit_type : str, optional
        Fit type: 'real' (default, impedance.py compatible), 'imag', or 'complex'
    weighting : str, optional
        Weighting scheme: 'uniform', 'sqrt', 'proportional', 'square'
        Default: 'proportional' (Lin-KK standard, w=1/|Z|)

    Returns
    -------
    M : int
        Optimal number of Voigt elements
    mu : float
        Final μ metric value (1.0 = no overfit, <0.85 = stop criterion reached)
    Z_fit : ndarray of complex
        Fitted impedance [Ω]
    res_real : ndarray of float
        Real part residuals normalized by |Z| (fractional, not %)
    res_imag : ndarray of float
        Imaginary part residuals normalized by |Z| (fractional, not %)
    L_value : float or None
        Fitted inductance [H] if include_L=True

    Notes
    -----
    Algorithm:
    1. Start with M=3 Voigt elements
    2. Fit using pseudoinverse (allows negative R for overfit detection)
    3. Calculate μ = 1 - Σ|R_negative|/Σ|R_positive|
    4. If μ > threshold, increase M and repeat
    5. When μ ≤ threshold, optimal M found

    References
    ----------
    Schönleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20–27 (2014)
    """
    from ..fitting.voigt_chain import (
        find_optimal_M_mu
    )

    # Find optimal M using μ metric
    M, mu, tau, elements, L_value = find_optimal_M_mu(
        frequencies, Z,
        mu_threshold=mu_threshold,
        max_M=max_M,
        extend_decades=0.0,  # Lin-KK standard: no extension
        include_Rs=True,
        include_L=include_L,
        fit_type=fit_type,
        allow_negative=True,  # Required for μ metric
        weighting=weighting
    )

    # Extract parameters from elements array
    # elements = [R_s, R_1, ..., R_M, L (optional)]
    R_s = elements[0]
    R_i_end = -1 if include_L else len(elements)
    R_i = elements[1:R_i_end]

    # Calculate Z_fit from fitted parameters
    omega = 2 * np.pi * frequencies
    Z_fit = np.full_like(frequencies, R_s, dtype=complex)

    for i, (r, t) in enumerate(zip(R_i, tau)):
        # K element: Z = R / (1 + jωτ)
        Z_fit += r / (1 + 1j * omega * t)

    # Add inductance contribution
    if include_L and L_value is not None:
        Z_fit += 1j * omega * L_value

    # Calculate residuals normalized by |Z| (Lin-KK standard)
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    res_real = (Z.real - Z_fit.real) / Z_mag_safe
    res_imag = (Z.imag - Z_fit.imag) / Z_mag_safe

    logger.info(f"Lin-KK native: M={M}, μ={mu:.4f}")
    logger.info(f"  Mean |res_real|: {np.mean(np.abs(res_real))*100:.2f}%")
    logger.info(f"  Mean |res_imag|: {np.mean(np.abs(res_imag))*100:.2f}%")

    return M, mu, Z_fit, res_real, res_imag, L_value


def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50
) -> Tuple[
    Optional[int],
    Optional[float],
    Optional[NDArray[np.complex128]],
    Optional[Tuple[NDArray[np.float64], NDArray[np.float64]]],
    Optional[plt.Figure]
]:
    """
    Perform Kramers-Kronig validation test on EIS data.

    Uses native Lin-KK implementation (Schönleber et al. 2014).

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ω]
    mu_threshold : float, optional
        Threshold for μ metric (default: 0.85)
    max_M : int, optional
        Maximum number of Voigt elements (default: 50)

    Returns
    -------
    M : int or None
        Number of RC elements in KK representation
    mu : float or None
        μ metric value from KK test
    Z_fit : ndarray of complex or None
        Fitted impedance [Ω]
    residuals : tuple or None
        Tuple (res_real, res_imag) of residuals (fractions, not %)
    fig : matplotlib.figure.Figure or None
        Visualization figure
    """
    logger.info("="*60)
    logger.info("Kramers-Kronig validation")
    logger.info("="*60)

    try:
        M, mu, Z_fit, res_real, res_imag, L_value = lin_kk_native(
            frequencies, Z,
            mu_threshold=mu_threshold,
            max_M=max_M,
            include_L=True,
            fit_type='real',
            weighting='proportional'
        )
        if L_value is not None:
            logger.info(f"Inductance L: {L_value:.3e} H")
    except Exception as e:
        logger.error(f"KK validation error: {e}", exc_info=True)
        return None, None, None, None, None

    if Z_fit is None:
        return None, None, None, None, None

    # Residuals are normalized by |Z| (fractions), convert to % for evaluation
    mean_res_real = np.mean(np.abs(res_real)) * 100
    mean_res_imag = np.mean(np.abs(res_imag)) * 100

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Fit comparison (Nyquist plot)
    ax1 = axes[0]
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=4)
    ax1.plot(Z_fit.real, -Z_fit.imag, '-', label='KK fit', linewidth=2)
    ax1.set_xlabel("Z' [Ω]")
    ax1.set_ylabel("-Z'' [Ω]")
    ax1.set_title(f"Kramers-Kronig fit in real domain (M={M})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='datalim')

    # Residuals plot (res_real/res_imag are fractions, convert to %)
    ax2 = axes[1]
    ax2.semilogx(frequencies, res_real * 100, 'o', label='Real', markersize=4)
    ax2.semilogx(frequencies, res_imag * 100, 's', label='Imaginary', markersize=4)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(y=5, color='r', linestyle=':', alpha=0.5)
    ax2.axhline(y=-5, color='r', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Residuals [%]")
    ax2.set_title(f"KK transformation residuals (μ={mu:.3f})")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Quality assessment
    if mean_res_real < 5 and mean_res_imag < 5:
        logger.info("Data quality is good (residuals < 5%)")
    else:
        logger.warning("Data may contain artifacts (residuals >= 5%)")

    return M, mu, Z_fit, (res_real, res_imag), fig
