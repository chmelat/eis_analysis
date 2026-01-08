"""
Kramers-Kronig validation for EIS data quality assessment.

Provides two implementations:
1. lin_kk_native() - Native implementation using Voigt chain (no external dependencies)
2. kramers_kronig_validation() - High-level wrapper with visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass
from typing import Tuple, Optional
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


@dataclass
class KKResult:
    """Result of Kramers-Kronig validation.

    Attributes
    ----------
    M : int
        Number of Voigt elements used
    mu : float
        Mu metric value (1.0 = no overfit, <0.85 = overfit)
    Z_fit : NDArray[np.complex128]
        Fitted impedance
    residuals_real : NDArray[np.float64]
        Real part residuals (normalized by |Z|)
    residuals_imag : NDArray[np.float64]
        Imaginary part residuals (normalized by |Z|)
    pseudo_chisqr : float
        Pseudo chi-squared (Boukamp 1995)
    noise_estimate : float
        Estimated noise in percent (Yrjana & Bobacka 2024)
    extend_decades : float
        Tau range extension in decades (0.0 = no extension)
    inductance : Optional[float]
        Fitted series inductance [H]
    figure : Optional[plt.Figure]
        Visualization figure
    """
    M: int
    mu: float
    Z_fit: NDArray[np.complex128]
    residuals_real: NDArray[np.float64]
    residuals_imag: NDArray[np.float64]
    pseudo_chisqr: float
    noise_estimate: float
    extend_decades: float = 0.0
    inductance: Optional[float] = None
    figure: Optional[plt.Figure] = None

    @property
    def mean_residual_real(self) -> float:
        """Mean absolute real residual in percent."""
        return float(np.mean(np.abs(self.residuals_real)) * 100)

    @property
    def mean_residual_imag(self) -> float:
        """Mean absolute imaginary residual in percent."""
        return float(np.mean(np.abs(self.residuals_imag)) * 100)

    @property
    def is_valid(self) -> bool:
        """Check if data passes KK validation (residuals < 5%)."""
        return self.mean_residual_real < 5 and self.mean_residual_imag < 5


def compute_pseudo_chisqr(
    Z_exp: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128]
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

    Based on Yrjana & Bobacka (2024).

    Parameters
    ----------
    chi2_ps : float
        Pseudo chi-squared value
    n_points : int
        Number of data points

    Returns
    -------
    float
        Estimated noise in percent

    References
    ----------
    Yrjana, V. and Bobacka, J. "Implementing Kramers-Kronig validity testing
    using pyimpspec." Electrochim. Acta 504, 144951 (2024)
    """
    return float(np.sqrt(chi2_ps * 5000 / n_points))


def reconstruct_impedance(
    frequencies: NDArray[np.float64],
    elements: NDArray[np.float64],
    tau: NDArray[np.float64],
    L_value: Optional[float],
    include_L: bool = True
) -> NDArray[np.complex128]:
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

    Returns
    -------
    Z_fit : array
        Reconstructed complex impedance
    """
    omega = 2 * np.pi * frequencies
    R_s = elements[0]
    R_i_end = -1 if include_L else len(elements)
    R_i = elements[1:R_i_end]

    Z_fit = np.full_like(frequencies, R_s, dtype=complex)
    for r, t in zip(R_i, tau):
        Z_fit += r / (1 + 1j * omega * t)

    if include_L and L_value is not None:
        Z_fit += 1j * omega * L_value

    return Z_fit


def find_optimal_extend_decades(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    M: int,
    search_range: Tuple[float, float] = (-1.0, 1.0),
    n_evaluations: int = 11,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'proportional'
) -> Tuple[float, float, NDArray[np.float64], NDArray[np.float64], Optional[float]]:
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
        Number of Voigt elements (from mu optimization)
    search_range : tuple
        Min and max extend_decades to search (default: -1.0 to 1.0)
    n_evaluations : int
        Number of grid points (default: 11)
    include_L : bool
        Include series inductance
    fit_type : str
        Fit type ('real', 'imag', 'complex')
    weighting : str
        Weighting scheme

    Returns
    -------
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
    """
    from ..fitting.voigt_chain import generate_tau_grid_fixed_M, estimate_R_linear

    candidates = np.linspace(search_range[0], search_range[1], n_evaluations)
    results = []

    for ext_dec in candidates:
        tau = generate_tau_grid_fixed_M(frequencies, M, extend_decades=ext_dec)
        elements, residual, L_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=True, include_L=include_L,
            fit_type=fit_type, allow_negative=True,
            weighting=weighting
        )

        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L)
        chi2 = compute_pseudo_chisqr(Z, Z_fit)
        results.append((ext_dec, chi2, tau, elements, L_value))

    # Find minimum chi^2
    best = min(results, key=lambda x: x[1])
    return best[0], best[1], best[2], best[3], best[4]


def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'proportional',
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (-1.0, 1.0)
) -> Tuple[int, float, NDArray[np.complex128], NDArray[np.float64], NDArray[np.float64], Optional[float], float, float]:
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
    auto_extend_decades : bool, optional
        If True, automatically optimize extend_decades to minimize chi^2
        (default: False)
    extend_decades_range : tuple, optional
        Search range for extend_decades optimization (default: (-1.0, 1.0))

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
    chi2_ps : float
        Pseudo chi-squared (Boukamp 1995)
    extend_decades : float
        Used extend_decades value (0.0 if not auto-optimized)

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

    # Find optimal M using μ metric (with extend_decades=0.0)
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

    # Track extend_decades value
    extend_decades = 0.0

    # Optionally optimize extend_decades
    if auto_extend_decades:
        logger.info(f"  Optimizing extend_decades in range {extend_decades_range}...")
        extend_decades, chi2_opt, tau, elements, L_value = find_optimal_extend_decades(
            frequencies, Z, M,
            search_range=extend_decades_range,
            n_evaluations=11,
            include_L=include_L,
            fit_type=fit_type,
            weighting=weighting
        )
        logger.info(f"  Optimal extend_decades: {extend_decades:.3f}")

    # Reconstruct Z_fit from fitted parameters
    Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L)

    # Calculate residuals normalized by |Z| (Lin-KK standard)
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    res_real = (Z.real - Z_fit.real) / Z_mag_safe
    res_imag = (Z.imag - Z_fit.imag) / Z_mag_safe

    # Compute pseudo chi-squared (Boukamp 1995)
    chi2_ps = compute_pseudo_chisqr(Z, Z_fit)
    noise_est = estimate_noise_percent(chi2_ps, len(Z))

    logger.info(f"Lin-KK native: M={M}, μ={mu:.4f}")
    logger.info(f"  Weighting: {weighting} (w=1/|Z|)")
    logger.info(f"  extend_decades: {extend_decades:.3f}")
    logger.info(f"  Mean |res_real|: {np.mean(np.abs(res_real))*100:.2f}%")
    logger.info(f"  Mean |res_imag|: {np.mean(np.abs(res_imag))*100:.2f}%")
    logger.info(f"  Pseudo chi^2: {chi2_ps:.2e}")
    logger.info(f"  Estimated noise: {noise_est:.2f}%")

    return M, mu, Z_fit, res_real, res_imag, L_value, chi2_ps, extend_decades


def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (-1.0, 1.0)
) -> Optional[KKResult]:
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
    auto_extend_decades : bool, optional
        If True, automatically optimize extend_decades to minimize chi^2
        (default: False)
    extend_decades_range : tuple of float, optional
        Search range for extend_decades optimization as (min, max).
        Only used when auto_extend_decades=True. Default: (-1.0, 1.0)

    Returns
    -------
    KKResult or None
        Result object containing:
        - M: Number of Voigt elements
        - mu: μ metric value
        - Z_fit: Fitted impedance
        - residuals_real, residuals_imag: Residuals (fractions)
        - pseudo_chisqr: Pseudo chi-squared (Boukamp 1995)
        - noise_estimate: Estimated noise in percent
        - extend_decades: Used tau range extension
        - inductance: Fitted inductance [H]
        - figure: Visualization figure

        Returns None if validation fails.
    """
    logger.info("="*60)
    logger.info("Kramers-Kronig validation")
    logger.info("="*60)

    try:
        # Weighting choice: 'proportional' (1/|Z|) vs 'modulus' (1/|Z|^2)
        #
        # We use 'proportional' (Lin-KK standard, Schönleber 2014) rather than
        # 'modulus' (Boukamp 1995 pseudo chi-squared weighting) because:
        # - Provides balanced relative weighting across entire spectrum
        # - Low-frequency region (high |Z|) often contains key electrochemical
        #   information (charge transfer, diffusion) that should not be de-emphasized
        # - Common artifacts (drift, non-stationarity) appear at low frequencies
        #   and would be masked with 1/|Z|^2 weighting
        #
        # Note: Pseudo chi-squared is still computed with 1/|Z|^2 per Boukamp (1995).
        # The difference only affects fitting, not the final chi^2 metric.
        M, mu, Z_fit, res_real, res_imag, L_value, chi2_ps, ext_dec = lin_kk_native(
            frequencies, Z,
            mu_threshold=mu_threshold,
            max_M=max_M,
            include_L=True,
            fit_type='real',
            weighting='proportional',
            auto_extend_decades=auto_extend_decades,
            extend_decades_range=extend_decades_range
        )
        if L_value is not None:
            logger.info(f"  Inductance L: {L_value:.3e} H")
    except Exception as e:
        logger.error(f"KK validation error: {e}", exc_info=True)
        return None

    if Z_fit is None:
        return None

    # Compute noise estimate from pseudo chi-squared
    noise_est = estimate_noise_percent(chi2_ps, len(Z))

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
    ax2.set_title(f"KK residuals (μ={mu:.3f}, χ²={chi2_ps:.2e}, noise~{noise_est:.1f}%)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Quality assessment
    if mean_res_real < 5 and mean_res_imag < 5:
        logger.info("Data quality is good (residuals < 5%)")
    else:
        logger.warning("Data may contain artifacts (residuals >= 5%)")

    return KKResult(
        M=M,
        mu=mu,
        Z_fit=Z_fit,
        residuals_real=res_real,
        residuals_imag=res_imag,
        pseudo_chisqr=chi2_ps,
        noise_estimate=noise_est,
        extend_decades=ext_dec,
        inductance=L_value,
        figure=fig
    )
