"""
Kramers-Kronig validation for EIS data quality assessment.

Clean design: No logging in core functions, all diagnostics returned as data.

Provides two implementations:
1. lin_kk_native() - Native implementation using Voigt chain (no external dependencies)
2. kramers_kronig_validation() - High-level wrapper with visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass, field
from typing import Tuple, Optional, List
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
    warnings: List[str] = field(default_factory=list)

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


@dataclass
class LinKKResult:
    """Result of Lin-KK native fitting.

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
        Estimated noise in percent
    extend_decades : float
        Tau range extension in decades
    inductance : Optional[float]
        Fitted series inductance [H]
    elements : NDArray[np.float64]
        Fitted elements [R_s, R_1, ..., R_M]
    tau : NDArray[np.float64]
        Time constants [s]
    weighting : str
        Weighting scheme used
    """
    M: int
    mu: float
    Z_fit: NDArray[np.complex128]
    residuals_real: NDArray[np.float64]
    residuals_imag: NDArray[np.float64]
    pseudo_chisqr: float
    noise_estimate: float
    extend_decades: float
    inductance: Optional[float]
    elements: NDArray[np.float64]
    tau: NDArray[np.float64]
    weighting: str = 'modulus'

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
    search_range: Tuple[float, float] = (0.0, 1.0),
    n_evaluations: int = 11,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'modulus'
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
        Number of Voigt elements
    search_range : tuple
        Min and max extend_decades to search
    n_evaluations : int
        Number of grid points
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
    min_chi2 = min(r[1] for r in results)
    tolerance = 0.001 * min_chi2
    near_optimal = [r for r in results if r[1] <= min_chi2 + tolerance]
    best = min(near_optimal, key=lambda x: abs(x[0]))
    return best[0], best[1], best[2], best[3], best[4]


def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'modulus',
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0)
) -> LinKKResult:
    """
    Native Lin-KK implementation using Voigt chain fitting.

    Implements the linear Kramers-Kronig test from Schönleber et al. (2014)
    without external dependencies.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85)
    max_M : int, optional
        Maximum number of Voigt elements to try (default: 50)
    include_L : bool, optional
        Include series inductance L (default: True)
    fit_type : str, optional
        Fit type: 'real', 'imag', or 'complex'
    weighting : str, optional
        Weighting scheme: 'uniform', 'sqrt', 'proportional', 'modulus'
    auto_extend_decades : bool, optional
        Automatically optimize extend_decades (default: False)
    extend_decades_range : tuple, optional
        Search range for extend_decades optimization

    Returns
    -------
    LinKKResult
        Dataclass containing M, mu, Z_fit, residuals, pseudo_chisqr,
        noise_estimate, extend_decades, inductance, elements, tau.

    References
    ----------
    Schönleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    from ..fitting.voigt_chain import find_optimal_M_mu

    # Find optimal M using mu metric
    M, mu, tau, elements, L_value = find_optimal_M_mu(
        frequencies, Z,
        mu_threshold=mu_threshold,
        max_M=max_M,
        extend_decades=0.0,
        include_Rs=True,
        include_L=include_L,
        fit_type=fit_type,
        allow_negative=True,
        weighting=weighting
    )

    extend_decades = 0.0

    # Optionally optimize extend_decades
    if auto_extend_decades:
        extend_decades, chi2_opt, tau, elements, L_value = find_optimal_extend_decades(
            frequencies, Z, M,
            search_range=extend_decades_range,
            n_evaluations=11,
            include_L=include_L,
            fit_type=fit_type,
            weighting=weighting
        )

    # Reconstruct Z_fit from fitted parameters
    Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L)

    # Calculate residuals normalized by |Z|
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    res_real = (Z.real - Z_fit.real) / Z_mag_safe
    res_imag = (Z.imag - Z_fit.imag) / Z_mag_safe

    # Compute pseudo chi-squared (Boukamp 1995)
    chi2_ps = compute_pseudo_chisqr(Z, Z_fit)
    noise_est = estimate_noise_percent(chi2_ps, len(Z))

    return LinKKResult(
        M=M,
        mu=mu,
        Z_fit=Z_fit,
        residuals_real=res_real,
        residuals_imag=res_imag,
        pseudo_chisqr=chi2_ps,
        noise_estimate=noise_est,
        extend_decades=extend_decades,
        inductance=L_value,
        elements=elements,
        tau=tau,
        weighting=weighting
    )


def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0)
) -> Optional[KKResult]:
    """
    Perform Kramers-Kronig validation test on EIS data.

    Uses native Lin-KK implementation (Schönleber et al. 2014).

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85)
    max_M : int, optional
        Maximum number of Voigt elements (default: 50)
    auto_extend_decades : bool, optional
        Automatically optimize extend_decades (default: False)
    extend_decades_range : tuple of float, optional
        Search range for extend_decades optimization

    Returns
    -------
    KKResult or None
        Result object with all diagnostics, or None if validation fails.
    """
    warnings = []

    try:
        lkk = lin_kk_native(
            frequencies, Z,
            mu_threshold=mu_threshold,
            max_M=max_M,
            include_L=True,
            fit_type='real',
            weighting='modulus',
            auto_extend_decades=auto_extend_decades,
            extend_decades_range=extend_decades_range
        )
    except Exception as e:
        logger.error(f"KK validation error: {e}")
        return None

    if lkk.Z_fit is None:
        return None

    # Quality assessment
    if not lkk.is_valid:
        warnings.append("Data may contain artifacts (residuals >= 5%)")

    # Generate interpolated frequencies for smooth curve
    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = reconstruct_impedance(freq_plot, lkk.elements, lkk.tau, lkk.inductance, include_L=True)

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Fit comparison (Nyquist plot)
    ax1 = axes[0]
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=4)
    ax1.plot(Z_fit_plot.real, -Z_fit_plot.imag, '-', label='KK fit', linewidth=2)
    ax1.set_xlabel("Z' [Ohm]")
    ax1.set_ylabel("-Z'' [Ohm]")
    ax1.set_title(f"Kramers-Kronig fit in real domain (M={lkk.M})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='datalim')

    # Residuals plot
    ax2 = axes[1]
    ax2.semilogx(frequencies, lkk.residuals_real * 100, 'o', label='Real', markersize=4)
    ax2.semilogx(frequencies, lkk.residuals_imag * 100, 's', label='Imaginary', markersize=4)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(y=5, color='r', linestyle=':', alpha=0.5)
    ax2.axhline(y=-5, color='r', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Residuals [%]")
    ax2.set_title(f"KK residuals (mu={lkk.mu:.3f}, chi^2={lkk.pseudo_chisqr:.2e}, noise~{lkk.noise_estimate:.1f}%)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return KKResult(
        M=lkk.M,
        mu=lkk.mu,
        Z_fit=lkk.Z_fit,
        residuals_real=lkk.residuals_real,
        residuals_imag=lkk.residuals_imag,
        pseudo_chisqr=lkk.pseudo_chisqr,
        noise_estimate=lkk.noise_estimate,
        extend_decades=lkk.extend_decades,
        inductance=lkk.inductance,
        figure=fig,
        warnings=warnings
    )
