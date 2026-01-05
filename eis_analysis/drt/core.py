"""
Core DRT (Distribution of Relaxation Times) calculation.

Refactored into smaller, testable functions:
- DRTResult: dataclass for results
- _estimate_r_inf(): R_inf estimation
- _validate_frequencies(): input validation
- _build_drt_matrices(): matrix construction
- _select_lambda(): regularization parameter selection
- _solve_nnls(): NNLS solver with validation
- _create_visualization(): plotting
- calculate_drt(): main orchestrator
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass
from typing import Tuple, Optional, List, Dict
from numpy.typing import NDArray
from scipy.signal import find_peaks
from scipy.optimize import nnls

from .gcv import find_optimal_lambda_gcv, find_optimal_lambda_hybrid
from ..rinf_estimation import estimate_rinf_with_inductance
from .peaks import gmm_peak_detection, GMM_AVAILABLE
from ..utils.compat import np_trapz
from ..fitting.config import DRT_PEAK_HEIGHT_THRESHOLD

logger = logging.getLogger(__name__)

# Constants
MIN_FREQUENCY_RANGE = 10
DRT_TOLERANCE = 1e-10
GAMMA_MAX_REASONABLE = 1e10
GAMMA_MIN_REASONABLE = 1e-15


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class DRTResult:
    """
    Container for DRT analysis results.

    Replaces the 5-tuple return value with a structured object.
    """
    tau: Optional[NDArray[np.float64]] = None
    gamma: Optional[NDArray[np.float64]] = None
    gamma_original: Optional[NDArray[np.float64]] = None  # Before normalization
    figure: Optional[plt.Figure] = None
    peaks: Optional[List[Dict]] = None
    figure_rinf: Optional[plt.Figure] = None

    # Additional diagnostics
    R_inf: Optional[float] = None
    R_pol: Optional[float] = None
    lambda_used: Optional[float] = None
    reconstruction_error: Optional[float] = None

    def as_tuple(self) -> Tuple:
        """Backward compatibility: return as original tuple format."""
        return (self.tau, self.gamma, self.figure, self.peaks, self.figure_rinf)

    @property
    def success(self) -> bool:
        """Check if DRT calculation was successful."""
        return self.tau is not None and self.gamma is not None


@dataclass
class DRTMatrices:
    """Container for DRT matrices."""
    A: NDArray[np.float64]      # System matrix [2N x n_tau]
    A_re: NDArray[np.float64]   # Real part of A
    A_im: NDArray[np.float64]   # Imaginary part of A
    b: NDArray[np.float64]      # Right-hand side
    L: NDArray[np.float64]      # Regularization matrix
    tau: NDArray[np.float64]    # Time constants
    d_ln_tau: float             # Log-spacing of tau


# =============================================================================
# Helper Functions
# =============================================================================

def _estimate_peak_resistance(tau: NDArray, gamma: NDArray,
                               peak_indices: NDArray,
                               tolerance: float = 0.1) -> List[float]:
    """
    Estimate resistance for each peak by integrating gamma in peak region.

    Parameters
    ----------
    tau : ndarray
        Time constants
    gamma : ndarray
        DRT spectrum
    peak_indices : ndarray
        Indices of peaks
    tolerance : float
        Threshold for peak boundaries (default: 0.1 = 10% of height)

    Returns
    -------
    list of float
        Resistance estimates for each peak [Ω]
    """
    if len(peak_indices) == 0:
        return []

    resistances = []
    ln_tau = np.log(tau)

    for peak_idx in peak_indices:
        peak_height = gamma[peak_idx]
        threshold = peak_height * tolerance

        # Find left boundary
        left = peak_idx
        while left > 0 and gamma[left] > threshold:
            left -= 1

        # Find right boundary
        right = peak_idx
        while right < len(gamma) - 1 and gamma[right] > threshold:
            right += 1

        # Integrate over peak
        if right > left:
            R_peak = np_trapz(gamma[left:right+1], ln_tau[left:right+1])
        else:
            R_peak = 0.0

        resistances.append(R_peak)

    return resistances


# =============================================================================
# R_inf Estimation
# =============================================================================

def _estimate_r_inf(frequencies: NDArray, Z: NDArray,
                    use_rl_fit: bool = False,
                    use_voigt_fit: bool = False,
                    r_inf_preset: Optional[float] = None
                    ) -> Tuple[float, Optional[plt.Figure], Optional[Dict]]:
    """
    Estimate high-frequency resistance R_inf.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance
    use_rl_fit : bool
        Use R+L fit for inductive data
    use_voigt_fit : bool
        Use Voigt fit for capacitive data
    r_inf_preset : float, optional
        Use this preset value instead of calculating

    Returns
    -------
    R_inf : float
        High-frequency resistance [Ω]
    fig_rl : Figure or None
        Diagnostic figure if fit was performed
    diag : dict or None
        Diagnostic info from fit
    """
    n_avg = min(5, max(1, len(frequencies) // 10))
    fig_rl = None
    diag = None

    if r_inf_preset is not None:
        logger.info("="*60)
        logger.info("R_inf: Using preset value from previous --ri-fit")
        logger.info("="*60)
        R_inf = r_inf_preset
        logger.info(f"R_inf = {R_inf:.3f} Ω (from preset)")
        logger.info("="*60)
        return R_inf, None, None

    logger.info("="*60)
    logger.info("R_inf estimation (high-frequency resistance)")
    logger.info("="*60)

    # Always calculate median for comparison
    high_freq_indices = np.argsort(frequencies)[-n_avg:]
    R_inf_median = np.median(Z.real[high_freq_indices])

    if use_rl_fit or use_voigt_fit:
        method_name = "Voigt fit" if use_voigt_fit else "Auto-detection (--ri-fit)"
        logger.info(f"Method: {method_name}")

        try:
            R_inf, L_fit, circuit_rl, diag, fig_rl = estimate_rinf_with_inductance(
                frequencies, Z, use_voigt_fit=use_voigt_fit, verbose=False, plot=True
            )

            behavior = diag.get('behavior', 'unknown')
            method = diag.get('method', 'unknown')

            if behavior == 'capacitive':
                logger.info("  Detected: Capacitive behavior (Im(Z) < 0)")
                if 'voigt_fit' in method:
                    logger.info("  Method: Voigt fit (3 params: R_inf, R_ct, C)")
                else:
                    logger.info("  Method: Re(Z) extrapolation to infinity")
            elif behavior == 'inductive':
                logger.info("  Detected: Inductive behavior (Im(Z) > 0)")
                logger.info("  Method: R-L fit (resistance + inductance)")
            else:
                logger.info(f"  Detected: {behavior}")

            logger.info(f"R_inf = {R_inf:.3f} Ω ({diag['n_points_used']} HF points)")

            if 'r_squared' in diag and diag['r_squared'] > 0:
                logger.info(f"  Quality: R² = {diag['r_squared']:.4f}")
            if 'L_nH' in diag and diag['L_nH'] > 0:
                logger.info(f"  Inductance: L = {diag['L_nH']:.2f} nH")

            if 'voigt' in method and 'R_ct' in diag and 'C_nF' in diag:
                logger.info(f"  Voigt params: R_ct = {diag['R_ct']:.2f} Ω, C = {diag['C_nF']:.2f} nF")
                if 'f_characteristic' in diag:
                    logger.info(f"  Characteristic freq: f_char = {diag['f_characteristic']/1e6:.3f} MHz")

            # Comparison with median
            diff_abs = R_inf - R_inf_median
            diff_pct = (diff_abs / R_inf_median * 100) if R_inf_median != 0 else 0
            logger.info(f"  Comparison: median = {R_inf_median:.3f} Ω (diff: {diff_abs:+.3f} Ω, {diff_pct:+.1f}%)")

            # Warnings
            if diag.get('warnings'):
                for warning in diag['warnings']:
                    logger.warning(f"  {warning}")

            if diag.get('L_nH', 0) > 500:
                logger.warning("="*60)
                logger.warning(f"WARNING: High inductance L = {diag['L_nH']:.1f} nH detected")
                logger.warning("  Possible causes: ground loops, long cables, bad grounding")
                logger.warning("="*60)

            logger.info("="*60)
            return R_inf, fig_rl, diag

        except Exception as e:
            logger.error(f"Auto-estimation failed: {e}")
            logger.warning("Fallback to median method")

    # Median method (default)
    logger.info("Method: Median of HF points (classic)")
    R_inf = R_inf_median
    logger.info(f"R_inf = {R_inf:.3f} Ω (median of {n_avg} highest frequencies)")
    logger.info("="*60)

    return R_inf, None, None


# =============================================================================
# Input Validation
# =============================================================================

def _validate_frequencies(frequencies: NDArray) -> None:
    """
    Validate frequency array for DRT analysis.

    Raises ValueError if validation fails.
    """
    if np.any(frequencies <= 0):
        raise ValueError("All frequencies must be positive for DRT analysis")

    f_max = frequencies.max()
    f_min = frequencies.min()

    if f_max <= 0 or f_min <= 0:
        raise ValueError(f"Invalid frequency range: min={f_min}, max={f_max}")

    if f_max <= f_min:
        raise ValueError(f"Max frequency ({f_max:.2e} Hz) must be > min ({f_min:.2e} Hz)")

    freq_range = f_max / f_min
    if freq_range < MIN_FREQUENCY_RANGE:
        logger.warning(f"Small frequency range for DRT: {freq_range:.1f}× (recommended >{MIN_FREQUENCY_RANGE}×)")
        logger.warning("DRT spectrum will have poor resolution")


# =============================================================================
# Matrix Construction
# =============================================================================

def _build_drt_matrices(frequencies: NDArray, Z: NDArray,
                        R_inf: float, n_tau: int = 100) -> DRTMatrices:
    """
    Build DRT system matrices A, b, and regularization matrix L.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance
    R_inf : float
        High-frequency resistance
    n_tau : int
        Number of tau points

    Returns
    -------
    DRTMatrices
        Container with A, A_re, A_im, b, L, tau, d_ln_tau
    """
    f_max = frequencies.max()
    f_min = frequencies.min()

    # Log-spaced tau grid
    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    omega = 2 * np.pi * frequencies

    # Check tau uniformity in log space
    d_ln_tau_array = np.diff(np.log(tau))
    if not np.allclose(d_ln_tau_array, d_ln_tau_array[0], rtol=DRT_TOLERANCE):
        logger.warning("Tau is not uniform in log space")
        logger.warning(f"d(ln tau) variation: {d_ln_tau_array.min():.3e} - {d_ln_tau_array.max():.3e}")

    d_ln_tau = np.mean(d_ln_tau_array)
    logger.debug(f"d(ln tau) = {d_ln_tau:.6f}")

    # Vectorized matrix construction
    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom

    A = np.vstack([A_re, A_im])

    # Condition number diagnostic
    cond_A = np.linalg.cond(A)
    logger.debug(f"Matrix A condition number: {cond_A:.2e}")
    if cond_A > 1e15:
        logger.warning("="*60)
        logger.warning(f"WARNING: Matrix A is ill-conditioned ({cond_A:.2e})")
        logger.warning("  Results may be numerically unstable.")
        logger.warning("  Recommendations: increase --lambda, reduce --n-tau, check data quality")
        logger.warning("="*60)
    elif cond_A > 1e12:
        logger.info(f"Matrix A has high condition number ({cond_A:.2e})")
        logger.info("  Regularization should prevent instability.")

    # Right-hand side
    b = np.concatenate([Z.real - R_inf, Z.imag])

    # Regularization matrix (2nd derivative)
    L = np.zeros((n_tau - 2, n_tau))
    for i in range(n_tau - 2):
        L[i, i] = 1
        L[i, i + 1] = -2
        L[i, i + 2] = 1

    return DRTMatrices(A=A, A_re=A_re, A_im=A_im, b=b, L=L, tau=tau, d_ln_tau=d_ln_tau)


# =============================================================================
# Lambda Selection
# =============================================================================

def _select_lambda(A: NDArray, b: NDArray, L: NDArray,
                   lambda_reg: Optional[float] = None,
                   auto_lambda: bool = False) -> float:
    """
    Select regularization parameter lambda.

    Parameters
    ----------
    A, b, L : ndarray
        DRT matrices
    lambda_reg : float, optional
        User-specified lambda
    auto_lambda : bool
        Use automatic selection (hybrid GCV + L-curve)

    Returns
    -------
    float
        Selected lambda value
    """
    if auto_lambda:
        logger.info("Mode: Automatic λ selection (hybrid GCV + L-curve)")
        try:
            lambda_reg, gcv_score, diag = find_optimal_lambda_hybrid(
                A, b, L,
                lambda_range=(1e-5, 1.0),
                n_search=20,
                lcurve_decades=1.5
            )
            if diag['method_used'] == 'lcurve_correction':
                logger.info(f"  L-curve correction: λ_gcv={diag['lambda_gcv']:.2e} → λ={lambda_reg:.2e}")
            return lambda_reg

        except Exception as e:
            logger.error(f"Hybrid optimization failed: {e}")
            logger.warning("Fallback to pure GCV...")
            try:
                lambda_reg, _ = find_optimal_lambda_gcv(A, b, L,
                                                        lambda_range=(1e-5, 1.0),
                                                        n_search=20)
                return lambda_reg
            except Exception as e2:
                logger.error(f"GCV also failed: {e2}")
                logger.warning("Using default λ = 0.1")
                return 0.1
    else:
        if lambda_reg is None:
            lambda_reg = 0.1
            logger.info(f"Using default λ = {lambda_reg}")
        else:
            logger.info(f"Using user-specified λ = {lambda_reg}")
        return lambda_reg


# =============================================================================
# NNLS Solver
# =============================================================================

def _solve_nnls(A: NDArray, b: NDArray, L: NDArray,
                lambda_reg: float, n_tau: int,
                Z: NDArray) -> Tuple[Optional[NDArray], bool]:
    """
    Solve regularized NNLS problem.

    Parameters
    ----------
    A, b, L : ndarray
        DRT matrices
    lambda_reg : float
        Regularization parameter
    n_tau : int
        Number of tau points
    Z : ndarray
        Original impedance (for inductive check)

    Returns
    -------
    gamma : ndarray or None
        DRT spectrum, None if failed
    success : bool
        Whether solution is valid
    """
    # Check for inductive data
    n_inductive = np.sum(Z.imag > 0)
    if n_inductive > 0:
        max_inductive = np.max(Z.imag)
        inductive_fraction = n_inductive / len(Z)

        if inductive_fraction > 0.1 or max_inductive > 50:
            logger.warning("="*60)
            logger.warning(f"WARNING: {n_inductive} points with inductive component ({inductive_fraction*100:.1f}%)")
            logger.warning(f"  Max Z'' = {max_inductive:.2f} Ω (positive imaginary)")
            logger.warning("  NNLS enforces γ(τ) ≥ 0, may be unsuitable for inductive processes.")
            logger.warning("="*60)
        else:
            logger.debug(f"Detected {n_inductive} inductive points (max Z''={max_inductive:.2f} Ω)")

    # Build regularized system
    A_reg = np.vstack([A, np.sqrt(lambda_reg) * L])
    b_reg = np.concatenate([b, np.zeros(n_tau - 2)])

    # Solve NNLS
    try:
        gamma, residual = nnls(A_reg, b_reg)
    except Exception as e:
        logger.error(f"NNLS solver error: {e}")
        logger.error("Matrix may be ill-conditioned. Try adjusting --lambda or --n-tau")
        return None, False

    # Validate solution
    if np.any(~np.isfinite(gamma)):
        logger.error("NNLS returned NaN or Inf values")
        return None, False

    if np.max(gamma) > GAMMA_MAX_REASONABLE:
        logger.warning(f"DRT contains very large values (max: {np.max(gamma):.2e} Ω)")
        logger.warning("Try increasing regularization --lambda")

    gamma_nonzero = gamma[gamma > 0]
    if len(gamma_nonzero) > 0 and np.min(gamma_nonzero) < GAMMA_MIN_REASONABLE:
        logger.debug(f"DRT contains very small nonzero values (min: {np.min(gamma_nonzero):.2e} Ω)")

    if np.sum(gamma) < 1e-10:
        logger.warning("DRT is nearly zero - data may have no relaxation structure")

    return gamma, True


# =============================================================================
# Visualization
# =============================================================================

def _create_visualization(tau: NDArray, gamma: NDArray,
                          gamma_original: Optional[NDArray],
                          Z: NDArray, Z_reconstructed: NDArray,
                          R_inf: float, lambda_reg: float,
                          normalize_rpol: bool,
                          peak_method: str,
                          peaks_result: Optional[List[Dict]],
                          bic_scores: Optional[List[float]]) -> plt.Figure:
    """
    Create DRT visualization figure.

    Parameters
    ----------
    tau : ndarray
        Time constants
    gamma : ndarray
        DRT spectrum (possibly normalized)
    gamma_original : ndarray or None
        Original gamma before normalization
    Z : ndarray
        Original impedance
    Z_reconstructed : ndarray
        Reconstructed impedance from DRT
    R_inf : float
        High-frequency resistance
    lambda_reg : float
        Used regularization parameter
    normalize_rpol : bool
        Whether gamma was normalized
    peak_method : str
        Peak detection method ('scipy' or 'gmm')
    peaks_result : list or None
        Detected peaks (GMM format)
    bic_scores : list or None
        BIC scores for GMM

    Returns
    -------
    Figure
        Matplotlib figure
    """
    use_gmm = (peak_method == 'gmm' and GMM_AVAILABLE and
               peaks_result is not None and len(peaks_result) > 0)

    if use_gmm:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        ax1, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        ax1, ax2 = axes[0], axes[1]

    # === DRT Spectrum ===
    ax1.semilogx(tau, gamma, 'b-', linewidth=2, label='DRT γ(τ)')
    ax1.fill_between(tau, 0, gamma, alpha=0.3)
    ax1.set_xlabel("τ [s]")

    if normalize_rpol:
        ax1.set_ylabel("γ(τ) / R_pol [-]")
        ax1.set_title(f"DRT normalized (λ = {lambda_reg})")
    else:
        ax1.set_ylabel("γ(τ) [Ω]")
        ax1.set_title(f"DRT (λ = {lambda_reg})")

    ax1.grid(True, alpha=0.3, which='both')

    # Mark peaks for scipy method
    if not use_gmm:
        peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD)
        if len(peaks_idx) > 0:
            ax1.plot(tau[peaks_idx], gamma[peaks_idx], 'ro', markersize=8,
                    label=f'{len(peaks_idx)} peaks', zorder=5)
            ax1.legend()

    # === Nyquist Comparison ===
    ax2.plot(Z.real, -Z.imag, 'o', label='Data', markersize=5)
    ax2.plot(Z_reconstructed.real, -Z_reconstructed.imag, '-',
             label='DRT reconstruction', linewidth=2)
    ax2.set_xlabel("Z' [Ω]")
    ax2.set_ylabel("-Z'' [Ω]")
    ax2.set_title("DRT fit verification")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # === GMM Visualization ===
    if use_gmm:
        log_tau = np.log10(tau)
        ax3.semilogx(tau, gamma, 'b-', linewidth=2, label='DRT γ(τ)', alpha=0.7)

        colors = plt.cm.tab10(np.linspace(0, 1, len(peaks_result)))
        for i, peak in enumerate(peaks_result):
            mu = np.log10(peak['tau_center'])
            sigma = peak['log_tau_std']

            gaussian_shape = np.exp(-0.5*((log_tau - mu)/sigma)**2)
            integral_log10 = sigma * np.sqrt(2 * np.pi)
            integral_ln = integral_log10 * np.log(10)
            height = peak['R_estimate'] / integral_ln
            gaussian_gamma = height * gaussian_shape

            ax3.fill_between(tau, 0, gaussian_gamma, alpha=0.4, color=colors[i],
                            label=f"Peak {i+1}: τ={peak['tau_center']:.2e}s")
            ax3.axvline(peak['tau_bounds'][0], color=colors[i], linestyle=':', alpha=0.8)
            ax3.axvline(peak['tau_bounds'][1], color=colors[i], linestyle=':', alpha=0.8)
            ax3.axvline(peak['tau_center'], color=colors[i], linestyle='--', alpha=0.8)

        ax3.set_xlabel("τ [s]")
        ax3.set_ylabel("γ(τ) [Ω]")
        ax3.set_title(f"GMM deconvolution ({len(peaks_result)} peaks)")
        ax3.legend(fontsize=8, loc='best')
        ax3.grid(True, alpha=0.3, which='both')

        # BIC plot
        if bic_scores:
            n_components_range = (1, 6)
            n_range = range(n_components_range[0], n_components_range[0] + len(bic_scores))
            valid_bic = [bic for bic in bic_scores if bic != np.inf]

            if valid_bic:
                ax4.plot(list(n_range), bic_scores, 'o-', linewidth=2, markersize=8)
                best_n = len(peaks_result)
                ax4.axvline(best_n, color='r', linestyle='--', alpha=0.7,
                           label=f'Optimal n={best_n}')
                ax4.set_xlabel("Number of components")
                ax4.set_ylabel("BIC")
                ax4.set_title("Model selection (lower BIC = better)")
                ax4.legend()
                ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# =============================================================================
# Peak Detection
# =============================================================================

def _detect_peaks(tau: NDArray, gamma: NDArray,
                  peak_method: str,
                  gmm_bic_threshold: float = 10.0) -> Tuple[Optional[List[Dict]], Optional[List[float]]]:
    """
    Detect peaks in DRT spectrum.

    Returns (peaks_result, bic_scores) for GMM or (None, None) for scipy.
    """
    logger.info("="*60)
    logger.info("Peak detection in DRT spectrum")
    logger.info("="*60)

    use_gmm = (peak_method == 'gmm' and GMM_AVAILABLE)

    if peak_method == 'gmm':
        if GMM_AVAILABLE:
            logger.info("Method: GMM (Gaussian Mixture Models)")
            logger.info("  → Robust detection with explicit peak boundaries")
        else:
            logger.warning("GMM unavailable (scikit-learn missing)")
            logger.warning("  → Fallback to scipy.signal.find_peaks")
            use_gmm = False

    if not use_gmm:
        logger.info("Method: scipy.signal.find_peaks")
        logger.info("  → Fast detection for well-separated peaks")

    logger.info("="*60)

    if use_gmm:
        peaks_result, gmm_model, bic_scores = gmm_peak_detection(
            tau, gamma,
            bic_threshold=gmm_bic_threshold
        )

        if len(peaks_result) == 0 or gmm_model is None:
            logger.warning("GMM detection failed, using scipy fallback")
            peaks_result = None
            bic_scores = None

        return peaks_result, bic_scores

    # Scipy method
    peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD)
    peak_resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)

    logger.info(f"Found {len(peaks_idx)} peaks:")
    for i, idx in enumerate(peaks_idx):
        R_peak = peak_resistances[i] if i < len(peak_resistances) else 0.0
        logger.info(f"  Peak {i+1}: τ = {tau[idx]:.2e} s (f = {1/tau[idx]:.2e} Hz), R ~ {R_peak:.2f} Ω")

    return None, None


# =============================================================================
# Main Orchestrator
# =============================================================================

def calculate_drt(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    n_tau: int = 100,
    lambda_reg: Optional[float] = None,
    auto_lambda: bool = False,
    normalize_rpol: bool = False,
    use_rl_fit: bool = False,
    use_voigt_fit: bool = False,
    peak_method: str = 'scipy',
    r_inf_preset: Optional[float] = None,
    gmm_bic_threshold: float = 10.0
) -> Tuple[
    Optional[NDArray[np.float64]],
    Optional[NDArray[np.float64]],
    Optional[plt.Figure],
    Optional[List[Dict]],
    Optional[plt.Figure]
]:
    """
    Calculate DRT (Distribution of Relaxation Times) using Tikhonov regularization.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    n_tau : int
        Number of tau points (default: 100)
    lambda_reg : float, optional
        Regularization parameter
    auto_lambda : bool
        Auto-select lambda using hybrid GCV + L-curve
    normalize_rpol : bool
        Normalize gamma by R_pol
    use_rl_fit : bool
        Use R+L fit for R_inf estimation
    use_voigt_fit : bool
        Use Voigt fit for R_inf estimation
    peak_method : str
        Peak detection method ('scipy' or 'gmm')
    r_inf_preset : float, optional
        Preset R_inf value

    Returns
    -------
    tuple
        (tau, gamma, fig, peaks_gmm, fig_rl)
        Returns (None, None, None, None, None) on error
    """
    # === Step 1: R_inf Estimation ===
    R_inf, fig_rl, diag_rl = _estimate_r_inf(
        frequencies, Z,
        use_rl_fit=use_rl_fit,
        use_voigt_fit=use_voigt_fit,
        r_inf_preset=r_inf_preset
    )

    # === Step 2: Input Validation ===
    logger.info("="*60)
    logger.info("DRT Analysis")
    logger.info("="*60)
    logger.info(f"Using R_inf = {R_inf:.3f} Ω")

    try:
        _validate_frequencies(frequencies)
    except ValueError as e:
        logger.error(str(e))
        return None, None, None, None, None

    # === Step 3: Build Matrices ===
    matrices = _build_drt_matrices(frequencies, Z, R_inf, n_tau)

    # === Step 4: Select Lambda ===
    lambda_reg = _select_lambda(matrices.A, matrices.b, matrices.L,
                                 lambda_reg, auto_lambda)
    logger.info(f"Regularization parameter λ = {lambda_reg}")

    # === Step 5: Solve NNLS ===
    gamma, success = _solve_nnls(matrices.A, matrices.b, matrices.L,
                                  lambda_reg, n_tau, Z)

    if not success:
        return None, None, None, None, None

    # === Step 6: R_pol Calculation & Normalization ===
    n_avg = min(5, max(1, len(frequencies) // 10))
    low_freq_indices = np.argsort(frequencies)[:n_avg]
    R_dc = np.median(Z.real[low_freq_indices])
    R_pol = R_dc - R_inf
    R_pol_from_gamma = np_trapz(gamma, np.log(matrices.tau))

    logger.info(f"R_pol (from data) = {R_pol:.2f} Ω")
    logger.info(f"R_pol (from DRT integral) = {R_pol_from_gamma:.2f} Ω")

    gamma_original = None
    if normalize_rpol:
        if R_pol_from_gamma > 1e-10:
            gamma_original = gamma.copy()
            gamma = gamma / R_pol_from_gamma
            logger.info("γ(τ) normalized by R_pol: ∫γ_norm d(ln τ) = 1")
        else:
            logger.warning("R_pol ≈ 0, normalization not possible")
            normalize_rpol = False

    # === Step 7: Peak Detection ===
    peaks_result, bic_scores = _detect_peaks(matrices.tau, gamma, peak_method, gmm_bic_threshold)

    # === Step 8: Reconstruction & Error ===
    gamma_for_recon = gamma_original if normalize_rpol else gamma
    Z_reconstructed = R_inf + (matrices.A_re + 1j * matrices.A_im) @ gamma_for_recon
    rel_error = np.mean(np.abs(Z - Z_reconstructed) / np.abs(Z)) * 100

    logger.info(f"Mean relative reconstruction error: {rel_error:.1f}%")

    if rel_error > 10.0:
        logger.warning("="*60)
        logger.warning("HIGH DRT RECONSTRUCTION ERROR!")
        logger.warning(f"Relative error {rel_error:.1f}% > 10% suggests DRT model")
        logger.warning("may not be suitable for this data.")
        logger.warning("Consider circuit fitting instead.")
        logger.warning("="*60)
    elif rel_error > 5.0:
        logger.warning(f"Elevated reconstruction error ({rel_error:.1f}%). Consider higher λ.")

    # === Step 9: Visualization ===
    fig = _create_visualization(
        matrices.tau, gamma, gamma_original,
        Z, Z_reconstructed,
        R_inf, lambda_reg,
        normalize_rpol, peak_method,
        peaks_result, bic_scores
    )

    return matrices.tau, gamma, fig, peaks_result, fig_rl
