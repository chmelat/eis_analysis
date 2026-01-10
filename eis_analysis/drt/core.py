"""
Core DRT (Distribution of Relaxation Times) calculation.

Clean design: No logging in core functions, all diagnostics returned as data.
CLI layer is responsible for user output.
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass, field
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
class RinfEstimate:
    """Result of R_inf estimation."""
    R_inf: float
    method: str  # 'preset', 'median', 'rl_fit', 'voigt_fit'
    R_inf_median: Optional[float] = None  # For comparison
    figure: Optional[plt.Figure] = None

    # Fit diagnostics (if applicable)
    behavior: Optional[str] = None  # 'capacitive', 'inductive', 'mixed'
    n_points_used: Optional[int] = None
    R_squared: Optional[float] = None
    L_nH: Optional[float] = None
    R_ct: Optional[float] = None
    C_nF: Optional[float] = None
    f_characteristic: Optional[float] = None

    warnings: List[str] = field(default_factory=list)


@dataclass
class LambdaSelection:
    """Result of lambda selection."""
    lambda_value: float
    method: str  # 'user', 'default', 'gcv', 'hybrid', 'fallback'
    lambda_gcv: Optional[float] = None  # If L-curve correction was applied
    gcv_score: Optional[float] = None


@dataclass
class NNLSSolution:
    """Result of NNLS solver."""
    gamma: Optional[NDArray[np.float64]]
    success: bool
    n_inductive_points: int = 0
    inductive_fraction: float = 0.0
    max_inductive_imag: float = 0.0
    gamma_max: Optional[float] = None
    gamma_min_nonzero: Optional[float] = None
    warnings: List[str] = field(default_factory=list)


@dataclass
class DRTDiagnostics:
    """Comprehensive DRT diagnostics."""
    # Frequency info
    freq_min: float
    freq_max: float
    freq_range_ratio: float
    n_points: int

    # Matrix info
    n_tau: int
    condition_number: float
    d_ln_tau: float

    # R_inf estimation
    rinf: RinfEstimate

    # Lambda selection
    lambda_sel: LambdaSelection

    # NNLS solution
    nnls: NNLSSolution

    # R_pol info
    R_pol_from_data: float
    R_pol_from_gamma: float
    normalized: bool

    # Reconstruction
    reconstruction_error_rel: float

    # Peak detection
    peak_method: str
    n_peaks: int
    scipy_peaks: Optional[List[Dict]] = None  # For scipy method


@dataclass
class DRTMatrices:
    """Container for DRT matrices."""
    A: NDArray[np.float64]
    A_re: NDArray[np.float64]
    A_im: NDArray[np.float64]
    b: NDArray[np.float64]
    L: NDArray[np.float64]
    tau: NDArray[np.float64]
    d_ln_tau: float
    condition_number: float


@dataclass
class DRTResult:
    """
    Complete DRT analysis result.

    All information previously logged is now available as structured data.
    """
    # Core results
    tau: Optional[NDArray[np.float64]] = None
    gamma: Optional[NDArray[np.float64]] = None
    gamma_original: Optional[NDArray[np.float64]] = None

    # Figures
    figure: Optional[plt.Figure] = None
    figure_rinf: Optional[plt.Figure] = None

    # GMM peaks (if GMM method used)
    peaks: Optional[List[Dict]] = None
    bic_scores: Optional[List[float]] = None

    # Key values for quick access
    R_inf: Optional[float] = None
    R_pol: Optional[float] = None
    lambda_used: Optional[float] = None
    reconstruction_error: Optional[float] = None

    # Full diagnostics
    diagnostics: Optional[DRTDiagnostics] = None

    @property
    def success(self) -> bool:
        """Check if DRT calculation was successful."""
        return self.tau is not None and self.gamma is not None

    @property
    def warnings(self) -> List[str]:
        """Collect all warnings from diagnostics."""
        if self.diagnostics is None:
            return []
        warnings = []
        warnings.extend(self.diagnostics.rinf.warnings)
        warnings.extend(self.diagnostics.nnls.warnings)
        return warnings


# =============================================================================
# Helper Functions
# =============================================================================

def _estimate_peak_resistance(tau: NDArray, gamma: NDArray,
                               peak_indices: NDArray,
                               tolerance: float = 0.1) -> List[float]:
    """
    Estimate resistance for each peak by integrating gamma in peak region.
    """
    if len(peak_indices) == 0:
        return []

    resistances = []
    ln_tau = np.log(tau)

    for peak_idx in peak_indices:
        peak_height = gamma[peak_idx]
        threshold = peak_height * tolerance

        left = peak_idx
        while left > 0 and gamma[left] > threshold:
            left -= 1

        right = peak_idx
        while right < len(gamma) - 1 and gamma[right] > threshold:
            right += 1

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
                    r_inf_preset: Optional[float] = None) -> RinfEstimate:
    """
    Estimate high-frequency resistance R_inf.

    Returns structured RinfEstimate with all diagnostics.
    """
    n_avg = min(5, max(1, len(frequencies) // 10))
    high_freq_indices = np.argsort(frequencies)[-n_avg:]
    R_inf_median = float(np.median(Z.real[high_freq_indices]))

    if r_inf_preset is not None:
        return RinfEstimate(
            R_inf=r_inf_preset,
            method='preset',
            R_inf_median=R_inf_median
        )

    if use_rl_fit or use_voigt_fit:
        try:
            R_inf, L_fit, circuit_rl, diag, fig_rl = estimate_rinf_with_inductance(
                frequencies, Z, use_voigt_fit=use_voigt_fit, verbose=False, plot=True
            )

            warnings = []
            if diag.get('warnings'):
                warnings.extend(diag['warnings'])
            if diag.get('L_nH', 0) > 500:
                warnings.append(f"High inductance L = {diag['L_nH']:.1f} nH detected")

            return RinfEstimate(
                R_inf=R_inf,
                method='voigt_fit' if use_voigt_fit else 'rl_fit',
                R_inf_median=R_inf_median,
                figure=fig_rl,
                behavior=diag.get('behavior'),
                n_points_used=diag.get('n_points_used'),
                R_squared=diag.get('r_squared'),
                L_nH=diag.get('L_nH'),
                R_ct=diag.get('R_ct'),
                C_nF=diag.get('C_nF'),
                f_characteristic=diag.get('f_characteristic'),
                warnings=warnings
            )
        except Exception as e:
            logger.debug(f"R_inf fit failed: {e}, using median fallback")

    # Median method (default)
    return RinfEstimate(
        R_inf=R_inf_median,
        method='median',
        R_inf_median=R_inf_median,
        n_points_used=n_avg
    )


# =============================================================================
# Input Validation
# =============================================================================

def _validate_frequencies(frequencies: NDArray) -> List[str]:
    """
    Validate frequency array for DRT analysis.

    Returns list of warnings. Raises ValueError for critical errors.
    """
    warnings = []

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
        warnings.append(f"Small frequency range: {freq_range:.1f}x (recommended >{MIN_FREQUENCY_RANGE}x)")

    return warnings


# =============================================================================
# Matrix Construction
# =============================================================================

def _build_drt_matrices(frequencies: NDArray, Z: NDArray,
                        R_inf: float, n_tau: int = 100) -> DRTMatrices:
    """
    Build DRT system matrices A, b, and regularization matrix L.
    """
    f_max = frequencies.max()
    f_min = frequencies.min()

    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    omega = 2 * np.pi * frequencies

    d_ln_tau_array = np.diff(np.log(tau))
    d_ln_tau = float(np.mean(d_ln_tau_array))

    # Vectorized matrix construction
    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom

    A = np.vstack([A_re, A_im])
    cond_A = float(np.linalg.cond(A))

    b = np.concatenate([Z.real - R_inf, Z.imag])

    # Regularization matrix (2nd derivative)
    L = np.zeros((n_tau - 2, n_tau))
    for i in range(n_tau - 2):
        L[i, i] = 1
        L[i, i + 1] = -2
        L[i, i + 2] = 1

    return DRTMatrices(
        A=A, A_re=A_re, A_im=A_im, b=b, L=L,
        tau=tau, d_ln_tau=d_ln_tau, condition_number=cond_A
    )


# =============================================================================
# Lambda Selection
# =============================================================================

def _select_lambda(A: NDArray, b: NDArray, L: NDArray,
                   lambda_reg: Optional[float] = None,
                   auto_lambda: bool = False) -> LambdaSelection:
    """
    Select regularization parameter lambda.
    """
    if auto_lambda:
        try:
            lambda_opt, gcv_score, diag = find_optimal_lambda_hybrid(
                A, b, L,
                lambda_range=(1e-5, 1.0),
                n_search=20,
                lcurve_decades=1.5
            )
            if diag['method_used'] == 'lcurve_correction':
                return LambdaSelection(
                    lambda_value=lambda_opt,
                    method='hybrid',
                    lambda_gcv=diag['lambda_gcv'],
                    gcv_score=gcv_score
                )
            return LambdaSelection(
                lambda_value=lambda_opt,
                method='gcv',
                gcv_score=gcv_score
            )
        except Exception:
            try:
                lambda_opt, gcv_score = find_optimal_lambda_gcv(
                    A, b, L, lambda_range=(1e-5, 1.0), n_search=20
                )
                return LambdaSelection(
                    lambda_value=lambda_opt,
                    method='gcv',
                    gcv_score=gcv_score
                )
            except Exception:
                return LambdaSelection(lambda_value=0.1, method='fallback')

    if lambda_reg is None:
        return LambdaSelection(lambda_value=0.1, method='default')

    return LambdaSelection(lambda_value=lambda_reg, method='user')


# =============================================================================
# NNLS Solver
# =============================================================================

def _solve_nnls(A: NDArray, b: NDArray, L: NDArray,
                lambda_reg: float, n_tau: int,
                Z: NDArray) -> NNLSSolution:
    """
    Solve regularized NNLS problem.
    """
    warnings = []

    # Check for inductive data
    n_inductive = int(np.sum(Z.imag > 0))
    inductive_fraction = n_inductive / len(Z) if len(Z) > 0 else 0.0
    max_inductive = float(np.max(Z.imag)) if n_inductive > 0 else 0.0

    if inductive_fraction > 0.1 or max_inductive > 50:
        warnings.append(
            f"{n_inductive} points with inductive component ({inductive_fraction*100:.1f}%), "
            f"max Z'' = {max_inductive:.2f} Ohm"
        )

    # Build regularized system
    A_reg = np.vstack([A, np.sqrt(lambda_reg) * L])
    b_reg = np.concatenate([b, np.zeros(n_tau - 2)])

    # Solve NNLS
    try:
        gamma, residual = nnls(A_reg, b_reg)
    except Exception as e:
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=[f"NNLS solver error: {e}"]
        )

    # Validate solution
    if np.any(~np.isfinite(gamma)):
        return NNLSSolution(
            gamma=None, success=False,
            n_inductive_points=n_inductive,
            inductive_fraction=inductive_fraction,
            max_inductive_imag=max_inductive,
            warnings=["NNLS returned NaN or Inf values"]
        )

    gamma_max = float(np.max(gamma))
    gamma_nonzero = gamma[gamma > 0]
    gamma_min_nonzero = float(np.min(gamma_nonzero)) if len(gamma_nonzero) > 0 else None

    if gamma_max > GAMMA_MAX_REASONABLE:
        warnings.append(f"DRT contains very large values (max: {gamma_max:.2e} Ohm)")

    if np.sum(gamma) < 1e-10:
        warnings.append("DRT is nearly zero - data may have no relaxation structure")

    return NNLSSolution(
        gamma=gamma,
        success=True,
        n_inductive_points=n_inductive,
        inductive_fraction=inductive_fraction,
        max_inductive_imag=max_inductive,
        gamma_max=gamma_max,
        gamma_min_nonzero=gamma_min_nonzero,
        warnings=warnings
    )


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
    ax1.semilogx(tau, gamma, 'b-', linewidth=2, label='DRT gamma(tau)')
    ax1.fill_between(tau, 0, gamma, alpha=0.3)
    ax1.set_xlabel("tau [s]")

    if normalize_rpol:
        ax1.set_ylabel("gamma(tau) / R_pol [-]")
        ax1.set_title(f"DRT normalized (lambda = {lambda_reg})")
    else:
        ax1.set_ylabel("gamma(tau) [Ohm]")
        ax1.set_title(f"DRT (lambda = {lambda_reg})")

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
    ax2.set_xlabel("Z' [Ohm]")
    ax2.set_ylabel("-Z'' [Ohm]")
    ax2.set_title("DRT fit verification")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # === GMM Visualization ===
    if use_gmm:
        log_tau = np.log10(tau)
        ax3.semilogx(tau, gamma, 'b-', linewidth=2, label='DRT gamma(tau)', alpha=0.7)

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
                            label=f"Peak {i+1}: tau={peak['tau_center']:.2e}s")
            ax3.axvline(peak['tau_bounds'][0], color=colors[i], linestyle=':', alpha=0.8)
            ax3.axvline(peak['tau_bounds'][1], color=colors[i], linestyle=':', alpha=0.8)
            ax3.axvline(peak['tau_center'], color=colors[i], linestyle='--', alpha=0.8)

        ax3.set_xlabel("tau [s]")
        ax3.set_ylabel("gamma(tau) [Ohm]")
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
                  gmm_bic_threshold: float = 10.0
                  ) -> Tuple[Optional[List[Dict]], Optional[List[float]], Optional[List[Dict]]]:
    """
    Detect peaks in DRT spectrum.

    Returns:
        (gmm_peaks, bic_scores, scipy_peaks)
    """
    use_gmm = (peak_method == 'gmm' and GMM_AVAILABLE)

    # Always calculate scipy peaks for diagnostics
    peaks_idx, _ = find_peaks(gamma, height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD)
    peak_resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)

    scipy_peaks = []
    for i, idx in enumerate(peaks_idx):
        R_peak = peak_resistances[i] if i < len(peak_resistances) else 0.0
        scipy_peaks.append({
            'index': int(idx),
            'tau': float(tau[idx]),
            'frequency': float(1/tau[idx]),
            'R_estimate': float(R_peak)
        })

    if use_gmm:
        peaks_result, gmm_model, bic_scores = gmm_peak_detection(
            tau, gamma, bic_threshold=gmm_bic_threshold
        )

        if len(peaks_result) == 0 or gmm_model is None:
            # GMM failed, scipy_peaks available as fallback
            return None, None, scipy_peaks

        return peaks_result, bic_scores, scipy_peaks

    return None, None, scipy_peaks


# =============================================================================
# Main Function
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
) -> DRTResult:
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
    gmm_bic_threshold : float
        BIC threshold for GMM peak detection

    Returns
    -------
    DRTResult
        Complete analysis result with all diagnostics
    """
    f_min, f_max = float(frequencies.min()), float(frequencies.max())
    freq_range_ratio = f_max / f_min

    # === Step 1: Validate frequencies ===
    try:
        freq_warnings = _validate_frequencies(frequencies)
    except ValueError as e:
        logger.error(str(e))
        return DRTResult()

    # === Step 2: R_inf Estimation ===
    rinf_est = _estimate_r_inf(
        frequencies, Z,
        use_rl_fit=use_rl_fit,
        use_voigt_fit=use_voigt_fit,
        r_inf_preset=r_inf_preset
    )
    R_inf = rinf_est.R_inf

    # === Step 3: Build Matrices ===
    matrices = _build_drt_matrices(frequencies, Z, R_inf, n_tau)

    # === Step 4: Select Lambda ===
    lambda_sel = _select_lambda(
        matrices.A, matrices.b, matrices.L,
        lambda_reg, auto_lambda
    )

    # === Step 5: Solve NNLS ===
    nnls_result = _solve_nnls(
        matrices.A, matrices.b, matrices.L,
        lambda_sel.lambda_value, n_tau, Z
    )

    if not nnls_result.success:
        return DRTResult(
            R_inf=R_inf,
            figure_rinf=rinf_est.figure,
            diagnostics=DRTDiagnostics(
                freq_min=f_min, freq_max=f_max,
                freq_range_ratio=freq_range_ratio,
                n_points=len(frequencies),
                n_tau=n_tau,
                condition_number=matrices.condition_number,
                d_ln_tau=matrices.d_ln_tau,
                rinf=rinf_est,
                lambda_sel=lambda_sel,
                nnls=nnls_result,
                R_pol_from_data=0.0,
                R_pol_from_gamma=0.0,
                normalized=False,
                reconstruction_error_rel=0.0,
                peak_method=peak_method,
                n_peaks=0
            )
        )

    gamma = nnls_result.gamma

    # === Step 6: R_pol Calculation & Normalization ===
    n_avg = min(5, max(1, len(frequencies) // 10))
    low_freq_indices = np.argsort(frequencies)[:n_avg]
    R_dc = float(np.median(Z.real[low_freq_indices]))
    R_pol_from_data = R_dc - R_inf
    R_pol_from_gamma = float(np_trapz(gamma, np.log(matrices.tau)))

    gamma_original = None
    normalized = False
    if normalize_rpol:
        if R_pol_from_gamma > 1e-10:
            gamma_original = gamma.copy()
            gamma = gamma / R_pol_from_gamma
            normalized = True

    # === Step 7: Peak Detection ===
    peaks_result, bic_scores, scipy_peaks = _detect_peaks(
        matrices.tau, gamma, peak_method, gmm_bic_threshold
    )

    n_peaks = len(peaks_result) if peaks_result else len(scipy_peaks) if scipy_peaks else 0

    # === Step 8: Reconstruction & Error ===
    gamma_for_recon = gamma_original if normalized else gamma
    Z_reconstructed = R_inf + (matrices.A_re + 1j * matrices.A_im) @ gamma_for_recon
    rel_error = float(np.mean(np.abs(Z - Z_reconstructed) / np.abs(Z)) * 100)

    # Add warning for high reconstruction error
    if rel_error > 10.0:
        nnls_result.warnings.append(
            f"High reconstruction error ({rel_error:.1f}%) - DRT model may not be suitable"
        )
    elif rel_error > 5.0:
        nnls_result.warnings.append(
            f"Elevated reconstruction error ({rel_error:.1f}%)"
        )

    # === Step 9: Visualization ===
    fig = _create_visualization(
        matrices.tau, gamma, gamma_original,
        Z, Z_reconstructed,
        R_inf, lambda_sel.lambda_value,
        normalized, peak_method,
        peaks_result, bic_scores
    )

    # === Build diagnostics ===
    diagnostics = DRTDiagnostics(
        freq_min=f_min,
        freq_max=f_max,
        freq_range_ratio=freq_range_ratio,
        n_points=len(frequencies),
        n_tau=n_tau,
        condition_number=matrices.condition_number,
        d_ln_tau=matrices.d_ln_tau,
        rinf=rinf_est,
        lambda_sel=lambda_sel,
        nnls=nnls_result,
        R_pol_from_data=R_pol_from_data,
        R_pol_from_gamma=R_pol_from_gamma,
        normalized=normalized,
        reconstruction_error_rel=rel_error,
        peak_method=peak_method,
        n_peaks=n_peaks,
        scipy_peaks=scipy_peaks
    )

    return DRTResult(
        tau=matrices.tau,
        gamma=gamma,
        gamma_original=gamma_original,
        figure=fig,
        figure_rinf=rinf_est.figure,
        peaks=peaks_result,
        bic_scores=bic_scores,
        R_inf=R_inf,
        R_pol=R_pol_from_gamma,
        lambda_used=lambda_sel.lambda_value,
        reconstruction_error=rel_error,
        diagnostics=diagnostics
    )
