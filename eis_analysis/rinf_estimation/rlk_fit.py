"""
R-L-K() linear fit for R_inf estimation.

Simple and robust approach: fit R_s - L - K(R, tau) model to highest decade
using linear least squares. Works for both inductive and capacitive data.

Model: Z(omega) = R_s + j*omega*L + R_k / (1 + j*omega*tau)

Clean design: No logging in core functions, all diagnostics returned as data.
"""

import numpy as np
import logging
from dataclasses import dataclass, field
from typing import Optional, List
from numpy.typing import NDArray

from .data_selection import select_highest_decade
from ..utils.impedance import sort_by_frequency

# Lazy import to avoid circular dependency
_estimate_R_linear = None


def _get_estimate_R_linear():
    """Lazy import to avoid circular dependency."""
    global _estimate_R_linear
    if _estimate_R_linear is None:
        from ..fitting.voigt_chain import estimate_R_linear
        _estimate_R_linear = estimate_R_linear
    return _estimate_R_linear


logger = logging.getLogger(__name__)


@dataclass
class RLKFitResult:
    """Result of R-L-K() linear fit."""
    # Core results
    R_inf: float
    L: float  # [H]
    R_k: float  # [Ohm]
    tau: float  # [s]
    Z_fit: NDArray[np.complex128]

    # Quality metrics
    R_squared: float
    rel_error: float  # [%]
    residual: float

    # Data info
    n_points_used: int
    freq_range: tuple  # (f_min, f_max) [Hz]
    behavior: str  # 'purely_capacitive', 'purely_inductive', 'mixed_with_crossing'
    method: str

    # Derived values
    L_nH: float  # Inductance [nH]
    tau_us: float  # Time constant [us]
    f_char_kHz: float  # Characteristic frequency [kHz]

    # Optional extras
    R_inf_hf: Optional[float] = None  # R_inf at highest frequency
    R_inf_poly: Optional[float] = None  # R_inf from polynomial fit
    f_zero_crossing: Optional[float] = None  # Frequency of Im(Z)=0 crossing
    poly_coeffs: Optional[list] = None

    # Status
    fit_success: bool = True
    warnings: List[str] = field(default_factory=list)


def _estimate_tau_from_data(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128]
) -> float:
    """
    Estimate tau from data based on where -Im(Z) is maximum.

    For a Voigt element, -Im(Z) peaks at omega = 1/tau.
    """
    im_abs = np.abs(Z.imag)
    idx_max = np.argmax(im_abs)
    is_at_boundary = (idx_max == 0) or (idx_max == len(frequencies) - 1)

    if is_at_boundary:
        f_geo = np.sqrt(frequencies.min() * frequencies.max())
        tau_est = 1.0 / (2 * np.pi * f_geo)
    else:
        f_peak = frequencies[idx_max]
        tau_est = 1.0 / (2 * np.pi * f_peak)

    return tau_est


def fit_rlk_model(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128]
) -> RLKFitResult:
    """
    Fit R-L-K() model to high-frequency EIS data using linear least squares.

    Model: Z(omega) = R_s + j*omega*L + R_k / (1 + j*omega*tau)

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (highest decade will be selected)
    Z : ndarray of complex
        Complex impedance [Ohm]

    Returns
    -------
    RLKFitResult
        Complete fit result with all diagnostics

    Notes
    -----
    The model handles both cases naturally:
    - Inductive data (Im(Z) > 0): L > 0, R_k ~ 0
    - Capacitive data (Im(Z) < 0): L ~ 0, R_k > 0
    """
    # Sort by frequency
    frequencies, Z = sort_by_frequency(frequencies, Z)

    # Select highest decade
    selection = select_highest_decade(frequencies, Z)
    f_high = selection.frequencies
    Z_high = selection.Z
    method_str = selection.method
    warnings = list(selection.warnings)

    # Detect behavior
    n_inductive = selection.n_inductive
    n_capacitive = selection.n_capacitive
    is_purely_capacitive = n_inductive == 0
    is_purely_inductive = n_capacitive == 0
    has_zero_crossing = n_inductive > 0 and n_capacitive > 0

    # Determine behavior string
    if is_purely_capacitive:
        behavior = "purely_capacitive"
    elif is_purely_inductive:
        behavior = "purely_inductive"
    elif has_zero_crossing:
        behavior = "mixed_with_crossing"
    else:
        behavior = "inductive" if n_inductive > n_capacitive else "capacitive"

    # BEST CASE: Im(Z) crosses zero - interpolate Re(Z) at crossing
    if has_zero_crossing:
        sort_idx = np.argsort(f_high)
        f_sorted = f_high[sort_idx]
        Z_sorted = Z_high[sort_idx]
        im_sorted = Z_sorted.imag

        for i in range(len(im_sorted) - 1):
            if im_sorted[i] * im_sorted[i+1] < 0:  # Sign change
                f1, f2 = f_sorted[i], f_sorted[i+1]
                im1, im2 = im_sorted[i], im_sorted[i+1]
                re1, re2 = Z_sorted[i].real, Z_sorted[i+1].real

                t = -im1 / (im2 - im1)
                f_zero = f1 + t * (f2 - f1)
                R_inf = re1 + t * (re2 - re1)

                Z_fit = np.full_like(Z_high, R_inf)
                ss_res = np.sum((Z_high.real - R_inf)**2)
                ss_tot = np.sum((Z_high.real - np.mean(Z_high.real))**2)
                R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
                rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(Z_high.real**2)) * 100

                tau_val = 1.0 / (2 * np.pi * f_zero)

                return RLKFitResult(
                    R_inf=R_inf,
                    L=0.0,
                    R_k=0.0,
                    tau=tau_val,
                    Z_fit=Z_fit,
                    R_squared=R_squared,
                    rel_error=rel_error,
                    residual=np.sqrt(ss_res),
                    n_points_used=len(f_high),
                    freq_range=(float(f_high.min()), float(f_high.max())),
                    behavior=behavior,
                    method=f'zero_crossing_{method_str}',
                    L_nH=0.0,
                    tau_us=tau_val * 1e6,
                    f_char_kHz=f_zero / 1e3,
                    f_zero_crossing=f_zero,
                    warnings=warnings
                )

    # For purely capacitive data: polynomial extrapolation
    if is_purely_capacitive:
        im_hf = Z_high.imag
        re_hf = Z_high.real
        coeffs = np.polyfit(im_hf, re_hf, 2)
        R_inf_poly = coeffs[2]

        idx_max = np.argmax(f_high)
        R_inf_hf = float(Z_high[idx_max].real)

        if R_inf_poly <= 0:
            R_inf = 0.0
            extrap_method = "zero_fallback"
        elif R_inf_poly > R_inf_hf:
            R_inf = R_inf_hf
            extrap_method = "hf_fallback"
        else:
            R_inf = float(R_inf_poly)
            extrap_method = "polynomial"

        tau_val = 1.0 / (2 * np.pi * f_high.max())
        Z_fit = np.full_like(Z_high, R_inf)

        ss_res = np.sum((Z_high.real - R_inf)**2)
        ss_tot = np.sum((Z_high.real - np.mean(Z_high.real))**2)
        R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(Z_high.real**2)) * 100

        return RLKFitResult(
            R_inf=R_inf,
            L=0.0,
            R_k=0.0,
            tau=tau_val,
            Z_fit=Z_fit,
            R_squared=R_squared,
            rel_error=rel_error,
            residual=np.sqrt(ss_res),
            n_points_used=len(f_high),
            freq_range=(float(f_high.min()), float(f_high.max())),
            behavior=behavior,
            method=f'capacitive_{extrap_method}_{method_str}',
            L_nH=0.0,
            tau_us=tau_val * 1e6,
            f_char_kHz=float(f_high.max()) / 1e3,
            R_inf_hf=R_inf_hf,
            R_inf_poly=float(R_inf_poly),
            poly_coeffs=coeffs.tolist(),
            warnings=warnings
        )

    # For inductive data: R-L-K model
    tau_est = _estimate_tau_from_data(f_high, Z_high)

    f_max = f_high.max()
    f_min = f_high.min()
    tau_min = 1.0 / (2 * np.pi * f_max)
    tau_max = 1.0 / (2 * np.pi * f_min)
    tau_val = float(np.clip(tau_est, tau_min, tau_max))

    tau = np.array([tau_val])

    # Linear fit
    estimate_R_linear = _get_estimate_R_linear()
    elements, residual, L_value = estimate_R_linear(
        f_high, Z_high, tau,
        include_Rs=True,
        include_L=True,
        fit_type='complex',
        allow_negative=False,
        weighting='modulus'
    )

    R_inf = float(elements[0])
    R_k = float(elements[1])
    L = float(L_value) if L_value is not None else 0.0

    # Compute fitted impedance
    omega = 2 * np.pi * f_high
    Z_fit = R_inf + 1j * omega * L + R_k / (1 + 1j * omega * tau_val)

    # Quality metrics
    Z_residual = Z_high - Z_fit
    ss_res = np.sum(np.abs(Z_residual)**2)
    ss_tot = np.sum(np.abs(Z_high - np.mean(Z_high))**2)
    R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(np.abs(Z_high)**2)) * 100

    return RLKFitResult(
        R_inf=R_inf,
        L=L,
        R_k=R_k,
        tau=tau_val,
        Z_fit=Z_fit,
        R_squared=R_squared,
        rel_error=rel_error,
        residual=float(residual),
        n_points_used=len(f_high),
        freq_range=(float(f_high.min()), float(f_high.max())),
        behavior=behavior,
        method=f'rlk_linear_{method_str}',
        L_nH=L * 1e9,
        tau_us=tau_val * 1e6,
        f_char_kHz=1 / (2 * np.pi * tau_val) / 1e3,
        warnings=warnings
    )


def estimate_rinf_with_inductance(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    max_L_nH: float = 1000.0,
    verbose: bool = True,  # Kept for backward compatibility, ignored
    plot: bool = False,
    save_plot: Optional[str] = None,
    use_voigt_fit: bool = False  # Kept for backward compatibility, ignored
) -> tuple:
    """
    Estimate R_inf using R-L-K() linear fit.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    max_L_nH : float
        Maximum reasonable inductance [nH] for warning (default: 1000)
    verbose : bool
        Ignored (kept for backward compatibility)
    plot : bool
        Create diagnostic plot (default: False)
    save_plot : str, optional
        Path to save plot
    use_voigt_fit : bool
        Ignored (kept for backward compatibility)

    Returns
    -------
    R_inf : float
        High-frequency resistance [Ohm]
    L : float
        Inductance [H]
    circuit : None
        Always None (legacy placeholder)
    diagnostics : dict
        Fitting diagnostics
    fig : matplotlib.Figure or None
        Diagnostic plot if plot=True
    """
    # Validate input
    if len(frequencies) < 3:
        diagnostics = {
            'n_points_used': len(frequencies),
            'fit_success': False,
            'warnings': ['Insufficient data (< 3 points)'],
            'method': 'fallback_insufficient_data',
            'behavior': 'unknown'
        }
        R_inf = float(np.median(Z.real)) if len(Z) > 0 else 0.0
        return R_inf, 0.0, None, diagnostics, None

    # Perform fit
    try:
        result = fit_rlk_model(frequencies, Z)
    except Exception as e:
        diagnostics = {
            'fit_success': False,
            'warnings': [str(e)],
            'method': 'fallback_fit_error',
            'behavior': 'unknown'
        }
        R_inf = float(np.median(Z.real))
        return R_inf, 0.0, None, diagnostics, None

    # Validate inductance
    if result.L_nH > max_L_nH:
        result.warnings.append(f"High inductance L = {result.L_nH:.1f} nH > {max_L_nH:.0f} nH")

    if result.L < 0:
        result.warnings.append(f"Negative inductance L = {result.L_nH:.1f} nH (non-physical)")

    # Build diagnostics dict for backward compatibility
    diagnostics = {
        'n_points_used': result.n_points_used,
        'freq_range': result.freq_range,
        'R_squared': result.R_squared,
        'r_squared': result.R_squared,  # Alias
        'rel_error': result.rel_error,
        'R_inf': result.R_inf,
        'L': result.L,
        'L_nH': result.L_nH,
        'R_k': result.R_k,
        'tau': result.tau,
        'tau_us': result.tau_us,
        'f_char_kHz': result.f_char_kHz,
        'method': result.method,
        'residual': result.residual,
        'warnings': result.warnings,
        'fit_success': result.fit_success,
        'behavior': result.behavior
    }

    # Add optional fields if present
    if result.R_inf_hf is not None:
        diagnostics['R_inf_hf'] = result.R_inf_hf
    if result.R_inf_poly is not None:
        diagnostics['R_inf_poly'] = result.R_inf_poly
    if result.f_zero_crossing is not None:
        diagnostics['f_zero_crossing'] = result.f_zero_crossing
    if result.poly_coeffs is not None:
        diagnostics['poly_coeffs'] = result.poly_coeffs

    # Plot if requested
    fig = None
    if plot:
        try:
            fig = _plot_rlk_fit(frequencies, Z, result, save_plot)
        except Exception as e:
            logger.debug(f"Visualization failed: {e}")

    return result.R_inf, result.L, None, diagnostics, fig


def _plot_rlk_fit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    result: RLKFitResult,
    save_path: Optional[str] = None
):
    """Create diagnostic plot for R-L-K() fit."""
    import matplotlib.pyplot as plt

    frequencies, Z = sort_by_frequency(frequencies, Z)

    f_max = frequencies.max()
    f_min_decade = f_max / 10.0
    mask = frequencies >= f_min_decade

    f_fit = frequencies[mask]
    Z_fit_data = Z[mask]

    omega_fit = 2 * np.pi * f_fit
    Z_fit = (result.R_inf + 1j * omega_fit * result.L +
             result.R_k / (1 + 1j * omega_fit * result.tau))

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Nyquist plot
    ax = axes[0, 0]
    ax.plot(Z_fit_data.real, -Z_fit_data.imag, 'o', color='#1f77b4',
            ms=6, alpha=0.8, label=f'HF data ({len(f_fit)} pts)',
            markeredgewidth=0.5, markeredgecolor='white')
    ax.plot(Z_fit.real, -Z_fit.imag, 'r-', lw=2, label='R-L-K() fit')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.axvline(result.R_inf, color='green', ls='--', lw=1.5,
               label=f'R_inf = {result.R_inf:.3f} Ohm')
    ax.set_xlabel("Z' [Ohm]")
    ax.set_ylabel("-Z'' [Ohm]")
    ax.set_title('Nyquist Plot (highest decade)')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. Re(Z) vs frequency
    ax = axes[0, 1]
    ax.semilogx(f_fit, Z_fit_data.real, 'o', color='#1f77b4',
                ms=6, alpha=0.8, label='HF data',
                markeredgewidth=0.5, markeredgecolor='white')
    ax.semilogx(f_fit, Z_fit.real, 'r-', lw=2, label='Fit')
    ax.axhline(result.R_inf, color='green', ls='--', lw=1.5,
               label=f'R_inf = {result.R_inf:.3f} Ohm')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("Re(Z) [Ohm]")
    ax.set_title('Re(Z) vs Frequency')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    # 3. Im(Z) vs frequency
    ax = axes[1, 0]
    ax.semilogx(f_fit, Z_fit_data.imag, 'o', color='#1f77b4',
                ms=6, alpha=0.8, label='HF data',
                markeredgewidth=0.5, markeredgecolor='white')
    ax.semilogx(f_fit, Z_fit.imag, 'r-', lw=2, label='Fit')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Im(Z) [Ohm]')
    ax.set_title('Im(Z) vs Frequency')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    # 4. Info panel
    ax = axes[1, 1]
    ax.axis('off')

    info_text = [
        "R-L-K() FIT RESULTS",
        "=" * 30,
        "",
        f"R_inf = {result.R_inf:.4f} Ohm",
        f"L     = {result.L_nH:.2f} nH",
        f"R_k   = {result.R_k:.4f} Ohm",
        f"tau   = {result.tau_us:.3f} us",
        "",
        "QUALITY METRICS",
        "-" * 20,
        f"R^2       = {result.R_squared:.4f}",
        f"Rel.error = {result.rel_error:.2f}%",
        f"Points    = {result.n_points_used}",
        "",
        f"Behavior: {result.behavior}",
        f"Method: {result.method}"
    ]

    if result.warnings:
        info_text.append("")
        info_text.append("WARNINGS:")
        for w in result.warnings:
            info_text.append(f"  * {w}")

    ax.text(0.1, 0.95, '\n'.join(info_text), transform=ax.transAxes,
            fontsize=10, fontfamily='monospace', verticalalignment='top')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')

    return fig


__all__ = ['fit_rlk_model', 'estimate_rinf_with_inductance', 'RLKFitResult']
