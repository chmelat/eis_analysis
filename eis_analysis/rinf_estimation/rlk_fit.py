"""
R-L-K() linear fit for R_inf estimation.

Simple and robust approach: fit R_s - L - K(R, τ) model to highest decade
using linear least squares. Works for both inductive and capacitive data.

Model: Z(ω) = R_s + jωL + R_k / (1 + jωτ)

Where:
- R_s = R_inf (series/solution resistance)
- L = parasitic inductance (can be ~0 for capacitive data)
- K(R_k, τ) = single Voigt element capturing initial arc

The τ is fixed to 1/(2π·f_max), making the problem linear in [R_s, R_k, L].
"""

import numpy as np
import logging
from typing import Tuple, Optional
from numpy.typing import NDArray

from .data_selection import select_highest_decade
from ..utils.impedance import sort_by_frequency

# Import estimate_R_linear lazily to avoid circular import
# (fitting.voigt_chain → fitting.__init__ → auto_suggest → drt → rinf_estimation)
_estimate_R_linear = None

def _get_estimate_R_linear():
    """Lazy import to avoid circular dependency."""
    global _estimate_R_linear
    if _estimate_R_linear is None:
        from ..fitting.voigt_chain import estimate_R_linear
        _estimate_R_linear = estimate_R_linear
    return _estimate_R_linear

logger = logging.getLogger(__name__)


def _estimate_tau_from_data(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128]
) -> float:
    """
    Estimate τ from data based on where -Im(Z) is maximum.

    For a Voigt element, -Im(Z) peaks at ω = 1/τ.
    If no clear peak exists (monotonic behavior), uses geometric mean frequency.
    """
    # Find frequency where |Im(Z)| is maximum
    im_abs = np.abs(Z.imag)
    idx_max = np.argmax(im_abs)

    # Check if maximum is at boundary (indicates monotonic behavior, no peak)
    is_at_boundary = (idx_max == 0) or (idx_max == len(frequencies) - 1)

    if is_at_boundary:
        # No clear peak - use geometric mean of frequency range
        # This gives τ in the middle of the decade (log scale)
        f_geo = np.sqrt(frequencies.min() * frequencies.max())
        tau_est = 1.0 / (2 * np.pi * f_geo)
        logger.debug(f"No Im(Z) peak found, using geometric mean f = {f_geo:.0f} Hz")
    else:
        # Clear peak - use peak frequency
        f_peak = frequencies[idx_max]
        tau_est = 1.0 / (2 * np.pi * f_peak)
        logger.debug(f"Im(Z) peak at f = {f_peak:.0f} Hz")

    return tau_est


def fit_rlk_model(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    verbose: bool = True
) -> Tuple[float, float, float, float, NDArray[np.complex128], dict]:
    """
    Fit R-L-K() model to high-frequency EIS data using linear least squares.

    Model: Z(ω) = R_s + jωL + R_k / (1 + jωτ)

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz] (highest decade will be selected)
    Z : ndarray of complex
        Complex impedance [Ω]
    verbose : bool
        Log fitting progress (default: True)

    Returns
    -------
    R_inf : float
        High-frequency resistance R_s [Ω]
    L : float
        Inductance [H]
    R_k : float
        Voigt element resistance [Ω]
    tau : float
        Voigt element time constant [s]
    Z_fit : ndarray of complex
        Fitted impedance for selected frequencies [Ω]
    diagnostics : dict
        Fitting diagnostics including:
        - 'n_points': number of points used
        - 'freq_range': (f_min, f_max) Hz
        - 'residual': fit residual norm
        - 'R_squared': coefficient of determination
        - 'L_nH': inductance in nH
        - 'tau_us': time constant in μs
        - 'f_char_kHz': characteristic frequency in kHz

    Notes
    -----
    The model handles both cases naturally:
    - Inductive data (Im(Z) > 0): L > 0, R_k ≈ 0
    - Capacitive data (Im(Z) < 0): L ≈ 0, R_k > 0

    Linear least squares is used (NNLS for R_s, R_k ≥ 0), making the fit:
    - Fast (analytical solution)
    - Robust (no initial guess needed)
    - Stable (no convergence issues)

    The τ is estimated from data (frequency of max |Im(Z)|) rather than
    fixed to f_max, giving better fit for capacitive data.

    Examples
    --------
    >>> R_inf, L, R_k, tau, Z_fit, diag = fit_rlk_model(freq, Z)
    >>> print(f"R_inf = {R_inf:.3f} Ω, L = {diag['L_nH']:.1f} nH")
    """
    # Sort by frequency
    frequencies, Z = sort_by_frequency(frequencies, Z)

    # Select highest decade
    f_high, Z_high, method_str, warnings = select_highest_decade(
        frequencies, Z, verbose=verbose
    )

    if verbose:
        logger.info("")
        logger.info("=" * 60)
        logger.info("R-L-K() LINEAR FIT")
        logger.info("=" * 60)
        logger.info(f"Data: {len(f_high)} points, "
                   f"{f_high.min()/1e3:.1f} - {f_high.max()/1e3:.1f} kHz")

    # Detect behavior: inductive (Im > 0) vs capacitive (Im < 0)
    n_inductive = np.sum(Z_high.imag > 0)
    n_capacitive = np.sum(Z_high.imag < 0)
    is_mostly_inductive = n_inductive > n_capacitive
    is_purely_capacitive = n_inductive == 0  # All points have Im(Z) < 0
    is_purely_inductive = n_capacitive == 0  # All points have Im(Z) > 0
    has_zero_crossing = n_inductive > 0 and n_capacitive > 0

    if verbose:
        if is_purely_capacitive:
            behavior = "purely_capacitive"
        elif is_purely_inductive:
            behavior = "purely_inductive"
        elif has_zero_crossing:
            behavior = "mixed_with_crossing"
        else:
            behavior = "inductive" if is_mostly_inductive else "capacitive"
        logger.info(f"Behavior: {behavior} ({n_inductive} inductive, {n_capacitive} capacitive points)")

    # BEST CASE: Im(Z) crosses zero - interpolate Re(Z) at crossing
    # This gives the most accurate R_inf
    if has_zero_crossing:
        if verbose:
            logger.info("Im(Z) crosses zero - interpolating R_inf at crossing")

        # Sort by frequency for interpolation
        sort_idx = np.argsort(f_high)
        f_sorted = f_high[sort_idx]
        Z_sorted = Z_high[sort_idx]

        # Find zero crossing (Im changes sign)
        im_sorted = Z_sorted.imag
        for i in range(len(im_sorted) - 1):
            if im_sorted[i] * im_sorted[i+1] < 0:  # Sign change
                # Linear interpolation
                f1, f2 = f_sorted[i], f_sorted[i+1]
                im1, im2 = im_sorted[i], im_sorted[i+1]
                re1, re2 = Z_sorted[i].real, Z_sorted[i+1].real

                # f at Im(Z) = 0
                t = -im1 / (im2 - im1)  # interpolation parameter
                f_zero = f1 + t * (f2 - f1)
                R_inf = re1 + t * (re2 - re1)

                if verbose:
                    logger.info(f"  Crossing at f = {f_zero/1e6:.3f} MHz")
                    logger.info(f"  R_inf = {R_inf:.4f} Ω (interpolated)")

                L = 0.0
                R_k = 0.0
                tau_val = 1.0 / (2 * np.pi * f_zero)
                Z_fit = np.full_like(Z_high, R_inf)

                # Quality metric
                ss_res = np.sum((Z_high.real - R_inf)**2)
                ss_tot = np.sum((Z_high.real - np.mean(Z_high.real))**2)
                R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
                rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(Z_high.real**2)) * 100

                diagnostics = {
                    'n_points_used': len(f_high),
                    'freq_range': (f_high.min(), f_high.max()),
                    'R_squared': R_squared,
                    'rel_error': rel_error,
                    'R_inf': R_inf,
                    'L': L,
                    'L_nH': 0.0,
                    'R_k': R_k,
                    'tau': tau_val,
                    'tau_us': tau_val * 1e6,
                    'f_char_kHz': f_zero / 1e3,
                    'method': f'zero_crossing_{method_str}',
                    'residual': np.sqrt(ss_res),
                    'warnings': warnings,
                    'fit_success': True,
                    'behavior': 'mixed_with_crossing',
                    'f_zero_crossing': f_zero
                }

                if verbose:
                    logger.info("")
                    logger.info(f"Quality: R² = {R_squared:.4f}, rel_error = {rel_error:.2f}%")
                    logger.info("=" * 60)

                return R_inf, L, R_k, tau_val, Z_fit, diagnostics

    # For purely capacitive data (no inductive points at all):
    # Use polynomial extrapolation in Nyquist space
    if is_purely_capacitive:
        if verbose:
            logger.info("Pure capacitive: polynomial extrapolation in Nyquist space")

        # Fit 2nd degree polynomial: Z_real = a*Z_imag² + b*Z_imag + c
        # At Z_imag = 0: R_inf = c
        im_hf = Z_high.imag
        re_hf = Z_high.real
        coeffs = np.polyfit(im_hf, re_hf, 2)
        R_inf_poly = coeffs[2]  # c = R_inf extrapolated

        # Fallback: if polynomial gives negative or unreasonable value,
        # use Re(Z) at highest frequency
        idx_max = np.argmax(f_high)
        R_inf_hf = Z_high[idx_max].real

        if R_inf_poly <= 0:
            R_inf = 0.0
            extrap_method = "zero_fallback"
            if verbose:
                logger.info(f"  Polynomial gave {R_inf_poly:.2f} Ω (negative), using R_inf = 0")
        elif R_inf_poly > R_inf_hf:
            R_inf = R_inf_hf
            extrap_method = "hf_fallback"
            if verbose:
                logger.info(f"  Polynomial gave {R_inf_poly:.2f} Ω (> R @ f_max), using HF fallback")
        else:
            R_inf = R_inf_poly
            extrap_method = "polynomial"

        L = 0.0
        R_k = 0.0
        tau_val = 1.0 / (2 * np.pi * f_high.max())
        Z_fit = np.full_like(Z_high, R_inf)

        ss_res = np.sum((Z_high.real - R_inf)**2)
        ss_tot = np.sum((Z_high.real - np.mean(Z_high.real))**2)
        R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(Z_high.real**2)) * 100

        diagnostics = {
            'n_points_used': len(f_high),
            'freq_range': (f_high.min(), f_high.max()),
            'R_squared': R_squared,
            'rel_error': rel_error,
            'R_inf': R_inf,
            'R_inf_hf': R_inf_hf,
            'R_inf_poly': R_inf_poly,
            'poly_coeffs': coeffs.tolist(),
            'L': L,
            'L_nH': 0.0,
            'R_k': R_k,
            'tau': tau_val,
            'tau_us': tau_val * 1e6,
            'f_char_kHz': f_high.max() / 1e3,
            'method': f'capacitive_{extrap_method}_{method_str}',
            'residual': np.sqrt(ss_res),
            'warnings': warnings,
            'fit_success': True,
            'behavior': 'purely_capacitive'
        }

        if verbose:
            logger.info("")
            logger.info("Results:")
            logger.info(f"  R_inf = {R_inf:.4f} Ω ({extrap_method})")
            logger.info(f"  R_inf @ f_max = {R_inf_hf:.4f} Ω")
            logger.info(f"  R_inf polynom = {R_inf_poly:.4f} Ω")
            logger.info("  L     = 0 nH (capacitive data)")
            logger.info("")
            logger.info(f"Quality: R² = {R_squared:.4f}, rel_error = {rel_error:.2f}%")
            logger.info("=" * 60)

        return R_inf, L, R_k, tau_val, Z_fit, diagnostics

    # For purely inductive data: use R-L-K model
    include_L = True  # Always include L for inductive data

    if verbose:
        logger.info("Model: R-L-K()")

    # Estimate τ from data (where |Im(Z)| is maximum)
    tau_est = _estimate_tau_from_data(f_high, Z_high)

    # Clamp τ to reasonable range within the frequency band
    f_max = f_high.max()
    f_min = f_high.min()
    tau_min = 1.0 / (2 * np.pi * f_max)
    tau_max = 1.0 / (2 * np.pi * f_min)
    tau_val = np.clip(tau_est, tau_min, tau_max)

    tau = np.array([tau_val])

    if verbose:
        logger.info(f"Estimated τ = {tau_val*1e6:.3f} μs "
                   f"(f_char = {1/(2*np.pi*tau_val)/1e3:.1f} kHz)")

    # Linear fit: [R_s, R_k] or [R_s, R_k, L]
    # estimate_R_linear returns elements in order: [R_s, R_1, ..., R_N, (L)]
    estimate_R_linear = _get_estimate_R_linear()
    elements, residual, L_value = estimate_R_linear(
        f_high, Z_high, tau,
        include_Rs=True,
        include_L=include_L,
        fit_type='complex',
        allow_negative=False,  # NNLS: R_s, R_k ≥ 0
        weighting='proportional'
    )

    # Extract parameters
    R_inf = elements[0]  # R_s
    R_k = elements[1]    # R from K element
    L = L_value if (include_L and L_value is not None) else 0.0
    tau_val = tau[0]

    # Compute fitted impedance
    omega = 2 * np.pi * f_high
    Z_fit = R_inf + 1j * omega * L + R_k / (1 + 1j * omega * tau_val)

    # Quality metrics
    Z_residual = Z_high - Z_fit
    ss_res = np.sum(np.abs(Z_residual)**2)
    ss_tot = np.sum(np.abs(Z_high - np.mean(Z_high))**2)
    R_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Relative error
    rel_error = np.sqrt(ss_res) / np.sqrt(np.sum(np.abs(Z_high)**2)) * 100

    diagnostics = {
        'n_points_used': len(f_high),
        'freq_range': (f_high.min(), f_high.max()),
        'R_squared': R_squared,
        'rel_error': rel_error,
        'R_inf': R_inf,
        'L': L,
        'L_nH': L * 1e9,
        'R_k': R_k,
        'tau': tau_val,
        'tau_us': tau_val * 1e6,
        'f_char_kHz': 1 / (2 * np.pi * tau_val) / 1e3,
        'method': f'rlk_linear_{method_str}',
        'residual': residual,
        'warnings': warnings,
        'fit_success': True
    }

    if verbose:
        logger.info("")
        logger.info("Results:")
        logger.info(f"  R_inf = {R_inf:.4f} Ω")
        logger.info(f"  L     = {L*1e9:.2f} nH")
        logger.info(f"  R_k   = {R_k:.4f} Ω")
        logger.info(f"  τ     = {tau_val*1e6:.3f} μs")
        logger.info("")
        logger.info(f"Quality: R² = {R_squared:.4f}, rel_error = {rel_error:.2f}%")
        logger.info("=" * 60)

    return R_inf, L, R_k, tau_val, Z_fit, diagnostics


def estimate_rinf_with_inductance(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    max_L_nH: float = 1000.0,
    verbose: bool = True,
    plot: bool = False,
    save_plot: Optional[str] = None
) -> Tuple[float, float, None, dict, Optional[object]]:
    """
    Estimate R_inf using R-L-K() linear fit.

    Uses linear least squares on R_s - L - K(R, τ) model.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ω]
    max_L_nH : float
        Maximum reasonable inductance [nH] for warning (default: 1000)
    verbose : bool
        Log fitting progress (default: True)
    plot : bool
        Create diagnostic plot (default: False)
    save_plot : str, optional
        Path to save plot

    Returns
    -------
    R_inf : float
        High-frequency resistance [Ω]
    L : float
        Inductance [H]
    circuit : None
        Always None (legacy placeholder)
    diagnostics : dict
        Fitting diagnostics
    fig : matplotlib.Figure or None
        Diagnostic plot if plot=True

    Examples
    --------
    >>> R_inf, L, _, diag, fig = estimate_rinf_with_inductance(freq, Z, verbose=True)
    >>> print(f"R_inf = {R_inf:.3f} Ω, L = {L*1e9:.1f} nH")
    """
    if verbose:
        logger.info("")
        logger.info("#" * 70)
        logger.info("##  R-L-K() ANALYSIS: Linear R_inf estimation from HF data  ##")
        logger.info("#" * 70)

    # Validate input
    if len(frequencies) < 3:
        diagnostics = {
            'n_points_used': len(frequencies),
            'fit_success': False,
            'warnings': ['Insufficient data (< 3 points)'],
            'method': 'fallback_insufficient_data'
        }
        R_inf = np.median(Z.real) if len(Z) > 0 else 0.0
        logger.warning(f"Insufficient data, fallback: R_inf = {R_inf:.3f} Ω")
        return R_inf, 0.0, None, diagnostics, None

    # Perform fit
    try:
        R_inf, L, R_k, tau, Z_fit, diagnostics = fit_rlk_model(
            frequencies, Z, verbose=verbose
        )
    except Exception as e:
        # Fallback on error
        logger.warning(f"R-L-K() fit failed: {e}")
        diagnostics = {
            'fit_success': False,
            'warnings': [str(e)],
            'method': 'fallback_fit_error'
        }
        R_inf = np.median(Z.real)
        logger.warning(f"Fallback: R_inf = {R_inf:.3f} Ω (median)")
        return R_inf, 0.0, None, diagnostics, None

    # Validate inductance
    L_nH = L * 1e9
    if L_nH > max_L_nH:
        msg = f"High inductance L = {L_nH:.1f} nH > {max_L_nH:.0f} nH"
        diagnostics['warnings'].append(msg)
        logger.warning(msg)

    if L < 0:
        msg = f"Negative inductance L = {L_nH:.1f} nH (non-physical)"
        diagnostics['warnings'].append(msg)
        logger.warning(msg)

    # Plot if requested
    fig = None
    if plot:
        try:
            fig = _plot_rlk_fit(frequencies, Z, R_inf, L, R_k, tau,
                               diagnostics, save_plot)
        except Exception as e:
            logger.warning(f"Visualization failed: {e}")

    return R_inf, L, None, diagnostics, fig


def _plot_rlk_fit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    R_inf: float,
    L: float,
    R_k: float,
    tau: float,
    diagnostics: dict,
    save_path: Optional[str] = None
):
    """Create diagnostic plot for R-L-K() fit."""
    import matplotlib.pyplot as plt

    # Sort data
    frequencies, Z = sort_by_frequency(frequencies, Z)

    # Get highest decade mask (same selection as in fit)
    f_max = frequencies.max()
    f_min_decade = f_max / 10.0
    mask = frequencies >= f_min_decade

    # Extract highest decade data
    f_fit = frequencies[mask]
    Z_fit_data = Z[mask]

    # Compute fit ONLY for the fit region
    omega_fit = 2 * np.pi * f_fit
    Z_fit = R_inf + 1j * omega_fit * L + R_k / (1 + 1j * omega_fit * tau)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Nyquist plot - only highest decade data
    ax = axes[0, 0]
    ax.plot(Z_fit_data.real, -Z_fit_data.imag, 'o', color='#1f77b4',
            ms=6, alpha=0.8, label=f'VF data ({len(f_fit)} pts)',
            markeredgewidth=0.5, markeredgecolor='white')
    ax.plot(Z_fit.real, -Z_fit.imag, 'r-', lw=2, label='R-L-K() fit')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.axvline(R_inf, color='green', ls='--', lw=1.5, label=f'R_inf = {R_inf:.3f} Ω')
    ax.set_xlabel("Z' [Ω]")
    ax.set_ylabel("-Z'' [Ω]")
    ax.set_title('Nyquist Plot (highest decade)')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. Re(Z) vs frequency - only highest decade
    ax = axes[0, 1]
    ax.semilogx(f_fit, Z_fit_data.real, 'o', color='#1f77b4',
                ms=6, alpha=0.8, label='VF data',
                markeredgewidth=0.5, markeredgecolor='white')
    ax.semilogx(f_fit, Z_fit.real, 'r-', lw=2, label='Fit')
    ax.axhline(R_inf, color='green', ls='--', lw=1.5, label=f'R_inf = {R_inf:.3f} Ω')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel("Re(Z) [Ω]")
    ax.set_title('Re(Z) vs Frequency')
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3, which='both')

    # 3. Im(Z) vs frequency - only highest decade
    ax = axes[1, 0]
    ax.semilogx(f_fit, Z_fit_data.imag, 'o', color='#1f77b4',
                ms=6, alpha=0.8, label='VF data',
                markeredgewidth=0.5, markeredgecolor='white')
    ax.semilogx(f_fit, Z_fit.imag, 'r-', lw=2, label='Fit')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Im(Z) [Ω]')
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
        f"R_inf = {R_inf:.4f} Ω",
        f"L     = {L*1e9:.2f} nH",
        f"R_k   = {R_k:.4f} Ω",
        f"τ     = {tau*1e6:.3f} μs",
        "",
        "QUALITY METRICS",
        "-" * 20,
        f"R²        = {diagnostics.get('R_squared', 0):.4f}",
        f"Rel.error = {diagnostics.get('rel_error', 0):.2f}%",
        f"Points    = {diagnostics.get('n_points', 0)}",
        "",
        f"Method: {diagnostics.get('method', 'unknown')}"
    ]

    if diagnostics.get('warnings'):
        info_text.append("")
        info_text.append("WARNINGS:")
        for w in diagnostics['warnings']:
            info_text.append(f"  • {w}")

    ax.text(0.1, 0.95, '\n'.join(info_text), transform=ax.transAxes,
            fontsize=10, fontfamily='monospace', verticalalignment='top')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        logger.info(f"Plot saved: {save_path}")

    return fig


__all__ = ['fit_rlk_model', 'estimate_rinf_with_inductance']
