"""
Visualization for DRT analysis.

Builds the DRT result figure: gamma(tau) spectrum, Nyquist reconstruction check,
and (for the GMM method) per-peak deconvolution and BIC model-selection panels.
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray
from scipy.signal import find_peaks

from ..drt.results import DRTResult
from ..fitting.config import DRT_PEAK_HEIGHT_THRESHOLD, GMM_N_COMPONENTS_RANGE


def plot_drt(Z: NDArray[np.complex128], result: DRTResult) -> plt.Figure:
    """
    Plot a DRT result: gamma(tau) spectrum and Nyquist reconstruction check.

    Parameters
    ----------
    Z : ndarray of complex
        Impedance the DRT was computed from [Ohm]
    result : DRTResult
        Successful output of `calculate_drt`

    Returns
    -------
    fig : Figure
        DRT spectrum (normalized when the result is) with the lambda-probe
        curves as thin overlays when they were computed, the tau grid past the
        slow end of the measured window shaded as extrapolated, and data vs.
        reconstruction in the Nyquist plane. For the GMM method with peaks,
        two more panels: per-peak Gaussians and the BIC model selection.
    """
    diag = result.diagnostics
    if (not result.success or diag is None or result.Z_reconstructed is None
            or diag.tau_window is None):
        raise ValueError("plot_drt needs a successful DRTResult")
    assert result.tau is not None and result.gamma is not None  # guaranteed by success

    tau, gamma, gamma_original = result.tau, result.gamma, result.gamma_original
    Z_reconstructed = result.Z_reconstructed
    lambda_reg = result.lambda_used
    normalize_rpol, peak_method = diag.normalized, diag.peak_method
    peaks_result, bic_scores = result.peaks, result.bic_scores
    tau_window = diag.tau_window
    # Lambda-probe curves are physical [Ohm]; match the displayed gamma
    probe_curves = None
    if diag.stability is not None:
        scale = diag.R_pol_from_gamma if normalize_rpol else 1.0
        probe_curves = [
            (p.lambda_value, p.gamma / scale)
            for p in diag.stability.probe_points
            if p.success and p.gamma is not None
        ]

    use_gmm = (peak_method == 'gmm' and
               peaks_result is not None and len(peaks_result) > 0)

    if use_gmm:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        ax1, ax2 = axes[0, 0], axes[0, 1]
        ax3, ax4 = axes[1, 0], axes[1, 1]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        ax1, ax2 = axes[0], axes[1]

    # === DRT Spectrum ===
    if probe_curves:
        for probe_lambda, probe_gamma in probe_curves:
            ax1.semilogx(tau, probe_gamma, '-', linewidth=1, alpha=0.35,
                         label=f'lambda = {probe_lambda:.1e}')
    ax1.semilogx(tau, gamma, 'b-', linewidth=2, label='DRT gamma(tau)')
    ax1.fill_between(tau, 0, gamma, alpha=0.3)
    if tau[-1] > tau_window[1] * (1 + 1e-9):
        ax1.axvspan(tau_window[1], tau[-1], color='gray', alpha=0.15,
                    label='past measured window (extrapolated)')
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

    if probe_curves:
        ax1.legend(fontsize=8)

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
        assert peaks_result is not None  # guaranteed by use_gmm
        log_tau = np.log10(tau)
        # GMM components below are scaled by R_estimate [Ohm], so this panel
        # must plot the unnormalized gamma even when ax1 shows gamma/R_pol.
        gamma_ohm = gamma_original if gamma_original is not None else gamma
        ax3.semilogx(tau, gamma_ohm, 'b-', linewidth=2, label='DRT gamma(tau)', alpha=0.7)

        colors = plt.get_cmap('tab10')(np.linspace(0, 1, len(peaks_result)))
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
            lo = GMM_N_COMPONENTS_RANGE[0]
            n_range = range(lo, lo + len(bic_scores))
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
