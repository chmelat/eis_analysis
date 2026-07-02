"""
Visualization for DRT analysis.

Builds the DRT result figure: gamma(tau) spectrum, Nyquist reconstruction check,
and (for the GMM method) per-peak deconvolution and BIC model-selection panels.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List, Dict
from numpy.typing import NDArray
from scipy.signal import find_peaks

from ..fitting.config import DRT_PEAK_HEIGHT_THRESHOLD


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
