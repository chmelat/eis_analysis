"""
Diagnostic visualization for R+L fitting.
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from typing import Optional
from numpy.typing import NDArray

logger = logging.getLogger(__name__)


def plot_rl_fit_diagnostics(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    R_inf: float,
    L: float,
    f_fit: NDArray[np.float64],
    Z_fit: NDArray[np.complex128],
    diagnostics: dict,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Create advanced R+L fit visualization with diagnostics.

    Parameters
    ----------
    frequencies : NDArray[np.float64]
        All frequencies [Hz]
    Z : NDArray[np.complex128]
        All impedance data [Ω]
    R_inf : float
        Fitted R_inf [Ω]
    L : float
        Fitted inductance [H]
    f_fit : NDArray[np.float64]
        Frequencies used for fit [Hz]
    Z_fit : NDArray[np.complex128]
        Fitted impedance [Ω]
    diagnostics : dict
        Diagnostic information from fit
    save_path : Optional[str]
        Path to save plot (None = display only)

    Returns
    -------
    fig : plt.Figure
        Matplotlib figure
    """

    # Create figure with 2x2 subplots
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # ===========================================================================
    # Subplot 1: Nyquist plot with R_inf (HF points only)
    # ===========================================================================
    ax1 = fig.add_subplot(gs[0, 0])

    # Mask for HF points used for fit
    mask_fit = np.isin(frequencies, f_fit)
    Z_fit_data = Z[mask_fit]

    # Only HF points used for fit
    ax1.plot(Z_fit_data.real, -Z_fit_data.imag, 'o', color='#1f77b4',
             label=f'HF points for fit ({len(f_fit)})',
             markersize=6, alpha=0.8, markeredgewidth=0.5, markeredgecolor='white')

    # R_inf vertical line
    ax1.axvline(R_inf, color='red', linestyle='--', linewidth=2,
                label=f'R_inf = {R_inf:.2f} Ω', alpha=0.8)

    ax1.set_xlabel("Z' [Ω]", fontsize=11, fontweight='bold')
    ax1.set_ylabel("-Z'' [Ω]", fontsize=11, fontweight='bold')
    ax1.set_title("Nyquist Plot", fontsize=12, fontweight='bold')
    ax1.legend(loc='best', framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Let matplotlib automatically determine axis range

    # ===========================================================================
    # Subplot 2: Re(Z) vs frequency (HF points only)
    # ===========================================================================
    ax2 = fig.add_subplot(gs[0, 1])

    # Only HF points
    ax2.semilogx(f_fit, Z_fit_data.real, 'o', color='#1f77b4',
                 label='HF points for fit', markersize=6, alpha=0.8,
                 markeredgewidth=0.5, markeredgecolor='white')

    # R_inf horizontal line
    ax2.axhline(R_inf, color='red', linestyle='--', linewidth=2,
                label=f'R_inf = {R_inf:.2f} Ω', alpha=0.8)

    # Fit line (if exists)
    if Z_fit is not None and len(Z_fit) > 0:
        ax2.semilogx(f_fit, Z_fit.real, '-', color='green',
                     linewidth=2.5, label='R+L fit', alpha=0.7)

    ax2.set_xlabel("Frequency [Hz]", fontsize=11, fontweight='bold')
    ax2.set_ylabel("Re(Z) [Ω]", fontsize=11, fontweight='bold')
    ax2.set_title("Re(Z) vs Frequency", fontsize=12, fontweight='bold')
    ax2.legend(loc='best', framealpha=0.9)
    ax2.grid(True, alpha=0.3, which='both', linestyle='--')

    # ===========================================================================
    # Subplot 3: Im(Z) vs frequency (HF points only)
    # ===========================================================================
    ax3 = fig.add_subplot(gs[1, 0])

    # Only HF points
    ax3.semilogx(f_fit, Z_fit_data.imag, 'o', color='#1f77b4',
                 label='HF points for fit', markersize=6, alpha=0.8,
                 markeredgewidth=0.5, markeredgecolor='white')

    # Fit line
    if Z_fit is not None and len(Z_fit) > 0:
        ax3.semilogx(f_fit, Z_fit.imag, '-', color='green',
                     linewidth=2.5, label=f'L fit: {diagnostics["L_nH"]:.1f} nH',
                     alpha=0.7)

    # Zero line
    ax3.axhline(0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    ax3.set_xlabel("Frequency [Hz]", fontsize=11, fontweight='bold')
    ax3.set_ylabel("Im(Z) [Ω]", fontsize=11, fontweight='bold')
    ax3.set_title("Im(Z) vs Frequency (inductance)", fontsize=12, fontweight='bold')
    ax3.legend(loc='best', framealpha=0.9)
    ax3.grid(True, alpha=0.3, which='both', linestyle='--')

    # ===========================================================================
    # Subplot 4: Diagnostics and residuals
    # ===========================================================================
    ax4 = fig.add_subplot(gs[1, 1])

    # Residuals (if fit exists)
    if Z_fit is not None and len(Z_fit) > 0:
        residuals_real = Z_fit_data.real - Z_fit.real
        residuals_imag = Z_fit_data.imag - Z_fit.imag
        residuals_mag = np.abs(Z_fit_data - Z_fit)

        ax4.semilogx(f_fit, residuals_real, 'o-', color='#d62728',
                     label="Res. Re(Z)", markersize=5, linewidth=1.5, alpha=0.7)
        ax4.semilogx(f_fit, residuals_imag, 's-', color='#9467bd',
                     label="Res. Im(Z)", markersize=5, linewidth=1.5, alpha=0.7)
        ax4.semilogx(f_fit, residuals_mag, '^-', color='black',
                     label="Res. |Z|", markersize=5, linewidth=1.5, alpha=0.7)

        ax4.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax4.set_xlabel("Frequency [Hz]", fontsize=11, fontweight='bold')
        ax4.set_ylabel("Residuals [Ω]", fontsize=11, fontweight='bold')
        ax4.set_title("Residuals (Data - Fit)", fontsize=12, fontweight='bold')
        ax4.legend(loc='best', framealpha=0.9)
        ax4.grid(True, alpha=0.3, which='both', linestyle='--')
    else:
        # If no fit, display diagnostics as text
        ax4.axis('off')

    # Text diagnostics in bottom right corner
    # Get new diagnostic parameters (if they exist)
    re_var = diagnostics.get('re_variation_pct', None)
    im_r2 = diagnostics.get('im_linearity_r2', None)

    diag_text = f"""
+===================================+
|   R_INF ESTIMATION DIAGNOSTICS    |
+===================================+
| R_inf:      {R_inf:8.3f} Ω         |
| L:          {diagnostics['L_nH']:8.2f} nH        |
| --------------------------------- |
| Fit quality:                      |
|   R2:       {diagnostics['R_squared']:8.6f}       |
|   RMSE:     {diagnostics['rmse']:8.4f} Ω         |
|   Rel err:  {diagnostics['rel_error']:8.2f} %         |
| --------------------------------- |
| Data quality:                     |"""

    if re_var is not None:
        re_status = '[OK]' if re_var < 10.0 else '[!] '
        diag_text += f"\n| {re_status} Re(Z) var: {re_var:6.1f}% (<10% OK)  |"

    if im_r2 is not None:
        im_status = '[OK]' if im_r2 > 0.8 else '[!] '
        diag_text += f"\n| {im_status} Im(Z) R²:  {im_r2:6.3f}  (>0.8 OK)  |"

    diag_text += f"""
| --------------------------------- |
| VF body:    {diagnostics['n_points_used']:8d}           |
| f_min:      {diagnostics['freq_range'][0]:8.1f} Hz       |
| f_max:      {diagnostics['freq_range'][1]:8.1f} Hz       |
| --------------------------------- |
| Metoda:     {diagnostics['method'][:23]:23s} |
| Fit OK:     {'[OK] Yes' if diagnostics['fit_success'] else '[!]  No':23s} |
+===================================+
"""

    if diagnostics['warnings']:
        diag_text += "\n[!] WARNINGS:\n"
        for i, warning in enumerate(diagnostics['warnings'][:3], 1):
            diag_text += f"  {i}. {warning[:50]}\n"

    ax4.text(0.05, 0.95, diag_text, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    # Main title - dynamic based on method
    behavior = diagnostics.get('behavior', 'unknown')
    if behavior == 'capacitive':
        method_name = 'Extrapolation'
    elif behavior == 'inductive':
        method_name = 'R-L fit'
    else:
        method_name = 'Estimation'

    fig.suptitle(f'R_inf Diagnostics ({method_name}): R_inf = {R_inf:.3f} Ω, L = {diagnostics["L_nH"]:.2f} nH',
                 fontsize=14, fontweight='bold', y=0.995)

    # Save or display
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        logger.info(f"Visualization saved: {save_path}")

    return fig
