"""
Visualization functions for EIS data.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Sequence
from numpy.typing import NDArray

from ..fitting.circuit import FitResult

PLOT_GRID_ALPHA = 0.3
# Points of the fitted curve: the data are often sparse (~10/decade), and a
# curve drawn only through them shows straight segments instead of arcs.
FIT_CURVE_POINTS = 300


def plot_circuit_fit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    result: FitResult,
    title: Optional[str] = None
) -> plt.Figure:
    """
    Plot a circuit fit: Nyquist plot with the fitted curve, and residuals.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies the fit ran on [Hz]
    Z : ndarray of complex
        Impedance the fit ran on [Ohm]
    result : FitResult
        Output of `fit_equivalent_circuit` (or `best_result` of the
        multistart / differential evolution fit)
    title : str, optional
        Nyquist panel title. Default "Circuit fit: <circuit>"; pass a short
        one for a long circuit (a Voigt chain), whose expression would
        stretch the figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure with Nyquist plot and residuals
    """
    circuit = result.circuit
    params = list(result.params_opt)
    freq_curve = np.logspace(np.log10(frequencies.min()), np.log10(frequencies.max()),
                             FIT_CURVE_POINTS)
    Z_fit = circuit.impedance(freq_curve, params)
    Z_fit_at_data = circuit.impedance(frequencies, params)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Nyquist plot
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=5)
    ax1.plot(Z_fit.real, -Z_fit.imag, '-', label='Fit', linewidth=2)
    ax1.set_xlabel("Z' [Ω]")
    ax1.set_ylabel("-Z'' [Ω]")
    ax1.set_title(title if title is not None else f"Circuit fit: {circuit}")
    ax1.legend()
    ax1.grid(True, alpha=PLOT_GRID_ALPHA)
    ax1.set_aspect('equal', adjustable='datalim')

    # Residuals normalized by |Z| (same as KK validation)
    Z_mag_safe = np.maximum(np.abs(Z), 1e-15)
    res_real = (Z.real - Z_fit_at_data.real) / Z_mag_safe * 100  # in %
    res_imag = (Z.imag - Z_fit_at_data.imag) / Z_mag_safe * 100  # in %
    mean_res_real = np.mean(np.abs(res_real))
    mean_res_imag = np.mean(np.abs(res_imag))

    ax2.semilogx(frequencies, res_real, 'o', label='Real', markersize=4, color='#1f77b4')
    ax2.semilogx(frequencies, res_imag, 's', label='Imaginary', markersize=4, color='#ff7f0e')
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(y=5, color='r', linestyle=':', alpha=0.5, label='5% threshold')
    ax2.axhline(y=-5, color='r', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Residuals [%]")
    ax2.set_title(f"Fit residuals (Re: {mean_res_real:.2f}%, Im: {mean_res_imag:.2f}%)")
    ax2.legend(loc='best')
    ax2.grid(True, alpha=PLOT_GRID_ALPHA, which='both')

    plt.tight_layout()
    return fig


def visualize_data(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    title: str = "EIS Spectrum"
) -> plt.Figure:
    """
    Plot Nyquist and Bode diagrams.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Nyquist diagram
    ax1 = axes[0]
    ax1.plot(Z.real, -Z.imag, 'o-', markersize=4)
    ax1.set_xlabel("Z' [Ω]")
    ax1.set_ylabel("-Z'' [Ω]")
    ax1.set_title("Nyquist Diagram")
    ax1.grid(True, alpha=0.3)

    # Let matplotlib automatically determine axis range

    # Bode diagram - magnitude
    ax2 = axes[1]
    ax2.loglog(frequencies, np.abs(Z), 'o-', markersize=4)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("|Z| [Ω]")
    ax2.set_title("Bode Diagram - Magnitude")
    ax2.grid(True, alpha=0.3, which='both')

    # Bode diagram - phase
    ax3 = axes[2]
    phase = np.angle(Z, deg=True)
    ax3.semilogx(frequencies, phase, 'o-', markersize=4)
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Phase [deg]")
    ax3.set_title("Bode Diagram - Phase")
    ax3.grid(True, alpha=0.3, which='both')

    plt.suptitle(title)
    plt.tight_layout()
    return fig


def visualize_ocv(
    ocv_data: dict,
    title: str = "OCV Curve"
) -> Optional[plt.Figure]:
    """
    Visualize Open Circuit Voltage (OCV) curve.

    Parameters
    ----------
    ocv_data : dict
        Dictionary with keys 'time', 'Vf', 'Vm' (numpy arrays)
        - time: time in seconds [s]
        - Vf: filtered voltage [V]
        - Vm: measured/raw voltage [V]
    title : str
        Plot title

    Returns
    -------
    fig : Figure or None
        Matplotlib figure, or None if ocv_data is None
    """
    if ocv_data is None:
        return None

    time_s = ocv_data['time']
    time_min = time_s / 60.0  # Convert to minutes
    Vf = ocv_data['Vf']
    Vm = ocv_data['Vm']

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot both Vf and Vm
    ax.plot(time_min, Vm, 'b-', alpha=0.5, linewidth=1, label='Vm (raw)')
    ax.plot(time_min, Vf, 'r-', linewidth=1.5, label='Vf (filtered)')

    ax.set_xlabel('Time [min]')
    ax.set_ylabel('Voltage [V]')
    ax.set_title(f'{title} - Open Circuit Voltage Stabilization')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=PLOT_GRID_ALPHA)

    # Add statistics
    delta_V = abs(Vf[-1] - Vf[0]) * 1000  # mV
    total_time = time_min[-1]
    ax.text(0.02, 0.98, f'Duration: {total_time:.1f} min\n'
            f'Vf final: {Vf[-1]*1000:.1f} mV\n'
            f'Delta Vf: {delta_V:.2f} mV',
            transform=ax.transAxes, fontsize=9,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    return fig


def plot_rinf_fit(result) -> plt.Figure:
    """
    Plot the R_inf estimate over its fit window.

    Parameters
    ----------
    result : RinfResult
        Output of `estimate_rinf`

    Returns
    -------
    fig : Figure
        Nyquist, Re(Z) and Im(Z) of the window with the R-L-(R|Q) fit (if it
        ran), the fitted R_s and the fallback HF upper bound.
    """
    f, Z = result.f_window, result.Z_window
    Z_fit = None
    if result.fit is not None:
        f_dense = np.logspace(np.log10(f.min()), np.log10(f.max()), 200)
        Z_fit = result.fit.circuit.impedance(f_dense, list(result.fit.params_opt))

    lines = [(result.R_inf_hf, 'gray',
              f'HF bound Re(Z) = {result.R_inf_hf:.4g} Ohm at {result.f_hf:.3g} Hz')]
    if result.R_inf_fit is not None:
        lines.append((result.R_inf_fit, 'green',
                      f'fit R_s = {result.R_inf_fit:.4g} +- {result.R_inf_stderr:.2g} Ohm'))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    ax = axes[0]
    ax.plot(Z.real, -Z.imag, 'o', label=f'data ({len(f)} pts)')
    if Z_fit is not None:
        ax.plot(Z_fit.real, -Z_fit.imag, 'r-', label='R-L-(R|Q) fit')
    for value, color, label in lines:
        ax.axvline(value, color=color, ls='--', label=label)
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel("Z' [Ohm]")
    ax.set_ylabel("-Z'' [Ohm]")
    ax.legend(fontsize=8)

    for ax, part, ylabel in ((axes[1], np.real, 'Re(Z) [Ohm]'),
                             (axes[2], np.imag, 'Im(Z) [Ohm]')):
        ax.semilogx(f, part(Z), 'o', label='data')
        if Z_fit is not None:
            ax.semilogx(f_dense, part(Z_fit), 'r-', label='fit')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8)
    for value, color, _ in lines:
        axes[1].axhline(value, color=color, ls='--')
    axes[2].axhline(0, color='gray', lw=0.5)

    for ax in axes:
        ax.grid(True, alpha=PLOT_GRID_ALPHA, which='both')
    used = 'fit' if result.method == 'rlq_fit' else 'upper bound'
    fig.suptitle(f'R_inf = {result.R_inf:.4g} Ohm ({used})')
    plt.tight_layout()
    return fig


def plot_local_exponent(result) -> plt.Figure:
    """
    Plot the local CPE exponent n(f).

    Parameters
    ----------
    result : LocalExponentResult
        Output of `local_exponent`

    Returns
    -------
    fig : Figure
        n(f) with its uncertainty at the determined points; the rest grey.
    """
    f, n, u, valid = result.frequencies, result.n, result.n_uncertainty, result.valid
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.errorbar(f[valid], n[valid], yerr=u[valid], fmt='o', ms=4, capsize=2,
                label='determined')
    rest = ~valid & np.isfinite(n)
    if rest.any():
        ax.semilogx(f[rest], n[rest], 'o', ms=3, color='lightgray',
                    label='not determined')
    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('local exponent n')
    ax.grid(True, alpha=PLOT_GRID_ALPHA, which='both')
    ax.legend(fontsize=8)
    title = 'Local CPE exponent n(f) = d ln Re(Y) / d ln(omega)'
    if valid.any():
        title += f'\nspan {result.span:.2f} (one CPE: constant)'
    ax.set_title(title)
    plt.tight_layout()
    return fig


def _residual_panel(ax, frequencies, result, threshold: float, title: str,
                    flagged_frequencies: Sequence[float]) -> None:
    """Real/imag residuals [%] with +-threshold lines and a red band per flagged frequency."""
    ax.semilogx(frequencies, result.residuals_real * 100, 'o', label='Real', markersize=4)
    ax.semilogx(frequencies, result.residuals_imag * 100, 's', label='Imaginary', markersize=4)
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    for y in (threshold, -threshold):
        ax.axhline(y=y, color='r', linestyle=':', alpha=0.5)
    for i, f in enumerate(flagged_frequencies):
        ax.axvline(x=f, color='red', linestyle='-', alpha=0.25, linewidth=3, zorder=0,
                   label='Flagged point' if i == 0 else None)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Residuals [%]")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=PLOT_GRID_ALPHA)


def plot_kk_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    result,
    flagged_frequencies: Sequence[float] = ()
) -> plt.Figure:
    """
    Plot a Kramers-Kronig validation: Nyquist fit and residuals.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies the validation ran on [Hz]
    Z : ndarray of complex
        Impedance the validation ran on [Ohm]
    result : KKResult
        Successful output of `kramers_kronig_validation`
    flagged_frequencies : sequence of float, optional
        Frequencies to mark in the residual panel (e.g. outliers)

    Returns
    -------
    fig : Figure
        Nyquist plot with a dense fit curve (left) and residuals with the
        +-KK_RESIDUAL_THRESHOLD lines that `is_valid` counts against (right).
    """
    from ..validation.kramers_kronig import KK_RESIDUAL_THRESHOLD, reconstruct_impedance

    freq_plot = np.logspace(np.log10(frequencies.min()), np.log10(frequencies.max()), 300)
    # lin_kk_native(include_L=False) leaves no L at the end of elements
    Z_fit_plot = reconstruct_impedance(freq_plot, result.elements, result.tau, result.inductance,
                                       include_L=result.inductance is not None,
                                       C_value=result.capacitance)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax1 = axes[0]
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=4)
    ax1.plot(Z_fit_plot.real, -Z_fit_plot.imag, '-', label='KK fit', linewidth=2)
    ax1.set_xlabel("Z' [Ohm]")
    ax1.set_ylabel("-Z'' [Ohm]")
    ax1.set_title(f"Kramers-Kronig fit in real domain (M={result.M})")
    ax1.legend()
    ax1.grid(True, alpha=PLOT_GRID_ALPHA)
    ax1.set_aspect('equal', adjustable='datalim')

    _residual_panel(axes[1], frequencies, result, KK_RESIDUAL_THRESHOLD,
                    f"KK residuals (stop mu={result.mu:.3f}, chi^2={result.pseudo_chisqr:.2e}, "
                    f"noise~{result.noise_estimate:.1f}%)", flagged_frequencies)

    plt.tight_layout()
    return fig


def plot_zhit_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    result,
    flagged_frequencies: Sequence[float] = ()
) -> plt.Figure:
    """
    Plot a Z-HIT validation: |Z| reconstruction and residuals.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies the validation ran on [Hz]
    Z : ndarray of complex
        Impedance the validation ran on [Ohm]
    result : ZHITResult
        Successful output of `zhit_validation`
    flagged_frequencies : sequence of float, optional
        Frequencies to mark in the residual panel (e.g. outliers)

    Returns
    -------
    fig : Figure
        Measured vs reconstructed |Z| (left) and complex residuals with
        +-result.quality_threshold lines (right).
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    ax1 = axes[0]
    ax1.loglog(frequencies, np.abs(Z), 'o', label='Measured', markersize=4)
    # The result keeps the caller's point order; a line needs sorted frequencies
    order = np.argsort(frequencies)
    ax1.loglog(frequencies[order], result.Z_mag_reconstructed[order], '-',
               label='Z-HIT reconstruction', linewidth=2, color='red')
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("|Z| [Ohm]")
    ax1.set_title("Z-HIT validation")
    ax1.legend()
    ax1.grid(True, alpha=PLOT_GRID_ALPHA, which='both')

    _residual_panel(axes[1], frequencies, result, result.quality_threshold,
                    f"Z-HIT residuals (χ²={result.pseudo_chisqr:.2e}, noise≤{result.noise_estimate:.1f}%)",
                    flagged_frequencies)

    plt.tight_layout()
    return fig


def plot_thd(
    frequencies: NDArray[np.float64],
    current_thd: Optional[NDArray[np.float64]],
    voltage_thd: Optional[NDArray[np.float64]],
    threshold: float
) -> plt.Figure:
    """
    Plot the per-point THD the instrument recorded against frequency.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz]
    current_thd, voltage_thd : ndarray of float or None
        THD per point as a fraction (see `LoadResult`); a None channel, or
        one without any value > 0, is left out
    threshold : float
        Linearity threshold drawn as a line (fraction)

    Returns
    -------
    fig : Figure
        THD [%] of each channel on a log-log axis with the threshold line.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    frequencies = np.asarray(frequencies, dtype=np.float64)
    # The data keep the file's point order; a line needs sorted frequencies
    order = np.argsort(frequencies)
    for thd, label, marker in ((current_thd, 'Current', 'o'), (voltage_thd, 'Voltage', 's')):
        if thd is None:
            continue
        # Missing values (NaN, or <= 0 as thd_check treats them) leave a gap;
        # a channel without any value gets no legend entry
        thd = np.asarray(thd, dtype=np.float64)
        percent = np.where(np.isfinite(thd) & (thd > 0), thd * 100, np.nan)[order]
        if np.isnan(percent).all():
            continue
        ax.loglog(frequencies[order], percent, marker=marker, linestyle='-',
                  linewidth=1, markersize=4, label=label)
    ax.axhline(y=threshold * 100, color='r', linestyle=':', alpha=0.5,
               label=f'Threshold {threshold * 100:.2g} %')
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("THD [%]")
    ax.set_title("Total harmonic distortion (instrument)")
    ax.legend()
    ax.grid(True, alpha=PLOT_GRID_ALPHA, which='both')
    plt.tight_layout()
    return fig
