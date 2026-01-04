"""
Visualization functions for EIS data.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional
from numpy.typing import NDArray

PLOT_GRID_ALPHA = 0.3


def plot_circuit_fit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128],
    circuit: Optional[object] = None,
    title: Optional[str] = None,
    figsize: tuple = (12, 5),
    Z_fit_at_data: Optional[NDArray[np.complex128]] = None
) -> plt.Figure:
    """
    Create Nyquist plot with measured data, circuit fit, and residuals.

    Parameters
    ----------
    frequencies : ndarray of float
        Measurement frequencies [Hz]
    Z : ndarray of complex
        Measured impedance data [Ω]
    Z_fit : ndarray of complex
        Fitted impedance for plotting (can be denser than data) [Ω]
    circuit : object, optional
        Circuit object (for title and residuals calculation)
    title : str, optional
        Custom plot title
    figsize : tuple, optional
        Figure size (default: (12, 5))
    Z_fit_at_data : ndarray of complex, optional
        Fitted impedance at original data frequencies [Ω].
        If None and circuit is provided, will be calculated from circuit.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure with Nyquist plot and residuals
    """
    # Calculate Z_fit at data frequencies for residuals
    if Z_fit_at_data is None and circuit is not None:
        if hasattr(circuit, 'impedance') and hasattr(circuit, 'get_all_params'):
            params = circuit.get_all_params()
            Z_fit_at_data = circuit.impedance(frequencies, params)

    # Create figure with 2 subplots if we have data for residuals
    if Z_fit_at_data is not None:
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        ax1 = axes[0]
        ax2 = axes[1]
    else:
        fig, ax1 = plt.subplots(figsize=(8, 6))
        ax2 = None

    # Nyquist plot
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=5)
    ax1.plot(Z_fit.real, -Z_fit.imag, '-', label='Fit', linewidth=2)
    ax1.set_xlabel("Z' [Ω]")
    ax1.set_ylabel("-Z'' [Ω]")

    if title:
        ax1.set_title(title)
    elif circuit is not None:
        ax1.set_title(f"Circuit fit: {circuit}")
    else:
        ax1.set_title("Circuit fit")

    ax1.legend()
    ax1.grid(True, alpha=PLOT_GRID_ALPHA)
    ax1.set_aspect('equal', adjustable='datalim')

    # Residuals plot
    if ax2 is not None and Z_fit_at_data is not None:
        # Calculate residuals normalized by |Z| (same as KK validation)
        Z_mag = np.abs(Z)
        Z_mag_safe = np.maximum(Z_mag, 1e-15)

        res_real = (Z.real - Z_fit_at_data.real) / Z_mag_safe * 100  # in %
        res_imag = (Z.imag - Z_fit_at_data.imag) / Z_mag_safe * 100  # in %

        # Calculate mean residuals for title
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
    Z: NDArray[np.complex128],
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
