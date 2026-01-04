#!/usr/bin/env python3
"""
Generátor syntetických EIS dat s Warburgovým chováním
======================================================

Simuluje elektrochemický systém vykazující difúzní impedanci (Warburg element)
a exportuje data v Gamry DTA formátu.

Typické použití:
    python3 generate_warburg_data.py output.DTA
    python3 generate_warburg_data.py output.DTA --noise 2.0 --ppd 10
    python3 generate_warburg_data.py output.DTA --circuit randles --plot

Autor: EIS Analysis Toolkit v2.0
Datum: 2025-12-13
"""

import argparse
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from typing import Optional


def warburg_impedance(freq: np.ndarray, sigma: float) -> np.ndarray:
    """
    Výpočet semi-infinite Warburg impedance.

    Z_W = σ/√ω - jσ/√ω = σ(1-j)/√ω

    Parameters
    ----------
    freq : np.ndarray
        Frekvence [Hz]
    sigma : float
        Warburg koeficient [Ω·s^(-1/2)]

    Returns
    -------
    np.ndarray
        Komplexní impedance [Ω]
    """
    omega = 2 * np.pi * freq
    Z_W = sigma / np.sqrt(omega) * (1 - 1j)
    return Z_W


def finite_warburg_impedance(freq: np.ndarray, R_W: float, tau_W: float) -> np.ndarray:
    """
    Výpočet finite-length Warburg impedance (bounded diffusion).

    Z_W = R_W * tanh(√(jωτ_W)) / √(jωτ_W)

    Parameters
    ----------
    freq : np.ndarray
        Frekvence [Hz]
    R_W : float
        Warburg odpor [Ω]
    tau_W : float
        Časová konstanta difúze [s]

    Returns
    -------
    np.ndarray
        Komplexní impedance [Ω]
    """
    omega = 2 * np.pi * freq
    arg = np.sqrt(1j * omega * tau_W)
    Z_W = R_W * np.tanh(arg) / arg
    return Z_W


def randles_circuit(
    freq: np.ndarray,
    R_s: float = 10.0,
    R_ct: float = 100.0,
    C_dl: float = 20e-6,
    sigma: float = 50.0,
    warburg_type: str = 'semi-infinite'
) -> np.ndarray:
    """
    Randles ekvivalentní obvod: R_s-p(R_ct,C_dl)-W

    Typický pro elektrochemické systémy s difúzí.

    Parameters
    ----------
    freq : np.ndarray
        Frekvence [Hz]
    R_s : float
        Solution resistance [Ω]
    R_ct : float
        Charge transfer resistance [Ω]
    C_dl : float
        Double layer capacitance [F]
    sigma : float
        Warburg koeficient [Ω·s^(-1/2)] pro semi-infinite
        nebo R_W [Ω] pro finite
    warburg_type : str
        'semi-infinite' nebo 'finite'

    Returns
    -------
    np.ndarray
        Komplexní impedance [Ω]
    """
    omega = 2 * np.pi * freq

    # Kapacitní impedance
    Z_C = 1 / (1j * omega * C_dl)

    # R_ct || C_dl
    Z_parallel = 1 / (1/R_ct + 1/Z_C)

    # Warburg
    if warburg_type == 'semi-infinite':
        Z_W = warburg_impedance(freq, sigma)
    else:  # finite
        tau_W = 1.0  # defaultní časová konstanta [s]
        Z_W = finite_warburg_impedance(freq, sigma, tau_W)

    # Celková impedance
    Z = R_s + Z_parallel + Z_W

    return Z


def complex_warburg_circuit(
    freq: np.ndarray,
    R_s: float = 15.0,
    R_ct1: float = 50.0,
    C_dl1: float = 10e-6,
    R_ct2: float = 150.0,
    C_dl2: float = 50e-6,
    sigma: float = 30.0
) -> np.ndarray:
    """
    Složitější obvod s dvěma časovými konstantami a Warburgem.

    R_s-p(R_ct1,C_dl1)-p(R_ct2,C_dl2)-W

    Typický pro systémy s více elektrodovými procesy a difúzí.

    Parameters
    ----------
    freq : np.ndarray
        Frekvence [Hz]
    R_s : float
        Solution resistance [Ω]
    R_ct1 : float
        První charge transfer resistance [Ω]
    C_dl1 : float
        První double layer capacitance [F]
    R_ct2 : float
        Druhý charge transfer resistance [Ω]
    C_dl2 : float
        Druhá double layer capacitance [F]
    sigma : float
        Warburg koeficient [Ω·s^(-1/2)]

    Returns
    -------
    np.ndarray
        Komplexní impedance [Ω]
    """
    omega = 2 * np.pi * freq

    # První R||C
    Z_C1 = 1 / (1j * omega * C_dl1)
    Z_par1 = 1 / (1/R_ct1 + 1/Z_C1)

    # Druhý R||C
    Z_C2 = 1 / (1j * omega * C_dl2)
    Z_par2 = 1 / (1/R_ct2 + 1/Z_C2)

    # Warburg
    Z_W = warburg_impedance(freq, sigma)

    # Celková impedance
    Z = R_s + Z_par1 + Z_par2 + Z_W

    return Z


def add_noise(Z: np.ndarray, noise_percent: float, seed: Optional[int] = None) -> np.ndarray:
    """
    Přidá gaussovský šum k impedančním datům.

    Parameters
    ----------
    Z : np.ndarray
        Komplexní impedance [Ω]
    noise_percent : float
        Úroveň šumu [%] relativně k |Z|
    seed : int, optional
        Random seed pro reprodukovatelnost

    Returns
    -------
    np.ndarray
        Zašuměná impedance [Ω]
    """
    if seed is not None:
        np.random.seed(seed)

    if noise_percent <= 0:
        return Z

    # Šum proporcionální k velikosti impedance
    noise_std = noise_percent / 100.0
    Z_mag = np.abs(Z)

    # Gaussovský šum na reálnou a imaginární část
    noise_real = np.random.normal(0, noise_std * Z_mag)
    noise_imag = np.random.normal(0, noise_std * Z_mag)

    Z_noisy = Z + noise_real + 1j * noise_imag

    return Z_noisy


def generate_frequencies(
    f_min: float = 0.01,
    f_max: float = 100000.0,
    points_per_decade: int = 10
) -> np.ndarray:
    """
    Generuje logaritmicky rozmístěné frekvence.

    Parameters
    ----------
    f_min : float
        Minimální frekvence [Hz]
    f_max : float
        Maximální frekvence [Hz]
    points_per_decade : int
        Počet bodů na dekádu

    Returns
    -------
    np.ndarray
        Frekvence [Hz], sestupně řazené (Gamry konvence)
    """
    n_decades = np.log10(f_max / f_min)
    n_points = int(n_decades * points_per_decade) + 1

    frequencies = np.logspace(np.log10(f_min), np.log10(f_max), n_points)

    # Gamry měří od vysokých k nízkým frekvencím
    frequencies = frequencies[::-1]

    return frequencies


def export_csv(
    filename: str,
    freq: np.ndarray,
    Z: np.ndarray
) -> None:
    """
    Exportuje data do CSV formátu.

    Parameters
    ----------
    filename : str
        Výstupní CSV soubor (např. 'data.csv')
    freq : np.ndarray
        Frekvence [Hz]
    Z : np.ndarray
        Komplexní impedance [Ω]
    """
    with open(filename, 'w') as f:
        f.write("frequency,Z_real,Z_imag\n")
        for i in range(len(freq)):
            f.write(f"{freq[i]:.6e},{Z[i].real:.6e},{Z[i].imag:.6e}\n")

    print(f"✓ CSV data exportována do: {filename}")


def export_gamry_dta(
    filename: str,
    freq: np.ndarray,
    Z: np.ndarray,
    circuit_description: str = "Randles circuit with Warburg",
    notes: str = ""
) -> None:
    """
    Exportuje data do Gamry DTA formátu.

    Parameters
    ----------
    filename : str
        Výstupní soubor (např. 'data.DTA')
    freq : np.ndarray
        Frekvence [Hz]
    Z : np.ndarray
        Komplexní impedance [Ω]
    circuit_description : str
        Popis obvodu
    notes : str
        Poznámky
    """
    Z_real = Z.real
    Z_imag = Z.imag
    Z_mag = np.abs(Z)
    Z_phase = np.angle(Z, deg=True)

    with open(filename, 'w') as f:
        # Header
        f.write("EXPLAIN\n")
        f.write("TAG\tSYNTHETIC\n")
        f.write("TITLE\tLABEL\tSimulated Warburg Data\tSynthetic\n")
        f.write(f"DATE\tLABEL\t{datetime.now().strftime('%m/%d/%Y')}\tDate\n")
        f.write(f"TIME\tLABEL\t{datetime.now().strftime('%H:%M:%S')}\tTime\n")
        f.write(f"NOTES\tLABEL\t{circuit_description}\tNotes\n")
        if notes:
            f.write(f"NOTES\tLABEL\t{notes}\tAdditional Notes\n")
        f.write("AREA\tQUANT\t1.0\tcm^2\tArea\n")
        f.write("PSTAT\tLABEL\tSimulated\tPotentiostat\n")
        f.write("TAG\tSYNTHETIC-END\n")
        f.write("ZCURVE\tTABLE\n")  # Changed from CURVE to ZCURVE for compatibility
        f.write("Pt\t#\tFreq\tZreal\tZimag\tZmod\tZphz\n")
        f.write("\t\t Hz\tOhm\tOhm\tOhm\t°\n")

        # Data - format: Pt, #, Freq, Zreal, Zimag, Zmod, Zphz
        for i in range(len(freq)):
            f.write(f"{i+1}\t{i+1}\t{freq[i]:.5e}\t{Z_real[i]:.6e}\t{Z_imag[i]:.6e}\t{Z_mag[i]:.6e}\t{Z_phase[i]:.4f}\n")

    print(f"✓ DTA data exportována do: {filename}")
    print(f"  Počet bodů: {len(freq)}")
    print(f"  Frekvenční rozsah: {freq[-1]:.2e} - {freq[0]:.2e} Hz")
    print(f"  |Z| rozsah: {Z_mag.min():.2f} - {Z_mag.max():.2f} Ω")


def plot_data(freq: np.ndarray, Z: np.ndarray, title: str = "Warburg Circuit Simulation") -> None:
    """
    Zobrazí Nyquist a Bode diagramy.

    Parameters
    ----------
    freq : np.ndarray
        Frekvence [Hz]
    Z : np.ndarray
        Komplexní impedance [Ω]
    title : str
        Titulek grafu
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Nyquist
    axes[0].plot(Z.real, -Z.imag, 'o-', markersize=4, linewidth=1)
    axes[0].set_xlabel("Z' [Ω]", fontsize=11)
    axes[0].set_ylabel("-Z'' [Ω]", fontsize=11)
    axes[0].set_title("Nyquist Plot", fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    axes[0].axis('equal')

    # Bode - magnitude
    axes[1].loglog(freq, np.abs(Z), 'o-', markersize=4, linewidth=1, color='C1')
    axes[1].set_xlabel("Frequency [Hz]", fontsize=11)
    axes[1].set_ylabel("|Z| [Ω]", fontsize=11)
    axes[1].set_title("Bode Plot - Magnitude", fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, which='both')

    # Bode - phase
    axes[2].semilogx(freq, np.angle(Z, deg=True), 'o-', markersize=4, linewidth=1, color='C2')
    axes[2].set_xlabel("Frequency [Hz]", fontsize=11)
    axes[2].set_ylabel("Phase [°]", fontsize=11)
    axes[2].set_title("Bode Plot - Phase", fontsize=12, fontweight='bold')
    axes[2].grid(True, alpha=0.3, which='both')
    axes[2].axhline(-45, color='gray', linestyle='--', alpha=0.5, label='Warburg (-45°)')
    axes[2].legend()

    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.show()


def main():
    """Hlavní funkce s CLI."""
    parser = argparse.ArgumentParser(
        description='Generátor syntetických EIS dat s Warburgovým chováním',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  # Základní použití (Randles obvod, 10 bodů/dekáda, 1% šum)
  python3 generate_warburg_data.py output.DTA

  # Vlastní šum a rozlišení
  python3 generate_warburg_data.py output.DTA --noise 2.5 --ppd 15

  # Složitější obvod
  python3 generate_warburg_data.py output.DTA --circuit complex --plot

  # Finite Warburg
  python3 generate_warburg_data.py output.DTA --warburg-type finite

  # Vlastní frekvenční rozsah
  python3 generate_warburg_data.py output.DTA --f-min 0.001 --f-max 1e6 --ppd 12

  # Bez šumu s vizualizací
  python3 generate_warburg_data.py output.DTA --noise 0 --plot

Obvody:
  randles  - R_s-p(R_ct,C_dl)-W (jednoduchý, typický pro elektrochemii)
  complex  - R_s-p(R_ct1,C_dl1)-p(R_ct2,C_dl2)-W (dva procesy + difúze)

Warburg typy:
  semi-infinite - Neomezená difúze (typický pro tlustou vrstvu elektrolytu)
  finite        - Omezená difúze (bounded, tenká vrstva)
        """
    )

    parser.add_argument('output', help='Výstupní soubor (např. warburg_data.DTA nebo warburg_data.csv)')

    parser.add_argument('--format', '-f', choices=['dta', 'csv', 'both'], default='dta',
                        help='Výstupní formát: dta (Gamry), csv, nebo both (default: dta)')

    parser.add_argument('--circuit', '-c', choices=['randles', 'complex'], default='randles',
                        help='Typ obvodu (default: randles)')

    parser.add_argument('--warburg-type', '-w', choices=['semi-infinite', 'finite'],
                        default='semi-infinite',
                        help='Typ Warburg elementu (default: semi-infinite)')

    parser.add_argument('--noise', '-n', type=float, default=1.0,
                        help='Úroveň šumu v procentech (default: 1.0)')

    parser.add_argument('--ppd', type=int, default=10,
                        help='Počet bodů na dekádu (points per decade, default: 10)')

    parser.add_argument('--f-min', type=float, default=0.01,
                        help='Minimální frekvence [Hz] (default: 0.01)')

    parser.add_argument('--f-max', type=float, default=1e5,
                        help='Maximální frekvence [Hz] (default: 100000)')

    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed pro reprodukovatelnost')

    parser.add_argument('--plot', '-p', action='store_true',
                        help='Zobrazit grafy (Nyquist, Bode)')

    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Verbose výstup')

    args = parser.parse_args()

    # Generování frekvencí
    if args.verbose:
        print(f"Generování frekvencí: {args.f_min} - {args.f_max} Hz, {args.ppd} bodů/dekáda")

    freq = generate_frequencies(args.f_min, args.f_max, args.ppd)

    # Výpočet impedance podle zvoleného obvodu
    if args.verbose:
        print(f"Výpočet impedance pro obvod: {args.circuit}")

    if args.circuit == 'randles':
        circuit_desc = f"Randles circuit: R_s-p(R_ct,C_dl)-W ({args.warburg_type})"
        Z = randles_circuit(
            freq,
            R_s=10.0,
            R_ct=100.0,
            C_dl=20e-6,
            sigma=50.0,
            warburg_type=args.warburg_type
        )
    else:  # complex
        circuit_desc = f"Complex circuit: R_s-p(R_ct1,C_dl1)-p(R_ct2,C_dl2)-W ({args.warburg_type})"
        Z = complex_warburg_circuit(
            freq,
            R_s=15.0,
            R_ct1=50.0,
            C_dl1=10e-6,
            R_ct2=150.0,
            C_dl2=50e-6,
            sigma=30.0
        )

    # Přidání šumu
    if args.noise > 0:
        if args.verbose:
            print(f"Přidávání {args.noise}% gaussovského šumu")
        Z_noisy = add_noise(Z, args.noise, seed=args.seed)
        notes = f"Noise level: {args.noise}%"
        if args.seed is not None:
            notes += f", Random seed: {args.seed}"
    else:
        Z_noisy = Z
        notes = "No noise added"

    # Export
    import os
    base_name = os.path.splitext(args.output)[0]

    if args.format == 'dta' or args.format == 'both':
        dta_file = args.output if args.output.endswith('.DTA') else f"{base_name}.DTA"
        export_gamry_dta(
            dta_file,
            freq,
            Z_noisy,
            circuit_description=circuit_desc,
            notes=notes
        )

    if args.format == 'csv' or args.format == 'both':
        csv_file = args.output if args.output.endswith('.csv') else f"{base_name}.csv"
        export_csv(csv_file, freq, Z_noisy)

    # Vizualizace
    if args.plot:
        if args.verbose:
            print("Zobrazování grafů...")
        plot_data(freq, Z_noisy, title=f"{circuit_desc} (noise: {args.noise}%)")

    # Statistiky
    if args.verbose:
        print("\nStatistiky dat:")
        print(f"  Z' rozsah: [{Z_noisy.real.min():.2f}, {Z_noisy.real.max():.2f}] Ω")
        print(f"  -Z'' rozsah: [{-Z_noisy.imag.max():.2f}, {-Z_noisy.imag.min():.2f}] Ω")
        print(f"  |Z| rozsah: [{np.abs(Z_noisy).min():.2f}, {np.abs(Z_noisy).max():.2f}] Ω")
        print(f"  Phase rozsah: [{np.angle(Z_noisy, deg=True).min():.1f}, {np.angle(Z_noisy, deg=True).max():.1f}] °")

        # Warburg charakteristika při nízkých frekvencích
        low_freq_idx = freq < 1.0  # frekvence < 1 Hz
        if np.any(low_freq_idx):
            phase_low = np.angle(Z_noisy[low_freq_idx], deg=True)
            avg_phase_low = np.mean(phase_low)
            print("\n  Warburg charakteristika (f < 1 Hz):")
            print(f"    Průměrná fáze: {avg_phase_low:.1f}° (ideální -45° pro čistý Warburg)")

    print("\n✓ Hotovo! Použij:")
    if args.format == 'dta' or args.format == 'both':
        output_file = args.output if args.output.endswith('.DTA') else f"{base_name}.DTA"
        print(f"  python3 eis_v2.py {output_file} --auto-lambda --peak-method gmm -v")
    if args.format == 'csv' or (args.format == 'both' and not args.output.endswith('.DTA')):
        output_file = args.output if args.output.endswith('.csv') else f"{base_name}.csv"
        print(f"  python3 eis_v2.py {output_file} --auto-lambda --peak-method gmm -v")


if __name__ == '__main__':
    main()
