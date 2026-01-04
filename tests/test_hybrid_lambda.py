#!/usr/bin/env python3
"""
Test hybridní metody výběru λ (GCV + L-curve) vs čisté GCV.

Testuje na různých typech syntetických dat:
1. Čistá Voigt data (2 RC elementy)
2. Voigt data s šumem
3. Warburg data (kde GCV typicky selhává)
"""

import numpy as np
import matplotlib.pyplot as plt
import logging

from eis_analysis.drt.gcv import (
    find_optimal_lambda_gcv,
    find_optimal_lambda_hybrid,
    compute_lcurve_point
)

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def generate_voigt_impedance(frequencies: np.ndarray,
                              R_inf: float,
                              voigt_elements: list) -> np.ndarray:
    """
    Generuje impedanci pro Voigt obvod.

    Parameters:
    - frequencies: frekvence [Hz]
    - R_inf: vysokofrekvenční odpor [Ω]
    - voigt_elements: list of (R, tau) tuples

    Returns:
    - Z: komplexní impedance
    """
    omega = 2 * np.pi * frequencies
    Z = np.full_like(omega, R_inf, dtype=complex)

    for R, tau in voigt_elements:
        # Voigt element: Z = R / (1 + jωτ)
        Z += R / (1 + 1j * omega * tau)

    return Z


def generate_warburg_impedance(frequencies: np.ndarray,
                                R_inf: float,
                                R_ct: float,
                                C_dl: float,
                                sigma_w: float) -> np.ndarray:
    """
    Generuje impedanci pro Randles obvod s Warburgem.

    Z = R_inf + (R_ct || C_dl) + W
    kde W = σ(1-j)/√ω
    """
    omega = 2 * np.pi * frequencies

    # Warburg: Z_W = sigma * (1 - j) / sqrt(omega)
    Z_W = sigma_w * (1 - 1j) / np.sqrt(omega)

    # R_ct || C_dl
    Z_ct = R_ct / (1 + 1j * omega * R_ct * C_dl)

    # Celková impedance
    Z = R_inf + Z_ct + Z_W

    return Z


def add_noise(Z: np.ndarray, noise_percent: float) -> np.ndarray:
    """Přidá Gaussovský šum proporcionální k |Z|."""
    noise_real = np.random.normal(0, noise_percent/100 * np.abs(Z))
    noise_imag = np.random.normal(0, noise_percent/100 * np.abs(Z))
    return Z + noise_real + 1j * noise_imag


def build_drt_matrices(frequencies: np.ndarray, Z: np.ndarray,
                       R_inf: float, n_tau: int = 100):
    """Sestaví matice A, b, L pro DRT problém."""
    omega = 2 * np.pi * frequencies

    f_max = frequencies.max()
    f_min = frequencies.min()
    tau_min = 1 / (2 * np.pi * f_max)
    tau_max = 1 / (2 * np.pi * f_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    d_ln_tau = np.mean(np.diff(np.log(tau)))

    omega_mesh, tau_mesh = np.meshgrid(omega, tau, indexing='ij')
    denom = 1 + (omega_mesh * tau_mesh)**2
    A_re = d_ln_tau / denom
    A_im = -omega_mesh * tau_mesh * d_ln_tau / denom
    A = np.vstack([A_re, A_im])

    b = np.concatenate([Z.real - R_inf, Z.imag])

    L = np.zeros((n_tau - 2, n_tau))
    for i in range(n_tau - 2):
        L[i, i] = 1
        L[i, i + 1] = -2
        L[i, i + 2] = 1

    return A, b, L, tau


def test_method_comparison(name: str, frequencies: np.ndarray, Z: np.ndarray,
                           R_inf: float, expected_behavior: str = "normal"):
    """
    Porovná GCV vs Hybrid metodu.

    Returns dict s výsledky.
    """
    print(f"\n{'='*60}")
    print(f"TEST: {name}")
    print(f"{'='*60}")

    A, b, L, tau = build_drt_matrices(frequencies, Z, R_inf)

    results = {'name': name, 'expected': expected_behavior}

    # Test GCV
    print("\n--- Čisté GCV ---")
    try:
        lambda_gcv, gcv_score = find_optimal_lambda_gcv(A, b, L)
        results['lambda_gcv'] = lambda_gcv
        results['gcv_score'] = gcv_score
    except Exception as e:
        print(f"GCV selhalo: {e}")
        results['lambda_gcv'] = None
        results['gcv_score'] = None

    # Test Hybrid
    print("\n--- Hybridní (GCV + L-curve) ---")
    try:
        lambda_hybrid, hybrid_score, diag = find_optimal_lambda_hybrid(A, b, L)
        results['lambda_hybrid'] = lambda_hybrid
        results['hybrid_score'] = hybrid_score
        results['method_used'] = diag['method_used']
        results['lambda_lcurve'] = diag['lambda_lcurve']
    except Exception as e:
        print(f"Hybrid selhalo: {e}")
        results['lambda_hybrid'] = None
        results['hybrid_score'] = None
        results['method_used'] = 'failed'

    # Porovnání
    if results['lambda_gcv'] and results['lambda_hybrid']:
        ratio = results['lambda_hybrid'] / results['lambda_gcv']
        print("\n--- Porovnání ---")
        print(f"λ_gcv    = {results['lambda_gcv']:.4e}")
        print(f"λ_hybrid = {results['lambda_hybrid']:.4e}")
        print(f"λ_lcurve = {results['lambda_lcurve']:.4e}")
        print(f"Ratio (hybrid/gcv) = {ratio:.2f}")
        print(f"Metoda použita: {results['method_used']}")
        results['ratio'] = ratio

    return results


def visualize_lcurve(name: str, frequencies: np.ndarray, Z: np.ndarray,
                     R_inf: float, ax: plt.Axes):
    """Vykreslí L-křivku pro daná data."""
    A, b, L, tau = build_drt_matrices(frequencies, Z, R_inf)

    lambda_range = np.logspace(-5, 0, 30)
    rho = []
    eta = []

    for lam in lambda_range:
        log_res, log_reg, _ = compute_lcurve_point(lam, A, b, L)
        rho.append(log_res)
        eta.append(log_reg)

    rho = np.array(rho)
    eta = np.array(eta)

    ax.plot(rho, eta, 'b.-', linewidth=1.5, markersize=4)
    ax.set_xlabel('log||Ax - b|| (reziduum)')
    ax.set_ylabel('log||Lx|| (hladkost)')
    ax.set_title(f'L-křivka: {name}')
    ax.grid(True, alpha=0.3)

    # Označ několik λ hodnot
    for i in [0, len(lambda_range)//4, len(lambda_range)//2,
              3*len(lambda_range)//4, -1]:
        ax.annotate(f'λ={lambda_range[i]:.1e}',
                   (rho[i], eta[i]), fontsize=7, alpha=0.7)


def main():
    np.random.seed(42)

    # Frekvence
    frequencies = np.logspace(5, -2, 71)  # 100kHz - 10mHz

    results_all = []

    # ========== TEST 1: Čistá Voigt data ==========
    print("\n" + "="*70)
    print("SCÉNÁŘ 1: Čistá Voigt data (2 RC elementy)")
    print("="*70)

    R_inf = 100
    voigt_elements = [
        (1000, 1e-3),   # R=1000Ω, τ=1ms (f=159Hz)
        (2000, 1e-1),   # R=2000Ω, τ=100ms (f=1.6Hz)
    ]
    Z_voigt = generate_voigt_impedance(frequencies, R_inf, voigt_elements)

    results = test_method_comparison(
        "Voigt (čistá data)",
        frequencies, Z_voigt, R_inf,
        expected_behavior="GCV a Hybrid by měly dát podobné výsledky"
    )
    results_all.append(results)

    # ========== TEST 2: Voigt s šumem ==========
    print("\n" + "="*70)
    print("SCÉNÁŘ 2: Voigt data s 2% šumem")
    print("="*70)

    Z_noisy = add_noise(Z_voigt, 2.0)

    results = test_method_comparison(
        "Voigt (2% šum)",
        frequencies, Z_noisy, R_inf,
        expected_behavior="Hybrid může preferovat vyšší λ kvůli NNLS"
    )
    results_all.append(results)

    # ========== TEST 3: Warburg data ==========
    print("\n" + "="*70)
    print("SCÉNÁŘ 3: Warburg data (Randles obvod)")
    print("="*70)

    R_inf_w = 50
    Z_warburg = generate_warburg_impedance(
        frequencies,
        R_inf=R_inf_w,
        R_ct=500,
        C_dl=1e-5,
        sigma_w=100
    )

    results = test_method_comparison(
        "Warburg (Randles)",
        frequencies, Z_warburg, R_inf_w,
        expected_behavior="GCV typicky underestimuje λ, Hybrid by měl korigovat"
    )
    results_all.append(results)

    # ========== TEST 4: Warburg s šumem ==========
    print("\n" + "="*70)
    print("SCÉNÁŘ 4: Warburg data s 1% šumem")
    print("="*70)

    Z_warburg_noisy = add_noise(Z_warburg, 1.0)

    results = test_method_comparison(
        "Warburg (1% šum)",
        frequencies, Z_warburg_noisy, R_inf_w,
        expected_behavior="L-curve korekce by měla být výrazná"
    )
    results_all.append(results)

    # ========== TEST 5: Komplexní 3-Voigt ==========
    print("\n" + "="*70)
    print("SCÉNÁŘ 5: Komplexní 3-Voigt s různými τ")
    print("="*70)

    R_inf_3 = 80
    voigt_3 = [
        (500, 1e-4),    # Rychlý proces (f=1.6kHz)
        (1500, 1e-2),   # Střední proces (f=16Hz)
        (3000, 1e0),    # Pomalý proces (f=0.16Hz)
    ]
    Z_3voigt = generate_voigt_impedance(frequencies, R_inf_3, voigt_3)
    Z_3voigt_noisy = add_noise(Z_3voigt, 1.5)

    results = test_method_comparison(
        "3-Voigt (1.5% šum)",
        frequencies, Z_3voigt_noisy, R_inf_3,
        expected_behavior="Komplexní spektrum - test robustnosti"
    )
    results_all.append(results)

    # ========== VIZUALIZACE ==========
    print("\n" + "="*70)
    print("Generuji vizualizace...")
    print("="*70)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # L-křivky pro různé scénáře
    visualize_lcurve("Voigt čistá", frequencies, Z_voigt, R_inf, axes[0, 0])
    visualize_lcurve("Voigt 2% šum", frequencies, Z_noisy, R_inf, axes[0, 1])
    visualize_lcurve("Warburg", frequencies, Z_warburg, R_inf_w, axes[0, 2])

    # Nyquist diagramy
    ax = axes[1, 0]
    ax.plot(Z_voigt.real, -Z_voigt.imag, 'b-', label='Voigt čistá')
    ax.plot(Z_noisy.real, -Z_noisy.imag, 'r.', alpha=0.5, label='Voigt + šum')
    ax.set_xlabel("Z' [Ω]")
    ax.set_ylabel("-Z'' [Ω]")
    ax.set_title("Nyquist: Voigt data")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    ax = axes[1, 1]
    ax.plot(Z_warburg.real, -Z_warburg.imag, 'b-', label='Warburg čistá')
    ax.plot(Z_warburg_noisy.real, -Z_warburg_noisy.imag, 'r.', alpha=0.5, label='Warburg + šum')
    ax.set_xlabel("Z' [Ω]")
    ax.set_ylabel("-Z'' [Ω]")
    ax.set_title("Nyquist: Warburg data")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    ax.plot(Z_3voigt.real, -Z_3voigt.imag, 'b-', label='3-Voigt čistá')
    ax.plot(Z_3voigt_noisy.real, -Z_3voigt_noisy.imag, 'r.', alpha=0.5, label='3-Voigt + šum')
    ax.set_xlabel("Z' [Ω]")
    ax.set_ylabel("-Z'' [Ω]")
    ax.set_title("Nyquist: 3-Voigt data")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('test_hybrid_lambda_results.png', dpi=150)
    print("Uloženo: test_hybrid_lambda_results.png")

    # ========== SOUHRN ==========
    print("\n" + "="*70)
    print("SOUHRN VÝSLEDKŮ")
    print("="*70)
    print(f"{'Scénář':<25} {'λ_gcv':>12} {'λ_hybrid':>12} {'Ratio':>8} {'Metoda':<20}")
    print("-"*70)

    for r in results_all:
        if r.get('lambda_gcv') and r.get('lambda_hybrid'):
            print(f"{r['name']:<25} {r['lambda_gcv']:>12.2e} {r['lambda_hybrid']:>12.2e} "
                  f"{r.get('ratio', 0):>8.2f} {r.get('method_used', 'N/A'):<20}")
        else:
            print(f"{r['name']:<25} {'FAILED':>12} {'FAILED':>12}")

    print("\n" + "="*70)
    print("INTERPRETACE:")
    print("="*70)
    print("- Ratio ≈ 1: GCV a L-curve jsou konzistentní")
    print("- Ratio > 1: L-curve preferuje vyšší λ (více smoothing)")
    print("- Ratio >> 1: GCV pravděpodobně underestimuje λ (NNLS efekt)")
    print("- method_used='lcurve_correction': L-curve výrazně korigovala GCV")
    print("="*70)


if __name__ == "__main__":
    main()
