#!/usr/bin/env python3
"""
Test Voigt chain initial guess estimation.

Tests:
1. generate_tau_grid() - správné pokrytí frekvenčního rozsahu
2. compute_voigt_matrix() - správná matice pro lineární problém
3. estimate_R_linear() - NNLS regression pro známá data
4. fit_voigt_chain_linear() - kompletní workflow
5. Recovery test - dokáže odhadnout parametry z čistých dat?
6. Noisy data test - robustnost vůči šumu

Author: EIS Analysis Toolkit
Date: 2025-12-20
"""

import numpy as np
import sys
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

print("="*70)
print("Test Voigt chain initial guess estimation")
print("="*70)

try:
    from eis_analysis.fitting import fit_equivalent_circuit
    from eis_analysis.fitting.voigt_chain import (
        generate_tau_grid,
        compute_voigt_matrix,
        estimate_R_linear,
        fit_voigt_chain_linear
    )

    # ========================================================================
    # Test 1: generate_tau_grid()
    # ========================================================================
    print("\n[Test 1] generate_tau_grid() - pokrytí frekvenčního rozsahu")
    print("-" * 70)

    freq1 = np.logspace(4, -1, 30)  # 10 kHz to 0.1 Hz
    tau1 = generate_tau_grid(freq1, n_per_decade=3, extend_decades=1.0)

    print(f"  Frequency range: [{freq1.min():.1e}, {freq1.max():.1e}] Hz")
    print(f"  Generated {len(tau1)} τ values")
    print(f"  τ range: [{tau1.min():.3e}, {tau1.max():.3e}] s")

    # Verify τ corresponds to frequency
    f_from_tau_min = 1 / (2 * np.pi * tau1.max())
    f_from_tau_max = 1 / (2 * np.pi * tau1.min())
    print(f"  Equivalent f range: [{f_from_tau_min:.1e}, {f_from_tau_max:.1e}] Hz")

    # Check extension
    expected_f_min = freq1.min() / 10  # 1 decade extension
    assert f_from_tau_min < expected_f_min * 1.5, \
        f"Extension failed: {f_from_tau_min:.2e} should be < {expected_f_min:.2e}"

    print(f"  ✓ Extension verified: covers down to {f_from_tau_min:.1e} Hz")

    # ========================================================================
    # Test 2: compute_voigt_matrix()
    # ========================================================================
    print("\n[Test 2] compute_voigt_matrix() - správná matice")
    print("-" * 70)

    freq2 = np.array([0.1, 1.0, 10.0, 100.0])
    tau2 = np.array([0.01, 0.1, 1.0])
    A2 = compute_voigt_matrix(freq2, tau2, include_Rs=True)

    print(f"  Frequencies: {freq2}")
    print(f"  τ values: {tau2}")
    print(f"  Matrix shape: {A2.shape}")
    print(f"  Expected: ({len(freq2)}, {len(tau2) + 1})")

    assert A2.shape == (len(freq2), len(tau2) + 1), "Matrix shape mismatch!"

    # First column should be all ones (R_s)
    assert np.allclose(A2[:, 0], 1.0), "First column should be 1.0 for R_s"
    print("  ✓ First column (R_s): all ones")

    # Check h(ω, τ) values are in [0, 1]
    assert np.all(A2[:, 1:] >= 0) and np.all(A2[:, 1:] <= 1), \
        "h(ω, τ) should be in [0, 1]"
    print("  ✓ h(ω, τ) values in [0, 1]")

    # ========================================================================
    # Test 3: estimate_R_linear() - známá data
    # ========================================================================
    print("\n[Test 3] estimate_R_linear() - NNLS regression")
    print("-" * 70)

    # Generate synthetic data from known Voigt chain
    # R_s = 100 Ω
    # Voigt 1: R=1000 Ω, τ=0.001 s → C=1e-6 F
    # Voigt 2: R=2000 Ω, τ=0.01 s  → C=5e-6 F
    freq3 = np.logspace(4, 0, 40)
    omega3 = 2 * np.pi * freq3

    R_s_true = 100.0
    R_true = np.array([1000.0, 2000.0])
    tau_true = np.array([0.001, 0.01])

    # Compute Z
    Z3 = np.zeros(len(freq3), dtype=complex)
    Z3 += R_s_true  # R_s contribution

    for R_i, tau_i in zip(R_true, tau_true):
        # Voigt: Z_i = R_i / (1 + jωτ_i)
        Z3 += R_i / (1 + 1j * omega3 * tau_i)

    print(f"  True R_s: {R_s_true} Ω")
    print(f"  True R_i: {R_true}")
    print(f"  True τ_i: {tau_true}")

    # Estimate using same τ values (include_L=False since synthetic data has no inductance)
    R_est, residual, L_value = estimate_R_linear(freq3, Z3, tau_true, include_Rs=True, include_L=False)

    print(f"\n  Estimated R_s: {R_est[0]:.3f} Ω")
    print(f"  Estimated R_i: {R_est[1:]}")
    print(f"  Residual: {residual:.3e}")
    print(f"  L (inductance): {L_value}")  # Should be None

    # Check recovery (should be nearly perfect for clean data with exact τ)
    assert np.allclose(R_est[0], R_s_true, rtol=0.01), \
        f"R_s mismatch: {R_est[0]:.1f} vs {R_s_true}"
    assert np.allclose(R_est[1:], R_true, rtol=0.01), \
        f"R_i mismatch: {R_est[1:]} vs {R_true}"

    print("  ✓ Perfect recovery from clean data!")

    # ========================================================================
    # Test 4: fit_voigt_chain_linear() - kompletní workflow
    # ========================================================================
    print("\n[Test 4] fit_voigt_chain_linear() - kompletní workflow")
    print("-" * 70)

    # Use same synthetic data as Test 3
    circuit4, params4 = fit_voigt_chain_linear(
        freq3, Z3,
        n_per_decade=3,
        extend_decades=1.0,
        prune_threshold=0.01
    )

    print(f"\n  Circuit: {circuit4}")
    print(f"  Number of parameters: {len(params4)}")
    print(f"  Parameters: {params4}")

    # Should have R_s and 2 main Voigt elements (might have extras from finer grid)
    assert len(params4) >= 5, "Should have at least R_s + 2 Voigt elements (5 params)"
    print("  ✓ Circuit built successfully")

    # ========================================================================
    # Test 5: Recovery test - noisy data
    # ========================================================================
    print("\n[Test 5] Recovery test - noisy data")
    print("-" * 70)

    # Add noise to synthetic data
    np.random.seed(42)
    noise_level = 0.02  # 2%
    noise_real = noise_level * np.random.randn(len(Z3))
    noise_imag = noise_level * np.random.randn(len(Z3))
    Z3_noisy = Z3 * (1 + noise_real + 1j * noise_imag)

    print(f"  Noise level: {noise_level * 100}%")

    # Estimate initial guess
    circuit5, params5 = fit_voigt_chain_linear(
        freq3, Z3_noisy,
        n_per_decade=2,  # Fewer points to avoid overfitting
        extend_decades=0.5,
        prune_threshold=0.05  # More aggressive pruning
    )

    print(f"\n  Circuit: {circuit5}")
    print(f"  Number of parameters: {len(params5)}")

    # Use as initial guess for final fit
    print("\n  Using as initial guess for final optimization...")
    result5, Z_fit5, fig5 = fit_equivalent_circuit(freq3, Z3_noisy, circuit5)

    print(f"\n  Final fit quality: {result5.quality}")
    print(f"  Final fit error: {result5.fit_error_rel:.2f}%")

    # Should achieve good fit
    assert result5.fit_error_rel < 5.0, \
        f"Fit error too high: {result5.fit_error_rel:.2f}%"

    print("  ✓ Good fit achieved with automatic initial guess!")

    # ========================================================================
    # Test 6: Realistic example - 3 time constants
    # ========================================================================
    print("\n[Test 6] Realistic example - 3 separated time constants")
    print("-" * 70)

    # Three well-separated processes
    freq6 = np.logspace(5, -2, 50)
    omega6 = 2 * np.pi * freq6

    R_s_true6 = 50.0
    R_true6 = np.array([500.0, 1500.0, 3000.0])
    tau_true6 = np.array([1e-4, 1e-2, 1.0])  # 3 decades apart

    print("  True parameters:")
    print(f"    R_s = {R_s_true6} Ω")
    for i, (R_i, tau_i) in enumerate(zip(R_true6, tau_true6)):
        C_i = tau_i / R_i
        f_peak = 1 / (2 * np.pi * tau_i)
        print(f"    Voigt {i+1}: R={R_i:.0f} Ω, C={C_i:.2e} F, τ={tau_i:.2e} s (f_peak={f_peak:.2e} Hz)")

    # Generate data
    Z6 = np.zeros(len(freq6), dtype=complex)
    Z6 += R_s_true6

    for R_i, tau_i in zip(R_true6, tau_true6):
        Z6 += R_i / (1 + 1j * omega6 * tau_i)

    # Add realistic noise
    np.random.seed(123)
    noise6 = 0.01
    Z6_noisy = Z6 * (1 + noise6 * np.random.randn(len(Z6)) + 1j * noise6 * np.random.randn(len(Z6)))

    print("\n  Estimating initial guess...")
    circuit6, params6 = fit_voigt_chain_linear(
        freq6, Z6_noisy,
        n_per_decade=3,
        extend_decades=1.0,
        prune_threshold=0.02
    )

    print("\n  Final optimization...")
    result6, Z_fit6, fig6 = fit_equivalent_circuit(freq6, Z6_noisy, circuit6)

    print("\n  Fit results:")
    print(f"    Quality: {result6.quality}")
    print(f"    Error: {result6.fit_error_rel:.2f}%")
    print(f"    Optimized parameters: {len(result6.params_opt)}")

    # Check if major peaks are recovered
    assert result6.fit_error_rel < 3.0, \
        f"Fit error too high: {result6.fit_error_rel:.2f}%"

    print("  ✓ Excellent fit achieved!")

    # ========================================================================
    # Shrnutí
    # ========================================================================
    print("\n" + "="*70)
    print("VŠECHNY TESTY PROŠLY!")
    print("="*70)
    print("\nImplementované funkce:")
    print("  ✓ generate_tau_grid() - logarithmická mřížka τ")
    print("  ✓ compute_voigt_matrix() - matice pro lineární problém")
    print("  ✓ estimate_R_linear() - NNLS regression pro R_i")
    print("  ✓ fit_voigt_chain_linear() - kompletní workflow")
    print("\nVýhody:")
    print("  + Lineární problém → stabilní, rychlá konvergence")
    print("  + Automatický initial guess → žádné manuální ladění")
    print("  + Física-informed τ grid → pokrývá celé spektrum")
    print("  + Pruning → eliminuje zbytečné parametry")
    print("  + NNLS → garantuje R_i ≥ 0")
    print("="*70)

    sys.exit(0)

except Exception as e:
    print(f"\n{'='*70}")
    print("CHYBA: Test selhal!")
    print(f"{'='*70}")
    print(f"Exception: {type(e).__name__}: {e}")
    import traceback
    print("\nTraceback:")
    traceback.print_exc()
    print(f"{'='*70}")
    sys.exit(1)
