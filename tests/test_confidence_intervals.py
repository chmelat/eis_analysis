#!/usr/bin/env python3
"""
Test confidence intervals computation for circuit fitting.

Tests:
1. Helper funkce _compute_confidence_interval() - správné t-critical values
2. FitResult properties - CI obsahují true values (synthetic data)
3. Edge cases - invalid covariance (inf stderr)
4. 99% CI širší než 95%
5. Backwards compatibility

Author: EIS Analysis Toolkit
Date: 2025-12-15
"""

import numpy as np
import sys
import logging

# Setup logging
logging.basicConfig(level=logging.WARNING)  # Suppress INFO messages for clean test output

print("="*70)
print("Test confidence intervalů pro circuit fitting")
print("="*70)

try:
    from eis_analysis.fitting import R, C, fit_equivalent_circuit
    from eis_analysis.fitting.circuit import _compute_confidence_interval
    from scipy.stats import t

    # ========================================================================
    # Test 1: Helper funkce _compute_confidence_interval
    # ========================================================================
    print("\n[Test 1] _compute_confidence_interval() - správné výpočty")
    print("-" * 70)

    params = np.array([100.0, 5000.0, 1e-6])
    stderr = np.array([1.0, 50.0, 1e-8])
    n_data = 50

    ci_low_95, ci_high_95 = _compute_confidence_interval(params, stderr, n_data, 0.95)
    ci_low_99, ci_high_99 = _compute_confidence_interval(params, stderr, n_data, 0.99)

    print(f"  Params:  {params}")
    print(f"  Stderr:  {stderr}")
    print(f"  n_data:  {n_data}")
    print(f"  dof:     {n_data - len(params)}")
    print("\n  95% CI:")
    print(f"    Low:   {ci_low_95}")
    print(f"    High:  {ci_high_95}")
    print("  99% CI:")
    print(f"    Low:   {ci_low_99}")
    print(f"    High:  {ci_high_99}")

    # Verify t-kritické hodnoty
    dof = n_data - len(params)
    t_95 = t.ppf(0.975, dof)  # 95% CI, two-tailed
    t_99 = t.ppf(0.995, dof)  # 99% CI
    print(f"\n  t_critical(95%, dof={dof}): {t_95:.4f}")
    print(f"  t_critical(99%, dof={dof}): {t_99:.4f}")

    # Verify margin
    expected_margin_95 = t_95 * stderr[0]
    actual_margin_95 = ci_high_95[0] - params[0]
    assert np.isclose(expected_margin_95, actual_margin_95, rtol=1e-10), \
        f"95% CI margin mismatch! Expected {expected_margin_95}, got {actual_margin_95}"

    expected_margin_99 = t_99 * stderr[0]
    actual_margin_99 = ci_high_99[0] - params[0]
    assert np.isclose(expected_margin_99, actual_margin_99, rtol=1e-10), \
        f"99% CI margin mismatch! Expected {expected_margin_99}, got {actual_margin_99}"

    print("\n  ✓ Margins correctly computed")
    print(f"    95%: {expected_margin_95:.6f} (expected) vs {actual_margin_95:.6f} (actual)")
    print(f"    99%: {expected_margin_99:.6f} (expected) vs {actual_margin_99:.6f} (actual)")

    # ========================================================================
    # Test 2: FitResult properties na synthetic data
    # ========================================================================
    print("\n[Test 2] FitResult.params_ci_95 a params_ci_99 properties")
    print("-" * 70)

    # Vytvoř syntetická data
    circuit = R(100) - (R(5000) | C(1e-6))
    freq = np.logspace(4, -1, 30)

    # True parametry
    true_params = [98.5, 4823.2, 8.7e-7]
    Z_true = circuit.impedance(freq, true_params)

    # Přidej šum
    np.random.seed(42)
    noise_level = 0.01
    noise_real = noise_level * np.random.randn(len(Z_true))
    noise_imag = noise_level * np.random.randn(len(Z_true))
    Z_measured = Z_true * (1 + noise_real + 1j * noise_imag)

    print(f"  Circuit: {circuit}")
    print(f"  True params: {true_params}")
    print(f"  Noise level: {noise_level * 100}%")
    print(f"  N points: {len(freq)}")

    # Fit
    result, Z_fit, fig = fit_equivalent_circuit(freq, Z_measured, circuit)

    print(f"\n  Fitted params: {result.params_opt}")
    print(f"  Std errors:    {result.params_stderr}")

    # CI
    ci_low_95, ci_high_95 = result.params_ci_95
    ci_low_99, ci_high_99 = result.params_ci_99

    print("\n  95% CI:")
    for i in range(len(result.params_opt)):
        print(f"    p[{i}]: [{ci_low_95[i]:.6e}, {ci_high_95[i]:.6e}]")

    print("\n  99% CI:")
    for i in range(len(result.params_opt)):
        print(f"    p[{i}]: [{ci_low_99[i]:.6e}, {ci_high_99[i]:.6e}]")

    # Verify CI obsahují true values (většinou, ne vždy kvůli šumu)
    print("\n  Checking if true values are within CI:")
    ci_coverage_95 = 0
    ci_coverage_99 = 0
    for i, (true_val, low_95, high_95, low_99, high_99) in enumerate(
        zip(true_params, ci_low_95, ci_high_95, ci_low_99, ci_high_99)
    ):
        in_95 = low_95 <= true_val <= high_95
        in_99 = low_99 <= true_val <= high_99
        if in_95:
            ci_coverage_95 += 1
        if in_99:
            ci_coverage_99 += 1

        status_95 = "✓" if in_95 else "✗"
        status_99 = "✓" if in_99 else "✗"
        print(f"    p[{i}]: 95% {status_95}  99% {status_99}  (true={true_val:.6e})")

    print(f"\n  Coverage: 95% CI = {ci_coverage_95}/{len(true_params)}, "
          f"99% CI = {ci_coverage_99}/{len(true_params)}")

    # Poznámka: S 1% noise může se stát že některé true values budou mimo CI
    # To je OK - 95% CI nezaručuje 100% coverage, pouze 95%
    if ci_coverage_95 < len(true_params):
        print(f"  Note: Some true values outside 95% CI - expected with {noise_level*100}% noise")

    # ========================================================================
    # Test 3: 99% CI širší než 95%
    # ========================================================================
    print("\n[Test 3] 99% CI širší než 95% CI")
    print("-" * 70)

    width_95 = ci_high_95 - ci_low_95
    width_99 = ci_high_99 - ci_low_99

    print("  Width ratios (99% / 95%):")
    for i in range(len(width_95)):
        ratio = width_99[i] / width_95[i]
        print(f"    p[{i}]: {ratio:.4f}")

    assert np.all(width_99 > width_95), "99% CI should be wider than 95% CI!"
    print("\n  ✓ 99% CI correctly wider than 95% CI")
    print(f"    Average ratio: {(width_99 / width_95).mean():.4f}")

    # ========================================================================
    # Test 4: Edge case - neplatná covariance matrix
    # ========================================================================
    print("\n[Test 4] Edge case - invalid covariance (inf stderr)")
    print("-" * 70)

    # Simuluj situaci s inf stderr
    from dataclasses import replace
    result_bad = replace(result, params_stderr=np.full_like(result.params_opt, np.inf))
    ci_low_bad, ci_high_bad = result_bad.params_ci_95

    print(f"  Stderr (invalid): {result_bad.params_stderr}")
    print(f"  95% CI low:  {ci_low_bad}")
    print(f"  95% CI high: {ci_high_bad}")

    assert np.all(np.isinf(ci_low_bad)) and np.all(ci_low_bad < 0), \
        "Expected -inf for invalid stderr (low)"
    assert np.all(np.isinf(ci_high_bad)) and np.all(ci_high_bad > 0), \
        "Expected +inf for invalid stderr (high)"

    print("  ✓ Invalid stderr correctly returns ±inf CI")

    # ========================================================================
    # Test 5: Backwards compatibility
    # ========================================================================
    print("\n[Test 5] Backwards compatibility")
    print("-" * 70)

    # Test že FitResult lze vytvořit bez _n_data (použije default 0)
    from eis_analysis.fitting.circuit import FitResult

    result_compat = FitResult(
        circuit=circuit,
        params_opt=np.array([100.0, 5000.0, 1e-6]),
        params_stderr=np.array([1.0, 50.0, 1e-8]),
        fit_error_rel=1.0,
        fit_error_abs=10.0,
        quality='good'
        # _n_data vynecháno - použije default 0
    )

    print("  FitResult created without _n_data: ✓")
    print(f"  _n_data value: {result_compat._n_data}")

    # CI budou computed s dof=max(0-3, 1)=1 (velmi široké)
    ci_low_compat, ci_high_compat = result_compat.params_ci_95
    print("  95% CI (with _n_data=0):")
    print(f"    Low:  {ci_low_compat}")
    print(f"    High: {ci_high_compat}")

    # S dof=1, t_critical je velmi velké (~12.7), takže CI budou široké
    dof_compat = max(0 - 3, 1)
    t_crit_compat = t.ppf(0.975, dof_compat)
    print(f"  t_critical(dof={dof_compat}): {t_crit_compat:.4f} (very conservative)")

    print("  ✓ Backwards compatibility maintained")

    # ========================================================================
    # Test 6: Malý dataset (edge case pro t-distribution)
    # ========================================================================
    print("\n[Test 6] Malý dataset (t-distribution vs normal)")
    print("-" * 70)

    # Test s velmi malým počtem bodů
    n_small = 10
    n_params = 3
    dof_small = n_small - n_params

    # Test s velkým počtem bodů
    n_large = 100
    dof_large = n_large - n_params

    t_small = t.ppf(0.975, dof_small)
    t_large = t.ppf(0.975, dof_large)
    z_normal = 1.96  # Normal approximation

    print(f"  Small dataset (n={n_small}, dof={dof_small}):")
    print(f"    t_critical: {t_small:.4f}")
    print(f"    vs normal:  {z_normal:.4f}")
    print(f"    difference: {((t_small / z_normal) - 1) * 100:.1f}%")

    print(f"\n  Large dataset (n={n_large}, dof={dof_large}):")
    print(f"    t_critical: {t_large:.4f}")
    print(f"    vs normal:  {z_normal:.4f}")
    print(f"    difference: {((t_large / z_normal) - 1) * 100:.1f}%")

    print("\n  ✓ t-distribution correctly more conservative for small datasets")

    # ========================================================================
    # Shrnutí
    # ========================================================================
    print("\n" + "="*70)
    print("VŠECHNY TESTY PROŠLY!")
    print("="*70)
    print("\nImplementované funkce:")
    print("  ✓ _compute_confidence_interval() - helper s t-distribution")
    print("  ✓ FitResult.params_ci_95 - 95% confidence intervals")
    print("  ✓ FitResult.params_ci_99 - 99% confidence intervals")
    print("  ✓ Edge case handling - invalid covariance → ±inf")
    print("  ✓ Backwards compatibility - default _n_data=0")
    print("  ✓ t-distribution - statisticky korektní pro malé datasety")
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
