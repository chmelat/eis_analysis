#!/usr/bin/env python3
"""
Test μ (mu) optimization for Voigt chain fitting.

Tests the new --voigt-auto-M feature that automatically determines
the optimal number of K elements using the μ metric from Lin-KK.
"""

import numpy as np
import matplotlib.pyplot as plt
from eis_analysis.fitting import fit_voigt_chain_linear, fit_equivalent_circuit

# Generate synthetic data: R_s + K(R1, τ1) + K(R2, τ2)
freq = np.logspace(4, -1, 40)  # 10 kHz to 0.1 Hz
omega = 2 * np.pi * freq

R_s_true = 100.0
R1_true = 500.0
tau1_true = 1e-4  # ~1.6 kHz
R2_true = 2000.0
tau2_true = 1e-3  # ~160 Hz

# True impedance
Z_true = (R_s_true +
          R1_true / (1 + 1j * omega * tau1_true) +
          R2_true / (1 + 1j * omega * tau2_true))

# Add 1% noise
np.random.seed(42)
noise_level = 0.01
noise = noise_level * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
Z_noisy = Z_true + noise

print("="*70)
print("Test μ optimization for Voigt chain")
print("="*70)
print(f"Synthetic data: R_s={R_s_true}, R1={R1_true}, τ1={tau1_true}, R2={R2_true}, τ2={tau2_true}")
print()

# ============================================================================
# Test 1: Default behavior (n_per_decade=3, no μ optimization)
# ============================================================================
print("Test 1: Default (n_per_decade=3, auto_optimize_M=False)")
print("-"*70)

circuit1, params1 = fit_voigt_chain_linear(
    freq, Z_noisy,
    n_per_decade=3,
    extend_decades=1.0,
    prune_threshold=0.01,
    weighting='uniform',
    allow_negative=False,
    include_L=False,
    auto_optimize_M=False  # Default behavior
)

print(f"Circuit: {circuit1}")
print(f"Number of parameters: {len(params1)}")
print()

# Count K elements
circuit_str = str(circuit1)
num_K_elements = circuit_str.count('K(')
print(f"Number of K elements: {num_K_elements}")

# Fit
result1, Z_fit1, fig1 = fit_equivalent_circuit(freq, Z_noisy, circuit1, weighting='proportional')
print(f"Fit error: {result1.fit_error_rel:.2f}% (rel), {result1.fit_error_abs:.2f} Ω (abs)")
print()
plt.close(fig1)

# ============================================================================
# Test 2: μ optimization (auto_optimize_M=True)
# ============================================================================
print("Test 2: μ optimization (auto_optimize_M=True)")
print("-"*70)

circuit2, params2 = fit_voigt_chain_linear(
    freq, Z_noisy,
    extend_decades=1.0,
    prune_threshold=0.01,
    weighting='uniform',
    allow_negative=False,
    include_L=False,
    auto_optimize_M=True,      # Enable μ optimization
    mu_threshold=0.85,         # Default
    max_M=50                   # Default
)

print(f"Circuit: {circuit2}")
print(f"Number of parameters: {len(params2)}")
print()

# Count K elements
num_K_elements2 = str(circuit2).count('K(')
print(f"Number of K elements: {num_K_elements2}")

# Fit
result2, Z_fit2, fig2 = fit_equivalent_circuit(freq, Z_noisy, circuit2, weighting='proportional')
print(f"Fit error: {result2.fit_error_rel:.2f}% (rel), {result2.fit_error_abs:.2f} Ω (abs)")
print()
plt.close(fig2)

# ============================================================================
# Test 3: μ optimization with stricter threshold
# ============================================================================
print("Test 3: μ optimization with stricter threshold (μ ≤ 0.95)")
print("-"*70)

circuit3, params3 = fit_voigt_chain_linear(
    freq, Z_noisy,
    extend_decades=1.0,
    prune_threshold=0.01,
    weighting='uniform',
    allow_negative=False,
    include_L=False,
    auto_optimize_M=True,
    mu_threshold=0.95,         # Stricter (fewer elements)
    max_M=50
)

print(f"Circuit: {circuit3}")
print(f"Number of parameters: {len(params3)}")
print()

# Count K elements
num_K_elements3 = str(circuit3).count('K(')
print(f"Number of K elements: {num_K_elements3}")

# Fit
result3, Z_fit3, fig3 = fit_equivalent_circuit(freq, Z_noisy, circuit3, weighting='proportional')
print(f"Fit error: {result3.fit_error_rel:.2f}% (rel), {result3.fit_error_abs:.2f} Ω (abs)")
print()
plt.close(fig3)

# ============================================================================
# Test 4: μ optimization with allow_negative (Lin-KK style)
# ============================================================================
print("Test 4: μ optimization + allow_negative (full Lin-KK)")
print("-"*70)

circuit4, params4 = fit_voigt_chain_linear(
    freq, Z_noisy,
    extend_decades=0.0,        # No extension (like Lin-KK)
    prune_threshold=0.0,       # No pruning (like Lin-KK)
    weighting='uniform',
    allow_negative=True,       # Lin-KK style
    include_L=False,
    auto_optimize_M=True,
    mu_threshold=0.85,
    max_M=50
)

print(f"Circuit: {circuit4}")
print(f"Number of parameters: {len(params4)}")

# Count K elements
num_K_elements4 = str(circuit4).count('K(')
print(f"Number of K elements: {num_K_elements4}")

# Check for negative R_i
params_array = np.array(params4)
R_params = params_array[1::2]  # R_s at 0, then R_1, τ_1, R_2, τ_2, ...
negative_R = R_params[R_params < 0]
print(f"Number of negative R_i: {len(negative_R)}")

# Fit
result4, Z_fit4, fig4 = fit_equivalent_circuit(freq, Z_noisy, circuit4, weighting='proportional')
print(f"Fit error: {result4.fit_error_rel:.2f}% (rel), {result4.fit_error_abs:.2f} Ω (abs)")
print()
plt.close(fig4)

# ============================================================================
# Comparison
# ============================================================================
print("="*70)
print("Comparison")
print("="*70)
print(f"Test 1 (default, n_per_decade=3):     {num_K_elements} elements, error: {result1.fit_error_rel:.2f}%")
print(f"Test 2 (auto μ, threshold=0.85):      {num_K_elements2} elements, error: {result2.fit_error_rel:.2f}%")
print(f"Test 3 (auto μ, threshold=0.95):      {num_K_elements3} elements, error: {result3.fit_error_rel:.2f}%")
print(f"Test 4 (auto μ + allow_negative):     {num_K_elements4} elements, error: {result4.fit_error_rel:.2f}%")
print()

# Check expectations
print("Validation:")
print("-"*70)

# μ optimization should find reasonable number of elements
if 2 <= num_K_elements2 <= 10:
    print(f"✓ PASS: Test 2 found reasonable number of elements ({num_K_elements2})")
else:
    print(f"? WARNING: Test 2 found unusual number of elements ({num_K_elements2})")

# Stricter threshold should give fewer elements
if num_K_elements3 <= num_K_elements2:
    print(f"✓ PASS: Stricter threshold gave fewer or equal elements ({num_K_elements3} ≤ {num_K_elements2})")
else:
    print("✗ FAIL: Stricter threshold should give fewer elements!")

# Fit quality should be good
if result2.fit_error_rel < 5.0:
    print(f"✓ PASS: μ optimization gives good fit ({result2.fit_error_rel:.2f}% < 5%)")
else:
    print(f"? WARNING: Fit error is high ({result2.fit_error_rel:.2f}%)")

print()

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Test 1: Default
ax1 = axes[0, 0]
ax1.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax1.plot(Z_fit1.real, -Z_fit1.imag, '-', label='Fit', linewidth=2)
ax1.set_xlabel("Z' [Ω]")
ax1.set_ylabel("-Z'' [Ω]")
ax1.set_title(f"Test 1: Default (M={num_K_elements}, error={result1.fit_error_rel:.2f}%)")
ax1.legend()
ax1.grid(True, alpha=0.3)

# Test 2: μ optimization
ax2 = axes[0, 1]
ax2.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax2.plot(Z_fit2.real, -Z_fit2.imag, '-', label='Fit', linewidth=2)
ax2.set_xlabel("Z' [Ω]")
ax2.set_ylabel("-Z'' [Ω]")
ax2.set_title(f"Test 2: μ opt (M={num_K_elements2}, error={result2.fit_error_rel:.2f}%)")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Test 3: Stricter threshold
ax3 = axes[1, 0]
ax3.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax3.plot(Z_fit3.real, -Z_fit3.imag, '-', label='Fit', linewidth=2)
ax3.set_xlabel("Z' [Ω]")
ax3.set_ylabel("-Z'' [Ω]")
ax3.set_title(f"Test 3: μ=0.95 (M={num_K_elements3}, error={result3.fit_error_rel:.2f}%)")
ax3.legend()
ax3.grid(True, alpha=0.3)

# Test 4: allow_negative
ax4 = axes[1, 1]
ax4.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax4.plot(Z_fit4.real, -Z_fit4.imag, '-', label='Fit', linewidth=2)
ax4.set_xlabel("Z' [Ω]")
ax4.set_ylabel("-Z'' [Ω]")
ax4.set_title(f"Test 4: Lin-KK (M={num_K_elements4}, error={result4.fit_error_rel:.2f}%)")
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/orangepi/Lang/python/eis_analysis/test_mu_optimization.png', dpi=150)
print("Plot saved: test_mu_optimization.png")
plt.close()

print()
print("="*70)
print("All tests completed!")
print("="*70)
