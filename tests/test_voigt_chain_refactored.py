#!/usr/bin/env python3
"""
Test refactored Voigt chain implementation with K elements.

Tests:
1. Basic K elements (replace R||C)
2. --voigt-allow-negative (pseudoinverse)
3. --voigt-include-L (add L element)
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
print("Test Voigt chain refactored implementation")
print("="*70)
print(f"Synthetic data: R_s={R_s_true}, R1={R1_true}, τ1={tau1_true}, R2={R2_true}, τ2={tau2_true}")
print()

# ============================================================================
# Test 1: Basic K elements (default behavior)
# ============================================================================
print("Test 1: Basic K elements (default, NNLS)")
print("-"*70)

circuit1, params1 = fit_voigt_chain_linear(
    freq, Z_noisy,
    n_per_decade=3,
    extend_decades=1.0,
    prune_threshold=0.01,
    weighting='uniform',
    allow_negative=False,  # Default
    include_L=False   # Default
)

print(f"Circuit: {circuit1}")
print(f"Parameters: {len(params1)}")
print()

# Fit
result1, Z_fit1, fig1 = fit_equivalent_circuit(freq, Z_noisy, circuit1, weighting='proportional')
print(f"Fit error: {result1.fit_error_rel:.2f}% (rel), {result1.fit_error_abs:.2f} Ω (abs)")

# Check if K elements are used
circuit_str = str(circuit1)
has_K = 'K(' in circuit_str
has_RC = '|' in circuit_str and 'C(' in circuit_str

print(f"Uses K elements: {has_K}")
print(f"Uses R||C: {has_RC}")

if has_K and not has_RC:
    print("✓ PASS: Circuit uses K elements (not R||C)")
else:
    print("✗ FAIL: Circuit should use K elements!")

print()

# ============================================================================
# Test 2: allow_negative=True (pseudoinverse, Lin-KK style)
# ============================================================================
print("Test 2: allow_negative=True (pseudoinverse)")
print("-"*70)

circuit2, params2 = fit_voigt_chain_linear(
    freq, Z_noisy,
    n_per_decade=5,  # Higher → more parameters → potential negative R_i
    extend_decades=0.0,  # No extension (like Lin-KK)
    prune_threshold=0.0,  # No pruning (like Lin-KK)
    weighting='uniform',
    allow_negative=True,  # Lin-KK style
    include_L=False
)

print(f"Circuit: {circuit2}")
print(f"Parameters: {len(params2)}")

# Check for negative R_i
params_array = np.array(params2)
# Assuming params = [R_s, R_1, τ_1, R_2, τ_2, ...]
R_params = params_array[1::2]  # Every other param starting from index 1
negative_R = R_params[R_params < 0]

print(f"Number of negative R_i: {len(negative_R)}")
if len(negative_R) > 0:
    print(f"  Negative values: {negative_R}")
    print("✓ PASS: allow_negative works (found negative R_i)")
else:
    print("  Note: No negative R_i found (may happen with good data)")

print()

# ============================================================================
# Test 3: include_L=True
# ============================================================================
print("Test 3: include_L=True")
print("-"*70)

circuit3, params3 = fit_voigt_chain_linear(
    freq, Z_noisy,
    n_per_decade=3,
    extend_decades=1.0,
    prune_threshold=0.01,
    weighting='uniform',
    allow_negative=False,
    include_L=True   # Add L
)

print(f"Circuit: {circuit3}")
print(f"Parameters: {len(params3)}")

# Check if L is in circuit
has_L = 'L(' in str(circuit3)
print(f"Has inductance L: {has_L}")

if has_L:
    print("✓ PASS: include_L works (L element added)")
    # Last parameter should be L
    L_value = params3[-1]
    print(f"  L value: {L_value:.3e} H")
else:
    print("✗ FAIL: L element not found in circuit!")

print()

# ============================================================================
# Test 4: Fit quality comparison
# ============================================================================
print("Test 4: Fit quality comparison")
print("-"*70)

# Fit test 2 circuit
result2, Z_fit2, fig2 = fit_equivalent_circuit(freq, Z_noisy, circuit2, weighting='proportional')
print(f"Test 1 (default):       {result1.fit_error_rel:.2f}%")
print(f"Test 2 (allow_negative): {result2.fit_error_rel:.2f}%")

# Fit test 3 circuit
result3, Z_fit3, fig3 = fit_equivalent_circuit(freq, Z_noisy, circuit3, weighting='proportional')
print(f"Test 3 (include_L): {result3.fit_error_rel:.2f}%")

print()

# ============================================================================
# Visualization
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Test 1: K elements
ax1 = axes[0]
ax1.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax1.plot(Z_fit1.real, -Z_fit1.imag, '-', label='Fit (K elements)', linewidth=2)
ax1.set_xlabel("Z' [Ω]")
ax1.set_ylabel("-Z'' [Ω]")
ax1.set_title(f"Test 1: K elements (error: {result1.fit_error_rel:.2f}%)")
ax1.legend()
ax1.grid(True, alpha=0.3)

# Test 2: allow_negative
ax2 = axes[1]
ax2.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax2.plot(Z_fit2.real, -Z_fit2.imag, '-', label='Fit (allow_negative)', linewidth=2)
ax2.set_xlabel("Z' [Ω]")
ax2.set_ylabel("-Z'' [Ω]")
ax2.set_title(f"Test 2: allow_negative (error: {result2.fit_error_rel:.2f}%)")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Test 3: include_L
ax3 = axes[2]
ax3.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data', markersize=4, alpha=0.5)
ax3.plot(Z_fit3.real, -Z_fit3.imag, '-', label='Fit (include_L)', linewidth=2)
ax3.set_xlabel("Z' [Ω]")
ax3.set_ylabel("-Z'' [Ω]")
ax3.set_title(f"Test 3: include_L (error: {result3.fit_error_rel:.2f}%)")
ax3.legend()
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/orangepi/Lang/python/eis_analysis/test_voigt_chain_refactored.png', dpi=150)
print("Plot saved: test_voigt_chain_refactored.png")
plt.close()

print()
print("="*70)
print("All tests completed!")
print("="*70)
