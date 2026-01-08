#!/usr/bin/env python3
"""
Test K element (Voigt with tau parametrization).

Verifies that K(R, τ) is equivalent to (R || C) where C = τ/R.
"""

import numpy as np
import matplotlib.pyplot as plt
from eis_analysis.fitting import R, C, K, fit_equivalent_circuit

# Generate test frequencies
freq = np.logspace(4, -1, 50)  # 10 kHz to 0.1 Hz

# Test parameters
R_val = 1000.0  # 1 kΩ
tau_val = 1e-4  # 100 μs
C_val = tau_val / R_val  # 100 nF

print("="*70)
print("Test K element (Voigt with tau parametrization)")
print("="*70)
print("Test parameters:")
print(f"  R = {R_val:.0f} Ω")
print(f"  τ = {tau_val:.3e} s")
print(f"  C = τ/R = {C_val:.3e} F")
print(f"  f_c = 1/(2πτ) = {1/(2*np.pi*tau_val):.2f} Hz")
print()

# Test 1: K element impedance
print("Test 1: K element impedance calculation")
print("-" * 70)
k_elem = K(R_val, tau_val)
Z_K = k_elem.impedance(freq, [R_val, tau_val])
print(f"K element: K(R={R_val:.0f}, τ={tau_val:.3e})")
print(f"  Z at f_min: {Z_K[0]:.2f} Ω")
print(f"  Z at f_max: {Z_K[-1]:.2f} Ω")
print()

# Test 2: Equivalent (R || C) circuit
print("Test 2: Equivalent (R || C) circuit")
print("-" * 70)
rc_circuit = R(R_val) | C(C_val)
Z_RC = rc_circuit.impedance(freq, [R_val, C_val])
print(f"R||C circuit: (R({R_val:.0f}) | C({C_val:.3e}))")
print(f"  Z at f_min: {Z_RC[0]:.2f} Ω")
print(f"  Z at f_max: {Z_RC[-1]:.2f} Ω")
print()

# Test 3: Verify equivalence
print("Test 3: Verify K(R, τ) ≡ (R || C) where C = τ/R")
print("-" * 70)
diff = np.abs(Z_K - Z_RC)
max_diff = np.max(diff)
rel_diff = max_diff / np.max(np.abs(Z_K))
print(f"  Max absolute difference: {max_diff:.3e} Ω")
print(f"  Max relative difference: {rel_diff:.3e}")
if max_diff < 1e-10:
    print("  ✓ PASS: K element is equivalent to (R || C)")
else:
    print("  ✗ FAIL: Difference too large!")
print()

# Test 4: K element properties
print("Test 4: K element properties")
print("-" * 70)
print(f"  k_elem.capacitance = {k_elem.capacitance:.3e} F")
print(f"  Expected C = {C_val:.3e} F")
print(f"  k_elem.characteristic_freq = {k_elem.characteristic_freq:.2f} Hz")
print(f"  Expected f_c = {1/(2*np.pi*tau_val):.2f} Hz")
print()

# Test 5: Conversion to (R || C)
print("Test 5: K.to_RC() conversion")
print("-" * 70)
rc_from_k = k_elem.to_RC()
print(f"  K element: {k_elem}")
print(f"  Converted: {rc_from_k}")
Z_RC_from_K = rc_from_k.impedance(freq, [R_val, C_val])
diff2 = np.abs(Z_K - Z_RC_from_K)
max_diff2 = np.max(diff2)
print(f"  Max difference: {max_diff2:.3e} Ω")
if max_diff2 < 1e-10:
    print("  ✓ PASS: to_RC() conversion works")
else:
    print("  ✗ FAIL: Conversion failed!")
print()

# Test 6: Series of K elements
print("Test 6: Series of K elements (simple Voigt chain)")
print("-" * 70)
circuit = R(100) - K(500, 1e-4) - K(2000, 1e-3)
print(f"  Circuit: {circuit}")
params = circuit.get_all_params()
print(f"  Parameters: {params}")
Z_chain = circuit.impedance(freq, params)
print(f"  Z at f_min: {Z_chain[0]:.2f} Ω")
print(f"  Z at f_max: {Z_chain[-1]:.2f} Ω")
print()

# Test 7: Circuit fitting with K element
print("Test 7: Fit data with K element")
print("-" * 70)

# Generate synthetic data: R_s + K(R, τ)
R_s_true = 100.0
R_true = 1000.0
tau_true = 2e-4
omega = 2 * np.pi * freq
Z_true = R_s_true + R_true / (1 + 1j * omega * tau_true)

# Add 1% noise
np.random.seed(42)
noise_level = 0.01
noise = noise_level * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
Z_noisy = Z_true + noise

# Fit with K element
circuit_fit = R(90) - K(900, 1.5e-4)  # Initial guess (slightly off)
result, Z_fit, fig = fit_equivalent_circuit(freq, Z_noisy, circuit_fit, weighting='modulus')

print(f"  True parameters: R_s={R_s_true}, R={R_true}, τ={tau_true:.3e}")
print(f"  Fitted parameters: {result.params_opt}")
print(f"  R_s_fit = {result.params_opt[0]:.2f} Ω (error: {abs(result.params_opt[0] - R_s_true)/R_s_true*100:.1f}%)")
print(f"  R_fit   = {result.params_opt[1]:.2f} Ω (error: {abs(result.params_opt[1] - R_true)/R_true*100:.1f}%)")
print(f"  τ_fit   = {result.params_opt[2]:.3e} s (error: {abs(result.params_opt[2] - tau_true)/tau_true*100:.1f}%)")
print(f"  Fit error (rel): {result.fit_error_rel:.2f}%")
print(f"  Fit error (abs): {result.fit_error_abs:.2f} Ω")

if result.fit_error_rel < 2.0:
    print("  ✓ PASS: Fit successful (<2% error)")
else:
    print("  ✗ FAIL: Fit error too large!")
print()

# Visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot 1: Nyquist comparison K vs (R||C)
ax1 = axes[0]
ax1.plot(Z_K.real, -Z_K.imag, 'o-', label='K(R, τ)', markersize=4)
ax1.plot(Z_RC.real, -Z_RC.imag, 'x--', label='R || C', markersize=6)
ax1.set_xlabel("Z' [Ω]")
ax1.set_ylabel("-Z'' [Ω]")
ax1.set_title("K element ≡ (R || C)")
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.axis('equal')

# Plot 2: Bode magnitude
ax2 = axes[1]
ax2.loglog(freq, np.abs(Z_K), 'o-', label='K(R, τ)', markersize=4)
ax2.loglog(freq, np.abs(Z_RC), 'x--', label='R || C', markersize=6)
ax2.axvline(1/(2*np.pi*tau_val), color='r', linestyle=':', alpha=0.5, label=f'f_c = {1/(2*np.pi*tau_val):.1f} Hz')
ax2.set_xlabel("Frequency [Hz]")
ax2.set_ylabel("|Z| [Ω]")
ax2.set_title("Bode magnitude")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Circuit fit
ax3 = axes[2]
ax3.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data (noisy)', markersize=4, alpha=0.5)
ax3.plot(Z_fit.real, -Z_fit.imag, '-', label='Fit (K element)', linewidth=2)
ax3.set_xlabel("Z' [Ω]")
ax3.set_ylabel("-Z'' [Ω]")
ax3.set_title(f"Fit with K element (error: {result.fit_error_rel:.2f}%)")
ax3.legend()
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/orangepi/Lang/python/eis_analysis/test_K_element.png', dpi=150)
print("Plot saved: test_K_element.png")
plt.close()

print()
print("="*70)
print("All tests completed!")
print("="*70)
