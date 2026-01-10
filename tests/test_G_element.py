#!/usr/bin/env python3
"""
Test G element (Gerischer for reaction-diffusion processes).

Verifies:
1. Impedance calculation matches analytical formula
2. Limiting behavior (low/high frequency)
3. Jacobian correctness (numerical vs analytical)
4. Circuit fitting with G element
5. Series and parallel combinations work
"""

import numpy as np
import matplotlib.pyplot as plt
from eis_analysis.fitting import R, G, fit_equivalent_circuit

# Generate test frequencies
freq = np.logspace(5, -2, 60)  # 100 kHz to 0.01 Hz

# Test parameters
sigma_val = 100.0  # Pre-factor
tau_val = 1e-3     # Reaction time constant (1 ms)

print("="*70)
print("Test G element (Gerischer for reaction-diffusion)")
print("="*70)
print("Test parameters:")
print(f"  sigma = {sigma_val:.0f} Ohm*s^(1/2)")
print(f"  tau = {tau_val:.3e} s")
print(f"  f_c = 1/(2*pi*tau) = {1/(2*np.pi*tau_val):.2f} Hz")
print()

# Test 1: G element impedance calculation
print("Test 1: G element impedance calculation")
print("-" * 70)
g_elem = G(sigma_val, tau_val)
Z_G = g_elem.impedance(freq, [sigma_val, tau_val])

# Calculate expected impedance
omega = 2 * np.pi * freq
Z_expected = sigma_val / np.sqrt(1 + 1j * omega * tau_val)

diff = np.abs(Z_G - Z_expected)
max_diff = np.max(diff)
print(f"G element: G(sigma={sigma_val:.0f}, tau={tau_val:.3e})")
print(f"  Z at f_min (0.01 Hz): {Z_G[-1]:.2f} Ohm")
print(f"  Z at f_max (100 kHz): {Z_G[0]:.2f} Ohm")
print(f"  Max difference from analytical: {max_diff:.3e} Ohm")
if max_diff < 1e-10:
    print("  PASS: Impedance calculation correct")
else:
    print("  FAIL: Impedance calculation error!")
print()

# Test 2: Limiting behavior
print("Test 2: Limiting behavior")
print("-" * 70)

# Low frequency limit: Z -> sigma
freq_low = np.array([1e-6])
Z_low = g_elem.impedance(freq_low, [sigma_val, tau_val])
print(f"  omega -> 0: Z = {Z_low[0]:.4f} Ohm (expected: sigma = {sigma_val:.0f})")
rel_err_low = np.abs(Z_low[0].real - sigma_val) / sigma_val
print(f"    Real part error: {rel_err_low*100:.4f}%")
print(f"    Imag part: {Z_low[0].imag:.6f} (should be ~0)")

# High frequency: |Z| -> 0
freq_high = np.array([1e8])
Z_high = g_elem.impedance(freq_high, [sigma_val, tau_val])
print(f"  omega -> inf: |Z| = {np.abs(Z_high[0]):.6f} Ohm (expected: ~0)")

if rel_err_low < 0.01 and np.abs(Z_high[0]) < 0.01 * sigma_val:
    print("  PASS: Limiting behavior correct")
else:
    print("  FAIL: Limiting behavior incorrect!")
print()

# Test 3: Gerischer vs Warburg comparison (informational)
print("Test 3: Gerischer vs Warburg comparison")
print("-" * 70)
print("  Gerischer: Z = sigma / sqrt(1 + j*omega*tau)")
print("  Warburg:   Z = sigma * (1 - j) / sqrt(omega)")
print()
print("  Key differences:")
print("    - Gerischer has finite DC resistance (Z(0) = sigma)")
print("    - Warburg has infinite DC impedance (Z(0) -> inf)")
print("    - Gerischer: reaction-diffusion coupling")
print("    - Warburg: pure semi-infinite diffusion")
print()

# Test 4: G element properties
print("Test 4: G element properties")
print("-" * 70)
print(f"  g_elem.sigma = {g_elem.sigma}")
print(f"  g_elem.tau = {g_elem.tau}")
print(f"  g_elem.characteristic_freq = {g_elem.characteristic_freq:.2f} Hz")
print(f"  Expected f_c = {1/(2*np.pi*tau_val):.2f} Hz")
print(f"  g_elem.get_param_labels() = {g_elem.get_param_labels()}")
print(f"  repr(g_elem) = {repr(g_elem)}")
print()

# Test 5: Series circuit: R - G
print("Test 5: Series circuit R - G")
print("-" * 70)
R_s = 10.0
circuit = R(R_s) - G(sigma_val, tau_val)
print(f"  Circuit: {circuit}")
params = circuit.get_all_params()
print(f"  Parameters: {params}")
labels = circuit.get_param_labels()
print(f"  Labels: {labels}")
Z_series = circuit.impedance(freq, params)
print(f"  Z at f_min: {Z_series[-1]:.2f} Ohm")
print(f"  Z at f_max: {Z_series[0]:.2f} Ohm")

# Verify series: Z_total = R_s + Z_G
Z_expected_series = R_s + Z_G
diff_series = np.abs(Z_series - Z_expected_series)
max_diff_series = np.max(diff_series)
print(f"  Max difference from expected: {max_diff_series:.3e}")
if max_diff_series < 1e-10:
    print("  PASS: Series circuit correct")
else:
    print("  FAIL: Series circuit error!")
print()

# Test 6: Parallel circuit: R | G
print("Test 6: Parallel circuit R | G")
print("-" * 70)
R_p = 200.0
circuit_par = R(R_p) | G(sigma_val, tau_val)
print(f"  Circuit: {circuit_par}")
Z_par = circuit_par.impedance(freq, circuit_par.get_all_params())
print(f"  Z at f_min: {Z_par[-1]:.2f} Ohm")
print(f"  Z at f_max: {Z_par[0]:.2f} Ohm")

# Verify parallel: 1/Z_total = 1/R_p + 1/Z_G
Z_R_p = R_p * np.ones_like(freq, dtype=complex)
Z_expected_par = 1 / (1/Z_R_p + 1/Z_G)
diff_par = np.abs(Z_par - Z_expected_par)
max_diff_par = np.max(diff_par)
print(f"  Max difference from expected: {max_diff_par:.3e}")
if max_diff_par < 1e-10:
    print("  PASS: Parallel circuit correct")
else:
    print("  FAIL: Parallel circuit error!")
print()

# Test 7: Jacobian verification (numerical vs analytical)
print("Test 7: Jacobian verification")
print("-" * 70)

from eis_analysis.fitting.jacobian import element_jacobian

Z_jac, dZ_jac = element_jacobian(g_elem, freq, [sigma_val, tau_val])

# Numerical Jacobian
eps = 1e-7
Z_sigma_plus = g_elem.impedance(freq, [sigma_val + eps, tau_val])
Z_sigma_minus = g_elem.impedance(freq, [sigma_val - eps, tau_val])
dZ_dsigma_num = (Z_sigma_plus - Z_sigma_minus) / (2 * eps)

Z_tau_plus = g_elem.impedance(freq, [sigma_val, tau_val + eps])
Z_tau_minus = g_elem.impedance(freq, [sigma_val, tau_val - eps])
dZ_dtau_num = (Z_tau_plus - Z_tau_minus) / (2 * eps)

# Compare
diff_sigma = np.abs(dZ_jac[:, 0] - dZ_dsigma_num)
diff_tau = np.abs(dZ_jac[:, 1] - dZ_dtau_num)
max_diff_sigma = np.max(diff_sigma) / np.max(np.abs(dZ_jac[:, 0]))
max_diff_tau = np.max(diff_tau) / np.max(np.abs(dZ_jac[:, 1]))

print(f"  dZ/dsigma: max relative error = {max_diff_sigma:.3e}")
print(f"  dZ/dtau:   max relative error = {max_diff_tau:.3e}")
if max_diff_sigma < 1e-5 and max_diff_tau < 1e-5:
    print("  PASS: Jacobian correct")
else:
    print("  FAIL: Jacobian error!")
print()

# Test 8: Circuit fitting with G element
print("Test 8: Fit data with G element")
print("-" * 70)

# Generate synthetic data: R_s + G(sigma, tau)
R_s_true = 15.0
sigma_true = 120.0
tau_true = 5e-4

omega = 2 * np.pi * freq
Z_true = R_s_true + sigma_true / np.sqrt(1 + 1j * omega * tau_true)

# Add 1% noise
np.random.seed(42)
noise_level = 0.01
noise = noise_level * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j * np.random.randn(len(freq)))
Z_noisy = Z_true + noise

# Fit with G element (initial guess slightly off)
circuit_fit = R(12) - G(100, 3e-4)
result, Z_fit, fig = fit_equivalent_circuit(freq, Z_noisy, circuit_fit, weighting='modulus', plot=False)

print(f"  True parameters: R_s={R_s_true}, sigma={sigma_true}, tau={tau_true:.3e}")
print(f"  Fitted parameters: {result.params_opt}")
R_s_fit, sigma_fit, tau_fit = result.params_opt
print(f"  R_s_fit   = {R_s_fit:.2f} Ohm (error: {abs(R_s_fit - R_s_true)/R_s_true*100:.1f}%)")
print(f"  sigma_fit = {sigma_fit:.2f} Ohm*s^(1/2) (error: {abs(sigma_fit - sigma_true)/sigma_true*100:.1f}%)")
print(f"  tau_fit   = {tau_fit:.3e} s (error: {abs(tau_fit - tau_true)/tau_true*100:.1f}%)")
print(f"  Fit error (rel): {result.fit_error_rel:.2f}%")
print(f"  Fit error (abs): {result.fit_error_abs:.2f} Ohm")

# Check fit quality
param_errors = [
    abs(R_s_fit - R_s_true) / R_s_true,
    abs(sigma_fit - sigma_true) / sigma_true,
    abs(tau_fit - tau_true) / tau_true
]
max_param_error = max(param_errors)

if result.fit_error_rel < 2.0 and max_param_error < 0.15:
    print("  PASS: Fit successful")
else:
    print("  FAIL: Fit error too large!")
print()

# Test 9: Fixed parameters
print("Test 9: Fixed parameters")
print("-" * 70)
g_fixed = G("100", 1e-3)  # sigma fixed
print(f"  G element with fixed sigma: {g_fixed}")
print(f"  fixed_params = {g_fixed.fixed_params}")
if g_fixed.fixed_params[0] == True and g_fixed.fixed_params[1] == False:
    print("  PASS: Fixed parameter handling correct")
else:
    print("  FAIL: Fixed parameter handling error!")

g_both_fixed = G("100", "1e-3")  # Both fixed
print(f"  G element with both fixed: {g_both_fixed}")
print(f"  fixed_params = {g_both_fixed.fixed_params}")
if g_both_fixed.fixed_params == [True, True]:
    print("  PASS: Both parameters fixed correctly")
else:
    print("  FAIL: Fixed parameter handling error!")
print()

# Visualization
print("Generating visualization...")
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Nyquist plot of G element
ax1 = axes[0, 0]
ax1.plot(Z_G.real, -Z_G.imag, 'b-', linewidth=2, label='G element')
ax1.plot(sigma_val, 0, 'ro', markersize=10, label=f'DC limit (sigma={sigma_val})')
ax1.set_xlabel("Z' [Ohm]")
ax1.set_ylabel("-Z'' [Ohm]")
ax1.set_title("Gerischer element Nyquist plot")
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, sigma_val * 1.1)

# Plot 2: Bode plot
ax2 = axes[0, 1]
ax2.loglog(freq, np.abs(Z_G), 'b-', linewidth=2, label='|Z|')
ax2.axhline(sigma_val, color='r', linestyle=':', alpha=0.5, label=f'DC limit = {sigma_val}')
ax2.axvline(1/(2*np.pi*tau_val), color='g', linestyle='--', alpha=0.5,
            label=f'f_c = {1/(2*np.pi*tau_val):.1f} Hz')
ax2.set_xlabel("Frequency [Hz]")
ax2.set_ylabel("|Z| [Ohm]")
ax2.set_title("Gerischer element Bode magnitude")
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Phase plot
ax3 = axes[1, 0]
phase = np.angle(Z_G, deg=True)
ax3.semilogx(freq, phase, 'b-', linewidth=2)
ax3.axhline(0, color='gray', linestyle=':', alpha=0.5)
ax3.axhline(-45, color='r', linestyle=':', alpha=0.5, label='-45 deg')
ax3.axvline(1/(2*np.pi*tau_val), color='g', linestyle='--', alpha=0.5)
ax3.set_xlabel("Frequency [Hz]")
ax3.set_ylabel("Phase [deg]")
ax3.set_title("Gerischer element phase")
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: Fit result
ax4 = axes[1, 1]
ax4.plot(Z_noisy.real, -Z_noisy.imag, 'o', label='Data (noisy)', markersize=4, alpha=0.5)
ax4.plot(Z_fit.real, -Z_fit.imag, 'r-', label='Fit (R-G)', linewidth=2)
ax4.set_xlabel("Z' [Ohm]")
ax4.set_ylabel("-Z'' [Ohm]")
ax4.set_title(f"Fit with G element (error: {result.fit_error_rel:.2f}%)")
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/orangepi/Lang/python/eis_analysis/test_G_element.png', dpi=150)
print("Plot saved: test_G_element.png")
plt.close()

print()
print("="*70)
print("All tests completed!")
print("="*70)
