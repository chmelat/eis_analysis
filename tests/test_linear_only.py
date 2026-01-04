#!/usr/bin/env python3
"""
Test --voigt-chain s pouze lineární optimalizací (bez nelineárního fitu).
"""

import numpy as np
from eis_analysis.fitting import fit_voigt_chain_linear

# Generate synthetic data
freq = np.logspace(4, -1, 40)
omega = 2 * np.pi * freq

R_s = 100.0
R1 = 500.0
tau1 = 1e-4
R2 = 2000.0
tau2 = 1e-3

Z_true = R_s + R1/(1 + 1j*omega*tau1) + R2/(1 + 1j*omega*tau2)

# Add noise
np.random.seed(42)
noise = 0.01 * np.abs(Z_true) * (np.random.randn(len(freq)) + 1j*np.random.randn(len(freq)))
Z = Z_true + noise

print("="*70)
print("Test: --voigt-chain s pouze lineárním fitem")
print("="*70)

# Test 1: NNLS (allow_negative=False)
print("\nTest 1: NNLS (allow_negative=False)")
print("-"*70)
circuit1, params1 = fit_voigt_chain_linear(
    freq, Z,
    n_per_decade=3,
    allow_negative=False
)
print(f"Circuit: {circuit1}")

# Compute Z_fit manually (what CLI now does)
Z_fit1 = circuit1.impedance(freq, params1)
error1 = np.sqrt(np.mean(np.abs(Z - Z_fit1)**2))
error1_rel = 100 * error1 / np.mean(np.abs(Z))
print(f"Linear fit error: {error1_rel:.2f}% (rel), {error1:.2f} Ω (abs)")

# Test 2: μ optimization + pseudoinverse
print("\nTest 2: μ optimization (allow_negative=True, auto)")
print("-"*70)
circuit2, params2 = fit_voigt_chain_linear(
    freq, Z,
    auto_optimize_M=True,
    mu_threshold=0.85
)
print(f"Circuit: {circuit2}")

# Compute Z_fit manually
Z_fit2 = circuit2.impedance(freq, params2)
error2 = np.sqrt(np.mean(np.abs(Z - Z_fit2)**2))
error2_rel = 100 * error2 / np.mean(np.abs(Z))
print(f"Linear fit error: {error2_rel:.2f}% (rel), {error2:.2f} Ω (abs)")

print("\n" + "="*70)
print("Summary")
print("="*70)
print(f"Test 1 (NNLS):           {error1_rel:.2f}%")
print(f"Test 2 (μ opt + pseudo): {error2_rel:.2f}%")

# Check for negative R
params2_array = np.array(params2)
R_params2 = params2_array[1::2]  # Skip tau values
neg_count = np.sum(R_params2 < 0)
print(f"\nNegative R values in Test 2: {neg_count}")
if neg_count > 0:
    print(f"  Values: {R_params2[R_params2 < 0]}")
    print("  → Tyto elementy by byly oříznuty při nelineárním fitu s bounds R≥0")
    print("  → S lineárním fitem zůstávají, což může být problém pro fyzikální interpretaci")
