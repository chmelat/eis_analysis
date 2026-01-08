#!/usr/bin/env python3
"""
Test script to compare different weighting types for circuit fitting.
"""

import matplotlib.pyplot as plt
from eis_analysis.io.data_loading import load_data
from eis_analysis.fitting.circuit_elements import R, C
from eis_analysis.fitting import fit_equivalent_circuit

# Load test data
print("Loading test data...")
freq, Z = load_data('real_gamry_example.DTA')

# Simple Voigt circuit for testing
circuit = R(100) - (R(5000) | C(1e-6))

# Test all weighting types
weighting_types = ['uniform', 'sqrt', 'proportional', 'modulus']
results = {}

print("\n" + "="*70)
print("Testing different weighting types")
print("="*70)

for wtype in weighting_types:
    print(f"\n{'='*70}")
    print(f"Testing: {wtype}")
    print(f"{'='*70}")

    result, Z_fit, fig = fit_equivalent_circuit(
        freq, Z, circuit, weighting=wtype
    )

    results[wtype] = {
        'params': result.params_opt,
        'error_rel': result.fit_error_rel,
        'error_abs': result.fit_error_abs,
        'quality': result.quality
    }

    plt.close(fig)  # Close individual figures

# Summary comparison
print("\n" + "="*70)
print("COMPARISON SUMMARY")
print("="*70)

print("\nParameters:")
print(f"{'Type':<15} {'R0 [Ω]':<15} {'R1 [Ω]':<15} {'C [F]':<15}")
print("-"*70)
for wtype in weighting_types:
    params = results[wtype]['params']
    print(f"{wtype:<15} {params[0]:<15.2e} {params[1]:<15.2e} {params[2]:<15.2e}")

print("\nFit Quality:")
print(f"{'Type':<15} {'Rel Error [%]':<20} {'Abs Error [Ω]':<20} {'Quality':<15}")
print("-"*70)
for wtype in weighting_types:
    r = results[wtype]
    print(f"{wtype:<15} {r['error_rel']:<20.4f} {r['error_abs']:<20.4f} {r['quality']:<15}")

print("\n" + "="*70)
print("Test completed successfully!")
print("="*70)
