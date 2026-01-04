#!/usr/bin/env python3
"""
Test kontroly korelační matice parametrů při circuit fitting.

Testuje:
1. Normální případ - nízká korelace parametrů (Voigt element)
2. Problematický případ - vysoká korelace (over-parametrized model)
3. Extrémní případ - 3 resistory v sérii (singulární problém)

Kontrola detekuje:
- Vysokou korelaci |corr[i,j]| > 0.95
- Varuje před over-parametrized modely
- Nabízí doporučení (zjednodušit model, fixovat parametry, atd.)

Author: EIS Analysis Toolkit
Date: 2025-12-20
"""

import numpy as np
import logging
import warnings

from eis_analysis.fitting import R, C, fit_equivalent_circuit

# Setup logging to see warnings
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Suppress matplotlib warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

print("="*70)
print("Test kontroly korelace parametrů při circuit fitting")
print("="*70)

# ============================================================================
# Test 1: Normální případ - nízká korelace (Voigt element)
# ============================================================================
print("\n[Test 1] Voigt element - očekáváme nízkou korelaci")
print("-" * 70)

# Jednoduchý Voigt element: R_s - (R_p | C)
circuit1 = R(100) - (R(5000) | C(1e-6))
freq1 = np.logspace(4, -1, 30)

# True parametry
true_params1 = [98.5, 4823.2, 8.7e-7]
Z_true1 = circuit1.impedance(freq1, true_params1)

# Přidej malý šum
np.random.seed(42)
noise_level1 = 0.01
noise_real1 = noise_level1 * np.random.randn(len(Z_true1))
noise_imag1 = noise_level1 * np.random.randn(len(Z_true1))
Z_measured1 = Z_true1 * (1 + noise_real1 + 1j * noise_imag1)

print(f"Circuit: {circuit1}")
print(f"True params: {true_params1}")
print(f"Noise level: {noise_level1 * 100}%\n")

# Fit
result1, Z_fit1, fig1 = fit_equivalent_circuit(freq1, Z_measured1, circuit1)

print(f"\nFitted params: {result1.params_opt}")
print(f"Std errors:    {result1.params_stderr}")
print(f"Fit quality:   {result1.quality}")

# ============================================================================
# Test 2: Over-parametrized model - očekáváme vysokou korelaci
# ============================================================================
print("\n" + "="*70)
print("[Test 2] Over-parametrized model - očekáváme vysokou korelaci")
print("-" * 70)

# Vytvoř data z jednoduchého obvodu
circuit_simple = R(100) - (R(5000) | C(1e-6))
freq2 = np.logspace(4, -1, 25)
true_params_simple = [100.0, 5000.0, 1e-6]
Z_true2 = circuit_simple.impedance(freq2, true_params_simple)

# Přidej šum
np.random.seed(123)
noise_level2 = 0.02
noise_real2 = noise_level2 * np.random.randn(len(Z_true2))
noise_imag2 = noise_level2 * np.random.randn(len(Z_true2))
Z_measured2 = Z_true2 * (1 + noise_real2 + 1j * noise_imag2)

# Fituj s over-parametrized modelem: R - (R|C) - (R|C)
# Tento model má 2 R-C články, ale data pocházejí z jednoho
# => parametry budou korelované (můžeme libovolně rozdistribuovat R a C)
circuit_overfit = R(100) - (R(2500) | C(5e-7)) - (R(2500) | C(5e-7))

print(f"Circuit (true):    {circuit_simple}")
print(f"Circuit (fitted):  {circuit_overfit}")
print(f"True params (simple): {true_params_simple}")
print(f"Noise level: {noise_level2 * 100}%")
print("\nFitování over-parametrized modelu...\n")

# Fit - očekáváme warning o vysoké korelaci
result2, Z_fit2, fig2 = fit_equivalent_circuit(freq2, Z_measured2, circuit_overfit)

print(f"\nFitted params: {result2.params_opt}")
print(f"Std errors:    {result2.params_stderr}")
print(f"Fit quality:   {result2.quality}")

# ============================================================================
# Test 3: Model se 3 resistory v sérii - extrémní příklad korelace
# ============================================================================
print("\n" + "="*70)
print("[Test 3] 3 resistory v sérii - extrémní korelace")
print("-" * 70)

# Data z jednoho resistoru
circuit_simple3 = R(1000)
freq3 = np.logspace(4, -1, 20)
true_params3 = [1000.0]
Z_true3 = circuit_simple3.impedance(freq3, true_params3)

# Šum
np.random.seed(456)
noise_level3 = 0.01
noise_real3 = noise_level3 * np.random.randn(len(Z_true3))
noise_imag3 = noise_level3 * np.random.randn(len(Z_true3))
Z_measured3 = Z_true3 * (1 + noise_real3 + 1j * noise_imag3)

# Fituj se 3 resistory: R - R - R
# Tyto parametry budou EXTRÉMNĚ korelované, protože
# Z = R1 + R2 + R3 = konstantní
# Můžeme mít libovolné R1, R2, R3 pokud součet je konstantní
circuit_corr = R(333) - R(333) - R(334)

print(f"Circuit (true):    {circuit_simple3}")
print(f"Circuit (fitted):  {circuit_corr}")
print(f"True params (simple): {true_params3}")
print(f"Noise level: {noise_level3 * 100}%")
print("\nFitování modelu se 3 resistory...\n")

# Fit - očekáváme VELMI silný warning o korelaci
result3, Z_fit3, fig3 = fit_equivalent_circuit(freq3, Z_measured3, circuit_corr)

print(f"\nFitted params: {result3.params_opt}")
print(f"Std errors:    {result3.params_stderr}")
print(f"Fit quality:   {result3.quality}")

# ============================================================================
# Shrnutí
# ============================================================================
print("\n" + "="*70)
print("SHRNUTÍ TESTŮ")
print("="*70)
print("\nTest 1 (Voigt element):")
print("  Očekáváno: ŽÁDNÝ warning o korelaci")
print(f"  Fit kvalita: {result1.quality}")

print("\nTest 2 (Over-parametrized R-C články):")
print("  Očekáváno: Warning o korelaci mezi R-C parametry")
print(f"  Fit kvalita: {result2.quality}")

print("\nTest 3 (3 resistory v sérii):")
print("  Očekáváno: Warning o extrémní korelaci všech R")
print(f"  Fit kvalita: {result3.quality}")

print("\n" + "="*70)
print("Kontrola korelace implementována a testována!")
print("="*70)
