# K Element - Voigt element s τ parametrizací

## Přehled

Element **K(R, τ)** je alternativní parametrizace Voigtova elementu (R||C), která používá časovou konstantu τ místo kapacity C.

**Impedance:**
```
Z_K(ω) = R / (1 + jωτ)
```

**Ekvivalence:**
```
K(R, τ) ≡ (R || C)  kde  C = τ/R
```

---

## Výhody K(R, τ) vs (R||C)

### 1. Fyzikální intuitivnost

```python
# Klasický Voigt: Jaká je charakteristická frekvence?
circuit = R(1000) | C(1e-7)  # τ = R×C = 1e-4 s → f = ? (musím počítat)

# K element: Přímý vztah k frekvenci
circuit = K(1000, 1e-4)      # τ = 1e-4 s → f = 1/(2πτ) = 1.59 kHz ✓
```

**Časová konstanta τ má jasný význam:**
- τ určuje charakteristickou frekvenci: **f_c = 1/(2πτ)**
- Při f = f_c: impedance je polovina maxima, fáze = -45°
- τ přímo odpovídá relaxačnímu času procesu

### 2. Lepší separace parametrů při fittingu

```python
# Problém s (R||C):
# Změna R → mění se f_c (pokud C zůstane konstantní)
R = 1000 Ω, C = 1e-7 F  →  τ = 1e-4 s, f_c = 1.6 kHz
R = 2000 Ω, C = 1e-7 F  →  τ = 2e-4 s, f_c = 800 Hz  (frekvence se změnila!)

# S K(R, τ):
# R a τ jsou nezávislé → lepší kondicionování
K(1000, 1e-4)  →  f_c = 1.6 kHz
K(2000, 1e-4)  →  f_c = 1.6 kHz  (frekvence zachována!)
```

**Parametry mají jasné role:**
- **R** kontroluje amplitudu (výšku semicircle v Nyquist)
- **τ** kontroluje frekvenci (pozici semicircle na frekvenční ose)

### 3. Konzistence s DRT analýzou

DRT (Distribution of Relaxation Times) používá τ jako primární parametr:

```
γ(τ) = distribuce časových konstant
```

**K element je přirozeným mostem mezi DRT a circuit fittingem:**
- DRT peak při τ_i → K(R_i, τ_i) element
- Integrace DRT píku → R_i = ∫ γ(τ) d(ln τ)
- Přímo použitelné pro konstrukci obvodu

### 4. Lin-KK kompatibilita

Lin-KK test (Schönleber et al. 2014) používá právě K(R, τ) parametrizaci:

```
Z_LinKK = R_0 + Σ K(R_k, τ_k)
```

K element umožňuje přímé použití výsledků z Lin-KK testu.

---

## Použití

### Základní syntaxe

```python
from eis_analysis.fitting import K, R

# Vytvoření K elementu
k = K(1000, 1e-4)  # R=1kΩ, τ=100μs

# Výchozí hodnoty
k = K()  # R=1kΩ, τ=100μs (výchozí)

# Fixované parametry (string)
k = K("1000", 1e-4)    # R fixní, τ volný
k = K("1000", "1e-4")  # Oba fixní
```

### Vlastnosti

```python
k = K(1000, 1e-4)

# Ekvivalentní kapacita
C = k.capacitance  # 1e-7 F (τ/R)

# Charakteristická frekvence
f_c = k.characteristic_freq  # 1591.5 Hz (1/(2πτ))

# Konverze na (R||C)
rc_circuit = k.to_RC()  # (R(1000) | C(1e-7))
```

### Výpočet impedance

```python
import numpy as np

freq = np.logspace(4, -1, 50)  # 10 kHz → 0.1 Hz

# Přímý výpočet
Z = k.impedance(freq, [1000, 1e-4])

# Nebo pomocí circuit objektu
circuit = R(100) - K(1000, 1e-4)
params = circuit.get_all_params()
Z = circuit.impedance(freq, params)
```

---

## Příklady

### Příklad 1: Jednoduchý Voigt element

```python
from eis_analysis.fitting import R, K, fit_equivalent_circuit
import numpy as np

# Obvod: R_s + K(R, τ)
circuit = R(100) - K(1000, 1e-4)

# Fitování dat
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Výsledné parametry
R_s, R, tau = result.params_opt
C = tau / R  # Ekvivalentní kapacita

print(f"R_s = {R_s:.2f} Ω")
print(f"R = {R:.2f} Ω")
print(f"τ = {tau:.3e} s")
print(f"f_c = {1/(2*np.pi*tau):.1f} Hz")
print(f"C = {C:.3e} F")
```

### Příklad 2: Voigt chain (série K elementů)

```python
# Obvod: R_s + K_1 + K_2 + K_3
# Každý K element modeluje jeden relaxační proces
circuit = (R(100) -
           K(500, 1e-5) -   # Rychlý proces (~3.2 kHz)
           K(1000, 1e-4) -  # Střední proces (~1.6 kHz)
           K(2000, 1e-3))   # Pomalý proces (~160 Hz)

# Fit
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit, weighting='proportional')

# Analýza procesů
params = result.params_opt
R_s = params[0]
for i in range(3):
    R_i = params[1 + 2*i]
    tau_i = params[2 + 2*i]
    f_i = 1/(2*np.pi*tau_i)
    print(f"Proces {i+1}: R={R_i:.1f} Ω, τ={tau_i:.3e} s, f={f_i:.1f} Hz")
```

### Příklad 3: Konverze z DRT výsledků

```python
# Předpokládejme DRT analýzu s 3 píky
drt_peaks = [
    {'tau': 1e-5, 'R': 500},   # Peak 1
    {'tau': 1e-4, 'R': 1000},  # Peak 2
    {'tau': 1e-3, 'R': 2000},  # Peak 3
]

# Sestavení obvodu z DRT píků
circuit = R(100)  # R_s
for peak in drt_peaks:
    circuit = circuit - K(peak['R'], peak['tau'])

print(circuit)
# Output: R(100) - K(R=500, τ=1e-05) - K(R=1000, τ=0.0001) - K(R=2000, τ=0.001)

# Použij jako initial guess pro fitting
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

### Příklad 4: Randles circuit s K elementem

```python
from eis_analysis.fitting import Q, W

# Randles circuit: R_s - (K || Q) - W
# K modeluje charge transfer proces
# Q modeluje double layer capacitor (non-ideal)
# W modeluje difúzi
circuit = R(10) - (K(100, 1e-3) | Q(1e-4, 0.85)) - W(50)

result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

### Příklad 5: Fixování parametrů

```python
# Scenario: Z DRT víme τ, ale chceme fitovat R
# Fixujeme τ, ale necháme R volný

circuit = R(100) - K(1000, "1e-4")  # τ fixní (string!)

result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

# Výsledek: fituje se pouze R_s a R (τ zůstane 1e-4)
R_s_fit, R_fit = result.params_opt  # Pouze 2 parametry
```

---

## Vztah k (R||C) elementu

### Matematická ekvivalence

```
K(R, τ):   Z = R / (1 + jωτ)
R||C:      Z = R / (1 + jωRC)
```

**Podmínka ekvivalence:**
```
τ = R × C
```

### Konverze

**Z K na (R||C):**
```python
k = K(1000, 1e-4)
C = k.capacitance  # C = τ/R = 1e-7 F
rc = k.to_RC()     # (R(1000) | C(1e-7))
```

**Z (R||C) na K:**
```python
rc = R(1000) | C(1e-7)
tau = 1000 * 1e-7  # τ = R×C = 1e-4 s
k = K(1000, 1e-4)
```

---

## Kdy použít K vs (R||C)?

### Použij K(R, τ) když:

+ **Máš DRT výsledky** - τ je přímo z DRT píků
+ **Znáš charakteristickou frekvenci** - snadný převod f → τ = 1/(2πf)
+ **Chceš separovat amplitudu a frekvenci** - lepší kondicionování fittingu
+ **Implementuješ Lin-KK test** - kompatibilita s literaturou
+ **Potřebuješ fyzikální interpretaci** - τ je relaxační čas

### Použij (R||C) když:

+ **Máš elektrochemický model** - C má přímý fyzikální význam (double layer, oxide capacitance)
+ **Kapacita je známá hodnota** - např. z geometrie elektrody
+ **Kompatibilita se starým kódem** - klasická parametrizace
+ **Didaktické účely** - R||C je standardní notation

---

## Příklady výpočtů

### Převody mezi parametry

```python
# Dáno: K(R=1000, τ=1e-4)
R = 1000  # Ω
tau = 1e-4  # s

# Výpočet:
C = tau / R                  # 1e-7 F = 100 nF
f_c = 1 / (2 * np.pi * tau)  # 1591.5 Hz
omega_c = 1 / tau            # 10000 rad/s

# Impedance při f_c:
Z_fc = R / (1 + 1j)          # R/√2 * (1 - j)
|Z_fc| = R / np.sqrt(2)      # 707.1 Ω
phase = -45°                 # -π/4
```

### Časové konstanty pro různé frekvence

```python
# Chci K element s charakteristickou frekvencí 1 kHz
f_c = 1000  # Hz
tau = 1 / (2 * np.pi * f_c)  # 1.59e-4 s

# Různé R hodnoty (všechny mají f_c = 1 kHz)
K(100, tau)    # C = 1.59 μF
K(1000, tau)   # C = 159 nF
K(10000, tau)  # C = 15.9 nF
```

### Voigt chain s dekádovým pokrytím

```python
# Pokrytí 4 dekády frekvencí (10 Hz - 100 kHz)
# 1 K element na dekádu

f_values = [10, 100, 1000, 10000]  # Hz
tau_values = [1/(2*np.pi*f) for f in f_values]

circuit = R(100)
for tau in tau_values:
    circuit = circuit - K(1000, tau)  # Všechny stejné R

# tau_values ≈ [0.016, 0.0016, 0.00016, 0.000016] s
```

---

## Reference

1. Schönleber, M. et al. "A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests." *Electrochimica Acta* 131, 20–27 (2014). doi: 10.1016/j.electacta.2014.01.034

2. Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance Data Validation." *J. Electrochem. Soc.* 142 (6), 1885-1894 (1995).

3. Ciucci, F. "Modeling electrochemical impedance spectroscopy." *Current Opinion in Electrochemistry* 13, 132-139 (2019).

---

## Implementační poznámky

### Numerical stability

K(R, τ) parametrizace má obecně **lepší numerickou stabilitu** než (R||C):

```python
# Problém s (R||C) při velkých R:
R = 1e6 Ω, C = 1e-10 F  →  τ = 1e-4 s
# Při fittingu: gradient wrt C je velmi malý (C << 1)

# S K(R, τ):
K(1e6, 1e-4)
# Gradient wrt τ je lépe škálovaný
```

### Bounds recommendation

Pro fitting s K elementem:

```python
# Typické bounds pro K(R, τ)
R_bounds = (0, 1e8)      # 0 Ω - 100 MΩ
tau_bounds = (1e-8, 1e2) # 10 ns - 100 s (pokrývá ~13 dekád frekvencí)

# Ekvivalentní bounds pro C:
# C_lower = tau_lower / R_upper = 1e-8 / 1e8 = 1e-16 F
# C_upper = tau_upper / R_lower = 1e2 / 0 → ∞ (problém!)
```

K parametrizace přirozeně vyhýbá problémům s boundsy.

---

## Changelog

**v4.2.0 (2025-12-20):**
- Přidán K element (Voigt s τ parametrizací)
- Vlastnosti: `.capacitance`, `.characteristic_freq`
- Metoda: `.to_RC()` pro konverzi na (R||C)
- Kompatibilní s operátory `-`, `|`, `*`
- Úplná dokumentace a testy
