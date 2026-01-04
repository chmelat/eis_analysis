# Analýza metody Lin-KK v knihovně impedance.py

## Přehled

Lin-KK test (Linear Kramers-Kronig) z knihovny `impedance.py` používá lineární regresi pro validaci EIS dat. Tato analýza je založena na zdrojovém kódu verze nainstalované v systému.

**Reference:**
Schönleber, M. et al. "A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20–27 (2014)
doi: 10.1016/j.electacta.2014.01.034

---

## 1. Model ekvivalentního obvodu

Lin-KK test fituje data pomocí následujícího modelu:

```
Z_fit = R_0 + L + Σ(k=1..M) [R_k / (1 + jωτ_k)]
```

Kde:
- `R_0` = seriový odpor (ohmic resistance)
- `L` = seriová indukčnost (zachycuje měřicí artefakty)
- `M` = počet RC elementů
- `R_k` = odpor k-tého RC elementu (fitovaný parametr)
- `τ_k` = časová konstanta k-tého RC elementu (PEVNÁ, nefitovaná)

**KLÍČOVÉ:** Časové konstanty τ_k jsou fixní a logaritmicky distribuované. Fitují se pouze odpory R_k a seriové komponenty.

---

## 2. Distribuce časových konstant

### Implementace (validation.py:112-124)

```python
def get_tc_distribution(f, M):
    """ Returns the distribution of time constants for the linKK method """

    t_max = 1/(2 * np.pi * np.min(f))
    t_min = 1/(2 * np.pi * np.max(f))
    ts = np.zeros(shape=(M,))
    ts[0] = t_min
    ts[-1] = t_max
    if M > 1:
        for k in range(2, M):
            ts[k-1] = 10**(np.log10(t_min) +
                           ((k-1)/(M-1))*np.log10(t_max/t_min))
    return ts
```

### Matematický vzorec

```
τ_1 = 1 / (2π * f_max)           (nejrychlejší proces)
τ_M = 1 / (2π * f_min)           (nejpomalejší proces)

Pro k = 2, 3, ..., M-1:
τ_k = 10^[log₁₀(τ_min) + (k-1)/(M-1) * log₁₀(τ_max/τ_min)]
```

**Vlastnosti:**
- Logaritmické rozložení v intervalu [τ_min, τ_max]
- Pokrývá přesně měřený frekvenční rozsah (žádné rozšíření!)
- Rovnoměrné rozložení v log(τ) prostoru
- M časových konstant → M RC elementů

**Rozdíl od --voigt-chain:**
- LinKK: **žádné** rozšíření rozsahu (extend_decades = 0)
- Voigt chain: rozšíření o +1 dekádu směrem k nízkým frekvencím

---

## 3. Určení počtu RC elementů M

### Iterativní algoritmus (validation.py:91-98)

```python
if c is not None:
    M = 0
    mu = 1
    while mu > c and M < max_M:
        M += 1
        ts = get_tc_distribution(f, M)
        elements, mu = fit_linKK(f, ts, M, Z, fit_type, add_cap)

        if M % 10 == 0:
            print(M, mu, rmse(eval_linKK(elements, ts, f), Z))
```

### Metrika kvality: μ (mu)

```python
def calc_mu(Rs):
    """ Calculates mu for use in LinKK """

    neg_sum = sum(abs(x) for x in Rs if x < 0)
    pos_sum = sum(abs(x) for x in Rs if x >= 0)

    return 1 - neg_sum/pos_sum
```

**Matematicky:**
```
μ = 1 - (Σ|R_k| pro R_k < 0) / (Σ|R_k| pro R_k ≥ 0)
```

**Interpretace:**
- μ → 1: Všechny R_k kladné → dobrý fit, fyzikálně korektní
- μ → 0: Velká záporná masa → overfit (příliš mnoho parametrů)
- μ < 0: Dominance záporných R_k → špatný fit

**Stop podmínka:**
- μ < c (výchozí c = 0.85)
- Když μ klesne pod práh, ukončí se iterace
- Výsledek: M = nejmenší počet RC elementů s přijatelným μ

**Proč záporné R_k znamenají overfit?**
- Fyzikálně by měly být všechny R_k ≥ 0 (odpor nemůže být záporný)
- Záporné R_k vznikají, když model má příliš mnoho volnosti
- Overfit: model fituje šum místo skutečného signálu

---

## 4. Lineární regrese

### Formulace problému (validation.py:127-263)

Lin-KK převádí nelineární problém na **lineární** tím, že fixuje τ_k a fituje pouze R_k.

**Normalizace:**
```
Z_norm = Z / |Z|
```
Všechna data jsou normalizována impedancí, aby se snížil vliv různých řádů velikosti.

### Design matice A

Pro fit_type = 'real' (výchozí):

```python
a_re = np.zeros((N_freq, M+2))

# První sloupec: R_0 (seriový odpor)
a_re[:, 0] = 1 / np.abs(Z)

# Sloupce 1..M: RC elementy
for i, tau in enumerate(ts):
    Z_RC = K([1, tau], f)  # R=1, τ=tau
    a_re[:, i+1] = Z_RC.real / np.abs(Z)

# Poslední sloupec: L (indukčnost) - pouze v imaginární části
```

**Element K (RC element):**
```python
def K(p, f):
    omega = 2 * np.pi * np.array(f)
    R, tau_k = p[0], p[1]
    Z = R / (1 + 1j * omega * tau_k)
    return Z
```

Pro R=1:
```
Z_RC(ω, τ_k) = 1 / (1 + jωτ_k)
Z_RC.real = 1 / (1 + (ωτ_k)²)
Z_RC.imag = -ωτ_k / (1 + (ωτ_k)²)
```

### Řešení pomocí pseudoinverze

**Reálná část (fit_type = 'real'):**
```python
elements = np.linalg.pinv(a_re).dot(Z.real / np.abs(Z))
```

**Matematicky:**
```
A_re @ x = b_re
x = (A_re)^+ @ b_re

Kde:
- A_re[k, 0] = 1/|Z_k|                          (pro R_0)
- A_re[k, i] = Re[1/(1+jω_k*τ_i)] / |Z_k|      (pro R_i, i=1..M)
- b_re[k] = Z_k.real / |Z_k|
- x = [R_0, R_1, R_2, ..., R_M]
```

**Pseudoinverze:**
```
A^+ = (A^T @ A)^(-1) @ A^T    (pro overdetermined systém, N > M+2)
```

**Po fitu reálné části:**
Dopočítá se L (indukčnost) z imaginární části:
```python
Z_fit_re = eval_linKK(elements, ts, f)
coefs = np.linalg.pinv(a_im).dot((Z.imag - Z_fit_re.imag) / np.abs(Z))
elements[-1] = coefs[-1]  # L
```

---

## 5. Tři režimy fittingu

### fit_type = 'real' (výchozí, validation.py:213-230)

1. Fituje pouze Re(Z/|Z|) → získá [R_0, R_1, ..., R_M]
2. Vypočítá Z_fit z těchto R_k
3. Fituje residuum Im(Z - Z_fit) → získá L (a volitelně C)

**Výhody:**
- Reálná část je obvykle méně zašuměná
- Rychlejší konvergence
- Stabilnější výsledky

### fit_type = 'imag' (validation.py:231-246)

1. Fituje pouze Im(Z/|Z|) → získá [R_1, ..., R_M, L]
2. Vypočítá Z_fit z RC elementů (bez R_0)
3. Fituje R_0 pomocí váhované regrese (Boukamp et al. 1995)

```python
ws = 1 / (Z.real**2 + Z.imag**2)  # Váhy
R_0 = Σ(w_i * (Z_i.real - Z_fit_i.real)) / Σ(w_i)
```

### fit_type = 'complex' (validation.py:247-253)

Fituje Re(Z) a Im(Z) současně pomocí kombinované least squares:

```python
# Minimalizuje: ||A_re @ x - b_re||² + ||A_im @ x - b_im||²
x = (A_re^T @ A_re + A_im^T @ A_im)^(-1) @
    (A_re^T @ b_re + A_im^T @ b_im)
```

**Matematicky (Eq 14 Schönleber):**
```
r = ||A' @ x - b'||² + ||A'' @ x - b''||²

∂r/∂x = 0  →  x = [(A')^T @ A' + (A'')^T @ A'']^(-1) @
                   [(A')^T @ b' + (A'')^T @ b'']
```

Kde ' = reálná část, '' = imaginární část.

---

## 6. Evaluace a rezidua

### Výpočet Z_fit (validation.py:266-280)

```python
def eval_linKK(elements, ts, f):
    """ Builds a circuit of RC elements to be used in LinKK """
    circuit_string = 's([R({},{}),'.format([elements[0]], f.tolist())

    for (Rk, tk) in zip(elements[1:], ts):
        circuit_string += f'K({[Rk, tk]},{f.tolist()}),'

    circuit_string += 'L({},{}),'.format([elements[-1]], f.tolist())
    if elements.size == (ts.size + 3):
        circuit_string += 'C({},{}),'.format([1/elements[-2]], f.tolist())

    circuit_string = circuit_string.strip(',')
    circuit_string += '])'

    return eval(circuit_string, circuit_elements)
```

Vytvoří circuit string: `s([R(R_0), K([R_1, τ_1]), ..., K([R_M, τ_M]), L(L)])`

### Výpočet reziduí (validation.py:283-297)

```python
def residuals_linKK(elements, ts, Z, f, residuals='real'):
    """ Calculates the residual between the data and a LinKK fit """

    err = (Z - eval_linKK(elements, ts, f))/np.abs(Z)

    if residuals == 'real':
        return err.real
    elif residuals == 'imag':
        return err.imag
```

**Matematicky:**
```
res_real = Re(Z - Z_fit) / |Z| × 100%
res_imag = Im(Z - Z_fit) / |Z| × 100%
```

**Interpretace:**
- |res| < 1% → data jsou KK-compliant, kvalitní
- |res| ≥ 1% → možné artefakty, nelinearita, nestabilita

---

## 7. Porovnání s --voigt-chain

| Aspekt | Lin-KK (impedance.py) | --voigt-chain (eis_analysis) |
|--------|----------------------|------------------------------|
| **Distribuce τ** | Logaritmická, f_min → f_max | Logaritmická, f_min/10 → f_max |
| **Rozšíření rozsahu** | **Ne** (přesně měřený rozsah) | **Ano** (+1 dekáda k nízkým f) |
| **Počet τ** | Určen iterativně (μ < 0.85) | Určen vzorcem (n_per_decade × dekády) |
| **Fit metoda** | Pseudoinverze (pinv) | NNLS (non-negative LS) |
| **Constraint R_i ≥ 0** | **Ne** (záporné R_i používá pro detekci overfitu) | **Ano** (fyzikální constraint) |
| **Normalizace** | Z/\|Z\| (proporcionální váhy) | Volitelné (uniform/sqrt/proportional/square) |
| **Seriové elementy** | R_0 + L (+ volitelně C) | Pouze R_s |
| **Fit co?** | real/imag/complex | Pouze real (Z') |
| **Pruning** | **Ne** (všechny RC elementy použity) | **Ano** (threshold = 1% max) |
| **Primární účel** | **Validace kvality dat** | **Automatický návrh obvodu** |

---

## 8. Klíčové rozdíly v přístupu

### Lin-KK (impedance.py)

**Filozofie:**
- Minimalistický rozsah τ (přesně f_min → f_max)
- Záporné R_k jsou **diagnostický nástroj** (ne chyba!)
- Iterativní hledání optimálního M pomocí μ metriky
-Fit imaginární části zachycuje indukčnost (měřicí artefakty)

**Výhody:**
+ Automatické určení optimálního M
+ Robustní detekce overfitu (μ metrika)
+ Zachycuje měřicí artefakty (L, C)
+ Dokumentovaná metoda (Schönleber 2014)

**Nevýhody:**
- Může vrátit záporné R_k (není to bug, ale feature!)
- Nepodporuje různé váhování
- Neprune malé příspěvky

### --voigt-chain (eis_analysis)

**Filozofie:**
- Rozšířený rozsah τ (zachycuje pomalé procesy mimo měřený rozsah)
- Fyzikální constraint R_i ≥ 0 (NNLS)
- Pruning malých příspěvků (zjednodušení modelu)
- Fit pouze Z' (stabilnější než Z'')

**Výhody:**
+ Fyzikálně korektní (R_i ≥ 0)
+ Flexibilní váhování
+ Automatické zjednodušení (pruning)
+ Rozšíření rozsahu → lepší zachycení pomalých procesů

**Nevýhody:**
- Fixní počet τ (není automaticky optimalizován)
- Může mít příliš mnoho parametrů (před pruningem)
- Nedetekuje měřicí artefakty (L)

---

## 9. Numerická ilustrace

### Příklad: 5 dekád frekvencí (10 kHz → 0.1 Hz)

**Lin-KK:**
```
f_min = 0.1 Hz,  f_max = 10000 Hz
τ_min = 1/(2π×10000) ≈ 1.59e-5 s
τ_max = 1/(2π×0.1) ≈ 1.59 s

Iterace:
M=1: μ = 0.95 > 0.85 → pokračuj
M=2: μ = 0.92 > 0.85 → pokračuj
M=3: μ = 0.89 > 0.85 → pokračuj
M=4: μ = 0.86 > 0.85 → pokračuj
M=5: μ = 0.83 < 0.85 → STOP

Výsledek: M = 5 RC elementů

τ distribuce:
τ_1 = 1.59e-5 s  (f = 10 kHz)
τ_2 = 1.59e-4 s  (f = 1 kHz)
τ_3 = 1.59e-3 s  (f = 100 Hz)
τ_4 = 1.59e-2 s  (f = 10 Hz)
τ_5 = 1.59e-1 s  (f = 1 Hz)
(Poznámka: τ_max by bylo 1.59 s pro 0.1 Hz, ale M=5 stop před tím)
```

**--voigt-chain:**
```
f_min_extended = 0.1 / 10^1 = 0.01 Hz  (rozšíření!)
τ_min = 1/(2π×10000) ≈ 1.59e-5 s
τ_max = 1/(2π×0.01) ≈ 15.9 s  (větší!)

n_decades = log₁₀(15.9 / 1.59e-5) = 6
n_tau = ceil(6 × 3) + 1 = 19 bodů

Po NNLS fit:
R_i rozsah: [0.01, 5000] Ω  (některé malé)

Po pruning (1%):
Odstraněno: 11 elementů
Zachováno: 8 elementů

Výsledek: 8 Voigt elementů (po pruningu)
```

**Porovnání:**
- Lin-KK: **5 RC elementů** (iterativně určeno)
- Voigt chain: **8 Voigt elementů** (19 → 8 po pruningu)
- Voigt chain má širší pokrytí (0.01 Hz vs 0.1 Hz)

---

## 10. Závěr

**Lin-KK metoda je elegantní implementace Schönleber et al. (2014):**

1. **Distribuce τ:** Logaritmická v rozsahu měřených frekvencí
2. **Počet M:** Iterativně určen pomocí μ metriky (poměr záporné/kladné masy)
3. **Řešení:** Lineární regrese s pseudoinverzí (numpy.linalg.pinv)
4. **Normalizace:** Z/|Z| pro vyvážení různých řádů velikosti
5. **Fit typy:** real/imag/complex (různé strategie pro robustnost)
6. **Měřicí artefakty:** Automaticky zachycuje L (a volitelně C)

**Hlavní rozdíl od --voigt-chain:**
- Lin-KK: **detekce kvality** → minimalistický, iterativní, diagnostický
- Voigt chain: **návrh obvodu** → extensivní, fyzikálně constrained, prunovaný

Obě metody používají lineární regresi, ale s různými cíli a filosofií!

**Co dělá empirický vzorec M ≈ log₁₀(Δf)/0.85?**
- To je **heuristika** z literatury, ne přímý vzorec z kódu!
- Kód skutečně určuje M iterativně pomocí μ metriky
- Vzorec je dobrá aproximace výsledku iterace (pro typická data)
- c = 0.85 je práh pro μ, ne dělitel v počtu dekád
