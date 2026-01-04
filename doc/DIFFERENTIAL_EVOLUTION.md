# Diferenciální evoluce pro globální optimalizaci EIS obvodů

## 1. Motivace

### 1.1 Problém lokálních minim

Fitování ekvivalentních obvodů na EIS data je **nelineární optimalizační problém** s potenciálně mnoha lokálními minimy:

```
RSS(θ)
   ^
   |     /\
   |    /  \     /\      /\
   |   /    \   /  \    /  \
   |  /      \_/    \__/    \
   | /                       \___  <- globální minimum
   +------------------------------> θ
```

Gradientní metody (jako `least_squares`) konvergují k nejbližšímu lokálnímu minimu, což závisí na initial guess.

### 1.2 Řešení: Globální optimalizace

**Diferenciální evoluce (DE)** je stochastická metoda, která:
- Prohledává celý parametrický prostor
- Nevyžaduje gradient (derivative-free)
- Efektivně uniká z lokálních minim
- Přirozeně respektuje bounds parametrů

### 1.3 Hybridní přístup

Samotná DE konverguje pomalu blízko minima. Proto používáme **hybridní strategii**:

```
1. DE pro globální prohledání -> najde oblast globálního minima
2. least_squares pro lokální refinement -> přesná konvergence + kovariance
```

Tato kombinace využívá silné stránky obou metod.


## 2. Matematické pozadí

### 2.1 Optimalizační problém

Hledáme parametry θ minimalizující vážený součet čtverců:

```
min_θ f(θ) = Σᵢ wᵢ² |Z_data(ωᵢ) - Z_model(ωᵢ, θ)|²
```

kde:
- θ = [θ₁, θ₂, ..., θₚ] jsou parametry obvodu
- wᵢ jsou váhy (uniform, sqrt, proportional)
- Z_model je impedanční funkce obvodu

### 2.2 Algoritmus diferenciální evoluce

DE pracuje s **populací** kandidátních řešení, která se vyvíjí přes generace.

**Inicializace:**
```
x_1^(0) = initial_guess    (z definice obvodu)
Pro i = 2, ..., NP:
    x_i^(0) = L + rand() * (U - L)
```
kde NP = popsize * p je velikost populace, L a U jsou dolní a horní bounds.

**Poznámka:** Jeden člen populace je inicializován hodnotami z definice obvodu (x0 seeding). Pokud máte dobrý odhad, DE ho využije jako výchozí bod pro jednoho "informovaného" jedince, zatímco zbytek populace zajišťuje globální exploraci.

**Mutace:**
Pro každého jedince xᵢ vytvoříme mutanta vᵢ kombinací jiných jedinců:
```
v_i = x_base + F * (x_r1 - x_r2)
```
kde F je mutation factor (typicky 0.5-1.0).

**Křížení (Crossover):**
Kombinujeme mutanta s rodičem:
```
u_i,j = v_i,j   pokud rand() < CR nebo j = j_rand
        x_i,j   jinak
```
kde CR je crossover probability.

**Selekce:**
Pokud je potomek lepší, nahradí rodiče:
```
x_i^(g+1) = u_i   pokud f(u_i) < f(x_i^(g))
            x_i   jinak
```

### 2.3 Strategie

Strategie definuje jak se vybírá **base vector** a kolik difference vectors se používá.

| Strategie | Base | Differences | Charakteristika |
|-----------|------|-------------|-----------------|
| `rand1bin` | náhodný | 1 | Nejlepší explorace, pomalejší |
| `best1bin` | nejlepší | 1 | Rychlá konvergence, může uváznout |
| `randtobest1bin` | náhodný + směr k best | 1 | Vyvážený kompromis |

**rand1bin:**
```
v_i = x_r1 + F * (x_r2 - x_r3)
```
Zcela náhodná explorace prostoru.

**best1bin:**
```
v_i = x_best + F * (x_r1 - x_r2)
```
Vždy směřuje k nejlepšímu jedinci - rychlá, ale riziko předčasné konvergence.

**randtobest1bin (doporučeno):**
```
v_i = x_i + F * (x_best - x_i) + F * (x_r1 - x_r2)
```
Kombinuje směr k nejlepšímu s náhodnou explorací.


## 3. Intuitivní vysvětlení

### 3.1 Populace jako hledači

Představte si populaci jako **skupinu hledačů** prohledávajících horský terén (kde výška = error):

```
     ^  error
     |    A     B
     |   /\    /\      C
     |  /  \  /  \    /\
     | /    \/    \__/  \___  <- globální minimum
     +------------------------> parametry
        ^    ^    ^
        |    |    |
     hledač1 hledač2 hledač3
```

- Každý hledač (jedinec) má pozici v parametrickém prostoru
- Hledači komunikují a sdílejí informace
- Kombinují své pozice pro vytvoření nových kandidátů
- Lepší pozice přežívají, horší jsou nahrazeny

### 3.2 Strategie jako chování hledačů

**rand1bin** - "Nezávislí průzkumníci":
- Každý hledač se pohybuje náhodně
- Dobrá explorace, ale pomalá konvergence
- Vhodné když netušíme kde je minimum

**best1bin** - "Následuj vůdce":
- Všichni směřují k aktuálně nejlepšímu
- Rychlá konvergence k lokálnímu minimu
- Riziko: celá skupina může uvíznout

**randtobest1bin** - "Vyvážený přístup":
- Kombinuje směr k vůdci s náhodným průzkumem
- Dobrý kompromis explorace/exploatace
- Doporučeno pro většinu problémů

### 3.3 Proč hybridní přístup?

DE je jako **horolezec s helikoptérou**:
1. Helikoptéra (DE) ho vysadí blízko nejvyššího vrcholu
2. Pak jde pěšky (least_squares) přímo na vrchol

Samotná helikoptéra nemá přesnost pro přistání na špičce.
Samotná chůze by trvala věčně a mohl by skončit na špatném kopci.


## 4. Implementace

### 4.1 Modul `diffevo.py`

```python
from eis_analysis.fitting import fit_circuit_diffevo, R, C

circuit = R(100) - (R(5000) | C(1e-6))

result, Z_fit, fig = fit_circuit_diffevo(
    circuit, freq, Z,
    strategy=1,      # 1=randtobest1bin, 2=best1bin, 3=rand1bin
    popsize=15,      # Populace = 15 * počet parametrů
    maxiter=1000,    # Maximální počet generací
    tol=0.01,        # Tolerance pro konvergenci
    workers=1,       # Paralelizace (-1 = všechna CPU)
    weighting='sqrt',
    verbose=True,
    use_analytic_jacobian=True  # Analytický Jacobián pro refinement
)

print(f"DE error: {result.de_error:.3f}%")
print(f"Final error: {result.final_error:.3f}%")
print(f"Improvement: {result.improvement:.1f}%")
```

### 4.2 CLI použití

```bash
# Základní použití
./eis.py --optimizer de --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA

# S vlastními parametry
./eis.py --optimizer de \
         --de-strategy 1 \
         --de-popsize 20 \
         --de-maxiter 500 \
         --de-workers -1 \
         --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

### 4.3 Algoritmus implementace

```
1. Inicializace
   - Načtení bounds z definice obvodu (PARAMETER_BOUNDS)
   - Vytvoření cost function pro DE (skalární)
   - Vytvoření residual function pro least_squares (vektorová)

2. Globální hledání (DE)
   - differential_evolution(cost_function, bounds, ...)
   - polish=False (vlastní refinement)
   - Výstup: přibližná pozice globálního minima

3. Lokální refinement (least_squares)
   - least_squares(residual_function, x0=de_result, ...)
   - Přesná konvergence ke skutečnému minimu
   - Výpočet Jacobiánu pro kovarianci

4. Výpočet kovariance
   - SVD-based robustní výpočet z Jacobiánu
   - Standard errors, confidence intervals
   - Condition number pro diagnostiku

5. Výstup
   - FitResult s optimálními parametry
   - DiffEvoResult s DE statistikami
   - Vizualizace fitu
```


## 5. Parametry a jejich vliv

### 5.1 Strategie (`--de-strategy`)

| Hodnota | Název | Kdy použít |
|---------|-------|------------|
| 1 | randtobest1bin | **Default.** Vyvážený, většina problémů |
| 2 | best1bin | Rychlá konvergence, jednodušší problémy |
| 3 | rand1bin | Maximální explorace, komplexní landscape |

**Doporučení:** Začněte s 1, pokud nenachází globální minimum, zkuste 3.

### 5.2 Velikost populace (`--de-popsize`)

```
Skutečná velikost = popsize * počet_parametrů
```

| popsize | Efekt | Vhodné pro |
|---------|-------|------------|
| 5-10 | Rychlé, může minout globální | Jednoduché obvody (2-3 parametry) |
| 15 | **Default.** Dobrý kompromis | Většina obvodů |
| 20-30 | Důkladnější prohledání | Komplexní obvody (>5 parametrů) |
| 50+ | Velmi důkladné, pomalé | Extrémně složité problémy |

**Intuice:** Větší populace = více hledačů = lepší pokrytí prostoru, ale pomalejší.

### 5.3 Maximální iterace (`--de-maxiter`)

| maxiter | Typické použití |
|---------|-----------------|
| 100-200 | Rychlý test, jednoduché obvody |
| 500-1000 | **Default.** Většina problémů |
| 2000+ | Velmi komplexní obvody |

**Poznámka:** DE má early stopping při konvergenci (kontrolováno `--de-tol`).

### 5.4 Tolerance (`--de-tol`)

Relativní změna v nejlepší hodnotě pro detekci konvergence:

| tol | Efekt |
|-----|-------|
| 0.1 | Velmi rychlé zastavení, hrubý výsledek |
| 0.01 | **Default.** Dobrý kompromis |
| 0.001 | Přesnější, ale refinement stejně doladí |

### 5.5 Paralelizace (`--de-workers`)

| workers | Efekt |
|---------|-------|
| 1 | Sekvenční (default) |
| 4 | 4 paralelní workers |
| -1 | Všechna dostupná CPU |

**Poznámka:** Pro malé populace režie převáží zisk. Doporučeno pro popsize > 15.

### 5.6 Analytický Jacobián (`use_analytic_jacobian`)

| Hodnota | Efekt |
|---------|-------|
| True | **Default.** Analytické derivace pro least_squares refinement |
| False | Numerická aproximace (finite differences) |

Analytický Jacobián je použit pouze pro least_squares refinement (krok 2), ne pro DE samotnou (ta je derivative-free).

**Podporované elementy:**
- R, C, L (základní elementy)
- Q (Constant Phase Element)
- W (Warburg)
- Wo (Warburg open/bounded)
- K (Voigt element)


## 6. Srovnání metod

### 6.1 DE vs Multi-start

| Aspekt | Differential Evolution | Multi-start |
|--------|------------------------|-------------|
| **Přístup** | Populační, evoluční | Opakované lokální hledání |
| **Explorace** | Systematická, celý prostor | Perturbace kolem initial |
| **Initial guess** | Volitelný (x0 seeding) | Vyžaduje rozumný odhad |
| **Rychlost** | Pomalejší (hodně evaluací) | Rychlejší (méně evaluací) |
| **Robustnost** | Vyšší pro komplexní landscape | Závisí na initial guess |
| **Paralelizace** | Nativní (populace) | Možná (nezávislé starty) |

**Kdy použít DE:**
- Nemáte dobrý initial guess (nebo ho chcete jen jako hint)
- Obvod má mnoho parametrů (>5)
- Podezření na mnoho lokálních minim
- Máte výpočetní čas

**Kdy použít Multi-start:**
- Máte rozumný initial guess
- Obvod je jednodušší (2-4 parametry)
- Potřebujete rychlý výsledek
- Chcete využít kovarianci pro perturbace

### 6.2 Počet evaluací

Pro obvod s p parametry:

**Multi-start (n restartů):**
```
Evaluace ≈ n * (lokální optimalizace)
         ≈ n * 50-200 evaluací
         ≈ 500-2000 pro n=10
```

**Differential Evolution:**
```
Evaluace = maxiter * popsize * p
         ≈ 1000 * 15 * p
         ≈ 15000 * p pro default nastavení
```

DE vyžaduje více evaluací, ale systematičtěji prohledává prostor.

### 6.3 Praktické srovnání

Pro typický EIS obvod R-(R|C)-(R|Q):

| Metoda | Parametry | Typický čas | Úspěšnost nalezení globálního |
|--------|-----------|-------------|-------------------------------|
| Single fit | 6 | <1s | ~60% (závisí na guess) |
| Multi-start (10) | 6 | 2-5s | ~85% |
| Multi-start (20) | 6 | 4-10s | ~95% |
| DE (default) | 6 | 10-30s | ~99% |
| DE + more popsize | 6 | 30-60s | ~99.9% |


## 7. Diagnostika a interpretace výstupu

### 7.1 Výstupní formát

```
==================================================
Differential Evolution optimization
==================================================
  Strategy: randtobest1bin (option 1)
  Population: 15 * n_params
  Max iterations: 1000
  Tolerance: 0.01
  Workers: 1
  Weighting: sqrt
  Parameters: 3

Running differential evolution...
  DE converged: True
  DE iterations: 127
  DE evaluations: 5715
  DE error: 1.342%

Refining with least_squares (analytic Jacobian)...
  Refined error: 1.320%
  Improvement: +1.6%

==================================================
Differential Evolution results
==================================================
  Strategy: randtobest1bin
  Total evaluations: 5736
  DE error: 1.342% -> Refined: 1.320%

Fit results:
  Parameters:
    R0    = 9.99e+01 +/- 1.02e+00  [95% CI: 9.76e+01, 1.02e+02]
    R1    = 5.00e+03 +/- 9.42e+00  [95% CI: 4.98e+03, 5.02e+03]
    C0    = 1.00e-05 +/- 5.05e-08  [95% CI: 9.91e-06, 1.01e-05]
  Fit error: 1.32% (rel), 28.72 Ohm (abs)
  Quality: Good (<10.0%)
==================================================
```

### 7.2 Klíčové indikátory

**DE converged:**
- `True` - DE dosáhla tolerance před maxiter
- `False` - Vyčerpány iterace, ale refinement může stále pomoci

**DE error vs Refined error:**
- Velký rozdíl (>10%) - DE našla oblast, ale ne přesné minimum
- Malý rozdíl (<2%) - DE už byla blízko

**Improvement:**
- Pozitivní - refinement zlepšil výsledek
- Negativní/nulový - DE už byla optimální (vzácné)

### 7.3 Řešení problémů

**Problém: DE nekonverguje (error stále vysoký)**

Možné příčiny a řešení:
1. **Špatný model** - Obvod neodpovídá datům
   - Zkuste jiný obvod
2. **Nedostatečná explorace**
   - Zvyšte popsize: `--de-popsize 25`
   - Zvyšte maxiter: `--de-maxiter 2000`
   - Zkuste strategii 3: `--de-strategy 3`
3. **Bounds příliš široké**
   - Zkontrolujte rozsahy parametrů

**Problém: DE je příliš pomalá**

Řešení:
1. Snižte popsize: `--de-popsize 10`
2. Snižte maxiter: `--de-maxiter 200`
3. Použijte paralelizaci: `--de-workers -1`
4. Zvyšte tol: `--de-tol 0.05`

**Problém: Velká nejistota parametrů**

To je problém modelu, ne DE. Možné příčiny:
1. Overparametrized model (příliš mnoho parametrů)
2. Korelované parametry
3. Nedostatečný frekvenční rozsah dat


## 8. Matematické detaily

### 8.1 Bounds handling

DE přirozeně respektuje bounds - všichni jedinci jsou vždy uvnitř:

```python
# Inicializace uvnitř bounds
x = lower + rand() * (upper - lower)

# Po mutaci/křížení: clipping
x = np.clip(x, lower, upper)
```

Bounds jsou generovány z `PARAMETER_BOUNDS` v `bounds.py`:

```python
PARAMETER_BOUNDS = {
    'R': (1e-4, 1e10),   # 0.1 mOhm - 10 GOhm
    'C': (1e-15, 1e-1),  # 1 fF - 100 mF
    'Q': (1e-12, 1e-1),  # Q koeficient
    'n': (0.4, 1.0),     # Q exponent
    ...
}
```

### 8.2 Cost function vs Residual function

**Pro DE (skalární):**
```python
def cost_function(params):
    Z_pred = circuit.impedance(freq, params)
    residuals_real = (Z.real - Z_pred.real) * weights
    residuals_imag = (Z.imag - Z_pred.imag) * weights
    return np.sum(residuals_real**2 + residuals_imag**2)
```

**Pro least_squares (vektorová):**
```python
def residual_function(params):
    Z_pred = circuit.impedance(freq, params)
    return np.concatenate([
        (Z.real - Z_pred.real) * weights,
        (Z.imag - Z_pred.imag) * weights
    ])
```

DE minimalizuje skalár, least_squares minimalizuje normu vektoru.

### 8.3 Konvergence least_squares po DE

Po DE máme dobrý starting point, least_squares typicky konverguje za:
- 5-20 iterací (vs 50-200 od špatného startu)
- Přesnější výsledek (kvadratická konvergence blízko minima)
- Jacobián pro výpočet kovariance


## 9. Příklady použití

### 9.1 Jednoduchý RC obvod

```bash
./eis.py --optimizer de --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

Pro 3 parametry, default nastavení obvykle stačí.

### 9.2 Komplexní obvod s Q

```bash
./eis.py --optimizer de \
         --de-popsize 20 \
         --de-maxiter 1500 \
         --circuit 'R(1)-(R(1000)|Q(1e-5,0.9))-(R(10000)|Q(1e-4,0.8))' \
         data.DTA
```

Více parametrů vyžaduje větší populaci a více iterací.

### 9.3 Rychlý test s paralelizací

```bash
./eis.py --optimizer de \
         --de-popsize 10 \
         --de-maxiter 200 \
         --de-workers -1 \
         --circuit 'R(100)-(R(5000)|C(1e-6))' \
         data.DTA
```

Rychlejší výsledek za cenu menší spolehlivosti.

### 9.4 Python API

```python
from eis_analysis import fit_circuit_diffevo, R, C, Q
import numpy as np

# Načtení dat
freq = np.logspace(-2, 5, 61)
Z = ...  # Naměřená impedance

# Definice obvodu (hodnoty = počáteční odhad, ale DE je nepoužije)
circuit = R(1) - (R(1000) | Q(1e-5, 0.9))

# DE optimalizace
result, Z_fit, fig = fit_circuit_diffevo(
    circuit, freq, Z,
    strategy=1,
    popsize=20,
    maxiter=1000,
    workers=-1,
    verbose=True
)

# Výsledky
print(f"Optimální parametry: {result.best_result.params_opt}")
print(f"Standard errors: {result.best_result.params_stderr}")
print(f"Fit error: {result.final_error:.3f}%")
print(f"DE evaluations: {result.n_evaluations}")
```


## 10. Závěr

### 10.1 Silné stránky DE

1. **Globální optimalizace** - Systematicky prohledává celý prostor
2. **Nevyžaduje gradient** - Funguje pro libovolné impedanční funkce
3. **Robustnost** - Efektivně uniká z lokálních minim
4. **Bounds nativně** - Přirozeně respektuje fyzikální omezení
5. **Paralelizovatelná** - Snadné škálování na více CPU
6. **x0 seeding** - Využije initial guess jako startovní bod jednoho jedince

### 10.2 Slabé stránky DE

1. **Pomalejší** - Více evaluací než lokální metody
2. **Stochastická** - Různé běhy mohou dát mírně různé výsledky
3. **Hrubá konvergence** - Vyžaduje refinement pro přesné hodnoty

### 10.3 Doporučení

| Situace | Doporučená metoda |
|---------|-------------------|
| Neznámý obvod, nejistý initial guess | **DE** (využije guess jako x0) |
| Mnoho parametrů (>5) | **DE** s větší populací |
| Známý přibližný obvod, rychlý výsledek | Multi-start |
| Rychlý výsledek potřeba | Multi-start nebo DE s malými parametry |
| Maximální spolehlivost | DE s konzervativním nastavením |

### 10.4 Hybridní workflow

Pro nejlepší výsledky doporučujeme:

```
1. DE pro nalezení globálního minima
2. least_squares pro přesné parametry + kovariance
3. Kontrola diagnostiky (condition number, stderr)
4. Případně multi-start pro verifikaci
```


---

*Dokument vytvořen pro EIS Analysis Toolkit*
*Autor: Claude Code*
*Poslední aktualizace: 2026-01-01*
