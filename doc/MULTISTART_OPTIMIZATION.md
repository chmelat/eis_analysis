# Multi-start optimalizace pro fitování EIS obvodů

## 1. Motivace

### 1.1 Problém lokálních minim

Fitování ekvivalentních obvodů je nelineární optimalizace s potenciálně mnoha lokálními minimy:

```
RSS(θ)
   ^
   |     /\
   |    /  \     /\
   |   /    \   /  \    <- lokální minima
   |  /      \_/    \
   | /                \___  <- globální minimum
   +-----------------------> θ
       ^
       |
    initial guess
```

Gradientní metody (least_squares) konvergují k **nejbližšímu** minimu - výsledek závisí na initial guess.

### 1.2 Řešení: Multi-start

Jednoduchá myšlenka: **spustit optimalizaci vícekrát** z různých počátečních bodů a vybrat nejlepší výsledek.

```
Start 1: guess_1 -> local_min_1 (error 5.2%)
Start 2: guess_2 -> local_min_2 (error 3.8%)
Start 3: guess_3 -> local_min_3 (error 1.2%)  <- nejlepší
Start 4: guess_4 -> local_min_1 (error 5.2%)
...
```

### 1.3 Klíčová inovace: Adaptivní perturbace

Naivní přístup: náhodné perturbace všech parametrů stejně.

**Náš přístup:** Využíváme **kovarianční matici** z prvního fitu:
- Parametry s vysokou nejistotou -> větší perturbace
- Parametry s nízkou nejistotou -> menší perturbace
- Korelované parametry -> perturbovány společně

Toto je **informovanější** prohledávání prostoru.


## 2. Matematické pozadí

### 2.1 Kovarianční matice

Po fitu získáme kovarianční matici parametrů:

```
        [ σ₁²    σ₁₂   σ₁₃  ]
Cov =   [ σ₁₂   σ₂²    σ₂₃  ]
        [ σ₁₃   σ₂₃   σ₃²   ]
```

kde:
- σᵢ² = variance parametru i (diagonála)
- σᵢⱼ = kovariance mezi parametry i a j (mimo-diagonální)

**Interpretace:**
- Velká σᵢ² -> parametr je špatně určen daty -> má smysl ho perturbovat více
- Velká |σᵢⱼ| -> parametry jsou korelované -> měly by se perturbovat společně

### 2.2 Korelace parametrů v EIS

**Příklad: Voigt element R || C**

```
Z(ω) = R / (1 + jωRC)
```

Pro nízké frekvence (ω << 1/RC):
```
Z ≈ R - jωR²C
```

Imaginární část závisí na součinu R²C. Pokud zvýšíme R a snížíme C tak, že R²C zůstane konstantní, fit se téměř nezmění.

**Důsledek:** R a C jsou **antikorelované** (ρ < 0).

Pokud perturbujeme R a C nezávisle, většina perturbací bude špatná. Ale pokud je perturbujeme **společně** ve směru korelace, zachováme R²C a najdeme alternativní minima.

### 2.3 Cholesky dekompozice

Pro generování korelovaných perturbací používáme **Cholesky dekompozici**:

```
Cov = L @ L^T
```

kde L je dolní trojúhelníková matice.

**Algoritmus:**
```python
z = randn(n_params)        # Nezávislé N(0,1)
perturbation = scale * L @ z  # Korelované N(0, scale²*Cov)
params_new = params + perturbation
```

**Proč to funguje:**

Pokud z ~ N(0, I), pak:
```
E[L@z] = 0
Cov[L@z] = L @ I @ L^T = L @ L^T = Cov
```

Tedy L @ z má přesně tu kovarianci, kterou chceme.

### 2.4 Geometrická interpretace

Kovarianční matice definuje **elipsoid nejistoty** v parametrickém prostoru:

```
        θ₂
         ^
         |    .  .
         |  .      .    <- elipsoid 2σ
         | .   *    .   <- optimum
         |  .      .
         |    .  .
         +---------------> θ₁
```

- Směr hlavní osy elipsoidu = směr největší nejistoty
- Délka os = standardní odchylky ve směru os

Cholesky perturbace generuje body uvnitř tohoto elipsoidu, s větší hustotou blízko středu.


## 3. Strategie perturbace

### 3.1 Covariance-based (preferováno)

```python
def perturb_from_covariance(params, cov, scale=2.0, bounds=None):
    L = np.linalg.cholesky(cov)  # Cov = L @ L^T
    z = np.random.randn(n_params)
    perturbation = scale * L @ z
    return clip(params + perturbation, bounds)
```

**Výhody:**
- Zachovává korelace mezi parametry
- Perturbuje více tam, kde je větší nejistota
- Efektivnější než náhodné prohledávání

**Kdy se použije:**
- Kovarianční matice je k dispozici
- Condition number < 10¹⁰ (well-conditioned)

### 3.2 Stderr-based (fallback)

```python
def perturb_from_stderr(params, stderr, scale=2.0, bounds=None):
    perturbation = scale * stderr * np.random.randn(n_params)
    return clip(params + perturbation, bounds)
```

Jednodušší varianta - ignoruje korelace, ale stále perturbuje více tam, kde je větší nejistota.

**Kdy se použije:**
- Kovarianční matice je ill-conditioned
- Cholesky dekompozice selhala

### 3.3 Log-uniform (last resort)

```python
def perturb_log_uniform(params, factor=3.0, bounds=None):
    multiplier = exp(uniform(-log(factor), log(factor)))
    return clip(params * multiplier, bounds)
```

Multiplikativní perturbace - každý parametr je vynásoben náhodným faktorem v rozsahu [1/factor, factor].

**Kdy se použije:**
- Standard errors jsou nekonečné nebo nulové
- Žádná informace o nejistotě není k dispozici

**Proč log-uniform:**
EIS parametry typicky pokrývají mnoho řádů (R: 1-10⁶ Ω, C: 10⁻¹²-10⁻³ F). Aditivní perturbace by byla buď příliš malá pro velké parametry, nebo příliš velká pro malé. Multiplikativní perturbace je scale-invariant.


## 4. Algoritmus

### 4.1 Schéma

```
1. Initial fit
   - Spuštění least_squares z hodnot v definici obvodu
   - Získání params_opt, cov, stderr

2. Generování perturbací
   - Pro i = 2, ..., n_restarts:
     - Pokud cov je well-conditioned:
         params_i = perturb_from_covariance(params_opt, cov, scale)
     - Jinak pokud stderr jsou finite:
         params_i = perturb_from_stderr(params_opt, stderr, scale)
     - Jinak:
         params_i = perturb_log_uniform(params_opt, factor=3.0)

3. Opakované fity
   - Pro každou perturbaci:
     - Spuštění least_squares z params_i
     - Uložení výsledku (nebo "fail" při chybě)

4. Výběr nejlepšího
   - Vrácení výsledku s nejnižším fit_error_rel
```

### 4.2 Hierarchie strategií

```
                  +------------------+
                  | Máme kovarianční |
                  | matici?          |
                  +--------+---------+
                           |
              +------------+------------+
              | ANO                     | NE
              v                         v
    +------------------+      +------------------+
    | Je well-         |      | Máme konečné    |
    | conditioned?     |      | stderr?         |
    +--------+---------+      +--------+---------+
             |                         |
    +--------+--------+       +--------+--------+
    | ANO             | NE    | ANO             | NE
    v                 v       v                 v
+----------+    +----------+  +----------+   +----------+
| Cholesky |    | Stderr   |  | Stderr   |   | Log-     |
| perturb  |    | perturb  |  | perturb  |   | uniform  |
+----------+    +----------+  +----------+   +----------+
```


## 5. Intuitivní vysvětlení

### 5.1 Analogie: Hledání nejvyššího kopce v mlze

Představte si, že hledáte nejvyšší kopec v horách zahalených mlhou:

```
          ?
         /\     ?
        /  \   /\     ?
       /    \_/  \   /\
      /           \_/  \___
     /                      \
```

**Naivní přístup (náhodný restart):**
- Teleportujete se na náhodné místo
- Jdete nahoru, dokud nedosáhnete vrcholu
- Opakujete z jiného náhodného místa

**Adaptivní přístup (náš multistart):**
- Najdete jeden vrchol
- Prozkoumáte okolí: "Kde jsem si nejméně jistý pozicí?"
- Teleportujete se více tím směrem
- Pokud je tam vyšší kopec, najdete ho efektivněji

### 5.2 Proč kovariance pomáhá

Představte si fit Voigt elementu R-(R₁|C₁):

```
Parametr    Hodnota     Stderr      Interpretace
R₀          100 Ω       ±2 Ω        Dobře určeno (2%)
R₁          5000 Ω      ±500 Ω      Méně jisté (10%)
C₁          1e-6 F      ±5e-7 F     Nejisté (50%)
```

**Naivní perturbace (±10% všechno):**
- R₀: 100 ± 10 Ω -> moc velká změna (data říkají ±2)
- C₁: 1e-6 ± 1e-7 F -> moc malá změna (data říkají ±5e-7)

**Adaptivní perturbace (scale=2σ):**
- R₀: 100 ± 4 Ω -> přiměřená
- R₁: 5000 ± 1000 Ω -> přiměřená
- C₁: 1e-6 ± 1e-6 F -> přiměřená

Perturbujeme **tam, kde data dovolují** - kde je cost function plochá.

### 5.3 Korelace = společný pohyb

Pokud R₁ a C₁ jsou antikorelované (τ = R₁C₁ je dobře určeno):

```
        C₁
         ^
         |   x x x
         |  x     x    <- elipsoid nejistoty
         | x   *   x   <- optimum
         |  x     x
         |   x x x
         +---------------> R₁
```

Cholesky perturbace generuje body podél elipsy, ne v čtverci:
- Pokud R₁ roste, C₁ klesá (zachovává τ)
- Toto jsou směry, kde se fit mění málo
- Tam mohou být alternativní minima


## 6. Implementace

### 6.1 Python API

```python
from eis_analysis.fitting import fit_circuit_multistart, R, C, Q

# Definice obvodu s initial guess
circuit = R(100) - (R(5000) | C(1e-6))

# Multi-start optimalizace
result, Z_fit, fig = fit_circuit_multistart(
    circuit, freq, Z,
    n_restarts=10,      # Počet startů
    scale=2.0,          # Perturbace = 2 sigma
    weighting='sqrt',   # Vážení
    parallel=False,     # Sekvenční běh
    max_workers=4,      # Workers pro paralelní
    verbose=True,       # Logovat průběh
    use_analytic_jacobian=True  # Analytický Jacobián (rychlejší, přesnější)
)

# Výsledky
print(f"Best error: {result.best_result.fit_error_rel:.2f}%")
print(f"Improvement: {result.improvement:.1f}%")
print(f"Successful: {result.n_successful}/{result.n_starts}")
```

### 6.2 CLI použití

```bash
# Základní multi-start s 10 restarty
./eis.py --multistart 10 --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA

# S vlastní škálou perturbace
./eis.py --multistart 20 --multistart-scale 3.0 \
         --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

### 6.3 Výstup MultistartResult

```python
@dataclass
class MultistartResult:
    best_result: FitResult      # Nejlepší výsledek
    all_results: List[FitResult]  # Všechny úspěšné výsledky
    n_starts: int               # Počet startů
    n_successful: int           # Počet úspěšných
    improvement: float          # Zlepšení oproti initial [%]
```


## 7. Parametry a jejich vliv

### 7.1 Počet restartů (`--multistart`)

| n_restarts | Efekt | Vhodné pro |
|------------|-------|------------|
| 5 | Rychlý test | Jednoduché obvody |
| 10 | **Doporučeno** | Většina případů |
| 20-30 | Důkladnější | Komplexní obvody |
| 50+ | Velmi důkladné | Podezření na mnoho minim |

**Intuice:** Pravděpodobnost nalezení globálního minima roste s počtem restartů, ale s klesajícím výnosem.

### 7.2 Škála perturbace (`--multistart-scale`)

| scale | Efekt |
|-------|-------|
| 1.0 | Perturbace = 1σ, blízko originálu |
| 2.0 | **Default.** Perturbace = 2σ, dobrý kompromis |
| 3.0 | Širší explorace |
| 5.0+ | Velmi široká, může být nestabilní |

**Matematický význam:**
- scale=2.0 znamená, že 95% perturbací je v rozsahu ±2σ
- Větší scale = větší explorace = větší šance najít vzdálená minima
- Příliš velká scale = mnohé perturbace mimo rozumný rozsah

### 7.3 Vážení (`--weighting`)

| Typ | Efekt |
|-----|-------|
| uniform | Velké |Z| dominuje |
| sqrt | **Default.** Kompromis |
| proportional | Všechny frekvence stejně |

Vážení ovlivňuje jak fit, tak kovarianci (a tedy perturbace).

### 7.4 Analytický Jacobián (`use_analytic_jacobian`)

| Hodnota | Efekt |
|---------|-------|
| True | **Default.** Analytické derivace, rychlejší a přesnější |
| False | Numerická aproximace (finite differences) |

**Podporované elementy pro analytický Jacobián:**
- R, C, L (základní elementy)
- Q (Constant Phase Element)
- W (Warburg)
- Wo (Warburg open/bounded)
- K (Voigt element)

Pro nepodporované elementy systém automaticky přepne na numerický Jacobián.


## 8. Formát výstupu

### 8.1 Příklad výstupu

```
==================================================
Multi-start optimization
==================================================
  Restarts: 10, scale: 2.0 sigma
  Jacobian: analytic

Progress (error %):
   1-10: 1.40, 1.40, 1.35, 1.40, [1.32], 1.40, 1.38, 1.40, 1.40, 1.40
==================================================
Multi-start summary
==================================================
  Successful fits: 10/10
  Initial error:   1.401%
  Best error:      1.320% (start #5)
  Improvement:     +5.8%
==================================================
Best fit results
==================================================

Fit results:
  Parameters:
    R0    = 9.99e+01 +/- 1.02e+00  [95% CI: 9.76e+01, 1.02e+02]
    R1    = 5.00e+03 +/- 9.42e+00  [95% CI: 4.98e+03, 5.02e+03]
    C0    = 1.00e-05 +/- 5.05e-08  [95% CI: 9.91e-06, 1.01e-05]
  Fit error: 1.32% (rel), 28.72 Ohm (abs)
  Quality: Good (<10.0%)
==================================================
```

### 8.2 Interpretace

**Progress řádek:**
- Čísla = fit error [%] pro každý start
- `[1.32]` = nejlepší výsledek (zvýrazněno)
- `fail` = optimalizace selhala

**Summary:**
- `Initial error` = chyba prvního fitu (z initial guess)
- `Best error` = nejnižší nalezená chyba
- `Improvement` = relativní zlepšení

**Kdy je improvement významný:**
- < 1% - marginální, initial guess byl dobrý
- 1-10% - významné zlepšení
- > 10% - initial guess byl daleko od optima


## 9. Srovnání s jinými metodami

### 9.1 Multi-start vs Single fit

| Aspekt | Single fit | Multi-start |
|--------|-----------|-------------|
| Čas | 1× | n× |
| Robustnost | Závisí na guess | Výrazně vyšší |
| Nalezení globálního | ~60% | ~85-99% |

### 9.2 Multi-start vs Differential Evolution

| Aspekt | Multi-start | DE |
|--------|-------------|-----|
| Initial guess | Vyžaduje | Nevyžaduje |
| Rychlost | Rychlejší | Pomalejší |
| Explorace | Lokální okolo guess | Globální |
| Využití informace | Kovariance | Žádná |
| Evaluací | ~500-2000 | ~15000+ |

**Kdy multi-start:**
- Máte rozumný initial guess
- Potřebujete rychlý výsledek
- Obvod je jednodušší (2-4 parametry)

**Kdy DE:**
- Nemáte initial guess
- Obvod je komplexní (>5 parametrů)
- Podezření na mnoho lokálních minim

### 9.3 Kombinace metod

Pro maximální robustnost:

```bash
# 1. DE pro globální prohledání
./eis.py --optimizer de --circuit '...' data.DTA

# 2. Multi-start pro verifikaci
./eis.py --multistart 20 --circuit '...' data.DTA
```

Pokud obě metody dají stejný výsledek, je vysoká pravděpodobnost, že je to globální minimum.


## 10. Diagnostika a řešení problémů

### 10.1 Všechny restarty dávají stejný výsledek

```
Progress (error %):
   1-10: 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40, 1.40
```

**Možné příčiny:**
1. Initial guess je již u globálního minima (dobře!)
2. Perturbace jsou příliš malé
3. Cost function je konvexní (jediné minimum)

**Řešení:**
- Zvýšit scale: `--multistart-scale 4.0`
- Zkontrolovat stderr: pokud jsou velmi malé, perturbace budou malé

### 10.2 Mnoho "fail" výsledků

```
Progress (error %):
   1-10: 1.40, fail, fail, 1.35, fail, fail, fail, 1.40, fail, fail
```

**Možné příčiny:**
1. Perturbace mimo fyzikální rozsah
2. Scale příliš velká
3. Bounds příliš úzké

**Řešení:**
- Snížit scale: `--multistart-scale 1.5`
- Zkontrolovat bounds parametrů

### 10.3 Velký rozptyl výsledků

```
Progress (error %):
   1-10: 5.40, 12.30, 2.35, 8.40, 1.32, 15.40, 3.38, 9.40, 6.40, 4.40
```

**Interpretace:**
- Cost function má mnoho lokálních minim
- Initial guess pravděpodobně nebyl u globálního

**Akce:**
- Nejlepší výsledek (1.32%) je pravděpodobně správný
- Zvážit více restartů pro větší jistotu
- Zvážit použití DE pro systematičtější prohledání


## 11. Matematické detaily

### 11.1 Výpočet kovariance

Z fitu získáme Jacobián J a residuály r. Kovariance:

```
s² = ||r||² / (n - p)     # residual variance
Cov = s² * (J^T @ J)^(-1)  # kovariance parametrů
```

kde n = počet datových bodů, p = počet parametrů.

### 11.2 Regularizace Cholesky

Pro numerickou stabilitu přidáváme malou regularizaci:

```python
cov_reg = cov + 1e-10 * np.eye(n_params)
L = np.linalg.cholesky(cov_reg)
```

Toto zajistí, že matice je pozitivně definitní i při numerických chybách.

### 11.3 Bounds clipping

Po perturbaci zajistíme, že parametry jsou v rozumném rozsahu:

```python
perturbed = np.clip(perturbed, lower_bounds, upper_bounds)
perturbed = np.maximum(perturbed, 1e-15)  # Zajistit pozitivní
```

### 11.4 Statistika úspěšnosti

Pro n restartů s pravděpodobností p nalezení globálního minima:

```
P(najdeme globální aspoň jednou) = 1 - (1-p)^n
```

| p | n=5 | n=10 | n=20 |
|---|-----|------|------|
| 0.3 | 83% | 97% | 99.9% |
| 0.5 | 97% | 99.9% | ~100% |
| 0.7 | 99.8% | ~100% | ~100% |


## 12. Příklady použití

### 12.1 Jednoduchý RC obvod

```bash
./eis.py --multistart 10 --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

Pro 3 parametry obvykle stačí 10 restartů.

### 12.2 Komplexní obvod s Q

```bash
./eis.py --multistart 20 --multistart-scale 2.5 \
         --circuit 'R(1)-(R(1000)|Q(1e-5,0.9))-(R(10000)|Q(1e-4,0.8))' \
         data.DTA
```

Více parametrů = více restartů + mírně větší scale.

### 12.3 Python API s analýzou

```python
from eis_analysis.fitting import fit_circuit_multistart, R, C

circuit = R(100) - (R(5000) | C(1e-6))

result, Z_fit, fig = fit_circuit_multistart(
    circuit, freq, Z,
    n_restarts=20,
    scale=2.0,
    verbose=True
)

# Analýza rozptylu výsledků
errors = [r.fit_error_rel for r in result.all_results]
print(f"Error range: {min(errors):.2f}% - {max(errors):.2f}%")
print(f"Error std: {np.std(errors):.2f}%")

# Pokud je rozptyl velký, máme mnoho lokálních minim
if np.std(errors) > 1.0:
    print("Warning: High variance suggests multiple local minima")
```


## 13. Závěr

### 13.1 Silné stránky

1. **Adaptivní perturbace** - Využívá kovarianci pro inteligentní prohledávání
2. **Zachování korelací** - Cholesky dekompozice respektuje strukturu problému
3. **Rychlost** - Méně evaluací než globální metody
4. **Robustnost** - Hierarchie fallback strategií
5. **Diagnostika** - Přehledný výstup s improvement statistikou

### 13.2 Slabé stránky

1. **Vyžaduje initial guess** - Na rozdíl od DE
2. **Lokální prohledávání** - Perturbace jsou kolem initial, ne globální
3. **Závislost na prvním fitu** - Pokud je špatný, perturbace budou špatné

### 13.3 Doporučení

| Situace | Doporučené nastavení |
|---------|---------------------|
| Běžný fit | `--multistart 10` |
| Komplexní obvod | `--multistart 20 --multistart-scale 2.5` |
| Nejistý initial guess | `--multistart 30 --multistart-scale 3.0` |
| Verifikace DE výsledku | `--multistart 10` (s params z DE) |

### 13.4 Best practices

1. **Začněte s rozumným initial guess** - Multi-start vylepšuje, ne nahrazuje
2. **Sledujte improvement** - Pokud je < 1%, guess byl dobrý
3. **Kontrolujte rozptyl** - Velký rozptyl = mnoho minim = zvažte DE
4. **Kombinujte s DE** - Pro maximální jistotu použijte obě metody


---

*Dokument vytvořen pro EIS Analysis Toolkit*
*Autor: Claude Code*
*Poslední aktualizace: 2026-01-01*
