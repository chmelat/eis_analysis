# Matematická analýza nelineárního fitování EIS obvodů

## 1. Formulace problému

### 1.1 Optimalizační problém

Fitování ekvivalentního obvodu na EIS data je nelineární nejmenších čtverců:

```
min_θ ||r(θ)||² = min_θ Σᵢ wᵢ² |Z_data(ωᵢ) - Z_model(ωᵢ, θ)|²
```

kde:
- θ = [θ₁, θ₂, ..., θₚ] jsou parametry obvodu (R, C, Q, n, τ, ...)
- wᵢ jsou váhy (uniform, sqrt, proportional)
- Z_model je impedanční funkce obvodu

### 1.2 Architektura modulu `eis_analysis.fitting`

Modulární struktura:

| Soubor | Funkce |
|--------|--------|
| `circuit.py` | Hlavní fit funkce, FitResult dataclass |
| `circuit_elements.py` | Definice elementů (R, C, Q, W, ...) |
| `circuit_builder.py` | Operátory pro skládání obvodů |
| `bounds.py` | Fyzikálně rozumné limity parametrů |
| `covariance.py` | SVD-based výpočet kovariance |
| `diagnostics.py` | Diagnostika kvality fitu, vážení |
| `jacobian.py` | Analytické derivace pro Jacobián |
| `multistart.py` | Multi-start optimalizace |
| `diffevo.py` | Diferenciální evoluce (globální optimalizace) |

### 1.3 Dostupné optimalizační metody

| Metoda | Funkce | Typ | Použití |
|--------|--------|-----|---------|
| Single fit | `fit_equivalent_circuit()` | Lokální | Rychlý fit s dobrým initial guess |
| Multi-start | `fit_circuit_multistart()` | Semi-globální | Více lokálních fitů s perturbacemi |
| Differential Evolution | `fit_circuit_diffevo()` | Globální | Systematické prohledání prostoru |

**Detailní dokumentace:**
- Multi-start: viz `MULTISTART_OPTIMIZATION.md`
- Differential Evolution: viz `DIFFERENTIAL_EVOLUTION.md`


## 2. Matematické problémy

### 2.1 Nekonvexnost účelové funkce

Impedanční funkce jsou **vysoce nelineární** a účelová funkce má **mnohočetná lokální minima**.

**Příklad: Voigt element (R || C)**

```
Z(ω) = R / (1 + jωRC)
```

Pro parametry θ = [R, C]:
- Gradient: ∂Z/∂R = 1/(1 + jωRC) - jωRC/(1 + jωRC)²
- Hessian má složitou strukturu

**Vizualizace problému:**

```
RSS (θ)
   ^
   |     /\
   |    /  \     /\
   |   /    \   /  \    <- lokální minima
   |  /      \_/    \
   | /                \
   |/                  \___  <- globální minimum
   +-----------------------> θ
```

**Řešení v implementaci:**
- Multi-start: Opakované lokální fity z perturbovaných startů
- Differential Evolution: Populační globální optimalizace

### 2.2 Špatná podmíněnost (Ill-conditioning)

#### 2.2.1 Rozdílné škály parametrů

Typický EIS obvod má parametry s **extrémně rozdílnými řády**:

| Parametr | Typická hodnota | Řád |
|----------|-----------------|-----|
| R (odpor) | 1-10⁶ Ω | 10⁰-10⁶ |
| C (kapacita) | 10⁻¹²-10⁻³ F | 10⁻¹²-10⁻³ |
| Q (Q) | 10⁻⁶-10⁻² | 10⁻⁶-10⁻² |
| n (exponent) | 0.5-1.0 | 10⁰ |
| τ (časová konstanta) | 10⁻⁶-10² s | 10⁻⁶-10² |
| L (indukčnost) | 10⁻⁹-10⁻⁶ H | 10⁻⁹-10⁻⁶ |

**Důsledek:** Jakobián J má sloupce s velmi rozdílnými normami.

Condition number κ(JᵀJ) může být 10¹⁰-10²⁰.

**Řešení v implementaci:**
```python
x_scale = np.maximum(np.abs(initial_guess), 1e-10)
```

#### 2.2.2 Korelované parametry

Některé parametry jsou silně korelované:

**Příklad 1: R a C ve Voigt elementu**

Pro nízké frekvence (ω << 1/RC):
```
Z ≈ R - jωR²C
```

Změna R a kompenzační změna C dává podobnou impedanci -> silná korelace.

**Příklad 2: Q a n v Q**

```
Z_Q = 1 / (Q·(jω)ⁿ)
```

Zvýšení Q a snížení n může dát podobný výsledek -> korelace.

**Důsledek:**
- Kovarianční matice je téměř singulární
- Standardní chyby parametrů jsou velké (nebo nekonečné)
- Více kombinací parametrů dává podobně dobrý fit

### 2.3 Singulární body a asymptoty

Některé impedanční funkce mají **singularity**:

| Element | Funkce | Singularita |
|---------|--------|-------------|
| Kapacitor | Z_C = 1/(jωC) | C -> 0: \|Z\| -> ∞ |
| Q | Z_Q = 1/(Q·(jω)ⁿ) | n -> 0: Z -> 1/Q |
| Warburg | Z_W = σ(1-j)/√ω | ω -> 0: \|Z\| -> ∞ |

**Důsledek:** Gradient může "explodovat" blízko singularit.

### 2.4 Bounds a constraints

**Fyzikální omezení (`bounds.py`):**

```python
PARAMETER_BOUNDS = {
    'R': (1e-4, 1e10),    # 0.1 mOhm - 10 GOhm
    'C': (1e-15, 1e-1),   # 1 fF - 100 mF
    'L': (1e-12, 1e-4),   # 1 pH - 100 uH (parasitní)
    'Q': (1e-12, 1e-1),   # Q koeficient
    'n': (0.4, 1.0),      # Q exponent
    'σ': (1e-2, 1e5),     # Warburg
    'τ': (1e-9, 1e4),     # 1 ns - 10000 s
    'R_W': (1e-2, 1e8),   # Warburg bounded - odpor
    'τ_W': (1e-6, 1e4),   # Warburg bounded - difuzní čas
}
```

**Výhody tight bounds:**
- Prevence konvergence k nefyzikálním hodnotám
- Zlepšení numerické stability
- Pomoc při špatném initial guess

**Automatická kontrola blízkosti k bounds:**
Systém varuje, pokud parametr konverguje blízko hranice.


## 3. Analýza Jakobiánu

### 3.1 Struktura Jakobiánu

Pro komplexní impedanci s N frekvenčními body a P parametry:

```
J ∈ ℝ^(2N × P)

J = [ ∂Re(Z)/∂θ₁  ∂Re(Z)/∂θ₂  ...  ∂Re(Z)/∂θₚ ]
    [ ∂Im(Z)/∂θ₁  ∂Im(Z)/∂θ₂  ...  ∂Im(Z)/∂θₚ ]
```

### 3.2 Analytické derivace pro základní elementy

**Resistor R:**
```
Z = R
∂Z/∂R = 1
```

**Kapacitor C:**
```
Z = 1/(jωC) = -j/(ωC)
∂Z/∂C = j/(ωC²) = -Z/C
```

**Q (Q, n):**
```
Z = 1/(Q·(jω)ⁿ)
∂Z/∂Q = -Z/Q
∂Z/∂n = -Z·ln(jω)
```

**Voigt K(R, τ):**
```
Z = R/(1 + jωτ)
∂Z/∂R = 1/(1 + jωτ) = Z/R
∂Z/∂τ = -jωR/(1 + jωτ)² = -jωZ²/R
```

### 3.3 Analytické vs numerické derivace

Implementace podporuje **analytické derivace** (modul `jacobian.py`) s fallback na numerické.

| Aspekt | Analytické | Numerické |
|--------|------------|-----------|
| Implementace | Explicitní vzorce | Automatická (scipy) |
| Rychlost | 1 evaluace | P+1 evaluací |
| Přesnost | Exaktní | O(h) nebo O(h²) |
| Stabilita | Robustní | Citlivé na h |

**Podporované elementy pro analytický Jacobián:**
- R, C, L (základní elementy)
- Q (Constant Phase Element) - derivace podle Q a n
- W (Warburg)
- Wo (Warburg open/bounded)
- K (Voigt element) - derivace podle R a τ

**Aktivace:**
```python
result, Z_fit, fig = fit_equivalent_circuit(
    freq, Z, circuit,
    use_analytic_jacobian=True  # Default
)
```

Pro nepodporované custom elementy systém automaticky přepne na numerický Jacobián.


## 4. Řešení problému initial guess

### 4.1 Citlivost na počáteční odhad

Nelineární optimalizace je **silně závislá** na initial guess (pro lokální metody).

**Příklad:**
Pro Randles obvod R₀ - (R₁ || C₁):
- Dobrý guess [100, 5000, 1e-6] -> konverguje k [98, 4823, 8.7e-7]
- Špatný guess [1, 1, 1] -> konverguje k lokálnímu minimu

### 4.2 Implementované přístupy

| Metoda | Initial guess | Jak řeší problém |
|--------|---------------|------------------|
| Single fit | Vyžaduje dobrý | Neřeší |
| Multi-start | Vyžaduje přibližný | Perturbace kolem initial |
| Differential Evolution | Využívá jako x0 seed | Jeden člen populace, zbytek náhodně |

**Multi-start:** Perturbuje initial guess pomocí kovariance z prvního fitu.

**Differential Evolution:** Inicializuje jednoho člena populace hodnotami z obvodu, zbytek náhodně v bounds. Tím kombinuje informovaný start s globální explorací.


## 5. Identifikovatelnost parametrů

### 5.1 Strukturální identifikovatelnost

Některé parametrizace jsou **strukturálně neidentifikovatelné**:

**Příklad: Dva Voigt elementy v sérii**
```
Z = R₁/(1+jωτ₁) + R₂/(1+jωτ₂)
```

Pokud τ₁ ≈ τ₂, pak [R₁, τ₁, R₂, τ₂] a [R₂, τ₂, R₁, τ₁] jsou ekvivalentní.

### 5.2 Diagnostika (`covariance.py`, `diagnostics.py`)

| Metrika | Threshold | Indikuje |
|---------|-----------|----------|
| Condition number | > 10¹⁰ | Ill-conditioned |
| Korelace \|ρᵢⱼ\| | > 0.95 | Silná korelace |
| Relativní stderr | > 100% | Neidentifikovatelný parametr |
| Numerický rank | < p | Rank-deficientní |

### 5.3 SVD-based výpočet kovariance

```python
# SVD dekompozice: J = U @ S @ V^T
U, S, Vt = np.linalg.svd(jacobian, full_matrices=False)

# Condition number
condition_number = S[0] / S[-1]

# Regularizovaná inverze
S_inv_sq = np.where(S > threshold, 1/(S*S), 1/(threshold*threshold))
JtJ_inv = V @ np.diag(S_inv_sq) @ Vt

# Kovariance
cov = residual_variance * JtJ_inv
```

**Výhody SVD přístupu:**
- Numericky stabilní i pro ill-conditioned problémy
- Automatická regularizace malých singulárních hodnot
- Detekce rank-deficientních problémů

### 5.4 Confidence intervals

Používá t-distribuci s (n_data - n_params) stupni volnosti:

```python
t_critical = t.ppf(1 - alpha/2, dof)
CI = params ± t_critical * stderr
```


## 6. Vážení dat

### 6.1 Typy vážení

| Typ | Váha wᵢ | Efekt |
|-----|---------|-------|
| uniform | 1 | Velké \|Z\| dominuje |
| sqrt | 1/√\|Z\| | Kompromis (doporučeno) |
| proportional | 1/\|Z\| | Všechny frekvence stejně |
| square | \|Z\|² | Extrémně velké \|Z\| dominuje |

**Poznámka:** Váhy jsou normalizovány tak, aby měly průměr 1.

### 6.2 Vliv na fit

**Uniform vážení:**
- Nízké frekvence (\|Z\| ~ kΩ) dominují
- Vysokofrekvenční detaily (\|Z\| ~ Ω) jsou ignorovány

**Proportional vážení:**
- Všechny frekvence přispívají podobně
- Lepší pro široký dynamický rozsah

### 6.3 Vliv na kovarianční matici

Kovariance závisí na vážení:
```
cov(θ) = σ² (JᵀWWJ)⁻¹
```

kde W = diag(w₁, w₂, ..., wₙ).


## 7. Přehled optimalizačních metod

### 7.1 Single fit (`fit_equivalent_circuit`)

```python
from eis_analysis.fitting import fit_equivalent_circuit, R, C

circuit = R(100) - (R(5000) | C(1e-6))
result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
```

**Algoritmus:** scipy.optimize.least_squares (Trust Region Reflective)

**Kdy použít:** Máte dobrý initial guess, rychlý výsledek.

### 7.2 Multi-start (`fit_circuit_multistart`)

```python
from eis_analysis.fitting import fit_circuit_multistart, R, C

circuit = R(100) - (R(5000) | C(1e-6))
result, Z_fit, fig = fit_circuit_multistart(
    circuit, freq, Z,
    n_restarts=10,
    scale=2.0
)
```

**Algoritmus:**
1. Initial fit z hodnot v obvodu
2. Generování perturbací pomocí kovariance (Cholesky)
3. Opakované lokální fity
4. Výběr nejlepšího

**Kdy použít:** Přibližný initial guess, střední složitost obvodu.

**Detaily:** viz `MULTISTART_OPTIMIZATION.md`

### 7.3 Differential Evolution (`fit_circuit_diffevo`)

```python
from eis_analysis.fitting import fit_circuit_diffevo, R, C

circuit = R(100) - (R(5000) | C(1e-6))
result, Z_fit, fig = fit_circuit_diffevo(
    circuit, freq, Z,
    strategy=1,      # randtobest1bin
    popsize=15,
    maxiter=1000
)
```

**Algoritmus:**
1. Inicializace populace (x0 = initial guess, zbytek náhodně)
2. Evoluční cyklus (mutace, křížení, selekce)
3. least_squares refinement pro přesné parametry
4. Výpočet kovariance

**Kdy použít:** Nejistý initial guess, komplexní obvod, mnoho parametrů.

**Detaily:** viz `DIFFERENTIAL_EVOLUTION.md`

### 7.4 Srovnání metod

| Aspekt | Single fit | Multi-start | DE |
|--------|-----------|-------------|-----|
| Initial guess | Vyžaduje dobrý | Vyžaduje přibližný | Využívá jako x0 |
| Rychlost | Nejrychlejší | Střední | Nejpomalejší |
| Robustnost | Nízká | Střední | Vysoká |
| Globální minimum | ~60% | ~85-95% | ~99% |
| Evaluací | ~50-200 | ~500-2000 | ~15000+ |


## 8. Lineární přístup (Voigt chain)

### 8.1 Speciální případ

Pro Voigt řetězec s fixovanými τᵢ:

```
Z(ω) = R_s + jωL + Σᵢ Rᵢ/(1 + jωτᵢ)
```

Pokud τᵢ jsou fixované, problém je **lineární v [R_s, R₁, R₂, ..., L]**.

### 8.2 Výhody a nevýhody

**Výhody:**
- Globální optimum garantováno
- Rychlé (NNLS nebo SVD)
- Žádný initial guess potřeba

**Nevýhody:**
- Vyžaduje znalost τᵢ předem (např. z DRT)
- Neplatí pro obecné obvody (Q, Warburg, atd.)

### 8.3 Hybridní přístup

Možná strategie:
1. DRT analýza pro identifikaci τᵢ
2. Lineární fit pro odhad Rᵢ
3. Nelineární fit pro doladění a přidání Q/Warburg


## 9. Implementované funkce

### 9.1 Core funkce

| Funkce | Modul | Popis |
|--------|-------|-------|
| `fit_equivalent_circuit` | circuit.py | Základní nelineární fit |
| `fit_circuit_multistart` | multistart.py | Multi-start optimalizace |
| `fit_circuit_diffevo` | diffevo.py | Diferenciální evoluce |
| `compute_covariance_matrix` | covariance.py | SVD-based kovariance |
| `generate_simple_bounds` | bounds.py | Fyzikální bounds |
| `compute_weights` | diagnostics.py | Váhové schémata |
| `compute_fit_metrics` | diagnostics.py | Error metriky |
| `make_jacobian_function` | jacobian.py | Analytický Jacobián pro least_squares |

### 9.2 Dataclasses

**FitResult:**
```python
@dataclass
class FitResult:
    circuit: Circuit
    params_opt: NDArray[np.float64]
    params_stderr: NDArray[np.float64]
    fit_error_rel: float
    fit_error_abs: float
    quality: str
    condition_number: float
    is_well_conditioned: bool
    cov: Optional[NDArray[np.float64]]
```

**MultistartResult:**
```python
@dataclass
class MultistartResult:
    best_result: FitResult
    all_results: List[FitResult]
    n_starts: int
    n_successful: int
    improvement: float
```

**DiffEvoResult:**
```python
@dataclass
class DiffEvoResult:
    best_result: FitResult
    de_result: OptimizeResult
    de_error: float
    final_error: float
    n_evaluations: int
    strategy: str
    improvement: float
```


## 10. CLI použití

### 10.1 Single fit

```bash
./eis.py --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

### 10.2 Multi-start

```bash
./eis.py --multistart 10 --multistart-scale 2.0 \
         --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

### 10.3 Differential Evolution

```bash
./eis.py --optimizer de \
         --de-strategy 1 --de-popsize 15 --de-maxiter 1000 \
         --circuit 'R(100)-(R(5000)|C(1e-6))' data.DTA
```

### 10.4 Přehled CLI parametrů

| Parametr | Default | Popis |
|----------|---------|-------|
| `--circuit` | - | Definice obvodu |
| `--weighting` | sqrt | Váhové schéma |
| `--multistart` | 0 | Počet restartů (0 = disabled) |
| `--multistart-scale` | 2.0 | Perturbace v sigma |
| `--optimizer` | multistart | Optimizér (multistart/de) |
| `--de-strategy` | 1 | DE strategie (1/2/3) |
| `--de-popsize` | 15 | DE populace |
| `--de-maxiter` | 1000 | DE max iterací |
| `--de-tol` | 0.01 | DE tolerance |
| `--de-workers` | 1 | DE paralelizace |


## 11. Závěr

### 11.1 Současný stav implementace

| Aspekt | Stav | Hodnocení |
|--------|------|-----------|
| Lokální optimalizace | least_squares (TRF) | Standard |
| Multi-start | Adaptive, covariance-based | Robustní |
| Globální optimalizace | Differential Evolution | Implementováno |
| Derivace | Analytické (s numeric fallback) | Robustní |
| Bounds | Fyzikálně rozumné | Dobré |
| Covariance | SVD-based, regularizované | Robustní |
| Diagnostika | Kompletní | Dobré |

### 11.2 Silné stránky

1. **Tři úrovně optimalizace** - Single fit, Multi-start, DE
2. **Robustní SVD-based kovariance** - Stabilní i pro ill-conditioned
3. **Adaptivní multi-start** - Využívá kovarianci pro perturbace
4. **DE s x0 seeding** - Kombinuje informovaný start s globální explorací
5. **Flexibilní API** - Operator overloading pro definici obvodů
6. **Fyzikální bounds** - Prevence nefyzikálních výsledků
7. **Kompletní diagnostika** - Condition number, korelace, CI

### 11.3 Doporučení pro použití

| Situace | Doporučená metoda |
|---------|-------------------|
| Dobrý initial guess, rychlý výsledek | Single fit |
| Přibližný guess, střední složitost | Multi-start (10-20) |
| Nejistý guess, komplexní obvod | DE |
| Maximální spolehlivost | DE + Multi-start verifikace |

### 11.4 Možná budoucí rozšíření

1. **Automatický initial guess** - Odhad z charakteristik dat
2. **Bayesiánský přístup** - MCMC pro posterior distribuce
3. **Model selection** - Automatická volba obvodu (AIC, BIC)
4. **Další elementy pro analytický Jacobián** - Rozšíření podpory


---

*Dokument vytvořen pro EIS Analysis Toolkit*
*Autor: Claude Code*
*Poslední aktualizace: 2026-01-01*
