# Automatický výběr regularizačního parametru λ (GCV + L-curve)

## Přehled

DRT analýza řeší **ill-posed** inverzní problém regularizací (Tikhonov 2. řádu)
s parametrem λ, který určuje kompromis mezi věrností datům a hladkostí výsledné
distribuce γ(τ). Tato implementace volí λ **automaticky a datově řízeně**, takže
odpadá subjektivní ruční tuning.

**Škála λ.** Minimalizuje se

```
(1/n) ‖Aγ − b‖²  +  λ ∫ (d²γ / d ln τ²)² d ln τ        (n = počet reziduí)
```

Obě části jsou v Ω², takže λ je bezrozměrná a nezávisí na hustotě τ-mřížky
(`--n-tau`), šířce frekvenčního okna, počtu bodů ani měřítku impedance - jen
na poměru šumu a struktury dat. Matice `L` je proto druhá diference dělená Δ²
(Δ = krok v ln τ), násobená √Δ (integrál) a √n (průměr rezidua). Holá
diference `[1, -2, 1]` (do v0.46) dávala λ úměrnou Δ³: automatická λ pak
rostla asi 10× s každým zdvojením mřížky. Naměřené automatické λ leží mezi
1e-9 a 1e-3 (syntetika s šumem 0.1 / 0.5 / 2 %: 1.6e-7 / 7.7e-6 / 2.2e-4);
rozsah hledání je `DRT_LAMBDA_RANGE = (1e-10, 1e-2)`, pevná výchozí hodnota
`DRT_LAMBDA_DEFAULT = 1e-6` (obojí `fitting/config.py`).

Automatický výběr je **výchozí chování** — spustí se vždy, když uživatel nezadá
λ ručně přes `--lambda`. Kombinuje dvě metody:

| Metoda | Role | Funkce v `drt/gcv.py` |
|--------|------|------------------------|
| **GCV** (Generalized Cross-Validation) | rychlý počáteční odhad | `compute_gcv_score`, `find_optimal_lambda_gcv` |
| **L-curve** | robustní korekce pro NNLS | `compute_lcurve_point`, `find_lcurve_corner` |
| **Hybrid** (GCV → L-curve) | **reálně používaný** výběr | `find_optimal_lambda_hybrid` |

`calculate_drt(..., auto_lambda=True)` volá `find_optimal_lambda_hybrid`;
při jeho selhání se použije čisté GCV jako fallback (`drt/core.py`).

## Proč hybrid, a ne čisté GCV?

GCV je odvozeno pro **lineární** least-squares. DRT ale řeší **NNLS**
(non-negative least squares, γ ≥ 0), což je nelineární constraint. Kvůli němu
má GCV pro DRT systematickou tendenci **podhodnocovat λ** (volí příliš malou
regularizaci → zašuměná γ s falešnými píky).

L-curve (graf log‖Ax−b‖ vs. log‖Lx‖) je vůči NNLS robustnější — „roh" křivky
odpovídá dobrému kompromisu nezávisle na linearitě. Je ale dražší na celém
rozsahu, proto se používá hybrid:

1. **GCV** dá rychlý hrubý odhad λ_gcv.
2. **L-curve** se prohledá jen v úzkém okolí λ_gcv (±1.5 dekády).
3. **Rozhodnutí** podle vzájemného poměru obou odhadů (viz níže).

## Použití

### CLI

Vstupní bod je `eis.py` (spouštěj `python3`).

```bash
# Automatický výběr λ (default — nic se nezadává)
python3 eis.py data.DTA

# Ruční λ (přebije automatiku)
python3 eis.py data.DTA --lambda 1e-6      # nebo -l 1e-6

# Vyšší rozlišení tau gridu
python3 eis.py data.DTA --n-tau 150        # nebo -n 150

# Bez vstupního souboru = demo se syntetickými daty (též auto λ)
python3 eis.py
```

### Relevantní přepínače

| Přepínač | Popis | Default |
|----------|-------|---------|
| `--lambda` / `-l` | Ruční λ. **Bez něj** se použije automatický výběr (GCV + L-curve). | None (= auto) |
| `--n-tau` / `-n` | Počet tau bodů (rozlišení DRT) | 100 |

> Žádný samostatný „zapínač" automatiky neexistuje — automatika je default a
> vypíná se právě a jen zadáním `--lambda`.

### Python API

```python
from eis_analysis.drt import calculate_drt
import numpy as np

freq = np.logspace(-2, 5, 61)
Z = ...  # naměřená impedance (complex)

# Automatický výběr λ (hybrid GCV + L-curve)
result = calculate_drt(freq, Z, n_tau=100, auto_lambda=True)

# Ruční λ
result = calculate_drt(freq, Z, n_tau=100, lambda_reg=1e-6)
```

## Implementační detaily

### 1. `compute_gcv_score(lambda_val, A, b, L)`

Vypočítá GCV score pro danou λ.

```
GCV(λ) = n · ‖b − A·x(λ)‖² / (n − trace K(λ))²
```

kde `K(λ) = A · (AᵀA + λ·LᵀL)⁻¹ · Aᵀ` je influence (projekční) matice a
`n − trace K = trace(I − K)` je efektivní počet stupňů volnosti.

Postup:
1. Vyřeš NNLS na rozšířeném systému `[A; √λ·L] x = [b; 0]`.
2. Spočítej **neregularizované** reziduum `r = b − A·x` a `‖r‖²`.
3. Spočítej `trace K` **přímou inverzí** (ne SVD):
   ```python
   M = A.T @ A + lambda_val * L.T @ L
   K = A @ np.linalg.solve(M, A.T)   # K = A·M⁻¹·Aᵀ
   trace_K = np.trace(K)
   ```
4. `GCV = n · ‖r‖² / (n − trace_K)²`.

Robustní error handling: při selhání NNLS, singulární M nebo
`n − trace_K ≤ 10⁻¹⁰` vrací `inf`.

> **Pozn.:** Pro NNLS je GCV pouze **aproximace** — vzorec předpokládá lineární
> řešení, kdežto `x ≥ 0` je nelineární constraint. Proto hybrid s L-curve.

### 2. L-curve

- **`compute_lcurve_point(λ, A, b, L)`** → `(log₁₀‖Ax−b‖, log₁₀‖Lx‖, x)`.
  Reziduum vs. regularizační norma pro dané λ (NNLS řešení).
- **`compute_lcurve_curvature(ρ, η)`** → křivost v log-log prostoru
  `κ = (ρ′η″ − ρ″η′) / (ρ′² + η′²)^{3/2}`, derivace přes `np.gradient`.
- **`find_lcurve_corner(λ, ρ, η)`** → roh = bod **maximální (kladné)** křivosti
  na vnitřku rozsahu.

  **Znaménková konvence:** L-křivka procházená rostoucím λ jde z malého rezidua /
  velké normy řešení do velkého rezidua / malé normy (ρ roste, η klesá). Roh je
  konvexní směrem k počátku — levotočivá (CCW) zatáčka s **kladnou** křivostí,
  proto `argmax`, ne `argmin`. Tuto konvenci hlídá `tests/test_lcurve_corner.py`.

### 3. `find_optimal_lambda_gcv(A, b, L, lambda_range=DRT_LAMBDA_RANGE, n_search=20)`

Čistě GCV, dvoufázové prohledání (fallback hybridu):

- **Fáze 1 (hrubá):** `n_search` bodů v log-prostoru `[1e-10, 1e-2]`, najdi
  minimum GCV.
- **Fáze 2 (jemná):** `n_search` bodů v okolí minima. Při minimu na okraji se
  rozsah rozšíří o dekádu ven.

Celkem 40 evaluací (2 × 20). Vrací `(λ_optimal, gcv_optimal)`.

### 4. `find_optimal_lambda_hybrid(..., lcurve_decades=1.5)` — reálně používaný

1. **GCV initial guess:** minimum GCV přes `n_search` bodů → λ_gcv.
2. **L-curve korekce:** L-křivka v rozsahu ±`lcurve_decades` (1.5) dekády kolem
   λ_gcv; najdi roh → λ_lcurve. Pokud roh padne na okraj okna, nastaví se
   `corner_at_edge` (varování — skutečný roh může ležet mimo úzké okno).
3. **Rozhodnutí** podle `ratio = λ_lcurve / λ_gcv`:

   | ratio | Volba | `method_used` |
   |-------|-------|---------------|
   | 0.1 < ratio < 10 | λ_lcurve (konzistentní) | `lcurve` |
   | ratio ≥ 10 | λ_lcurve (GCV podhodnotilo kvůli NNLS) | `lcurve_correction` |
   | ratio ≤ 0.1 | √(λ_gcv·λ_lcurve) (geom. průměr, nekonzistentní) | `geometric_mean` |

Vrací `(λ_optimal, gcv_score, diagnostics)`; `diagnostics` obsahuje `lambda_gcv`,
`lambda_lcurve`, `method_used`, `curvature`, `rho`, `eta`, `corner_at_edge`.

### 5. Integrace v `drt/core.py`

`_select_lambda(A, b, L, lambda_reg, auto_lambda)` vrací `LambdaSelection`:

```python
@dataclass
class LambdaSelection:
    lambda_value: float
    method: str            # 'user' | 'default' | 'gcv' | 'hybrid' | 'fallback'
    lambda_gcv: Optional[float] = None   # jen při L-curve korekci
    gcv_score: Optional[float] = None
    corner_at_edge: bool = False         # roh L-křivky na okraji okna (F7)
    lambda_at_edge: bool = False         # zvolené λ na mezi GCV rozsahu (F3/F7)
```

Když `auto_lambda=True`, volá se `find_optimal_lambda_hybrid`. Reportovaná
`method` je **`hybrid`** pouze pokud L-curve výrazně korigovala
(`method_used == 'lcurve_correction'`), jinak **`gcv`** (L-curve souhlasila s
GCV). `lambda_gcv` se ukládá/zobrazuje jen pro `hybrid`.

Detekce okrajů (náprava F3/F7): pokud λ_opt nebo λ_gcv narazí na mez rozsahu
`[1e-10, 1e-2]`, nebo je roh na okraji okna, nastaví se `lambda_at_edge` —
signál, že optimizér chce extrémnější regularizaci, než rozsah dovoluje
(typicky problém s daty / modelem). `calculate_drt(..., auto_lambda=...)` celý
výběr orchestruje (`drt/core.py`).

## Matematické pozadí

### Odvození GCV

Standardní leave-one-out cross-validation `CV(λ) = (1/n) Σᵢ [bᵢ − x̂₋ᵢ(λ)]²`
vyžaduje n řešení (pomalé). GCV ho aproximuje pomocí stopy influence matice:

```
GCV(λ) = n · ‖b − A·x(λ)‖² / (n − trace K)²
```

s předpokladem `[bᵢ − x̂₋ᵢ]² ≈ [bᵢ − x̂ᵢ]² / (1 − Kᵢᵢ)²` a uniformní `Kᵢᵢ ≈
trace(K)/n`. Výhoda: jediné řešení místo n.

### Influence matice K(λ)

```
K(λ) = A · (AᵀA + λ·LᵀL)⁻¹ · Aᵀ
```

- `Kᵢᵢ` měří vliv i-tého bodu na vlastní predikci, `trace K` = efektivní dimenze
  modelu.
- λ → 0: `trace K → rank(A)` (bez regularizace); λ → ∞: `trace K → 0`
  (plná regularizace).

Stopa se počítá **přímou inverzí** přes `np.linalg.solve(M, A.T)` (kód i tato
dokumentace; dřívější verze chybně uváděla SVD).

### L-curve a křivost

L-křivka `(ρ, η) = (log‖Ax−b‖, log‖Lx‖)` má v log-log prostoru tvar „L"; roh
(maximum kladné křivosti) odpovídá optimálnímu kompromisu reziduum/hladkost.
Křivost se počítá centrálními diferencemi (`np.gradient`).

## Diagnostika a interpretace výstupu

CLI vypisuje v sekci `DRT Analysis` obě fáze hledání i výsledné λ:

```
Lambda: Hybrid GCV + L-curve
  lambda = 5.80e-10  (GCV 6.95e-10 -> L-curve corner 5.80e-10, ratio 0.83)
```

Když `ratio` vypadne z pásma konsenzu, přibude varování — buď o korekci
nahoru (`lcurve_correction`, typické pro NNLS), nebo o použití geometrického
průměru (`geometric_mean`). Samotný průběh hledání (Fáze 1 / Fáze 2) je na
úrovni DEBUG, tedy viditelný jen s `-v`.

`Lambda: GCV (L-curve correction failed)` znamená fallback po selhání
hybridního hledání, ne shodu obou kritérií.

**Klíčové indikátory:**
- `lambda_at_edge` / „λ na mezi rozsahu" — optimum je mimo `[1e-10, 1e-2]`;
  zkontroluj kvalitu dat (KK validace) nebo zadej λ ručně.
- `corner_at_edge` — roh L-křivky na okraji okna; skutečný roh může ležet dál.
- Velmi šumová data → GCV/L-curve mohou preferovat vyšší λ (over-smoothing).

## Kdy použít ruční λ

Automatika je dobrá baseline (objektivní, reprodukovatelná, vhodná pro batch a
publikace). Ruční `--lambda` má smysl pro fine-tuning u známého systému, při
specifických požadavcích na zvýraznění rysů, nebo když je auto λ na okraji
rozsahu.

## Validace

Relevantní testy:

```bash
python3 -m pytest tests/test_hybrid_lambda.py tests/test_lcurve_corner.py \
                  tests/test_drt_spikiness.py -q
```

- `tests/test_hybrid_lambda.py` — GCV score (konečný/kladný; GCV-vybrané λ není
  horší než okraje rozsahu), hybridní výběr, integrace s `_build_drt_matrices`.
- `tests/test_lcurve_corner.py` — znaménková konvence rohu (max kladné křivosti).
- `tests/test_drt_spikiness.py` — detekce degenerované / zašuměné γ.

## Reference

1. **Golub, Heath, Wahba (1979)** — *Generalized Cross-Validation as a Method
   for Choosing a Good Ridge Parameter*, Technometrics 21, 215-223.
2. **Wahba (1985)** — *A comparison of GCV and GML…*, Annals of Statistics 13,
   1378-1402.
3. **Hansen, P.C. (1998)** — *Rank-Deficient and Discrete Ill-Posed Problems*,
   SIAM. (L-curve a regularizační metody.)
4. **Saccoccio, Wan, Chen, Ciucci (2014)** — *Optimal regularization in DRT…*,
   Electrochimica Acta 147, 470-482.
5. **Maradesa, Py, et al. (2023)** — *Selecting the regularization parameter in
   the distribution of relaxation times*, J. Electrochem. Soc. 170, 030502.

---

*Dokument pro EIS Analysis Toolkit*
*Poslední aktualizace: 2026-06-30 (GCV + L-curve hybrid)*
