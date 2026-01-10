# GMM Peak Detection v DRT spektrech

**Verze:** 0.9.4
**Datum:** 2026-01-05
**Implementace:** `eis_analysis/drt/peaks.py`

## Přehled

Gaussian Mixture Models (GMM) poskytují robustnější a objektivnější metodu pro detekci píků v DRT (Distribution of Relaxation Times) spektrech než tradiční `scipy.signal.find_peaks`.

### Klíčové výhody

1. **Explicitní hranice píků**: τ ∈ [μ-2σ, μ+2σ] poskytuje 95% confidence interval
2. **Automatická dekonvoluce**: Separuje překrývající se píky
3. **Objektivní výběr**: BIC (Bayesian Information Criterion) automaticky určí počet píků
4. **Kvantifikace nejistoty**: Šířka píků v log-prostoru (σ)
5. **Přesnější odhady R_i**: Integrace gaussovských komponent místo ad-hoc hranice

## Použití

### Základní použití

```bash
# GMM detekce píků (místo výchozí scipy metody)
eis data.DTA --peak-method gmm
```

### S automatickou analýzou

```bash
# Kompletní automatizace: GMM + auto-lambda + auto-circuit
eis data.DTA --peak-method gmm --auto-lambda --auto-circuit

# S verbose výstupem pro diagnostiku
eis data.DTA --peak-method gmm --auto-lambda -vv
```

### Uložení výsledků

```bash
# Uložit grafy a skrýt display
eis data.DTA --peak-method gmm --save results --no-show
```

## Výstup

### Konzole

GMM vypíše detailní informace o každém detekovaném píku:

```
==================================================
GMM detekce píků
==================================================
Hledám optimální počet komponent v rozsahu (1, 6)
✓ Optimální počet píků: 2 (BIC=56.59)

Detekované píky (seřazeno podle τ):
  Pík 1:
    τ = 9.11e-04 s (f = 1.75e+02 Hz)
    Hranice τ: [3.25e-04, 2.55e-03] s
    Šířka (σ): 0.224 dekád
    Váha: 0.314
    R ~ 1876.67 Ω
  Pík 2:
    τ = 4.92e-02 s (f = 3.24e+00 Hz)
    Hranice τ: [1.72e-02, 1.41e-01] s
    Šířka (σ): 0.228 dekád
    Váha: 0.686
    R ~ 4178.99 Ω
```

### Grafy (layout 2×2)

Při použití GMM metody se vytvoří 4 subploty:

**Top-left: DRT spektrum**
- Standardní γ(τ) graf
- Základní vizualizace distribuce

**Top-right: Rekonstrukce Nyquist**
- Porovnání dat vs DRT rekonstrukce
- Kontrola kvality fitu

**Bottom-left: GMM dekonvoluce**
- Překryv jednotlivých gaussovských komponent (barevně odlišené)
- Vertikální čáry označují:
  - Tečkované: hranice píků (μ±2σ)
  - Přerušované: střed píku (μ)
- Legenda s τ hodnotami pro každý pík

**Bottom-right: BIC křivka**
- Model selection diagnostika
- X-osa: počet komponent
- Y-osa: BIC hodnota (nižší = lepší)
- Červená hvězda: optimální volba

## Jak GMM hledá komponenty

GMM (Gaussian Mixture Model) používá iterativní EM (Expectation-Maximization) algoritmus pro nalezení optimálních parametrů směsi Gaussiánů.

### Inicializace (k-means)

GMM nepoužívá náhodnou inicializaci, ale **k-means clustering**:

```
Data (log10(τ) vážené podle gamma):
••••    •••••••••    ••••••

K-means rozdělí do n shluků:
[----shluk 1----]  [--shluk 2--]

Počáteční odhad:
  μ₁ = centroid shluku 1
  μ₂ = centroid shluku 2
  σ₁, σ₂ = variance v každém shluku
  w₁, w₂ = podíl bodů v každém shluku
```

Parametr `random_state=42` zajišťuje reprodukovatelnost výsledků.

### EM iterace

**E-step (Expectation):** Pro každý bod vypočti pravděpodobnost příslušnosti ke každé komponentě:

```
P(bod → komponenta A) = w_A × N(bod|μ_A,σ_A) / Σ[w_i × N(bod|μ_i,σ_i)]
```

**M-step (Maximization):** Aktualizuj parametry podle vážených příslušností:

```
μ_nový = vážený průměr bodů podle pravděpodobností příslušnosti
σ_nový = vážená variance
w_nový = průměrná pravděpodobnost příslušnosti
```

### Průběh pro 2 komponenty

```
Iterace 0 (k-means init):
  Komponenta A: μ=-3.0, σ=0.5
  Komponenta B: μ=-1.0, σ=0.4

Iterace 1:
  E: Přepočítej příslušnosti bodů
  M: A: μ=-3.02, σ=0.48  B: μ=-1.05, σ=0.42

Iterace 2:
  E: Přepočítej příslušnosti
  M: A: μ=-3.01, σ=0.47  B: μ=-1.04, σ=0.41

... (typicky 5-20 iterací)

Konvergence:
  A: μ=-3.01, σ=0.47, w=0.35
  B: μ=-1.04, σ=0.41, w=0.65
```

### Klíčové vlastnosti

- **Měkké přiřazení:** Každý bod má pravděpodobnost příslušnosti ke všem komponentám (ne jen k jedné)
- **Dekonvoluce překryvu:** I překrývající se peaky se správně separují
- **Automatické vyvážení:** Váhy w se naučí z dat

## Technické detaily

### Algoritmus

1. **Transformace**: DRT se fituje v log₁₀(τ) prostoru pro lepší separaci píků
2. **Váhování**: Body jsou replikovány podle γ(τ) amplitudy (nové v sklearn 1.8+)
3. **BIC optimalizace**: Testuje n=1 až n=6 komponent
4. **Konzervativní výběr**: Preferuje jednodušší model pokud BIC improvement < 10
5. **Parametry píků**:
   - μ (mean): pozice píku v log₁₀(τ)
   - σ (std): šířka v log₁₀(τ) prostoru
   - w (weight): relativní váha komponenty
   - R_estimate: numerický integrál γ(τ) přes hranice píku (viz níže)

### BIC formula

```
BIC(n) = -2·log(L) + k·log(N)
```

kde:
- L = likelihood modelu s n komponentami
- k = počet parametrů (3n pro GMM: μ, σ, w pro každou komponentu)
- N = počet datových bodů

Nižší BIC = lepší kompromis mezi fitem a složitostí.

### Výpočet R_estimate a gamma_max

GMM přímo nepočítá odpor peaku ani jeho výšku. Tyto hodnoty se stanovují dodatečně:

**1. R_estimate - numerická integrace z DRT spektra**

Po detekci peaku (μ, σ) se odpor R_estimate počítá **numericky z vypočteného DRT spektra**:

```python
# Hranice peaku = 95% confidence interval
tau_lower = 10**(mu - 2*sigma)
tau_upper = 10**(mu + 2*sigma)

# Numerická integrace skutečných hodnot gamma
R_estimate = trapz(gamma[idx_lower:idx_upper+1], ln_tau[idx_lower:idx_upper+1])
```

Proč numericky ze skutečného DRT:
- GMM předpokládá symetrické Gaussiány, ale reálné peaky mohou být asymetrické
- Integrál ∫γ(τ)d(ln τ) přímo odpovídá odporu daného elektrochemického procesu
- Zachovává se skutečná plocha pod peakem, ne idealizovaná

**2. gamma_max - zpětný výpočet pro vizualizaci**

Výška peaku se dopočítává až při vizualizaci tak, aby rekonstruovaný Gaussián měl **stejnou plochu** jako integrovaný DRT peak:

```python
integral_log10 = sigma * sqrt(2 * pi)    # integrál Gaussiánu v log10 prostoru
integral_ln = integral_log10 * log(10)    # převod na ln prostor
gamma_max = R_estimate / integral_ln      # výška peaku
```

**Vzorec:**

```
gamma_max = R_estimate / (σ × √(2π) × ln(10))
```

Fyzikální význam: Gaussián s touto výškou má stejnou plochu (odpor R) jako skutečný DRT peak.

### Srovnání s scipy.find_peaks

| Vlastnost | `scipy.find_peaks` | GMM |
|-----------|-------------------|-----|
| Překrývající se píky | ❌ Jen jeden max | ✅ Dekonvoluce |
| Hranice píků | ❌ Nedefinované | ✅ μ±2σ (95% CI) |
| Šířka píků | ❌ | ✅ σ v log-prostoru |
| Odhad R_i | Manuální integrace | ✅ Z plochy gaussiánu |
| Objektivita | ⚠️ Hard thresholds | ✅ BIC |
| Rychlost | ✅✅ Velmi rychlé | ✅ Rychlé (~1-2s) |
| Závislosti | scipy | sklearn |

## Integrace s auto-circuit

Při použití `--auto-circuit`, GMM píky se automaticky použijí pro:

1. **Počet Voigtových článků**: Jeden článek R||C pro každý GMM pík
2. **Initial guess pro R_i**: Použije R_estimate z GMM
3. **Initial guess pro τ_i**: Použije τ_center z GMM
4. **Initial guess pro C_i**: Vypočte C_i = τ_i / R_i

Výsledek: Výrazně lepší konvergence circuit fitu.

## Požadavky

### Python balíčky

```bash
# Nutné pro GMM
pip install scikit-learn --break-system-packages

# Nebo systemově (Debian/Ubuntu)
sudo apt install python3-sklearn
```

### Fallback

Pokud sklearn není dostupný:
- GMM metoda se automaticky přepne na `scipy.find_peaks`
- Uživatel dostane warning zprávu
- Analýza pokračuje s fallback metodou

## Kdy použít GMM vs scipy

### ✅ Použij GMM pokud:
- Máš překrývající se píky
- Potřebuješ přesné hranice pro fyzikální interpretaci
- Chceš objektivní výběr počtu píků
- Potřebuješ publikovatelné výsledky s kvantifikací nejistoty

### 🔧 Použij scipy pokud:
- Píky jsou dobře separované
- Potřebuješ rychlou exploraci
- sklearn není dostupný
- Stačí ti kvalitativní analýza

## Příklady výstupů

### Příklad 1: Dva dobře separované píky

```
Input: R1=1000Ω, τ1=1ms, R2=5000Ω, τ2=50ms
GMM output:
  Pík 1: τ=0.91ms [0.33-2.6ms]  R~1877Ω
  Pík 2: τ=49.2ms [17-141ms]    R~4179Ω
BIC: n=2 optimal (BIC=56.59)
```

### Příklad 2: Překrývající se píky

```
Input: R1=2000Ω, τ1=5ms, R2=3000Ω, τ2=8ms (částečný překryv)
scipy: Detekuje 1 široký pík při τ~6ms
GMM: Separuje na 2 píky při τ1=4.8ms a τ2=8.3ms
```

## Limitace a budoucí vylepšení

### Současné limitace

1. **Gaussovský předpoklad**: Reálné píky mohou být asymetrické
2. **Pouze kladné γ**: Zatím nepodporuje záporné píky (induktivní smyčky)
3. **Fixní rozsah**: n=1 až n=6 komponent (může být nedostatečné pro složité systémy)

### Sloučení malých peaků do jednoho Gaussiánu

GMM může aproximovat několik blízkých malých peaků **jedním širokým Gaussiánem**:

```
Skutečné DRT spektrum:           GMM aproximace (n=1):

γ(τ)                             γ(τ)
  │    ▲                           │
  │   ███  ▲   ▲                   │      ████████
  │  █████ ██ ███                  │    ████████████
  │ ███████████████                │  ████████████████
  └─────────────────               └─────────────────
      3 malé peaky                 1 široký Gaussián

   R₁ + R₂ + R₃  ≈  R_estimate (zachováno)
```

**Kdy k tomu dochází:**
- Peaky jsou blízko u sebe v log(τ) prostoru
- Přidání další komponenty nezlepší BIC o více než threshold
- GMM vidí skupinu peaků jako jeden široký shluk

**Důsledky:**

| Parametr | Výsledek |
|----------|----------|
| R_estimate | ≈ součet R jednotlivých peaků (fyzikálně správné) |
| τ_center | vážený průměr pozic peaků |
| σ (šířka) | větší než u jednotlivých peaků |
| Počet procesů | podhodnocen |

**Kdy je to problém:**
- Potřebujete identifikovat jednotlivé elektrochemické procesy
- Různé procesy mají různou teplotní závislost nebo jinou charakteristiku

**Kdy je to žádoucí:**
- Šum vytváří falešné "picky"
- Blízké časové konstanty reprezentují jeden distribuovaný proces (CPE)

**Jak předejít sloučení:**

```bash
# Snížit BIC threshold = citlivější detekce
eis data.dta --gmm-bic-threshold 5
```

### Co GMM neověřuje

GMM pracuje **pouze s DRT spektrem** a nemá přístup k původním impedančním datům:

```
Tok dat v analýze:

Z(ω) měřená
    ↓
[Tikhonov regularizace] ← validace: ||A·γ - Z||²
    ↓
γ(τ) DRT spektrum
    ↓
[GMM peak detection] ← BIC vybírá "statisticky nejlepší" model
    ↓
Detekované peaky (μ, σ, R_estimate)
    ↓
[suggest_circuit] → návrh obvodu
    ↓
[CNLS fit] ← ZDE probíhá validace proti Z(ω)
    ↓
Fitované parametry + χ²
```

**Co GMM dělá:**
- Fituje směs Gaussiánů na distribuci pozic v log(τ) vážených podle γ
- Likelihood = pravděpodobnost, že pozorovaná distribuce pochází z modelu
- BIC penalizuje složitost modelu

**Co GMM nedělá:**
- Nerekonstruuje impedanci Z(ω) z detekovaných peaků
- Neporovnává rekonstrukci s měřenými daty
- Neověřuje fyzikální smysluplnost výsledků

**Důsledek:**
GMM může vrátit statisticky optimální model, který neodpovídá fyzikální realitě. Validace probíhá až při CNLS fitu navrženého obvodu, kde se porovnává model s měřenou impedancí.

### Plánovaná vylepšení

- [ ] **Log-normální směs**: Lepší pro asymetrické píky
- [ ] **Signed GMM**: Podpora záporných γ pro induktivní procesy
- [ ] **Bayesian GMM**: Dirichlet Process pro automatický počet komponent bez horní meze
- [ ] **Uncertainty propagation**: Chyby R_i a C_i z GMM kovariancí
- [ ] **L-curve visualization**: Alternativní diagnostika k BIC

## Reference

### Teorie GMM

- Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer. Chapter 9.
- Murphy, K. P. (2012). *Machine Learning: A Probabilistic Perspective*. MIT Press. Chapter 11.

### BIC v kontextu regularizace

- Schwarz, G. (1978). "Estimating the dimension of a model". *Annals of Statistics* 6(2), 461-464.
- Burnham, K. P., & Anderson, D. R. (2004). *Multimodel Inference*. Springer.

### DRT aplikace

- Saccoccio, M. et al. (2014). "Optimal Regularization in Distribution of Relaxation Times applied to Electrochemical Impedance Spectroscopy". *Electrochimica Acta* 147, 470-482.

## Troubleshooting

### "ModuleNotFoundError: No module named 'sklearn'"

```bash
pip install scikit-learn --break-system-packages
```

### "GMM fit selhal pro n=X"

Obvykle neškodné - některé hodnoty n mohou selhat kvůli numerické nestabilitě. GMM automaticky přeskočí tyto hodnoty.

### "Všechny GMM fity selhaly"

Vzácné, ale možné při:
- Velmi malém počtu bodů
- Extrémně šumných datech
- Téměř nulové γ(τ)

Řešení: Automatický fallback na scipy.find_peaks.

### GMM najde příliš mnoho/málo píků

Zkus:
1. Změnit regularizaci λ (méně píků → zvětši λ)
2. Použít `--auto-lambda` pro optimální λ
3. Zvýšit `--n-tau` pro lepší rozlišení

## Příspěvky

Máte-li nápady na vylepšení GMM implementace:
- Otevřete issue na GitHubu
- Navrhněte pull request s vylepšením
- Sdílejte příklady dat kde GMM selhává

---

**Autor:** EIS Analysis Toolkit
**Licence:** MIT
