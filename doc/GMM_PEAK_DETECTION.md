# GMM Peak Detection v DRT spektrech

**Verze:** 0.16.6
**Datum:** 2026-06-30
**Implementace:** `eis_analysis/drt/peaks.py`

## Přehled

Vážená Gaussovská směs (GMM) poskytuje robustnější a objektivnější metodu pro
detekci píků v DRT (Distribution of Relaxation Times) spektrech než tradiční
`scipy.signal.find_peaks`.

Implementace je **vlastní vážený EM** (`_weighted_gaussian_mixture_1d`) postavený
jen na `numpy` a `scipy.special.logsumexp` — **nevyžaduje scikit-learn ani žádnou
další závislost** nad rámec core balíčků projektu.

### Klíčové výhody

1. **Explicitní hranice píků**: τ ∈ [10^(μ−2σ), 10^(μ+2σ)] (95% interval)
2. **Automatická dekonvoluce**: separuje překrývající se píky
3. **Objektivní výběr**: BIC (Bayesian Information Criterion) určí počet píků
4. **Kvantifikace nejistoty**: šířka píku v log-prostoru (σ)
5. **Konzistentní odhad R_i**: rozklad R_pol podle vah komponent (Σ R_i = R_pol)

## Použití

### CLI

Vstupní bod je `eis.py` (spouštěj `python3`).

```bash
# GMM detekce píků (místo výchozí scipy metody)
python3 eis.py data.DTA --peak-method gmm

# Konzervativnější / citlivější výběr počtu píků (BIC práh, default 10)
python3 eis.py data.DTA --peak-method gmm --gmm-bic-threshold 5

# Vyšší rozlišení tau gridu
python3 eis.py data.DTA --peak-method gmm --n-tau 150

# Uložit grafy a skrýt display
python3 eis.py data.DTA --peak-method gmm --save results --no-show
```

> Výběr regularizace λ je **automatický by default** (GCV + L-curve, viz
> [GCV_IMPLEMENTATION.md](GCV_IMPLEMENTATION.md)); samostatný „auto" přepínač
> neexistuje, ruční λ se zadá přes `--lambda`.

### Relevantní přepínače

| Přepínač | Popis | Default |
|----------|-------|---------|
| `--peak-method` | `scipy` (rychlé) nebo `gmm` (robustní) | `scipy` |
| `--gmm-bic-threshold` | Min. BIC zlepšení pro přidání komponenty (vyšší = méně píků) | 10.0 |
| `--n-tau` / `-n` | Počet tau bodů (rozlišení DRT) | 100 |

### Python API

```python
from eis_analysis.drt import gmm_peak_detection, calculate_drt

# Přímo na DRT spektru
peaks, gmm_model, bic_scores = gmm_peak_detection(tau, gamma, bic_threshold=10.0)

# Nebo přes plnou DRT analýzu
result = calculate_drt(freq, Z, peak_method='gmm')
```

## Výstup

### Konzole

```
============================================================
GMM detekce píků
============================================================
Hledám optimální počet komponent v rozsahu (1, 6)
n=1 komponenty: BIC=...
n=2 komponenty: BIC=56.59
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

Při GMM metodě se vykreslí 4 subploty (plot fce v `drt/core.py`):

- **Top-left: DRT spektrum** — γ(τ).
- **Top-right: Nyquist rekonstrukce** — data vs. DRT rekonstrukce.
- **Bottom-left: GMM dekonvoluce** — jednotlivé Gaussovské komponenty (barevně);
  tečkované čáry = hranice μ±2σ, přerušované = střed μ. Výška každého Gaussiánu
  je škálovaná tak, aby měl plochu rovnou R_estimate píku (viz `gamma_max` níže).
- **Bottom-right: BIC křivka** — BIC vs. počet komponent; svislá čára značí
  zvolené n (nižší BIC = lepší).

## Jak GMM hledá komponenty

GMM fituje směs Gaussiánů na **vážená** data v log₁₀(τ) prostoru pomocí
iterativního EM (Expectation-Maximization). DRT γ(τ) je vážená hustota, ne
náhodný vzorek — proto váhy vstupují přímo do EM (žádná replikace bodů).

### Inicializace (kvantilová, deterministická)

Není použito k-means ani náhodný start. Středy komponent se inicializují na
rovnoměrně rozmístěných **vážených kvantilech** dat:

```
x = log10(τ) pro biny s γ > 0,  váhy w = γ
cdf = kumulativní Σw / Σw
μ_k = interpolace x na kvantilech (k+0.5)/K,  k = 0..K-1
σ²_k = (vážený celkový rozptyl) / K
w_k  = 1/K
```

Pro 1D je kvantilový init stabilní a reprodukovatelný bez `random_state`.

### Vážený EM

**E-krok:** spočti log-responsibility (přes `logsumexp`, numericky stabilní):

```
P(bod → komponenta k) ∝ w_k · N(x | μ_k, σ²_k)
```

**M-krok (vážený):** aktualizuj parametry s váhami `w·resp`:

```
N_k = Σ w·resp_k
π_k = N_k / Σw
μ_k = Σ(w·resp_k·x) / N_k
σ²_k = Σ(w·resp_k·(x−μ_k)²) / N_k  + reg_covar
```

Iteruje do konvergence log-likelihood (typicky jednotky až desítky iterací).

### Klíčové vlastnosti

- **Měkké přiřazení:** každý bin přispívá všem komponentám podle responsibility.
- **Dekonvoluce překryvu:** překrývající se píky se separují.
- **Vážení γ:** váhy = γ normalizované na Σw = n_data (počet měření), takže
  log-likelihood má škálu měření a `bic_threshold` je srovnatelný napříč
  datasety (náprava auditu F1).
- **Podlaha rozptylu (`var_floor`):** (½ mediánu rozestupu mřížky)² brání tomu,
  aby se komponenta smrskla na jediný bin (singularita / nereálně úzký pík).

## Technické detaily

### Algoritmus

1. **Transformace**: směs se fituje v log₁₀(τ) prostoru.
2. **Vážení**: body váženy přímo γ (vážený EM, **bez replikace bodů** — audit F1);
   váhy normalizovány na Σw = n_data.
3. **BIC optimalizace**: testuje n = 1 až 6 komponent.
4. **Konzervativní výběr (Occam)**: přidává komponentu jen pokud BIC zlepšení >
   `bic_threshold`, jinak zastaví.
5. **Parametry píku**: μ (log₁₀τ), σ (šířka v dekádách), w (váha), R_estimate.

### BIC

```
BIC(n) = −2·log(L) + k·log(N)
```

kde:
- `L` = (vážená) likelihood modelu s n komponentami,
- `k = 3n − 1` (parametry 1D směsi: μ(n) + σ²(n) + π(n−1)),
- `N` = počet skutečných **měření (frekvencí)**.

Nižší BIC = lepší kompromis fit/složitost. Výběr používá early-stopping
(přidávej komponenty, dokud zlepšení > `bic_threshold`); na okrajích rozsahu
(n=1 nebo n=6) se vypíše varování.

### Výpočet R_estimate

GMM nepočítá odpor píku přímo. R_i se odvodí z DRT spektra rozkladem celkového
polarizačního odporu podle vah komponent (náprava auditu F10):

```python
# Celkový R_pol obdélníkovým pravidlem (konzistentní s jádrem DRT)
d_ln_tau = mean(diff(ln(tau)))
R_pol = sum(gamma) * d_ln_tau

# Dílčí odpor = podíl podle váhy komponenty (Σ weight = 1)
R_estimate = weight * R_pol      # => Σ R_i = R_pol přesně
```

Proč rozklad podle vah a ne integrace γ přes ±2σ okno každého píku: okna
překrývajících se píků by integrovala **totéž** γ vícekrát (double-counting),
takže Σ R_i by neodpovídalo R_pol. Rozklad podle vah to řeší přesně.

### gamma_max (jen pro vizualizaci)

Výška píku se dopočítává **pouze pro graf** dekonvoluce tak, aby rekonstruovaný
Gaussián měl stejnou plochu jako R_estimate:

```
gamma_max = R_estimate / (σ · √(2π) · ln(10))
```

(plocha Gaussiánu v log₁₀ prostoru `σ·√(2π)` převedená na ln prostor faktorem
`ln 10`). Není to vstup do žádného výpočtu, jen měřítko křivky v grafu.

### Srovnání s scipy.find_peaks

| Vlastnost | `scipy.find_peaks` | GMM (vážený EM) |
|-----------|-------------------|-----------------|
| Překrývající se píky | ❌ jen jeden max | ✅ dekonvoluce |
| Hranice píků | ❌ nedefinované | ✅ μ±2σ (95% interval) |
| Šířka píků | ❌ | ✅ σ v log-prostoru |
| Odhad R_i | per-pík integrace | ✅ rozklad R_pol podle vah |
| Objektivita | ⚠️ pevné prahy | ✅ BIC |
| Rychlost | ✅✅ velmi rychlé | ✅ rychlé (~1-2 s) |
| Závislosti | numpy/scipy | numpy/scipy |

## Integrace s návrhem obvodu

GMM píky slouží jako podklad pro návrh ekvivalentního obvodu. V DRT výstupu se
automaticky vypíše **Voigt report** (`analyze_voigt_elements` /
`format_voigt_report` z `fitting/auto_suggest.py`):

1. **Počet R‖C článků**: jeden na každý GMM pík.
2. **Initial guess R_i**: R_estimate z GMM.
3. **Initial guess τ_i**: τ_center z GMM.
4. **Initial guess C_i**: C_i = τ_i / R_i.

Samotný fit se spustí zadáním obvodu přes `--circuit` (např.
`--circuit 'R()-(R()|Q())-(R()|Q())'`). Lepší initial guess z GMM zlepšuje
konvergenci fitu.

## Požadavky

GMM vyžaduje pouze core závislosti projektu (`numpy`, `scipy`) — **žádná extra
instalace** (scikit-learn není potřeba).

### Fallback

Pokud vážený EM selže numericky (např. téměř nulová γ, extrémně málo binů),
`gmm_peak_detection` vrátí prázdný seznam a DRT analýza se automaticky vrátí ke
`scipy.signal.find_peaks` (ty se počítají vždy a slouží jako fallback).

## Kdy použít GMM vs scipy

### ✅ Použij GMM pokud:
- Máš překrývající se píky
- Potřebuješ přesné hranice pro fyzikální interpretaci
- Chceš objektivní výběr počtu píků
- Potřebuješ publikovatelné výsledky s kvantifikací šířky

### 🔧 Použij scipy pokud:
- Píky jsou dobře separované
- Potřebuješ nejrychlejší exploraci
- Stačí kvalitativní analýza

## Příklady výstupů

### Příklad 1: Dva dobře separované píky

```
Vstup: R1=1000Ω, τ1=1ms, R2=5000Ω, τ2=50ms
GMM:   Pík 1: τ=0.91ms [0.33-2.6ms]  R~1877Ω
       Pík 2: τ=49.2ms [17-141ms]    R~4179Ω
BIC: n=2 optimální (BIC=56.59)
```

### Příklad 2: Překrývající se píky

```
Vstup: R1=2000Ω, τ1=5ms, R2=3000Ω, τ2=8ms (částečný překryv)
scipy: detekuje 1 široký pík při τ~6ms
GMM:   separuje na 2 píky při τ1≈4.8ms a τ2≈8.3ms
```

## Limitace

1. **Gaussovský předpoklad**: reálné píky mohou být asymetrické.
2. **Pouze kladné γ**: záporné píky (induktivní smyčky) zatím nepodporovány.
3. **Fixní rozsah**: n = 1 až 6 komponent.

### Sloučení malých peaků do jednoho Gaussiánu

GMM může aproximovat několik blízkých malých peaků **jedním širokým Gaussiánem**,
pokud přidání další komponenty nezlepší BIC nad práh. Důsledky:

| Parametr | Výsledek |
|----------|----------|
| R_estimate (součet) | ≈ součet R sloučených peaků (R_pol zachováno) |
| τ_center | vážený průměr pozic |
| σ (šířka) | větší než u jednotlivých peaků |
| Počet procesů | podhodnocen |

**Žádoucí**, když jde o šum nebo jeden distribuovaný proces (CPE); **problém**,
když potřebuješ rozlišit jednotlivé procesy. Citlivější detekci vynutíš nižším
prahem: `--gmm-bic-threshold 5`.

### Co GMM neověřuje

GMM pracuje **pouze s DRT spektrem** a nemá přístup k původní impedanci:

```
Z(ω) → [Tikhonov] → γ(τ) → [GMM] → píky (μ, σ, R_estimate)
                                      → [Voigt report / --circuit] → [CNLS fit] ← validace proti Z(ω)
```

GMM tedy vrací statisticky nejlepší rozklad γ, ne fyzikálně ověřený model —
validace proti měřené impedanci probíhá až při CNLS fitu navrženého obvodu.

## Reference

- Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer,
  kap. 9.
- Murphy, K. P. (2012). *Machine Learning: A Probabilistic Perspective*. MIT
  Press, kap. 11.
- Schwarz, G. (1978). *Estimating the dimension of a model*. Annals of
  Statistics 6(2), 461-464.
- Saccoccio, M. et al. (2014). *Optimal Regularization in DRT applied to EIS*.
  Electrochimica Acta 147, 470-482.

## Troubleshooting

**„GMM fit selhal pro n=X"** — neškodné; některé n mohou selhat kvůli numerické
nestabilitě, GMM je přeskočí (BIC = inf).

**„Všechny GMM fity selhaly"** — vzácné (velmi málo bodů, extrémní šum, téměř
nulová γ). Automatický fallback na `scipy.find_peaks`.

**GMM najde příliš mnoho/málo píků** — uprav `--gmm-bic-threshold` (vyšší = méně
píků), případně regularizaci λ (ruční `--lambda`) nebo rozlišení `--n-tau`.

---

**Autor:** EIS Analysis Toolkit
**Licence:** MIT
