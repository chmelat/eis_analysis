# GMM Peak Detection v DRT spektrech

**Verze:** 1.6.0
**Datum:** 2025-12-12
**Implementace:** `eis_analysis_1_6.py`

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
python3 eis_analysis_1_6.py data.DTA --peak-method gmm
```

### S automatickou analýzou

```bash
# Kompletní automatizace: GMM + auto-lambda + auto-circuit
python3 eis_analysis_1_6.py data.DTA --peak-method gmm --auto-lambda --auto-circuit

# S verbose výstupem pro diagnostiku
python3 eis_analysis_1_6.py data.DTA --peak-method gmm --auto-lambda -vv
```

### Uložení výsledků

```bash
# Uložit grafy a skrýt display
python3 eis_analysis_1_6.py data.DTA --peak-method gmm --save results --no-show
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
   - R_estimate: plocha gaussiánu = w × γ_max × σ × √(2π) × ln(10)

### BIC formula

```
BIC(n) = -2·log(L) + k·log(N)
```

kde:
- L = likelihood modelu s n komponentami
- k = počet parametrů (3n pro GMM: μ, σ, w pro každou komponentu)
- N = počet datových bodů

Nižší BIC = lepší kompromis mezi fitem a složitostí.

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

### Plánovaná vylepšení (v1.7+)

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

**Autor:** EIS analýza toolkit v1.6
**Licence:** MIT (nebo dle projektu)
**Kontakt:** [vaše email/GitHub]
