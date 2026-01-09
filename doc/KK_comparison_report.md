# Technicka zprava: Srovnani implementaci Kramers-Kronig validace

**Projekt**: eis_analysis vs pyimpspec
**Datum**: 2026-01-09 (aktualizace)
**Autor**: Analyza provedena pomoci Claude Code

---

## 1. Uvod

Tato zprava poskytuje podrobne srovnani implementaci Kramers-Kronig (KK) validace ve dvou Python knihovnach pro analyzu elektochemicke impedancni spektroskopie (EIS):

- **eis_analysis** - lokalni projekt v aktivnim vyvoji
- **pyimpspec** - etablovana open-source knihovna

Kramers-Kronig validace je zakladni nastroj pro overeni kvality EIS dat. Vyuziva kauzality a linearity systemu k detekci experimentalnich artefaktu, sumu a nestacionarnich jevu.

---

## 2. Teoreticky zaklad

Obe implementace vychazi z metody **Lin-KK** (Linear Kramers-Kronig), ktera aproximuje impedancni spektrum radou Voigtovych elementu s pevnymi casovymi konstantami:

```
Z(w) = R_s + sum_{k=1}^{M} R_k / (1 + jw*tau_k) + jw*L
```

Kde:
- R_s: seriovy odpor
- R_k: odpor k-teho Voigtova elementu
- tau_k: casova konstanta (fixni, logaritmicky rozlozena)
- L: seriova induktance
- M: pocet elementu (optimalizovany parametr)

Hlavni reference:
- Boukamp, B.A. (1995): J. Electrochem. Soc., 142, 1885-1894
- Schonleber, M. et al. (2014): Electrochim. Acta, 131, 20-27

---

## 3. Srovnani implementaci

### 3.1 Prehled

| Charakteristika | eis_analysis | pyimpspec |
|-----------------|--------------|-----------|
| Zakladni metoda | Lin-KK | Lin-KK |
| Varianty testu | 1 | 7 |
| Metody volby M | 1 | 6 |
| Reprezentace | Impedance | Impedance + Admitance |
| Pseudo chi^2 | Ano | Ano |
| Odhad sumu | Ano | Ano |
| Optimalizace tau rozsahu | Automaticka (--auto-extend) | Automaticka |

### 3.2 Varianty testu

**eis_analysis** implementuje jednu variantu:
- Real-part fitting s extrakci induktance z imaginarnich residuii

**pyimpspec** nabizi 7 variant:

| Varianta | Popis | Solver |
|----------|-------|--------|
| complex | Simultanni fit realne a imaginarni casti | lstsq |
| real | Fit realne casti | lstsq |
| imaginary | Fit imaginarni casti | lstsq |
| complex-inv | Komplexni fit | Matricova inverze |
| real-inv | Fit realne casti | Matricova inverze |
| imaginary-inv | Fit imaginarni casti | Matricova inverze |
| cnls | Komplexni nelinearni fit | lmfit |

**Poznamka k volbe varianty:** Pro ucely **validace** je real-part fitting metodologicky spravny pristup (Schonleber 2014). KK relace svazuji realnou a imaginarni cast - pokud fitujeme pouze Z', imaginarni cast Z'' by mela automaticky sedet pro KK-kompatibilni data. Rezidua v Z'' pak slouzi jako diagnostika KK compliance. Complex fitting muze maskovat KK nekonzistence tim, ze optimalizuje obe casti soucasne.

### 3.3 Vahove funkce

Vahova funkce ovlivnuje relativni dulezitost bodu v ruznych castech spektra:

| Knihovna | Vahova funkce | Nazev |
|----------|---------------|-------|
| eis_analysis | w = 1/\|Z\| | Modulus (Lin-KK standard) |
| pyimpspec | w = 1/\|Z\|^2 | Proportional (Boukamp 1995) |

**Poznamka k terminologii (v0.11.0):** eis_analysis pouziva nazvy dle bezne EIS konvence:
- `modulus` = w = 1/\|Z\| (vyvazene vazeni pres cele spektrum)
- `proportional` = w = 1/\|Z\|^2 (Boukamp 1995)

Default pro KK validaci je `modulus` (1/\|Z\|), coz:
- Poskytuje vyvazenou relativni vahu pres cele spektrum
- Nizkofrekencni oblast (vysoka \|Z\|) obsahuje klicove elektrochemicke informace (prenos naboje, difuze)
- Bezne artefakty (drift, nestacionarita) se objevuji na nizkych frekvencich a nebyly by s 1/\|Z\|^2 dostatecne zvyrazneny

### 3.4 Automaticka volba poctu RC elementu

**eis_analysis** pouziva jednu metodu:

**Mu criterion** (Schonleber 2014):
```
mu = 1 - (sum|R_k < 0|) / (sum|R_k >= 0|)
```
- Iterativni zvysovani M od 3
- Zastaveni pri mu <= 0.85
- Negativni R_k indikuji overfit

**pyimpspec** kombinuje 6 metod:

| # | Metoda | Princip |
|---|--------|---------|
| 1 | Mu criterion | Pomer negativnich/pozitivnich R_k |
| 2 | Norma parametru | Minimalizace \|\|theta\|\| / M |
| 3 | Norma krivosti | Minimalizace krivosti Z_fit |
| 4 | Zmeny znamenka krivosti | Pocet oscilaci ve fitu |
| 5 | Vzdalenost mezi zmenami | Prumerna vzdalenost oscilaci |
| 6 | Variance log tau | Kubicka aproximace trendu |

Kombinovany pristup: kazda metoda vraci score 0-1, vyber nejlepsiho kandidata.

### 3.5 Distribuce casovych konstant

**eis_analysis**:
```python
tau_min = 1 / (2*pi*f_max)
tau_max = 1 / (2*pi*f_min) * 10^extend_decades

# extend_decades pracuje v CASOVE domene - rozsiruje tau_max
# smerem k vyssim hodnotam (tj. nizsim frekvencim)
# Default: 0.0, s --auto-extend: automaticka optimalizace v rozsahu (0, max)
```

**pyimpspec**:
```python
F_ext = 10^log_F_ext  # faktor rozsireni
tau_min = 1 / (w_max * F_ext)
tau_max = F_ext / w_min
# Automaticka optimalizace log_F_ext v rozsahu [-1, 1]
```

Obe knihovny podporuji automatickou optimalizaci rozsireni tau rozsahu:
- **eis_analysis**: `--auto-extend` flag s `--extend-decades-max` parametrem
- **pyimpspec**: automaticky pri `log_F_ext=0` a `num_F_ext_evaluations>1`

### 3.6 Vypocet residuii

Obe knihovny pouzivaji stejnou formuli dle Schonleber (2014):

```
res_real = (Z_real_exp - Z_real_fit) / |Z_exp|
res_imag = (Z_imag_exp - Z_imag_fit) / |Z_exp|
```

### 3.7 Metriky kvality

**eis_analysis**:
- Mu metrika s prahem 0.85
- Vizualni reference +-5% pro rezidua
- Pseudo chi^2 (Boukamp 1995): `sum(w_i * [dRe^2 + dIm^2])` kde `w_i = 1/|Z|^2`
- Odhad sumu: `noise% = sqrt(chi^2 * 5000 / N)` (Yrjana & Bobacka 2024)
- KKResult dataclass s property `is_valid` (True pokud rezidua < 5%)

**pyimpspec**:
- Pseudo chi^2 (Boukamp 1995): `sum(w_i * [dRe^2 + dIm^2])`
- Odhad sumu: `SD = sqrt(chi^2 * 5000 / N)`
- Statisticke testy normality residuii:
  - Lilliefors test
  - Shapiro-Wilk test
  - Kolmogorov-Smirnov test

**Poznamka:** eis_analysis neimplementuje statisticke testy normality residuii (Lilliefors, Shapiro-Wilk), coz je oblast pro potencialni rozsireni.

### 3.8 Admitancni reprezentace

**eis_analysis**: Nepodporovano

**pyimpspec**: Plne podporovano
- Uzitecne pro systemy s negativni diferencialni rezistenci
- Prepina mezi seriovym a paralelnim zapojenim elementu
- Automaticka detekce vhodne reprezentace

---

## 4. Implementacni detaily

### 4.1 Struktura kodu

**eis_analysis**:
```
eis_analysis/
  validation/
    kramers_kronig.py      # Hlavni KK funkce, KKResult dataclass
                           # compute_pseudo_chisqr(), estimate_noise_percent()
                           # find_optimal_extend_decades(), reconstruct_impedance()
  fitting/
    voigt_chain/
      __init__.py          # Public API export
      fitting.py           # estimate_R_linear(), fit_voigt_chain_linear()
      mu_optimization.py   # calc_mu(), find_optimal_M_mu()
      tau_grid.py          # generate_tau_grid(), generate_tau_grid_fixed_M()
      solvers.py           # compute_voigt_matrix()
      validation.py        # Vstupni data validace
```

**pyimpspec**:
```
pyimpspec/
  analysis/
    kramers_kronig/
      single.py            # Hlavni API
      result.py            # Vysledkova trida
      exploratory.py       # F_ext optimalizace
      least_squares.py     # LS solver
      matrix_inversion.py  # Legacy solver
      cnls.py              # CNLS solver
      algorithms/
        method_1.py        # Mu criterion
        method_2.py        # Norma parametru
        method_3.py        # Norma krivosti
        method_4.py        # Zmeny znamenka
        method_5.py        # Vzdalenosti
        method_6.py        # Log tau variance
```

### 4.2 API design

**eis_analysis** - strukturovane API s dataclass vysledkem:
```python
def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0)
) -> Optional[KKResult]

@dataclass
class KKResult:
    M: int                              # Pocet Voigt elementu
    mu: float                           # Mu metrika
    Z_fit: NDArray[np.complex128]       # Fitovana impedance
    residuals_real: NDArray[np.float64] # Rezidua realne casti
    residuals_imag: NDArray[np.float64] # Rezidua imaginarni casti
    pseudo_chisqr: float                # Pseudo chi^2 (Boukamp 1995)
    noise_estimate: float               # Odhad sumu v %
    extend_decades: float               # Pouzite rozsireni tau
    inductance: Optional[float]         # Fitovana induktance [H]
    figure: Optional[plt.Figure]        # Vizualizace
    # Properties: mean_residual_real, mean_residual_imag, is_valid
```

**pyimpspec** - rozsahle API:
```python
def perform_kramers_kronig_test(
    data: DataSet,
    test: str = "real",
    num_RC: int = 0,
    add_capacitance: bool = True,
    add_inductance: bool = True,
    admittance: Optional[bool] = None,
    min_log_F_ext: float = -1.0,
    max_log_F_ext: float = 1.0,
    log_F_ext: float = 0.0,
    num_F_ext_evaluations: int = 20,
    rapid_F_ext_evaluations: bool = True,
    cnls_method: str = "leastsq",
    max_nfev: int = 0,
    timeout: int = 60,
    num_procs: int = -1,
) -> KramersKronigResult
```

### 4.3 Zavislosti

**eis_analysis**:
- numpy
- scipy (optimize.nnls)
- matplotlib

**pyimpspec**:
- numpy
- scipy
- lmfit
- pandas
- matplotlib

---

## 5. Hodnoceni

### 5.1 Silne stranky eis_analysis

1. **Jednoduchost**: Srozumitelna implementace, snadne pochopeni algoritmu
2. **Rychlost**: Jedna metoda = rychlejsi vypocet
3. **Minimalismus**: Mene zavislosti, snadnejsi integrace
4. **Dokumentace**: Dobra matematicka dokumentace v kodu
5. **Flexibilita**: Snadna modifikace pro specificke potreby
6. **Strukturovany vystup**: KKResult dataclass s convenience properties
7. **Kvantitativni metriky**: Pseudo chi^2 a odhad sumu dle literatury

### 5.2 Silne stranky pyimpspec

1. **Komplexnost**: Vice variant pro ruzne typy dat
2. **Robustnost**: 6 metod pro volbu M zvysuje spolehlivost
3. **Statistika**: Formalni testy normality residuii
4. **Admitance**: Podpora pro problematicka data
5. **Automatizace**: Optimalizace F_ext bez uzivatelske intervence
6. **Ekosystem**: Soucasti vetsiho analytickeho frameworku

### 5.3 Slabe stranky eis_analysis

1. Chybi admitancni reprezentace (pro systemy s negativni diferencialni rezistenci)
2. Chybi statisticke testy normality residuii (Lilliefors, Shapiro-Wilk)
3. Jedna metoda volby M (pouze mu criterion, pyimpspec ma 6 metod)

Poznamka: Real-part fitting (jedina varianta) je pro validaci metodologicky spravny - neni to slabina, ale vlastnost.

### 5.4 Slabe stranky pyimpspec

1. Vyssi komplexita API
2. Vice zavislosti
3. Obtiznejsi customizace
4. Pomalejsi (vice vypoctu)

---

## 6. Doporuceni pro eis_analysis

### 6.1 Implementovano (od v0.10.0)

Nasledujici funkce byly implementovany od puvodni verze tohoto reportu:

1. **Pseudo chi^2** - Implementovano v `compute_pseudo_chisqr()` (Boukamp 1995)
2. **Odhad sumu** - Implementovano v `estimate_noise_percent()` (Yrjana & Bobacka 2024)
3. **Automaticka optimalizace extend_decades** - CLI flag `--auto-extend`
4. **Strukturovany vystup** - KKResult dataclass s convenience properties

### 6.2 Stredni priorita (potencialni rozsireni)

1. **Statisticke testy residuii**
   - Shapiro-Wilk test pro normalitu
   - Objektivni kriterium kvality dat (nad ramec vizualni kontroly)

2. **Dalsi metody volby M**
   - Method 3 (norma krivosti) jako doplnek k mu criterion
   - Zvyseni robustnosti automaticke volby pro specificke typy dat

### 6.3 Nizka priorita

3. **Admitancni reprezentace**
   - Pouze pro specialni pripady (negativni dif. rezistence)
   - Vyzaduje vyraznejsi zmeny v architekture
   - Typicka EIS data tuto funkcionalitu nevyzaduji

### 6.4 Neni potreba

- **Complex fitting** - Real-part fitting je pro validaci metodologicky spravny. Complex fitting by maskoval KK nekonzistence a snizil diagnostickou hodnotu testu.

---

## 7. Zaver

Obe implementace jsou zalozeny na solidnim teoretickem zakladu (Lin-KK) a produkuji srovnatelne vysledky pro standardni EIS data. Od verze 0.10.0 eis_analysis implementuje vsechny klicove metriky (pseudo chi^2, odhad sumu) a automatickou optimalizaci.

Hlavni rozdily:

1. **Rozsah funkcionalit**: pyimpspec ma vice variant testu (7 vs 1) a metod volby M (6 vs 1)
2. **Filozofie designu**: eis_analysis preferuje jednoduchost a metodologickou spravnost
3. **Statisticke testy**: pyimpspec ma testy normality residuii (Lilliefors, Shapiro-Wilk)
4. **Admitancni reprezentace**: pouze pyimpspec (pro specialni pripady)

Pro beznou KK validaci jsou obe knihovny funkcne ekvivalentni. Implementace v eis_analysis pouziva metodologicky spravny pristup (real-part fitting), ktery umoznuje diagnostiku KK compliance pres rezidua v imaginarni casti.

Pro pokrocile analyzy (statisticke testy normality, data s negativni diferencialni rezistenci) nabizi pyimpspec vice nastroju. Pro typicke EIS analyzy je eis_analysis plne dostatecna.

---

## Reference

1. Boukamp, B.A. (1995). A Linear Kronig-Kramers Transform Test for Immittance Data Validation. J. Electrochem. Soc., 142, 1885-1894.

2. Schonleber, M., Klotz, D., and Ivers-Tiffee, E. (2014). A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests. Electrochim. Acta, 131, 20-27.

3. Plank, C., Ruther, T., and Danzer, M.A. (2022). A novel approach for model selection in linear Kramers-Kronig analysis. IWIS Workshop, 1-6.

4. Yrjana, V. and Bobacka, J. (2024). Implementing Kramers-Kronig validity testing using pyimpspec. Electrochim. Acta, 504, 144951.
