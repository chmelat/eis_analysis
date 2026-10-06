# Lin-KK implementace v eis_analysis

## Prehled

Modul `eis_analysis.validation.kramers_kronig` poskytuje nativni implementaci Lin-KK testu (Linear Kramers-Kronig) pro validaci kvality EIS dat. Intuitivni vysvetleni principu je v [KK_INTUITION.md](KK_INTUITION.md). Implementace je zalozena na metode Schonleber et al. (2014) a krome numpy nevyzaduje externi zavislosti. Pro `fit_type='real'` i `'complex'` dava stejne M, mu i rezidua jako referencni `impedance.py` (`impedance.validation.linKK`), az na odchylky popsane v sekci 3 (start hledani M) a 5 (rozsireni tau).

**Reference:**
- Schonleber, M. et al. "A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
- Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
- Yrjana, V. and Bobacka, J. "Implementing Kramers-Kronig validity testing using pyimpspec." Electrochim. Acta 504, 144951 (2024)

---

## 1. Teoreticky zaklad

Lin-KK test fituje impedancni data pomoci rady Voigtovych elementu s pevnymi casovymi konstantami:

```
Z_fit(w) = R_s + sum_{k=1}^{M} R_k / (1 + jw*tau_k) + jw*L [+ 1/(jw*C)]
```

Kde:
- `R_s` = seriovy odpor [Ohm]
- `R_k` = odpor k-teho Voigtova elementu [Ohm], muze vyjit zaporny
- `tau_k` = casova konstanta (fixni, logaritmicky rozlozena) [s]
- `L` = seriova induktance [H]
- `C` = seriova kapacita [F], jen s `include_C=True` (CLI `--kk-series-c`)
- `M` = pocet elementu (volen pomoci mu metriky, sekce 3)

**Klicove:** Casove konstanty tau_k jsou fixni. Fitovany jsou pouze odpory R_k a seriove cleny (R_s, L, pripadne C), takze fit je linearni.

Kazdy clen modelu splnuje KK relace. Kdyz model data nedokaze popsat, data KK relace porusuji (drift, nestacionarita, nelinearita, artefakty pristroje).

Seriova kapacita je KK-kompatibilni, ale ma nulovou realnou cast, takze ji Voigtuv retezec neumi reprezentovat. U blokujicich systemu (dvouelektrodove cely, pasivni vrstvy) bez ni rostou imaginarni rezidua k nizkym frekvencim, i kdyz data KK splnuji.

---

## 2. Distribuce casovych konstant

### Zakladni rozsah

```python
tau_min = 1 / (2*pi*f_max)
tau_max = 1 / (2*pi*f_min)
```

Casove konstanty jsou logaritmicky rozlozeny v intervalu [tau_min, tau_max]:

```
tau_k = 10^[log10(tau_min) + (k-1)/(M-1) * log10(tau_max/tau_min)]
```

### Rozsireni rozsahu (extend_decades)

Parametr `extend_decades` rozsiruje rozsah tau **jen smerem k nizkym frekvencim** (k delsim tau):

```python
tau_max_extended = tau_max * 10^(extend_decades)
tau_min          = beze zmeny
```

- `extend_decades = 0.0`: zadne rozsireni (standardni Lin-KK)
- `extend_decades = 0.4`: tau_max posunuto o 0.4 dekady vys
- zaporne hodnoty se orezou na 0

Pocet elementu M se rozsirenim nemeni, jen se zvetsi rozestup mezi tau_k. Vysokofrekvencni induktivni chvost pokryva clen L, ne rozsireni.

---

## 3. Mu metrika a volba poctu elementu

### Algoritmus

1. Najdi `M_lower`: prvni M, od ktereho u 3 po sobe jdoucich M (`CHI2_PLATEAU_RUN`) lezi log10(pseudo chi^2) do 0.3 dekady od minima pres nasledujicich 8 M (`CHI2_PLATEAU_DECADES`, `CHI2_PLATEAU_WINDOW`). Jedine M by mohlo padnout do nahodneho poklesu chi^2 na hrubem gridu.
2. Zacni s M = M_lower (puvodni Lin-KK zacina s M=3)
3. Fituj pomoci pseudoinverze (povol zaporne R_k)
4. Vypocitej mu metriku
5. Pokud mu > threshold a M < max_M, zvys M o 1 a opakuj
6. Kdyz mu <= threshold, zastav: toto M je vysledek

Krok 1 je odchylka od puvodniho Lin-KK. Mu neni v M monotonni: na hrube mrizce tau se relaxace mezi dvema body mrizky fituje stridavymi znamenky R_k, mu na chvili klesne pod prah a puvodni algoritmus zastavi s nedostatecne rozlisenym modelem (3.8 % rezidua na presnem 2-RC spektru). Pseudo chi^2 ukazuje, od ktereho M je mrizka dost husta. Okno je lokalni (8 M), protoze na driftujicich datech chi^2 dal pomalu klesa, jak zaporne R_k absorbuji drift.

`max_M` (vychozi 50) se omezuje poctem bodu N: fit si musi nechat aspon jeden stupen volnosti (M <= N - 2). Pod 5 body `lin_kk_native` vyhodi `ValueError`.

### Mu metrika (Schonleber 2014)

```
mu = 1 - (sum|R_k| pro R_k < 0) / (sum|R_k| pro R_k >= 0)
```

**Interpretace:**
- mu = 1: zadne zaporne R_k. Model muze byt i nedostatecne rozliseny, mu samo o kvalite fitu nic nerika.
- mu klesa pod prah (vychozi 0.85): zaporne R_k zacinaji mit vyznamnou vahu. Pridavani dalsich elementu by vedlo k preuceni.
- Vysledne `mu` je hodnota, pri ktere hledani zastavilo, takze je obvykle pod prahem. Nad prahem je jen pri dosazeni max_M (s varovanim).

**Proc zaporne R_k znamenaji preuceni?**
- Fyzikalne by mely byt vsechny R_k >= 0 (odpor nemuze byt zaporny)
- Zaporne R_k vznikaji, kdyz model ma prilis mnoho volnosti
- Model zacina fitovat sum misto skutecneho signalu

**Mu neni kriterium kvality dat.** Kvalitu dat posuzuji rezidua (sekce 4).

---

## 4. Metriky kvality

### Pseudo chi-squared (Boukamp 1995)

```
chi2_ps = sum_i w_i * [(Z'_exp - Z'_fit)^2 + (Z''_exp - Z''_fit)^2]
```

kde `w_i = 1/|Z_i|^2` je Boukampova vaha.

### Odhad sumu (Yrjana & Bobacka 2024)

```
noise_estimate = sqrt(chi2_ps * 5000 / N) [%]
```

kde N je pocet bodu. Je to horni odhad: do chi^2 prispiva i nedokonalost modelu a KK poruseni. Na datech, ktera KK porusuji, proto vyjde vysoky (4.8 % u `real_gamry_example.DTA`, sekce 9) a neni to sum mereni.

### Rezidua

```
res_real = (Z'_exp - Z'_fit) / |Z_exp|
res_imag = (Z''_exp - Z''_fit) / |Z_exp|
```

Rezidua dobreho fitu jsou nahodne rozptylena kolem nuly, radove na urovni sumu mereni (u bezneho potenciostatu obvykle pod 1 %). Systematicky prubeh (hrb, trend k jednomu konci spektra) je znakem KK poruseni i pri malych hodnotach.

### Kriterium platnosti (`is_valid`)

Pro kazdy bod se bere vetsi slozka `max(|res_real|, |res_imag|)`. Data jsou platna, kdyz nanejvys 5 % bodu (`KK_MAX_FRACTION_ABOVE`) ma tuto hodnotu nad 5 % (`KK_RESIDUAL_THRESHOLD`). Mez 5 % odpovida teckovanym caram v grafu reziduii. Povoleny podil toleruje nekolik krajnich bodu, kde Lin-KK rezidua prirozene rostou (pri N = 72 jsou to 3 body).

Obe hodnoty jsou empiricke, hrube "zjevne rozbite", ne test na urovni sumu. Prumer pres vsechny body (kriterium do v0.46.0) schoval lokalni poruseni: 10 z 70 bodu s rezidui 20 % da prumer 2.9 %.

### Popisek kvality v CLI

CLI k verdiktu vypisuje popisek podle `max(prumer |res_real|, prumer |res_imag|)`:
- < 0.5 %: excellent
- < 1 %: good
- < 2.5 %: acceptable
- < 5 %: marginal (check for drift/nonlinearity)
- neplatna data (`is_valid` False): poor

Rozhoduje verdikt: neplatna data jsou vzdy "poor", platna nejhure "marginal".

### Vazeni pri fittingu

Pouzivame `modulus` vazeni (w = 1/|Z|) podle Schonleber (2014), nikoli `proportional` (w = 1/|Z|^2) podle Boukamp (1995):

- **1/|Z|** - vyrovnane relativni vazeni pres cele spektrum
- **1/|Z|^2** - silny duraz na body s nizkym |Z| (obvykle vysoke frekvence), nizkofrekvencni oblast ma maly vliv

Pro typicka EIS data je 1/|Z| vhodnejsi, protoze nizkofrekvencni oblast casto obsahuje klicove elektrochemicke informace (prenos naboje, difuze) a bezne artefakty (drift, nestacionarita).

**Pozn.:** Pseudo chi^2 se pocita s vahou 1/|Z|^2 podle Boukampa - rozdil je pouze ve fittingu, ne ve vysledne metrice.

---

## 5. Automaticka optimalizace extend_decades

### Motivace

Prilis uzky rozsah tau vede k vysokym reziduim na nizkofrekvencnim konci, kdyz relaxace pokracuje za nejnizsi merenou frekvenci (kapacitni chvost). Prilis siroky rozsah pri pevnem M zhrubne mrizku a fit muze kompenzovat stridavymi zapornymi R_k, tedy preucenim, ktere ma mu metrika zachytit.

### Implementace

M se urci nejdriv bez rozsireni (sekce 3). Pak `find_optimal_extend_decades()` pri tomto M projde mrizku 11 hodnot v `search_range` (vychozi (0.0, 1.0)) a vybere tu s nejnizsim pseudo chi^2:

- Kandidat, jehoz vlastni mu klesne pod mu, pri kterem hledani M zastavilo, se zahodi (`min_mu`). Rozsireni tedy nesmi obejit pojistku proti preuceni.
- Mezi kandidaty do 0.1 % od minima chi^2 vyhrava nejmensi rozsireni.
- Kdyz zadny kandidat neprojde, funkce vrati `None`, zustane nerozsirena mrizka a do `warnings` se prida varovani.

### CLI pouziti

V CLI je optimalizace **zapnuta ve vychozim stavu**:

```bash
eis data.DTA                           # s optimalizaci, rozsah 0 az 1 dekada
eis data.DTA --extend-decades-max 2.0  # rozsah 0 az 2 dekady
eis data.DTA --no-auto-extend          # bez rozsireni
```

Zvolene rozsireni je v souhrnnem radku, napr. `extend_decades=0.60` pro `example/EISPOT-test1.DTA` (sekce 9).

---

## 6. API reference

Uplne signatury a popis parametru jsou v docstringech. Zde je prehled.

### KKResult dataclass

Vysledek `kramers_kronig_validation()`:

```python
@dataclass
class KKResult:
    M: int                           # Pocet Voigtovych elementu (0 pri chybe)
    M_lower: int                     # Kde hledani M zacalo (plato chi^2)
    mu: float                        # Mu, pri kterem hledani zastavilo
    Z_fit: Optional[NDArray]         # Fitovana impedance
    residuals_real: Optional[NDArray]  # Rezidua realne casti (podil |Z|)
    residuals_imag: Optional[NDArray]  # Rezidua imaginarni casti (podil |Z|)
    pseudo_chisqr: float             # Pseudo chi^2 (Boukamp 1995)
    noise_estimate: float            # Horni odhad sumu [%]
    extend_decades: float            # Pouzite rozsireni tau
    inductance: Optional[float]      # Seriova induktance [H]
    capacitance: Optional[float]     # Seriova kapacita [F] (jen s include_C)
    elements: Optional[NDArray]      # Fitovane prvky [R_s, R_1, ..., R_M, L]
    tau: Optional[NDArray]           # Casove konstanty [s]
    weighting: str                   # Vazeni fitu ('modulus')
    warnings: List[str]              # Varovani (max_M, odmitnute rozsireni)
    error: Optional[str]             # Chybova zprava, kdyz validace selhala

    # Odvozene vlastnosti
    success: bool                    # Validace probehla (Z_fit existuje)
    mean_residual_real: float        # Prumer |res_real| [%]
    mean_residual_imag: float        # Prumer |res_imag| [%]
    n_above_threshold: int           # Body s max(|res_real|, |res_imag|) > 5 %
    is_valid: bool                   # Nanejvys 5 % bodu nad 5 % (sekce 4)
```

### kramers_kronig_validation()

```python
def kramers_kronig_validation(
    frequencies, Z,
    mu_threshold=0.85,
    max_M=50,
    auto_extend_decades=True,
    extend_decades_range=(0.0, 1.0),
    include_C=False
) -> KKResult
```

Vysokourovnova funkce: vzdy `fit_type='real'`, `weighting='modulus'` a L v modelu. Pri chybe nevyhazuje vyjimku, ale vrati `KKResult` s `error` a `success == False`.

### lin_kk_native()

```python
def lin_kk_native(
    frequencies, Z,
    mu_threshold=0.85,
    max_M=50,
    include_L=True,
    include_C=False,
    fit_type='real',            # 'real', 'imag' nebo 'complex'
    weighting='modulus',        # 'uniform', 'sqrt', 'modulus', 'proportional'
    auto_extend_decades=False,
    extend_decades_range=(0.0, 1.0)
) -> KKResult
```

Samotny fit s plnou kontrolou parametru. Vraci stejny `KKResult` jako `kramers_kronig_validation()`, jen chybu nevraci v `error`, ale vyhodi ji (pod 5 body `ValueError`). `LinKKResult` je alias `KKResult`, zachovany kvuli kompatibilite.

### Pomocne funkce

```python
compute_pseudo_chisqr(Z_exp, Z_fit) -> float
estimate_noise_percent(chi2_ps, n_points) -> float
reconstruct_impedance(frequencies, elements, tau, L_value, include_L=True, C_value=None) -> NDArray
find_optimal_extend_decades(frequencies, Z, M, ...) -> Optional[Tuple[...]]
```

Konstanty kriteria platnosti: `KK_RESIDUAL_THRESHOLD`, `KK_MAX_FRACTION_ABOVE` (v `eis_analysis.validation.kramers_kronig`).

### Graf

Validace sama zadny graf nevytvari. Graf fitu a reziduii vykresli `eis_analysis.visualization.plot_kk_validation(frequencies, Z, result, flagged_frequencies=())`; `flagged_frequencies` oznaci body cervenymi pasy (CLI tam predava body z kontroly jednotlivych bodu).

---

## 7. Priklady pouziti

### Python API

```python
from eis_analysis.validation import kramers_kronig_validation

result = kramers_kronig_validation(frequencies, Z)
if not result.success:
    print(f"KK validace selhala: {result.error}")
else:
    print(f"M = {result.M}, mu = {result.mu:.3f}, extend = {result.extend_decades:.2f}")
    print(f"Estimated noise (upper bound): {result.noise_estimate:.2f}%")
    print(f"Points above threshold: {result.n_above_threshold}/{len(frequencies)}")
    print(f"Valid: {result.is_valid}")
    for w in result.warnings:
        print(f"Warning: {w}")

    from eis_analysis.visualization import plot_kk_validation
    fig = plot_kk_validation(frequencies, Z, result)

# Dvouelektrodova cela nebo jiny blokujici system
result = kramers_kronig_validation(frequencies, Z, include_C=True)
print(f"Series C: {result.capacitance:.2e} F")
```

### CLI

```bash
# Pouze KK validace (bez DRT, fitu obvodu a Z-HIT)
eis data.DTA --no-drt --no-fit --no-zhit

# Se seriovou kapacitou (blokujici chovani)
eis data.DTA --no-drt --no-fit --no-zhit --kk-series-c

# Bez zobrazeni grafu, graf ulozen do souboru
eis data.DTA --no-drt --no-fit --no-zhit --no-show --save vysledek
```

---

## 8. Porovnani metod fitu

`kramers_kronig_validation()` vzdy pouziva `fit_type='real'`. Ostatni typy jsou dostupne pres `lin_kk_native()`.

### fit_type = 'real' (vychozi)

1. Fituje pouze realnou cast (vazenou 1/|Z|) -> ziska [R_s, R_1, ..., R_M]
2. Z techto parametru predpovi imaginarni cast
3. Ze zbytku Im(Z_exp - Z_fit) fituje L (a C, je-li zapnuto)

**Vyhody pro validaci:**
- Imaginarni cast je dopocitana z realne pres KK relace, takze imaginarni rezidua primo ukazuji KK nekonzistenci
- Pokud data splnuji KK relace, imaginarni cast by mela automaticky sedet

### fit_type = 'imag'

Zrcadlove: fituje imaginarni cast, R_s se dopocita z realne. Realna rezidua slouzi jako diagnostika.

### fit_type = 'complex'

Fituje Re(Z) a Im(Z) soucasne. Poruseni se rozlozi do obou slozek, takze maximum reziduii je nizsi. Pro validaci proto neni doporuceno.

Na `example/real_gamry_example.DTA`, ktery KK porusuje (sekce 9), poruseni zachyti vsechny tri typy:

| fit_type | M | prumer \|res_real\| / \|res_imag\| | max bod | bodu > 5 % |
|---|---|---|---|---|
| real | 22 | 0.37 / 3.82 % | 20.4 % | 19 z 72 |
| imag | 21 | 5.12 / 0.16 % | 22.4 % | 14 z 72 |
| complex | 21 | 1.89 / 2.51 % | 11.2 % | 20 z 72 |

---

## 9. Interpretace vysledku

Priklady jsou skutecny vystup CLI (`--no-drt --no-fit --no-zhit`) na spektrech z adresare `example/`.

### Priklad: Dobra data (`EISPOT-test1.DTA`)

```
KK: M=19 (from M=4, chi^2 plateau), mu=0.8456 (Lin-KK stop, threshold 0.85), extend_decades=0.60
  Mean |res_real|: 0.08%
  Mean |res_imag|: 0.41%
  Pseudo chi^2: 1.65e-03
  Estimated noise (upper bound): 0.34%
Data quality: excellent (0/72 points above 5.0%, allowed 5%; max mean |res|=0.41%)
```

- Zadny bod nad 5 %, nejvetsi reziduum 1.1 %
- Odhad sumu 0.34 % odpovida beznemu potenciostatu

### Priklad: Blokujici system (`EISPOT-M136113-4.DTA`, ZrO2 na Zr, dvouelektrodove)

Bez `--kk-series-c` data projdou (mean |res_imag| 0.59 %), ale imaginarni rezidua systematicky rostou k nejnizsim frekvencim (4.3 % v poslednim bodu 1.6 mHz). Se seriovou kapacitou:

```
KK: M=41 (from M=16, chi^2 plateau), mu=0.8084 (Lin-KK stop, threshold 0.85), extend_decades=0.00
  Mean |res_real|: 0.05%
  Mean |res_imag|: 0.15%
  Pseudo chi^2: 3.22e-04
  Estimated noise (upper bound): 0.13%
  Series C: 3.53e-05 F
Data quality: excellent (0/89 points above 5.0%, allowed 5%; max mean |res|=0.15%)
```

- Nizkofrekvencni trend zmizel, maximum klesne na 0.37 %: byl to chybejici seriovy clen, ne KK poruseni

### Priklad: Problematicka data (`real_gamry_example.DTA`)

```
KK: M=22 (from M=8, chi^2 plateau), mu=0.8477 (Lin-KK stop, threshold 0.85), extend_decades=0.00
  Mean |res_real|: 0.37%
  Mean |res_imag|: 3.82%
  Pseudo chi^2: 3.37e-01
  Estimated noise (upper bound): 4.84%
! Data quality: poor (19/72 points above 5.0%, allowed 5%; max mean |res|=3.82%)
```

- Imaginarni rezidua tvori hladky hrb jednoho znamenka mezi 0.03 a 4 Hz s maximem 20 % u 0.25 Hz. Realna cast sedi do 3 %.
- Prumer (3.82 %) je pod 5 %. Poruseni ukazuje az pocet bodu nad mezi.
- Odhad sumu 4.84 % neni sum mereni, ale dusledek poruseni (sekce 4)
- Mozne priciny: nestacionarita, artefakty, nelinearita

---

## 10. Literatura

1. Schonleber, M., Klotz, D., and Ivers-Tiffee, E. (2014). A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests. Electrochim. Acta, 131, 20-27.

2. Boukamp, B.A. (1995). A Linear Kronig-Kramers Transform Test for Immittance Data Validation. J. Electrochem. Soc., 142, 1885-1894.

3. Yrjana, V. and Bobacka, J. (2024). Implementing Kramers-Kronig validity testing using pyimpspec. Electrochim. Acta, 504, 144951.
