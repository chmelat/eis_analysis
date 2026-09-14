# Zahner Analysis - co si z toho vzít

Analýza manuálu **Zahner Analysis** (35 stran, 11/2023) v kontextu tohoto projektu:
které myšlenky stojí za převzetí, které už máme (a lépe) a které patří jinam.

Dokument je z větší části popis GUI, ale kapitoly 2.2 (*How Fitter Works*) a 2.3
(*Impedance electrical elements*) jsou svým způsobem malá učebnice.

## Kalibrace: v čem je projekt napřed

Automatická volba λ (GCV, L-curve), DRT, GMM detekce peaků, Lin-KK s μ-kritériem,
DE a multistart místo ručně zadávaných start values, analytický jakobián,
kovariance a konfidenční intervaly, AIC/BIC, reziduální diagnostika.
Zahner tohle buď nemá, nebo to nechává na uživateli.

Z-HIT máme také. Tafel, Butler-Volmer, CV a solární články jsou mimo záběr projektu.

## Přehled návrhů

| # | Návrh | Přínos | Práce | Dotčené moduly |
|---|-------|--------|-------|----------------|
| 1 | Significance - **HOTOVO v 0.33.0** (plot zatím ne) | vysoký | malá | `fitting/diagnostics.py`, `circuit.py`, `diffevo.py` |
| 2 | Young-Göhr element - **HOTOVO v 0.37.0** | vysoký | střední | `circuit_elements/composite.py`, `analysis/oxide.py` |
| 3 | Blokující difuze (coth) | střední | malá | `circuit_elements/distributed.py` |
| 4 | Log-polární residuum (clog) | střední | střední | `fitting/diagnostics.py`, `circuit.py` |
| 5 | Series Fit (série spekter) | vysoký | velká | CLI, handlery, vizualizace |
| 6 | Rovina komplexní kapacity, C(f) | střední | malá | `visualization/plots.py` |
| 7 | Fit na Z-HIT rekonstrukci - **HOTOVO v 0.34.0** | malý | triviální | `cli/parser.py`, `handlers/validation.py` |
| 8 | Element beroucí podobvod (PE) | koncepční | velká | `circuit_elements/base.py` |

---

## 1. Significance - nejostřejší myšlenka dokumentu

> **Stav: implementováno ve verzi 0.33.0.** `compute_significance()`
> v `eis_analysis/fitting/diagnostics.py`, pole `FitResult.params_significance`,
> sloupec `S=` u každého parametru ve výpisu fitu. Významy symbolů a odvození
> níže platí; podrobnosti pro uživatele jsou v `doc/WEIGHTING_AND_STATISTICS.md`
> §3.7. Significance plot S_i(f) implementován **není** - viz konec sekce.

Zahner k parametru hlásí **dvě různá čísla**: `error` (přesnost) a `significance`

```
S_i = max_n ( d|Z_n| * P_i / (dP_i * |Z_n|) )
```

Je to citlivost, ne nejistota, a odpovídá na jinou otázku:
*je ten prvek v tomhle frekvenčním okně vůbec k něčemu?*

### Význam symbolů

| Symbol | Význam |
|--------|--------|
| `Z` | **celková impedance modelu** (`Z_theo = f(ω, P_1, ..., P_k)`), ne impedance jednoho prvku |
| `\|Z_n\|` | její **modul** při n-té experimentální frekvenci ω_n; fáze do vzorce nevstupuje |
| `P_i` | i-tý **skalární fitovací parametr** přenosové funkce, ne prvek |
| `n` | index experimentálního vzorku (maximum jde přes měřené frekvence) |

Dokument je u `Z` explicitní: *"reflects its influence on the **network**
impedance"*. Proto má rezistor 338 nΩ v sérii s obloukem o kΩ nízkou
significance - jeho změna neposune výslednou impedanci sítě, i když sám o sobě
je dobře definovaný.

`P_i` je parametr, ne element: rezistor přispívá jedním, CPE dvěma (V a α),
Young-Göhr třemi. "Significance prvku" ve výsledkovém okně je ve skutečnosti
significance jednotlivých parametrů.

Vyhodnocuje se **v nafitovaném optimu** a **jen na měřených frekvencích** -
significance je tedy vlastnost modelu *v daném měřicím okně*, ne modelu obecně.

### Proč zrovna tenhle podíl

Přeuspořádáním je vidět logaritmická derivace:

```
d|Z|/|Z|  /  dP/P  =  d ln|Z| / d ln P
```

Bezrozměrné, takže R [Ω] a C [F] leží na téže škále - jednotky se vykrátí.
Obyčejné `d|Z|/dP` by porovnatelné nebylo.

Pro sériový rezistor je `dZ/dR = 1`, tedy `d|Z|/dR = Re(Z)/|Z| = cos φ` a

```
S = R * cos φ / |Z|   <= 1,   rovno 1 právě když Z = R
```

S je tedy zhruba **největší podíl |Z|, za který ten parametr odpovídá**.
Odtud plyne celá kalibrace dokumentu:

- rezistor dominující části spektra dá S ~ 1
- S << 0.01 znamená, že prvek lze z modelu vypustit
  (jejich příklad: R = 338 nΩ, S = 0.002)
- nelineárně vstupující parametry (exponent α u CPE) mohou 1 přesáhnout

Poslední bod je vidět přímo: `|Z_CPE| = (1/ω_0 V)*(ω/ω_0)^(-α)`, takže
`ln|Z| = -α*ln(ω/ω_0) + konst.` a

```
d ln|Z| / d ln α = -α * ln(ω/ω_0)
```

což roste se vzdáleností od normalizační frekvence bez omezení - tři dekády
od ω_0 dají při α ~ 0.9 hodnotu kolem 6. Není to chyba, jen u nelineárně
vstupujícího parametru mizí interpretace "podíl na |Z|".

**Proč to chceme:** `covariance.py` dá stderr a condition number. Velký stderr ale
míchá dohromady dvě různé situace - *parametr je nepodstatný* a *parametr je
zaměněný s jiným*. Significance ty dvě věci rozplete.

**Cena:** téměř nulová. `jacobian.py` už analytický J počítá, zbytek je

```
d|Z|/dP = (Re Z * dReZ/dP + Im Z * dImZ/dP) / |Z|
```

a jeden `max` přes vzorky. Zapadá to přímo do porovnávání `--circuit` variant
a do `auto_suggest` - automatické "tenhle prvek lze odebrat" je přesně
deklarovaný cíl *Automation*.

**Odchylka od zdroje:** Zahnerův vzorec bere maximum **znaménkové** veličiny,
bez vnějších svislic. Rozdíl nastane u parametru, jehož zvýšení |Z| *snižuje* -
takovému Zahner vyhodnotí nízkou significance, zatímco `max |...|` vysokou.
Implementace používá absolutní hodnotu

```
S_i = max_n | d ln|Z_n| / d ln P_i |
```

protože ptát se chceme na *velikost* vlivu, ne na jeho směr, a práh
"S << 0.01 -> lze vypustit" dává smysl jen pro nezáporné S. Je to vědomý
odklon od zdroje, ne jeho reprodukce, a je uvedený v docstringu
`compute_significance()`.

**Bonus:** significance plot S_i(f), tedy křivka na frekvenci pro každý prvek -
"kde který prvek řídí spektrum". Vizuálně komplementární k DRT peakům,
nic podobného v projektu není.

## 2. Young-Göhr element - trefa do oxidové domény

> **Stav: implementováno ve verzi 0.37.0.** Element `YG(C, p, tau)`
> v `eis_analysis/fitting/circuit_elements/composite.py`, analytický jakobián
> v `jacobian.py`, napojení na oxidovou analýzu v `analysis/oxide.py`.
> Uživatelský popis je v `doc/CIRCUIT_PARSER.md` (sekce YG) a
> v `doc/OXIDE_ANALYSIS_GUIDE.md`. Dvě věci nad rámec zdroje:
>
> - `e^(1/p)` přeteče float64 už pro `p < 1/709`, tedy **uvnitř** rozumného
>   rozsahu `p`. Implementace proto vzorec nikdy nepočítá jak je vytištěný -
>   `1/p` se přičítá až v log prostoru (`_yg_log_terms`).
> - `R_dc = p·τ·(e^(1/p) − 1)/C` je sice korektní limita ω→0, ale leží
>   `e^(1/p)` pod kapacitním rohem, tedy desítky dekád mimo jakékoli okno.
>   Do heuristiky „největší R = bariéra" **nevstupuje**; reportuje se zvlášť
>   spolu s frekvencí, pod kterou by ho bylo vidět.
>
> Navíc oproti CPE cestě: `δ = p·d`, hloubka průniku vodivosti v nm.

```
Z_Y = p/(jωC) * ln[ (1 + jωτ*e^(1/p)) / (1 + jωτ) ]
```

Model pasivní vrstvy s exponenciálně klesající vodivostí od povrchu.

| Parametr | Význam | Jednotka |
|----------|--------|----------|
| C | kapacita vrstvy (vysokofrekvenční limita) | F |
| p | relativní hloubka průniku vodivosti δ/d | - |
| τ | časová konstanta v místě nejvyšší vodivosti | s |

Dokument to říká přímo:

> If the origin of the CPE is known, it should be substituted by a more accurate
> model, like a Porous Electrode or the Young-Göhr element.

Uváděné aplikace: *oxide layers on metal electrodes like Fe, Al, Ti, and Ta* a
*organic coatings under soaking*.

**Proč to chceme:** `analysis/oxide.py` dnes fituje `Q`, pak přes Hsu-Mansfeld
**a** Brug dopočítá C_eff dvěma soupeřícími vzorci, hlídá jejich divergenci
a teprve z toho určí tloušťku. Young-Göhr fituje `C` přímo a `p` je fyzikálně
tloušťkový poměr - CPE→C konverzi nepotřebuje vůbec. Pro fitovaný oxid je to
koncepčně čistší cesta k témuž číslu a zároveň nezávislá kontrola stávající.

Souvisí s existujícím `CC` (Cole-Cole): obojí je "nahraď CPE něčím, co má fyziku".

Dokument uvádí i CPE aproximaci fáze Young-Göhra:

```
φ = -90° * (1 - q)     kde  q = 1 / (ln(ω*τ) + 1/p)
```

což je použitelné jako sanity check při implementaci.

## 3. Chybějící difuzní element: blokující terminace (coth)

Zahner rozlišuje čtyři difuzní případy, projekt má tři:

| Zahner | vzorec | v projektu |
|--------|--------|------------|
| Warburg | W/sqrt(jω) | `W` |
| Nernst (konečná délka, konstantní koncentrace) | W/sqrt(jω) * tanh sqrt(jω/k) | `Wo` (identické, τ = 1/k) |
| **Finite diffusion (blokující)** | W/sqrt(jω) * **coth** sqrt(jω/k) | **chybí** |
| Spherical diffusion | W/sqrt(jω + k) | `GE` |
| Homogeneous reaction (Gerischer) | W*/sqrt(k + jω) | `GE` |

Blokující varianta se pro ω→0 chová jako **kondenzátor** (ne rezistor jako Nernst) -
standardní model interkalační elektrody a blokujícího rozhraní obecně.
Implementačně ~30 řádků zrcadlících `Wo`, včetně analytické derivace.

**Zjištění zdarma:** spherical diffusion a Gerischer jsou **matematicky týž element**:

```
Z_R = W/sqrt(jω + k) = (W/sqrt(k)) / sqrt(1 + jω/k)      s τ = 1/k
```

`GE` tedy už umí i mikroelektrody, mikroelektrodová pole a bodovou korozi -
jen to nikde nestojí. Řádek do docstringu.

Ověření ekvivalence `Wo` = Zahnerův Nernst:

```
Wo  = R_W * tanh(sqrt(jωτ)) / sqrt(jωτ)
Z_N = W/sqrt(jω) * tanh(sqrt(jω/k_N))
```

s τ = 1/k_N a W = R_W * sqrt(k_N) jsou identické.

## 4. Log-polární residuum (clog) - jiná filozofie fitu

Zahner neminimalizuje kartézské residuum, ale komplexní logaritmus:

```
E = Σ_n | wgt * clog( Z_theo(ω_n) / Z_exp(ω_n) ) |² / N

wgt(real, imag) = real²/weight + imag² * weight
```

Reálná část clog je odchylka v ln|Z|, imaginární odchylka ve fázi. Parametr
`weight` (default **2.2222**) je ladicí knoflík na váhu fáze. Zdůvodnění v
dokumentu: logaritmické škálování dává stejnou váhu malým i velkým hodnotám
napříč dekádami - což u EIS platí.

Ukončovací kritéria: absolutní E < 0.001, nebo relativní zlepšení poslední
iterace < 1e-7.

**Vztah k projektu:** fitujeme `(ΔRe, ΔIm) * w` s `w = 1/|Z|`. To je aproximace
prvního řádu téhož (relativní chyba), ale ne totéž - chybová plocha má jiný tvar
a fázi nelze vážit nezávisle na modulu.

Návrh: `--weighting log-polar` + `--phase-weight`. Zapadá to k existujícím
`doc/WEIGHTING_AND_STATISTICS.md` a `doc/LEVM_WEIGHTING_REPORT.md` a dá se to
srovnat s LEVM váhováním na stejných datech.

Dokument nabízí i kalibrační škálu fitovací chyby:

- < 3 % - dobrý výsledek
- \> 10 % - přehodnotit model

Pozor: jejich procento je definované jejich metrikou, převzít lze myšlenku
kalibrované škály, ne konkrétní čísla.

Dále explicitní varování, které stojí za převzetí do dokumentace:

> The fitting error should not be used to estimate parameter uncertainties
> because the exact value differs strongly for different impedance element
> parameters and different models.

## 5. Series Fit - největší praktická mezera

Zahner fituje **sérii spekter jedním modelem najednou**:

1. fitne se první spektrum (ručně, pečlivě)
2. výsledek n-tého fitu jsou start values pro (n+1)-ní
3. výstupem je graf parametrů proti sériové proměnné

Sériová proměnná je typicky čas, ale uživatel si může přidat libovolný sloupec
(potenciál, teplota, cyklus). Výsledek se ukládá jako jeden objekt (`.zsfx`).

**Stav projektu:** CLI bere jeden soubor (`input, nargs='?'`).

Filozoficky je to nejsilnější bod dokumentu: **EIS se skoro nikdy neměří jednou.**
Měří se trend - nasákavost povlaku, růst oxidu, degradace. Reálný výstupní
artefakt tedy není "report jednoho fitu", ale trajektorie R_ct(t), C(t), α(t).

Kaskádování start values navíc drasticky zlevňuje globální optimalizaci:
plné DE / multistart je potřeba jen u prvního spektra.

Je to největší kus práce ze všech bodů (CLI, handlery, vizualizace,
`*Result` dataclassy), ale i největší změna v tom, co s nástrojem jde dělat.

## 6. Reprezentace: komplexní kapacita a C(f)

Dokument otevírá právě touhle myšlenkou: různé reprezentace zaostřují na různé
části spektra a dobrý fit musí sedět **ve všech**. Zahner jich nabízí 12,
projekt kreslí Nyquist + Bode.

Dvě nejužitečnější pro naši doménu:

**Komplexní kapacita** C* = 1/(jωZ), roviny C''(C') a C'(f), C''(f):

```
C' = Y''/(2πf)        C'' = -i * Y'/(2πf)
```

Docstring našeho `CC` elementu sám říká, že *the arc is depressed in the complex
capacitance plane, which is where the physics of an oxide or polymer film lives* -
a ten graf neexistuje. Element fitujeme v rovině, kterou neumíme zobrazit.

**Paralelní kapacita C_p(f)** jako model-free předanalýza (1/Z = 1/R + jωC):
plateau při nízkých frekvencích dá kapacitu vrstvy bez jakéhokoli fitu.
Nezávislá kontrola oxidové tloušťky proti fitované cestě, prakticky zadarmo.

Implementačně: rozšíření `visualize_data` o další panely, žádná nová matematika.

## 7. Fit na Z-HIT rekonstruovaných datech

> **Stav: implementováno ve verzi 0.34.0.** Přepínač `--fit-on
> {original,zhit,all}`; `apply_zhit_reconstruction()` v
> `eis_analysis/cli/handlers/validation.py`. Hodnota `all` jde nad rámec
> Zahnera - propíše rekonstrukci i do R_inf, DRT a oxidové analýzy, protože
> drift nízkofrekvenčního modulu kazí DRT stejně jako fit. Uživatelský popis
> je v README, sekce "Z-HIT as a correction".

Zahner v Toolboxu nabízí fitovat proti *Original / Smoothed / Z-HIT*.

Z-HIT tedy není jen validátor ("byl systém stacionární?"), ale i **oprava**:
u driftujícího nízkofrekvenčního konce se |Z| zrekonstruuje z fáze a fituje se
rekonstrukce. Typický případ z dokumentu: povlak nasákávající vodu, jehož
impedance během měření nízkých frekvencí klesá.

U nás byla celá instalatérská část hotová (`zhit_validation` vrací
`ZHITResult.Z_fit`), chyběl jen přepínač. Nejlevnější položka v seznamu -
matematika žádná nová, jen protažení už spočítaného pole přes frekvenční filtr
do fitu.

## 8. Architektonická myšlenka: element, který bere podobvod

Porous Electrode (Göhr / De Levie) a Surface Relaxation Impedance nejsou
v Zahnerově modelu obyčejné prvky - jsou to **transformátory topologie**.
PE vezme libovolný podobvod `q` (impedance stěny póru) a zabalí ho do
přenosového vedení s `Z_p` (elektrolyt v póru) a `Z_s` (pevná fáze):

```
Z_S = Z_par + Z* * [C + (1-C)*2ps + S*(p²*q_n + s²*q_0)]
                 / [S*(1 + q_n*q_0) + C*(q_n + q_0)]

Z*  = sqrt((Z_p + Z_s) * Z_q)      Z_par = Z_p*Z_s/(Z_p + Z_s)
C   = cosh((Z_p + Z_s)/Z*)         S     = sinh((Z_p + Z_s)/Z*)
p   = Z_p/(Z_p + Z_s)              s     = 1 - p
q_0 = Z*/Z_0                       q_n   = Z*/Z_n
```

Zdroj: H. Göhr et al., *Kinetic Properties of Smooth and Porous Lead / Lead
Sulfate Electrodes*, 34th I.S.E. Meeting, Erlangen (1983). Dokument dodává, že
podle L. Bay a Key West, *Solar Energy Materials and Solar Cells* 87 (2005)
613-628 platí model i pro prostorově distribuovanou elektronovou a iontovou
vodivost ve formě žebříkové sítě (pak se Z_o a Z_n vynechají).

**Proč je to zajímavé:** naše elementová algebra (`-`, `|`) dnes zná jen listy
a spojení. Element, jehož *parametrem je jiný obvod*, je skutečný krok v návrhu
DSL a otevře i další distribuované modely. Zjednodušená De Levie forma
(uniformní póry, Z_o a Z_n vynechány) je přitom docela krotký vzorec.

Aplikace podle dokumentu: baterie, palivové články, **povrchové vrstvy
korodujících elektrod, nehomogenní oxidové vrstvy** - opět naše doména.

Surface Relaxation Impedance je stejná architektonická třída (aplikuje se na
přenos náboje, ne samostatně), ale exotičtější - beru ji jen jako důkaz, že se
ta abstrakce vyplatí.

---

## Co nebrat

**User Element (JavaScript + complex.js).** V Pythonu je odpověď "podědit
`CircuitElement`", což je lepší než skriptovací okno. Z CLI stringu by to
navíc znamenalo `eval`. Vyřešeno, jen elegantněji.

**Start values, fixování parametrů, "zafixuj pár, fitni, odfixuj".**
Máme DE, multistart, `bounds.py` a fixování přes string. Tady jsme napřed.

**Tafel, Butler-Volmer, CV, fill factor.** Jiná metoda, jiný projekt.

Jediné, co bych si z non-EIS části vzal, je koncept **metadat datasetu**
(referenční elektroda, iR-drop korekce, plocha elektrody), která se propagují
do normalizovaných os (Ω*cm²). U nás `--area` teče jen do oxidové analýzy.

---

## Doporučené pořadí

1. ~~**Significance**~~ - HOTOVO v 0.33.0 (zbývá jen plot S_i(f))
2. **Blokující difuze (coth)** - triviální doplnění mřížky, kterou máme skoro celou
3. ~~**Young-Göhr**~~ - HOTOVO v 0.37.0
4. ~~**Fit na Z-HIT rekonstrukci**~~ - HOTOVO v 0.34.0
5. **Rovina komplexní kapacity** - vizualizace pro `CC`, které chybí
6. **Log-polární residuum** - navazuje na existující srovnání váhování
7. **Series Fit** - největší práce, největší změna použitelnosti
8. **Element beroucí podobvod** - architektonický krok, až bude poptávka

---

## Zdroj

Zahner Analysis, uživatelský manuál, 11/2023, 35 stran.
Reference citované dokumentem a relevantní pro nás:

- C. H. Hsu, F. Mansfeld, *Corrosion* **57**/9 (2001) 747-748 - CPE → C_eff (už používáme)
- Brug et al., *J. Electroanal. Chem.* **176** (1984) 275 - alternativní CPE → C_eff (už používáme)
- H. Göhr et al., 34th I.S.E. Meeting, Erlangen (1983) - porézní elektroda
- L. Bay, Key West, *Sol. Energy Mater. Sol. Cells* **87** (2005) 613-628 - žebříková síť
- *J. Electrochem. Soc.* **150**/3 (2003) A292-A300 - komplexní kapacita u superkondenzátorů
