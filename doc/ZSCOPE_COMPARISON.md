# Srovnani: eis_analysis vs ZScope

Datum: 2026-08-30
Porovnavane verze: eis_analysis v0.25.4, ZScope v2.2.0 (adresar `tmp/ZScope-main`)

## Co je predmetem srovnani

Adresar `tmp/ZScope-main` **neobsahuje zdrojovy kod**. Je to snapshot GitHub
Pages vetve projektu ZScope: sablona HTML5UP "Future Imperfect"
(`index.html`, `single.html`, CSS/JS/webfonty), screenshoty, `README.md`,
`CITATION.cff` a adresar `benchmarks/`. Celkem 94 souboru / 12 MB, z toho
zhruba 11 MB obrazky a fonty. Zadny Python.

ZScope je zamerne closed-source ("Source code is not publicly distributed"),
distribuovany jako instalator pro Windows. Autor: Tecush Mohammadi,
DOI 10.5281/zenodo.20357547.

**Licence:** jediny licencni soubor ve snapshotu je `LICENSE.txt` s textem
Creative Commons Attribution 3.0 -- to je licence pouzite webove sablony
HTML5UP, ne projektu. README ma odznak "MIT License" odkazujici na "LICENSE",
takze licence benchmarkovych dat neni jednoznacne urcena.

**Dusledek pro toto srovnani:** algoritmy ZScope nelze overit ctenim kodu.
Srovnani stoji na deklarovanych vlastnostech z README a na dodanych
benchmarkovych datech, ktera jsou jedina nezavisle overitelna cast.

---

## Prehled

| | eis_analysis | ZScope |
|---|---|---|
| Typ | CLI + Python knihovna | desktopova GUI aplikace |
| Platformy | Linux, macOS, Windows (Python >= 3.9) | pouze Windows; macOS/Linux "planned" |
| Zdrojovy kod | verejny (MIT) | neverejny |
| Rozsah | ~13 900 radku, modularni balik | neznamy |
| Testy | 31 testovacich souboru (pytest) | neverejne; publikovan jen vysledny benchmark |
| Zavislosti | numpy, scipy, matplotlib | irelevantni (standalone instalator) |
| Instalace | `pip install -e .` | instalator, bez Pythonu |
| Cena | zdarma | zdarma |

Filozoficky rozdil je stejneho druhu jako u [ADAPT-ECM](ADAPT_ECM_COMPARISON.md):
eis_analysis mira na skriptovatelnost, reprodukovatelnost a pouzitelnost jako
knihovna; ZScope na interaktivni praci s jednim spektrem pred ocima uzivatele.

---

## Metody

| Metoda | eis_analysis | ZScope |
|---|---|---|
| Lin-KK validace | ano (auto `extend_decades`, mu metrika, volitelna seriova C) | ano (Voigt baze, modulus vahy, residual mapa) |
| Z-HIT validace | ano (vc. optimalizace offsetu) | ne |
| DRT | Tikhonov, auto-lambda pres GCV + L-curve | ano (metoda volby lambda nedokumentovana) |
| Stabilita DRT piku | `--lambda-probe` (stable/marginal/artifact) | ne; publikovan rozlisovaci limit tau2/tau1 |
| Detekce piku | scipy nebo GMM s BIC | nedokumentovano |
| Voigt chain (auto-M) | ano, mu kriterium | ne (Voigt jen uvnitr KK) |
| Odhad R_inf | ano (RLK fit) | ne |
| Fitovani obvodu | DE (default) / multistart / single + TRF | LHS multi-start -> TRF, DE jako eskalace |
| Analyticky Jacobian | ano | ano (adjoint sensitivity) |
| Vahovani | uniform / sqrt / modulus / proportional | modulus + robustni soft-L1 |
| Nejistoty | kovariancni matice, 95% CI, condition number | **MCMC (emcee)**, kredibilni intervaly, R-hat, autokorelace, robustni sandwich kovariance |
| Vyber modelu | ne | **AIC/BIC** |
| Warm-start pro serie | ne | ano |
| Analyza oxidicke vrstvy | ano (tloustka <-> permitivita) | ne |

Obe implementace se v jadru fitovani prekryvaji vice, nez by README ZScope
naznacoval: analyticky Jacobian, DE i modulus vahovani mame take. Skutecne
rozdily jsou MCMC, AIC/BIC a warm-start na jejich strane; Z-HIT, auto-lambda
DRT, Voigt chain s auto-M a odhad R_inf na nasi.

---

## Prvky obvodu

**eis_analysis** (parser vyrazu, bez preddefinovanych obvodu):
R, C, L, G (vodivost), Q (CPE), W, Wo, K (Voigt v parametrizaci R-tau),
GE (Gerischer),
CC (Cole-Cole v rovine permitivity). Syntaxe `R(100) - (R(5000) | C(1e-6))`.

**ZScope** (kanvas + presety): R, C, L, CPE, W, FLW, FSW, Gerischer,
prenosove vedeni (porezni elektrody), plus **uzivatelske prvky definovane
v GUI** libovolnym vyrazem Z(omega) nebo Y(omega), exportovatelne jako JSON.

Chybi nam FLW/FSW (konecna difuze) a prenosove vedeni. Naopak CC element
(depresovany oblouk v rovine komplexni kapacity) v jejich seznamu neni.

---

## Import dat

| | eis_analysis | ZScope |
|---|---|---|
| Formaty | Gamry DTA (nativni parser, ZCURVE, metadata), CSV | Gamry, BioLogic, Autolab, Zahner, Ivium, CHI, PalmSens, PAR + genericky parser |
| CSV | auto oddelovac (`,` `;` tab), US i EU desetinna carka | totez, vc. locale detekce pri vlozeni ze schranky |
| Davkove nacteni | ne (resi se shellem: `for f in *.DTA`) | ano, s revizi mapovani sloupcu |
| Vlozeni ze schranky | ne | ano (Ctrl+V) |
| Perzistence relace | ne | `.zscope` projekt (ZIP + JSON manifest) |
| Sledovani puvodu dat | ne | ano (auto / rucni mapovani / rucni editace) |

---

## Overeni na benchmarkovych datech ZScope

Adresar `benchmarks/zview_test_data/` obsahuje 15 CSV se znamou pravdou
(4 obvody + jedna degenerovana sada, kazdy ve variante clean / 2 % / 5 % sum).
Nas loader je cte bez uprav.

**Metodika overeni:** ground truth ziskan fitem na `clean` variante (na
bezsumnych datech DE dosahuje presne parametry), nasledne fit na zasumenych
datech a porovnani. DE se `seed=0`, vychozi nastaveni.

| Obvod | topologie | nejhorsi param @2 % | nejhorsi param @5 % |
|---|---|---|---|
| Randles | `R-(R\|C)` | C 0,27 % | R_s 0,54 % |
| Randles+Warburg | `R-((R-W)\|C)` | C 0,40 % | C 0,89 % |
| CPE Randles | `R-(R\|Q)` | Q 1,75 % | Q 4,35 % |
| Two-TimeConstants | `R-(R\|C)-(R\|C)` | viz nize | C2 5,68 % |

Publikovane hodnoty ZScope pro 2 % sum: C_dl 0,04 %, sigma_W 0,04 %,
Q 0,99 %, C2 1,48 %. **Jsme ve stejnem radu**, u CPE mirne horsi, jinak
srovnatelne. Doby fitu 0,6-2,0 s vs. jejich 0,7-2,6 s (jine HW, jen orientacne).

**Nalezena slabina:** u `Two_TimeConstants` @2 % nasel DE reseni s prohozenymi
RC vetvemi (R1<->R2, C1<->C2). RMSE 2,04 % je v poradku, jde o degeneraci
znaceni, ne o spatny fit -- ale ukazuje, ze u dvou blizkych casovych konstant
nemame zadne poradove omezeni ani detekci teto permutace.

---

## Nesrovnalosti v materialech ZScope

Zaznamenano pro uplnost, protoze cast srovnani stoji na jejich cislech:

1. **`CITATION.cff` uvadi `version: v1.2.0`, `date-released: 2026-06-02`**,
   zatimco README i prilozeny BibTeX uvadeji 2.2.0 / 2026. Citace generovana
   z CFF tedy odkaze na spatnou verzi.
2. **README slibuje 5 referencnich obvodu, 3 modely sumu a >19 000 fitu**
   na jeden benchmarkovy beh; dodany `benchmarks/` obsahuje 4 obvody x 3 stavy
   = 12 fitu, jeden seed (`n_iter` 1-4). Pata sada
   (`Two_TimeConstants_Degenerate`) je v datech, ale ne ve vysledcich.
3. **"analysis scripts are available in `benchmarks/`"** -- v tomto snapshotu
   zadne skripty nejsou. Nezavisle zopakovani benchmarku tedy neni mozne.
4. **Dve ruzne definice chi2.** `benchmarks/reports/README.md` uvadi pravidlo
   "chi2_reduced ~ 1,0 = idealni fit" a hned pod nim hodnoty 1e-4 az 1e-3
   (zjevne nenormalizovane na sigma). README projektu uvadi "chi2_calibrated"
   0,992 / 1,006. Ktera definice patri ke ktere tabulce, dokumentovano neni.

Sama dokumentace ZScope je jinak nadprumerne poctiva v tom, ze publikuje
i limity metod (rozlisovaci limit DRT, pokryti CI pri kontaminovanych datech,
korelace Q-alfa u CPE) -- to je prakticky sdelnejsi nez vetsina srovnatelnych
nastroju.

---

## Prevzato: regresni benchmark

Referencni obvody jsou od v0.25.5 soucasti testu jako
`tests/test_zscope_benchmark.py` (20 testu, ~20 s). Pro kazdy ze ctyr obvodu
se overuje:

1. **presna rekonstrukce z bezsumnych dat** (< 0,1 %) -- nejsilnejsi
   regresni signal, bez statisticke vule,
2. **recovery parametru** pri 2 % a 5 % proporcionalnim sumu v kalibrovanem
   pasmu,
3. **dosazeni sumove podlahy** -- relativni chyba fitu mezi 1,0x a 1,4x
   urovne sumu (namereno stabilne 1,19-1,20x).

Data se v testu **generuji znovu**, nevendoruji se. Duvody: nejasna licence
(viz vyse) a to, ze test ukazujici do neverzovaneho `tmp/` by se vsude jinde
nez na jednom stroji preskocil. Ground truth byl overen proti jejich CSV na
max|dZ|/|Z| ~ 5e-7, coz je zaokrouhleni jejich 7-mistneho exportu.

Sumovy model testu (proporcionalni Gauss, sigma = uroven x |Z| zvlast na
realnou a imaginarni cast) neni bitove shodny s jejich generatorem -- ten je
uzsi a nedokumentovany -- proto jsou tolerance kalibrovany z vlastnich behu,
ne prevzaty z jejich tabulek. Permutace RC vetvi je osetrena serazenim podle
tau pred porovnanim.

---

### AIC/BIC pro vyber obvodu

Od v0.26.0 lze `--circuit` zadat vicekrat; kandidati se nafituji na tychz
datech se stejnym vazenim a seradi se v jedne tabulce podle BIC. Doplneno
`compute_information_criteria` v `fitting/diagnostics.py` a verejne pole
`FitResult.n_free_params`.

Voigtuv retez do zebricku zamerne nevstupuje: jeho prvky `K(R, tau)` lezi na
pevne mrizce tau, ktera se neoptimalizuje, takze naivni pocet parametru jeho
slozitost nadhodnocuje o M (pri M = 7 asi 35 bodu BIC, dost na prevraceni
poradi) a ferovy pocet by stejne byl jen dolni odhad -- mrizka i
`prune_threshold` jsou datove informovane.

---

## Co stoji za prevzeti

| Vlastnost | Priorita | Slozitost | Poznamka |
|---|---|---|---|
| Robustni loss (soft-L1) + sandwich kovariance | vysoka | stredni | jejich mereni: pokryti CI 89 % -> 94 % pri 5 % kontaminovanych bodu |
| Dokumentovat rozlisovaci limit DRT | stredni | nizka | tau2/tau1 ~ 15 x (sum v %); mame `--lambda-probe`, ale zadne cislo |
| FLW / FSW (konecna difuze) | stredni | nizka | doplneni do `circuit_elements/distributed.py` |
| Warm-start pro serie mereni | nizka | stredni | dava smysl az s davkovym rezimem |

Mimo zamer projektu: MCMC posteriory, interaktivni kanvas a `.zscope`
projektove soubory smeruji k jinemu typu nastroje, nez je skriptovatelna
CLI knihovna.

---

## Reference

**ZScope:**
- GitHub: https://github.com/Tecush/ZScope
- DOI: https://doi.org/10.5281/zenodo.20357547
- Mohammadi, T. (2026). *ZScope: Publication-Grade Electrochemical Impedance
  Spectroscopy Analysis Platform*, v2.2.0

**eis_analysis:**
- GitHub: https://github.com/chmelat/eis_analysis
