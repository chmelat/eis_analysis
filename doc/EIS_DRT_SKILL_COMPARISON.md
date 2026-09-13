# Srovnani: eis_analysis vs eis-drt-batch-analysis

Datum: 2026-09-13
Porovnavane verze: eis_analysis v0.35.0, eis-drt-batch-analysis v1.0.0
(2026-09-11, adresar `tmp/eis-drt-batch-analysis-skill-main`)

## Co je predmetem srovnani

Verejny zdrojovy kod (GPL-3.0-or-later), autor pod prezdivkou R3ze959. Snapshot
obsahuje kompletni Python kod (~12 000 radku ve 39 skriptech), osm referencnich
dokumentu, testy a licencni soubory tretich stran. Algoritmy lze tedy overit
ctenim kodu, na rozdil od [ZScope](ZSCOPE_COMPARISON.md).

Dokumentace repozitare je cinsky; `SKILL.md` a `references/` jsou anglicky.
Overeno bylo: `SKILL.md`, `references/methodology.md`, oba README, seznam
importu, seznam funkci `advanced_drt.py` a CLI argumenty obou vstupnich bodu.
Radek po radku byl cten `advanced_drt.py` jen castecne a `batch_drt.py`
(2 910 radku) neni cten cely -- tvrzeni o detailech jeho vnitrni logiky nize
vychazi z dokumentace a signatur, ne z auditu kodu.

**Zasadni rozdil v zamereni.** eis_analysis je **interaktivni nastroj na jedno
spektrum**: nacti soubor, zvaliduj, sprav DRT, nafituj obvod, ukaz grafy,
vypis vysledky do konzole. eis-drt-batch-analysis je **davkovy pipeline a
auditni vrstva**: nacti adresar, kazde spektrum zpracuj, vsechno zapis do
strukturovaneho stromu vystupu s hashi zdroju a otiskem behu, a nic
neinterpretuj. Je to zaroven Codex Skill -- balicek instrukci pro LLM agenta --
takze velka cast repozitare je text pro model, ne kod.

Ten rozdil vysvetluje vetsinu neshod nize.

---

## Prehled

| | eis_analysis | eis-drt-batch-analysis |
|---|---|---|
| Typ | CLI + Python knihovna | davkove CLI + Codex Skill |
| Licence | MIT | **GPL-3.0-or-later** |
| Rozsah kodu | ~16 500 radku, 63 modulu | ~12 000 radku, 39 skriptu |
| Testy | 575 testu (pytest) | 329 testu (unittest) |
| Zavislosti | numpy, scipy, matplotlib | + pyimpspec, CVXOPT, pandas, openpyxl |
| Python | >= 3.9 | **pouze 3.11** (vlastni `.venv`) |
| Instalace | `pip install -e .` | `bootstrap.py` do izolovaneho venv (~477 MiB) |
| Overene platformy | Linux (vyvoj), prenositelne | **jen macOS arm64**; Win/Linux neoveren |
| Vstup | Gamry `.DTA`, CSV | text + pyimpspec: `.dta .mpt .z .idf .ids .dfr .P00 .i2b .pssession .ods .xlsx` |
| Davka | ne, jedno spektrum na spusteni | ano, adresare + manifest |
| Export dat | **zadny** (jen grafy a konzole) | CSV, JSON, Origin-ready, dlouhe tabulky |
| Vlastni numerika DRT | ano, cela | ne, primarni DRT je pyimpspec |

---

## Puvod numeriky: kdo co pocita sam

Toto je klic k celemu srovnani.

**eis-drt-batch-analysis** nepocita primarni DRT sam. `calculate_drt_tr_rbf`,
`calculate_drt_tr_nnls`, `calculate_drt_lm` (Loewner), Lin-KK, Z-HIT i vsechny
parsery vendor formatu jsou volani do **pyimpspec 5.1.3**. Vlastni je az to,
co je nad tim: vyber lambdy, rozhodovani o serioveho L, klasifikace piku,
audit stability, export a reprodukovatelnost. Vlastni numerika je jen ve
`advanced_drt.py` -- znamenkove GDRT, difuzni DDT kandidati a "model-and-reduce"
-- a ta stoji na `scipy.optimize.lsq_linear` a `least_squares`.

**eis_analysis** nema pyimpspec ani zadnou EIS zavislost. DRT, Lin-KK, Z-HIT,
fitovani obvodu, Voigt retezec -- vse je napsane v projektu proti holemu
numpy/scipy.

Dusledky:

- Jejich licence **musi** byt GPL, protoze pyimpspec je GPL-3.0. Kazdy, kdo by
  chtel jejich kod pouzit, dedi GPL. Nas MIT je v tomhle nesrovnatelne volnejsi
  a je to duvod, proc odtud nelze kopirovat kod, jen napady.
- Oni dostali TR-RBF, Loewner a sadu vendor parseru zdarma; my mame kontrolu
  nad kazdym radkem a zadnou 477 MiB instalaci.
- Overitelnost je asymetricka: jejich jadro je overene tim, ze je overeny
  pyimpspec; nase jadro je overene jen nasimi testy.

---

## DRT: srovnani po castech

| | eis_analysis | eis-drt-batch-analysis |
|---|---|---|
| Primarni metoda | vlastni Tikhonov + NNLS | pyimpspec TR-RBF |
| Baze | po castech konstantni (kolokace), `n_tau=100` | gaussovske RBF, FWHM koeficient 0.5 |
| Penalizace | **2. derivace** | **1. derivace** (prepinatelne) |
| Nezapornost | ano (NNLS) | ano |
| Rozsah tau | z merici okna, `1/(2*pi*f)` | tez, ale RBF kolokace na frekvencich |
| Vyber lambda | hybrid GCV + L-krivka, mrizka 1e-5..1, 20 bodu | konsenzus: mGCV ze 3 startu, shoda do 0.25 dekady; jinak roh L-krivky; mrizka 1e-7..1e-1, 13 bodu |
| Skalovani | zadne | deleni `median(abs(Z))` pred solverem, zpetne obnoveni |
| Seriova indukcnost v DRT | **ne**, jen varovani na indukcni body | adaptivni souboj modelu bez L / se seriovym L |
| Kontrola nezavislym algoritmem | ne | ano, TR-NNLS jako druhy nazor |
| Stabilita v lambda | 4 sondy: `10^(+-0.5)`, `10^(+-1)` | 2 sondy: `/10`, `*10` |
| Trida piku | 3 (`stable`/`marginal`/`artifact`) | 7 (vc. `boundary-sensitive`, `unsupported-outside-window`, `inductive-overlap`, `split-merge-ambiguous`) |
| Plocha piku | integral povodi, delici body v udolich | totez, plus zvlast "podporena" cast v mericim okne |
| Detekce piku | `scipy.find_peaks` nebo vazeny GMM s BIC | lokalni maxima |
| Pasy nejistoty | ne | ano (pyimpspec podminena posterior 0.5-99.5 %) |
| Znamenkove DRT (RL vetve) | ne | ano, vlastni regularizovane GDRT |
| Loewner RC/RL | ne | ano (pyimpspec) |
| Difuzni distribuce (DDT) | ne | ano, 3 kandidati + Warburg baseline |

### Kde se nezavisle shodujeme

Dve veci stoji za zminku, protoze k nim oba projekty dosly nezavisle:

1. **Plocha povodi, ne vyska piku.** Obe implementace deli osu tau v udolich
   mezi piky a integruji `gamma` pres `ln(tau)`. Jejich `methodology.md` to
   zduvodnuje stejne jako nas docstring v `estimation.py:32`: vyska je
   ovlivnena sirenim, prekryvem a regularizaci.
2. **Lambda neni jedno cislo, ale rozhodovaci cesta.** Oba projekty zaznamenavaji
   GCV i L-krivku, oba detekuji, ze optimum sedi na okraji rozsahu, a oba
   odmitaji okrajovou hodnotu prohlasit za optimum.

### Kde jsou vecne napred

**Seriova indukcnost uvnitr DRT.** Nas `linear_system.py:163-171` indukcni
body jen spocita a vyda varovani; pak je posle do NNLS, kde je nezaporna RC
distribuce nemuze reprezentovat. Oni resi dva modely (bez L a se seriovym L)
na trech lambdach a vyberou podle rekonstrukce plus regrese `-Im(Z)` proti
`omega` ve vysokych dekadach (pozadavek: kladne L a R^2 >= 0.90). To je
konkretni, prenositelny recept -- a nas `rinf_estimation/rlk_fit.py` uz kus te
matematiky ma, jen ji nepropojuje do DRT.

**Klasifikace okrajovych piku.** Jejich trida `boundary-sensitive` -- pik blize
nez 0.7 dekady k okraji mericiho okna -- resi slabinu, kterou nas kod mel:
`_estimate_peak_resistance` vede krajni povodi az na konec pole, takze
posledni pik pohlti i to, co tam navrsila regularizace.

Jejich druha trida `unsupported-outside-window` je pro nas bezpredmetna. Nase
mrizka tau je **presne** merici okno (`linear_system.py:59-62`:
`tau_min = 1/(2*pi*f_max)`, `tau_max = 1/(2*pi*f_min)`), takze pik mimo okno u
nas vzniknout nemuze. Misto nej ale nastane jina vec, kterou oni pojmenovanou
nemaji: nezaporne NNLS nema kam dat odezvu s casovou konstantou za oknem, tak
ji navrsi do krajniho binu. `scipy.find_peaks` tam pik vratit nedokaze
(krajni index nema souseda), takze to zustane bez povsimnuti a jen nafoukne
sousedni pik.

Oboji je opraveno ve v0.36.0: `boundary_sensitive` na kazdem piku a
`edge_bin_R_fraction` jako podil R_pol v krajnim binu.

**Tiche orezani sondy lambda.** V nasem `stability.py:203` se sondy orezavaji
na `[1e-6, 1]` a pak deduplikuji -- bez varovani. Pri `lambda* = 0.5` se obe
horni sondy slozi na 1.0, zbydou tri body misto ctyr a verdikt `stable` se
vyda pres uzsi rozsah, nez rika napoveda `--lambda-probe`. Jejich pravidlo
("orezany sweep neni ekvivalent plne +-1 dekady a je explicitne
`boundary-sensitive`") je jednoradkova oprava, kterou stoji za to prevzit.

### Kde jsme napred my

**Hustota sondovani lambda.** Ctyri sondy vcetne pulkrokovych `10^(+-0.5)`
rozlisi pik, ktery zmizi uz pri pulce dekady, od piku, ktery prezije celou.
Jejich dve sondy tohle nerozlisi.

**Detekce piku.** Vazena smes gaussianu s vyberem poctu komponent pres BIC
(`drt/peaks.py`) je vecne robustnejsi nez hledani lokalnich maxim, hlavne u
prekryvajicich se procesu. Oni na prekryv reaguji az ex post tridou
`split-merge-ambiguous`.

**Penalizace 2. derivace.** Pro hladke relaxacni spektrum je to obvykle
vhodnejsi volba nez jejich vychozi 1. derivace; ta ma sklon davat schodovite
reseni. Oni maji prepinac `--derivative-order`, my mame pevne 2.

---

## Validace: Lin-KK a Z-HIT

Oba projekty maji obe metody. Rozdil je v roli a v implementaci.

| | eis_analysis | eis-drt-batch-analysis |
|---|---|---|
| Lin-KK | vlastni implementace, vychozi zapnuta | pyimpspec, vychozi **vypnuta** (`--kk off`) |
| Automatika M | `mu` kriterium, prah 0.85 | pyimpspec |
| Rozsireni tau | `--auto-extend` minimalizuje pseudo chi^2 | neni |
| Seriova C | `--kk-series-c` (Schonleber add_cap) | neni |
| Z-HIT | vlastni, vychozi zapnuty | pyimpspec, vychozi vypnuty |
| Fit na rekonstrukci | `--fit-on zhit/all` | neni |
| Hranice tvrzeni | zminena v dokumentaci | explicitni pravidlo v `methodology.md` |

Nase Lin-KK je propracovanejsi jako **metoda** (auto-extend, seriova C, moznost
fitovat na Z-HIT rekonstrukci). Jejich je propracovanejsi jako **rezim**: KK a
Z-HIT maji oddelene stavy od kvality vstupu, nizke reziduum nesmi prepsat
varovani o vstupnich datech, a vysledek "nelze posoudit" je platny vystup, ne
prochazi/neprochazi.

Jejich formulace "KK pass neni dukaz linearity" s citaci Urquidi-Macdonald 1990
je presne to, co nas [VALIDATION_METHOD_COMPARISON.md](VALIDATION_METHOD_COMPARISON.md)
rika taky. Tady neshoda neni.

---

## Co ma jen jeden z nich

### Jen eis_analysis

- **Fitovani ekvivalentnich obvodu.** Cely modul `fitting/` -- parser vyrazu
  `R()-(R()|Q())`, analyticky jakobian, diferencialni evoluce, multistart,
  kovariance a konfidencni intervaly, AIC/BIC zebricek kandidatu,
  `auto_suggest` na navrh obvodu z dat. Oni to **zamerne nedelaji**: SKILL.md
  explicitne zakazuje "physical equivalent-circuit fitting".
- **Voigt retezec** linearni regresi, s `mu` auto-M a variantami real/imag/complex.
- **Odhad R_inf** vcetne RLK fitu (`rinf_estimation/`).
- **Oxidova analyza** (`analysis/oxide.py`) -- tloustka a permitivita vrstvy,
  vcetne DQ elementu z poslednich commitu.
- **Detekce odlehlych bodu** (`validation/outliers.py`).
- **OCV vizualizace**.

### Jen eis-drt-batch-analysis

- **Davkove zpracovani** adresaru s manifestem jako allowlistem.
- **Strukturovany vystup** do sedmi adresaru (`00_overview` az
  `06_reproducibility`) se stabilnim `spectrum_uid` na spojovani vysledku.
- **Otisk behu** -- hash zdroju, zamcene verze zavislosti, identita
  interpretu; `--resume` jede jen pri presne shode, `--export-only`
  pregeneruje vystupy bez prepocitavani.
- **Znamenkove GDRT** (kladna gamma = RC, zaporna = RL) pro indukcni smycky,
  ktere seriove L nepokryva.
- **Loewner RC/RL** a **difuzni DDT kandidati** (blocking `coth(s)/s`,
  transmissive `tanh(s)/s`, Gerischer) s explicitnimi prijimacimi prahy.
- **Model-and-reduce** -- odecteni rusiveho prvku s kontrolou, ze se piky
  po odecteni neposunuly (limit 0.35 dekady / 50 % plochy).
- **Kontrola experimentalni platnosti z metadat** -- opakovana mereni,
  amplitudova linearita, ustaleni po klidu. Rozlisuje "nepodarilo se overit"
  od "neproslo".
- **Screening driftu behem sweepu** -- oddeluje monotonni trend residui podle
  poradi sberu od nahodneho sumu.
- **3D vodopadove grafy a heatmapy** pres napetove uzly, na vyzadani.

---

## Davkove zpracovani a export: nejvetsi rozdil

Nas projekt **neexportuje zadna data**. Zadny `to_csv`, zadny `json.dump`,
zadny `savetxt` nikde v `eis_analysis/`. Vystup je matplotlib okno a text v
konzoli. Kdo chce cisla dal zpracovat, musi pouzit Python API a napsat si
export sam.

To je pro interaktivni praci s jednim spektrem v poradku a zamerne to sedi s
nasi CLI filozofii (`doc/CLI_OUTPUT_UNIFICATION.md`: knihovna vraci
`*Result` dataclassy, tiskne jen `cli/handlers/`). Ale znamena to, ze
srovnani deseti spekter mezi sebou je u nas rucni prace.

Jejich reseni je druhy extrem: grafy jsou az posledni krok, hlavni produkt jsou
CSV a JSON, a prezentacni grafy se generuji **jen na vyzadani** (`--trend-plots
off` je vychozi). Diagnosticke grafy automaticke jsou; selhani vykresleni
nikdy nesmi smazat ciselny vysledek -- JSON se zapisuje pred rendrovanim.

Ta posledni vec je dobry napad nezavisly na davce: **checkpoint ciselneho
vysledku pred vykreslenim**. U nas vyjimka v `visualization/` shodi cely beh
a spoctena DRT se ztrati.

---

## Nesrovnalosti a slaba mista v jejich releasu

- **Overeno jen na macOS arm64.** README to prizna: Windows a Linux maji
  vstupni body, ale 1.0 na nich overen nebyl. Pri 477 MiB instalaci a vazbe na
  presne CPython 3.11 to neni maly zavazek.
- **Cislo 329 testu neni srovnatelne s nasimi 575.** Jejich testy jsou z velke
  casti kontrakty exportu, privacy scan a chovani resume; vlastni numeriky
  testuji mene, protoze numerika je z pyimpspec.
- **Dokumentace je rozdvojena.** Cinsky README pro cloveka, anglicky SKILL.md
  pro model. Obsah se casto prekryva -- stejny problem s jedinym zdrojem
  pravdy, jaky nas CLAUDE.md zakazuje.
- **SKILL.md je z velke casti prompt, ne specifikace.** Instrukce typu
  "Wait for the answer; no answer is not permission" jsou pravidla pro LLM
  agenta. Kdo pouziva jen Python CLI, cte je zbytecne.
- **Verze v SKILL.md nesedi.** Text uvadi "Release 1.0.0", ale dal se odkazuje
  na chovani "0.2.1" a "0.2.2" jako na aktualni pravidla.

---

## Co si odtud vzit

Serazeno podle pomeru uzitku k praci:

1. **Varovat, kdyz se sondy lambda orezou.** ~~Hotovo ve v0.36.0.~~
   `stability.py:203` orezavalo tise.
   Jedno pole navic v `StabilityDiagnostics` a jeden radek do `warnings`.
   Bez teto zmeny muze `--lambda-probe` hlasit `stable` na zaklade tri sond
   v uzsim rozsahu, nez slibuje napoveda.
2. **Oznacit piky u okraje mericiho okna.** ~~Hotovo ve v0.36.0.~~ Jejich
   pravidlo 0.7 dekady je konzervativni heuristika, ne konstanta -- ale nemit
   zadne takove oznaceni je horsi nez mit heuristiku a priznat ji. K tomu
   navic hlaseni podilu R_pol v krajnim binu, ktere zachyti navrseni odezvy
   zpoza okna i tam, kde zadny pik detekovan neni.
3. **Checkpoint DRT vysledku pred vykreslenim.** Vyjimka v grafech nema
   zahodit spoctenou analyzu.
4. **Seriove L jako soubezny model v DRT.** Nejvetsi vecny prinos, ale i
   nejvic prace: vyzaduje rozsirit navrhovou matici a rozhodovaci logiku.
   Kus matematiky uz je v `rinf_estimation/rlk_fit.py`.
5. **Export do CSV/JSON.** Ne cely jejich sedmiadresarovy strom -- ale jeden
   prepinac `--export-csv`, ktery zapise tau, gamma, piky a residua, by
   odstranil nejvetsi prakticke omezeni naseho CLI.

Nebrat si: davkovou architekturu, otisky behu a resume (resi problem, ktery
nemame), Codex Skill vrstvu, ani zavislost na pyimpspec -- ta by nam vymenila
MIT za GPL a pridala 400 MiB instalace za funkce, ktere z 80 % uz mame vlastni.

---

## Zaver

Nejsou to konkurenti, jsou to dve poloviny jineho problemu. eis_analysis je
**hlubsi na jednom spektru** -- fitovani obvodu, Voigt retezec, oxidova
analyza a propracovanejsi Lin-KK nemaji u nich protejsek a mit ho nemaji
zamerne. eis-drt-batch-analysis je **siroky pres davku** -- reprodukovatelnost,
export, znamenkove a difuzni vetve a oddeleni "spocitano" od "overeno".

Nejcennejsi na jejich repozitari neni kod (ten je stejne GPL a z velke casti
pyimpspec), ale `references/methodology.md`. Je to poctivy, ocitovany vycet
toho, co DRT **netvrdi** -- a ctyri z peti bodu v seznamu vyse jsou prime
dusledky jeho cteni.
