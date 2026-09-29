# Srovnani: eis_analysis vs eis-drt-batch-analysis

Datum: 2026-09-28 (aktualizace; puvodni verze 2026-09-13)
Porovnavane verze: eis_analysis v0.46.0, eis-drt-batch-analysis **v1.2.0rc1**
(2026-09-28, predbezne vydani, adresar `tmp/eis-drt-batch-analysis-skill-main`).
Puvodni srovnani bylo proti eis_analysis v0.35.0 a eis-drt-batch-analysis v1.0.0.

## Co je predmetem srovnani

Verejny zdrojovy kod (GPL-3.0-or-later), autor pod prezdivkou R3ze959. Snapshot
obsahuje kompletni Python kod (~15 300 radku v 50 skriptech), dvanact
referencnich dokumentu, testy a licencni soubory tretich stran. Algoritmy lze
tedy overit ctenim kodu, na rozdil od [ZScope](ZSCOPE_COMPARISON.md).

Dokumentace repozitare je cinsky; `SKILL.md` a `references/` jsou anglicky.
Pri aktualizaci bylo overeno: oba changelogy a release notes, `SKILL.md`,
`references/lambda-and-peak-interpretation.md`, cely novy `dual_lambda.py`,
CLI argumenty `batch_drt.py` s vychozimi hodnotami a jeho klasifikace piku.
`batch_drt.py` (2 918 radku) a `advanced_drt.py` (1 343 radku) nejsou ctene
cele -- tvrzeni o detailech jejich vnitrni logiky, ktera nejsou vyse
vyjmenovana, vychazi z dokumentace a signatur, ne z auditu kodu.

**Zasadni rozdil v zamereni.** eis_analysis je **interaktivni nastroj na jedno
spektrum**: nacti soubor, zvaliduj, sprav DRT, nafituj obvod, ukaz grafy,
vypis vysledky do konzole. eis-drt-batch-analysis je **davkovy pipeline a
auditni vrstva**: nacti adresar, kazde spektrum zpracuj, vsechno zapis do
strukturovaneho stromu vystupu s hashi zdroju a otiskem behu, a nic
neinterpretuj. Je to zaroven Codex Skill -- balicek instrukci pro LLM agenta --
takze velka cast repozitare je text pro model, ne kod.

Ten rozdil vysvetluje vetsinu neshod nize.

---

## Co se zmenilo od minuleho srovnani

### eis-drt-batch-analysis v1.0.0 -> v1.2.0rc1

- **Novy vychozi vyber lambda `gcv-lcurve`** (`scripts/dual_lambda.py`).
  Obycejne GCV a L-krivka se skutecnou penalizaci `x^T M x` se pocitaji
  nezavisle na stejne mrizce (1e-7..1e-1, **31 bodu**, drive 13). Prednost ma
  vnitrni minimum GCV, L-krivka je zaloha. Obe kandidatni krivky i jejich
  rekonstrukce se exportuji. Drivejsi konsenzus mGCV je uz jen explicitni
  legacy volba (`--lambda-policy consensus`).
- **Znamenkove GDRT a predikcni souboj pozitivni/znamenkove vetve jsou nove
  vypnute** (`--signed-gdrt off --predictive-fit off`). Drive se spoustely
  automaticky na indukcni nebo spatne nafitovana spektra. Loewner zustava
  zapnuty, DDT dal startuje na difuzni spoustec.
- **Hlavni `SKILL.md` zkracen o ~30 %** (137 radku); interakce s uzivatelem a
  grafy se presunuly do novych referenci nactenych jen v pripade potreby.
- **Nove zobrazovaci volby**: 2D vrstvene DRT, parove Nyquist/Bode srovnani,
  heatmapy `plot_drt_map.py` s volitelnou (jen zobrazovaci) interpolaci.
- **Metodicke texty** o vyberu modelu, mazani bodu, Gaussove rozkladu piku a
  zastaveni iteraci. Spolecne DRT s vice RL/RLC vetvemi je popsane jako
  postup, ale **neni implementovane** v CLI -- a text to otevrene rika.
- 390 testu (drive 329). Zavislosti beze zmeny, porad overeno jen na macOS
  arm64.

Primarni resic se nezmenil: pyimpspec TR-RBF, 1. derivace, FWHM koeficient
0.5, skalovani `median(|Z|)`, okraj 0.7 dekady, kriterium seriove L
R^2 >= 0.90.

### eis_analysis v0.35.0 -> v0.46.0 (jen to, co se tyka srovnani)

- v0.36.0: `boundary_sensitive` na pikach, `edge_bin_R_fraction`, varovani
  pri orezani sond lambda -- body 1 a 2 z puvodniho seznamu "Co si odtud vzit".
- v0.38.0: **vazeni datoveho clenu DRT**, vychozi `1/sqrt|Z|`
  (`--drt-weighting`).
- v0.39.0: **volitelne prodlouzeni mrizky tau** za pomaly konec okna
  (`--tau-extend DECADES|auto`), piky s `outside_window`,
  `R_pol_extrapolated_fraction`.
- v0.40.0-v0.41.0: **seriova indukcnost L v modelu DRT**
  (`--drt-inductance auto`, bod 4 puvodniho seznamu); R_inf fitem
  `R-L-(R|Q)`.
- v0.41.1: dostatek iteraci NNLS pro vazena spektra s velkym dynamickym
  rozsahem.

---

## Prehled

| | eis_analysis | eis-drt-batch-analysis |
|---|---|---|
| Typ | CLI + Python knihovna | davkove CLI + Codex Skill |
| Licence | MIT | **GPL-3.0-or-later** |
| Rozsah kodu | ~17 300 radku, 66 modulu | ~15 300 radku, 50 skriptu |
| Testy | 612 testu (592 ve vychozim behu, pytest) | 390 testu (unittest) |
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
co je nad tim: vyber lambdy, rozhodovani o seriovem L, klasifikace piku,
audit stability, export a reprodukovatelnost. Vlastni numerika je ve
`advanced_drt.py` -- znamenkove GDRT, difuzni DDT kandidati a "model-and-reduce"
-- a ta stoji na `scipy.optimize.lsq_linear` a `least_squares`.

Novy `dual_lambda.py` je hybrid: sklada matice TR-RBF privatnimi funkcemi
pyimpspec (`tr._assemble_A_matrix`, `tr._solve_qp_cvxopt`, ...) a GCV pocita
sam. Proto odmitne bezet s jinou verzi pyimpspec nez 5.1.3 -- rozumna pojistka,
ale zaroven dukaz, ze jejich kod je na vnitrnostech pyimpspec zavisly vic nez
drive.

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
| Primarni metoda | vlastni Tikhonov + NNLS | pyimpspec TR-RBF (kvadraticky program, CVXOPT) |
| Baze | po castech konstantni (kolokace), `n_tau=100` | gaussovske RBF, FWHM koeficient 0.5 |
| Penalizace | **2. derivace** | **1. derivace** (prepinatelne) |
| Nezapornost | ano (NNLS) | ano |
| Rozsah tau | merici okno `1/(2*pi*f)`; volitelne prodlouzeni za pomaly konec (`--tau-extend`, `auto`) | RBF centra na `1/f`, export orezan na merici okno |
| Vyber lambda | GCV na mrizce 1e-10..1e-2 (20 bodu), pak L-krivka +-1.5 dekady kolem nej; plati vetsi z obou lambd | GCV a L-krivka nezavisle na mrizce 1e-7..1e-1 (31 bodu); prednost ma vnitrni minimum GCV |
| Hlaseni neshody GCV/L-krivka | obe hodnoty v `LambdaSelection`; varovani, kdyz je roh L-krivky > 1 dekadu pod GCV | `requires_review` pri rozdilu > 1 dekada, plochem minimu GCV nebo okraji |
| Vazeni / skalovani | vazeni `1/sqrt|Z|` (vychozi), prenormovane na `||w*Z|| = ||Z||` | deleni `median(|Z|)` pred resicem -- skalovani, ne vazeni |
| Seriova indukcnost v DRT | ano, neregularizovany sloupec `j*omega*L`; `auto` pri Im(Z) > 0 v horni dekade | souboj modelu bez L / se L; kladne L a R^2 >= 0.90 regrese `-Im(Z)` proti `omega` |
| Kontrola nezavislym algoritmem | ne | ano, TR-NNLS jako druhy nazor |
| Stabilita v lambda | 4 sondy: `10^(+-0.5)`, `10^(+-1)`, varovani pri orezani | 2 sondy: `/10`, `*10`, orezani = `boundary-sensitive` |
| Trida piku | 3 verdikty (`stable`/`marginal`/`artifact`) + priznaky `boundary_sensitive`, `outside_window` | 8 trid (`robust-to-lambda`, `tentative`, `unstable`, `minor-low-signal`, `boundary-sensitive`, `unsupported-outside-window`, `inductive-overlap`, `split-merge-ambiguous`) |
| Plocha piku | integral povodi, delici body v udolich | totez, plus zvlast "podporena" cast v mericim okne |
| Detekce piku | `scipy.find_peaks` nebo vazeny GMM s BIC | lokalni maxima |
| Pasy nejistoty | ne | ano (pyimpspec podminena posterior 0.5-99.5 %) |
| Znamenkove DRT (RL vetve) | ne | ano, vlastni GDRT; **od v1.2 vypnute ve vychozim stavu** |
| Loewner RC/RL | ne | ano (pyimpspec), zapnuty |
| Difuzni distribuce (DDT) | ne | ano, 3 kandidati + Warburg baseline, na spoustec |

### Kde se nezavisle shodujeme

Tri veci stoji za zminku, protoze k nim oba projekty dosly nezavisle:

1. **Plocha povodi, ne vyska piku.** Obe implementace deli osu tau v udolich
   mezi piky a integruji `gamma` pres `ln(tau)`. Jejich `methodology.md` to
   zduvodnuje stejne jako nas docstring v `estimation.py`: vyska je
   ovlivnena sirenim, prekryvem a regularizaci.
2. **Lambda neni jedno cislo, ale rozhodovaci cesta.** Oba projekty pocitaji
   GCV i L-krivku, oba detekuji, ze optimum sedi na okraji rozsahu, a oba
   odmitaji okrajovou hodnotu prohlasit za optimum.
3. **Seriove L jako neregularizovana rusiva promenna, lambda vybrana na
   systemu vcetne L.** Jejich `dual_lambda.py` sklada matice s indukcnim
   sloupcem a penalizaci na nem nulovou ("nuisance columns have zero
   penalty"); nase `extension._solve_on_grid` vybira lambdu na systemu, ktery
   L obsahuje, se stejnym zduvodnenim (jinak neregularizovana promenna
   pohlti, co penalizace vytlaci z gamma).

A jedna vec, kde se shodujeme v tom, co **neni** mozne: jejich reference
explicitne rika, ze GCV je jen aproximace linearnim vyhlazovacem, ne presna
krizova validace nezaporneho odhadu. Nas `gcv.py` ma v komentari totez.

### Kde jsou vecne napred

**Vyber lambda: dve nezavisla kriteria misto jednoho retezeneho.** Nas hybrid
hleda L-krivku jen v okne +-1.5 dekady kolem minima GCV, takze L-krivka neni
nezavisly nazor -- vetsi neshodu nez 1.5 dekady nemuze principialne ukazat.
Jejich pravidlo "Never average the two lambdas merely to obtain a compromise"
jsme prevzali: drivejsi geometricky prumer pri rohu hluboko pod GCV je pryc.
K tomu hlasi **ploche minimum GCV** (body do 1 % minima pres vic nez dekadu),
ktere nas kod nerozpozna -- s jejich prahem by se ale sepnulo i na bezne
zasumenych datech (real_gamry 1.5 dekady, syntetika s 2 % sumu 1.0 dekady).

Poznamka k preferenci: oni davaji prednost GCV, my od 2026-09-29 vetsi z obou
lambd. Mereni ukazalo, ze na NNLS-DRT obe kriteria chybuji stejnym smerem
(lambda prilis mala), takze vetsi hodnota je mensi chyba. Nezavisla L-krivka
pres cely rozsah (jejich cesta) by nam nepomohla: na `real_gamry_example`
nasla falesny roh u 0.18.

**GCV je u nich vnitrne konzistentni.** Jejich skore bere reziduum i stopu z
tehoz linearniho vyhlazovace. Nase `compute_gcv_score` bere reziduum z NNLS
reseni a stopu z linearniho, tedy smes dvou ruznych odhadu. Neni to chyba --
obe varianty jsou aproximace -- ale hodnota nasi lambdy neni srovnatelna s
"ucebnicovym" GCV a nas docstring to rika jen castecne.

**Rozhodnuti o seriovem L souteznim modelu.** Tohle drive byla nase nejvetsi
mezera a ve v0.41.0 je zavrena. Zbyva rozdil v rozhodovacim pravidle: nase
`auto` zapne L podle znamenka Im(Z) v horni dekade (levne, ale jen spoustec),
oni resi oba modely a porovnaji rekonstrukci plus regresi v horni dekade.
Nase pravidlo neodhali indukcni prispevek, ktery Im(Z) jeste neprevrati do
kladnych hodnot. Obrana je levna: pri zavedeni ve v0.41.0 bylo zmereno, ze
vynucene L na datech bez indukcniho konce posune R_pol o < 0.1 %, takze
`--drt-inductance on` je pri podezreni bezpecna volba.

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

**Vazeni datoveho clenu.** Jejich skalovani `median(|Z|)` je jedna konstanta
pro cele spektrum: zlepsi podminenost, ale nezmeni, ktere body rozhoduji.
Nase vazeni `1/sqrt|Z|` (v0.38.0) je skutecne vazeni po frekvencich a na
syntetice rozlisi tri oblouky 50/200/2000 Ohm, ktere nevazene reseni slije do
dvou.

**Odezva za mericim oknem.** Puvodne jsme meli jen hlaseni podilu R_pol v
krajnim binu. Od v0.39.0 muze mrizka tau pokracovat za pomaly konec okna
(`--tau-extend auto` vybere nejmensi prodlouzeni, ktere navrseni odstrani) a
piky za oknem nesou `outside_window` a vzdalenost od okna. Jejich
`unsupported-outside-window` ma tedy u nas protejsek -- a navic
`R_pol_extrapolated_fraction`, ktery rika, kolik polarizacniho odporu lezi
mimo data. Oni RBF mimo okno maji, ale export orezou a mnozstvi mimo okno
nekvantifikuji.

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
rika taky. Tady neshoda neni. Ve v1.2 se na teto casti nic nezmenilo.

---

## Co ma jen jeden z nich

### Jen eis_analysis

- **Fitovani ekvivalentnich obvodu.** Cely modul `fitting/` -- parser vyrazu
  `R()-(R()|Q())`, analyticky jakobian, diferencialni evoluce, multistart,
  kovariance a konfidencni intervaly, AIC/BIC zebricek kandidatu,
  `auto_suggest` na navrh obvodu z dat, od v0.37.0 i Young-Gohruv prvek `YG`.
  Oni to **zamerne nedelaji**: SKILL.md explicitne zakazuje "physical
  equivalent-circuit fitting". Jejich nove texty o "podminene RC-DRT" se
  spolecnymi RL/RLC vetvemi se tomu blizi, ale zustavaji navodem bez
  implementace.
- **Voigt retezec** linearni regresi, s `mu` auto-M a variantami real/imag/complex.
- **Odhad R_inf** fitem `R-L-(R|Q)` s kontrolou urcitelnosti (`rinf_estimation/`).
- **Oxidova analyza** (`analysis/oxide.py`) -- tloustka a permitivita vrstvy,
  vcetne prvku DQ a YG.
- **Detekce odlehlych bodu** (`validation/outliers.py`).
- **OCV vizualizace**.

### Jen eis-drt-batch-analysis

- **Davkove zpracovani** adresaru s manifestem jako allowlistem.
- **Strukturovany vystup** do sedmi adresaru (`00_overview` az
  `06_reproducibility`) se stabilnim `spectrum_uid` na spojovani vysledku.
- **Otisk behu** -- hash zdroju, zamcene verze zavislosti, identita
  interpretu; `--resume` jede jen pri presne shode, `--export-only`
  pregeneruje vystupy bez prepocitavani.
- **Export obou kandidatu lambda** (GCV a L-krivka) vcetne rekonstrukci, aby
  uzivatel mohl volbu prezkoumat; pri doruceni lze ukazat dva nahledy vedle
  sebe se stejnymi osami.
- **Znamenkove GDRT** (kladna gamma = RC, zaporna = RL) pro indukcni smycky,
  ktere seriove L nepokryva -- nove jen na vyzadani.
- **Loewner RC/RL** a **difuzni DDT kandidati** (blocking `coth(s)/s`,
  transmissive `tanh(s)/s`, Gerischer) s explicitnimi prijimacimi prahy.
- **Model-and-reduce** -- odecteni rusiveho prvku s kontrolou, ze se piky
  po odecteni neposunuly (limit 0.35 dekady / 50 % plochy).
- **Kontrola experimentalni platnosti z metadat** -- opakovana mereni,
  amplitudova linearita, ustaleni po klidu. Rozlisuje "nepodarilo se overit"
  od "neproslo".
- **Screening driftu behem sweepu** -- oddeluje monotonni trend residui podle
  poradi sberu od nahodneho sumu.
- **Prezentacni grafy na vyzadani**: 3D vodopady, heatmapy pres napetove
  uzly, 2D vrstvene DRT, parove srovnani EIS.

---

## Davkove zpracovani a export: nejvetsi rozdil

Nas projekt **porad neexportuje zadna data**. Zadny `to_csv`, zadny
`json.dump`, zadny `savetxt` nikde v `eis_analysis/`. Vystup je matplotlib
okno a text v konzoli. Kdo chce cisla dal zpracovat, musi pouzit Python API a
napsat si export sam.

To je pro interaktivni praci s jednim spektrem v poradku a zamerne to sedi s
nasi CLI filozofii (`doc/CLI_OUTPUT_UNIFICATION.md`: knihovna vraci
`*Result` dataclassy, tiskne jen `cli/handlers/`). Ale znamena to, ze
srovnani deseti spekter mezi sebou je u nas rucni prace.

Jejich reseni je druhy extrem: grafy jsou az posledni krok, hlavni produkt jsou
CSV a JSON, a prezentacni grafy se generuji **jen na vyzadani** (`--trend-plots
off` je vychozi). Diagnosticke grafy automaticke jsou; selhani vykresleni
nikdy nesmi smazat ciselny vysledek -- JSON se zapisuje pred rendrovanim.

Ta posledni vec je dobry napad nezavisly na davce: **oddelit ciselny vysledek
od vykresleni**. U nas `calculate_drt` vola `_create_visualization` uvnitr
vypoctu (`drt/core.py`, krok 9) bez osetreni vyjimky, takze chyba v
`visualization/` shodi cely beh a spoctena DRT se ztrati.

---

## Nesrovnalosti a slaba mista v jejich releasu

- **Overeno jen na macOS arm64.** README to porad prizna: Windows a Linux maji
  vstupni body, ale 1.2.0rc1 na nich overen nebyl. Pri 477 MiB instalaci a
  vazbe na presne CPython 3.11 to neni maly zavazek.
- **Predbezne vydani s nejasnym stavem.** Korenovy changelog a GitHub release
  mluvi o zverejneni, `RELEASE_NOTES.md` v baliku o "locally prepared, not
  uploaded or publicly published". Release notes to vysvetluji (zmrazeny balik
  se po zverejneni nemenil), ale ctenar baliku dostane zastaraly stav.
- **Vazba na privatni API pyimpspec.** `dual_lambda.py` vola funkce s
  podtrzitkem (`_assemble_A_matrix`, `_solve_qp_cvxopt`, ...) a proto
  vynucuje presne pyimpspec 5.1.3. Kazda aktualizace pyimpspec znamena
  rucni revizi.
- **Cislo 390 testu neni srovnatelne s nasimi 612.** Jejich testy jsou z velke
  casti kontrakty exportu, privacy scan, chovani resume a nove kontroly
  grafu; vlastni numeriky testuji mene, protoze numerika je z pyimpspec.
- **Dokumentace je rozdvojena.** Cinsky README pro cloveka, anglicky SKILL.md
  pro model. Obsah se casto prekryva -- stejny problem s jedinym zdrojem
  pravdy, jaky nas CLAUDE.md zakazuje. Zkraceni SKILL.md ho neodstranilo,
  jen presunulo text do dvanacti referenci.
- **SKILL.md je z velke casti prompt, ne specifikace.** Nova sekce "Clarify
  only decisions that change the work" je cela o tom, kdy se ma agent ptat
  uzivatele. Kdo pouziva jen Python CLI, cte ji zbytecne.
- **Verze v textech porad nesedi.** SKILL.md uvadi "Release candidate
  1.2.0rc1", ale dal pise "0.2.1 preserves ..." a "0.2.2 validates ..." jako
  aktualni pravidla; `references/eis-only.md` ma v nadpisu "(0.2.2)".

---

## Co si odtud vzit

Serazeno podle pomeru uzitku k praci:

1. **Varovat, kdyz se sondy lambda orezou.** ~~Hotovo ve v0.36.0.~~
2. **Oznacit piky u okraje mericiho okna.** ~~Hotovo ve v0.36.0~~, doplneno
   prodlouzenim mrizky tau a `outside_window` ve v0.39.0.
3. **Nahradit tichy geometricky prumer v hybridnim vyberu lambda varovanim.**
   ~~Hotovo~~: plati vetsi z obou lambd, roh L-krivky > 1 dekadu pod GCV jde
   do `warnings`. Mereni pritom ukazalo, ze geometricky prumer se na 21
   spektrech nespustil ani jednou; skutecna slabina je jinde a zustava
   otevrena: pri sumu 0.1-0.5 % je auto-lambda o 1-3 dekady nizko a DRT ma
   postranni laloky (~10 Ohm, 4-5 piku misto 2), ktere sondy stability jen
   oznaci jako `marginal`. Ploche minimum GCV s prahem 1 % by se sepnulo i na
   bezne zasumenych datech, nevyplati se.
4. **Oddelit vykresleni od vypoctu DRT.** Vyjimka v grafech nema zahodit
   spoctenou analyzu. Stale otevrene.
5. **Seriove L jako soubezny model v DRT.** ~~Hotovo ve v0.41.0~~
   (`--drt-inductance`), s jinym rozhodovacim pravidlem nez jejich, viz vyse.
6. **Export do CSV/JSON.** Ne cely jejich sedmiadresarovy strom -- ale jeden
   prepinac `--export-csv`, ktery zapise tau, gamma, piky a residua, by
   odstranil nejvetsi prakticke omezeni naseho CLI. Stale otevrene.

Nebrat si: davkovou architekturu, otisky behu a resume (resi problem, ktery
nemame), Codex Skill vrstvu, ani zavislost na pyimpspec -- ta by nam vymenila
MIT za GPL a pridala 400 MiB instalace za funkce, ktere z 80 % uz mame vlastni.
Ani jejich preferenci GCV pred L-krivkou: na NNLS-DRT GCV lambdu podcenuje
jeste vic nez L-krivka (zmereno, viz bod 3).

---

## Zaver

Nejsou to konkurenti, jsou to dve poloviny jineho problemu. eis_analysis je
**hlubsi na jednom spektru** -- fitovani obvodu, Voigt retezec, oxidova
analyza, propracovanejsi Lin-KK a nove i vazena DRT se seriovym L a
prodlouzenou mrizkou tau. eis-drt-batch-analysis je **siroky pres davku** --
reprodukovatelnost, export, znamenkove a difuzni vetve a oddeleni "spocitano"
od "overeno".

Verze 1.2 posunula jejich projekt smerem, ktery je nam blizsi: obycejna
nezaporna DRT se seriovym L je znovu hlavni cesta a exoticke vetve jsou na
vyzadani. Nejcennejsi zustava `references/methodology.md` a nove
`references/lambda-and-peak-interpretation.md` -- poctivy vycet toho, co DRT
a vyber lambda **netvrdi**. Tri z puvodnich peti doporuceni jsou od minula
hotova; nove pribylo jedno (bod 3), prime z cteni jejich `dual_lambda.py`,
a je take hotove.
