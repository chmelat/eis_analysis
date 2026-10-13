# Cold agents: vysledky

Plan a metodika: `doc/COLD_AGENT_PLAN.md`. Piskoviste `~/coldbox/`, pravda,
generator a doslovne debriefy v `~/coldbox-private/` (mimo repozitar).

## Kolo 1 (2026-10-13)

Commit piskoviste: ad2b9dd (eis_analysis 0.59.1). Model: Sonnet 5.5, kazda
uloha jednim agentem. Cena: 1-1.5 min a 0.21-0.36 USD na agenta, kolo
vcetne pilotu a zamitnutych behu asi 1.5 USD.

### Ulohy a vysledek proti pravde

Synteticka data (numpy, uzavrene vzorce), CSV bez komentaru, sum
Gaussovsky sigma = p*|Z|. Tolerance 5 %, tau 0.15 dekady.

| Id | Data (pravda) | Agent | Cesta | Proti pravde |
|---|---|---|---|---|
| T1 | R_s - (R_f \| C_f), sum 1 %; eps_r 25, A 0.5 cm^2: d = 50.0 nm | 50.1 +- 0.3 nm (stat. 0.09 nm + rozptyl vah) | CLI `--analyze-oxide`, pak API | ok |
| T2 | R_s 15 - (120 \| 1 uF) - (300 \| 300 uF), sum 0.5 %: tau 1.2e-4 a 0.09 s | 2 procesy; tau 1.20e-4 / 8.95e-2 s; R 120.0 / 299.4; R_s 15.01 | CLI: DRT, 4 kandidati, BIC | ok |
| T3 | mereni_2 = T2; mereni_1 = T2 s Re(Z) +15 %/dekadu pod 1 Hz | mereni_2 cele ok; mereni_1 ok nad ~1 Hz, pod drift | CLI: KK/Z-HIT, porovnani, orez | ok, vcetne hranice |
| T4 (pilot) | R_s 8 - (C_dl 2e-5 \| (R_ct 60 - Ws(150, 2 s))), sum 1 % | 8.00 / 2.01e-5 / 59.8 / 149.5 / 1.97; Ws | CLI | ok |

Izolace: bod 0 debriefu u vsech "knihovnu neznal"; prectene soubory jen
z `~/coldbox/`.

### Nalezy

Overeno reprodukci nebo ctenim dokumentace (princip 6). Kategorie podle
planu: A tichy spatny vysledek, B nespravna dokumentace, C tren
v API/chybejici informace, D chyba agenta.

**Z vlastniho reseni uloh (krok 1, pred agenty):**

| # | Kat. | Nalez | Narazil agent? |
|---|---|---|---|
| F1 | B | `load_data('x.csv')` -> "No ZCURVE section found"; PYTHON_API ("Both loaders...", "Explicit CSV loading") naznacuje, ze load_data CSV zvladne | ne (T1 sel rovnou na load_csv_data) |
| F2 | C | poradi argumentu `fit_equivalent_circuit(f, Z, circuit)` vs `fit_circuit_multistart(circuit, f, Z)` | ne |
| F3 | B | `calculate_drt` s vychozim `peak_method='scipy'`: `result.peaks is None`, PYTHON_API slibuje "List of dicts"; piky v `diagnostics.scipy_peaks` | ne (T2 pres CLI) |
| F4 | C | ciste 2xRC: "Auto-lambda at search-range edge (1.45e-08)" pri dokumentovanem rozsahu (1e-10, 1e-2) - hrana L-krivky; plus "sparse/spiky" | ano (T2): citoval, overil fitem, vysledku neublizilo |
| F5 | A/C | API: `zhit_validation(...).is_valid=True` pro mereni_1 (prumer 4.66 % < 5 %, lokalne az 25 %). CLI totez oznaci "marginal" s podezrelymi body | ne (T3 pres CLI, kde verdikt sedi) |
| F6 | C | `ZHITResult.quality` 0-1 proti `quality_threshold=5.0` "[%]" | ne |

**Od agentu:**

| # | Kat. | Nalez | Uloha | Overeni |
|---|---|---|---|---|
| G1 | B | `doc/OXIDE_ANALYSIS_GUIDE.md:34,469`: `freq, Z = load_data(...)` -> TypeError (loader vraci LoadResult) | T1 | grep |
| G2 | C | `OXIDE_ANALYSIS_GUIDE.md` (5x), `GCV_IMPLEMENTATION.md`, `GMM_PEAK_DETECTION.md` volaji `python3 eis.py`; README `eis` | T1 | grep |
| G3 | C | tloustka bez nejistoty: `OxideAnalysisResult` ji nema, pruvodce nerika, jak ji dostat (agent propagoval stderr C rucne) | T1 | kod |
| G4 | C | CLI "Suggested circuit: R(R_inf) - (R(R1) \| C(C1)) - ..." symbolicky, do `--circuit` nejde vlozit | T2 | reprodukce |
| G5 | C | tau clanku (R\|C) neni ve vypisu ani v JSON (kovariance v JSON je, agent ji prehledl: D) | T2 | JSON |
| G6 | C | urovne "Data quality: poor/marginal/acceptable/good/excellent" nemaji v README definici | T3 | grep |
| G7 | C | chybi navod, jak z validace urcit platny frekvencni rozsah a cemu verit pri neshode KK a Z-HIT (`--f-min` KK neomezuje, README:724) | T3 | README |
| G8 | D | `-q` skryje i vysledky KK; `--help` rika "hide INFO messages" | T3 | bez akce |

### Poznatky o metode

- **Agenti volili CLI** (T2-T4 jen CLI, T1 hlavne). Pasti Python API z kroku 1
  (F1, F2, F3, F6) proto zustaly netknute, F5 se v CLI neprojevi. Kolo 2
  potrebuje ulohy, ktere API vyzaduji (napr. "napis funkci, ktera ...").
- **Agenti ctou vyberove** (zacatek README, tabulky, grep), ne celou
  dokumentaci. Riziko v planu ("LLM cte celou a doslova") neplati.
- **Filtr Sonnetu** (`[reasoning_extraction]`) zamitl zadani, kdyz osnova
  zadala popis vlastnich uvah ("kde jsi tapal nebo hadal, co jsi cekal").
  Osnova preformulovana na zpravu o pouzitelnosti dokumentace (nejasna
  mista, chybove hlasky, nesoulad); obsah pro triage stejny, prosla 4x.
- **Stav Claude Code mezi behy:** `~/.claude/projects/<cwd>/` (prazdna
  pamet, odlozeny vystup nastroju) - `run_agent.sh` ho pred behem maze.
- Umele ulohy kola 1 byly pro agenta snadne (vse v toleranci). Hodnota kola
  je v nalezech G1-G7, ne ve vysledcich.

## Kolo 2 (2026-10-13): totez s Haiku

Stejna data, zadani a piskoviste (ad2b9dd) jako kolo 1, model Haiku 5.5.
Cena 0.03-0.05 USD na agenta (Sonnet 0.21-0.36), cas 2.5-4 min (Sonnet
1-1.5), 16-25 tahu (Sonnet 7-10).

| Id | Agent | Proti pravde |
|---|---|---|
| T1 | 50.1 +- 0.2 nm | ok |
| T2 | 2 procesy; tau 1.20e-4 / 8.95e-2 s; R 120.0 / 299.4; R_s 15.0 +- 0.1 | ok |
| T3 | mereni_2 cele ok; mereni_1 jen 1 Hz - 1 MHz | ok, vcetne hranice |
| T4 | 8.000 / 20.07 uF / 59.84 / 149.5 / 1.969 s; Ws proti Wo a Q podle BIC | ok |

Izolace: bod 0 u vsech "neznal"; jen `~/coldbox/`.

Haiku pouzival Python API vic nez Sonnet (T1-T3), proto narazil na pasti
z kroku 1: F1 (T1, `load_data` na CSV), F3 (T2, `drt.peaks` je None,
spravne nasel `scipy_peaks`), F5/F6 (T3: `is_valid` True pri `quality`
0.067 "by ctenare zmatlo"), F4 (T2: rozlisil hranu L-krivky od hrany
rozsahu podle GCV_IMPLEMENTATION.md), G1 (T1: rozbaleni LoadResult).

### Nove nalezy

| # | Kat. | Nalez | Uloha | Overeni |
|---|---|---|---|---|
| H1 | B | PYTHON_API Basic Workflow importuje `compute_information_criteria`, `analyze_residuals`, `ResidualDiagnostics` z `eis_analysis` -> ImportError (jsou jen v `eis_analysis.fitting`) | T3, T2 | import |
| H2 | B | PYTHON_API (validation): `kramers_kronig_validation(..., include_L=True, allow_negative=True)` -> TypeError; signatura ma jen `mu_threshold, max_M, auto_extend_decades, extend_decades_range, include_C` | T3 | inspect |
| H3 | B | PYTHON_API (drt, ~419): `from eis_analysis.drt import ... GMM_AVAILABLE` - v modulu neni | T3 | hasattr |
| H4 | C | CLI "Analysis quality: GOOD" (DRT) nedefinovano; u mereni_1 GOOD vedle "Data quality: poor" - dva podobne stitky, ruzny vyznam | T3 | grep README |
| H5 | C | jmena parametru ve vypisu (`R0`, `R_W0`, `τ_W0`) proti dokumentaci (`R_W`, `tau_W`) | T4 | - |
| H6 | C | `params_stderr`: nerika, zda 1 sigma a zda je kovariance skalovana rozptylem rezidui | T1 | - |

### Srovnani modelu

- Oba modely vyresily vsechny ulohy. Haiku za sedminu ceny, dvakrat
  pomaleji, s vic kontrolami (vahy, multistart, lambda-probe, hybrid
  spektrum v T3) a delsimi zpravami.
- Haiku nasel vic nalezu typu B (H1-H3), protoze pouzil Python API; Sonnet
  zustal u CLI. Rozdil je v ceste, ne v schopnosti: pasti API ceka kazdy,
  kdo pujde pres API.
- Pro dalsi kola staci Haiku (levny, dukladny); ulohy musi vest pres API,
  pokud ma kolo testovat API.

## Stav nalezu

| Nalez | Stav |
|---|---|
| G1, F1, F3, H1, H2, H3 (kategorie B) | opraveno v dokumentaci (PYTHON_API, OXIDE_ANALYSIS_GUIDE; navic klice piku GMM, scikit-learn v instalaci, `len(result.peaks)` v prikladech). Kontrola: vsech 55 importu z PYTHON_API a pruvodce oxidem projde |
| F5 (A/C) | dokumentovano: `ZHITResult.is_valid` zustava prumer, docstring a PYTHON_API rikaji, ze lokalni drift zredi (studie nize) |
| F6 (C) | opraveno: pole `ZHITResult.quality` (0-1) odstraneno; `quality_threshold` je prah prumeru pro `is_valid` |
| G2, G3, G6, H5, H6 (C) | opraveno v dokumentaci: `eis` misto `python3 eis.py` (pruvodce oxid, GCV, GMM); nejistota tloustky z fitu (OXIDE_ANALYSIS_GUIDE, 4); tabulka stitku "Data quality" a rozdil verdiktu KK/Z-HIT (README, Data quality); jmena parametru, stderr a CI (README, Circuit fitting) |
| G4, H4, F4 (C) | opraveno v CLI: navrzeny obvod s hodnotami (vlozitelny do `--circuit`, s R_inf a seriovym L z DRT); "Analysis quality" -> "Suggestion quality"; hrana lambda rozdelena na dva nezavisle priznaky a dve varovani (mez rozsahu hledani, okraj okna L-krivky). "sparse/spiky" na cistem 2xRC zustava (kalibrace prahu je samostatna uloha) |
| F2, G5, G7 (C) | otevreno |

### F5: studie pravidla platnosti Z-HIT

30 spekter z kazde rodiny stress testu (KK-konzistentni, sum 0-3 %), ciste
a s driftem (|Z| nebo Re(Z) x (1 + 5-15 %/dekadu) pod nahodnou f_x), plus
17 souboru z `example/`. Podil spekter uznanych za validni:

| Pravidlo | ciste, sum 0 / 0.1 / 1 % | drift \|Z\|, sum 0.1 % | example/ |
|---|---|---|---|
| prumer \|res_mag\| < 5 % (stav) | 33/33, 49/49, 42/44 | 46/49 | real_gamry ne, ostatni ano |
| <= 5 % bodu nad 5 % (jako KK) | 33/33, 46/49, 29/44 | 2/49 | odmitne 8 cistych s 1 % sumem |
| klouzavy median 1 dekada < 5 % | 33/33, 48/49, 39/44 | 2/49 | real_gamry ne, ostatni ano |
| Lin-KK `is_valid` | 26/33, 31/49, 21/44 | 15/49 | - |

Z-HIT zesiluje sum (pri 1 % sumu body az 8 %), proto pravidlo na body
odmita cista data; klouzavy median sum potlaci a drift ponecha. Prakticky
dopad je maly: CLI drift ukaze stitkem "marginal (check for drift)"
a hlasenim podezrelych bodu, na `is_valid` nic v knihovne nestavi, verdikt
na trech merenych DTA se nemeni. Zvoleno: jen dokumentace; vyhlazene
pravidlo az pro mereni s driftem, ktery Lin-KK propusti (kandidat: ZrO2,
`EISPOT-M136113-4.DTA` ma vyhlazene maximum 4.76 %).
