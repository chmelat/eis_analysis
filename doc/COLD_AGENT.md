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

## Kolo 3 (2026-10-13): Python API, Haiku

Piskoviste z f1b1f33 (po opravach B, F5/F6, G2-G6, H4-H6, F4). Ulohy
vyzaduji funkci v `reseni.py` nad Python API (ne CLI). Cena 0.05-0.09 USD,
3.5-4.5 min, 21-40 tahu na agenta. Kazdy `reseni.py` jsem spustil sam.

| Id | Data (pravda) | Vysledek `reseni.py` | Proti pravde |
|---|---|---|---|
| T5 | 3 vrstvy, eps_r 18, 0.25 cm^2: 30 / 50 / 80 nm | 29.95 +- 0.06 / 50.10 +- 0.08 / 79.96 +- 0.13 nm | ok |
| T6 | 2 a 3 RC (tau 1e-5 - 0.2 s), R_s 10 | 2 a 3 procesy, vsechna tau a R do 1 %, R_s 10.01 / 9.94 | ok |
| T7 | X ciste; Y drift \|Z\| +10 %/dek pod 0.5 Hz | X cele; Y od 0.63 Hz | ok, hranice 0.1 dekady |
| T8 | Randles s Wo (R_W 80, tau 5 s) | Wo, dBIC 985, parametry do 1 % | ok |

Opravene pasti uz nebrzdi: T5 zkusil `load_data` na CSV a rozbaleni
LoadResult, oboje oznacil za "ocekavane" podle dokumentace; T6 sel rovnou na
`diagnostics.scipy_peaks`; T8 importoval `compute_information_criteria`
z `eis_analysis.fitting` se spravnou trojici; T7 rozhodoval podle Lin-KK
a reziduí Z-HIT, ne podle Z-HIT `is_valid`.

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| K7 | A/C | `kramers_kronig_validation` na spektru s kapacitnim NF koncem (reflexni Wo) bez `include_C`: `is_valid` False, 80/81 bodu nad 5 %, `noise_estimate` 34.7 %, `warnings` prazdne. Napoveda `--kk-series-c` zije jen v CLI handleru; v API zadna. S `include_C=True` sum 1.23 %. Sum 34.7 % jde dal do `compare_fit_to_noise` | T8 |
| K6 | C | `find_outliers` na T7/Y: Lin-KK ma 6 bodu nad 5 % (verdikt poor), ale hlaseni podezrelych bodu z nich neukaze zadny - pojistka 4x median (median 2 %, prah 8 %). README rika, ze na merenych spektrech pojistka "nic nemeni"; tady hlaseni a verdikt nesouhlasi | T7 |
| K1 | B | PYTHON_API (fitting, ~835): `R(1000, fixed=True)` -> TypeError; plati `R("1000")` | T5 |
| K8 | B | PYTHON_API Example 6: `Z_fit = circuit.impedance(frequencies, params)` - `circuit`, `params` nedefinovane (`chain.circuit`, `chain.initial_params`) | T8 |
| K2 | C | `params_stderr` pevneho parametru je 0.0; dokumentace zmini `inf` (neidentifikovatelny) a `null` (export), pevny ne | T5 |
| K3 | C | PYTHON_API (rinf): "R_inf_range ... (misses: docstring)" - nesrozumitelne | T6 |
| K4 | C | vysledky `lambda_probe` (`diagnostics.stability.peak_stability`) v PYTHON_API nepopsane | T6 |
| K5 | C | `find_outliers`, `OutlierReport`, `OutlierPoint` chybi v PYTHON_API | T7 |
| K9 | C | Ws a Wo: nazvy ("finite-length" / "finite-space") nerikaji hranici, navic do v0.52.0 byla tanh forma `Wo` (podruhe, po kole 2) | T8 |

## Kolo 4 (2026-10-14): T1-T8 znovu, nove parametry, Haiku

Piskoviste z cefa1f8 (po G5, F2, G7, K7 a opakovani Lin-KK se seriovym C),
stare behy smazany. Zadani beze zmeny krome cisel (T1, T5); vsechny
parametry a seedy nove (`~/coldbox-private/gen_r4.py`), protoze README nese
T3 z kola 1 jako priklad. T3 a T7 prohodily spatny soubor, T7 driftuje ve
fazi (5 deg/dek pod 2 Hz) misto v \|Z\| - spousti `series_c_passes=False`.
Cena 0.41 USD celkem (0.03-0.10), 2-5 min, 16-25 tahu. Bod 0 vsude "ne",
zadna zminka o vyvoji. `reseni.py` T5-T8 jsem spustil sam.

| Id | Data (pravda) | Vysledek | Proti pravde |
|---|---|---|---|
| T1 | eps_r 9.5, 0.8 cm^2: d = 35.0 nm | 35.0 +- 0.06 nm (stat.), +- 0.5 nm metodicka | ok |
| T2 | R_s 22 - (60 \| 0.4 uF) - (450 \| 150 uF) | 2 procesy, tau 24.0 us / 67.6 ms, R 60.1 / 450.1, R_s 22.00 | ok |
| T3 | m1 ciste; m2 Re +25 %/dek pod 0.3 Hz | m1 cele; m2 od 1.6 Hz (beh Lin-KK), ze srovnani s m1 0.4 Hz | ok, hranice o 0.7 dek. vys (konzervativne) |
| T4 | Randles Ws (R_W 220, tau 0.8 s), jen CLI | Ws (Wo dAIC 1007), 5 parametru do 0.5 % | ok |
| T5 | 3 vrstvy, eps_r 12, 0.4 cm^2: 25 / 60 / 110 nm | 24.96 / 59.90 / 109.92 nm | ok |
| T6 | 2 a 3 RC (tau 3e-6 - 6e-2 s), R_s 18 | 2 a 3 procesy, tau i R do 1.5 %, R_s 17.96; tau_stderr z `arc_time_constants` | ok |
| T7 | Y ciste; X faze -5 deg/dek pod 2 Hz | Y cele; X od 0.4 Hz (tabulka), z grafu odhad 1-2 Hz | ok, tabulkova hranice o 0.7 dek. niz |
| T8 | Randles Wo (R_W 120, tau 3 s) | Wo, dBIC 1028, 5 parametru do 0.5 % | ok |

Nove funkce se osvedcily: T7 precetl `series_c_passes=False` a varovani
"drift or another violation" spravne a `find_outliers` mu sam dal rezidua
behu s C; T8 pouzil `include_C=True` na zaklade varovani; T6 vzal nejistoty
tau z `arc_time_constants`; T3 a T7 postupovali podle "Finding the valid
frequency range".

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| L1 | B | README "Finding the valid frequency range", krok 3: "the run in the Lin-KK residuals reaches the true boundary" neplati obecne. T3: souvisly beh Lin-KK konci na 0.4 Hz, za mezerou body nad 5 % az do 1.26 Hz (agent vzal 1.6 Hz), Z-HIT 0.16 Hz, pravda 0.3 Hz. T7 (faze 5 deg/dek): Z-HIT konci na 0.13 Hz, Lin-KK s C nema beh, pravda 2 Hz. Konec behu je dolni mez rezu, ne zacatek driftu (0.1 dek. pod az 1.2 dek. nad vyssim koncem) | T3, T7 |
| L2 | A/C | PYTHON_API `compare_fit_to_noise`: priklad predava `kk.noise_estimate`; na kapacitnim konci je to sum Voigtova retezce (37.7 %), ne dat. T8 dostal pomer 0.03 ("vypadalo dobre") a zachranilo ho jen varovani KK. Priklad ma brat `kk.with_series_c`, kdyz `kk.capacitive_end` | T8 |
| L3 | B | PYTHON_API (~960, Voigt elements): `elem['id']`, `elem['tau']` -> `TypeError: 'VoigtElement' object is not subscriptable`; plati atributy (`elem.tau`) | T6 |
| L4 | B | README (Voigt chain, ~608 a `--voigt-chain`, ~853): "then nonlinear refinement"; CLI vypise "Using linear result (no nonlinear optimization)" a zadny nelinearni krok neudela | T2 |
| L5 | B | PYTHON_API `arc.stderr_note` uvadi jen 'no covariance' a 'n0 at its upper bound'; kod vraci i 'covariance not usable' | T2 |
| L6 | C | README (time constants, ~509): "every (R\|C) and (R\|Q) arc" - u Randlese `(R - Ws \| C)` je `time_constants` prazdne (zamerne, C vedle delsi vetve), README to nerika | T4 |
| L7 | C | OXIDE_ANALYSIS_GUIDE 4 (nejistota): `param_labels.index('C0')` jen pro jediny C; `OxideAnalysisResult` neuvadi label dominantniho prvku, ve vicecelnem obvodu ho uzivatel nenajde | T1 |
| L8 | C | Tloustka (oxid) neni v `--save` JSON/CSV, jen v konzoli | T1 |
| L9 | C | PYTHON_API (DE): `seed` a `archive_check` chybi v prikladu i popisu, i kdyz jsou v podpisu (multistart `rng` je) | T8 |
| L10 | C | PYTHON_API (`find_optimal_lambda_gcv`, ~469): "searches DRT_LAMBDA_RANGE = (1e-10, 1e-2)", ale `calculate_drt` na T6/B vrati GCV i vybranou lambda 5.97e-11, pod spodni mezi (varovani "at search-range edge" je) | T6 |
| L11 | C | Stitky parametru jsou pozice v obvodu; dva clanky (R\|C) si pri fitu muzou vymenit hodnoty (T3: R1 = 452 misto 60) - porovnavat podle tau. Nezmineno | T3 |

Nereprodukovano: T2 hlasil degenerovany lokalni 3RC fit (cond 3e31,
stderr inf) bez varovani; vlastni 3RC fit na T2 da cond 1e30 a varovani
"Rank-deficient Jacobian ... Non-identifiable parameter(s)". Drobnosti
bez akce (D nebo opakovane zname): sparse/spiky DRT na cistych datech (F4),
tichy default plochy 1 cm^2 u CSV (dokumentovano, vypis "assuming"),
`initial_guess` v poradi `param_labels`, `warnings` vs `all_warnings`.

## Kolo 5, krok 1 (2026-10-14): slozitejsi obvody, Haiku

Piskoviste z 0b54101 (v0.60.0), behy kola 4 smazany. Ulohy maji kratkou
cestu, ktera vede ke spatne odpovedi (pocet piku DRT, R_inf z VF konce,
R_pol z DRT), a spravnou, ke ktere vede varovani nebo fit
(`~/coldbox-private/gen_r5.py`). Zadani T9 je zneni T2, T10 a T12 jen
jina otazka. Cena 0.18 USD celkem, 4-6 min, 21-27 tahu. Bod 0 vsude "ne".

| Id | Data (pravda) | Kratka cesta | Vysledek agenta | Proti pravde |
|---|---|---|---|---|
| T9 | R_s 15 - ZARC(200, 1 ms, 0.85) - ZARC(150, 6.3 ms, 0.80), 0.8 dek. od sebe | DRT: 1 pik, bez varovani | 2 procesy (BIC: 1 oblouk +474), tau 0.99 / 5.9 ms, R 193 / 158, R_s 14.96; sam: korelace R1-R2 -0.999 | ok |
| T10 | L 1 uH - R_s 5 - (30, tau 0.3 us) - (300, 0.3 s) | Re(Z) na 1 MHz 11.6, HF median 17.5 | R_s 5.0 +- 0.4 (fit s L 5.20 +- 0.16), median i horni mez odmitl | ok |
| T12 | R_s 12 - (150, 1 ms) - (800, 5 s), okno 0.1 Hz - 100 kHz | DRT R_pol 425 (s varovanim) | 2 procesy, R_pol 966 +- 14, R2 jako extrapolace za okno | ok |

Vsichni tri vzali kratkou cestu, videli, ze nesedi, a dosli fitem
a porovnanim kandidatu ke spravne odpovedi.

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| M1 | A | Vychozi R_inf (median 5 VF bodu) muze lezet nad Re(Z) na f_max: T10 17.48 Ohm proti 11.62 na 1 MHz (pravda 5). Bez varovani; DRT bezi s nim, "Suggestion quality: GOOD", navrzeny obvod `R(17.48) - ...`. U pasivniho systemu je Re(Z(f_max)) horni mez R_inf (tak ji bere `--ri-fit`) | T10 |
| M2 | C | DRT slouci dva procesy 0.6-1 dekadu od sebe do jednoho piku bez varovani (T9: 1 pik, BIC chce 2 s dBIC 474); `MODEL_SELECTION_AIC_BIC.md` radi "let DRT tell you how many processes there are". Jen `--lambda-probe` rekne "not stable" | T9 |
| M3 | C | Navrzeny obvod se vypise i pri vlastnim "Inconsistency": T12 `R(425.5)` s kvalitou ACCEPTABLE vedle "sum(R_i) = 425.5 vs R_pol = 180.2" (pravda 950); T10 dedi horni mez R_inf 11.62 jako odpor | T12, T10 |
| M4 | C | `--save` do neexistujiciho adresare: `!! Error saving figure` u kazdeho grafu, nic se neulozi, exit 0 | T12 |
| M5 | C | `--drt-inductance auto` prida L jen pri Im(Z) > 0 v horni dekade; T10 ma L 1 uH zakryte rychlym obloukem (fit s L: dBIC 459) - DRT L neuvede a README to nerozlisuje od "zadna indukcnost" | T10 |
| M6 | C | `--voigt-chain --voigt-auto-M` vypise "Allow negative R_i: True" bez `--voigt-allow-negative`; README rika, ze bez nej se zaporne R_i odstrani (iterace mu je pripousti, vysledek ne) | T9 |
| M7 | C | `RINF_ESTIMATION.md` odkazuje na `archive/AUDIT_ri_fit_2026-09-25.md` a `doc/STRESS_TEST.md` jako na cesty, ne odkazy; v dokumentaci, kam README vede, nejsou | T10 |
| M8 | C | "R_pol (from data)" neni v README ani PYTHON_API definovano (T12: 180 Ohm pri pravde 950) | T12 |

Drobnosti bez akce: knihovni varovani jmenuje `include_C` (CLI hned pod nim
pise `--kk-series-c`), `KK: M=23` neni pocet procesu, `parse_circuit_expression`
neni top-level (dokumentovano).

## Kolo 5, krok 2 (2026-10-14): T11, povlak s pory, Haiku

Piskoviste z d0db5e0 (v0.60.1, po opravach M1-M8), behy kroku 1 smazany.
Zadani popisuje fyziku (povlak s pory, koroze na dne poru s
semi-nekonecnou difuzi), ne obvod (`~/coldbox-private/gen_r5b.py`).
Cena 0.10 USD, 7.5 min, 30 tahu. Bod 0 "ne".

| Id | Data (pravda) | Vysledek agenta | Proti pravde |
|---|---|---|---|
| T11 | R_s 20 - (Q_c 2e-8, 0.92 \| (R_po 3000 - (Q_dl 5e-6, 0.80 \| (R_ct 2e4 - W sigma 2000)))), sum 1 % | seriovy model `R_s - (R_p \| Q_c) - ((R_ct - W) \| Q_dl)`: R_s 19.88, R_po 2940, Q_c 2.03e-8 / 0.920, Q_dl 5.08e-6 / 0.797, R_ct 20080, sigma 1991; C misto Q a model bez W zamitl (dBIC > 530), overil multistartem, vahami a f_min | ok, vsech 8 do 2 % |

Vnoreny model povlaku agent nezkusil. Pri oddelenych casovych konstantach
(26 us proti ~0.1 s) a R_po << R_ct se s nim seriovy prakticky kryje;
jinde by parametry nesedely.

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| M9 | B | README (~501): "With `--circuit` or `--de`" - prepinac `--de` neexistuje (je `--optimizer de`) | T11 |
| M10 | B | CLI blok DE vypise `Weighting: <vaha> (w=1/\|Z\|)` pro kazdou vahu (`cli/handlers/fitting.py:139` ma vzorec natvrdo), napr. `uniform (w=1/\|Z\|)`; jinde CLI pouziva spravnou tabulku (fitting.py:531, 585) | T11 |
| M11 | C | Dokumentace nema vnoreny model povlaku s pory (`R_s - (Q_c \| (R_po - (Q_dl \| (R_ct - W))))`) ani rozdil proti seriovemu; agent sestavil seriovy | T11 |
| M12 | C | Tabulka porovnani kandidatu zkracuje vyrazy obvodu (`R(20)-(R(2800)\|Q(1e-8, 0.9))-((R(30...`); `DE error` pred zpresnenim vedle `Fit error` bez vysvetleni | T11 |

Bez akce: `LoadResult` caveat (PYTHON_API hned pod popisem ukazuje
`result.warnings` - chyba agenta), helper pro C_eff obecneho Q (mimo
analyzu oxidu) - pozadavek, ne chyba.

## Kolo 6 (2026-10-15): T13, ZARC s ridkymi vadnymi body, Haiku

Piskoviste z d2ed640 (v0.60.3). Uloha od uzivatele: dva blizke oblouky,
maly R_s, vetsina bodu presna, asi 6 % bodu s vetsi chybou
(`~/coldbox-private/gen_r6.py`, seed 6112 vybran z 20 jako nejzretelnejsi).
Zadani je zneni T2 s "vse s nejistotami". Cena 0.09 USD, 9 min, 24 tahu.
Bod 0 "ne".

| Id | Data (pravda) | Kratka cesta | Vysledek agenta | Proti pravde |
|---|---|---|---|---|
| T13 | R_s 2 - ZARC(100, 1 ms, 0.85) - ZARC(80, 3.16 ms, 0.80), 0.5 dek. od sebe, sum 0.1 %; 5 z 81 bodu (20 kHz, 6.3 kHz, 31.6 Hz, 15.8 Hz, 0.03 Hz) navic s chybou sigma 2.5 % na Re i Im | KK "excellent", sum 0.70 %, podezrely jen bod 0.03 Hz ("often sample drift"); DE fit 2 ZARC: R 141.5 +- 7.7 / 37.9 +- 7.8, Fit vs. noise 1.0 | 2 procesy (BIC), R_s 2.004 +- 0.004, R 142 +- 8 / 38 +- 8, tau 1.18 / 5.7 ms; zavislost na f_min si vsiml (f >= 1 Hz: ~110 / 70) a rozdeleni R oznacil za nespolehlive, vadne body nehledal | pocet procesu a R_s ok; R1/R2 a tau2 spatne, 5 sigma mimo uvedenou nejistotu |

Spravna cesta: rezidua fitu nad 5x median oznaci presne pet vadnych bodu,
fit bez nich da R 104.9 +- 3.8 / 75.1 +- 3.8, tau 1.02 / 3.31 ms. Fit
startovany z pravdy skonci ve stejnem vychylenem minimu jako DE: chyba je
v datech, ne v hledani minima.

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| N1 | A | `KKResult.noise_estimate` (RMS pres pseudo chi^2) s peti vadnymi body 0.70 %, bez nich 0.10 %. Proto KK "excellent" a `Fit vs. noise` 1.0; se sumem 0.10 % by fit (0.45 %) byl 4.5x nad sumem a CLI by varovalo. Odhad z medianu rezidui KK je 0.20 % | T13 |
| N2 | A/C | `find_outliers`: rezidua Lin-KK radi pet vadnych bodu presne na prvni mista (2-5 %, 9-20x median 0.24 %), ale absolutni prah 5 % je propusti; oznacen jen bod 0.03 Hz pres Z-HIT, a hlaseni ho vysvetli driftem | T13 |
| N3 | C | Varovani "Residuals are not random" nabizi jen chybejici prvek nebo malo prvku, izolovane vadne body ne. Agent zkusil treti oblouk (BIC +10) | T13 |
| N4 | C | Dokumentace neukazuje, jak vyradit jednotlive body a fitovat znovu; CLI ma jen `--f-min/--f-max`, soft-L1 (`ROBUST_LOSS_SOFT_L1.md`) je neimplementovany navrh | T13 |
| N5 | C | CLI nema `--seed` pro DE (API `seed=` ano) | T13 |
| N6 | B | PYTHON_API "Fit result export": ukazka vola `parse_circuit_expression` bez importu (odkud, rika az text pod ni) | T13 |
| N7 | D | S `--no-kk` radek `Fit vs. noise` beze slova chybi; README (Are the residuals noise?) to rika ("With `--no-kk` the line is absent") - chyba agenta | T13 |

Bez akce: vztah `RinfResult` a `RinfEstimate`, R_inf z DRT (2.038, median
VF bodu na ocasu ZARC s n < 1) proti R_s z fitu (2.004).

## Kolo 6, krok 2 (2026-10-15): T14, totez s N1 + N2 + napovedou, Haiku

Overeni, zda N1, N2 a napoveda ve varovani `Residuals are not random`
(508c06a) staci, nebo je treba N3 (izolovane body z rezidui fitu a novy fit
bez nich). Jine parametry nez T13, ktere README ted cituje
(`~/coldbox-private/gen_r6b.py`, seed 6201 vybran ze 40 podle vychyleni
fitu). Piskoviste z 508c06a, poprve bez CHANGELOG (prani uzivatele: zapisuje
nalezy predchozich kol). Cena 0.05 USD, 4 min, 25 tahu. Bod 0 "ne".

| Id | Data (pravda) | Kratka cesta | Vysledek agenta | Proti pravde |
|---|---|---|---|---|
| T14 | R_s 3 - ZARC(120, 2 ms, 0.82) - ZARC(90, 7.96 ms, 0.78), 0.6 dek., 97 bodu, sum 0.15 %; 6 bodu (38 kHz, 14.7 kHz, 17.8 Hz, 0.26, 0.1, 0.046 Hz) s chybou sigma 2.5 % | KK: varovani N1 (sum 0.72 proti 0.29 % z medianu); kontrola po bodech 5 bodu, z nich 2 vadne (3 sousedy oznacil Z-HIT); runs test fit nezamitne (p 0.024), napoveda se nevypise; fit 150/60 | 2 procesy (BIC), R_s 3.00, R 149.9 +- 10.6 / 59.8 +- 10.7, tau2 11.8 ms; varovani N1 citoval, "Podezrele body z kontroly po bodech jsem nevyrazoval" | pocet procesu a R_s ok; R1/R2 a tau2 spatne jako u T13 |

N3 z vetve na tychz datech a agentove vlastnim prikazu (multistart, 2x R|Q):
5 ze 6 vadnych bodu, zadny dobry, novy fit 121.5 +- 5.5 / 88.5 +- 5.5.
Uzivatel N3 presto zamitl: ladene a overene na jednom modelu vadnych bodu
(stejnem jako T13 a T14), slozite. Prototyp soft-L1 dava 107/73 (T13)
a 128/82 (T14) a spatny model nezamaskuje; navrh, cisla a duvody
v `doc/ROBUST_LOSS_SOFT_L1.md`, 4.1-4.2.

### Nove nalezy (overeno ctenim dokumentace a vystupu)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| P1 | B | PYTHON_API Basic Workflow: po `fit_circuit_multistart` blok vysledku cte `result` z predchoziho jednoducheho fitu misto `ms_result.best_result` | T14 |
| P2 | C | "Suggestion quality: GOOD" u navrhu s jednim prvkem, ktery fit odmitne (fit/sum 2.4); dokumentace nerika, ze kvalita meri jen shodu souctu R s R_pol | T14 |
| P3 | C | `--ri-fit` vypise L = 2.5e-15 nH bez nejistoty; RINF_ESTIMATION nerika, jak cist L ~ 0 | T14 |
| P4 | C | Vyznamnost S = 6.1 u exponentu n; tabulka ve WEIGHTING_AND_STATISTICS 3.7 hodnoty nad 1 nevysvetluje | T14 |

## Kolo 6, krok 3 (2026-10-15): T13 a T14 s Opus

Tytez ulohy a data jako v krocich 1 a 2, model Opus misto Haiku. Piskoviste
z 508c06a (stejny kod jako HEAD, starsi `ROBUST_LOSS_SOFT_L1.md` bez
pribehu T13), bez CHANGELOG. Pro T14 navic z kopie README v piskovisti
vymazan popis pripadu T13 (tri pasaze s cisly: dva oblouky 0.5 dek., 5 z 81
bodu, 142/38 proti 100/80, sum 0.70 proti 0.10 %); obecny mechanismus
(varovani N1, kriterium 8x median, rho1 ~ 0 se zamitnutym runs testem jsou
izolovane body) zustal. Puvodni kopie `~/coldbox-private/r6/README.sandbox.orig.md`.

| Id | README | Cena, cas, tahy | Vysledek agenta | Proti pravde |
|---|---|---|---|---|
| T13 | s pripadem T13 | 1.16 USD, 6 min, 22 | vyradil 5 bodu z kontroly po bodech vlastnim filtrem CSV, Lin-KK pak sum 0.10 % a zadny dalsi bod; R 105 +- 4 / 75 +- 4, tau 1.021 / 3.31 ms, R_s 2.0008; citlivost na vahu (proportional 121/59) uvedl | ok, vse v toleranci |
| T14 | bez pripadu T13 | 1.00 USD, 6 min, 18 | vyradil jen 2 body oznacene Lin-KK (17.8 Hz, 14.7 kHz), 3 body jen ze Z-HIT vyhodnotil jako sousedy; R 123 +- 11 / 87 +- 11, tau 2.02 / 8.2 ms, R_s 3.000; systematicky rozptyl R1 123-150 podle vahy a vyrazeni uvedl | ok, vse v toleranci |

Bod 0 u obou "ne". Agent T13 si prikladu v README vsiml a napsal, ze ho
nepouzil; T14 bez prikladu dopadl stejne, takze k vyrazeni bodu vedl obecny
text (varovani N1, seznam N2, veta o izolovanych bodech), ne priklad. Oba
cetli README cele, T14 nic jineho. Haiku na tez dokumentaci (kroky 1-2)
varovani citoval a body nevyradil: dokumentace ted staci silnemu ctenari,
vysledek ale zavisi na rucnim zasahu, ktery CLI neumi (N4) - argument pro
soft-L1 (`doc/ROBUST_LOSS_SOFT_L1.md`).

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| O1 | C | Jeden vadny bod Z-HIT rozmaze do sousedu (T14: KK oznaci 17.8 Hz, Z-HIT 14.7 a 21.5 Hz; KK 14.7 kHz, Z-HIT 17.8 kHz); README u kontroly po bodech to nerika, agent to odvodil sam | T14 |
| O2 | C | Radek parametru pise `[at upper bound - CI not meaningful]` (`cli/handlers/fitting.py:345`) i pro hodnotu o dekadu pod mezi (Q2 = 1.37, mez 10; `classify_bound_status` meri 1 dekadu), varovani pod nim totez nazve "near upper bound" | T14 |
| O3 | C | `_fit.json`: pole `covariance` bez popisku a README nerika jeho poradi (je podle `parameters`); chybi navod na nejistotu odvozene veliciny (R1 + R2 pri korelaci -1), kterou si oba agenti spocitali sami | T13, T14 |
| O4 | C | Pojmenovani oblouku zavisi na vaze: s `modulus` je R1 rychly proces, se `sqrt`/`uniform` pomaly (L11; doklad pro `doc/BLOCK_ORDER_PLAN.md`) | T13, T14 |
| O5 | C | `proportional`: README "emphasizes high-frequency", PYTHON_API "emphasizes low-impedance region" (w = 1/\|Z\|^2 je to druhe; vysoke frekvence jen u typickeho spektra) | T13 |
| O6 | C | Titulek grafu DRT ukazuje lambda v plne presnosti (`visualization/drt.py:83, 86`, `lambda = 0.00011396899106368631`) | T13, T14 |

Opakovane: CLI nema vyrazeni jednotlivych bodu (N4, oba). Nereprodukovano
jako nalez: "95% CI u tau" (README rika CI u parametru, u casovych konstant
jen stderr - nejednoznacne cteni), zkracene obvody v tabulce porovnani (M12:
cele se vypisuji pod tabulkou - chyba agenta).

## Kolo 6, krok 4 (2026-10-15): T13 znovu s Haiku a `--exclude-flagged`

Puvodni T13 (seed 6112, stejna data jako v kroku 1), Haiku. Piskoviste
z 878b498: `--exclude-flagged`, kontrola po bodech ho doporucuje, README bez
prikladu T13/T14, bez CHANGELOG a `ROBUST_LOSS_SOFT_L1.md` (`setup.sh`,
`EXCLUDED`). Cena 0.11 USD, 7.4 min, 29 tahu. Bod 0 "ne".

| | Haiku, krok 1 | Haiku, krok 4 | Opus, krok 3 | Pravda |
|---|---|---|---|---|
| R1 / R2 [Ohm] | 142 +- 8 / 38 +- 8 | 104.9 +- 3.8 / 75.1 +- 3.8 | 105 / 75 | 100 / 80 |
| tau1 / tau2 [ms] | 1.18 / 5.7 | 1.02 / 3.31 | 1.02 / 3.31 | 1.00 / 3.16 |
| R_s [Ohm] | 2.004 | 2.001 | 2.001 | 2.0 |

Postup: fit na vsech bodech dal 142/38 jako v kroku 1; pak `--exclude-flagged
kk` (z doporuceni kontroly po bodech) jako hlavni vysledek, BIC kandidatu na
ponechanych bodech (2 CPE proti 1 o 764), stabilita pres `--f-min 0.1`,
v API `frequencies_flagged_by`. Rozdeleni R1/R2 oznacil za zavisle na
vyrazenych bodech (systematika +-35 Ohm), bod 0.03 Hz v nejnizsi dekade
overil oriznutim a doporucil opakovane mereni. Varovani bez nastroje (kroky
1-2) Haiku necetl jako vyzvu k cinu; jmenovany prepinac ano, za desetinu
ceny Opus.

### Nove nalezy

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| Q1 | C | README u `--exclude-flagged` nerikal, ze `kk`/`both` nejde s `--no-kk` (agent narazil na chybu parseru) | T13 |

Bez akce (chyba agenta): `drt.peaks` je pri `peak_method='scipy'` None -
priklad v PYTHON_API pouziva `gmm` a pro scipy `diagnostics.scipy_peaks`
(F3); DE 1.64 % -> 0.45 % po doladeni bez varovani - radek `Improvement (SSR)`
to vypisuje. Opakovane: "Suggestion quality: GOOD" u navrhu, ktery fit
odmitne (P2).

## Kolo 6, krok 5 (2026-10-16): T14 znovu s Haiku a `--exclude-flagged`

Puvodni T14 (seed 6201, stejna data jako v kroku 2), Haiku. Piskoviste
z 58b510f: `--exclude-flagged`, `--seed`, opravy P2-P4, O3, O4; README bez
prikladu T13/T14, bez CHANGELOG a `ROBUST_LOSS_SOFT_L1.md`. Cena 0.10 USD,
4.5 min, 26 tahu. Bod 0 "ne".

| | Haiku, krok 2 | Haiku, krok 5 | Opus, krok 3 | Pravda |
|---|---|---|---|---|
| R1 / R2 [Ohm] | 149.9 +- 10.6 / 59.8 +- 10.7 | 123.2 +- 10.6 / 86.6 +- 10.7 | 123 / 87 | 120 / 90 |
| tau1 / tau2 [ms] | 2.26 / 11.8 | 2.02 / 8.2 | 2.02 / 8.2 | 2.00 / 7.96 |
| R_s [Ohm] | 3.00 | 2.999 | 3.000 | 3.0 |

Postup: vychozi beh (varovani N1 a napoveda `--exclude-flagged kk`),
porovnani 4 kandidatu (BIC: 2 CPE, jeden oblouk o 327 hur), treti oblouk
(S = 0.005, BIC +14.6), pak `--exclude-flagged kk` jako hlavni vysledek
a beh na vsech bodech jako citlivost; `--lambda-probe`. Nejistotu R1 + R2
(0.2 Ohm pri korelaci -1) spocital z `covariance` v JSON (O3), `--seed`
pouzil a overil bitovou shodu dvou behu. Vysledek v toleranci, stejny jako
Opus. Chyby agenta: posun rozdeleni udava "60-70 Ohm" (je 27); vychylky
~2 % v 0.04-0.3 Hz videl (tri vadne body, ktere kontrola po bodech
nevypise), nevyradil je.

### Nove nalezy (overeno)

| # | Kat. | Nalez | Uloha |
|---|---|---|---|
| R1 | C | Varovani `Residuals are not random` vypise trend (p 0.84 / 0.34) a strukturu (power 0.06), oba nevyznamne, a nerekne, ze ho spustil runs test (Re p = 1e-3); agent to cetl jako rozpor (beh se 3 oblouky) | T14 |
| R2 | C | CLI vypise `Suggestion quality: GOOD` a `Suggested circuit: R(3.032) - (R(209.6) \| C(1.403e-05))` bez poznamky; fit ho zamitne (BIC o 327). README to od P2 rika, vystup ne | T14 |
| R3 | C | README u `--exclude-flagged` rika, ze par bodu posune rozdeleni R o desitky procent, ale ne, kterou variantu reportovat | T14 |
| R4 | C | README u `--lambda-probe` nerika, ze stabilni pik muze byt slouceni dvou procesu (STABLE neurcuje pocet procesu) | T14 |

Bez akce: DE 2.0 % -> 0.45 % po doladeni bez varovani (opakovane z kroku
4; radek `Improvement (SSR)` to vypisuje).

## Stav nalezu

| Nalez | Stav |
|---|---|
| G1, F1, F3, H1, H2, H3 (kategorie B) | opraveno v dokumentaci (PYTHON_API, OXIDE_ANALYSIS_GUIDE; navic klice piku GMM, scikit-learn v instalaci, `len(result.peaks)` v prikladech). Kontrola: vsech 55 importu z PYTHON_API a pruvodce oxidem projde |
| F5 (A/C) | dokumentovano: `ZHITResult.is_valid` zustava prumer, docstring a PYTHON_API rikaji, ze lokalni drift zredi (studie nize) |
| F6 (C) | opraveno: pole `ZHITResult.quality` (0-1) odstraneno; `quality_threshold` je prah prumeru pro `is_valid` |
| G2, G3, G6, H5, H6 (C) | opraveno v dokumentaci: `eis` misto `python3 eis.py` (pruvodce oxid, GCV, GMM); nejistota tloustky z fitu (OXIDE_ANALYSIS_GUIDE, 4); tabulka stitku "Data quality" a rozdil verdiktu KK/Z-HIT (README, Data quality); jmena parametru, stderr a CI (README, Circuit fitting) |
| G4, H4, F4 (C) | opraveno v CLI: navrzeny obvod s hodnotami (vlozitelny do `--circuit`, s R_inf a seriovym L z DRT); "Analysis quality" -> "Suggestion quality"; hrana lambda rozdelena na dva nezavisle priznaky a dve varovani (mez rozsahu hledani, okraj okna L-krivky). "sparse/spiky" na cistem 2xRC zustava (kalibrace prahu je samostatna uloha) |
| G5 (C) | opraveno: `arc_time_constants` (fitting) - tau, stderr metodou delta z kovariance a f pro kazdy clanek (R\|C) a (R\|Q); CLI je vypise pod parametry, `--save` do JSON (`time_constants`). Stderr overen Monte Carlo (50 realizaci sumu: pomer 0.86-1.12) |
| F2 (C) | opraveno: `fit_circuit_multistart` a `fit_circuit_diffevo` berou `(frequencies, Z, circuit, ...)` jako `fit_equivalent_circuit`; bez zpetne kompatibility (rozhodnuti uzivatele 2026-10-14) |
| G7 (C) | opraveno: README, Data quality, "Finding the valid frequency range" - postup (beh od kraje, hranice z grafu rezidui Lin-KK, Z-HIT konci driv, `--f-min`) s prikladem T3 a tabulkou neshod KK/Z-HIT; tabulka podezrelych bodu v CLI razena podle frekvence. Tvrzeni overena na T3 a 7 syntetickych spektrech (drift Re, \|Z\|, faze; kapacitni konec; VF artefakt; spicka) |
| K7 (A/C) | opraveno: `lin_kk_native` pri selhanem testu na kapacitnim NF konci bez `include_C` nastavi `KKResult.capacitive_end` a vyhradu ve `warnings` (API i CLI); CLI pak fit se sumem Lin-KK neporovnava a napovi `--kk-series-c`. Test, ktery projde, priznak nenastavi (nanejvys par okrajovych bodu) |
| K1-K5, K8, K9 (kolo 3) | opraveno v dokumentaci: pevny parametr `R("1000")`, jeho stderr 0 a CI; `R_inf_range` bez "(misses: docstring)"; vysledky `lambda_probe` (`diagnostics.stability`); `find_outliers` a `OutlierReport` v PYTHON_API; priklad 6 bez nedefinovanych promennych; Ws/Wo s "transmissive"/"reflective" v tabulkach README a PYTHON_API |
| K6 (kolo 3) | uzavreno dokumentaci: hlaseni podezrelych bodu muze mit mene KK bodu nez verdikt (pojistka 4x median, zduvodnena na merenych spektrech v outliers.py); verdikt dava `is_valid`, seznam jen lokalizuje. Popsano v PYTHON_API, chovani nezmeneno (rozhodnuti uzivatele 2026-10-13) |
| L1-L5 (kolo 4) | opraveno v dokumentaci: README "Finding the valid frequency range" - konec behu je dolni mez rezu, ne zacatek driftu (namereno na 4 driftech: 0.1 dek. pod az 1.2 dek. nad vyssim koncem), rez na nebo nad vyssim koncem, opakovane mereni jako nejlepsi lokalizace; PYTHON_API `compare_fit_to_noise` bere sum z `kk.with_series_c` na kapacitnim konci (T8: pomer 1.05 misto 0.027); priklad Voigt elementu pres atributy; README `--voigt-chain` bez "nonlinear refinement"; `stderr_note` i s "covariance not usable" |
| L6-L11 (kolo 4) | opraveno: README - casove konstanty jen pro par R/G s C/Q (Randles bez radku, R1*C0 rucne), jmena parametru jsou pozice, oblouky srovnavat podle tau, lambda az o dekadu za DRT_LAMBDA_RANGE (i PYTHON_API); `OxideAnalysisResult.element_label` (parametr dominantniho prvku ve fitu, navod na nejistotu v OXIDE guide uz bez pevneho C0); `--save` s `--analyze-oxide` zapise `<prefix>_oxide.json` (`oxide_result_record`, `save_oxide_result`); priklad DE se `seed` a `archive_check` |
| M1 (kolo 5) | opraveno: vychozi R_inf (DRT i `calculate_rpol` Voigt navrhu, `hf_r_inf`) omezen shora Re(Z) v nejvyssim kapacitnim bode, kdyz median lezi o vic nez 5 % nad obema nejvyssimi kapacitnimi body (jeden zakmit na f_max ho nespusti), s varovanim; T10 17.48 -> 11.62 Ohm. Na plochem konci se nezapne (0/300 pri 0.1-3 % sumu); v `example/` u dvou souboru s otevrenym obloukem nahore |
| M2-M8 (kolo 5) | opraveno: M3 navrh s "Inconsistency" (sum R_i vs R_pol za faktorem 2) ma kvalitu poor (T12: drive ACCEPTABLE); M4 `--save` vytvori chybejici adresar, nevytvoritelny skonci chybou s exit 1; dokumentace: M2 pocet piku DRT je dolni mez (procesy do ~1 dekady splynou, T9), MODEL_SELECTION radi kandidata s prvkem navic; M5 L pod obloukem otevrenym na f_max auto nenajde; M6 hledani M v `--voigt-auto-M` pripousti zaporna R_i vzdy; M7 RINF_ESTIMATION odkazuje na audit a STRESS_TEST.md jako odkazy; M8 definice "R_pol (from data)" a ze je dolni mez pri oblouku otevrenem na f_min |
| M9-M12 (kolo 5) | opraveno: M9 README bez `--de`; M10 blok DE vypisuje vzorec zvolene vahy (sdilena tabulka `_WEIGHTING_LABELS`); M11 CIRCUIT_PARSER ukazuje vnoreny model povlaku a jeho rozdil proti seriovemu (T11: BIC pro vnoreny o 3.1); M12 zkracene vyrazy v tabulce porovnani se pod ni vypisou cele podle `-c` (cast "DE error" byla chyba agenta) |
| N1, N2 (kolo 6) | opraveno (Unreleased): `find_outliers` oznaci Lin-KK bod nad 8x medianem rezidui i pod `max_residual`, mimo ridkou mrizku `sparse_fit` M > N/2 (T13: presne 5 bodu; kalibrace na stress spektrech v README); `KKResult.noise_estimate_median` a varovani, kdyz `noise_estimate` je vic nez 2x vetsi (T13: 0.70 vs 0.20 %); poznamka CLI k NF bodum zminuje i nesoulad Lin-KK na konci spektra |
| N3 (kolo 6) | detekce izolovanych bodu z rezidui fitu a novy fit bez nich (N3) zkousena a zamitnuta (kolo 6, krok 2; vetev smazana); v `main` napoveda ve varovani `Residuals are not random`. Reseni: soft-L1 (`--robust-loss`) podle `doc/ROBUST_LOSS_SOFT_L1.md`, neimplementovano |
| P1 (kolo 6, krok 2) | opraveno v dokumentaci: PYTHON_API Basic Workflow bere po multistartu `result = ms_result.best_result` |
| P2-P4 (kolo 6, krok 2) | opraveno v dokumentaci: P2 README a PYTHON_API - kvalita navrhu hodnoti jen bilanci sum(R_i) vs R_pol a pocet varovani, ne fit; P3 RINF_ESTIMATION - L = 0 je dolni mez (zadna indukcnost v okne), opravena mez `L >= 0` (psalo 1 pH); P4 README a WEIGHTING_AND_STATISTICS 3.7 - S > 1 cist jen proti prahu 0.01, velikost neradi parametry |
| O1 (kolo 6, krok 3) | zdokumentovano (README, Per-point residual check: Z-HIT rozmaze vadny bod do sousedu, proto `--exclude-flagged` bere vychozi body od Lin-KK) |
| O5, O6 (kolo 6, krok 3) | opraveno: O5 README u `proportional` "low-\|Z\| points" (vysoke frekvence jen u typickeho spektra), PYTHON_API stejne; O6 titulek grafu DRT pise lambda na 3 platne cislice |
| O2 (kolo 6, krok 3) | opraveno: radek parametru pise `[near ... bound]` jako varovani pod nim (kriterium 1 dekada beze zmeny, hodnotu meze uvadi varovani) |
| O3, O4 (kolo 6, krok 3) | opraveno v dokumentaci: O3 PYTHON_API (Fit result export) - poradi a jednotky `covariance`, nejistota R1 + R2 z mimodiagonalniho clenu; O4 README - poradi oblouku se muze zmenit i s vahou a behem DE, porovnavat podle tau (kanonicke poradi dal `BLOCK_ORDER_PLAN.md`, neimplementovano) |
| Q1 (kolo 6, krok 4) | opraveno v README (popis `--exclude-flagged`) |
| N4 (kolo 6) | CLI: `--exclude-flagged {kk,both,any}` vyradi body z kontroly po bodech ze vseho za validaci (Unreleased); kontrola po bodech a napoveda u rezidui ho jmenuji. Robustni fit (soft-L1) odlozen 2026-10-15: prepinac to pokryva (`ROBUST_LOSS_SOFT_L1.md`, 4.3) |
| N6 (kolo 6) | opraveno v dokumentaci: ukazka "Fit result export" importuje `parse_circuit_expression` |
| N5 (kolo 6) | opraveno: CLI `--seed N` pro DE i multistart (N7 je chyba agenta) |
| R1 (kolo 6, krok 5) | opraveno: varovani `Residuals are not random` ma prvni radek `runs test: p = ... (fails below 0.01)` - spoustec; README a RESIDUAL_DIAGNOSTICS s nim. Vyklad "runs propadne, trend a struktura nevyznamne = izolovane body" v CLI ne (chybi prah nevyznamneho trendu) |
| R2 (kolo 6, krok 5) | opraveno: pod navrzenym obvodem poznamka - kvalita meri jen sum(R_i) vs R_pol, blizke piky splynou, prvky jsou idealni C; fitovat navrh, navrh s prvkem navic a s Q, porovnat BIC a Fit vs. noise. Vypocet kvality beze zmeny |
| R3 (kolo 6, krok 5) | opraveno v README (Per-point residual check): reportovat fit bez bodu s jeho +/- a vynechanymi frekvencemi; rozdil variant neni nejistota ani rozpeti; posun nad stderr = nevypsane body mohou tahnout dal; robustni jsou veliciny, ktere posun nemeni |
| R4 (kolo 6, krok 5) | opraveno v README (`--lambda-probe`): STABLE vylucuje artefakt regularizace, ne slouceni dvou procesu; na T14 zustane 1 pik i pri lambda 1e-8 |

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
