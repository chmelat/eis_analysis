# Kontrakt v hlavicce modulu a cold review proti nemu

Inspirace: driftlet 1.2.2 az 1.3.0 (`doc/DRIFTLET_COMPARISON1.md`, bod 6).
Tam kazdy modul dostal kontrakt a cerstvy agent ho cetl jen proti nemu;
skoro kazde z deseti takovych review naslo tichy spatny vysledek.

## Co je kontrakt

Sekce `CONTRACT` v docstringu modulu. Neni to docstring jedne funkce, ale
overitelny popis chovani celeho modulu:

- **Co modul cte**, a jen to: kdo ho vola, co mu predava (tvar, jednotky,
  poradi, co volajici zarucuji), co pouziva z jinych modulu.
- **Co vraci:** pole a typy, jednotky, kdy vyhodi vyjimku a kdy misto toho
  hlasi problem ve vysledku.
- **Co musi platit:** vztahy mezi vstupem a vystupem (invariance,
  presnost) a hlavne okrajove pripady: prazdna data, NaN, jeden bod,
  nulovy sloupec, fixni parametry, hodnota na mezi.

Kontrakt popisuje **zamyslene** chovani (ze statistiky, fyziky,
docstringu), ne prepis kodu. Kde autor pri psani kontraktu vi, ze kod dela
neco jineho, kontrakt pise, co ma platit, a odchylku najde review.

**Pravidlo (CLAUDE.md, Code Organization):** modul s kontraktem se meni
spolu s kontraktem a testem v jednom commitu. Bez toho kontrakt zastara
jako kazda dokumentace.

## Cold review

Cerstvy agent bez kontextu projektu dostane modul s kontraktem a ukol:
"najdi mista, kde kod kontraktu neodpovida, nebo kde kontrakt mlci
o pripadu, ktery kod resi pochybne". Necte diff ani historii; cte cely
modul a porovnava, co kod dela, s tim, co hlavicka slibuje.

Funguje proto, ze kontrakt pojmenuje i okrajove pripady, na ktere autor
pri psani kodu nemyslel, a agent pak hleda cilene ("kontrakt rika X, plati
to i pro Y?").

**Piskoviste:** `~/coldreview/<modul>/` mimo repozitar, kopie modulu
(a modulu, ktere kontrakt jmenuje, pokud bez nich nejde spustit). Agent
smi cist jen tam a spoustet tam vlastni skripty
(`~/venvs/eis/bin/python`). Izolace je mekka (instrukce v zadani), proto
debrief obsahuje seznam prectenych souboru.

**Zadani** (jednotne, menit jen cestu): nalezy jen potvrzene behem
konkretniho vstupu; ke kazdemu klauzule kontraktu (nebo "contract
silent"), misto v kodu, co kod dela misto toho, reprodukcni skript
s vystupem, zavaznost. Na konci debrief: prectene soubory, spustene
skripty, co nesel otestovat.

**Zavaznost:**

- **A** tichy spatny vysledek, ktery volajici pouzije (konecny stderr tam,
  kde ma byt inf, spatna kovariance)
- **B** spatny nebo nejednoznacny kontrakt
- **C** vada bez vlivu na cisla (hlaska, robustnost, styl)

**Overeni a triage** (jako princip 6 v `doc/COLD_AGENT_PLAN.md`): kazdy
nalez reprodukovat vlastnim skriptem; co se nereprodukuje, je chyba
agenta. Overit i dopad pres volajici: nastane ten vstup v knihovne vubec?
Nalezy se zapisuji sem, opravy az po schvaleni, kazda samostatne a s
testem.

## Vztah k ostatnim kontrolam

| Kontrola | Hleda | Cte |
|---|---|---|
| Cold agents (`doc/COLD_AGENT.md`) | chyby API a dokumentace | jen dokumentaci, pouziva knihovnu zvenku |
| `/code-review` | chyby ve zmene | diff |
| Stress test (`doc/STRESS_TEST.md`) | chyby algoritmu na nahodnych spektrech | vysledky, jen pres invarianty |
| **Cold review proti kontraktu** | chyby kodu proti zamyslenemu chovani | cely modul zevnitr |

Kandidati (male, matematicke moduly, jejichz vystup se pouziva dal):
`fitting/covariance.py`, `fitting/series_resistance.py`,
`rinf_estimation/`, `validation/` (KK).

## Kolo 1: fitting/covariance.py (2026-10-20)

Kontrakt: docstring `eis_analysis/fitting/covariance.py`. Reviewer: Opus,
jeden beh (2.5 min, 64k tokenu, 6 volani nastroju), piskoviste
`~/coldreview/covariance/`. Debrief: precten jen `covariance.py`, tri
vlastni skripty v piskovisti. Vsechny nalezy reprodukovany nezavisle
(`~/coldreview/verify/verify.py` proti modulu v repozitari).

Kontrakt byl psan podle zamysleneho chovani; F1 autor kontraktu cekal
(klauzule dof <= 0 je zamerne jina nez kod), ostatni ne.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | dof <= 0: stderr inf | `dof = max(n_res - n_free, 1)`: s^2 = RSS/1, konecny stderr, `dof` 1, bez varovani (J 2x2: stderr 0.22, 0.11) | fit s velmi malo body pres API (hlavni fit pocet bodu nehlida; okno R_inf ma >= 5 bodu, tj. 10 rezidui na 5 parametru) |
| F2 | A | rank deficient: parametr v null space ma stderr inf | siroky J (n_res < n_free): `svd(full_matrices=False)` vrati jen n_res radku Vt, null space (n_free - n_res smeru) chybi, vsechny stderr konecne, varovani "Non-identifiable parameter(s) []" | stejne jako F1 |
| F3 | A | NaN/inf v reziduich: cov None, stderr inf, varovani | rezidua se nekontroluji: NaN dava cov NaN, stderr NaN, rank 2, wc True, bez varovani | jen pokud optimalizator vrati nekonecna rezidua; CI je pak (-inf, inf), cov a priznaky ale lzou |
| F4 | C | totez pro jakobian | spravny tvar vysledku jen proto, ze LAPACK na NaN/inf selze; varovani "SVD did not converge", navic RuntimeWarning a `logger.error` | jako F3 |
| F5 | B/C | "never raises on numerical trouble" | prazdna rezidua: `IndexError` na `S[-1]` | volajici nikdy (vzdy >= 1 bod) |
| F6 | B | inf radek neidentifikovatelneho vs. 0 radek fixniho | kontrakt si odporuje; kod dava na prusecik 0 (rozumne) | jen kontrakt |
| F7 | C | log CI: preteceni -> (0, inf) | test preteceni jen na f, ne na p*f: pri velkem p (1e10) RuntimeWarning a dolni mez 3e-298 misto 0 | parametr nedeterminovany na 700 % a vic; cisla prakticky stejna |
| N1 | - | zaporna variance -> varovani | vetev je nedosazitelna: diagonala je soucet nezapornych clenu i v plovouci carce | mrtvy kod |
| N2 | C | kontrakt mlci | `confidence_level = 1.0`: t = inf, stderr 0 dava NaN misto (p, p) | volajici pouzivaji 0.95 a 0.99 |

Overeno, ze plati (reviewer a reprodukce): invariance vuci jednotkam
(relativni chyba 3.6e-15 pri skalovani sloupce 1e6), plny rank presne
s^2 (J^T J)^-1, pseudo-inverze identifikovatelnych parametru rovna
presnemu profilu, nulovy sloupec neidentifikovatelny, indexy v plnem
prostoru, fixni parametry 0, vsechny klauzule CI krome F7.

Navrh oprav:

1. **Provedeno 2026-10-20:** F1 + F2 + F5 jednou vetvi na zacatku
   `compute_covariance_matrix`: `n_res <= n_free` -> stderr a cov volnych
   parametru inf, skutecny dof (<= 0), condition_number inf, rank 0,
   varovani; test `test_too_few_residuals_gives_inf_stderr` (2x2, 2x3,
   prazdna rezidua). Ctyri starsi testy pouzivaly jakobian 2x2
   (dof 0, drzeny na 1); dostaly nulovy radek, ktery J^T J nemeni.

2. **Provedeno 2026-10-20 (po /code-review oprav 1):** F3 + F4 kontrolou
   `np.isfinite` jakobianu a rezidui uplne na zacatku (pred vetvi pro malo
   bodu): cov None, kazdy stderr inf, varovani "Non-finite jacobian or
   residuals"; test `test_non_finite_input_gives_no_covariance`. N1: mrtva
   vetev zaporne variance smazana i s klauzuli. Ve stejne zmene: jeden
   pomocnik `_full` (rozbaleni do vsech parametru) misto tri kopii,
   `_failed` pro selhani; s `fixed_params=None` maji vsechny cesty pole
   delky n_free (dve mely n_params); `CovarianceResult.dof` bez vychozi
   hodnoty 1. Soubor 513 -> 471 radku.

Zbyle body (neprovedeno, bez dopadu na bezna data):

- F6: kontrakt - fixni ma prednost (0 i na prusecicich).
- F7: test preteceni na `log p + z`, ne na `z`.
- N2: kontrakt - `confidence_level` v (0, 1).
- `cli/handlers/voigt_fit.py` drzi `_dof = max(2N - k, 1)`, stejne jako
  driv kovariance (nalez /code-review); Voigt chain je okrajovy, neopraveno.
- Parametr `n_params` funkce `compute_covariance_matrix` uz nic neurcuje
  (velikost dava `fixed_params`, bez nej n_free); ponechan kvuli API.

## Kolo 2: fitting/time_constants.py (2026-10-21)

Kontrakt: docstring `eis_analysis/fitting/time_constants.py` (K zamerne
mimo soucet odporu). Reviewer: Opus, jeden beh (5 min, 90k tokenu, 16
volani nastroju), piskoviste `~/coldreview/time_constants/` (modul,
`circuit_builder.py`, `circuit_elements/`; `FitResult` si reviewer nahradil
`SimpleNamespace`). Debrief: cetl jen piskoviste, 11 vlastnich skriptu
v `review/`. Nalezy reprodukovany nezavisle se skutecnym `FitResult` proti
modulu v repozitari (skript v scratchpadu sezeni), dosah overen pres
`covariance.py`, `bounds.py` a `block_order.py`.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | blok, ktery zjevne neni PSD: None | tau kontroluje jen `var_log < 0`, soucet jen celkovy rozptyl a diagonalu: korelace 3 dava konecny tau_stderr, blok [[1, 2], [2, 1]] stderr 2.45 a correlation 2.0 | zadny: `cov` z `covariance.py` je s^2 V S^-2 V^T, PSD z konstrukce |
| F2 | A | zaporny rozptyl jen zaokrouhlenim: stderr 0 | tolerance je jen u souctu; tau pri korelaci presne -1 vraci None ('covariance not usable') v 780 z 2000 nahodnych pripadu, soucet 0 ve 2000 z 2000 | prakticky zadny: potrebuje var(ln tau) ~1e-16 clenu; nekonzistence obou cest |
| F3 | A | nekonecny parametr: zadne tau | n = inf projde `n > 0`: tau 1 s, stderr 0 | zadny: n v (0.3, 1) z `PARAMETER_BOUNDS` |
| F4 | B | correlation "Sigma_12 scaled by the parts" | kod bere (D Sigma D)_12, korelaci R1 a 1/G2 (spravne); doslova Sigma_12 je cov(R1, G2), opacne znamenko | jen kontrakt |
| F5 | B | vse fixni: stderr 0, vs. inf v bloku: None | pri vsem fixnim a cele cov inf None | zadny: `covariance.py` dava fixnim radkum 0 i ve vetvi s malo rezidui (overeno) |
| F6 | B | `_arc_resistor`: None, kdyz dva prvky nejsou par | delku nekontroluje, (R\|C\|Q) da (0, +1) | zadny: `walk` i `block_order._block_arc` kontroluji 2 prvky predem |
| F7 | B | nekonecne hodnoty: zadny soucet | G = inf: 1/G = 0 konecne, soucet hlasen, rozptyl clenu G vynulovan | zadny: G v (0, 1e4) |
| F8 | C | bez RuntimeWarning | `_arc_tau` primo: overflow / divide warning | zadny: `block_order` vola v `np.errstate` |
| F9 | C | - | subnormalni fixni C: gradient inf, inf*0 = NaN, poznamka 'covariance not usable' misto duvodu | zadny (meze C) |
| F10 | C | correlation v [-1, 1] | pasmo tolerance 1e-9*trace pusti -1.0000000004 | zadny (PSD) |
| F11 | C | contract silent | kratsi `bound_status` = vnitrni body, polozka None da poznamku "near its None bound" | zadny: `build_bound_status` dava plnou delku |

Overeno, ze plati (reviewer a reprodukce): nezavisle R a C, soucet
nekorelovanych clenu, korelace -1 se stejnymi castmi dava 0, ekvivalence
R/G s transformovanou kovarianci a prohozeni prvku oblouku (1e-16),
gradient proti centralnim diferencim (R|Q, Q|G, C|G), `cov None`, mez
jen u parametru oblouku, fixni = 0, inf radek mimo blok, `param_labels
None`, vstupy nezmeneny, topologie (vnoreny Series, Series v Parallel,
(R|C|Q), Parallel v Parallel, K, Randles), extremy pres verejne API bez
vyjimky a varovani.

Zadny nalez neni dosazitelny z fitu knihovny: A plati jen pro vstup, ktery
volajici nikdy nepredaji.

Provedeno 2026-10-21 (jedna zmena s testem `tests/test_time_constants.py`):

- Kod: F2 - `_delta_variance` s toleranci zaokrouhleni pro tau i soucet
  (driv jen soucet); F3 + F7 - parametr oblouku nebo clenu souctu, ktery
  neni konecny, nedava tau ani soucet.
- Kontrakt: F1 + F10 - `cov` je PSD jako zaruka volajiciho (klauzule
  o ne-PSD bloku vypustena); F4 - korelace jako (D Sigma D)_12; F5 - inf
  radek jen u volneho parametru, fixni vzdy 0; F6 + F8 - predpoklady
  privatnich funkci (dvouprvkovy Parallel, `np.errstate` u volajiciho);
  F11 - `bound_status` ma polozku pro kazdy parametr.
- F9 bez akce (subnormalni C, mimo meze).
