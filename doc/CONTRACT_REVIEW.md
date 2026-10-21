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

## Kolo 3: fitting/series_resistance.py (2026-10-21)

Kontrakt: docstring `eis_analysis/fitting/series_resistance.py`. Reviewer:
Opus, jeden beh (7 min, 100k tokenu, 30 volani nastroju), piskoviste
`~/coldreview/series_resistance/` = cely balik (slouceny refit vola
`_fit_local`); cist smel jen modul, zbytek jako cernou skrinku. Debrief:
precten jen `series_resistance.py`, zbytek pres `help`/`inspect`; 16
vlastnich skriptu v `review/`. Nalezy reprodukovany nezavisle vlastnimi
spektry (1 % sum, 0.1 Hz - 100 kHz) proti modulu v repozitari.

Autor kontraktu predem cekal F3 (kontrakt slibuje vyrazeni kandidata
s nekonecnym RSS, kod ho nekontroluje).

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | (A) -> B | hlavicka: profil Rs plochy, protoze data urcuji jen Rs + R_k | dva bliske oblouky (K 40/40 Ohm, tau 1e-4 a 1.6e-4): po odebrani K1 ho prevezme K2, slouceny Rs zustane 9.98 (ne Rs + R_k), dchi2 0.02-1.0, CI Rs 95 % [1e-4, 10.1] misto [9.83, 10.1] | bezny: obvod s obloukem navic. Profil Rs je ale opravdu plochy: K1 s tau 1e-9 s (v mezich) je rychly rezistor, Rs 10 -> 0.5 Ohm za dchi2 0.89-1.04 (overeno). Rozsireni je spravne, vysvetleni ve varovani ne: pricinou je nadbytecny oblouk |
| F2 | (A) -> B | remiza: prvni v poradi vyrazu; `arc_index` R oblouku | K a (R\|C) maji stejnou impedanci; po odebrani 50Ohm (R\|C) prevezme jeho roli K a dchi2 se shoduje s odebranim K na 10 cislic; jmenovan a rozsiren R2 50 Ohm | bezny (nas beh: 1 z 3 seminek, v serii (R\|C)(R\|C) 1 z 6). Dva bloky se stejnym tau v serii jsou jeden oblouk, takze i R2 v [0, 50] je plochy profil |
| F3 | A | kandidat s nekonecnym RSS vypadne, vyhrava nejmensi konecne dchi2 | NaN dchi2 prvniho kandidata se stane `best`, `dchi2 < nan` pak nikdy neplati: None misto platneho rozdeleni | jen simulovane (prvni refit vrati NaN); prirozeny vstup nenalezen |
| F4 | B | "with ... a G ... it returns None" | pocita jen `type(node) is R`; G vedle jednoho volneho R se ignoruje a R0 se rozsiri | obvod R-G-(R\|C), neobvykly |
| F5 | B/C | zaporne dchi2: fit neni minimum | ve skutecne degenerovanem pripade je dchi2 +-1e-5 az 1e-13 a znamenko nahodne: varovani "fits better (dchi2 = -0.00), so this fit is not the least-squares minimum" ve 2 z 6 seminek | bezny, presne pripad, pro ktery modul vznikl; cisla nezmeni |
| F6 | C | `series_split` se prepise | varovani se pri opakovanem volani pridava znovu | zadny: kazdy volajici vola jednou (multistart fituje pres `_fit_local`) |
| F7 | C | - | `merged.circuit` muze obsahovat prazdny Series (vypis `R(89.74) -  - (...)`) | jen vypis verejneho pole |
| - | C | - | `param_labels` jako numpy pole: `ValueError` | zadny (List[str]) |

Overeno, ze plati (reviewer a reprodukce): fit se nezmeni (snimek vsech
poli), kde hleda (vnoreny Series, Series v Parallel, fixni R, top-level R
nebo Parallel), kandidati vcetne (C|R), (G|C) a K, `merged_index`
a `arc_index` pred i za Rs, dchi2 meri objektiv fitu (uniform, sqrt,
modulus, proportional), invariance vuci skale vah, dof <= 0, NaN v datech
a parametrech bez vyjimky, orez na horni mez Rs, monotonie `widens`,
`diagnostics` None, CI ve varovani = `params_ci_95`.

F1 a F2 jako A neobstaly: rozsireni CI odpovida profilu obvodu, ktery ma
oblouk navic. Chyba je ve vysvetleni (hlavicka i varovani znaji jen
mechanismus "oblouk je konstantni odpor").

Provedeno 2026-10-21 (jedna zmena s testem v `tests/test_series_resistance.py`):

- Kod: F3 - kandidat s dchi2, ktere neni konecne, vypadne; F4 - G mezi
  cleny hlavniho retezce dava None; F5 - `SPLIT_TIE_DCHI2 = 0.01`:
  v jejim okoli nuly varovani rika "fits as well", "not the least-squares
  minimum" az pod -0.01.
- F1 + F2: vypocet beze zmeny, nove zneni varovani ("the fit without the
  arc of R1 ...; that arc can act as a plain resistance in series, or
  another arc takes its place") a druhy mechanismus v hlavicce; kontrakt:
  dva oblouky, ktere si mohou prevzit roli, daji sloucene fity shodne do
  konvergence refitu a vyhrat muze kterykoli.
- F6, F7 a poznamka C bez akce.
- Po /code-review: jen volne G dava None (pevne je konstanta jako pevne R);
  refit s parametry, ktere nejsou konecne, nevypisuje RuntimeWarning
  (`np.errstate`); v pasmu remizy varovani neuvadi cislo, mimo nej
  `.3g` (driv -0.008 i -0.012 jako "-0.01" s opacnym zaverem); varovani
  rika primo "the circuit has an arc too many", README radi oblouk vypustit
  misto hlaseni souctu; test vybira kandidata s NaN podle startu, ne podle
  poradi volani.

## Kolo 4: drt/significance.py (2026-10-21)

Kontrakt: docstring `eis_analysis/drt/significance.py`. Reviewer: Opus,
jeden beh (6.5 min, 92k tokenu, 22 volani nastroju), piskoviste
`~/coldreview/significance/` = cely balik; cist smel jen modul, zbytek
jako cernou skrinku. Debrief: precten jen `significance.py`, zbytek pres
`help`/`inspect`; 13 vlastnich skriptu v `review/`. F1 reprodukovan
nezavisle vlastnim spektrem pres `calculate_drt` proti modulu v repozitari.

Autor kontraktu predem bez pevneho ocekavani; podezrival ploche vrcholy
(find_peaks bere prostredni index, `_peak_basins` argmin) a stejne vysoka
maxima (netestovana, kontrakt to uvadi zamerne). Ani jedno nalez nedalo.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | hlavicka: "gamma must not increase from the neighbour's apex through the candidate's basin"; kontrakt: usek do nejnizsiho bodu povodi vpravo od j (nejvzdalenejsi ze stejnych minim) | povodi z `_peak_basins` je polootevrene a bod udoli patri pravemu sousedovi (u useku nul prvni nula). S vyssim sousedem vlevo usek konci bod pred udolim, u udoli hned za vrcholem primo na vrcholu j; refit pak posune maximum na prvni neomezeny bod a rameno nevynuti. S vyssim sousedem vpravo usek udoli i nuly obsahuje: test neni symetricky | bezny. Vlastni spektrum R10-(R100\|C1e-6)-(R5\|C4e-4)-(R3\|C0.1), 71 bodu 1e5-1e-2 Hz, 0.2 % sum, `calculate_drt` default: skutecny proces 5 Ohm (tau 2 ms) zamitnut ve 3 z 10 seminek (dchi2 12-14), 3 Ohm (tau 0.3 s) v 1 z 10 (10.2); s usekem pres nuly dchi2 69-148. I prijate piky jen 19-30 misto ~100 |
| F2 | C | sigma: noise_sigma, floor na eps | `max(nan, eps)` = NaN: vsechna dchi2 NaN, vse vyznamne; zaporne noise_sigma tise eps (dchi2 ~1e28) | zadny: volajici predavaji vypocitane konecne sigma (gamma z NNLS je konecne) |
| F3 | C | kontrakt mlci | `threshold` NaN: vse vyznamne | zadny: volajici pouzivaji konstantu |
| - | B | "A gamma with no interior maximum (empty, ...)" | prazdna gamma je nedosazitelna: `_build_drt_matrices` pro n_tau < 2 vyhodi ValueError | jen kontrakt |

Prototyp opravy F1 (usek se sestupem konci v nejvzdalenejsim minimu pred
dalsim kandidatem, zrcadlove k useku se vzestupem; monkeypatch, repozitar
beze zmeny): 7 syntetik (3RC, 1RC, 1 ZARC n 0.8, 2 ZARC 0.5 a 1 dek.,
2RC 1:100 a 1:1000) x sum 0.1 / 1 / 3 % x 10 seminek, `calculate_drt`
default:

| Spektrum, sum | chybejici ted | chybejici s opravou |
|---|---|---|
| 3RC 0.1 % | 2 | 0 |
| 3RC 1 % | 16 | 12 |
| 2RC 1:100 0.1 % | 5 | 0 |
| ostatni | beze zmeny | beze zmeny |

Falesne piky: 1 (1RC 1 %) pred i po, jinak 0. Kalibrace prahu 16
(`DRT_PEAK_DCHI2_MIN`, 248 behu) byla merena s chybou F1; kalibracni
skript je na druhem systemu.

Overeno, ze plati (reviewer): vstupy nezmeneny, zarovnani a usporadani
poli, nejvyssi a stejne vysoky kandidat vzdy vyznamny (inf), monotonie
v prahu (62 prahu, +-inf), dchi2 ~ 1/sigma^2, nezavislost na jednotce Z
(c = 1e-3 .. 1e6, s L i bez, rel. 1e-11), nulova/konstantni/monotonni
gamma bez vyjimky, chi2_full = Boukamp, sigma = sqrt(chi2/2N), refit na
uz monotonnim useku reprodukuje chi2_full na 1e-15, volba q (vyssi sedlo,
pri shode levy), nnls na limitu iteraci -> NaN a pik zachovan, zaporne
dchi2 (-0.0018).

Provedeno 2026-10-21 (jedna zmena s testem v `tests/test_drt_peak_significance.py`):

- F1: usek konci v nejvzdalenejsim minimu pred hranici; hranice je dalsi
  kandidat, ktery je aspon tak vysoky jako j nebo neni zamitnuty (jinak
  konec mrizky). Kandidati se testuji od nejnizsiho, takze mensi maxima
  jsou rozhodnuta driv a zamitnuta lezi v udoli. Testuji se oba vyssi
  sousede (napred ten za vyssim sedlem, druhy jen kdyz prvni neni pod
  prahem) a plati mensi dchi2. `_peak_basins` modul uz nepouziva.
- Cesta k tomu (vse mereno, `example/` pres CLI, 3RC a 7 syntetik):
  1. hranice = sousedni kandidat: `example_eis_data.csv` 1 -> 2 piky,
     novy (tau 4.8e-5 s, dchi2 17.2) je hrot mezi nulami testovany proti
     levemu HF hrotu (pravidlo "pri shode sedel levy"), proti vyssimu
     pravemu je ramenem (-0.1); dosud maskovano F1.
  2. pri shode sedel oba sousede: `example/` beze zmeny, 3RC 1 % 16 -> 12
     chybejicich. /code-review: shoda jen pri presne rovnosti floatu
     (skok) a drobny lalok v udoli zastavi usek jako driv (overeno: u 3RC
     1 % je za kazdym zamitnutym skutecnym pikem lalok vysky 0.4-4 proti
     12-14).
  3. hranice = dalsi vyssi kandidat: 12 -> 8, ale do useku spadne i mensi
     skutecny pik a jeho struktura se pripise kandidatovi
     (`example_eis_data`: falesny pik 0.38 s s dchi2 105).
  4. hranice = dalsi vyssi nebo vyznamny (vysledna): 3RC 0.1 % 2 -> 0,
     1 % 16 -> 0; 2RC 1:100 0.1 % 5 -> 0; zadny falesny pik navic (1RC
     1 %: 1 pred i po). `example_eis_data` 1 -> 2: proces R0||Q0 (tau
     0.08 s, R 95 kOhm z 100), ktery se se zavedenim testu ztratil;
     diffusion_Ws_noise1 1 -> 2 a diffusion_RpWs 2 -> 3: dalsi mod
     konecneho Warburgu (tau/9, tau/25). Cena: 151 refitu jako drive
     (druhy soused jen u prijatych).
- Kontrakt: hranice, poradi a oba sousede popsane primo; NaN jen kdyz
  zadny refit neuspel; monotonie v prahu a presne skalovani 1/sigma^2
  uz neplati (zamitnuti mensich maxim zavisi na prahu a sigma), uvedeno
  vyslovne; symetrie mimo plochy vrchol sude sirky; F2 + F3 jako zaruka
  volajiciho; prazdna gamma z klauzule vyskrtnuta.
- /code-review dale: komentar o konci mrizky, hlavicka bez "basin",
  docstring testu (zamitnuty byl proces 3 Ohm). Selhani refitu jednoho
  souseda: zamitne jen mensi dchi2 druheho, ktere je pod prahem, coz
  plati pro minimum tak jako tak.
- Kalibrace `DRT_PEAK_DCHI2_MIN` (248 behu) byla merena pred opravou;
  zopakovano 2026-10-21 novym `tests/calibrate_peak_threshold.py` (puvodni
  skript je na druhem systemu; sada podle popisu v `config.py` + 3RC z F1):
  pri 16 chybi 77 pred opravou 82 (3RC 0.1 % 3 -> 0, 1 % 19 -> 17),
  falesne 0 pred i po; nejvyssi falesne dchi2 5.5 -> 7.7. Prah 16 zustava.

## Kolo 5: rinf_estimation/estimate.py (2026-10-21)

Kontrakt: docstring `eis_analysis/rinf_estimation/estimate.py`. Reviewer:
Opus, jeden beh (5 min, 89k tokenu, 19 volani nastroju), piskoviste
`~/coldreview/rinf_estimation/` = cely balik (fit pres
`fit_equivalent_circuit`); cist smel jen modul, zbytek jako cernou
skrinku. Debrief: precten jen `estimate.py`, `hf_bound` a fit pres
`help`/volani; 16 vlastnich skriptu v `review/`. Vsechny nalezy
reprodukovany nezavisle proti modulu v repozitari (skript v scratchpadu
sezeni).

Autor kontraktu predem cekal F3 (`rel = stderr / R_fit` na Python floatech,
dolni mez R_s je 0) a prijeti nezkonvergovaneho fitu jen podle stderr
(kontrakt zamerne mlci; reviewer to nehlasil). `fit.series_split` kontrakt
vyslovne nepouziva (otevrene rozhodnuti, menilo by --ri-fit).

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | poradi bodu meni jen poradi f_window / Z_window | `hf_bound` bere z bodu se stejnou nejvyssi frekvenci ten, ktery je v poradi prvni: R_inf_hf, R_inf mezi, R_inf_upper i R_inf_range zavisi na poradi vstupu (otevreny oblouk, f_max dvakrat s Re +1 Ohm: R_inf 109.47 nebo 110.47 podle permutace). Kontrakt nerika, ktery z duplikatu plati; mez je kazdy, nejtesnejsi je minimum | `hf_bound` sdili i vychozi R_inf DRT (`hf_r_inf`). Jen pri opakovane nejvyssi kapacitni frekvenci; `example/` duplikaty nema |
| F2 | B | zaloha: varovani jmenuje mez a duvod | bez bodu s Im <= 0 vraci `hf_bound` f_max a Re(Z) tam se hlasi jako "HF upper bound" bez poznamky, R_inf_range (0, Re) skutecne R_s neobsahuje (Rs 10, Re 9.0..9.9, Im > 0: rozsah (0, 9) pri 4 bodech, (0, 9.26) s fitem pri 8) | jen cele spektrum induktivni (kabel, zkrat); `example/` zadne |
| F3 | C | jen ValueError; R_s = 0 -> zaloha | `rel = stderr / R_fit` na Python floatech: R_s = 0 dava ZeroDivisionError (fit podvrzeny mockem) | zadny pozorovany: 14 spekter se skutecnym R_s <= 0, fit vzdy > 0 (nejmene 2.5e-34) |
| F4 | C | jen ValueError; \|Z\| = 0 -> zaloha | \|Z\| ~1e-320: `min|Z| / 1e6` podtece na 0, `_window_range` deli nulou (ZeroDivisionError) | zadny: nefyzikalni |
| F5 | B/C | invariance vuci jednotce, bez rozsahu | do c = 1e-16 rel. 1e-15; od c <= 1e-17 R_inf +1.3 %, L -94 az -100 % (absolutni tolerance ve fitu); c = 1e-200 zaloha, c >= 1e200 fit selze a zaloha | zadny: Z pod ~1e-15 Ohm |
| F6 | C | okno f >= f_max / 100 | deleni v plovouci carce: 1.1 / 100 = 0.011000000000000001, bod 0.011 Hz vypadne (taky horni dekada `_window_spread`) | zanedbatelny: o jeden bod na libovolne hranici; mocniny 10 presne |
| - | C | text "Ignored k non-finite point(s)" | pocita i konecne body s f <= 0 (text predepsal kontrakt) | jen hlaska |
| - | C | - | R_hf = -0.0: R_inf -0.0, hlaska "Re(Z) = -0 Ohm" | jen hlaska |
| - | B | sum: "RMS ... scaled by 2N / (2N - 5)" | skaluje se rozptyl (RMS odmocninou) | jen kontrakt |

Overeno, ze plati (reviewer): ValueError pro tvar, prazdny vstup, vse NaN,
jen f <= 0; seznamy, 2-D i 0-d vstup; vyrazeni NaN/inf/f <= 0 s poctem;
vsechny vetve zalohy (< 5 bodu, \|Z\| = 0 v okne i mimo nej, cele Z nulove,
stderr NaN/inf, zaporne R_hf, poznamka o induktivnim vrcholu); invarianty
(konecne a >= 0, rozsah obsahuje R_inf, R_inf <= R_inf_upper u zalohy, L 0
mimo fit, rlq_fit -> stderr <= 5 %) na vsech pripadech bez vyjimky;
R_inf_upper s fitem presne podle vzorce (huste, ridke spektrum s pravidlem
3 nejblizsich, f_hf v induktivnim okne); R_inf_window_spread i jeho nulove
vetve; invariance c = 1e-16 .. 1e100; permutace (rel. 1e-15); vstupy
nezmeneny; determinismus; f az 1e300 / 1e-300.
