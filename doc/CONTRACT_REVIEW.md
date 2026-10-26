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

Provedeno 2026-10-21 (jedna zmena s testy v `tests/test_rinf_estimation.py`;
F2 jen varovani, rozhodnuti uzivatele):

- F1: `hf_bound` (`utils/impedance.py`, sdileny s `hf_r_inf` v DRT) radi
  `lexsort` podle klesajici f a pri shode rostouciho Re: z bodu na stejne
  frekvenci plati nejnizsi Re. Bez duplikatu beze zmeny.
- F2: bez bodu s Im <= 0 varovani "No point with Im(Z) <= 0: Re(Z) at
  f_max is not an upper bound of R_inf ..." u kazde metody; R_inf_upper
  a R_inf_range se pocitaji jako dosud, kontrakt rika, ze R_s drzet
  nemusi.
- F3: `rel` je inf pro R_s <= 0 (zaloha s varovanim); `_window_spread`
  prijme jen R_s > 0.
- F4, F5, F6: jen kontrakt (|Z| nad ~1e-300 Ohm, invariance pro |Z|
  ~1e-15 .. 1e100 Ohm, hranice okna v plovouci carce). Hlaska o vyrazenych
  bodech jmenuje i f <= 0; kontrakt u sumu rika rozptyl, ne RMS. -0.0 bez
  akce.
- Overeno: 3 nove testy bez opravy padaji (ZeroDivisionError, poradi,
  chybejici varovani); `pytest tests/` 863 passed, stress smoke proti
  baseline prosel.
- /code-review (10 nalezu): `hf_r_inf` (dve horni kapacitni body) a
  `hf_median` (shoda na hranici n bodu) razeny stejne jako `hf_bound`,
  jinak DRT R_inf dal zavisel na poradi (test: 3 body na f_max, median vs.
  mez podle permutace); `hf_bound` a spol. prevadi f na float (uint by
  pretekl). `RinfResult.hf_is_bound`: bez bodu s Im <= 0 zalozni hlaska,
  CLI (`--ri-fit`, `--local-exponent`) i graf nerikaji "upper bound";
  poznamka o preskocenem induktivnim bodu i pri shode f s f_hf; CLI
  netiskne procento pri R_s = 0. `_window_spread`: R_s na mezi 0
  v uzsim okne je chyba modelu, vraci cely posun R_fit (moje F3 oprava
  vracela 0 a zuzila rozsah). Docstring Raises. Bez akce: horni konec
  rozsahu prijateho fitu na induktivnim spektru dal oriznut R_inf_upper
  (rozhodnuti u F2: jen varovani, kontrakt rika, ze R_s drzet nemusi).
- Druhe /code-review (8 nalezu): vetev R_s <= 0 v `_window_spread` je
  mrtvy kod (least_squares 'trf' drzi parametr uvnitr mezi: R_s na mezi
  vraci 1e-11 az 2e-16, ne 0); `_window_spread` vracen do puvodniho
  stavu. Otevrene: uzsi fit s R_s na mezi dava sirku 0, i kdyz je to
  znamka chyby modelu; rozsireni by zmenilo R_inf_range, rozhodnout az po
  mereni na stress spektrech. `hf_median` bere pri shode nizsi Re jako
  `hf_bound` (sdilene `_top_down`); graf (titulek, R_inf_upper, "Re(Z) <
  0") a varovani `--local-exponent` bez "bound" pri hf_is_bound False;
  test induktivniho spektra vynucuje zalohu se 4 i 8 body. Bez akce:
  rozsah oriznuty R_inf_upper bez meze (rozhodnuti u F2), text varovani
  pri R_s presne 0 (nedosazitelne).

## Kolo 6: analysis/local_exponent.py (2026-10-22)

Kontrakt: docstring `eis_analysis/analysis/local_exponent.py`. Reviewer:
Opus, jeden beh (11 min, 113k tokenu, 33 volani nastroju), piskoviste
`~/coldreview/local_exponent/` = cely balik; cist smel jen modul, zbytek
jako cernou skrinku. Debrief: precten jen `local_exponent.py`,
`estimate_rinf` a `validate_eis_data` pres `help`; 15 vlastnich skriptu
v `review/`. F1-F4 reprodukovany nezavisle (skript v scratchpadu sezeni)
proti modulu v repozitari, F5 a F7 skripty reviewera proti repozitari.

Autor kontraktu predem cekal F1 (deleni nulou mimo `np.errstate`), F2
(okno s jedinou frekvenci), F3 (nestabilni `argsort`) a F5 (nejistota jen
z koncu rozsahu). F4 ne.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | C | bez RuntimeWarning; Z - R_inf presne 0 nema exponent | `1 / (Z - r_inf - j w L)` mimo `np.errstate`: "divide by zero" a "invalid value"; n v tom bode NaN spravne | CLI: zalozni R_inf = Re Z(f_max) pri Im presne 0 a L = 0 (overeno pres `estimate_rinf`, metoda hf_bound) |
| F2 | C | okno s jedinou frekvenci: NaN, bez RuntimeWarning | `sum(dx**2)` = 0, 0/0 = NaN: 3 varovani na bod (R_inf a oba konce) | 4 opakovani jedne frekvence a nic do +-0.5 dekady; neobvykle |
| F3 | B | body se stejnou f v poradi vstupu | `np.argsort` bez `kind='stable'`: dva behy po 61 bodech, 29 z 61 paru prohozeno (amd64, numpy 2.5.4; pri 16 bodech 0) | dva sweepy v jednom souboru; n kazdeho bodu spravne pro svou f, prirazeni k behu ne |
| F4 | B | invariance vuci skale frekvenci (s L / c) | hrana okna `<= window / 2` padne na presne logaritmicke mrizce (10 b/dek) na bod; zaokrouhleni log10 rozhodne, zda je uvnitr: f x 1e3 meni n az o 0.019-0.028, 5-6 bodu platnych v obou behu s \|dn\| nad nejistotou; bez skalovani 60 oken s 11 body a 13 s 10 | CLI ne: frekvence zaokrouhlene na 6 cislic, invariance plati (1e-15); synteticka data a API |
| F5 | B | nejistota pokryje kazde R v `R_inf_range` | n jen na koncich rozsahu; bod j s Re Z_j uvnitr rozsahu z okna na jednom konci vypadne a pro R tesne pod Re Z_j jde ln Re Y_j do -inf: n 0.800 +- 0.014 platne, n(102.4999) = 0.689 | konstruovany pripad; 110 nahodnych CLI spekter pres `estimate_rinf` bez poruseni |
| F6 | C | jen ValueError | `R_inf=None`/`"10"`/komplexni: TypeError; `R_inf_range` s 1 prvkem IndexError, se 3 prvky treti tise ignorovan; 2-D vstup "need at least 4 points, got 1"; komplexni f ComplexWarning; L 1e305 overflow | jen API |
| F7 | C | kontrakt mlci | dve dvojice opakovani na 1000 a 1000(1 + 1e-9) Hz bez sumu, R_inf = 0 s vychozim rozsahem nulove sirky: n = -999499, nejistota 0, platne, varovani o rozpeti 999500 | jen konstruovany vstup (presne shodne hodnoty, sirka rozsahu 0) |
| - | B | - | docstring tridy (NaN jen pro Re Y <= 0 a malo bodu), Raises bez R_inf/L/rozsahu, "R_inf +-LOCAL_EXPONENT_RINF_REL" bez x \|R_inf\| | jen dokumentace |

Overeno, ze plati (reviewer): ValueError pro prazdny vstup, ruzne delky,
f <= 0, NaN/inf v Z, malo bodu, R_inf/L/rozsah nekonecne, rozsah bez R_inf
(i obraceny); vychozi rozsah (i R_inf = 0 a zaporne); vstupy nezmeneny
(pole i seznam); determinismus; permutace bez duplikatu (2e-16); skala Z,
R_inf, L a rozsahu (1e-15); jeden CPE || C s presnym R_inf a L: n na
4e-13 (R_inf 0 .. 1e4, L do 1e-5, zaporne R_inf), bez varovani; NaN
nejistoty kde n NaN na kterekoli strane; nejistota >= 0; valid -> n
konecne; souhrnne hodnoty jen z platnych bodu, NaN s varovanim bez nich;
varovani o rozpeti 0.31 ano, 0.06 a 0.009 ne.

Provedeno 2026-10-22 (minimalni varianta, rozhodnuti uzivatele: nalezy
jsou malickosti; jedna zmena s testem v `tests/test_local_exponent.py`):

- Kod: F1 + F2 - deleni, logaritmus i regrese v `np.errstate`; F3 -
  `argsort(kind='stable')`. Cisla beze zmeny.
- Kontrakt: F4 - invariance vuci skale f neplati pro bod na hrane okna
  (presna logaritmicka mrizka), merene frekvence to nepostihne; F5 -
  nejistota pokryva konce `R_inf_range`, uvnitr se n(R) nevyhodnocuje
  (mechanismus uveden); F6 - typy vstupu jako zaruka volajiciho; F7 -
  u vychoziho rozsahu nulove sirky plati jen stderr regrese.
- Bez akce: tolerance hrany okna (F4), docstring tridy a Raises.
- Test bez opravy pada (RuntimeWarning, poradi dvojic).
- /code-review (6 nalezu): F2 nebyl opraven jen potlacenim varovani -
  prumer shodnych floatu se zaokrouhli o ~1e-16, dx pak neni 0 a okno
  jedne frekvence (5x 123.4 Hz) dalo n = 0 +- 0 platne (i pred zmenou);
  okno s `ptp(x) == 0` se preskoci, `np.errstate` jen kolem deleni
  a logaritmu. Test pokryva i tento pripad. Kontrakt: zalozni R_inf je
  Re(Z) mereneho bodu (hf_bound), ne nutne nejvyssi f; klauzule
  o skale frekvenci prepsana. Bez akce: test poradi zavisi na tom, ze
  vychozi `argsort` duplikaty prohazuje (na amd64 ano, ARM neovereno).
- Druhe /code-review (7 nalezu): `LOCAL_EXPONENT_MIN_POINTS` pocital body,
  ne frekvence; dva sweepy na mrizce 0.4 dekady daly okno se 2 frekvencemi
  platne (n 0.8 +- 5e-5, bez chyby modelu). Okno ted potrebuje 4 ruzne
  frekvence; to nahrazuje kontrolu `ptp` a pokryva i F7 (1000 a
  1000(1 + 1e-9) Hz jsou dve frekvence). Kontrakt: poradi meni n o
  zaokrouhleni (~1e-16) u duplikatu; docstring tridy jmenuje vsechny
  priciny NaN. Bez akce: f nad ~1e307 (preteceni omega), `over='ignore'`
  skryje nefyzikalni L (1e305), test poradi na ARM.

## Kolo 7: drt/gcv.py (2026-10-23)

Kontrakt: docstring `eis_analysis/drt/gcv.py`. Reviewer: Opus, jeden beh
(8 min, 105k tokenu, 25 volani nastroju), piskoviste `~/coldreview/gcv/`
= cely balik; cist smel jen modul, zbytek jako cernou skrinku. Debrief:
precten jen `gcv.py`; realne systemy (A, b, L) zachytil obalenim
`find_optimal_lambda_hybrid` pri `calculate_drt` (6 systemu, 16x100 az
100x138, i s indukcnim sloupcem a prodlouzenim tau); 20 skriptu
v `review/`. F1-F4, F6-F8 reprodukovany nezavisle (skript v scratchpadu
sezeni) proti modulu v repozitari, druha varianta F1 a F5 skripty
reviewera (balik v piskovisti shodny s repozitarem, `diff -rq`).

Autor kontraktu predem cekal F1 (roh bez kladne krivosti), F2 (`argmax`
na NaN), F3 (vsechna skore inf bez ValueError, i pri spatne delce b), F4
(jemna mrizka bez hrubeho minima), F5 (spatne podminena M), F6 (bez
validace) a F7 (deleni v `np.where`). Necekal F8 a to, ze F1 nastane
pres `calculate_drt` (cisty rezistor). Cekana podlaha 1e-15 invarianci
neporusila.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | krivka bez kladne vnitrni krivosti nema roh a lambdu nezvysi | hybrid bere `argmax` krivosti bez ohledu na znamenko: diagonalni 10x6 A, krivost v okne vsude zaporna (max -0.0118), lambda 5.1e-9 -> 1.1e-7 (22x), `lcurve`, `corner_at_edge` | CLI: cisty rezistor (Z = 100 Ohm), b = 0, krivost vsude 0: `lcurve`, `corner_at_edge` True; gamma 0, cisla beze zmeny, spatne jen `method_used` a priznak; 40 nahodnych spekter bez dalsiho pripadu |
| F2 | A | selhany bod (inf) roh nerozhodne | `np.gradient` udela z inf inf/NaN, `np.argmax` vrati prvni NaN: cista krivka roh 4, inf v bodu 7 -> roh 6 s krivosti NaN | jen pri NNLS maxiter (50 k) v jednom bode okna; simulovano (synth: lambda 6.6e-10 -> 4.9e-9, `lcurve`, 7x), na realnych datech nenastalo |
| F3 | A | ValueError, kdyz zadna lambda neda konecne skore; nesouhlas tvaru je ValueError | `argmin` nad samymi inf vrati (1e-11, inf); spatna delka b: ValueError z nnls se chyti jako inf, stejny vysledek; NaN v A totez | jen API; N = 3 az 6 frekvenci pres `calculate_drt` skore konecna |
| F4 | A | zjemneni nevrati skore nad hrubym minimem | minimum jen z jemne mrizky, ta pri sudem `n_search` hrube minimum neobsahuje | CLI: 1 ze 40 nahodnych spekter, pomer skore nejhur 1.00034 |
| F5 | A | konecne skore je vzorec; kde M neumozni presny trace K, inf | `np.linalg.solve` bez kontroly podminenosti: ext3 pri lambda 1e-22 m - tr K 23.6 misto 13.5 (reference 60 cislic), skore 3x nizsi; duplicitni neregularizovany sloupec cond 2e22, solve nevyhodi | jen API: od lambda 1e-12 souhlas na 6+ cislic na vsech systemech, volajici zacina na 1e-11 |
| F6 | A | ValueError pro spatny `lambda_range`, `n_search` < 3, `lcurve_decades` | zadna kontrola: obraceny rozsah, `n_search` 2, `lcurve_decades` 0/nan/inf vrati cisla, (1e-10, inf) vrati NaN | jen API; volajici predava pevne platne hodnoty |
| F7 | C | bez RuntimeWarning | `np.where` deli i tam, kde je jmenovatel 0 ("invalid value encountered in divide"), vysledek 0 spravne | CLI: cisty rezistor; dalsi varovani jen z API (zaporna lambda, inf, F2) |
| F8 | C | nikdy NaN; vyjimky jen ValueError a LinAlgError | A a b radu 1e160: A^T A pretece, skore NaN; `find_lcurve_corner` na prazdnem vstupu IndexError | jen API, extremni vstup |
| - | B | kontrakt mlci | lambda = 0 v `compute_gcv_score` a `compute_lcurve_point` vrati cisla, ackoli Problem definuje x jen pro lambda > 0 | jen API |

Overeno, ze plati (reviewer): b x c pro c 1e-6 az 1e8 na 7 systemech
(lambda, `method_used`, priznaky beze zmeny, skore / c^2 na 1e-6, rho
a eta posunute o log10 c na 6e-13, krivost na 2e-8); permutace radku;
determinismus; vstupy nezmenene; zadna varovani na beznych systemech;
(inf, inf, None) a inf pro selhane NNLS (NaN v A, inf v b, NaN/zaporna/
nekonecna lambda); trace na 6+ cislic od lambda 1e-12; krivost nulova pod
3 body a pri nulovem jmenovateli; roh uprostred bez vnitrniho bodu,
nejnizsi index pri shode; rozsahy jemneho hledani (uvnitr i na obou
koncich) a lambda v [lo/10, hi x 10]; okno L-krivky (stred orezany do
rozsahu, okno orezane), max, pravidlo `method_used`, vzorce
`corner_at_edge` a `corner_below_gcv`; skore pri vyhre GCV; rozpocet
NNLS 50 x sloupce.

Provedeno 2026-10-23 (minimalni varianta, rozhodnuti uzivatele; F4 bez
opravy; jedna zmena s testem v `tests/test_hybrid_lambda.py`):

- Kod: F1 - roh bez kladne krivosti neni roh: lambda z GCV,
  `lambda_lcurve` None, `method_used` 'gcv', priznaky False; F2 -
  krivost do dvou bodu od selhaneho bodu NaN, roh jen z konecnych
  krivosti; F3 - ValueError, kdyz zadna lambda hrubeho hledani neda
  konecne skore (pokryva i spatnou delku b a NaN v A); F7 - derivace
  a deleni v `np.errstate`. Cisla na beznych spektrech beze zmeny.
- Kontrakt: F4 - minimum jen z jemne mrizky, muze byt nepatrne nad hrubym
  (1 ze 40, 3.4e-4); F5 - podminenost se nekontroluje, presnost
  namerena od 1e-12; F6, F8 a lambda = 0 - zaruky volajiciho, nic se
  neslibuje mimo ne; roh bez konecne krivosti = prvni vnitrni bod s NaN.
- Test bez opravy pada.
- /code-review (8 nalezu): maska selhanych bodu `mode='same'` mela pro 3-4
  body delku 5 (IndexError, jen API s `n_search` 3-4) -> `mode='full'`
  orezany; jemna mrizka bez konecneho skore vratila (lambda, inf) -> vrati
  hrube minimum. Bez akce: poradi kontroly znamenka a debug logu, znamenko
  ve verejne `find_lcurve_corner` (kontrakt to necha volajicimu), dvoji
  hrube hledani v zalozni ceste, delka modulu 571 radku, b = 0 hlasi
  `lambda_at_edge` (tak i pred opravou), prah bez tolerance (podlaha
  1e-15 pres data nedosazitelna).

## Kolo 8: validation/noise.py (2026-10-23)

Kontrakt: docstring `eis_analysis/validation/noise.py`. Reviewer: Opus, jeden
beh (6 min, 86k tokenu, 18 volani nastroju), piskoviste `~/coldreview/noise/`
= cely balik; cist smel jen modul, zbytek jako cernou skrinku. Debrief:
precten jen `noise.py` (dva radky zdroje jinde jen z tracebacku); volajici
overeni pres Lin-KK, Z-HIT, `kk_grid`, DRT a `compare_fit_to_noise`
(20 nahodnych spekter, Monte Carlo sumoveho modelu); 9 skriptu v `review/`.
F1-F7, F9-F11 reprodukovany nezavisle (skript v scratchpadu sezeni) proti
modulu v repozitari, F8 uvahou (odhad sigma lezi pod sigma v polovine
behu). Neoveril tvrzeni o `NOISE_FLOOR_PERCENT` v `outliers` (cerna
skrinka, import uvnitr funkce); plati (`outliers.py:260`).

Autor kontraktu predem cekal F1, F5 (broadcast misto ValueError), F2 (NaN
misto inf pri Z_exp = Z_fit = 0), F3, F4 (n_points a zaporne chi2) a F9.
Necekal F6, F7 (prehnana robustnost medianu), F8, F10 a ze F2 rozbije
`kk_grid` a Lin-KK na samych nulach.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | A | ruzne tvary jsou ValueError, ne broadcast | `compute_pseudo_chisqr` tvary nekontroluje: (3,) proti (1,) 1.36, proti (3,1) 7.09 (vnejsi soucin), proti skalaru 0 3.0 | jen API; Lin-KK, Z-HIT, DRT predavaji stejne delky, `compare_fit_to_noise` je kontroluje |
| F2 | B | bod se Z_exp = 0 da inf; 0 pro shodna pole | Z_exp = Z_fit = 0 v bode: weight inf x 0 = NaN; shodna pole s nulou NaN; kontrakt si sam odporuje | volajici jen se Z = 0 vsude: `lin_kk_native` chi2 a sum NaN, `success` True bez varovani; `find_optimal_extend_decades` ValueError "min() iterable argument is empty"; jeden bod Z = 0 da inf (Z_fit ~1e-30) |
| F3 | A | ValueError, kdyz n_points neni kladne cele cislo | kontroluje jen `<= 0`: n = inf da 0.0 %, 2.5 a True cisla, NaN da nan | jen API; volajici predavaji `len(Z)` nebo pocet masky |
| F4 | B | zaporne chi2_ps je ValueError | `np.sqrt` vrati nan (RuntimeWarning), -inf take | jen API |
| F5 | A | ruzne delky residui jsou ValueError | `np.hypot` broadcastuje: delka 4 proti 1 da 1.20; 2-D vstup se zplosti | jen API; KKResult predava stejne dlouha residua |
| F6 | A | nan, kde je residuum NaN | `hypot(inf, nan)` = inf, median mensiny inf konecny: 1.20 a bod s NaN se tise vynecha | jen API; Lin-KK NaN/inf v Z odmita (`estimate_R_linear`) |
| F7 | B | mensina velkych residui posune median nejvys o sousedni poradovou statistiku | k velkych residui z 81 posune odhad o k poradovych statistik ciste sady (k = 1, 2, 5, 20: #42, 43, 46, 61 proti #41), o k/2 vuci medianu nedotcenych bodu; omezene do poloviny bodu | kazdy volajici s odlehlymi body |
| F8 | B | kazdy odhad je horni mez sumu | horni mez jen ve stredni hodnote: pri Z_fit = Z_true pod sigma v 51 % (stredni kvadrat) a 49 % (median) behu; Lin-KK 7 %, median 15-20 %, Z-HIT 1 %, DRT 20-30 %; stredni hodnoty volajicich 1.17-1.39 % pri 1 % | vsichni volajici |
| F9 | B | `noise_estimate_on` se vsemi body = `noise_estimate` | se Z = 0 v jednom bode `noise_estimate` inf, `noise_estimate_on` 69.8 % (Lin-KK normuje residua s podlahou 1e-15) | Lin-KK se Z = 0 v bode |
| F10 | C | invariance skalovanim pro \|Z\| 1e-150..1e150 | ctverce residui subnormalni: pri relativnim residuu 1e-10 a skale 1e-150 relativni zmena 1e-5 az 3.5e-6 | jen API; impedance 1e-6..1e13 |
| F11 | C | kontrakt mlci | seznam: `compute_pseudo_chisqr` AttributeError, median ho prijme; skalar v medianu TypeError | jen API |

Overeno, ze plati (reviewer): 0 pro prazdna pole; inf pro Z_exp = 0 s
nenulovym Z_fit; NaN a inf v kterekoli poli nikdy nedaji konecne cislo;
konjugace presne, komplexni skala 1e+-150 na 1e-15 pri realnych
residuich; `estimate_noise_percent` inf -> inf, NaN -> nan, n <= 0
ValueError; median NaN -> nan, prazdny -> nan, Rayleighova konstanta
(Monte Carlo 1.004 % pri 1 %, stredni kvadrat 0.999 %); volajici na 20
spektrech: Lin-KK a Z-HIT `noise_estimate` = vzorec presne,
`noise_estimate_on` se vsemi body = `noise_estimate` na 2e-16 pri Z != 0,
median = funkce nad ulozenymi residui, Z-HIT Z_fit nese merenou fazi na
1e-16; DRT `noise_sigma_used` = vzorec / 100 na 2e-18;
`compare_fit_to_noise` ValueError pri ruznych delkach a prazdnem vstupu,
nekonecny vstup hlasi jako neporovnany; `kk_grid` s jednim Z = 0 dostane
inf; `NOISE_DOMINATED_RATIO` 2.0.

Mimo modul (vedlejsi nalezy reviewera): `compare_fit_to_noise` vraci pri
nekonecnem vstupu `adequate` True; Lin-KK a Z-HIT neodmitaji Z = 0, na
samych nulach `lin_kk_native` `success` True se sumem NaN,
`find_optimal_extend_decades` nesrozumitelna ValueError (wrapper
`kramers_kronig_validation` ji hlasi jako `error`).

Provedeno 2026-10-23 (minimalni varianta, rozhodnuti uzivatele; bez zmeny
kodu modulu; test kontraktu v `tests/test_kramers_kronig.py`):

- F2 nejdriv opraven v kodu (inf i pri Z_fit = 0). /code-review (8 nalezu)
  ukazal regresi: na samych nulach `kramers_kronig_validation` misto
  `error` vratil `is_valid` True se sumem inf, `lin_kk_native` s
  auto-extend misto ValueError totez; jinde (jeden bod Z = 0, DRT) beze
  zmeny. Oprava vracena: Z = 0 je neplatny vstup pro vsechna relativni
  residua a patri odmitnout na vstupu Lin-KK a Z-HIT (samostatna oprava,
  ne v tomto kole). Ostatni nalezy review (RuntimeWarning, CHANGELOG,
  test) odpadly s ni.
- Kontrakt: F2 - Z_exp = 0 da inf, NaN kde je i Z_fit = 0, odmitnout
  Z = 0 je na volajicim; F1, F5, F11 - ndarray stejneho tvaru (median
  1-D) je zaruka volajiciho, tvary se nekontroluji; F3, F4 - kladne cele
  n_points a chi2_ps >= 0 zaruka volajiciho, ValueError jen pro
  n_points <= 0; F6 - (inf, NaN) se pocita jako velke residuum; F7 - k
  velkych residui posune median asi o k/2 poradovych statistik, omezene
  pro k < N/2; F8 - horni mez ve stredni hodnote, jednotlivy odhad kolem
  sigma; F9 - rovnost s `noise_estimate_on` jen pro Z != 0 (i v docstringu
  `KKResult.noise_estimate_on`); F10 - invariance, dokud |Z| i
  |Z_exp - Z_fit| lezi v 1e-150..1e150.
- Test pinuje inf (na nem stavi `kk_grid`) a NaN pri Z_fit = 0.

Navazujici oprava 2026-10-23 (rozhodnuti uzivatele, vlastni commit):
Z = 0 odmita sdileny `validate_eis_data` ValueError (Lin-KK vcetne
`find_optimal_extend_decades`, wrapper: `error`; fit Voigtovym retezcem;
`local_exponent`, kontrakt doplnen, ackoli Y tam je konecne),
`zhit_validation` vraci prazdny vysledek s duvodem v logu. Podlaha
max(|Z|, 1e-15) v residuich Lin-KK smazana: rovnost `noise_estimate_on`
s `noise_estimate` (F9) plati pro kazde |Z| > 0. Test fallbacku `kk_grid`
se Z = 0 smazan (nedosazitelny). /code-review (10 nalezu): vedle techto
opravena mezera; bez akce traceback v Z-HIT (CLI ho nevypisuje), 0.0
v prazdnem vysledku Z-HIT (stavajici cesta selhani), podteceni |Z| <
1e-154 (mimo rozsah kontraktu). Druhy /code-review (9 nalezu): kontrola
ve sdilenem validatoru odmita Z = 0 i tam, kde se |Z| nedeli (fit
Voigtovym retezcem s vahami `uniform`, `local_exponent`: driv 60 z 61
bodu platnych) - ponechano, rozhodnuti uzivatele (spektrum se Z = 0 je
poskozene; vychozi fit s vahami 1/|Z| ho driv tise pretahl vahou 1e15).
Z-HIT odmita i NaN a inf (driv `success` True se sumem NaN); upresnen
kontrakt (rozsah rovnosti `noise_estimate_on`, sigma v DRT jen vychozi,
v `local_exponent` konecne 1/(Z - R_inf - j omega L)). Zbyva: DRT Z = 0
neodmita (test vyznamnosti piku pak nedefinovany). Treti /code-review
(8 nalezu): Z-HIT propoustel NaN a zaporne frekvence -> vlastni kontrola
nahrazena `validate_eis_data`; ostatni bez akce (rozhodnuto vyse,
pretekajici |Z| u 1e308, podlahy vah pro |Z| < 1e-15, opakovana kontrola
O(N) v iteracich).

DRT 2026-10-23 (rozhodnuti uzivatele, vlastni commit): `calculate_drt`
vola `validate_eis_data` misto vlastni kontroly NaN/inf, Z = 0 dava
prazdny vysledek (driv pseudo chi^2 inf v testu vyznamnosti piku,
chyba rekonstrukce inf %). Kontrakty `noise.py` a `drt/significance.py`
doplneny. /code-review (7 nalezu): validator bere za nulu i |Z|^2 = 0
(podteceni pod ~1e-154, driv chyba rekonstrukce 1.8e201 %); `f_min` az po
validaci (prazdny vstup padal nezachycenou vyjimkou); zalomeni kontraktu.
Bez akce: dvoji kontrola frekvenci, prazdny `DRTResult` bez duvodu
(stavajici cesta selhani); presun validatoru z `fitting/voigt_chain/`
do `utils/` nabidnut jako samostatny refactor.
/code-review rozsahu b3f7fe6..3389852 (9 nalezu): Z-HIT validuje pred
razenim (delsi Z se tise orezal, kratsi IndexError), `calculate_drt`
vola validator pred `_validate_frequencies` (prazdny vstup dostane jeho
hlasku), prah podteceni v CHANGELOGu ~2e-162, komentar `kk_grid`,
mrtva ochrana `+ 1e-30` ve vahach R_s smazana. Bez akce: kontrola ve
validatoru a v loaderech (varianta A), traceback v Z-HIT (CLI ho
nevypisuje), preteceni pro |Z| > 1e154.
/code-review a050aa1 (6 nalezu): `+ 1e-30` nebyl mrtvy - pro |Z| mezi
~1.6e-162 a ~1.5e-154 je |Z|^2 subnormalni a 1/|Z|^2 pretece; validator
proto odmita |Z|^2 pod nejmensim normalnim floatem (|Z| < ~1.5e-154).
Z-HIT loguje traceback jen u jine nez ValueError; nedosazitelne kontroly
NaN a f <= 0 v `_validate_frequencies` smazany; komentar `kk_grid`
upresnen. Bez akce: varovani o malem rozsahu frekvenci z
`_validate_frequencies` se zahazuje (tak i driv).
/code-review 5886ae7 (9 nalezu), posledni kolo na toto tema (rozhodnuti
uzivatele): kontrakt `local_exponent` a docstring Z-HIT doplneny o
|Z| < ~1.5e-154, hlaska odvozena z prahu, Z-HIT bez tracebacku jen pro
odmitnuty vstup (ne pro kazdou ValueError), `_validate_frequencies` ma
zpet pojistku NaN a f <= 0. Znama mez (docstring `validate_eis_data`):
mimo |Z| ~1e-150..1e150 vypocty nejsou invariantni vuci meritku (soucty
1/|Z|^2, podlaha 1e-15 v `compute_weights`, complex64); oprava by byla
prepis vzorcu na relativni residua, samostatny ukol.

## Kolo 9: validation/kk_grid.py (2026-10-26)

Kontrakt: docstring `eis_analysis/validation/kk_grid.py`. Reviewer: Opus,
jeden beh (10 min, 129k tokenu, 47 volani nastroju), piskoviste
`~/coldreview/kk_grid/` = cely balik; cist smel jen modul, zbytek jako
cernou skrinku. Debrief: precten jen `kk_grid.py` (jeden grep pres balik
ukazal jednotlive radky jinych modulu, zdroj jinde jen z varovani);
volajici overeni pres `lin_kk_native`, `kramers_kronig_validation`
a volani KK grafu (120 nahodnych spekter); 20 skriptu v `review/`.
Pozor: z piskoviste `~/venvs/eis/bin/python` importuje editovatelnou
instalaci z repozitare, reviewer to obesel `review/_setup.py` (sys.path).
F1-F3 a F5-F10 reprodukovany nezavisle (skript v scratchpadu sezeni)
proti modulu v repozitari, F7 jen na `find_optimal_extend_decades`, F4
skriptem reviewera (balik v piskovisti shodny s repozitarem, `diff -r`).

Kontrakt puvodne u `_chi2_lower_M` sliboval plny beh CHI2_PLATEAU_RUN
pred `max_M` a jinak M = 3; mereni pred review (15 KK bodu, 0.1 %: plny
beh start 3, sum 12.8 %, neplatne; kod start 13, 0.8 %, platne) -> kontrakt
prepsan podle kodu (beh zkraceny u `max_M`), rozhodnuti uzivatele.

Autor kontraktu predem cekal F1 (prazdne `elements` IndexError), F3
(NaN kandidat nevypadne, `min()` nad prazdnym seznamem) a ze vetev None
v `lin_kk_native` je pri rozsahu od 0 nedosazitelna (potvrzeno: 120
nahodnych spekter, i rozsahy (0.5, 1) a (1, 2)). Necekal F2, F4-F10.
Cekany rozpor docstringu (min_chi2 = chi^2 viteze, ne minimum) reviewer
nehlasil; `n_evaluations` = 0 (None bez min_mu) nezkousel.

| Id | Zav. | Klauzule | Co kod dela | Dosah v knihovne |
|---|---|---|---|---|
| F1 | C | ValueError, kdyz `elements` nema 1 + len(tau) + [include_L] polozek, i prazdne | `elements[0]` pred kontrolou: prazdne `elements` IndexError | jen API; Lin-KK a graf predavaji vysledek fitu |
| F2 | C | totez | s `include_L` a jednou polozkou je `elements[1:-1]` prazdne a sedi na tau = []: [5] misto dvou polozek projde; cisla jako pri spravnem volani (slot L se nepouziva) | jen API |
| F3 | B | kandidat s nekonecnym chi^2 nebo mu nesoutezi; None, kdyz zadny nezbude | konecnost se nekontroluje, NaN < min_mu je False: NaN u prvniho kandidata -> `min()` vrati NaN, mnozina shody prazdna, ValueError "min() arg is an empty sequence"; NaN jinde se tise preskoci (zavisi na poradi); same inf -> (0.0, inf) misto None | Z x 1e200 (mimo znamou mez 1e-150..1e150): `kramers_kronig_validation` `error` s touto hlaskou a ~430 RuntimeWarning; NaN/inf ve fyzikalnim rozsahu jen simulovane |
| F4 | A | kontrakt mlci | okno hleda minimum jen dopredu: na chi^2 stoupajicim s M (fit 'imag' na vysokem M) je kazde M "na plateau". 40 bodu, 1 RC, 0.1 %: l(M) minimum -4.0 u M = 19, u M = 37 +1.1; start 33, `lin_kk_native` imag M 33, mu 0.000, sum 0.75 % (real 0.11 %, complex 0.10 %) | jen `lin_kk_native(fit_type='imag')` pres API; CLI a `kramers_kronig_validation` fituji 'real', tam start vzdy pod minimem chi^2; vlastni 1-RC spektrum 0 z 20 seminek |
| F5 | B | f x c necha extenzi, chi^2 a mu, tau, L, C / c; pro kazdou vahu | 'imag' s C: f x 1e3 extenze 0.0 -> 0.2, chi^2 1.8e-4 -> 7.7e-5, C x c -0.51 -> 94 (podminenost fitu v `estimate_R_linear`); 'real' s C plati do 1e3, rozpadne se az u 1e6 | jen API ('imag'); CLI 'real' v realnem rozsahu f invariantni |
| F6 | B | Z x c bez rozsahu | nad ~1e150 pretece |Z|^2 v `noise.py`: 1e152 extenze 0.4 -> 0.5, start M 20 -> 27, 1e153 ValueError (F3) | mimo znamou mez (docstring `validate_eis_data`); kontrakt ji nejmenuje |
| F7 | B | permutace a skala meni vysledek jen o zaokrouhleni | presna data (cisty R, R-L) maji chi^2 1e-31 az 1e-28: relativni shoda a log10 rozhoduji zaokrouhleni, extenze pres poradi bodu 0.0 / 0.1 / 0.2 / 0.6; `kramers_kronig_validation` cisty R podle reviewera M 6 / 16 / 6 | jen presna syntetika; sum ~1e-14 %, verdikt beze zmeny |
| F8 | C | zadne RuntimeWarning na tom, co volajici predavaji; f <= 0 je ValueError | mrizka tau se stavi pred validaci v `estimate_R_linear`: f = 0 da "divide by zero" a "invalid value" pred ValueError, f < 0 "invalid value in log10" | jen API: loadery (Gamry, CSV) body s f <= 0 zahodi; `kramers_kronig_validation` s f = 0 varovani a `error` |
| F9 | B | n_evaluations hodnot pres rozsah vcetne obou koncu, n_evaluations >= 1 | `np.linspace(lo, hi, 1)` = [lo]: n = 1 fituje jen dolni mez | jen API; Lin-KK predava 11 |
| F10 | C | kontrakt mlci | degenerovane frekvence projdou: vsechny f stejne (mrizka tau v jednom bode) vrati (0.0, 10.3), M = 40 > N = 20 vrati vysledek | jen API; `lin_kk_native` omezuje `max_M` na N - 2 |

Overeno, ze plati (reviewer): `reconstruct_impedance` (slot L se
nepouziva, linearita v (elements, L, 1/C) na 2e-16, shoda s modelem, tvar
0-D i 2-D, complex128, vstupy nezmenene); mrizka tau (konce, zaporna
extenze = 0); vraceny tuple reprodukuje sve chi^2; pravidlo shody
a nejmensi extenze; determinismus a nezmenene vstupy i pres
`kramers_kronig_validation`; permutace na 30 spektrech se sumem; Z x c
1e-6 az 1e6 pro vsechny vahy a typy fitu; `M_lower` = vysledek
`_chi2_lower_M`; kandidat s extenzi 0 v 120 behach vzdy ponechan (mu
presne shodne); volani KK grafu reprodukuje Z_fit bitove pro vsechny
kombinace L / C / typu fitu; ValueError pro neplatnou vahu, typ fitu
a data.

Provedeno 2026-10-26 (rozhodnuti uzivatele: oprava jen F3, F5-F10
kontraktem; F1, F2 a F4 bez akce; test v `tests/test_kramers_kronig.py`):

- Kod: F3 - kandidat bez konecneho chi^2 nebo mu vypadne; zadny
  nezbyde -> None, varovani `lin_kk_native` jmenuje oba duvody. Z x 1e200
  (nad znamou mezi) dava misto `error` sum NaN a `is_valid` False.
  Cisla na merenych datech beze zmeny.
- Kontrakt: F5 - invariance vuci f x c slibena pro 'real' a 'complex'
  (s C 'real' do 1e3), pro 'imag' ne; F6 - Z x c jen v mezi
  1e-150..1e150; F7 - "na zaokrouhleni" plati pro cisla, ne pro volby na
  chi^2 v urovni zaokrouhleni; F8 - f <= 0 da RuntimeWarning pred
  ValueError; F9 - n_evaluations = 1 jen dolni mez; F10 - stejne f
  a M > N se nekontroluji.
- Bez akce a kontrakt je nepopisuje: F4 (stoupajici chi^2 u 'imag').
  F1, F2 nakonec opraveny po /code-review (CR4, nize).
- Test bez opravy pada.
- /code-review (9 nalezu), provedeno: CR1 - nad mezi |Z| chi^2 nejdriv
  pretece na 0 (vaha 1/inf; Z x 1e154 a 1e155 na 2-RC: test projde se
  sumem 0 %, tak i pred opravou), NaN az od ~1e160: CHANGELOG a kontrakt
  opraveny; CR3 - None pri rozsahu od 0 i pri stop fitu bez konecneho
  chi^2 (kontrakt); CR4 - F1, F2 opraveny v kodu (pocet prvku pred
  indexaci, ValueError), test rozsiren; CR5, CR6 - docstring `min_mu`
  a navratovych hodnot (chi^2 viteze shody, ne minimum); CR8 - kontrola
  konecnosti mu zbytecna (`calc_mu` je NaN jen pri nekonecnem R_k, pak je
  nekonecne i chi^2), vracena na `mu < min_mu`. Bez akce: CR2 (pricina je
  vzorec v `noise.py`, prepis na relativni residua uz zapsan jako
  samostatny ukol u kola 8), CR7 (varovani spojuje dve priciny, rozliseni
  by menilo navratovy typ kvuli pripadu mimo rozsah), CR9 (pokryti testu,
  rozpocet testu).
- Druhy /code-review (7 nalezu), provedeno: prahy preteceni v CHANGELOGu
  a kontraktu platily jen pro jedno spektrum (2-RC se sumem 0.1 %: x 1e153
  sum 0.051 % misto 0.093 %, x 1e154 0 % platne, od x 1e155 NaN; presne
  totez spektrum 0 % az do x 1e156) -> popsano obecne: nad ~1.3e154 ve
  slozce Z pretece ctverec, chi^2 klesne (az 0 %, test projde) nebo je
  NaN podle dat; text "konecne mu" po CR8 opraven; zalomeni kontraktu.
  Bez akce: pricina v `noise.py` a slouceni pricin ve varovani (jako
  v prvnim review), NaN u pozdejsiho kandidata v testu (volitelne,
  nevybrano).
