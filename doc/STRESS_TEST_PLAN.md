# Plan: stress test s invarianty

Stav: schvaleny plan (2026-10-06), revize po kritickem cteni (2026-10-07).
Krok 2 implementovan a roztriden (2026-10-07): vysledky a zname limity
v `doc/STRESS_TEST.md`. Rozsireni 2026-10-08: rodina `anomalous` (Wa, Wat,
model ZrO2 vrstvy) a mapa n(f) (`local_exponent`), viz nize.
Krok 3 rozpracovan (2026-10-08): D, E, F, I, M v `tests/stress_consistency.py`,
prahy PROVIZORNI. Prvni plny beh roztriden: chyby testu opraveny (piky DRT
z `scipy_peaks`, serie C v Lin-KK podle pravidla CLI, Edc proti intervalu
do DC, Erec jen na uzavrenem VF konci), knihovna opravena (Z-HIT bez
unwrap, mez R_inf se sumem, n(f) s `R_inf_range`), konstantni sum omezen
na SNR >= 10. Plny beh 2026-10-09: prahy zkalibrovany (2x rezerva),
vysledky a zname limity 5-9 v `doc/STRESS_TEST.md`; `--ri-fit` zustava
volitelny (Ifit 0.80 / 0.64 pri 1 / 3 % sumu). Overovaci beh hotov.
Kontrolni bod 2026-10-09: z kroku 4 jen G (DE ze startu az 2 dekady od pravdy,
kazdy paty pripad; 300/300 proslo), H, J a L vynechany; krok 5 hotov
(`stress_baseline.json`, `--check`, `--update-baseline`, smoke test
s markerem `stress`). Zbyva: oprava `R_inf_range` (znamy limit 7).
Vychozi verze: eis_analysis v0.56.7.

## Kontext

Inspirace: driftlet (`doc/DRIFTLET_COMPARISON.md`). Jeho `npm run stress`
generuje tisice nahodnych realistickych pripadu a vysledky nesrovnava
s ocekavanou hodnotou, ale s **invarianty**, ktere musi platit vzdy.
Vetsina jeho numerickych oprav vznikla jako odpoved na selhani v nem.

U nas existuje jen ZScope benchmark (`tests/test_zscope_benchmark.py`,
4 pevne obvody x 3 urovne sumu, marker `slow`). Chybi pokryti mimo tyto
ctyri obvody, extremy (mOhm baterie, GOhm oxidy, blokujici elektrody, L)
a kontrola vlastnosti, ktere nezavisi na pravde. Uz ted je znamy pripad,
ktery by invariant chytil: absolutni meze `RESISTANCE_RANGE`
(`fitting/bounds.py:15`) rozbiji fit plosne normovanych dat (mereni B
v driftlet srovnani). Pozdeji (2026-10-07) uzavreno jako zamer: meze jsou
fyzikalni v zakladnich jednotkach, `doc/STRESS_TEST.md` limit 2.

**Cil:** nastroj, ktery na nahodnych, ale fyzikalne vernych spektrech
overi invarianty vsech hlavnich casti (DRT, Lin-KK, Z-HIT, R_inf, fit,
Voigt, oxid, CLI), vypise uspesnost po rodinach a invariantech a kazde
selhani umi zopakovat samostatne. Selhani se tridi na chybu (oprava +
regresni test) nebo znamy limit (zdokumentovany s cisly).

## Forma (rozhodnuto)

- **Skript** `tests/stress.py`: plny beh az ~1 hodina, max 4 procesy
  (`ProcessPoolExecutor(max_workers=4)`, viz pamet max-four-processes),
  tabulka uspesnosti, JSON souhrn, opakovani jednoho pripadu.
- **Maly pytest** `tests/test_stress_smoke.py` s markerem `stress`
  (vyrazen z vychoziho behu): prvnich N seedu kazde rodiny musi projit
  vsemi invarianty krome tech, ktere baseline pro presne tento pripad
  vede jako selhani.

Soubory (limit 500 radku na soubor):

| Soubor | Obsah |
|---|---|
| `tests/stress_cases.py` | generatory rodin, sit frekvenci, model sumu |
| `tests/stress_invariants.py` | jedna funkce na invariant, vraci `(ok, detail)` |
| `tests/stress.py` | runner: CLI, paralelizace, tabulka, JSON |
| `tests/stress_baseline.json` | mnozina selhani + miry statistickych invariantu |
| `tests/test_stress_smoke.py` | pytest nad malym vzorkem |
| `doc/STRESS_TEST.md` | vysledky, zname limity s cisly |

Pouzit: `parse_circuit_expression` (`cli/utils.py`, jako ZScope test)
pro stavbu pravdiveho obvodu z retezce, `Circuit.impedance`
(`fitting/circuit_builder.py`) pro pravdive Z. Model sumu prevzit ze
ZScope testu (proporcionalni, sigma = level * |Z| na Re i Im).

## Reprodukovatelnost

Bez nasledujicich tri bodu `--index` nezopakuje vysledek plneho behu
a invariant K nic neoveri.

1. **`rng` parametr v multistartu (zmena knihovny, rozhodnuto).**
   `fit_circuit_multistart` dnes perturbuje pres globalni `np.random`
   (`fitting/multistart.py:101,106,133,154`) a nema seed. Ve workeru
   `ProcessPoolExecutor` se globalni stav prenasi mezi pripady, takze
   vysledek zavisi na tom, co v procesu bezelo predtim. Pridat
   `rng` (None, int seed nebo `np.random.Generator`, cokoli bere
   `np.random.default_rng`; int jako `seed` u DE) do
   `fit_circuit_multistart` a do `perturb_from_covariance`,
   `perturb_from_stderr`, `perturb_log_uniform` (prvni dve jsou verejne
   exporty z `fitting/__init__.py`). `None` -> cerstva entropie:
   chovani zustava nahodne, ale bez globalniho stavu. CLI se nemeni.
   Regresni test: stejny seed -> shodny vysledek. CHANGELOG.
2. **Oddelene proudy nahody.** Kazdy ucel ma vlastni generator
   `np.random.default_rng([family_id, index, purpose])`, purpose:
   0 obvod, 1 sum, 2 start F2, 3 multistart, 4 seed DE. Pridani losovani
   do generatoru tak neposune sum ani starty a baseline zustane platny.
   `family_id` je pevne cislo v tabulce rodin (ne `hash(str)`, ten se
   kvuli `PYTHONHASHSEED` mezi procesy lisi); nova rodina dostane nove
   cislo, existujici se nikdy neprecisluji.
3. **Jedno vlakno BLAS.** `OMP_NUM_THREADS=1` nastavit na zacatku
   `tests/stress.py` pred importem numpy. 4 procesy x vicevlaknovy
   OpenBLAS pretizi jadra a vicevlaknove redukce nejsou bitove shodne.

## Generator pripadu

`python3 tests/stress.py --family oxide --index 37` zopakuje presne jeden
pripad. Kazdy pripad = (vyraz obvodu s pravdivymi hodnotami, f, Z, sum,
`n_arcs`, `min_sep`, `min_frac`).

Sit: f_max = 10^U(3, 7), f_min = 10^U(-3, 1), 5-15 bodu/dekadu.
Sum: {0, 0.1 %, 1 %, 3 %} proporcionalni; 20 % pripadu konstantni
(sigma = level * max|Z|).

| Rodina | Obvod (nahodne parametry) |
|---|---|
| `rc` | Rs - 1..4 x (R\|C), R log-uniformne 1e-3..1e9 Ohm v ramci pripadu do 3 dekad; 30 % s L (1e-8..1e-6 H, kabelova/parazitni) v serii |
| `cpe` | Rs - 1..3 x (R\|Q), n v [0.6, 0.98] |
| `diffusion` | Randles s W, Ws nebo Wo: Rs-(Q\|(R-W)) atd. |
| `blocking` | Rs-(R\|Q)-C nebo -Q (n v [0.85, 0.98]): kapacitni NF konec |
| `oxide` | Rs 1-100 - (R 1e6..1e9 \| Q 1e-11..1e-8, n 0.8-0.98) [+ druhy oblouk]; \|Z\| pres mnoho dekad (Zr oxidy) |
| `anomalous` | pulka Randles Rs-(Q\|(R-Wa nebo Wat)), pulka ZrO2 vrstva Rs-(G\|Wa\|Q\|C) nebo Rs-(Wat\|Q\|C); gamma 0.5-0.95 |

`anomalous` (id 6, pridana 2026-10-08) je samostatna rodina, ne rozsireni
`diffusion`: pridani Wa/Wat do losovani `diffusion` by zmenilo obvody
existujicich indexu a seedy citovane v `doc/STRESS_TEST.md` by prestaly
platit. gamma <= 0.95 drzi pravdu mimo gamma = 1 (presne Wo/Ws, horni mez).

Model ZrO2 vrstvy (fit realnych spekter M136, CHANGELOG `Wa`) se losuje
pres charakteristicke frekvence, ne pres hodnoty prvku, aby kazdy prvek
v okne neco urcoval: tri casy s rozestupy 0.3-3 dek uvnitr okna (jako
oblouky) davaji 1/omega prechodu Q -> C (Q omega^n = C omega), prechodu
Wa -> Q (\|Y_Wa\| = Q omega^n, s presnou admitanci prvku z knihovny, ne
s VF asymptotou: tau_W muze lezet jen 0.3 dek nad prechodem) a tau_W;
gamma < n jako ve vsech fitech, jinak by Q prevzal NF konec. Rs 1-100 Ohm, C 1e-9..1e-7 F (vrstva 0.1-10 um, plocha mm^2-cm^2),
n 0.6-0.9 (fit: 0.73-0.79), G 0.1-1x \|Y_Wa(f_min)\| (vodivost, ktera NF
konec ohyba, ale blokaci neschova). Vetev Wat (pulka ZrO2 pripadu) je bez
G: DC cestu ma Wat sam (model `(Wat|Q|C)` z CHANGELOGu).

Idealni n = 1 pokryvaji prvky C v `rc` a `blocking`.

### Poloha tau v okne

tau kazdeho oblouku lezi uvnitr okna [1/(2 pi f_max), 1/(2 pi f_min)]
s rezervou 0.5 dekady na kazde strane. Z tau a R se dopocita C = tau/R,
u CPE Q = tau^n / R. Rozestupy, ktere se do okna nevejdou, se losuji
znovu (nejuzsi okno 2 dek - 2 x 0.5 = 1 dek pojme 4 oblouky po 0.3 dek).
Oblouky mimo okno nejsou predmetem testu: shodily by E a F, aniz by slo
o chybu.

### Odstup pravdy od mezi

Kazdy pravdivy parametr lezi aspon 1 dekadu uvnitr `PARAMETER_BOUNDS`
(log-skala), u linearnich (n) mimo poslednich 1 % rozsahu. To je presne
prah `classify_bound_status` (`fitting/bounds.py:201`); bliz by fit
hlasil varovani z konstrukce. Pripad, ktery to porusi (napr. baterie
R = 1e-3, tau = 1e3 s -> C = 1e6 F), se losuje znovu; pocet zahozenych
se hlasi po rodinach. Zamerne k mezim vede jen kategorie "meze"
invariantu B.

### Pocet oblouku a jejich rozestup

Rodina urcuje typ prvku, ne jejich pocet: pocet oblouku se v ramci rodiny
losuje (`rc` 1-4, `cpe` 1-3, `oxide` 1-2; `diffusion` a `blocking` maji
pevnou strukturu). Jeden oblouk je jen nejjednodussi pripad. Vic oblouku
je nutnych, protoze vetsina realnych problemu je v jejich interakci:
rozliseni blizkych piku v DRT, preteceni R mezi sousednimi piky, lokalni
minima fitu s prohozenymi oblouky, neidentifikovatelnost v stderr.

- **Rozestup** sousednich oblouku Delta = |log10(tau_i / tau_j)| se losuje
  U(0.3, 3) dekady, takze tesne dvojice vznikaji zamerne (zadna pevna
  spodni mez). tau oblouku = R*C, u CPE (R*Q)^(1/n).
- Kazdy pripad si nese `n_arcs` a `min_sep` (nejmensi rozestup; u jednoho
  oblouku nedefinovan). `min_sep` se zarazuje do trid:
  **tesny** < 1 dek, **stredni** 1-2 dek, **volny** > 2 dek.

### Podil nejmensiho oblouku

Identifikovatelnost nezavisi jen na rozestupu: R se v ramci pripadu lisi
az o 3 dekady, a maly oblouk vedle velkeho (nebo s R pod urovni sumu)
nelze urcit ani pri volnem rozestupu. Druha osa trideni je proto
`min_frac` = min R_i / sum R_i: **slaby** < 0.05, **normalni** >= 0.05.
Prah 0.05 je vychozi odhad (oblouk pod 5 % R_pol se pri 3 % sumu ztraci
v sumu Re Z), upresnit z prvniho behu.

## Invarianty

Presne invarianty (rtol ~1e-6, plati vzdy):

- **A Robustnost.** Zadna vyjimka, zadne NaN/Inf ve vysledcich
  `calculate_drt`, `kramers_kronig_validation`, `zhit_validation`,
  `estimate_rinf`, `fit_equivalent_circuit`, `fit_voigt_chain_linear`,
  `analyze_oxide_layer`.
- **B Jednotky.** Z -> kZ, k v {1e-3, 1e3}: R-parametry a R_inf x k,
  C a Q / k, n a tau beze zmeny; DRT: stejne lambda, gamma x k, stejne
  tau piku (vzor: `test_weighted_drt_is_scale_invariant`,
  `tests/test_drt_weighting.py:130`); Lin-KK: stejne relativni rezidua
  a mu; Z-HIT: stejne relativni rezidua. U fitu se skaluje start i meze
  (PARAMETER_BOUNDS x k^power): absolutni meze ridi trust region scipy
  podle vzdalenosti i daleko od nich, takze s nimi B z principu neplati.
  Jejich vliv meri zvlast **Babs** (fit s absolutnimi mezemi, rozdil =
  status "meze", nikdy selhani). Fit v A/B/C/K je jeden LM, ne multistart
  (ten orezava perturbace na absolutni meze; patri do F2). Na DE se B neaplikuje:
  DE vzorkuje v absolutnich mezich, invariance tam z principu neplati.
- **C Poradi.** Obracene a zamichane poradi bodu dava stejny vysledek
  (nebo jasnou chybu, nikdy tichy rozdil). KK, Z-HIT a R_inf si data
  tridi samy, C je levna pojistka.
- **K Determinismus.** Dva behy se stejnym vstupem a stejnymi seedy
  (multistart purpose 3, DE purpose 4) jsou bitove shodne.

Konzistence (tolerance kalibrovat z prvniho behu s ~2x rezervou, jako
ZScope):

- **D Kramers-Kronig.** Kazde generovane spektrum je KK-konzistentni
  z konstrukce, takze Lin-KK musi projit (rezidua <= c * sum + podlaha);
  `blocking` s `include_C=True`. Z-HIT totez s vlastni podlahou
  (~1 % aproximacni chyba, viz DRT_PEAK_SIGNIFICANCE_PLAN).
- **E DRT.** gamma >= 0 (konecnost hlida A); u `rc`
  chyba rekonstrukce <= c * sum + podlaha; u uzavrenych spekter (faze na
  f_min > -5 deg) R_inf + R_pol ~ Re Z(f_min); u `rc` bez sumu, tridy
  `min_frac` normalni, s tau >= 1 dek od sebe ma kazde pravdive tau pik
  do 0.15 dek. Kontrola piku zavisi na tom, ktere piky DRT hlasi, a to
  zmeni DRT_PEAK_SIGNIFICANCE_PLAN: pokud ten pujde driv, kalibrovat az
  po nem, jinak pocitat s druhou kalibraci.
  (Puvodni "kazdy pik uvnitr okna nebo oznaceny" vypusten 2026-10-08:
  priznaky se pocitaji prave ze vzdalenosti od okna,
  `_flag_boundary_peaks`, takze nemuze selhat.
  Puvodni "soucet `R_estimate` <= R_pol" vypusten: u GMM je
  `R_estimate = weight * R_pol` a vahy davaji soucet 1, `drt/peaks.py:300`,
  takze nemuze selhat.)
- **F Fit (LM).** cost = sum |(Z - Z_model) * w|^2 s
  `w = compute_weights(Z, weighting)`, stejnou funkci vah a stejnym
  `weighting` jako fit (`fitting/circuit.py:328`), ne vlastni vahy.
  - F1 pevny bod: bez sumu, start v pravde -> pravda (rtol 1e-6).
  - F2 ne hur nez pravda: start = pravda x U(0.3, 3),
    `fit_circuit_multistart` -> cost(fit) <= cost(pravda) * (1 + 1e-6)
    + podlaha. Podlaha je nutna pro nulovy sum, kde cost(pravda) ~ 0;
    hodnotu kalibrovat z prvniho behu. Selhani = lokalni minimum; hlasi
    se mira, ne automaticky chyba. Pocita se zvlast po tridach rozestupu
    a `min_frac`: ve tride **tesny** a **slaby** je fit z principu spatne
    urceny, takze tam selhani znaci neidentifikovatelnost, ne chybu
    programu, a nesmi zaplavit souhrn ostatnich trid.
  - F3 pokryti intervalu (agregat pres pripady s proporcionalnim sumem
    a `is_well_conditioned`): podil parametru, jejichz pravda lezi
    v `FitResult.params_ci_95`, ma byt ~0.95 (pass 0.85-0.99). Testuje
    se interval, ktery knihovna uzivateli ukazuje (`t_critical`, pro
    kladne parametry v log-skale, `log_scale_ci_mask`), ne vlastni
    2 * stderr. Konstantni sum je vyrazen: vaha `modulus` mu neodpovida
    a pokryti 0.95 tam neplati.
  - F4 parametr na mezi (rtol 1e-6) => varovani
    v `diagnostics.bounds_warnings`. Jen tento smer: varovani prichazi
    uz do 1 dekady od meze, opacny smer by jen kopiroval
    `classify_bound_status`.
- **G DE** (kazdy paty pripad, 20 %): `fit_circuit_diffevo` se seedem (purpose 4)
  ze startu az 2 dekady od pravdy (purpose 6, `de_start`),
  cost <= cost(pravda) * (1 + 1e-3) + podlaha. Horsi konec je selhani.
- **H Voigt.** `fit_voigt_chain_linear`: Rs + sum R ~ Re Z(f_min)
  u uzavrenych spekter.
- **I R_inf.** Jen pripady bez L. 0 <= R_inf <= Re Z(f_max) * (1 + tol);
  kdyz je VF konec uzavreny (faze na f_max > -5 deg, obdoba NF konce
  v E), |R_inf - Rs| / Rs mala.
- **J Oxid.** (Z / a, plocha a) pro a v {1, 2} cm^2 -> stejna tloustka.
- **M Mapa n(f).** `local_exponent` dostava A, B, C, K jako ostatni
  analyzy (n, nejistota a maska `valid` se pri Z -> kZ nemeni), s
  odectenym pravdivym Rs a L (x k): testuje se mapa, ne `estimate_rinf`.
  Fitovane R_inf se pri B a C hybe v ramci sve stderr a n u f_max, vzor
  NaN i prah `valid` by ho nasledovaly (code review 2026-10-08: polovina
  pripadu falesne selhala). Konzistence pres cestu CLI (R_inf a L
  z `estimate_rinf`): u pripadu se sumem musi body oznacene `valid`
  souhlasit s mapou bezsumoveho Z s pravdivym Rs a L,
  |n - n_ref| <= 2 * `n_uncertainty`. Hlasi se mira (jako F3), ne
  selhani po pripadech; overuje prahy 0.02 a 0.1, kalibrovane zatim jen
  na ctyrech spektrech M136.
- **L CLI end-to-end** (podvzorek): `eis.py case.csv --no-show` skonci
  s kodem 0 a bez tracebacku.

## Runner

    python3 tests/stress.py                     # plny beh, vsechny rodiny
    python3 tests/stress.py --n 20 --family rc  # mensi beh
    python3 tests/stress.py --family oxide --index 37 -v   # jeden pripad
    python3 tests/stress.py --update-baseline   # zapsat stress_baseline.json
    python3 tests/stress.py --check             # selze na novem selhani

Vystup: tabulka invariant x radek `rodina/n_arcs` (napr. `cpe/1`,
`cpe/2`, `cpe/3`), aby se selhani u tri prekryvajicich se ZARC
nerozpustilo v prumeru s jednoduchymi pripady. Invarianty citlive na
identifikovatelnost (E rozliseni piku, F2, F3) maji navic sloupce po
tridach **tesny / stredni / volny** a **slaby / normalni**. Dale seznam
selhanych `(rodina, index, invariant, detail)`, pocet zahozenych
pripadu (odstup od mezi), celkovy cas. Selhani jednoho invariantu
nezastavi ostatni.

### Baseline

`stress_baseline.json` neuklada pocty, ale:

- **mnozinu selhani** `(rodina, index, invariant)` pro per-case
  invarianty. `--check` selze na kazdem selhani, ktere v baseline neni
  (pocty by neodhalily vymenu: jeden pripad opraven, jiny rozbit). Selhani
  z baseline, ktere zmizelo, se jen vypise (kandidat na `--update-baseline`).
- **miry** statistickych invariantu (F2 po tridach, F3 pokryti) s prahem;
  `--check` selze, kdyz mira klesne pod prah.

Ze stejne mnoziny cte smoke test, ktere invarianty ma u prvnich N seedu
preskocit.

### Rozpocet casu

Zmereno na RK3588, 1 vlakno, 3 oblouky, 71 bodu, 1 % sum: DRT
(auto-lambda) 0.46 s, Lin-KK 0.11 s, multistart 0.14 s, Z-HIT a R_inf
zanedbatelne, DE 21 s. Pri ~6 behach na pripad (original, 2x B, 2x C, K)
~4-5 s/pripad -> 1250 pripadu (250 na rodinu) ~25 min na 4 procesech;
DE u 10 % pripadu dalsich ~11 min. DE dominuje, proto B, C a K na DE
nebezi. Presne pocty doladit v kroku 4.

## Postup

1. (hotovo) Plan ulozen v `doc/STRESS_TEST_PLAN.md`.
2. `rng` parametr v multistartu (viz Reprodukovatelnost 1) s regresnim
   testem a CHANGELOG. Pak generator + runner + invarianty A, B, C, K
   na DRT, Lin-KK, Z-HIT, R_inf, LM fit. Zmerit cas, spustit, roztridit
   selhani.
3. Rodina `anomalous` a analyza n(f) s A, B, C, K; plny beh a trideni.
   Pak invarianty D, E, F, I, M. Kalibrovat tolerance a prah `min_frac`
   z prvniho behu.
   **Kontrolni bod:** zastavit a s uzivatelem podle realnych selhani
   rozhodnout, zda kroky 4 a 5 maji smysl v plnem rozsahu.
4. G, H, J, L. Doladit pocty na rodinu na rozpocet ~1 h.
5. `stress_baseline.json`, `test_stress_smoke.py`, marker `stress`
   v `pyproject.toml` (`markers`, `addopts = '-m "not slow and not stress"'`),
   sekce Testing v `CLAUDE.md`, `doc/STRESS_TEST.md` se znamymi limity.

Trideni selhani (kazde zvlast, s uzivatelem):
- **chyba** -> oprava v prislusnem modulu, regresni unit test s presnym
  spektrem z daneho seedu v beznem `tests/`, CHANGELOG, commit
  `fix(module): ...`;
- **znamy limit** -> zapis do `doc/STRESS_TEST.md` s cislem a seedem
  (driftlet styl: "impedance rozlisena na ~1e-4 |Z|", ne "muze byt
  nepresne").

Prvni ocekavany kandidat: invariant B na absolutnich mezich
`PARAMETER_BOUNDS`. Rozhodnuti (meze relativni k datum vs. zdokumentovany
limit) az s cisly z behu.

Commity, bump verze a tagy jen na pokyn (pamet review-before-commit).

## Overeni

- `python3 tests/stress.py --n 3` dobehne, vytiskne tabulku.
- `--family X --index N` reprodukuje bitove stejny vysledek jako v plnem
  behu (i kdyz v plnem behu bezel v procesu po jinych pripadech).
- Umele zanesena chyba (napr. docasne `R_inf * 1.1` v `estimate_rinf`)
  shodi invariant I; po vraceni zase projde.
- `--check` po umele vymene (jeden pripad opraven, jiny rozbit) selze.
- `python3 -m pytest tests/` beze zmeny poctu a casu (stress vyrazen);
  `python3 -m pytest tests/ -m stress` projde.
- `ruff check eis_analysis/ tests/`, `mypy eis_analysis/` ciste.
- Plny beh zmeren a cas zapsan do `doc/STRESS_TEST.md`.
