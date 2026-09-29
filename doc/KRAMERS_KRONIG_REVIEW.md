# Kriticka analyza modulu kramers_kronig

Datum: 2026-09-28, verze v0.44.1

**Stav:** body 1.2 a 1.3 jsou opraveny ve v0.45.0 (commit 4df6cb6), bod 1.1
ve v0.46.0 (commit 77e6e35). Z bodu 1.4 je vyresena
duplicita cisla 5 (konstanta `KK_RESIDUAL_THRESHOLD`); prumer jako kriterium
zustava.

Rozsah: `eis_analysis/validation/kramers_kronig.py` (590 radku) a funkce,
ktere vola: `fitting/voigt_chain/mu_optimization.py` (`find_optimal_M_mu`,
`calc_mu`), `fitting/voigt_chain/tau_grid.py` (`generate_tau_grid_fixed_M`),
`fitting/voigt_chain/fitting.py` (`estimate_R_linear`). Dale
`cli/handlers/validation.py` a `tests/test_kramers_kronig.py`.

Nalezy oznacene "overeno" byly reprodukovany skriptem na syntetickych
Voigtovych spektrech (presne KK-kompatibilnich, bez sumu) a na
`example/EISPOT-test1.DTA`. Referenci byla `impedance.py` 1.7.1
(`impedance.validation.linKK`).

## Celkovy dojem

Implementace je vernou kopii referencni Lin-KK: pro `fit_type='real'` i
`'complex'` dava stejne M, mu i rezidua jako `impedance.py` (overeno).
Hlavni problemy proto nejsou chyby prepisu, ale vedecke: s vychozim
nastavenim ukazuje test na dokonale KK datech rezidua v rade procent a
optimalizace `extend_decades` obchazi pojistku proti preuceni.

## 1. Chyby s dopadem na vysledek

### 1.1 Mu kriterium s prahem 0.85 zastavuje prilis brzy (overeno)

mu(M) neni monotonni. Na hrubem tau-gridu pinv aproximuje tau lezici mezi
body gridu stridave kladnymi a zapornymi R, takze mu kratce spadne pod prah
a hledani skonci.

| Spektrum (presne, bez sumu) | Prah | M | prumer \|res\| re/im | max \|res\| |
|---|---|---|---|---|
| 2-RC z `test_kramers_kronig.py` (10^-1 .. 10^5 Hz, 60 bodu) | 0.85 (vychozi CLI) | 7 | 3.8 / 4.0 % | 14 % |
| totez | 0.7 | 17 | 0.14 / 0.15 % | 0.6 % |
| 2-RC, tau = 1 ms a 0.5 s, 10^-2 .. 10^5 Hz, 71 bodu | 0.85 | 6 | 1.8 / 2.2 % | -9 % (imag, 0.01 Hz) |
| totez, `fit_type='complex'` | 0.85 | 7 | 11.0 / 10.9 % | - |

Prubeh mu pro posledni spektrum (real fit): M = 3..5 -> mu ~ 1.0, M = 6 ->
0.81 (stop), M = 7..9 -> 0.52 .. 0.65, M = 10 -> 0.86.

Testy se tomu vyhybaji pouzitim `mu_threshold=0.7`
(`tests/test_kramers_kronig.py:183`, `:210`); CLI pouziva 0.85.

Jde o znamou slabinu puvodni metody (Schonleber 2014), kterou resi pyimpspec
(Yrjana & Bobacka 2024) - na ten modul uz odkazuje u odhadu sumu.

Navrh minimalni opravy: nezastavovat na prvnim poklesu mu. Projit M az do
`max_M` a vzit prvni M, kde mu < prah a zaroven se pseudo chi^2 uz vyrazne
nezlepsuje. Pred zmenou vychoziho chovani zmerit na realnych datech.

**Stav: opraveno ve v0.46.0.** Hledani mu zacina na `M_lower`: prvni M,
jehoz log10(chi^2) je do 0.3 dekady od minima v nasledujicich 8 M
(`CHI2_PLATEAU_DECADES`, `CHI2_PLATEAU_WINDOW`). M je navic omezeno na
N - 2; pod 5 body Lin-KK skonci chybou. Mereni pri oprave:

- Realne soubory beze zmeny (EISPOT: M=19, ext 0.6; real_gamry: M=22).
- Presne 2-RC: 0.15 % misto 3.8 %; ZARC s 0.3 % sumem: 0.30 % misto 19 %;
  2-RC s 1 % sumem: 0.8-1.1 % misto 2.4 %.
- Jediny vadny bod (spike 2-20 %): velikost rezidua beze zmeny.
- Okno misto minima pres vsechna M: na datech s driftem chi^2 pomalu klesa,
  jak zaporna R drift pohlcuji, a globalni minimum posunulo start na
  M ~ 45-49 ve 2 ze 6 seedu (10% drift pak jen 2.4 % misto 9 %).

Dusledek pro drift (R_ct roste o 0-20 % behem mereni, sum 0.2 %): drive
"detekovan" 20% prumernym reziduem, ktere ale davalo i validni spektrum.
Nyni se drift projevi lokalne a v odhadu sumu (rozsah pres 6 seedu):

| drift | prumer \|res\| | max \|res\| | odhad sumu |
|---|---|---|---|
| 0 % | 0.18-0.23 % | 0.5-1.5 % | 0.19-0.26 % |
| 2 % | 0.32-0.41 % | 1.6-3.1 % | 0.36-0.50 % |
| 10 % | 1.04-1.14 % | 7.6-9.1 % | 1.52-1.66 % |
| 20 % | 1.91-2.00 % | 14.0-15.4 % | 2.88-3.01 % |

Prumer zustava pod 5 %, takze `is_valid` takova data prohlasi za validni.
Drift od ~10 % zachyti per-point report (`--max-residual` 5 %), 2% drift
jen zvyseny odhad sumu. To je nalez 1.4.

### 1.2 `auto_extend_decades` obchazi mu (overeno)

`lin_kk_native` vybere M na gridu s `extend_decades=0`. Potom
`find_optimal_extend_decades` fituje pri tomtez M na jinem gridu a vybira
pouze podle pseudo chi^2. Mu na vyslednem gridu se neprepocita.

Vynucene rozsireni na spektru z 1.1 (radek 3):

| ext [dekady] | hlasene mu | skutecne mu vysledneho modelu |
|---|---|---|
| 0.3 | 0.808 | -0.23 |
| 0.6 | 0.808 | -4.8 |
| 1.0 | 0.808 | -43 |

Skutecne mu znamena silne oscilujici zaporna R - presne to preuceni, pred
kterym ma mu chranit. Na `EISPOT-test1.DTA` je to naopak: hlasene mu 0.846,
skutecne 1.000.

Upresneni (overeno pri oprave): na 5 spektrech (EISPOT, dve presna 2-RC,
LF chvost tau = 200 s se sumem i bez) vybral chi^2 pokazde kandidata s
mu >= stop mu - preuceni kandidati fituji hure. Slo tedy o chybejici
zaruku, ne o pozorovany spatny vysledek; projevila se jen pri vynucenem
rozsahu rozsireni.

**Stav: opraveno ve v0.45.0.** `find_optimal_extend_decades` ma parametr
`min_mu` a zahodi kandidaty s mu pod nim; `lin_kk_native` preda stop mu a
kdyz nezbude zadny kandidat, ponecha grid bez rozsireni s varovanim. Stop
mu je tak dolni mez mu vraceneho modelu. Vysledky na overenych spektrech
se nezmenily.

### 1.3 Varovani z hledani M se zahazuji (overeno)

`lin_kk_native` (`kramers_kronig.py:424-425`) prevezme z `MuOptimization`
M, mu, tau, prvky, L a C, ale ne `warnings` ani `reached_max_M`.
`LinKKResult` pro ne nema pole. Varovani "Reached max_M ... model may still
be overfit" se k uzivateli nikdy nedostane.

Priklad: pri N = 3 bodech bezi smycka az do M = 50, protoze `max_M` neni
omezeno poctem bodu. Uloha je silne podurcena a uzivatel se to nedozvi.

**Stav: opraveno ve v0.45.0.** Varovani jdou do `LinKKResult.warnings`,
`KKResult.warnings` a CLI je vypisuje. Hlaska "Data may contain artifacts"
z `KKResult.warnings` zmizela, protoze opakovala `is_valid`. Omezeni M
poctem bodu reseno neni.

### 1.4 `is_valid` (prumer |res| < 5 %) je prilis volne a nekonzistentni

- Prumer schova lokalni poruseni: 10 z 70 bodu s rezidui 20 % da prumer
  ~2.9 % a data "projdou".
- Graf kresli +-5 % jako bodovou mez, `is_valid` ale testuje prumer.
- Cislo 5 neni nikde zduvodneno (pravidlo projektu o magickych cislech) a
  opakuje se na ctyrech mistech: `KKResult.is_valid`, `LinKKResult.is_valid`,
  CLI handler a graf.

### 1.5 Rozsireni gridu jde jen k nizkym frekvencim

Docstring `kramers_kronig_validation` (`kramers_kronig.py:497-500`) slibuje
reseni "capacitive/inductive tails". `generate_tau_grid_fixed_M` ale
rozsiruje pouze smerem k nizkym frekvencim. Vysokofrekvencni induktivni
chvost pokryva jen clen L. Docstring je treba opravit.

## 2. Architektura a robustnost

### 2.1 Vizualizace v jadru knihovny (overeno)

Odporuje pravidlu CLAUDE.md "Visualization separated from algorithms".

- `kramers_kronig_validation` vytvori figure pri kazdem volani; po trech
  volanich zustanou tri otevrene figure (unik pameti pri davkovem
  zpracovani).
- Uz samotny `import eis_analysis.validation` natahne `matplotlib.pyplot`.
- `KKResult` neobsahuje `elements` ani `tau`, takze si volajici graf mimo
  modul nevykresli.

### 2.2 `except Exception` polyka programatorske chyby (overeno)

`kramers_kronig.py:531`. Predani listu misto ndarray vrati
`KKResult(error="'<=' not supported between instances of 'list' and 'int'")`,
zaloguje se jen na urovni debug a CLI to poda jako "KK validation failed".
Stacilo by chytat `ValueError` a `np.linalg.LinAlgError`.

### 2.3 Duplicita `KKResult` / `LinKKResult`

Deset stejnych poli a tri stejne property (`mean_residual_real`,
`mean_residual_imag`, `is_valid`). Cistsi by bylo
`KKResult(fit: Optional[LinKKResult], warnings, error)`. Zaroven by to
vyresilo chybejici tau a prvky (2.1) i zahozena varovani (1.3).

## 3. Drobnosti

- **Mrtvy kod:** `if lkk.Z_fit is None` (`kramers_kronig.py:535`) nemuze
  nastat, `Z_fit` je v `LinKKResult` povinne pole.
- **Tuple misto dataclassy:** `find_optimal_extend_decades` vraci 6-tuple,
  projekt predepisuje `*Result` dataclassy.
- **Nezdokumentovane konstanty:** tolerance `0.001` (`:344`),
  `n_evaluations=11` (`:434`), `5000` v `estimate_noise_percent`
  (= 100^2 / 2 za predpokladu stejneho relativniho sumu v obou slozkach;
  patri do docstringu spolu s tim, ze jde o horni odhad).
- **`estimate_noise_percent(chi2, 0)`** skonci `ZeroDivisionError`.
- **Nekonzistentni ochrana |Z| = 0:** rezidua maji floor 1e-15,
  `compute_pseudo_chisqr` ne.
- **Krehke API `reconstruct_impedance`:**
  - L je ulozene uvnitr `elements`, C mimo ne;
  - `zip(R_i, tau)` pri nesouhlasnych delkach potichu orizne data;
  - `include_L=True` u pole bez L potichu zahodi posledni R.
- `kramers_kronig_validation` natvrdo pouziva `fit_type='real'`,
  `weighting='modulus'`, `include_L=True`; neni to v parametrech ani v
  docstringu.

## 4. Co je v poradku

- Vernost referenci (`impedance.py`) je overena.
- Vzorec pro mu, normalizace rezidui pres |Z| a pseudo chi^2 podle Boukampa
  jsou spravne.
- Seriova kapacita (`include_C`, ekvivalent `add_cap`) funguje.
- Hranice knihovna / CLI je, az na graf, dodrzena.

## 5. Doporucene poradi

1. **1.2 a 1.3** - male, jasne zmeny: prepocitat mu po rozsireni gridu (a
   odmitnout rozsireni, ktere mu shodi pod prah) a propagovat varovani.
2. **1.1** - zmena vyberu M; nutne zmerit na realnych datech pred zmenou
   vychoziho chovani.
3. **1.4 a 1.5** - kriterium platnosti a docstring.
4. **2.1-2.3** - refaktor: graf do `visualization/`, jeden result typ,
   uzsi `except`.
5. Sekce 3 prubezne.
