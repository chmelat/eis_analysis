# Audit modulu R_inf estimace (`--ri-fit`)

**Datum:** 2026-09-25 | **Verze:** 0.37.0 | **Typ:** kritické code-level review s numerickým ověřením (Claude Code)

**Rozsah:** `eis_analysis/rinf_estimation/` (`rlk_fit.py`, `data_selection.py`),
CLI handler `cli/handlers/rinf.py`, napojení do DRT (`drt/estimation.py:_estimate_r_inf`,
`cli/handlers/drt.py`), dokumentace `doc/RINF_ESTIMATION.md` a README (`--ri-fit`),
testy `tests/test_rinf_estimation.py` a `tests/test_cli_integration.py`.

Audit nemění produkční kód. Všechna čísla níže jsou naměřená. Pro každý nález
existoval regresní test (`xfail(strict=True)`), který chybu reprodukoval. Testy
jsou v commitu `6a7f8f1` a revert `c6eb751` je zase odstranil.
Obnovit je lze příkazem `git show 6a7f8f1:tests/test_rinf_estimation.py`.

---

## Verdikt

`--ri-fit` se prezentuje jako "R-L-K() fit", ale samotný fit běží jen v jedné
ze tří větví a v té je rozbitý. Výsledek určuje znaménko Im(Z) v nejvyšší
dekádě:

| Větev | Kdy | Hodnocení |
|---|---|---|
| průchod nulou | Im(Z) v dekádě mění znaménko | **systematicky zkreslená**, horší než výchozí medián |
| polynom Re(Im) | celá dekáda kapacitní | jediná užitečná část modulu, ale s magickou konstantou a nestabilní na šum |
| R-L-K lineárně | celá dekáda indukční | **rozbitá**: pevné tau, vrací i R_inf = 0 při R^2 = 0.994 |

Hlavní závěr: v indukčním případě, pro který je `--ri-fit` v dokumentaci
doporučen, vrací horší odhad než výchozí medián (+10 % proti +0.6 %,
+89 % proti +38 %). Chybu přitom nijak nesignalizuje: `fit_success=True`,
žádné varování, u rozbité větve navíc R^2 blízko 1.

Model R-L-(R|Q) nafitovaný přes horní dvě dekády stávajícím
`fit_equivalent_circuit` měl ve všech testovaných případech chybu pod 1 %
(šum 1 %: bias <= 0.2 %, směrodatná odchylka <= 1.6 %). Směrodatná chyba R_inf
z Jacobiánu přitom spolehlivě označila případy, kdy R_inf z dat určit nejde.
Z toho vychází hlavní návrh N1+N2.

---

## Metodika

**Syntetická spektra** se známým Rs, 10 bodů/dekádu. Kde není uvedeno jinak,
Rs = 10 Ohm, R = 100 Ohm a spektrum sahá do 1 MHz (případy B a D do 100 kHz):

| Případ | Model | fáze(f_max) |
|---|---|---|
| flat | Rs + RC, oblouk při 16 Hz, bez L | -0.1 deg |
| pure R-L | Rs + jwL, L = 0.1 uH | +3.6 deg |
| A1 | L = 0.1 uH, oblouk při 16 MHz (nad f_max) | -2.9 deg |
| A2 | L = 1 uH, oblouk při 1.6 kHz | +31.5 deg |
| A3 | L = 10 uH, oblouk při 160 kHz | +75.2 deg |
| C | L = 10 uH, RC tau = 10 us (16 kHz) | +80.7 deg |
| C2 | jako C, ZARC n = 0.7 | +77.7 deg |
| B n | kapacitní ZARC, n = 1 / 0.8 / 0.7, tau = 100 us | -9 / -17 / -20 deg |
| B2 | Rs + Warburg 30/sqrt(jw) | -0.2 deg |
| D | Rs = 0.05 Ohm, ZARC n = 0.6, tau = 1 ms | -52 deg |

Šum: multiplikativní komplexní gaussovský 1 %, 30-50 seedů (`default_rng(s)`).

**Reálná a vzorová data:** `example/EISPOT-test1.DTA`, `example/real_gamry_example.DTA`
(skutečné Rs neznámé), `example/example_eis_data.csv` (syntetický,
Rs = 10 Ohm, R0 = 1e5, CPE n = 0.6, šum 1 %).

**Srovnávané odhady:** výchozí HF medián, současný `--ri-fit` a čtyři
kandidáti na náhradu (viz "Srovnání metod").

---

## Nálezy

### R1 - Větev "průchod nulou" je strukturálně zkreslená (Vysoká)

`rlk_fit.py:176-218`. Hodnota Re(Z) v bodě Im = 0 není Rs. Obsahuje ještě
zbytek oblouku, Rk/(1 + w^2 tau^2), a indukčnost jen posune průchod k frekvenci,
kde se oblouk ještě neuzavřel. Model R-L-K má právě tuto chybu odstranit, ale
větev s průchodem nulou se vyhodnotí dřív a fit úplně přeskočí.

| Případ | medián (default) | `--ri-fit` |
|---|---|---|
| A2 (oddělený oblouk) | +0.0 % | +0.1 % |
| C | +0.6 % | **+10.0 %** |
| C2 | +37.5 % | **+89.0 %** |
| A3 | +59.8 % | **+101.7 %** |

Chyba je systematická, ne šumová: při šumu 1 % má C bias +10.1 % a směrodatnou
odchylku 0.7 %. `fit_success=True`, `warnings=[]`. Dokumentace
(`RINF_ESTIMATION.md`, "This is the strongest case: the real-axis intercept is
interpolated, not extrapolated") i README tuto větev popisují jako nejspolehlivější.

### R2 - Větev R-L-K s pevným tau vrací nesmysly s vysokým R^2 (Vysoká)

`rlk_fit.py:289-295`. Tau se odhaduje z maxima |Im(Z)|. V čistě indukční
dekádě ale |Im| roste s frekvencí, maximum leží vždy na okraji a tau se
nastaví na geometrický střed dekády. Prvek K pak neumí popsat ocas oblouku
a NNLS kompenzuje chybu na úkor Rs.

| Spektrum (Rs = 10 Ohm, dekáda čistě indukční) | R_inf | R^2 |
|---|---|---|
| L = 100 uH, RC tau = 10 us | **0.000** (-100 %) | 0.994 |
| L = 100 uH, ZARC 0.7 tau = 10 us | **0.000** (-100 %) | 0.995 |
| L = 100 uH, RC tau = 100 us | 8.03 (-20 %) | 1.000 |
| L = 10 uH, RC tau = 100 us | 8.59 (-14 %) | 0.998 |

Jediné varování je "High inductance", které se týká skutečné hodnoty L, ne
chyby R_inf. Pokud se tau neurčí z maxima, ale projde se (1D scan minimalizující
reziduum), vyjde u C -0.4 % a u A3 +5.3 %. Samotný model R-L-K je tedy
v pořádku, chybná je volba tau.

### R3 - Magická konstanta `R_inf = 1.0` Ohm je dosažitelná (Vysoká)

`rlk_fit.py:247-249`. Pokud polynom extrapoluje k R_inf <= 0, kód vrátí
pevně 1.0 Ohm bez ohledu na měřítko dat a bez varování. Mezi 3000 náhodnými
kapacitními spektry (Rs 1e-3 až 10 Ohm, R 1 až 1e5 Ohm, n 0.4 až 1, šum 0/1 %)
se tato větev spustila ve **201 případech (6.7 %)**, například:

- Rs = 3.41 Ohm, R = 755 Ohm, n = 0.55 -> R_inf = 1.0 (-71 %), `warnings=[]`
- Rs = 0.0043 Ohm, R = 104 Ohm, n = 0.49 -> R_inf = 1.0 (+23000 %), `warnings=[]`

V `RINF_ESTIMATION.md` je tato pojistka popsána ("clamp to 1 Ohm"), ale
magická hodnota není nikde zdůvodněná, což porušuje pravidla v CLAUDE.md.

### R4 - Tiché fallbacky (Vysoká)

a) **`hf_fallback` bez varování** (`rlk_fit.py:250-253`). Na
`example/real_gamry_example.DTA` (Im(f_max) = -1370 Ohm, Re = 826 Ohm, fáze -59 deg)
vrací `--ri-fit` 825.9 Ohm, tedy prostě Re(Z) na f_max. `fit_success=True` a žádné
varování, ačkoli se oblouk zjevně ani neblíží reálné ose. Stejně je to
u šumového případu (L = 0.1 uH, oblouk 160 kHz, šum 0.5 %): +25 % bez varování.

b) **Fallback po chybě fitu bere medián celého spektra** (`rlk_fit.py:389-393`,
podobně `:377-383`). `np.median(Z.real)` přes všechny body padne doprostřed
oblouku. U spektra Rs = 10 Ohm + RC 100 Ohm vrací 30.2 Ohm (+202 %). DRT
(`estimation.py:150-152`) přitom počítá medián jen z nejvyšších frekvencí.
Obsahují-li data NaN, fallback vrátí `nan`.

c) **Handler ignoruje `fit_success`** (`cli/handlers/rinf.py:57`). Při
selhání fitu CLI vypíše:
```
R_inf = 30.211 Ohm (0 HF points)
  For comparison: median = 10.000 Ohm (diff: +20.211 Ohm, +202.1%)
  forced
```
Chybná hodnota se předá do DRT jako `r_inf_preset`, přestože správný medián
je o řádek níž. Varování je jen holý text výjimky bez kontextu.

### R5 - Polynom Re(Im) je zkreslený pro CPE a nestabilní na šum (Střední)

`rlk_fit.py:244`. Pro oblouk s CPE platí na HF konci Re - Rs ~ |Im|^(1/n),
což kvadratický polynom nepostihne. Bez šumu:

| n | 1.0 | 0.9 | 0.8 | 0.7 | D (n = 0.6, Rs = 0.05) |
|---|---|---|---|---|---|
| chyba | +0.1 % | +0.4 % | +1.4 % | +4.0 % | +26.7 % |

Se šumem 1 % je případ D **+526 % +- 770 %**, tedy odhad je čistě náhodný.
Přesto je tato větev ze současných metod pro kapacitní data nejlepší:
v případech B/D je medián horší o jeden až dva řády.

### R6 - Chybí test identifikovatelnosti (Střední)

`--ri-fit` vždy vrátí číslo a nikdy neřekne, že R_inf z dat určit nelze.
Na `example/example_eis_data.csv` (Rs = 10 Ohm, fáze na f_max -48 deg) se
mýlí všechny metody: medián +3262 %, `--ri-fit` +258 % (s R^2 = -5.9), žádné
varování. README a `RINF_ESTIMATION.md` přitom doporučují `--ri-fit` právě pro
případ "an arc that is not closed", kde extrapolace z horní dekády principiálně
selhat může.

Fáze na f_max jako detektor nestačí. Případ A1 (oblouk nad f_max) má fázi jen
-2.9 deg a chybu +996 %, a naopak C má fázi +81 deg a lze ho spolehlivě vyřešit.
Funkční detektor je až směrodatná chyba R_inf z nelineárního fitu (návrh N2).

### R7 - R^2 ve dvou větvích nic neměří (Střední)

`rlk_fit.py:193-196, 261-264`. Ve větvích průchodu nulou a polynomu se R^2
počítá pro konstantní "model" Re = R_inf, takže je z principu záporný:
-2.58 (`EISPOT-test1.DTA`), -2.73 (`real_gamry_example.DTA`), -5.92
(`example_eis_data.csv`), až -7.8 (syntetický Warburg). Handler ho skryje
(`if R_squared > 0`), ale diagnostický graf ho ukazuje jako "QUALITY METRICS"
a dokumentace ho uvádí jako "quality over the fitted window". V rozbité větvi
R-L-K (R2) je naopak R^2 = 0.994 u výsledku -100 %. Metrika tedy mate
v obou směrech.

### R8 - Dokumentace popisuje chování, které kód nemá (Střední)

`doc/RINF_ESTIMATION.md` a README (`--ri-fit`):
- průchod nulou je popsán jako "strongest case" (viz R1),
- varování "`L < 0` - non-physical negative inductance" nemůže nastat (R10),
- "If that window happens to be empty, it falls back to the whole dataset"
  také nemůže nastat (R10),
- `R_squared` jako "quality over the fitted window" (R7),
- "`estimate_rinf_with_inductance()` never raises" neplatí: vstup jako
  Python `list` vyhodí `AttributeError` (`Z.real`) ještě před `try`. Pole
  nestejné délky naopak výjimku nevyhodí a tiše skončí ve fallbacku (medián
  celého spektra, v testu 107.4 Ohm při Rs = 10 Ohm),
- u pojistky 1 Ohm chybí zdůvodnění (R3).

### R9 - Okrajové případy (Nízká)

a) **Bod s Im == 0 přesně** (`rlk_fit.py:183`). Podmínka `im[i] * im[i+1] < 0`
bod ležící na nule nepovažuje za průchod a spektrum propadne do rozbité větve
R-L-K. Na datech s (12 - 1j, 11 + 0j, 10.5 + 1j) vyjde 10.04 místo 11.0.

b) **NaN v horní dekádě** projde `np.polyfit` a výsledek je `R_inf = nan`,
`fit_success=True`, `warnings=[]`. CLI je chráněné (`io/data_loading.py` NaN
odfiltruje), veřejné API (`eis_analysis.estimate_rinf_with_inductance`) ne.

c) **1-2 body v horní dekádě v indukční větvi.** Pro 3 neznámé (Rs, R_k, L)
a 2-4 reálné rovnice kontrola chybí. Kapacitní větev ji od opravy
z 2026-07-26 má (`len(im_hf) < 3`), indukční ne.

### R10 - Mrtvý kód (Nízká)

- `data_selection.py:172-180`: `n_in_decade == 0` nemůže nastat, f_max leží v dekádě vždy.
- `rlk_fit.py:399-400`: varování `L < 0`. NNLS (`allow_negative=False`) vynucuje L >= 0.
- `rlk_fit.py:172-173`: `else` v určení `behavior`. Tři předchozí podmínky pokrývají všechny kombinace.
- `rlk_fit.py:177-179`: druhé třídění dat, která už setřídil `sort_by_frequency` (`:149`).

### R11 - Duplicity a nekonzistence (Nízká)

- Dva prahy pro indukčnost: 1000 nH (`rlk_fit.py:348`) a 500 nH
  (`drt/estimation.py:164`). V knihovní cestě může vzniknout dvojí varování,
  v CLI cestě se druhé nikdy neuplatní (`preset` má přednost). Ani jeden práh
  není zdůvodněný.
- Výpočet HF mediánu je zkopírovaný v `cli/handlers/rinf.py:65-67` i `drt/estimation.py:150-152`.
- `_plot_rlk_fit` (`rlk_fit.py:424-429`) vybírá dekádu vlastním kódem místo `select_highest_decade`.
- `use_rl_fit=args.ri_fit` (`cli/handlers/drt.py:306`) v CLI nemá účinek.
  Pokud handler selže, DRT stejný fit tiše zopakuje.
- Lazy import přes globální proměnnou (`rlk_fit.py:21-31`) nahradí prostý
  lokální import ve funkci.
- `method` je složený řetězec (`capacitive_hf_fallback_highest_decade_11points`),
  takže se musí parsovat přes `in`/`startswith`. Počet bodů je přitom už v `n_points_used`.

### R12 - Pokrytí testy (Nízká)

- `tests/test_rinf_estimation.py` (4 testy) hlídá jen regresi polynomu pro málo
  bodů. Větve průchodu nulou a R-L-K nemají **žádný test přesnosti**,
  `test_inductive_data_unaffected` kontroluje jen název metody.
- `tests/test_cli_integration.py:339` (`test_rinf_estimation`) přijme i `None`,
  takže nic neověřuje.

---

## Srovnání metod

Chyba R_inf v %, data bez šumu:

| Případ | medián | `--ri-fit` | R-L-K, tau scan | R-L-2K | Lin-KK R_s | **R-L-ZARC, 2 dek.** |
|---|---|---|---|---|---|---|
| flat | +0.0 | +0.0 | -0.0 | -0.0 | -6.4 | **-0.0** |
| pure R-L | +0.0 | -0.0 | -0.0 | -0.0 | +0.0 | **+0.0** |
| A2 | +0.0 | +0.1 | -0.0 | -0.0 | +32.2 | **+0.0** |
| A3 | +59.8 | +101.7 | +5.3 | -1.1 | +91.9 | **-0.0** |
| C | +0.6 | +10.0 | -0.4 | -0.0 | +203.4 | **+0.0** |
| C2 | +37.5 | +89.0 | +41.4 | +22.5 | +10.1 | **+0.0** |
| B n = 1 | +0.6 | +0.1 | -0.1 | -0.0 | +57.1 | **-0.0** |
| B n = 0.8 | +18.4 | +1.4 | +13.1 | +6.2 | +3.8 | **-0.0** |
| B n = 0.7 | +37.5 | +4.0 | +26.4 | +12.1 | +10.1 | **-0.0** |
| B2 Warburg | +0.3 | -0.0 | +0.3 | +0.2 | +0.3 | **+0.0** |
| D | +3289 | +26.7 | +2375 | +1415 | +1155 | **-0.0** |
| A1 (oblouk nad f_max) | +998 | +996 | +0.6 | -40.2 | +992 | **-20.4** |

Šum 1 % (30 seedů), bias / směrodatná odchylka v %:

| Případ | medián | `--ri-fit` | R-L-K, tau scan | R-L-2K | R-L-ZARC, 2 dek. |
|---|---|---|---|---|---|
| flat | -0.0 / 0.4 | +0.1 / 0.6 | -0.1 / 0.3 | -35.8 / 37.1 | **-0.0 / 0.2** |
| C | +0.8 / 1.3 | +10.1 / 0.7 | -0.1 / 0.4 | -6.7 / 23.7 | **-0.2 / 0.4** |
| C2 | +37.5 / 1.9 | +89.0 / 1.1 | +41.3 / 0.5 | +23.2 / 5.0 | **-0.2 / 1.6** |
| B n = 0.8 | +18.2 / 0.9 | +1.2 / 2.3 | +13.6 / 0.9 | +6.1 / 0.5 | **-0.1 / 0.6** |
| D | +3280 / 34 | +526 / 770 | +2373 / 15 | +1390 / 90 | -22 / 43 |

Kandidáti:
- **R-L-K, tau scan:** současný model, tau hledané 1D scanem přes 120 hodnot (+-2 dekády kolem okna). Opraví R1/R2 pro RC oblouky, na CPE nestačí.
- **R-L-2K:** dva prvky K, scan přes páry tau. Nestabilní na šum (flat: +-37 %). **Zamítnuto.**
- **Lin-KK R_s:** `lin_kk_native(fit_type='complex', include_L=True)` na celém spektru, výsledkem `elements[0]`. Soustavně zkreslené (C +203 %). **Zamítnuto.**
- **R-L-ZARC:** Rs + jwL + R/(1 + (jw tau)^n), vážení modulem, meze L >= 0 a 0.3 <= n <= 1, multistart 9 x 3. Při fitu přes horní dvě dekády je ve všech případech nejlepší.

**Výhrada:** syntetická data jsou generovaná z modelu ZARC, takže kandidát
R-L-ZARC tu fituje "svůj" model a výsledek je optimistický. Na
`example_eis_data.csv` (dva překrývající se CPE, Rs/R0 = 1e-4) dává +153 %
a ostatní metody selhávají také. Právě proto je nutný návrh N2.

**Reálná data** (skutečné Rs neznámé, uvedeno pro rozptyl metod):

| Soubor | fáze(f_max) | medián | `--ri-fit` | R-L-K, tau scan | R-L-ZARC 2 dek. (sigma/R) |
|---|---|---|---|---|---|
| EISPOT-test1.DTA | -32 deg | 1.446 | 1.180 (polynom) | 1.365 | 1.082 (0.9 %) |
| real_gamry_example.DTA | -59 deg | 1402 | 825.9 (tichý hf_fallback) | 213.5 | ~0 (sigma astronomická -> neurčitelné) |
| example_eis_data.csv (Rs = 10) | -48 deg | 336.2 | 35.8 | 288.7 | 25.3 (14.4 %) |

U `EISPOT-test1.DTA` se metody rozcházejí o 30 % (1.08 až 1.45 Ohm) a pravdu
z dostupných dat určit nelze. Doporučuji ověřit proti nezávislému fitu celého
spektra.

**Dopad na DRT.** R_inf se odečítá před řešením, takže jeho chyba přejde do R_pol
zhruba 1:1 v Ohmech (případ C: R_inf 10.0 -> R_pol 100.5, R_inf 16.0 -> R_pol 94.6).
Relativní dopad na R_pol roste s poměrem Rs/R_pol. Posun ve tvaru gamma(tau)
na krátkých tau jsem v tomto auditu nekvantifikoval.

---

## Návrhy na vylepšení

Návrhy jsou seřazené podle přínosu. Varianta A je doporučená, varianta B je
nejmenší změna, která odstraní vážné nálezy. Obě staví na existujícím kódu.

### Varianta A (doporučená)

**N1 - Jedna metoda místo tří větví: R-L-(R|Q) přes horní dvě dekády.**
Stávající `fit_equivalent_circuit` (`fitting/circuit.py:271`) s obvodem
`R - L - (R | Q)` na datech s `f >= f_max/100` vrátil pro C2 i B n = 0.7 přesně
10.0 Ohm. Vlastní fitter tedy psát netřeba, stačí vhodný počáteční odhad
(Rs ~ min Re, L ~ Im(f_max)/w_max je-li Im > 0, jinak malé kladné,
tau ~ 1/w v okolí maxima -Im) a případně `fit_circuit_multistart`
(`fitting/multistart.py:158`). R_inf = `params_opt[0]`.
Zmizí tak R1, R2, R3, R5, R7 a R9a, protože větve přestanou existovat.
Šířka okna dvou dekád je kompromis. Při jedné dekádě měl ZARC na případu A1
-6.8 % místo -20.4 %, ale na datech se šumem byl výrazně rozptýlenější
(flat: +-2.1 % proti +-0.2 %, D: +-53 % proti +-43 %). Pokud má horní dekáda
málo bodů, rozšíření okna pomáhá i u R9c.

**N2 - Identifikovatelnost z `params_stderr`.** `FitResult` už vrací
`params_stderr` i `condition_number`. Navrhuji vracet `R_inf_stderr` a při
`stderr/R_inf > 5 %` přidat varování, že R_inf nelze z horní části spektra
spolehlivě určit, a nastavit příznak, podle kterého se CLI rozhodne. Naměřené
sigma/R (šum 1 %):

| Případ | sigma/R | skutečná chyba |
|---|---|---|
| flat, pure R-L, A2, C, C2, B, B2 | 0.2-1.1 % | <= 1.6 % |
| A3 | 2.4 % | -2.9 % |
| EISPOT-test1.DTA | 0.9 % | ? |
| example_eis_data.csv | **14.4 %** | +153 % |
| D | **77 %** | -48 % |
| A1, real_gamry_example.DTA | **> 1e20 %** | -100 % / ? |

Práh 5 % odděluje obě skupiny s rezervou 2x nahoru i dolů. U nesprávně
zvoleného modelu (CSV: 14 % proti skutečným 153 %) sigma chybu podhodnocuje,
je to tedy **příznak, ne chybový interval**, a tak ho je třeba dokumentovat.
Případ A1 (oblouk celý nad f_max) nepozná žádná metoda z dat samotných.
Dokumentace by to měla říct výslovně.

**N3 - Jednotná pravidla pro fallbacky (R4).**
- Každý fallback přidá do `warnings` záznam s tím, co se stalo a proč.
- Jediný fallback je HF medián, který už používá DRT. Výpočet vytáhnout do
  jedné funkce a volat ho z `rinf_estimation` i z `drt/estimation.py`
  (odstraní i duplicitu z R11).
- Handler (`cli/handlers/rinf.py`) při `fit_success=False` nebo varování
  z N2 vypíše obě hodnoty a do DRT předá medián, ne neúspěšný fit, nebo
  alespoň varování zopakuje vedle hodnoty. Rozhodnout by měl autor.
- Zrušit konstantu 1.0 Ohm (R3).

**N4 - Validace vstupu na veřejném API (R9b).** Nekonečné hodnoty a NaN
odmítnout nebo odfiltrovat na vstupu `estimate_rinf_with_inductance` stejně
jako v `io/data_loading.py`. Vstup převést přes `np.asarray` a kontrolovat shodnou délku `frequencies` a `Z`
(dnes nestejné délky tiše skončí ve fallbacku, viz R8).

**N5 - Úklid (R10, R11).** Odstranit mrtvé větve, jeden zdůvodněný práh pro L,
`_plot_rlk_fit` přes `select_highest_decade`, lokální import, `method` jako
krátký enum (`'rlc_fit' | 'median_fallback' | ...`) a počet bodů jen
v `n_points_used`. Po N1 se bude `RLKFitResult` (dnes 21 polí, polovina
specifických pro jednotlivé větve) dát zmenšit. Změna názvu funkce
`estimate_rinf_with_inductance` by rozbila veřejné API, proto jen alias
a deprecation.

**N6 - Testy (R12).** Jako základ regresní sady poslouží testy z commitu
`6a7f8f1`: 10 případů s tolerancí 2 % bez šumu, doplněných o test N2
(A1/D/CSV musí nést varování a flat/C nesmí). CLI test
`test_rinf_estimation` musí vyžadovat `R_inf is not None` a známou toleranci.
Náhodná data vždy se seedem (viz `generate_synthetic_data`).

**N7 - Dokumentace (R8).** Přepsat `RINF_ESTIMATION.md` a popis `--ri-fit`
v README podle nové metody. Výslovně uvést, kdy `--ri-fit` pomoct nemůže
(oblouk nad f_max, fáze na f_max desítky stupňů s CPE), a vysvětlit význam
`R_inf_stderr`.

### Varianta B (nejmenší změna)

Pokud je N1 příliš velký zásah:
1. Zrušit větev průchodu nulou a pro jakákoli data s indukční částí použít
   R-L-K se **scanem tau** (1D, ~120 hodnot, NNLS uvnitř; vše existuje
   v `estimate_R_linear`). C: +10 % -> -0.4 %, A3: +102 % -> +5.3 %.
2. Polynom ponechat pro čistě kapacitní data, ale bez konstanty 1.0 Ohm
   a s varováním u každé pojistky.
3. N3 (fallbacky) a N4 beze změny.

Tato varianta odstraní R1, R2, R3, R4 a R9a. **Nevyřeší** C2 (+41 %), R5 ani R6
a identifikovatelnost bude dál chybět.

---

## V pořádku (ověřeno)

- Řešení R-L-K přes `estimate_R_linear` je lineárně správné: pro čisté R-L
  vrací Rs i L přesně (10.000 Ohm, 100.0 nH).
- Pro čistě kapacitní data je polynom Re(Im) ze současných metod nejlepší:
  B n = 0.8 +1.4 % proti +18.4 % u mediánu. Oprava pro málo bodů
  z 2026-07-26 funguje.
- `RLKFitResult.failed()` dává všem cestám jednotný tvar výsledku.
- Knihovní část nic neloguje a vše vrací v datech, v souladu
  s `CLI_OUTPUT_UNIFICATION.md`.
- CLI loader odfiltruje NaN a nekonečna, takže R9b se v CLI neprojeví.

---

## Priority

| Priorita | Nálezy | Návrh |
|---|---|---|
| 1 | R1, R2 (chybné výsledky bez signálu) | N1, nebo B.1 |
| 2 | R3, R4 (tiché fallbacky, magická konstanta) | N3 |
| 3 | R6 (identifikovatelnost) | N2 |
| 4 | R8, R12 (dokumentace, testy) | N6, N7 |
| 5 | R5, R7, R9-R11 | po N1 zčásti zmizí, zbytek N4, N5 |
