# Plan: sum Lin-KK se seriovym C pro fit/KK i u testu, ktery projde

Datum: 2026-11-02. Stav: **realizovano 2026-11-02, B zamitnuta** (viz Vysledek);
omezeni zapsano do README (Fit vs. noise).
Puvod: cold agent kolo 8, nalez Z1 (`doc/COLD_AGENT.md`, "Z1 podrobne").

## Cil

Zjistit merenim, jestli ma `Fit vs. noise` (sloupec fit/KK v porovnani
obvodu) brat sum Lin-KK se seriovym C i tehdy, kdyz zakladni test projde,
pokud je NF konec kapacitni. Vystupem je rozhodnuti s cisly: prevzit, nebo
zamitnout a zapsat omezeni do README. Zadna zmena `eis_analysis/` behem
experimentu.

## Problem

`kramers_kronig_validation` pousti beh se `include_C` jen pri **neuspesnem**
testu a zapornem sklonu -Z'' v nejnizsi dekade (`capacitive_end`). Komentar
v kodu predpoklada, ze test, ktery projde, ma sum blizko sumu dat. Na
bateriich NMC811 (Zenodo 15422339) to neplati: polonekonecna difuze
pokracuje pod f_min, realny fit Lin-KK ji nevidi a test projde se sumem
0.82-0.87 %, kdezto se seriovym C 0.14-0.18 %. fit/KK pak u obvodu bez
R||L ukazuje 0.8 ("na urovni sumu") misto 4.0-5.1, prestoze BIC ho odmita
o 167-185. Varovaly jen testy tvaru rezidui.

## Soucasny stav (k porovnani)

- `validation/kramers_kronig.py`: `capacitive_end = True` a
  `with_series_c = lin_kk_native(..., include_C=True)` jen pri
  `not result.is_valid and low_frequency_slope(...) < 0`.
- `cli/handlers/fit_noise.py: kk_noise_for_fit`: pri `capacitive_end` bere
  sum z `with_series_c`, jinak ze zakladniho behu (`noise_estimate_on`).
- `FIT_NOISE_RATIO_WARN = 2.0`; kalibrace 2026-10-13 (`fitting/config.py`):
  spravny obvod 0.60-1.02, chybejici oblouk 9-219, Q -> C 1.7-41, jen
  uzavrene oblouky.
- README (Fit vs. noise): "Lin-KK noise is only an upper bound - it contains
  drift - but so does the circuit fit, and the ratio cancels it."

## Varianty

- **A (dnes):** viz vyse.
- **B (navrh):** pri `low_frequency_slope < 0` pustit beh se `include_C`
  vzdy; pro fit/KK brat **mensi** ze dvou sumu (zakladni, se seriovym C).
  Oba modely jsou KK konzistentni, oba odhady jsou horni mezi sumu, mensi
  je tesnejsi horni mez. Priznak `capacitive_end` a vse, co na nem stoji
  (per-point kontrola v `outliers.py`, varovani, popisek grafu, vypis KK),
  zustava beze zmeny - meni se jen reference fit/KK.
- **C (jen pro srovnani, levne):** mensi z Lin-KK (A) a Z-HIT. Z-HIT je
  presny pro konstantni fazi, ale na ostrych obloucich nafouknuty (1 % na
  bezsumovem RC), takze sam o sobe neni kandidat; sloupec ukaze, jestli
  by stacil.

Merit se bude jen reference sumu: skript spocita oba behy Lin-KK
(`kramers_kronig_validation` s `include_C` False/True) a Z-HIT primo, kod
knihovny se nemeni.

## Mereni

Skript v `tmp/fitnoise_c/` (mimo git), max 4 procesy, cisla do
`tmp/fitnoise_c/results.json` a tabulek.

### S1: stress pripady (spravny a chybny obvod)

`tests/stress_cases.py`, rodiny `diffusion`, `blocking`, `anomalous`,
`oxide` (kapacitni konce) a jako kontrola `rc`, `cpe`; 40 pripadu na
rodinu, sum podle generatoru (0, 0.1, 1, 3 %). Na kazdem pripadu:

- **spravny obvod:** pravdivy vyraz fitovany `fit_equivalent_circuit`
  z pravdivych hodnot (LM, bez DE - zajima nas reference, ne optimalizator);
- **chybne obvody:** (1) bez nejmensiho oblouku, (2) nejvetsi Q -> C,
  (3) u Warburgu (W/Ws/Wo/Wa/Wat) prvek nahrazeny Q; start z pravdy
  zbytku, LM.

Pro A, B, C: pomer fit/KK. Metriky po rodine a urovni sumu:

- **falesny poplach:** podil spravnych obvodu s pomerem > 2 (A: kalibrace
  ~0);
- **detekce:** podil chybnych obvodu s pomerem > 2;
- rozdeleni pomeru spravneho obvodu (median, 95. percentil, max);
- u kolika pripadu B referenci meni (sklon < 0 a test projde) a o kolik.

### S2: drift

Komentar v kodu varuje: na uzavrenem konci seriove C jen pohlti drift. B
se spousti jen pri zapornem sklonu, ale fazovy drift ohne -Z'' nahoru
i na uzavrenem oblouku. Spektra: `drifting_zarc`
(`tests/test_kramers_kronig.py`, drift R 0 / 5 / 10 / 20 %) a ctyri drifty
z README "Finding the valid frequency range" (Re 25 %/dek., |Z|, faze
5 deg/dek. na uzavrenem oblouku, blokujici konec). Spravny (nedriftujici)
obvod fitovany na driftujici data; pomer A, B, C. Otazka: kolikrat B
zmeni "drift se vykrati" (pomer ~1) na pomer > 2. To muze byt i zadouci
(model data nepopisuje), ale je to zmena slibu v README a musi se
rozhodnout vedome.

### S3: realna data

- NMC811, vsech 54 clanku: sklon, sum A/B/C; pomer pro tri obvody
  (`R-L-(R|Q)-W`, `R-(R|L)-(R|Q)-W`, s druhym obloukem). DE je pomale -
  obvody fitovat jen na 10 clancich (los), sum na vsech 54.
- `example/` (EISPOT DTA, csv s difuzi `diffusion_*`): kde se reference
  zmeni, porovnat s vykladem v README (EISPOT-test1: dnes 9.3x sumu).

## Kriteria rozhodnuti

Prevzit B, pokud soucasne:

1. **S1, spravne obvody:** falesny poplach B v zadne rodine a urovni
   sumu <= 1 % a nikde vic nez A + 1 pripad; max pomer spravneho obvodu
   na kapacitnich koncich <= 2 pri sumu >= 0.1 %.
2. **S1, chybne obvody:** detekce na kapacitnich koncich vyssi nez u A
   (alespon u jednoho typu chyby o >= 10 procentnich bodu), jinde stejna.
3. **S2:** chovani pri driftu popsatelne jednou vetou v README; pokud B
   meni pomer spravneho obvodu nad 2 jen u driftu, ktery Lin-KK + C
   sam hlasi (test neprojde), je to prijatelne.
4. **S3:** NMC811 jako v kole 8 (agentuv obvod > 2, obvody s R||L < 2);
   `example/` bez noveho poplachu, ktery by README nevysvetlilo.

Jinak zamitnout a do README (Fit vs. noise) zapsat omezeni: pri difuzi
pokracujici za f_min je reference Lin-KK nafouknuta a fit/KK < 2 nic
nedokazuje; rozhoduji testy tvaru rezidui (runs test).

## Pri prevzeti (az po rozhodnuti, samostatny krok)

- Kontrakt a docstring `kramers_kronig_validation` (novy beh se C i pri
  uspechu; pole pro nej, aby `capacitive_end` zustal jen pro selhani),
  `kk_noise_for_fit`, README (Fit vs. noise, `--kk-series-c`),
  test v `tests/test_fit_noise_adequacy.py` (Warburg za f_min: obvod bez
  prvku > 2, spravny <= 2), CHANGELOG, verze.
- Cena: jeden beh Lin-KK navic jen na spektrech se zapornym NF sklonem.
- Pred commitem `/code-review`, pak stress `--check`.

## Odhad casu

S1: 6 rodin x 40 pripadu x (1 + 3) LM fity + 2x Lin-KK + Z-HIT, ~1-2 s na
pripad: ~10 min na 4 procesech. S2: sekundy. S3: sum 54x pod minutu, DE
10 clanku x 3 obvody ~10 min. Pred dlouhym behem se zeptat (review skriptu).

## Vysledek (2026-11-02)

Skript `tmp/fitnoise_c/measure.py` (mimo git), behy S1 20 s, S2 10 s,
S3 ~10 min. Rozhodnuti uzivatele: **B zamitnuta**, omezeni v README.

**S1 (240 stress pripadu, 40 na rodinu).** Spravny obvod: falesny poplach
A/B/C 0/0/0 ve vsech rodinach a urovnich sumu (max pomer 1.81, oxide 3 %).
B snizila referenci jen v 5 pripadech: stress ohony jsou tak silne, ze
zakladni Lin-KK selze a series C bezi uz dnes (`capacitive_end`); rezim NMC
(test projde, sum nafouknuty) stress skoro nepokryva. Detekce chybnych
obvodu A = B krome diffusion q_to_c 15 -> 17 ze 40 (+5 p. b.):
**kriterium 2 nesplneno.** C = A vsude (Z-HIT nikdy nizsi).

**S2 (drift, 0.2 % sum, 2 seedy).** Drift R oblouku, Re a |Z| na
uzavrenem ZARC: A = B 0.64-0.72 (drift se vykrati). Blokujici konec:
A = B (series C bezi uz dnes). Fazovy drift na uzavrenem ZARC ohne -Z''
nahoru (sklon -0.35 az -0.73), Lin-KK projde:

| drift faze | A | B | Lin-KK + C |
|---|---|---|---|
| 1 deg | 1.26-1.40 | 1.62-1.65 | projde, 0.19 % |
| 2 deg | 1.15-1.25 | 2.56-2.57 | projde, 0.21 %, C = 0.57 F |
| 5 deg | 1.00-3.26 | 4.25-4.38 | projde, 0.30 % |

Series C drift pohlti; ani ono, ani Z-HIT (`is_valid`, max |res_mag|
4.4 % pri 5 deg) drift nehlasi. Spravny obvod by B oznacila a vedla
k pridani C, ktere modeluje drift: **kriterium 3 nesplneno.**

**S3.** NMC811, 54 clanku: sklon -0.51 az -0.53, Lin-KK projde 54/54,
`capacitive_end` 0; reference A 0.78-1.12 %, B 0.12-0.23 %, C 0.24-0.40 %.
DE (seed 1) na 10 vylosovanych clancich:

| obvod | A | B | C |
|---|---|---|---|
| R-L-(R\|Q)-W (agent T16) | 0.64-0.88 | 4.3-5.3 | 2.2-2.7 |
| R-(R\|L)-(R\|Q)-W | 0.24-0.31 | 1.57-2.07 | 0.76-0.95 |
| R-(R\|L)-(R\|Q)-(R\|Q)-W | 0.13-0.18 | 0.87-1.18 | 0.45-0.54 |

`example/`: Warburgovy soubory 0.92-0.99 ve vsech variantach;
EISPOT-test1 5.9 (A) / 8.7 (B), uz dnes nad 2. **Kriterium 4 splneno.**

**Zaver.** Uvnitr okna nelze odlisit ohon pokracujici pod f_min od NF
fazoveho driftu: oba ohnou -Z'' nahoru a series C pohlti oba. B ma pravdu
na bateriich a myli se na driftu. README nyni rika, ze pomer < 2 na
kapacitnim konci, ktery Lin-KK projde, nic nedokazuje, rozhoduji testy
rezidui, a kdo vi, ze jde o skutecny ohon, pusti `--kk-series-c` (overeno:
ID11, agentuv obvod, pomer 4.5). C (mensi z Lin-KK a Z-HIT) by NMC
oznacila (2.2-2.7) a na drift nereagovala, ale je podlozena jen NMC;
samostatne mereni na spektrech, kde Z-HIT nadhodnocuje, neprobehlo.
