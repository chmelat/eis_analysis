# Kriticka analyza modulu Z-HIT validace

Datum: 2026-09-28, verze v0.43.1

Rozsah: `eis_analysis/validation/zhit.py` (477 radku), pouziti v
`eis_analysis/cli/handlers/validation.py`. Modul nema vlastni testovaci soubor,
testuje se jen neprimo pres `tests/test_outliers.py`.

Nalezy oznacene "overeno" byly reprodukovany skriptem na syntetickych datech
bez sumu (71 bodu, 10^-2 .. 10^5 Hz): `R+RC` (R0 = 10 Ohm, R = 100 Ohm,
tau = 1e-2 s nebo 1 s) a `R+CPE-RC` (n = 0.8).

## Celkovy dojem

Samotna rekonstrukce je spravne: znamenko korekce 2. radu gamma = -pi/6 sedi
(s +pi/6 vychazi chyba 17 %, bez korekce 8 %, s -pi/6 0.9 %). Problem je v tom,
co se deje kolem ni. Vychozi kotva v jedinem bode prenasi lokalni chybu
aproximace na cele spektrum, `optimize_offset` v jedne ceste nefunguje, vstup se
nekontroluje a real/imag rezidua vypadaji jako nezavisla kontrola, prestoze
nejsou.

## 1. Chyby s dopadem na vysledek (overeno)

### 1.1 Kotva v jedinem bode udela z idealniho RC "acceptable"

**Stav:** opraveno ve v0.44.0, posun je median pres cele spektrum.

`zhit.py:376-382`, `zhit.py:285`. Aproximace Z-HIT 2. radu ma sama o sobe chybu,
na idealnim R+RC az 2.8 % na vrcholu faze. Vychozi `ref_freq` je geometricky
stred rozsahu a ten casto lezi prave na charakteristicke frekvenci. Lokalni
chyba se pak jako konstantni posun prenese na cele spektrum:

| Data (bez sumu) | ref = geom. stred | ref na okraji | `optimize_offset` |
|---|---|---|---|
| R+RC, tau = 1e-2 | **1.72 %** ("acceptable") | 0.63 % | 0.70 % |
| R+RC, tau = 1 | 0.76 % | 0.60 % | 0.60 % |
| R+CPE-RC, n = 0.8 | 0.77 % | 0.29 % | 0.34 % |

Stejne se to projevi, kdyz sum nese bod na referencni frekvenci. Jeden bod
posunuty o 5 % zvedne median rezidui na 3.3 %, s `optimize_offset=True` je
median 0.29 %. Vychozi nastaveni je tedy to horsi z obou.

Navrh: vychozi posun urcovat pres vazeny fit (`optimize_offset=True`) nebo
medianem rozdilu ln|Z_exp| - ln|Z_recon|, ne jednim bodem.

### 1.2 `optimize_offset` se ignoruje pri `use_second_order=False`

**Stav:** opraveno ve v0.44.0, cesta 1. radu (`use_second_order`) odstranena.

`zhit.py:261-262`. Cesta 1. radu vraci vysledek driv, nez se offset spocita.
Overeno: pri `optimize_offset=True` zustane chyba na kotve 0.1, nic se
neoptimalizuje. Obdobne `optimize_offset=True` bez `ln_Z_exp` tise spadne zpet
na pevny bod, bez varovani.

### 1.3 Vlastni chyba metody je radove stejna jako prahy kvality

**Stav:** analyzovano 2026-09-28 (v0.44.0), neopraveno. Pricina je v aproximaci
samotne, ne v implementaci. Nize je podrobny rozbor.

I po odecteni konstanty zbyva na idealnim RC prumerne 0.6-0.9 % (maximum 2.8 %
u f = 16 Hz). Podle `_quality_label` je to "good", nikdy "excellent".

#### 1.3.1 Odvozeni: Z-HIT je useknuta asymptoticka rada

V promenne x = ln(omega) plati pro minimalne-fazovou impedanci Bodeho vztah
mezi A(x) = ln|Z| a fazi phi(x):

    phi(x0) = (1/pi) * integral[ A'(x) * ln coth(|x - x0| / 2) dx ]

Fourierova transformace jadra je `pi * tanh(pi*k/2) / k`, takze ve spektralni
oblasti `phi^ = i * tanh(pi*k/2) * A^` a obracene

    A^ = -i * coth(pi*k/2) * phi^

Laurentuv rozvoj `coth z = 1/z + z/3 - z^3/45 + 2 z^5/945 - ...` se `z = pi*k/2`
a prevod mocnin k na derivace (`ik <-> d/dx`) dava:

    ln|Z(x)| = C + (2/pi) * integral[phi dx]
                 - (pi/6)       * phi'(x)
                 - (pi^3/360)   * phi'''(x)
                 - (pi^5/15120) * phi^(5)(x)
                 - ...

Prvni dva cleny presne odpovidaji implementaci (gamma = -pi/6), coz odvozeni
potvrzuje. Kod konci za phi'.

Rozvoj coth konverguje jen pro |z| < pi, tj. |k| < 2. Spektrum faze relaxace
tam nekonci: u Debye je phi'(x) ~ sech(x), jehoz transformace klesa jen jako
exp(-pi |k| / 2). Rada je proto **asymptoticka, ne konvergentni**. Dalsi cleny
nejdriv pomahaji a pak zhorsuji, s optimem po 2-3 clenech.

#### 1.3.2 Mereni: chyba useknuti bez vlivu diskretizace

Derivace faze byly spocteny na huste siti (2000 bodu/dekadu, 10^-4 .. 10^7 Hz)
a chyba byla vyhodnocena na 10^-2 .. 10^5 Hz, tedy daleko od okraju husteho
rozsahu. Posun je median, jako ve v0.44.0. Hodnoty jsou mean / max |res| v %:

| obvod | 2 cleny (kod) | + phi''' | + phi^(5) |
|---|---|---|---|
| R+RC (Debye) | 0.71 / 3.19 | 0.31 / 1.86 | 0.47 / 3.22 |
| R+CPE, n = 0.9 | 0.49 / 1.96 | 0.17 / 0.85 | 0.17 / 1.10 |
| R+CPE, n = 0.8 | 0.33 / 1.16 | 0.09 / 0.37 | 0.06 / 0.38 |
| R+CPE, n = 0.6 | 0.12 / 0.35 | 0.02 / 0.06 | 0.02 / 0.12 |
| 2x RC, tau 1e-2 a 1e-3 | 0.79 / 3.24 | 0.50 / 2.45 | 0.79 / 4.61 |
| RC + Warburg | 0.68 / 3.23 | 0.29 / 1.88 | 0.42 / 3.22 |

Obvody: R0 = 10 Ohm, R = 100 Ohm, tau = 1e-2 s, Warburg `30/sqrt(j omega)`.

Zavery:

- Chyba roste s tim, jak ostre se faze ohyba. Debye je nejhorsi bezny pripad,
  CPE s nizsim n je nekolikrat lepsi. Realna data s n ~ 0.8-0.9 lezi v pasmu
  0.3-0.5 %.
- Maximum chyby lezi na relaxaci (vrchol |phi'''|), ne na okraji spektra.
- Pridani phi''' by na cistych datech chybu snizilo 2-4x. Clen phi^(5) ji u
  ostrych relaxaci opet zvysi, coz je typicke chovani asymptoticke rady.
- Pri fitu volneho koeficientu u phi''' vychazi -0.052 az -0.076 misto
  -pi^3/360 = -0.086. Fit pohlcuje vyssi cleny, konstantni "lepsi gamma_3"
  neexistuje.

#### 1.3.3 Diskretizace neni pricina

Implementace (`cumulative_trapezoid` + `np.gradient`) na R+RC (Debye),
10^-2 .. 10^5 Hz:

| body/dekadu | mean % | max % |
|---|---|---|
| 5 | 0.44 | 1.68 |
| 10 | 0.63 | 2.80 |
| 20 | 0.69 | 3.09 |
| 50 | 0.70 | 3.18 |
| presne derivace (1.3.2) | 0.71 | 3.19 |

Pri 10 bodech/dekadu se implementace od presneho vzorce prakticky nelisi.
Hrubsi sit chybu dokonce snizuje, protoze diferencni chyba centralni diference
(+h^2/6 * phi''') a lichobeznikove integrace pusobi proti chybe useknuti. Je to
nahodna kompenzace, ne vlastnost, na kterou by se dalo spolehat.

#### 1.3.4 Vysvetleni chyby na okraji v kodu je chybne

CLI varovani pro `--fit-on all` (`cli/handlers/validation.py`), README (sekce
`--fit-on`) i docstring `test_reconstruction_is_least_accurate_at_the_high_frequency_edge`
pripisuji ~1 % v nejvyssi dekade jednostrannym diferencim `np.gradient` na
okraji. Na tomtez referencnim spektru (`tests/test_zhit_fit_on.py`,
1 mHz .. 100 kHz) vychazi:

| | implementace | presny 2-clenny vzorec (bez `np.gradient`) |
|---|---|---|
| nejvyssi dekada | 1.03 % | 1.13 % |
| nejnizsi dekada | 0.08 % | 0.08 % |

HF relaxace toho spektra, `R(100) | Q(1e-6, 0.9)`, ma
tau0 = (R*Q)^(1/n) = 3.6e-5 s, tedy f_peak = 4.4 kHz, prave v nejvyssi dekade.
Chyba tam pochazi z useknuti rady, ne z okraje.

Vliv okraje samotneho (Debye, 71 bodu, 10^-2 .. 10^5 Hz, posledni tri body,
chyba v %):

| f_peak | implementace | presny vzorec |
|---|---|---|
| 1e3 Hz | -0.54, -0.59, -0.01 | -0.55, -0.61, -0.60 |
| 1e4 Hz | 1.99, 2.53, 1.92 | 2.23, 2.85, 3.16 |
| 3e4 Hz | -1.45, -0.72, -2.75 | -1.60, -0.79, 0.02 |

Jednostranna diference meni jen posledni bod, a to obema smery (jednou chybu
zmensi, jindy zvetsi) az o ~2.8 %, a jen kdyz relaxace lezi primo u okraje.
Kdyz je relaxace daleko (f_peak = 10 Hz), maji okrajove body chybu ~0.01 %.

Varovani tedy spravne rika, ze v HF oblasti muze byt chyba ~1 %, ale z
nespravneho duvodu. Vypisuje se pokazde, i kdyz na okraji zadna relaxace neni,
a naopak nevaruje pred stejne velkou chybou uprostred spektra.

#### 1.3.5 Proc nepridat clen phi'''

Treti derivace zesiluje sum zhruba jako 1/h^3 (h = krok v ln omega, 0.23 pri
10 bodech/dekadu). Debye, 10 bodu/dekadu, 1% komplexni sum, 20 seedu, prumerna
chyba rekonstrukce vuci *skutecnemu* |Z|:

| rad | chyba |
|---|---|
| 1. (jen integral) | 6.02 % |
| 2. (-pi/6 phi', soucasny kod) | 1.69 % |
| 3. (+ -pi^3/360 phi''') | 2.26 % |

Na realnych datech je tedy soucasne useknuti za phi' optimalni. Clen phi'''
by pomohl jen na datech s velmi malym sumem nebo po vyhlazeni faze, a to by
byla dalsi volba k ladeni.

#### 1.3.6 Dusledky

1. **Podlaha chyby na K-K kompatibilnich datech** je 0.1-0.7 % prumerne a az
   ~3 % v jednotlivych bodech u relaxace. Zavisi na tvaru dat, ne na sumu.
   Idealni RC proto podle `_quality_label` (0.5 / 1 / 2.5 / 5 %, kalibrovano
   pro Lin-KK) nikdy nedostane "excellent".
2. **Prah pro odlehle body je bezpecny.** `find_outliers` pouziva 5 %, takze
   cista data false positive nedaji. Rezerva u vrcholu Debye je ale jen ~1.8 %.
3. **`--fit-on`**: ostre relaxace se v rekonstrukci systematicky zkresli az o
   ~3 % v |Z| kolem f_peak. Na fitu se to zatim projevuje jen 0.12 % v
   odporech (`test_reconstruction_error_floor_on_undisturbed_data`).

#### 1.3.7 Navrhy

1. **Stav:** CLI varovani a docstring testu opraveny ve v0.44.1, README
   zatim ne.
   Opravit vysvetleni v CLI varovani, README a docstringu testu: chyba sedi
   tam, kde se faze ohyba, ne na okraji. Varovani pro `--fit-on all` bud
   podminit velikosti |pi^3/360 * phi'''| v HF oblasti, nebo ho preformulovat
   obecne.
2. Z-HIT by mel mit vlastni prahy kvality, nebo aspon vypisovat odhad podlahy
   metody. Levny odhad je velikost prvniho vynechaneho clenu
   |pi^3/360 * phi'''| v kazdem bode. Z derivace zasumene faze je ale sam
   zasumeny, takze by se musel pocitat z vyhlazene faze nebo jen jako median
   pres oblast.
3. Presna alternativa je aplikovat nasobitel `-i*coth(pi*k/2)` pres FFT v
   ln omega. Vyzaduje ale extrapolaci faze mimo mereny rozsah, coz je presne
   to, cemu se Ehm et al. lokalnim rozvojem vyhybaji. Nedoporucuji.

## 2. Stredne vazne

### 2.1 Real/imag rezidua nenesou novou informaci

`zhit.py:414-421`. `Z_fit` prebira zmerenou fazi, takze
`res_real = res_mag * cos(phi)` presne (overeno `np.allclose`). Graf
"Real/Imag", vypis v CLI i `pseudo_chisqr` jsou jen jinak zobrazene reziduum
modulu, ne nezavisla kontrola.

Navic graf ukazuje real/imag s pevnymi carami +-5 %, zatimco `is_valid`
rozhoduje podle modulu a podle `quality_threshold`.

### 2.2 Chybi validace vstupu (overeno)

| Vstup | Chovani |
|---|---|
| 0 bodu | `IndexError` na `zhit.py:377`, mimo `try`, spadne |
| 1 bod | `success=False` s nesrozumitelnou hlaskou z `np.gradient` |
| 2 body | "projde" s reziduem 0.005 % (nesmysl) |
| duplicitni frekvence | `success=True`, NaN v reziduich, jen RuntimeWarning |
| f <= 0 nebo \|Z\| = 0 | `log` -> NaN/inf |

Siroke `except Exception` (`zhit.py:394`) NaN nezachyti, protoze numpy
nevyhazuje vyjimku, jen varuje. Cesta s prazdnym vysledkem tak prakticky chrani
jen pred pripady, ktere by stejne skoncily chybou.

### 2.3 Poruseni projektovych pravidel

- Vykreslovani primo ve validacni funkci (`zhit.py:432-461`). Kazde volani
  vytvori figuru, coz vadi v davkach a bez displeje. CLAUDE.md vyzaduje oddelit
  vizualizaci od algoritmu.
- `logger.error` v knihovnim modulu misto pole `warnings` ve vysledku.
- `ZHITResult` nema `warnings: List[str]`.

## 3. Drobnosti

- **Docstring vs. kod** (`zhit.py:247-249`): uvadi "first_order - gamma * dphi"
  s "gamma ~ 0.2-0.5", kod dela `+ gamma` s pevnym gamma = -pi/6. Kod je
  spravne, plati jen dokumentace. Zapis "(2/pi) * H[phi]" (`zhit.py:13`, `:347`)
  neni presny, jde o integral, ne Hilbertovu transformaci.
- **Jednotky**: `residuals_mag` je v %, `residuals_real/imag` jako zlomek.
  Vlastnosti `mean_residual_*` to musi obchazet.
- **Nejblizsi bod k `ref_freq`** se hleda linearne (`zhit.py:380`), na
  log-skale by mela byt logaritmicka vzdalenost. V testovanych pripadech to
  vyslo stejne.
- **`estimate_noise_percent`** (Yrjana & Bobacka 2024) je odvozeny pro Lin-KK
  fit. Pro Z-HIT, kde se nic nefituje a rezidua jsou projekce modulu, neni
  pouziti zduvodnene.
- **`np.clip(weights, 0, 1)`** (`zhit.py:185`) nic nedela, Gaussovo okno je z
  definice v [0, 1].
- **Testy**: chybi test na idealnim RC (absolutni presnost), test
  `optimize_offset` v obou cestach a okrajove vstupy.

## 4. Priorita oprav

1. Vychozi offset pres vazeny fit nebo median misto jednoho bodu (1.1).
2. Aplikovat offset i v ceste 1. radu (1.2).
3. Validace vstupu (2.2).
4. Rozhodnout, co s real/imag rezidui v grafu a vypisu (2.1).
5. Opravit vysvetleni chyby na okraji v CLI, README a testu (1.3.4).
6. Srovnat prahy kvality s podlahou metody, nebo ji vypisovat (1.3.6, 1.3.7).
