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

I po odecteni konstanty zbyva na idealnim RC prumerne 0.6-0.9 % (maximum 2.8 %
u f = 16 Hz). Podle `_quality_label` je to "good", nikdy "excellent". Bud prahy
neodpovidaji teto implementaci, nebo chybi dalsi clen rozvoje (Ehm et al. 2001
uvadeji i vyssi liche derivace faze; presne koeficienty jsem neoveroval).

Komentar v `cli/handlers/validation.py` (vetev `--fit-on all`) uvadi 0.08 % v
nejnizsi dekade a 1.0 % v nejvyssi. S merenim vyse nesedi, ale jde o jina
referencni data, takze jen k provereni.

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
5. Srovnat prahy kvality s vlastni chybou metody (1.3).
