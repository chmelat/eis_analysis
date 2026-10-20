# Plan: kanonicke poradi zamenitelnych bloku po fitu

Stav: **ve fronte (2026-10-13, revize 2026-10-20), neimplementovano**.
Revize overena na v0.61.0 (97c86b8). Podnet: revize eis-fitting-protocol
(`doc/EIS_FITTING_PROTOCOL_COMPARISON.md`, nalez 1 a bod 3 "Co si odtud
vzit").

## Problem

Seriove zapojene bloky stejne struktury, napr. dva `(R|Q)` v
`R - (R|Q) - (R|Q)`, jsou zamenitelne: prohozeni jejich hodnot nemeni
impedanci. Fit tak muze vratit tytez oblouky v libovolnem poradi a ktery
z nich se jmenuje `R1|Q0` a ktery `R2|Q1` (labely se cisluji od nuly
zvlast pro kazdy typ parametru), je nahoda.

### DE: poradi je nahodne

Mereni 2026-10-13 (`fit_circuit_diffevo`, 0.5 % sum, 8 seedu na scenar,
stejne spektrum, stejny pocatecni odhad; pocet seedu, kdy rychlejsi oblouk
skoncil v prvnim bloku):

| Spektrum | odhad ve spravnem poradi | prohozeny odhad | oba bloky stejne |
|---|---|---|---|
| A: oblouky daleko (tau ~2e-2 s vs ~2e3 s) | 5/8 | 3/8 | 4/8 |
| B: oblouky blizko (pomer tau ~25) | 4/8 | 2/8 | 4/8 |

Konecna chyba fitu byla ve vsech behech stejna (0.617 %, resp. 0.630 %):
tentyz minimum, jen prohozene stitky. Pocatecni odhad poradi skoro
neovlivni - DE populace pokryva cele meze a x0 je jeden clen z mnoha.
Protoze vychozi `seed=None`, muze prohodit stitky i opakovany fit tehoz
spektra.

### LM: poradi drzi vstup

Mereni 2026-10-20 (`fit_equivalent_circuit`, `R-(R|Q)-(R|Q)`, tau 1e-4 s
a 0.5 s, 0.5 % sum, odhad = pravda x U(0.3, 3), 10 behu na scenar):
s rychlym obloukem v odhadu prvnim vysel rychly prvni 10/10, s pomalym
prvnim 0/10. Lokalni optimalizator z bloku neutece; poradi vysledku je
poradi, ktere zapsal uzivatel.

Pro LM tedy zmena neopravuje nahodu, ale **meni chovani**: kdo dnes zapise
pomaly oblouk prvni, dostane ho po zmene druhy. Rozhodnuti 2026-10-20:
presto vzdy rostouci tau (varianta a), protoze je objektivni, nezavisi na
vstupu a shoduje se s `auto_suggest`. Zamitnuta varianta b (poradi podle
tau v pocatecnim odhadu): respektuje zapis, ale pridava logiku a pri
shodnych tau v odhadu stejne konci u varianty a.

### Proc na tom zalezi

Pro jeden fit to nevadi. Vadi to pri srovnani fitu rady spekter (oxidove
serie, export do CSV/JSON a dalsi zpracovani): "R1" v jednom souboru muze
byt jiny fyzikalni proces nez "R1" ve druhem.

Cold-agent kolo 6 (T13, T14, nalez O4): na tomtez spektru byl s vahou
`modulus` R1 rychly proces, se `sqrt`/`uniform` pomaly - stitky meni i
volba vahy. README to od 2026-10-16 rika (porovnavat podle tau).

Uz se to projevilo dvakrat:

- `tests/test_zscope_benchmark.py` ma vlastni `_canonical()`, ktery vetve
  `Two_TimeConstants` radi podle tau pred srovnanim, protoze fit je na
  ZScope 2 % datech vratil prohozene (CHANGELOG v0.25.5).
- Stress test, invariant F3 (`doc/ROBUST_LOSS_SOFT_L1.md`, sekce 6,
  97c86b8): 2-3 % parametru lezi 5-1000 stderr od pravdy, temer vsechno
  prohozene bloky, ktere F3 porovnava podle pozice. Serazeni bloku podle
  tau zvedne pokryti 95% intervalu z 0.903 na 0.918.

## Navrh

Po optimalizaci, **pred** vypoctem kovariance, stavu mezi a varovani,
preradit hodnoty zamenitelnych bloku tak, aby tau rostlo zleva doprava
(nejrychlejsi proces prvni = nejvyssi frekvence, jak se cte Nyquist).
Stejne poradi uz pouziva `auto_suggest.py` (piky DRT razene podle
`peak_tau`), takze navrzeny obvod a jeho fit budou cist bloky stejne.
Stejne poradi pak vyjde i ve vystupu `arc_time_constants`.

Co poradi znamena: **poradi ve spektru, ne identitu procesu**. V serii
spekter, kde se tau dvou procesu krizi (teplota, starnuti oxidu), si
stitky v miste krizeni prohodi. Kdo potrebuje sledovat proces, porovnava
podle kapacity (C, resp. C_eff); oxidovy modul to tak dela
(`_element_size` v `analysis/oxide_elements.py`). Patri to do README.

### 1. Ktere bloky jsou zamenitelne

V kazdem seriovem **retezci** (rekurzivne, nejen na nejvyssi urovni) najit
skupiny clenu, ktere:

- maji **stejnou strukturu**: stejny strom typu prvku a zapojeni, bez
  ohledu na hodnoty (rekurzivni tuple `(typ uzlu, typy potomku)`);
- maji **stejnou masku fixovanych parametru** a u fixovanych **stejne
  hodnoty**. Blok s fixovanym `R = 100` a blok s volnym R zamenitelne
  nejsou - prohozeni by zmenilo model;
- maji **stejne meze** (`lower`, `upper`) po parametrech.
  `fit_equivalent_circuit(bounds=...)` je verejne API a meze bere po
  indexech; prohozeni bloku s ruznymi mezemi by mohlo dat hodnotu mimo meze
  jejiho noveho bloku. Bez `bounds` se meze generuji z labelu
  (`generate_simple_bounds`) a stejna struktura ma stejne meze sama;
  `log_search_bounds` v DE zavisi jen na mezich a labelech (overeno).

Retezec je Series i se Series vnorenou primo v ni: builder zprava
neslucuje, `(R|Q) - (R - (R|Q))` je `Series([P, Series([R, P])])` a oba P
jsou v serii. `_collect_arcs` v `time_constants.py` retezce takhle uz
spojuje - prevzit tu logiku, ne psat druhou.

Strukturalni klic citlivy na poradi prvku v Parallel bere `(R|Q)` a
`(Q|R)` jako ruzne bloky, takze se neseradi. Pro prvni verzi prijatelne
(builder ani `auto_suggest` takove smisene zapisy nevytvari).

Poznamka k pravidlu "fyzikalni kriterium, ne strukturalni": rozpoznani
symetrie je nutne strukturalni (symetrie je vlastnost vyrazu), ale klic
razeni je fyzikalni (tau). Strukturalni cast nerozhoduje o fyzikalnim
vyznamu prvku, jen o tom, co lze prohodit beze zmeny modelu.

### 2. Klic razeni: tau bloku

Prvni verze jen analyticky, pro bloky, kde je tau jednoznacne:

| Blok | tau |
|---|---|
| `(R\|C)`, `(G\|C)` | `R*C`, `C/G` |
| `(R\|Q)`, `(G\|Q)` | `(R*Q)^(1/n)` (maximum -Im v `omega = 1/tau`), R = 1/G |
| prvek K | jeho vlastni `tau` |

Vzorce **nepsat znovu**: `_arc_resistor` a `_arc` v `time_constants.py` je
uz pocitaji, vcetne G a ochran (n <= 0, nekladne R/C/Q, preteceni,
nekonecne tau). Vytahnout z `_arc` vypocet tau z parametru do jedne funkce
bez kovariance a pouzit ji v `_arc` i pri razeni.

Skupina, v niz nektery blok nema tau (tabulka nepokryva, nebo ho ochrany
odmitnou: `(R|Wo)`, `(R|DQ)`, vnorene podobvody, n <= 0), **zustava beze
zmeny**. Rozsireni pozdeji, az bude potreba: tau jako frekvence maxima
-Im(Z_bloku) na siroke logaritmicke mrizce.

Shodna tau (na relativni toleranci) neprehazovat - stabilni razeni.

### 3. Kde to udelat (jedno misto)

Nova funkce, napr. `canonical_block_order(circuit, params_full,
fixed_params, lower, upper) -> permutation` (index plneho vektoru
parametru), v `fitting/` - nejspis primo v `time_constants.py`, kde je
detekce oblouku a retezcu. Permutace volnych indexu plyne z plne, protoze
fixovane masky bloku jsou shodne.

**`fit_equivalent_circuit`** (`circuit.py`): hned po `least_squares`,
pred krokem 4 (bound status). Preradit `params_opt` a **sloupce
`opt_result.jac`** (ve volnem prostoru). `least_squares_normalized` vraci
`jac` uz v puvodnich promennych, takze staci permutace sloupcu. Rezidua se
nemeni (stejna impedance). Varovani z kroku 3 (`diag_warnings`: orezany
pocatecni odhad, konvergence) se tykaji vstupu nebo optimalizatoru, labely
v nich zustavaji spravne.

**`fit_circuit_diffevo`** (`diffevo.py`): hned po vyberu `chosen`
(blok `if chosen is not None: ... else: ...`, ~r. 505-513), **pred**
`selection_warnings` (~r. 522) - ta tiskne labely a hodnoty alternativy.
Preradit:

- `params_opt_free` i `params_opt`: analyticky Jakobian i
  `final_residuals` se v kroku 3 pocitaji z `params_opt_free`;
- sloupce `final_ls.jac` (numericky Jakobian, ve volnem prostoru);
- `sel.alternatives` (cele plne vektory parametru) stejnou funkci, jinak
  se alternativa ve varovani lisi jen stitky.

Detekce dvojznacnosti v `de_archive.py` porovnava spektra, ne parametry,
takze ji prohozene bloky nematou (overeno: "Permuted series blocks give an
identical spectrum and never pass (b)").

Multistart vola `fit_equivalent_circuit`, takze dostane kanonicke poradi
automaticky. Parametry startu mezi sebou nesrovnava (vybira podle RSS),
jiny prinos tam neni. Voigt (`cli/handlers/voigt_fit.py`) je linearni fit
bez `fit_equivalent_circuit`, razeni se tam nevola; K lezi na tau mrizce
uz vzestupne (overeno).

Vse, co se pocita az potom - `bound_status`, `params_at_bounds`,
`bounds_warnings`, `param_warnings`, kovariance, CI, signifikance,
`circuit.update_params`, `arc_time_constants`, `arc_resistance_sums` (soucet
i korelace jsou vuci prohozeni symetricke) - pak vyjde v kanonickem poradi
samo, bez permutovani matic a prepisovani textu varovani.

### 4. Co se nemeni

- `DiffEvoResult.diagnostics.initial_guess` zustava v poradi zadani
  (je to vstup uzivatele).
- Vyraz obvodu (`str(result.circuit)`) zustava strukturne stejny; meni se
  jen hodnoty v blocich.
- Vychozi `seed=None` zustava.

### 5. Hlaseni

Knihovna: `logger.debug`, zadny novy field ve `FitResult`.

CLI: jedna radka poznamky, kdyz se poradi bloku v pocatecnim odhadu lisi
od kanonickeho. U DE tiskne CLI pocatecni odhad hned nad vysledkem
(`cli/handlers/fit_output.py`, "Initial guess:"); u LM ho netiskne, ale
odhad jsou hodnoty v `--circuit`, ktere uzivatel sam zapsal, a protoze LM
jejich poradi dnes drzi, bez poznamky by prohozeni vypadalo jako chyba.
Handler na to field nepotrebuje: zavola `canonical_block_order` na
pocatecni odhad (DE: `diag.initial_guess`; LM: parametry obvodu pred
fitem) a meze fitu, a kdyz nevrati identitu, vypise napr. "Blocks
reported in ascending tau, not in the order of the initial guess".

## Testy

Test budget: zadne nove pomale testy. Jeden novy soubor,
`tests/test_block_order.py`:

1. **Cista funkce, parametrizovane**: `R-(R|Q)-(R|Q)` v "pomalem" poradi
   -> serazeno podle tau a impedance pred a po se shoduje na strojovou
   presnost; jeden blok s fixovanym R -> beze zmeny; ruzne meze -> beze
   zmeny; `(R|Q)` a `(R|C)` -> beze zmeny; `(R|Q)` a `(R|Wo)` -> beze
   zmeny; `(R|Q) - (R - (R|Q))` -> serazeno.
2. **Konzistence kovariance (LM, rychle)**: `fit_equivalent_circuit`
   s prohozenym a se spravnym odhadem -> stejne `params_opt`,
   `params_stderr` i `cov` v ramci tolerance fitu. Permutace sloupcu
   Jakobianu je nejrizikovejsi misto zmeny.

DE regresni test se seedy (puvodni bod 6) vynechan: ~1 min, zavisly na
seedech, ktere se pri kazde zmene DE posunou, a determinismus plyne
z konstrukce.

Pak odstranit `_canonical()` z `tests/test_zscope_benchmark.py`: pravda
`R(10)-(R(200)|C(1e-5))-(R(100)|C(1e-3))` ma tau 2e-3 s a 0.1 s, je tedy
uz vzestupne. Spustit `-m slow` a overit, ze projde bez nej.

**Stress test** (`tests/stress.py --check`) pred vydanim. Ocekavani: F3
stoupne z ~0.903 na ~0.918 (emulace v 97c86b8; generator zapisuje pravdu
vzestupne v 860/860 pripadech, F3 tedy nepotrebuje zmenu); F1 (start
v pravde) se nezmeni. Posun F3 je akceptacni kriterium; pak
`--update-baseline` a upravit znamy limit 9 v `doc/STRESS_TEST.md`.
(Invariance vuci permutaci ve stress testu je permutace frekvencnich bodu,
s poradim bloku nesouvisi.)

## Dokumentace a vydani

- CHANGELOG: zmena chovani - poradi zamenitelnych bloku ve vysledku je
  nove deterministicke (rostouci tau), u DE i LM. U LM dosud vysledek drzel
  poradi pocatecniho odhadu; kdo zapisoval pomaly oblouk prvni nebo
  porovnava stare a nove vysledky podle indexu, uvidi prohozene hodnoty.
- `doc/PYTHON_API.md` a README (sekce fitovani): konvence a ze jde o
  poradi ve spektru, ne identitu procesu (krizeni tau v serii -> porovnavat
  podle C). Upravit vetu z 2026-10-16 o porovnavani podle tau.
- `doc/ROBUST_LOSS_SOFT_L1.md`, sekce 6: "planovane kanonicke poradi" ->
  zavedene, s namerenou hodnotou F3.
- Verze: minor bump (zmena vystupu), pred bumpem `git fetch`, tag az po
  dotazu.
- Zavrit odkaz v `doc/EIS_FITTING_PROTOCOL_COMPARISON.md`, bod 3.

## Rozhodnuto (2026-10-20)

1. Strukturalni rozpoznani zamenitelnosti je v poradku: struktura
   rozhoduje jen o tom, co lze prohodit, poradi urcuje fyzikalni tau.
2. Prvni verze: `(R|C)`, `(R|Q)`, G misto R (z `time_constants.py`
   zadarmo) a K; ostatni beze zmeny.
3. Hlaseni: `logger.debug` v knihovne, poznamka v CLI odvozena
   v handleru (oddil 5), bez noveho fieldu.
4. Klic: vzdy rostouci tau (varianta a), i za cenu zmeny chovani LM.
