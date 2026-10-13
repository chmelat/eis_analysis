# Plan: kanonicke poradi zamenitelnych bloku po fitu

Stav: **ve fronte (2026-10-13), neimplementovano**. Vychozi verze
eis_analysis v0.58.3. Podnet: revize eis-fitting-protocol
(`doc/EIS_FITTING_PROTOCOL_COMPARISON.md`, nalez 1 a bod 3 "Co si odtud
vzit").

## Problem

Seriove zapojene bloky stejne struktury, napr. dva `(R|Q)` v
`R - (R|Q) - (R|Q)`, jsou zamenitelne: prohozeni jejich hodnot nemeni
impedanci. Fit tak muze vratit tytez oblouky v libovolnem poradi a ktery
z nich se jmenuje `R1/Q1` a ktery `R2/Q2`, je nahoda.

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

Pro jeden fit to nevadi. Vadi to pri srovnani fitu rady spekter (oxidove
serie, export do CSV/JSON a dalsi zpracovani): "R1" v jednom souboru muze
byt jiny fyzikalni proces nez "R1" ve druhem.

Uz se to jednou projevilo: `tests/test_zscope_benchmark.py` ma vlastni
`_canonical()`, ktery vetve `Two_TimeConstants` radi podle tau pred
srovnanim, protoze fit je na ZScope 2 % datech vratil prohozene
(CHANGELOG v0.25.5).

## Navrh

Po optimalizaci, **pred** vypoctem kovariance, stavu mezi a varovani,
preradit hodnoty zamenitelnych bloku tak, aby tau rostlo zleva doprava
(nejrychlejsi proces prvni = nejvyssi frekvence, jak se cte Nyquist).
Stejne poradi uz pouziva `auto_suggest.py` (piky DRT razene podle
`peak_tau`), takze navrzeny obvod a jeho fit budou cist bloky stejne.

### 1. Ktere bloky jsou zamenitelne

V kazdem uzlu `Series` (rekurzivne, nejen na nejvyssi urovni) najit
skupiny primych potomku, ktere:

- maji **stejnou strukturu**: stejny strom typu prvku a zapojeni, bez
  ohledu na hodnoty (porovnani napr. pres `repr` s vymazanymi cisly nebo
  pres rekurzivni tuple `(typ uzlu, typy potomku)`);
- maji **stejnou masku fixovanych parametru** a u fixovanych **stejne
  hodnoty**. Blok s fixovanym `R = 100` a blok s volnym R zamenitelne
  nejsou - prohozeni by zmenilo model.

Meze se generuji z labelu (`generate_simple_bounds(param_labels)`), takze
bloky se stejnou strukturou maji stejne meze. Pri implementaci overit, ze
to plati i pro meze odvozene jinak (napr. `log_search_bounds`).

Poznamka k pravidlu "fyzikalni kriterium, ne strukturalni": rozpoznani
symetrie je nutne strukturalni (symetrie je vlastnost vyrazu), ale klic
razeni je fyzikalni (tau). Strukturalni cast nerozhoduje o fyzikalnim
vyznamu prvku, jen o tom, co lze prohodit beze zmeny modelu.

### 2. Klic razeni: tau bloku

Prvni verze jen analyticky, pro bloky, kde je tau jednoznacne:

| Blok | tau |
|---|---|
| `(R\|C)` | `R*C` |
| `(R\|Q)` | `(R*Q)^(1/n)` (maximum -Im v `omega = 1/tau`) |
| prvek K | jeho vlastni `tau` |

Skupina, v niz nektery blok nema tau podle tabulky (napr. `(R|Wo)`,
`(R|DQ)`, vnorene podobvody), **zustava beze zmeny**. Rozsireni pozdeji,
az bude potreba: tau jako frekvence maxima -Im(Z_bloku) na siroke
logaritmicke mrizce - obecne a fyzikalni, ale drazsi a s okrajovymi pripady
(blok bez maxima).

Shodna tau (na relativni toleranci) neprehazovat - stabilni razeni.

### 3. Kde to udelat (jedno misto)

Nova funkce, napr. `canonical_block_order(circuit, params_full,
fixed_params) -> permutation` (index plneho vektoru parametru), v
`fitting/` (maly modul nebo `circuit_builder.py`). Vola se ve dvou cestach,
ktere stavi `FitResult`:

- `fit_equivalent_circuit` (`circuit.py`), hned po `least_squares`,
  pred krokem 5 (kovariance). Prepermutovat `params_opt` a **sloupce
  `opt_result.jac`** (ve volnem prostoru; permutace volnych indexu plyne
  z plne permutace, protoze fixovane masky bloku jsou shodne). Rezidua se
  nemeni (stejna impedance).
- `fit_circuit_diffevo` (`diffevo.py`), po vyberu `params_opt` a pred
  krokem 3 (kovariance). Analyticky Jakobian se pocita az v bode, takze se
  spocita z prerazenych parametru; numericky `final_ls.jac` prepermutovat
  po sloupcich.

Multistart vola `fit_equivalent_circuit`, takze dostane kanonicke poradi
automaticky (a jeho srovnani startu bude konzistentni). Cesta v
`cli/handlers/fitting.py:677` (Voigt) zamenitelne bloky v tomto smyslu
nema - overit a pripadne vynechat.

Vse, co se pocita az potom - `bound_status`, `params_at_bounds`,
`bounds_warnings`, `param_warnings`, kovariance, CI, signifikance,
`circuit.update_params` - pak vyjde v kanonickem poradi samo, bez
permutovani matic a prepisovani textu varovani.

Pred implementaci zkontrolovat varovani vznikla **pred** permutaci
(`diag_warnings` v `diffevo.py`, archiv DE v `de_archive.py`): pokud
obsahuji indexovane labely (`R1`, `Q2`), byla by po prerazeni zavadejici.
`ambiguous_alternatives` nese cele vektory parametru alternativ - ty
prerazit stejnou funkci, jinak se alternativa lisi jen stitky.

### 4. Co se nemeni

- `DiffEvoResult.diagnostics.initial_guess` zustava v poradi zadani
  (je to vstup uzivatele). CLI ho tiskne vedle vysledku - zvazit, zda
  pridat poznamku, kdyz se poradi zmenilo.
- Vyraz obvodu (`str(result.circuit)`) zustava strukturne stejny; meni se
  jen hodnoty v blocich.
- Vychozi `seed=None` zustava.

### 5. Hlaseni

Ze se poradi zmenilo, je informace o behu, ne o vysledku. Navrh:
`logger.debug` v knihovne; zadny novy field ve `FitResult`, dokud ho nikdo
nepotrebuje (YAGNI). Pokud ano, pak podle `doc/CLI_OUTPUT_UNIFICATION.md`
jako field vysledku, ne `logger.info`.

## Testy

Jeden novy soubor, napr. `tests/test_block_order.py`:

1. **Jednotkove**: obvod `R-(R|Q)-(R|Q)` s parametry v "pomalem" poradi ->
   permutace je seradi podle tau; impedance pred a po se shoduje na
   strojovou presnost.
2. **Fixovany parametr**: jeden blok s fixovanym R -> zadna permutace.
3. **Ruzna struktura**: `(R|Q)` a `(R|C)` v serii -> zadna permutace;
   `(R|Q)` a `(R|Wo)` -> zadna permutace.
4. **Vnoreni**: zamenitelne bloky uvnitr vnorene `Series` se seradi.
5. **Konzistence kovariance**: fit s prohozenym odhadem a fit se spravnym
   odhadem -> stejne `params_opt`, `params_stderr` i `cov` (v ramci
   tolerance fitu).
6. **Regrese z mereni**: spektrum B z tabulky vyse, `fit_circuit_diffevo`
   s nekolika seedy, ktere driv vracely pomalejsi oblouk prvni -> vzdy
   rychlejsi prvni. Seedy vybrat tak, aby test bezel do ~1 min.

Pak odstranit `_canonical()` z `tests/test_zscope_benchmark.py`
(benchmark je `slow`, spustit `-m slow`) a overit, ze projde bez nej.

Stress test (`tests/stress.py --check`) spustit pred vydanim: invariance
vuci jednotkam a permutaci se nesmi zhorsit.

## Dokumentace a vydani

- CHANGELOG: zmena chovani - poradi zamenitelnych bloku ve vysledku je
  nove deterministicke (rostouci tau). Kdo porovnava stare a nove vysledky
  podle indexu, uvidi prohozene hodnoty.
- `doc/PYTHON_API.md` a README (sekce fitovani): jedna veta o konvenci.
- Verze: minor bump (zmena vystupu), tag az po dotazu.
- Zavrit odkaz v `doc/EIS_FITTING_PROTOCOL_COMPARISON.md`, bod 3.

## Otevrene otazky (rozhodnout pred implementaci)

1. Je strukturalni rozpoznani zamenitelnosti (oddil 1) v poradku vuci
   pravidlu "fyzika, ne struktura"?
2. Staci prvni verze jen s `(R|C)`, `(R|Q)` a K, ostatni beze zmeny?
3. Hlasit zmenu poradi uzivateli (poznamka v CLI vedle initial guess), nebo
   jen `logger.debug`?
