# Mutacni beh nad modulem s kontraktem

Inspirace: driftlet 2026-10-18 (`docs/reliability.md`, mutace `laws.js`
a `lawRows.js`), `doc/DRIFTLET_COMPARISON1.md` bod 8.

## Princip

Do modulu se zanese jedna mala chyba (mutant): otoceny operator, `<` misto
`<=`, jina konstanta, `min` misto `max`, negovana podminka, smazany prikaz.
Na mutantovi se pusti testy. Selze-li nektery, mutant je **zabity**:
testy tu chybu vidi. Projdou-li vsechny, mutant **prezil**. Opakuje se to
pro kazde takove misto v modulu, vzdy na ciste kopii.

Mutacni beh nemeri kod, ale testy: odpovida na otazku, zda by sada
zachytila chybu, kterou nekdo pri pristi zmene udela.

## Smysl vedle ostatnich kontrol

| Kontrola | Porovnava |
|---|---|
| Cold review (`doc/CONTRACT_REVIEW.md`) | kod s kontraktem |
| Testy | kod s ocekavanim v nekolika bodech |
| **Mutacni beh** | testy s kodem: hlidaji testy kazdou klauzuli kontraktu? |

Pravidlo z `CLAUDE.md` (modul s kontraktem se meni s kontraktem a testem)
zajisti test ke kazde zmene, ne ke kazde klauzuli. Po review vznikly testy
hlavne k nalezenym chybam; klauzule, ktere platily uz predtim, test mit
nemusi a drzi jen tim, ze na kod nikdo nesahl. Mutacni beh je najde.

Neumi najit chybu v kodu: je-li kod spatne a test tu spatnost zabetonuje,
mutant je zabity. Proto az po cold review a jen u modulu s kontraktem.

## Trideni prezivsich

1. **Ekvivalentni:** zadny test ho zachytit nemuze (`>` misto `>=` na
   float prahu, na kterem hodnota nelezi; zmena jen v logu). Zapsat, dal
   neresit.
2. **Prah nebo tolerance:** test nepinuje presnou hodnotu konstanty.
   Test jen tam, kde na hodnote neco zavisi.
3. **Dira:** klauzule kontraktu bez testu. Novy test napsany proti
   klauzuli, ne proti tomu, co kod dela.

Prezivsi se zkousi ve dvou stupnich: nejdriv testy modulu (rychle), pak
prezivsi proti cele vychozi sade. Mutant, ktereho zabije az vzdaleny
test, je viden; dira je az to, co projde vsim.

## Metoda

Skript `mutate.py` (zatim scratchpad sezeni 1669e97b, `mut/`): mutanty
z tokenu (operatory porovnani a aritmetiky, `and`/`or`, `True`/`False`,
`min`/`max`, `argmin`/`argmax`, `any`/`all`, cisla: 0 <-> 1, n -> n+1,
x -> 10x) a z AST (negace podminky `if`/`while`/`x if c else y`, smazani
prikazu-vyrazu). Mutant, ktery se neprelozi nebo dava stejny AST, se
zahodi. Kazdy proces ma vlastni kopii repozitare (bez `tmp/`, `.git`),
import se overi, ze vede do ni, a bytecode cache je vypnuta
(`PYTHONDONTWRITEBYTECODE`), aby mutant stejne velikosti nebezel ze stare
`.pyc`. Nemutovana sada musi projit. 4 procesy, `pytest -x`.

## Kolo 1: fitting/covariance.py (2026-10-25)

91 mutantu (40 operatoru, 28 cisel, 11 jmen, 11 negaci, 1 smazani).

**Stupen 1** (`tests/test_covariance.py`, `tests/test_confidence_intervals.py`,
26 testu, 3.3 s): 55 zabito, 36 prezilo, 0 timeout.

| Skupina | Mutanty (id, radek) | Pocet |
|---|---|---|
| 1 ekvivalentni | 15 (L248 mez nuloveho sloupce 1.0 -> 10.0), 19 (L253 `full_matrices=True`), 25, 26 (L259 jen `S[-1] = 0`, dava inf i tak), 33, 36, 38, 71 (strict/nestrict na float prahu), 73 (L452 hodnota zahozena `where`), 90 (L346 jen log) | 10 |
| 2 prah | 0 (`NULLSPACE_OVERLAP_TOL` 1e-3 -> 1e-2), 1 (`rcond` 1e-10 -> 1e-9), 30, 31 (L263 prah `rcond * S[0]` -> `/`, `S[1]`) | 4 |
| 3 dira | nize | 22 |

Diry ve stupni 1 (klauzule kontraktu, kterou testy modulu nehlidaji):

- **D1 hodnota kovariance pri plne hodnosti** ("cov is s^2 (J^T J)^-1"):
  47, 48 (L331 `1/(S*S)` -> `S*S`, `1`), 49 (L332 preskalovani). Testy
  pouzivaji J s ortogonalnimi sloupci, kde S = 1, a kontroluji jen
  konecnost; s referenci `s^2 inv(J^T J)` se nesrovnava.
- **D2 hodnoty identifikovatelnych parametru pri nedostatecne hodnosti**
  ("the others get their variance from the pseudo-inverse"): 40-44
  (L296-298, kazda aritmetika vetve). Testy kontroluji jen, ktere stderr
  jsou inf.
- **D3 cislo podmineni a nezavislost na jednotkach**: 14 (L248 skaluj jen
  sloupce s normou > 1), 22, 24 (L259 `S[-1]` -> `S[+1]`). Testy maji jen
  2 volne parametry (tam S[1] = S[-1]) a zadny nepreskaluje sloupec.
- **D4 varovani pri spatnem podmineni**: 84 (L321 negace). Zadny test
  nekontroluje text ani nepritomnost varovani.
- **D5 pole degenerovanych vysledku**: 9 (L220 `rank=0` pri dof <= 0),
  52, 53 (L370 `_failed`: rank, `is_well_conditioned`), 87 (L368
  `_failed` s `fixed_params`: delka stderr).
- **D6 okrajove pripady**: 13 (L248 nulovy sloupec, `> 0` -> `>= 0`),
  32 (L264 cely J nulovy), 56 (L426 CI s dof <= 0), 66 (L440 CI pro p = 0
  pod maskou), 77, 82 (L460 CI pro NaN stderr).

**Stupen 2** (36 prezivsich proti cele vychozi sade, 901 testu, 108 s na
beh, 3 procesy): 7 zabito, 29 prezilo.

| Mutant | Zabil ho |
|---|---|
| 14 (D3, skaluj jen normu > 1) | `test_rinf_estimation.py::test_rinf_does_not_depend_on_units` |
| 22, 24 (D3, `S[+1]`), 31 (prah z `S[1]`) | `test_bounds_diagnostics.py::test_bounds_warnings_consistent_with_bound_status` |
| 47, 48 (D1, L331) | `test_drt_recovery.py::test_default_r_inf_warning_joins_the_drt_warnings` |
| 49 (D1, L332) | `test_cli_integration.py::test_rinf_estimation` |

D3 je tedy videna. D1 je videna jen nahodou: zabily ji testy textu
varovani v DRT a CLI vystupu R_inf, ne test kovariance, a pri zmene techto
testu by zase byla slepa.

**Vysledek kola 1** (29 prezivsich cele sady):

| Skupina | Mutanty | Pocet |
|---|---|---|
| ekvivalentni | 15, 19, 25, 26, 33, 36, 38, 71, 73, 90 | 10 |
| prah | 0, 1, 30 | 3 |
| dira D2: hodnoty pri nedostatecne hodnosti | 40, 41, 42, 43, 44 | 5 |
| dira D4: varovani pri spatnem podmineni | 84 | 1 |
| dira D5: pole degenerovanych vysledku | 9, 52, 53, 87 | 4 |
| dira D6: okrajove pripady | 13, 32, 56, 66, 77, 82 | 6 |

Skore: 62 z 91 zabito (68 %), z neekvivalentnich 62 z 81 (77 %).

Zaver:

- **D2 je skutecna dira.** Vetev nedostatecne hodnosti hlasi uzivateli
  stderr identifikovatelnych parametru a kazda jeji aritmetika muze byt
  spatne (i `residual_variance / JtJ_pinv`), aniz by to cokoli v sade
  videlo. Existujici `test_partial_rank_deficiency_per_parameter` ji
  nevidel, protoze identifikovatelny sloupec je jednotkovy a kolmy na
  ostatni (S = 1, skala 1: kazda mutace da stejne cislo). **Uzavreno
  2026-10-25:** `test_partial_rank_deficiency_identifiable_variance_matches_pinv`
  (sloupce ruzne normy, c0 korelovany s c1, c2 = 5 c1; reference
  `s^2 pinv(J^T J)`, presna pro parametr kolmy na null space) zabije vsech
  5 mutantu 40-44.
- **D1 dostat vlastni test**, i kdyz je videna: `s^2 inv(J^T J)` pro J
  s neortogonalnimi sloupci ruzne normy, proti primemu vypoctu.
- D4, D5, D6 jsou male a levne: text varovani pri cond > 1e10 a jeho
  nepritomnost pri dobrem podmineni; `rank`, `is_well_conditioned`
  a delka stderr v `_failed` a pri dof <= 0 (i s `fixed_params`); nulovy
  sloupec a cely nulovy J; CI s dof <= 0, p = 0 pod maskou, NaN stderr.
  Vsechny jsou klauzule kontraktu, takze testy patri k nemu.
- Prahy (`NULLSPACE_OVERLAP_TOL`, `rcond`) netestovat na presnou hodnotu:
  overlap je v praxi ~1 nebo ~eps a konstanta ma zdokumentovany duvod.

**Poucky pro dalsi kola:**

- Plna sada spousti vlastni vlakna (multistart 4, DE `workers=2`) a BLAS
  pouziva vsechna jadra: 4 procesy zaplnily 8 jader. Skript nastavuje
  `OMP/OPENBLAS/MKL_NUM_THREADS=1`; pak 2 procesy vytizi ~2.5 jadra
  (spicky 3.3), 3 procesy jsou v poradku.
- Stupen 1 je levny (91 mutantu za ~1.5 min), stupen 2 na prezivsich
  ~25 min. Bez stupne 2 by D1 a D3 vypadaly jako diry.
