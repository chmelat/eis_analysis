# Audit diagnostiky fittingu

**Datum:** 2026-07-03 | **Verze:** 0.16.17 | **Typ:** kritické code-level review (Claude Code)

## Stav oprav (aktualizováno 2026-07-08)

| Nález | Stav | Verze | Commit |
|-------|------|-------|--------|
| D1, D2, D5 | Opraveno | 0.16.18 | 3b41082 |
| D3 | Opraveno | 0.16.21 | — |
| D4, D6-D10 | Otevřené | — | — |

**Rozsah:** `fitting/diagnostics.py` (váhy, metriky), `fitting/covariance.py`
(kovariance, CI), `fitting/bounds.py`, diagnostické části `fitting/circuit.py`
(FitDiagnostics, FitResult, kroky 4-8), konzumace v `cli/handlers/fitting.py`
a v multistart/diffevo.

## Verdikt

Statistické jádro je správné: s² = RSS/dof s váženými rezidui (škálově
invariantní vůči normalizaci vah), dof = 2N − p, SVD kovariance
s rank-deficient → inf (curve_fit konvence), cond(JᵀJ) = cond(J)²,
t-distribuce CI s konzistentním dof. Rizika leží v protichůdných a špatně
indexovaných diagnostických zprávách, tichém oříznutí initial guess
a zavádějícím stderr = 0 u Voigt chain.

## Nálezy

### D1 — Dvě protichůdná kritéria „parametr na mezi" (Střední)

Krok 4 (`circuit.py`, ř. 414-422) plní `FitDiagnostics.bounds_warnings`
kritériem *do 1 % od hodnoty meze* (`|p−b|/max(|b|,1e-10) < 0.01`);
`bound_status` (krok 8b) používá `classify_bound_status` (`bounds.py`,
ř. 107) s kritériem *1 dekáda od meze v logu* (pro meze > 6 dekád).
Liší se řádově: R = 5e-4 u dolní meze 1e-4 je podle `bound_status`
„lower" (CLI potlačí CI: „CI not meaningful"), ale `bounds_warnings`
mlčí. `_log_fit_result` vypisuje obojí → výstup si protiřečí.
**Náprava:** jediné kritérium (classify_bound_status) pro oba kanály;
krok 4 z něj odvodit.

### D2 — Indexy parametrů ve varováních nesedí při fixních parametrech (Střední)

`bounds_warnings` („Parameter {i} at lower bound", `circuit.py` ř. 419)
indexuje ve *free-space* (bez fixních), `param_warnings` (ř. 438) a labels
ve *full-space*. U `R("100") - (R()|C())` znamená „Parameter 0" v jedné
zprávě R0, v druhé R1. **Náprava:** všude full-space indexy, ideálně
s labely místo čísel.

### D3 — Tiché oříznutí initial guess do bounds (Střední)

`_prepare_optimization` (`circuit.py`, ř. 216-229) ořízne odhad mimo meze
a zaznamená `clipped_params` — ale ten se nikam nepropaguje (není ve
FitDiagnostics, žádný log). Fit tiše startuje z jiného bodu, než uživatel
zadal (např. R(0) → 1e-4). Stejný typ tichého předpokladu jako oxide O3.
**Náprava:** propsat do FitDiagnostics.warnings + log.

### D4 — Voigt chain: stderr = 0 (Střední)

`cli/handlers/fitting.py` ř. 302: `params_stderr=np.zeros(...)` znamená
„známo přesně"; `params_ci_95` vrátí CI ±0. Lineární fit odhad nejistoty
nemá — poctivá hodnota je inf/NaN (konvence dodržovaná jinde).

### D5 — Mrtvý kód `check_bounds_proximity` (Nízká)

`bounds.py` ř. 139: exportovaná, nikým nevolaná — třetí nezávislá
implementace téže kontroly (vedle kroku 4 a classify_bound_status).

### D6 — Hrubá CI při částečné neidentifikovatelnosti (Nízká)

`params_ci_95/99` (`circuit.py` ř. 120): jediný inf stderr → CI
(−inf, inf) pro *všechny* parametry, včetně dobře určených a fixních
(stderr = 0). CLI to maskuje tagem [fixed], API uživatel dostane nesmysl.

### D7 — Nekonzistentní kontrakt pro neznámý weighting (Nízká)

`fit_equivalent_circuit` raisuje ValueError; přímé volání
`compute_weights` (`diagnostics.py` ř. 49) jen warninguje a tiše použije
uniform.

### D8 — Nedokumentované magic konstanty (Nízká)

- 0.01 a 1e-10 v kroku 4 (`circuit.py` ř. 417),
- práh 10 p.b. pro poznámku weighted vs. unweighted (`diagnostics.py`
  ř. 99) — absolutní rozdíl navíc nezachytí řádový nepoměr u malých chyb
  (0.5 % vs. 8 % projde bez poznámky),
- `s > 2|p|` (`circuit.py` ř. 439),
- „acceptable" = `FIT_QUALITY_GOOD_ERROR * 2` — odvozený práh bez
  vlastního jména, duplikovaný v CLI (`fitting.py` ř. 164).

### D9 — Deklarace „No logging in core functions" neplatí (Nízká)

`circuit.py` ř. 4 tvrdí „all diagnostics returned as data" —
`compute_fit_metrics` loguje info, `compute_weights` warning,
`covariance.py` error. Opravit kód, nebo docstring.

### D10 — Edge case prázdná data (Nízká)

`compute_weights` dělí průměrem prázdného pole (RuntimeWarning, NaN);
projekt jinak edge-cases hlídá (CLAUDE.md).

### Kosmetika

- `rank` jako np.int64 v plné cestě, `int()` jen v rank-deficient větvi
  (`covariance.py` ř. 182 vs. 233).
- Lokální `warnings` v `all_warnings` stíní modul warnings
  (`circuit.py` ř. 146).

## V pořádku (ověřeno)

Vážená relativní chyba bez dvojího započtení váhy (pro modulus =
mean(|dZ|/|Z|)); multistart i diffevo používají tutéž metriku
z `compute_fit_metrics`; výběr DE-vs-refined podle cost (ne zobrazované
chyby) je zdokumentovaný záměr; potlačení CI u parametrů na mezi v CLI
je správně (Jacobiánové CI na mezi neplatí); testy diagnostiky existují
(covariance, CI, fit metrics, correlation, jacobian, multistart).

## Priority

1. D1 + D2 — sjednotit kritérium a indexování (přímý dopad na výstup).
2. D3 — zviditelnit oříznutí initial guess.
3. D4 — poctivé stderr u Voigt chain.
4. D5-D10 + kosmetika; testy k opravám průběžně.
