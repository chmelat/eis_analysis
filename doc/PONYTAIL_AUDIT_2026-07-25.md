# Ponytail Audit 2026-07-25 — over-engineering audit celého repa

Datum: 2026-07-25
Rozsah: celý strom projektu (61 souborů v `eis_analysis/`, `tests/`, `example/`,
`eis.py`, root). Pouze složitost a over-engineering; korektnost, bezpečnost
a výkon jsou mimo rozsah auditu.

Navazuje na `doc/PONYTAIL_AUDIT.md` (2026-07-19), jehož nálezy #1-#6 byly
aplikovány. Tento audit je nezávislý průchod se zaměřením na mrtvé cesty kódu,
vlečené příznaky a vrstvy s jediným volajícím.

**Stav:** mechanická skupina (#1, #3, #8-#11) aplikována 2026-07-25.
Skupina měnící veřejné signatury (#2, #4, #5, #7) zůstává otevřená.

## Celkové zhodnocení

Repo zůstává na svůj rozsah štíhlé — závislosti (numpy, scipy, matplotlib) jsou
všechny nosné, obě ABC vrstvy mají reálné vícenásobné implementace, dataclassy
mají živé producenty i konzumenty a re-exportní `__init__.py` jsou tenké.

Co audit našel, je jiný vzorec než minule: **vestigiální kód po refaktoringu**.
Několik funkcí bylo v minulosti nahrazeno novější implementací, ale stará verze
i parametry, které ji zapínaly, zůstaly v kódu jako mrtvé cesty. Největší nález
(#1) je celý 224řádkový soubor, který nikdo nevolá.

## Nálezy (seřazeno od největšího řezu)

| # | Tag | Nález | Náhrada | Kde |
|---|-----|-------|---------|-----|
| 1 | delete | **[APLIKOVÁNO 2026-07-25]** `plot_rl_fit_diagnostics()` — 224řádkový soubor s nulovým počtem volání. Nahrazen funkcí `_plot_rlk_fit()` (`rlk_fit.py:423`), která se skutečně používá pro `--ri-fit`. Jediné odkazy jsou dva re-exporty v `__init__.py` a ukázka v `PYTHON_API.md` | nic; smazat soubor, oba exportní bloky a sekci v dokumentaci | `eis_analysis/visualization/diagnostics.py` |
| 2 | delete | Příznak `use_voigt_fit` vlečený čtyřmi vrstvami (`compute_drt` -> `_estimate_r_inf` -> `estimate_rinf_with_inductance`), kde je dole označen `# Kept for backward compatibility, ignored`. Jediný efekt je kosmetický label `method='voigt_fit'`. CLI navíc natvrdo posílá `use_rl_fit=False, use_voigt_fit=args.ri_fit` — dva příznaky pro jedno chování | ponechat `use_rl_fit`, předávat do něj `args.ri_fit`, `use_voigt_fit` smazat všude | `drt/core.py:128`, `drt/estimation.py:84`, `rinf_estimation/rlk_fit.py:314`, `cli/handlers/drt.py:213` |
| 3 | delete | **[APLIKOVÁNO 2026-07-25]** `OnePerLineHelpFormatter` — 37 řádků, které se nikdy nevykonají. `parse_arguments()` nastavuje `usage='eis [input] [options]'`, takže `_format_usage()` vždy skončí v časné větvi `usage is not None`. Ověřeno proti živému výstupu `--help` | `formatter_class=argparse.RawDescriptionHelpFormatter` | `eis_analysis/cli/parser.py:18` |
| 4 | delete | Pole `RinfEstimate.R_ct` / `C_nF` / `f_characteristic` jsou vždy `None` — `estimate_rinf_with_inductance()` staví diagnostický dict z explicitního seznamu polí, který tyto klíče vůbec nemá. Obě podmíněné výpisové větve jsou tím nedosažitelné | nic; zrušit 3 pole a obě větve | `drt/results.py:28-30`, `cli/handlers/drt.py:58-61`, `cli/handlers/rinf.py:67-73` |
| 5 | yagni | 20řádkový diagnostický dict znovu poskládaný z `RLKFitResult` "for backward compatibility", včetně aliasu `r_squared`/`R_squared`. Třetí prvek návratové pětice je dokumentován jako "Always None (legacy placeholder)" a všichni tři volající ho zahazují | vracet `(result, fig)`; volající čtou atributy dataclassy | `rinf_estimation/rlk_fit.py:381-420` |
| 6 | shrink | `utils/compat.py` — 32řádkový modul, z toho 26 řádků docstring, kolem 4řádkového import-shimu s jediným volajícím | try/except přímo v `auto_suggest.py:15`; smazat modul a dva exportní řádky v `utils/__init__.py` | `eis_analysis/utils/compat.py` |
| 7 | delete | Parametr `verbose` na 4 veřejných vstupních bodech fittingu, dokumentovaný jako "Ignored". Tři volající ho stále zbytečně předávají | nic; úroveň logování je skutečný ovládací prvek | `fitting/circuit.py:254`, `fitting/diffevo.py:156`, `fitting/multistart.py:166`, `rinf_estimation/rlk_fit.py:311` |
| 8 | stdlib | **[APLIKOVÁNO 2026-07-25]** `_detect_delimiter()` — 12řádková smyčka sledující maximum | `return max(',', '\t', ';', key=header_line.count)` (čárka první zachová výchozí hodnotu pro případ nulových výskytů) | `eis_analysis/io/data_loading.py:456` |
| 9 | delete | **[APLIKOVÁNO 2026-07-25]** Parametr `refit_positive` — nikdy nečten, docstring přiznává "Unused parameter for backward compatibility" | nic | `fitting/voigt_chain/fitting.py:347` |
| 10 | delete | **[APLIKOVÁNO 2026-07-25]** Konstanty `DRT_TOLERANCE`, `GAMMA_MIN_REASONABLE` — komentář přiznává "Unused (pre-existing)". Nejde o importovatelné API (nejsou v `__all__`, nula externích referencí) | nic | `drt/core.py:65-66` |
| 11 | delete | **[APLIKOVÁNO 2026-07-25]** Nepoužité parametry `weighting` v `_prepare_optimization()` a `R_inf` v `_create_visualization()` — obě funkce privátní, každá s jediným volajícím | zrušit parametr i argument na místě volání | `fitting/circuit.py:188`, `drt/plotting.py:19` |

## Ověřeno jako čisté (nejde o nálezy)

- **Závislosti** — numpy, scipy i matplotlib jsou nosné, nic k odstranění.
- **ABC vrstvy** — `CompositeCircuit` -> Series/Parallel, `CircuitElement` ->
  basic/distributed/composite. Obě mají reálné vícenásobné implementace.
- **Dataclassy** — všechny (`LoadedData`, `OptimizationSetup`,
  `WeightedGMMResult`, `DiffEvoDiagnostics`, `MultistartDiagnostics`,
  `DataSelectionResult`, ...) mají živého producenta i konzumenta.
- **Re-exportní `__init__.py`** — tenké a oprávněné.
- **Ostatní modulové konstanty** — AST sweep nenašel žádnou další mrtvou mimo #10.

## Návrh postupu aplikace

Nálezy se dělí na dvě skupiny podle rizika:

**Mechanické, chování beze změny** (#1, #3, #8, #9, #10, #11) — mrtvý soubor,
mrtvý formatter, one-liner ze stdlib, nepoužité parametry a konstanty.
Zhruba -300 řádků.

**Mění veřejné signatury** (#2, #4, #5, #7) — `use_voigt_fit`, `verbose`,
návratová pětice `estimate_rinf_with_inductance()`. Podle CHANGELOGu jsou tyto
prvky drženy záměrně kvůli zpětné kompatibilitě, takže jde o rozhodnutí, zda
přijmout breaking change. Zhruba -100 řádků.

## Bilance

**net: -400 řádků, -0 závislostí.**
