# Ponytail Audit — over-engineering audit celého repa

Datum: 2026-07-19
Rozsah: celý strom projektu (eis_analysis/, eis.py, root). Pouze složitost a
over-engineering; korektnost, bezpečnost a výkon jsou mimo rozsah auditu.

## Celkové zhodnocení

Repo je na svůj rozsah štíhlé:

- žádné mrtvé funkce ani třídy (AST + grep sweep všech veřejných symbolů),
- závislosti minimální a všechny nutné (numpy, scipy, matplotlib),
- ABC vrstvy mají reálné vícenásobné implementace (CompositeCircuit -> Series,
  Parallel; CircuitElement -> basic/distributed/composite),
- re-exportní `__init__.py` soubory jsou tenké a oprávněné.

## Nálezy (seřazeno od největšího řezu)

| # | Tag | Nález | Náhrada | Kde |
|---|-----|-------|---------|-----|
| 1 | delete | **[APLIKOVÁNO 2026-07-19]** Mrtvé konstanty `DEFAULT_R0_GUESS`, `DEFAULT_Q_N_GUESS`, `DRT_PEAK_MIN_SPACING_DECADES` — 0 referencí mimo config (~45 řádků vč. docstringů) | nic; smazat i z `__all__` | `eis_analysis/fitting/config.py` |
| 2 | shrink | **[APLIKOVÁNO 2026-07-19]** 4 formatter třídy (`InfoFormatter`, `WarningFormatter`, `ErrorFormatter`, `DebugFormatter`) — stejná logika, jiný prefix | jedna třída `PrefixFormatter` s dict `{WARNING: "! ", ERROR: "!! ", DEBUG: "[DEBUG] "}` | `eis_analysis/cli/logging.py` |
| 3 | delete | **[APLIKOVÁNO 2026-07-19]** Duplicitní datové soubory v rootu: `example_eis_data.csv`, `real_gamry_example.DTA` — identické kopie verzované v `example/` | nic; příklady odkazovat na `example/` | `./` |
| 4 | delete | **[APLIKOVÁNO 2026-07-19]** `requirements.txt` — 3 řádky duplikující dependencies z `pyproject.toml` (porušení single source of truth); v README ani CI neodkazovaný | `pip install .` | `requirements.txt` |
| 5 | shrink | **[APLIKOVÁNO 2026-07-19]** `sort_by_frequency` — 45řádkový docstring na 2řádkovém wrapperu kolem `np.argsort` | zdůvodnění ascending pořadí ve 2 větách (-21 řádků) | `eis_analysis/utils/impedance.py` |
| 6 | delete | **[APLIKOVÁNO 2026-07-19]** Odkaz na neexistující `AUDIT_REPORT.md` v docstringu modulu | nic | `eis_analysis/fitting/config.py` |

## Mimo bodování

V rootu leží neverzované pracovní soubory: `README.pdf`,
`drt_clean_vs_noisy.png`, `mod1d.sh`, `mod2c.sh` a pět srovnávacích poznámek
`doc/*_COMPARISON*.md`. Nejsou v gitu, takže nejde o bloat repa — jen
rozhodnout, zda je commitovat, nebo přidat do `.gitignore`.

## Bilance

**net: -110 řádků, -0 závislostí.**

Audit je one-shot report. Všechny nálezy #1-#6 aplikovány 2026-07-19
(-68 řádků, -3 duplicitní soubory, chování beze změny).
