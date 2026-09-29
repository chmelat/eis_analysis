# Audit modulu `cli/parser.py`

**Datum:** 2026-07-08 | **Verze:** 0.16.22 | **Typ:** kritické code-level review (Claude Code)

## Stav oprav (aktualizováno 2026-07-08)

| Nález | Stav | Verze | Commit |
|-------|------|-------|--------|
| P1 | Opraveno | 0.16.23 | — |
| P2-P7 | Otevřené | — | — |

**Rozsah:** `cli/parser.py` (265 ř.), konzumace argumentů v `eis.py`,
`cli/logging.py`, `cli/data_handling.py`, `cli/handlers/*`; konzistence
s README a s test fixture `create_test_args`.

## Verdikt

Struktura je čistá: logické skupiny, choices pro enumy, `dest='lambda_reg'`
řeší kolizi s klíčovým slovem, help texty jsou po K2 auditu věcně správné
(mu jako stop hodnota). Rizika leží v tichém ignorování `--multistart`,
mrtvém formatteru a chybějící validaci číselných rozsahů.

## Nálezy

### P1 — `--multistart N` bez `--optimizer multistart` se tiše ignoruje (Střední)

Výchozí optimizer je 'de'. `eis data --circuit ... --multistart 4` spustí
Differential Evolution a N se zahodí bez varování (ověřeno živě). Navíc
si protiřečí dokumentace: help říká „default: 0 = disabled", ale
v multistart režimu 0 znamená „použij 16" (`handlers/fitting.py` ř. 396,
skrytá konstanta; README ř. 316 dokumentuje 16 správně). Záporné N
→ taky 16, bez varování.
**Náprava:** `--multistart N > 0` buď implikuje `--optimizer multistart`,
nebo aspoň warning „--multistart ignored with --optimizer de"; sjednotit
help s README (default 16); vyhodit sémantiku „0 = disabled".

### P2 — Mrtvé jádro `OnePerLineHelpFormatter` (Střední)

Parser předává `usage='eis [input] [options]'`, takže `_format_usage`
vždy skončí na early-return pro custom usage (ř. 26-28) — celá smyčka
„one option per line" (ř. 30-54) je nedosažitelná a docstring třídy
(„puts each option on a separate line") nepopisuje skutečné chování
(ověřeno: `--help` tiskne jednořádkové usage). Formatter tak reálně
slouží jen jako `RawDescriptionHelpFormatter` pro epilog.
**Náprava:** smazat mrtvou větev a podědit jen
`RawDescriptionHelpFormatter`, nebo odstranit custom `usage=` a mrtvý
kód tím oživit. Drobnost: literál `%` v usage stringu by shodil
`usage % dict(...)`.

### P3 — Žádná validace číselných rozsahů (Střední)

`--n-tau -5` → syrový numpy traceback z `linspace` hluboko v DRT
(ověřeno; exit 1, ale bez čistého „!!" formátu chyb, který CLI jinak
používá). Stejně neošetřené: `--n-tau 0`, `--de-popsize 0`,
`--de-maxiter 0`, záporné `--voigt-n-per-decade`, `--voigt-max-M`,
`--extend-decades-max`, `--gmm-bic-threshold`. Projekt přitom deklaruje
edge-case hygienu (CLAUDE.md). Souvisí: `main()` v `eis.py` chytá jen
`EISAnalysisError` a KeyboardInterrupt — cokoli jiného je traceback.
Pozn.: `--f-min > --f-max` končí čistou chybou jen proto, že filtr
vyprázdní data („No data remaining").
**Náprava:** malý `type=`-helper (positive_int/positive_float) na
dotčených volbách; případně catch-all v `main()` s odkazem na -v.

### P4 — `--area 1.0` nelze odlišit od defaultu (Nízká)

Default `1.0` v parseru + `if args.area == 1.0: # Default value was not
changed` v `handlers/oxide.py` ř. 48: uživatel s elektrodou přesně
1 cm², který zadá `--area 1.0` explicitně, dostane místo své hodnoty
area z DTA metadat. Už zmíněno okrajově v oxide auditu; správné místo
opravy je parser: `default=None`, default 1.0 doplnit v handleru.

### P5 — `parse_arguments()` nelze testovat; defaulty zdvojené v testech (Nízká)

Funkce nemá parametr `argv` a volá `parser.parse_args()` napřímo —
integrační testy proto ručně duplikují všech ~40 defaultů
(`create_test_args` v `test_cli_integration.py`). Změna defaultu
v parseru se do testů nepropíše a nikdo si toho nevšimne (drift už
dnes: help `--multistart` vs. README, viz P1).
**Náprava:** `def parse_arguments(argv=None)` → `parse_args(argv)`;
testy pak mohou stavět Namespace přes
`parse_arguments(['--voigt-chain', ...])`.

### P6 — `-v` je count, ale čte se jen jako bool (Kosmetika)

`action='count'` slibuje úrovně (-vv), `logging.py` ř. 88 testuje jen
`>= 1`. Kombinace `--quiet -v` je povolená a dává zvláštní výstup
(debug bez info) — nezdokumentováno.

### P7 — Epilog příklad zavádí (Kosmetika)

„`eis --ri-fit data.csv` — Analyze CSV data file": --ri-fit je robustní
odhad R_inf, s analýzou CSV nesouvisí; příklad směšuje dvě věci.

## V pořádku (ověřeno)

Skupiny voleb odpovídají handlerům; choices kryjí všechny enumy
(weighting, optimizer, peak-method, formát, DE strategie, voigt-fit-type);
`--lambda` přes `dest` bez kolize s keywordem; `BooleanOptionalAction`
kompatibilní s requires-python >= 3.9; `--ocv` má konzumenta
(`eis.py` ř. 96); mu help texty konzistentní s K2 opravou (obě volby);
`--f-min/--f-max` s prázdným výsledkem končí čistou EISAnalysisError;
defaulty parseru aktuálně sedí s `create_test_args` (namátkou ověřeno).

## Priority

1. P1 — tiché ignorování `--multistart` (přímý dopad na výsledky:
   uživatel dostane jiný optimizer, než si myslí).
2. P3 — validace rozsahů (UX: traceback místo chyby).
3. P2 + P5 — mrtvý kód a testovatelnost (levné, zvyšuje udržovatelnost).
4. P4, P6, P7 — průběžně.
