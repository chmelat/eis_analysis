# Code Review — EIS Analysis Toolkit

**Verze projektu:** 0.13.6 (2026-04-26)
**Konsolidováno:** 2026-04-26
**Zdrojové reporty:** `doc/archive/CODE_ANALYSIS_REPORT.md` (2026-01-09, v0.11.0),
`doc/archive/CODE_REVIEW_REPORT.md` (2026-01-10, v0.13.1),
`doc/archive/CODE_DESIGN_ANALYSIS.md` (2026-01-11, v0.13.2),
`doc/archive/CODE_REVIEW_2026.md` (2026-02-07, v0.13.3)

Tento dokument je konsolidovaným pohledem na stav projektu. Pro aktuální
audit a otevřené priority viz `doc/AUDIT_2026-04-24.md`.

---

## 1. Shrnutí

EIS Analysis Toolkit je modulární jednoautorský vědecký nástroj pro
analýzu elektrochemické impedanční spektroskopie. Architektura je čistá,
vědecké metody jsou korektně implementované, dokumentace je nadstandardní.

| Metrika | Hodnota |
|---------|---------|
| Python soubory | 49 |
| Produkční kód | 12 476 řádků |
| Testovací kód | ~2 500 řádků |
| Testy | 75 (100% passing) |
| Lint (ruff) | čistý |
| Mypy | 110 chyb (non-blocking v CI) |
| Cyklické závislosti | 0 |
| Závislosti core | 3 (numpy, scipy, matplotlib) |
| CI pipeline | ruff + pytest + mypy (Python 3.9, 3.12) |
| Licence | MIT |

**Hodnocení:** Nadprůměrný projekt. Hlavní zbývající rizika leží
v **udržitelnosti při růstu** (testovací pokrytí kritických modulů,
soulad s vlastním 500-řádkovým limitem).

---

## 2. Architektura

### 2.1 Struktura modulů

```
eis_analysis/                  12 476 řádků
  fitting/          ~5 600 (45%)  Circuit elements + 3 optimalizační strategie
  drt/              ~1 700 (14%)  DRT analýza, GCV, GMM detekce píků
  cli/              ~1 600 (13%)  CLI rozhraní, parsování, handlery
  validation/       ~1 050  (8%)  Lin-KK + Z-HIT validace
  io/                 750  (6%)  Gamry DTA + CSV parser
  rinf_estimation/    650  (5%)  Odhad R_inf s induktancí
  visualization/      440  (4%)  Nyquist, Bode, diagnostické grafy
  analysis/           395  (3%)  Doménová analýza (oxidová vrstva)
  utils/              182  (1%)  Sdílené utility
```

### 2.2 Vrstvená architektura

```
+------------------+
|     CLI Layer    |  eis.py, cli/
+------------------+
         |
         v
+------------------+
|   Core Modules   |  drt/, fitting/, validation/, analysis/
+------------------+
         |
         v
+------------------+
|   I/O + Utils    |  io/, utils/, visualization/
+------------------+
```

| Aspekt | Stav |
|---|---|
| Cyklické závislosti | 0 (acyklický graf) |
| Maximální hloubka importu | 3 |
| Cross-module coupling | Nízký |

Cyklická závislost mezi `rinf_estimation/` a `fitting/voigt_chain/` je
ošetřena lazy importem v `rinf_estimation/rlk_fit.py:26-35`.

### 2.3 Silné stránky architektury

- **Žádný side-effect v jádře.** Core funkce vrací strukturovaná data
  (dataclassy), nelogují, nepíší na stdout. CLI vrstva řídí veškerý
  uživatelský výstup. Umožňuje použití jako Python knihovny.
- **Dataclassy jako API kontrakty.** `DRTResult`, `FitResult`, `KKResult`,
  `ZHITResult` místo n-tic nebo slovníků; IDE autocomplete funguje.
- **Operátorové přetížení pro obvody.** `R(100) - (R(5000) | C(1e-6))`
  místo stringového parseru — typově bezpečné a rozšiřitelné.
- **Tři úrovně optimalizace:** single fit, multi-start, differential
  evolution. Volba podle scénáře.

---

## 3. Návrhové vzory

### 3.1 Dataclass Pattern

Konzistentně použit pro výsledky:

```python
@dataclass
class DRTResult:
    tau: Optional[NDArray[np.float64]] = None
    gamma: Optional[NDArray[np.float64]] = None
    lambda_used: Optional[float] = None
    # ...
```

15+ dataclassů napříč moduly. Immutable-friendly, dobře dokumentované.

### 3.2 Operator Overloading

```python
circuit = R(100) - (R(5000) | C(1e-6))
# __sub__ (sériová kombinace), __or__ (paralelní), __mul__ (škálování)
```

### 3.3 Handler Pattern (CLI)

```python
def run_kk_validation(freq, Z, args) -> Figure
def run_drt_analysis(freq, Z, args, R_inf, peak_method) -> DRTResult
def run_circuit_fitting(freq, Z, args) -> (FitResult, Figure)
```

Každý handler 30-100 řádků, jasné oddělení odpovědností.

### 3.4 Strategy Pattern (optimalizace)

```python
--optimizer single      # jeden lokální fit
--optimizer multistart  # multi-start (~85-95% global)
--optimizer de          # differential evolution (~99% global)
```

### 3.5 Composition over Inheritance

`CircuitElement` jako base s 9 konkrétními implementacemi (R, C, Q, L, W,
Wo, K, G, ...), `Series`/`Parallel` builder třídy pro kompozici.

| Strategie | Rychlost | Robustnost | Globální min. |
|---|---|---|---|
| Single fit | rychlá | nízká | ~60% |
| Multi-start | střední | střední | ~85-95% |
| Differential Evolution | pomalá | vysoká | ~99% |

---

## 4. Kvalita kódu

### 4.1 Linting a typy

| Kontrola | Stav |
|---|---|
| `ruff check eis_analysis/` | Čistý |
| `mypy eis_analysis/` | 110 chyb (non-blocking) |
| `pytest tests/` | 75/75 pass (~12 s) |

Mypy v CI od v0.13.5 jako `continue-on-error: true` job. Ze 153 původních
chyb 43 odpadlo díky `ignore_missing_imports`, 2 reálné None-unchecked
přístupy opraveny. Zbývajících 110 jsou postupné cleanupy.

### 4.2 Velikost souborů (limit 500 řádků z CLAUDE.md)

| Soubor | Řádků | Překročení |
|---|---:|---:|
| `cli/handlers.py` | 840 | +68% |
| `drt/core.py` | 810 | +62% |
| `io/data_loading.py` | 650 | +30% |
| `fitting/circuit_elements.py` | 647 | +29% |
| `validation/kramers_kronig.py` | 536 | +7% |
| `fitting/voigt_chain/fitting.py` | 532 | +6% |
| `rinf_estimation/rlk_fit.py` | 529 | +6% |
| `fitting/auto_suggest.py` | 528 | +6% |
| `fitting/circuit.py` | 500 | přesně |

Žádný z modulů není "god object" — všechny mají koherentní obsah, ale
`handlers.py` a `circuit_elements.py` jsou nejočividnější kandidáti na
rozpad. Detailní plán dekompozice viz `doc/AUDIT_2026-04-24.md`, Priorita 5.

### 4.3 Pozitivní aspekty

- Type hints konzistentně použité (`NDArray[np.float64]`, `Optional[...]`)
- NumPy-style docstringy
- Žádné TODO/FIXME/HACK komentáře v kódu
- Vlastní výjimka `EISAnalysisError` pro CLI chyby
- Kompatibilita s Python 3.9+ (bez walrus operátoru a novějších)
- `21× except Exception` — v CLI handlerech rozumné (uživatelská chyba),
  v interních funkcích (`drt/core.py:379`, `drt/peaks.py:98`) skrývá chyby.

---

## 5. Testování

### 5.1 Přehled

```
75 passed in ~12s
```

| Soubor | Testů | Pokrývá |
|---|---:|---|
| `test_cli_integration.py` | 16 | End-to-end workflow |
| `test_jacobian.py` | 15 | Analytic vs numerical Jacobian (přidáno v0.13.4) |
| `test_G_element.py` | 8 | Gerischer element |
| `test_voigt_chain.py` | 8 | Voigt chain fitting |
| `test_weighting.py` | 8 | Váhové funkce |
| `test_hybrid_lambda.py` | 6 | Volba regularizace |
| `test_K_element.py` | 5 | Voigt element |
| `test_confidence_intervals.py` | 5 | Intervaly spolehlivosti |
| `test_correlation_check.py` | 3 | Korelace parametrů |

### 5.2 Pokrytí podle modulů

| Modul | Vlastní testy | Integrační | Hodnocení |
|---|:---:|:---:|---|
| `fitting/` | Ano (5 souborů) | Ano | Dobré |
| `drt/` | Ano (lambda) | Ano | Přijatelné |
| `cli/` | -- | Ano (16 testů) | Přijatelné |
| `validation/` | -- | Částečně | **Nedostatečné** |
| `io/` | -- | Částečně | **Nedostatečné** |
| `rinf_estimation/` | -- | Ano | Přijatelné |
| `visualization/` | -- | -- | **Chybí** |
| `analysis/` | -- | -- | **Chybí** |

**Poměr test/produkce: ~20%**. Pro vědecký software je doporučovaný cíl
30-40%.

### 5.3 Kritické mezery

- **`io/data_loading.py` (650 řádků)** — žádné unit testy parseru DTA/CSV.
  Chybí testy pro poškozené soubory, chybějící sloupce, extrémní hodnoty.
- **`analysis/oxide.py` (375 řádků)** — žádné testy. Riziko tichých regresí.
- **`visualization/` (440 řádků)** — žádné smoke testy.

---

## 6. Vědecká korektnost

| Metoda | Implementace | Hodnocení |
|---|---|---|
| Lin-KK validace | Parametrická, Voigt model | Standardní přístup |
| Z-HIT validace | Neparametrická, Hilbert transform v log-omega | Moderní alternativa |
| DRT (Tikhonov) | Regularizace s auto-lambda | Korektní |
| GCV | Generalized Cross-Validation | Korektní |
| GMM peak detection | scikit-learn GMM | Robustní volba |
| Circuit fitting | L-BFGS-B, multi-start, DE | Tři úrovně robustnosti |
| Analytický Jacobian | Ručně odvozený pro každý element | Otestováno (v0.13.4) |

### Silné stránky

- Duální validace (Lin-KK + Z-HIT) snižuje riziko false positives.
- Automatická volba lambda (GCV + hybrid L-curve) eliminuje subjektivní
  volbu.
- Tři optimalizační strategie pokrývají různé scénáře.
- Srovnání s existujícími nástroji (pyEIS, eisfit.py) dokumentuje shodu
  výsledků.
- **Z-HIT review (v0.13.6)** — externí code review odhalil 6 problémů
  (sort permutace, np.unwrap, kvalitativní stupně, magic numbers, ...);
  5 z 6 opraveno. Detail viz `doc/ZHIT_AUDIT_2026-04-26.md`.

### Pokrytá rizika

- **Jacobian** dříve nebyl testován proti numerické aproximaci. v0.13.4
  přidal `test_jacobian.py` (15 testů systematicky srovnávajících analytic
  a numerický Jacobian).
- `term_classification.py` (klasifikace fyzikálních procesů) stále nemá
  unit testy.

---

## 7. Závislosti a balíčkování

```
numpy>=1.20
scipy>=1.7
matplotlib>=3.4
```

Minimální konzervativní sada. `scikit-learn` je volitelná závislost pro
GMM s graceful degradation. Moderní `pyproject.toml` (PEP 621), dynamická
verze z `eis_analysis/version.py`, entry point `eis` pro CLI, dev závislosti
oddělené (`.[dev]`).

---

## 8. Dokumentace

34+ markdown souborů v `doc/`, ~14 500 řádků (rozsáhlejší než kód —
neobvyklé a pozitivní). Single source of truth pravidla:

- Verze → `eis_analysis/version.py`
- Historie verzí → `CHANGELOG.md`
- API reference → docstrings + `PYTHON_API.md`
- Implementační detaily → README.md

| Kategorie | Souborů |
|---|---:|
| Algoritmické specifikace | 10 |
| Srovnání s jinými nástroji | 4 |
| Uživatelské příručky | 5 |
| Vývojářské příručky | 3 |
| API reference | 1 |
| Audit / review | 3 (`AUDIT_*`, `ZHIT_AUDIT_*`, `CODE_REVIEW.md`) |

**Otevřené:** 4 srovnávací dokumenty v `doc/` jsou stále untracked
(`EISFITPYTHON_COMPARISON.md`, `KK_IMPLEMENTATION_COMPARISON.md`,
`PYEIS_COMPARISON_REPORT.md`, `PYEIS_PREDEFINED_CIRCUITS.md`) — viz
`AUDIT_2026-04-24.md`.

---

## 9. Historie refaktoringů (v0.11 → v0.13.6)

Konsolidovaný přehled významných změn z předchozích reviews:

| Verze | Hlavní změna |
|---|---|
| **v0.12.0** | Rozdělení `eis.py` (1 107 → 146 řádků): nový `cli/` submodul (parser, handlers, data_handling, logging, utils). 87% redukce hlavního souboru. |
| **v0.12.1** | Sjednocení duplicitních `DRTResult` dataclassů. Jediný zdroj pravdy: `eis_analysis.drt.DRTResult`. |
| **v0.13.0** | Přesun loggingu z core do CLI vrstvy. Odstraněno ~240 logger volání napříč 6 moduly. Přidány diagnostické dataclassy (`DRTDiagnostics`, `RinfEstimate`, `LambdaSelection`, `FitDiagnostics`, `MultistartDiagnostics`, `DiffEvoDiagnostics`). |
| **v0.13.1** | Sjednocení error handlingu. `KKResult` rozšířen o `error: Optional[str]` + `success` property. Multi-start tracking selhání přes `failed_errors`. |
| **v0.13.4** | CI pipeline (GitHub Actions): ruff + pytest na Python 3.9 a 3.12. Přidán `test_jacobian.py` (analytic vs numerical Jacobian). |
| **v0.13.5** | Oprava prohozených `log_separator` argumentů (latentní bug). KK return type → `KKResult` (bez Optional). Mypy přidán do CI (non-blocking). |
| **v0.13.6** | Z-HIT validation fixes z externího review: un-sort výstupů, `np.unwrap` na fázi, stratifikované quality labels, data-driven `offset_center`. |

---

## 10. Soulad s CLAUDE.md

| Princip | Dodržení | Poznámka |
|---|---|---|
| Modularita | Plně | Každý modul má jasnou odpovědnost |
| Single Responsibility | Plně | Žádné god objects |
| DRY | Plně | Minimální duplikace |
| Testability | Převážně | Core testovatelný, viz mezery v sekci 5 |
| Documentation | Plně | NumPy-style docstrings, doc/ adresář |
| No emoji in docs | Plně | PDF export bezpečný |
| Version management | Plně | Single source v `version.py` |
| Git workflow | Plně | Deskriptivní commity, žádný force-push |
| **Limit 500 řádků** | **Částečně** | 8 souborů přes limit (sekce 4.2) |

---

## 11. Otevřené priority

Aktuální seznam otevřených úkolů a jejich stav je v
**`doc/AUDIT_2026-04-24.md`**. Stručný výtah:

- **Priorita 4** — konsolidace 4 review reportů ✅ **(hotovo, viz tento
  dokument)**
- **Priorita 5** — rozdělit soubory nad 500 řádků (`cli/handlers.py`,
  `drt/core.py`, ...).
- Rozhodnout o 4 untracked COMPARISON dokumentech v `doc/`.
- Postupný cleanup zbylých 110 mypy chyb.
- Doplnit testy pro `io/data_loading.py`, `analysis/oxide.py`,
  `visualization/`, `term_classification.py`.

---

## 12. Závěr

Projekt EIS Analysis je pro jednoautorský vědecký software **výborný**.
Architektura je promyšlená, vědecké metody jsou správně implementované,
dokumentace je nadstandardní.

Hlavní rizika neleží v kvalitě kódu, ale v **udržitelnosti při dalším
růstu**:

1. Testovací pokrytí (~20%) je nedostatečné pro kritické moduly
   (I/O parser, doménová analýza, vizualizace).
2. CI/CD od v0.13.4 existuje (ruff + pytest + mypy non-blocking) — kvalita
   už není jen na disciplíně jediného vývojáře.
3. Některé soubory překračují vlastní limity velikosti a zaslouží
   dekompozici.

Projekt má solidní základy pro další rozvoj. Investice do testovací
infrastruktury (rozšíření testů na netestované moduly + zpřísnění mypy
v CI) by výrazně snížila riziko regresí.
