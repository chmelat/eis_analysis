# Code Review Report: EIS Analysis Toolkit

**Datum:** 2026-01-09
**Verze:** 0.11.3
**Reviewer:** Claude Code (Opus 4.5)
**Rozsah:** Kritická analýza z hlediska čistoty designu a jasnosti zápisu

---

## Executive Summary

EIS Analysis Toolkit je **dobře strukturovaný** vědecký software s jasnou modulární architekturou. Kód je převážně čitelný a dobře dokumentovaný. Hlavní oblasti pro zlepšení jsou:

1. **Přerostlý CLI soubor** (1107 řádků) - potřeba refactoring
2. **Duplicitní datové struktury** - porušení DRY principu
3. **Verbose logging v core modulech** - míchání odpovědností

Celkové hodnocení: **B+** (dobré s prostorem pro zlepšení)

---

## 1. Metriky projektu

### 1.1 Rozsah kódu

| Kategorie | Počet souborů | Celkem řádků |
|-----------|---------------|--------------|
| Core package | 32 | 10,878 |
| Testy | 11 | 2,422 |
| CLI | 1 | 1,107 |
| **Celkem** | **44** | **14,407** |

### 1.2 Velikost modulů

| Soubor | Řádků | Stav |
|--------|-------|------|
| eis.py | 1,107 | KRITICKÉ - příliš velký |
| drt/core.py | 837 | Akceptovatelné |
| io/data_loading.py | 650 | Akceptovatelné |
| rinf_estimation/rlk_fit.py | 612 | Akceptovatelné |
| validation/kramers_kronig.py | 578 | Akceptovatelné |
| fitting/circuit_elements.py | 556 | Akceptovatelné |

**Pravidlo:** Soubory nad 500 řádků jsou kandidáti na refactoring. Soubory nad 1000 řádků vyžadují rozdělení.

---

## 2. Architektura

### 2.1 Modulární struktura (POZITIVNÍ)

```
eis_analysis/
├── io/                 # Data I/O - jasně oddělené
├── validation/         # KK a Z-HIT validace
├── drt/                # DRT algoritmy
├── fitting/            # Circuit fitting
│   └── voigt_chain/    # Správně vyčleněný submodul
├── rinf_estimation/    # R_inf odhad
├── analysis/           # Doménová logika
├── visualization/      # Grafy - oddělené od algoritmů
└── utils/              # Utility funkce
```

**Hodnocení:** Výborná separace concerns. Každý modul má jasnou odpovědnost.

### 2.2 Dependency flow

```
Data Loading → Validation → R_inf Estimation → DRT → Circuit Fitting → Visualization
     ↓              ↓              ↓             ↓           ↓              ↓
   io/        validation/    rinf_estimation/  drt/      fitting/    visualization/
```

**Hodnocení:** Správný unidirectional flow bez cyklických závislostí.

---

## 3. Identifikované problémy

### 3.1 KRITICKÉ: Přerostlý CLI soubor

**Soubor:** `eis.py` (1107 řádků)

**Popis:** Hlavní CLI soubor obsahuje příliš mnoho odpovědností:
- Custom logging formattery (řádky 75-107)
- Argument parsing (řádky 330-492)
- Data loading logika (řádky 499-607)
- Analytické funkce (řádky 614-1000)
- Main orchestrace (řádky 1034-1107)

**Porušuje:** Single Responsibility Principle

**Dopad:**
- Obtížná navigace a údržba
- Komplikované unit testování
- Vysoká kognitivní zátěž při čtení

**Doporučené řešení:**
```
eis_analysis/cli/
├── __init__.py
├── logging.py      # Custom formattery (33 řádků)
├── parser.py       # Argument parsing (163 řádků)
├── pipeline.py     # Orchestrace workflow (200 řádků)
├── handlers.py     # Analytické funkce (400 řádků)
└── utils.py        # Helper funkce (50 řádků)
```

**Priorita:** VYSOKÁ

---

### 3.2 VYSOKÁ: Duplicitní DRTResult dataclass

**Lokace:**
- `eis.py:142-148`
- `eis_analysis/drt/core.py:43-71`

**Popis:** Existují dvě verze stejné datové struktury s různými poli a konvencemi pojmenování.

```python
# eis.py - zjednodušená verze
@dataclass
class DRTResult:
    tau: Optional[NDArray]
    gamma: Optional[NDArray]
    peaks_gmm: Optional[list]      # <- "peaks_gmm"
    fig_drt: Optional[plt.Figure]  # <- "fig_drt"
    fig_rinf: Optional[plt.Figure]

# drt/core.py - plná verze
@dataclass
class DRTResult:
    tau: Optional[NDArray] = None
    gamma: Optional[NDArray] = None
    gamma_original: Optional[NDArray] = None
    figure: Optional[plt.Figure] = None  # <- "figure"
    peaks: Optional[List[Dict]] = None   # <- "peaks"
    figure_rinf: Optional[plt.Figure] = None
    R_inf: Optional[float] = None
    R_pol: Optional[float] = None
    lambda_used: Optional[float] = None
    reconstruction_error: Optional[float] = None
```

**Porušuje:** DRY (Don't Repeat Yourself)

**Dopad:**
- Nejednoznačnost která verze se má použít
- Nutnost manuálního mapování mezi verzemi
- Riziko nekonzistence při změnách

**Doporučené řešení:**
1. Ponechat pouze verzi v `drt/core.py`
2. Importovat v `eis.py`: `from eis_analysis.drt import DRTResult`
3. Sjednotit pojmenování polí

**Priorita:** VYSOKÁ

---

### 3.3 STŘEDNÍ: Verbose logging v core modulech

**Lokace:** `drt/core.py:181-261` (funkce `_estimate_r_inf`)

**Popis:** Core algoritmy obsahují rozsáhlý logging místo čistého vrácení dat.

```python
def _estimate_r_inf(...):
    # 38 volání loggeru v 80 řádcích:
    logger.info("="*60)
    logger.info("R_inf estimation (high-frequency resistance)")
    logger.info("="*60)
    logger.info(f"R_inf = {R_inf:.3f} Ω (median of {n_avg} highest frequencies)")
    # ... dalších 30+ volání
```

**Porušuje:** Separation of Concerns

**Dopad:**
- Míchání business logiky s prezentační vrstvou
- Obtížné použití jako knihovna (nežádoucí výstup)
- Komplikované testování

**Doporučené řešení:**
```python
# Core funkce vrací data
def _estimate_r_inf(...) -> RinfResult:
    result = RinfResult(R_inf=..., diagnostics=...)
    return result

# CLI vrstva loguje
def run_rinf_estimation(...):
    result = _estimate_r_inf(...)
    logger.info(f"R_inf = {result.R_inf:.3f} Ω")
```

**Priorita:** STŘEDNÍ

---

### 3.4 STŘEDNÍ: Magic numbers bez dokumentace

**Lokace:** `eis.py:518-520`, `drt/core.py:33-36`

**Popis:** Numerické konstanty bez vysvětlení jejich původu.

```python
# eis.py - syntetická data bez vysvětlení
frequencies, Z = generate_synthetic_data(
    Rs=10, R0=1e5, Q0=(1e-6, 0.6),   # Proč 0.6?
    R1=8e5, Q1=(3e-5, 0.43),          # Proč 0.43?
    noise=0.01                         # Proč 1%?
)

# drt/core.py - konstanty bez zdůvodnění
MIN_FREQUENCY_RANGE = 10      # Proč 10?
DRT_TOLERANCE = 1e-10         # Proč 1e-10?
GAMMA_MAX_REASONABLE = 1e10   # Odkud tato hodnota?
```

**Doporučené řešení:**
```python
# Dokumentované konstanty
MIN_FREQUENCY_RANGE = 10  # f_max/f_min ratio pro spolehlivou DRT analýzu
                          # Empiricky ověřeno: < 10 způsobuje špatné rozlišení

DRT_TOLERANCE = 1e-10  # Numerická tolerance pro NNLS solver
                       # Menší hodnoty způsobují numerickou nestabilitu
```

**Priorita:** STŘEDNÍ

---

### 3.5 STŘEDNÍ: Nekonzistentní error handling

**Lokace:** Napříč projektem

**Popis:** Některé funkce vracejí `None` při chybě, jiné vyhazují výjimky.

```python
# Pattern A: Vrací None
if result is None:
    return None

# Pattern B: Vyhazuje výjimku
if not os.path.exists(args.input):
    raise EISAnalysisError(f"File '{args.input}' does not exist!")

# Pattern C: Tichý fallback
except Exception as e:
    logger.warning("Fallback to median method")
```

**Doporučené řešení:**
- Výjimky pro neočekávané stavy
- `Optional` return pro legitimně prázdné výsledky
- Explicitní Result pattern pro operace které mohou selhat

**Priorita:** STŘEDNÍ

---

### 3.6 NÍZKÁ: Nestrukturované CLI argumenty

**Lokace:** `eis.py:330-492`

**Popis:** 163 řádků argumentů bez logického seskupení.

**Doporučené řešení:**
```python
# Použít argument groups
drt_group = parser.add_argument_group('DRT Analysis')
drt_group.add_argument('--lambda', ...)
drt_group.add_argument('--n-tau', ...)

kk_group = parser.add_argument_group('Kramers-Kronig Validation')
kk_group.add_argument('--mu-threshold', ...)
```

**Priorita:** NÍZKÁ

---

### 3.7 NÍZKÁ: Eval pro circuit parsing

**Lokace:** `eis.py:222-224`

**Popis:** Použití `eval()` pro parsování circuit expressions.

```python
circuit = eval(expr, {"__builtins__": {}}, safe_namespace)
```

**Hodnocení:** Riziko je minimalizované prázdným `__builtins__`, ale pro produkční kód by byl vhodnější vlastní parser.

**Priorita:** NÍZKÁ (aktuálně bezpečné)

---

## 4. Pozitivní aspekty

### 4.1 Elegantní operator overloading

**Lokace:** `fitting/circuit_elements.py:101-118`

```python
def __sub__(self, other):  # Series: Z_total = Z1 + Z2
    return Series([self, other])

def __or__(self, other):   # Parallel: 1/Z_total = 1/Z1 + 1/Z2
    return Parallel([self, other])
```

**Výsledek:** Intuitivní syntax `R(100) - (R(5000) | C(1e-6))`

### 4.2 Dataclasses pro strukturované výsledky

```python
@dataclass
class FitResult:
    circuit: Circuit
    params_opt: NDArray[np.float64]
    params_stderr: NDArray[np.float64]
    fit_error_rel: float
    fit_error_abs: float
    quality: str
    # ...
```

### 4.3 Type hints

Většina funkcí má type hints:
```python
def fit_equivalent_circuit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    circuit: Circuit,
    weighting: str = 'modulus',
) -> Tuple[FitResult, NDArray[np.complex128], Optional[plt.Figure]]:
```

### 4.4 Comprehensive docstrings

NumPy-style dokumentace s Parameters, Returns, Examples:
```python
def read_gamry_native(filename: str) -> Tuple[NDArray, NDArray]:
    """
    Native parser for Gamry .DTA files.

    Parameters
    ----------
    filename : str
        Path to .DTA file

    Returns
    -------
    frequencies : ndarray of float
        Frequency values [Hz]
    Z : ndarray of complex
        Complex impedance values [Ω]
    """
```

---

## 5. Skóre kvality

| Kritérium | Skóre | Poznámka |
|-----------|-------|----------|
| Modulární architektura | 9/10 | Výborná separace |
| Single Responsibility | 6/10 | CLI soubor porušuje |
| DRY | 7/10 | Duplicitní dataclasses |
| Čitelnost | 8/10 | Dobré pojmenování |
| Dokumentace | 9/10 | Comprehensive docstrings |
| Type safety | 8/10 | Dobré type hints |
| Error handling | 6/10 | Nekonzistentní |
| Testovatelnost | 7/10 | Core moduly OK, CLI těžké |

**Celkové skóre: 7.5/10**

---

## 6. Akční plán

### Fáze 1: Kritické (týden 1-2)
- [ ] Rozdělit `eis.py` na menší moduly
- [ ] Sjednotit `DRTResult` na jedno místo
- [ ] Odstranit duplicitní definice

### Fáze 2: Střední priorita (týden 3-4)
- [ ] Přesunout logging z core do CLI vrstvy
- [ ] Dokumentovat magic numbers
- [ ] Sjednotit error handling strategii

### Fáze 3: Nízká priorita (podle kapacity)
- [ ] Strukturovat CLI argumenty do groups
- [ ] Vyčistit zastaralé komentáře
- [ ] Přesunout inline importy na začátek souborů

---

## 7. Závěr

EIS Analysis Toolkit je **kvalitní vědecký software** s dobrou architekturou. Hlavní problémy jsou v CLI vrstvě, která přerostla a potřebuje refactoring. Core moduly (DRT, fitting, validation) jsou dobře navržené.

**Doporučení:** Prioritně řešit rozdělení `eis.py` a sjednocení datových struktur před přidáváním nových funkcí.

---

*Report vygenerován: 2026-01-09*
*Nástroj: Claude Code (Opus 4.5)*
