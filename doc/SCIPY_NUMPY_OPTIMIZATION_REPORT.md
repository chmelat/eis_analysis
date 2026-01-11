# Report: Příležitosti pro náhradu scipy/numpy funkcemi

**EIS Analysis Toolkit v0.13.2**

Datum: 2026-01-11

---

## Shrnutí

Kód je **velmi dobře napsaný** a již efektivně využívá scipy/numpy knihovny. Nalezeny byly **4 příležitosti** pro zlepšení, z toho 1 s vysokou prioritou.

| Priorita | Počet | Popis |
|----------|-------|-------|
| Vysoká | 1 | Manuální kumulativní integrace |
| Střední | 2 | Manuální derivace, konstrukce matice |
| Nízká | 1 | Vektorizace for cyklu |

---

## Již správně implementované vzory

Před nálezem nedostatků je důležité zdůraznit, že kód **již správně používá**:

| Funkce | Použití | Soubor |
|--------|---------|--------|
| `scipy.signal.find_peaks()` | Detekce píků v DRT | auto_suggest.py |
| `scipy.optimize.nnls()` | Non-negative least squares | solvers.py |
| `scipy.optimize.least_squares()` | Circuit fitting | circuit.py |
| `scipy.optimize.differential_evolution()` | Globální optimalizace | diffevo.py |
| `np.gradient()` | Numerická derivace | zhit.py:236 |
| `np.linalg.solve()` | Řešení lin. systémů | solvers.py |
| `np.linalg.lstsq()` | Least squares via SVD | solvers.py |
| `np.searchsorted()` | Efektivní hledání indexů | peaks.py |
| `np.polyfit()` | Polynomiální fit | rlk_fit.py |
| `np.trapz()` | Numerická integrace | peaks.py, auto_suggest.py |

---

## Nalezené příležitosti pro zlepšení

### 1. Manuální kumulativní integrace (VYSOKÁ priorita)

**Soubor:** `eis_analysis/validation/zhit.py`
**Řádky:** 220-229
**Úspora:** ~10 řádků

#### Současný kód:

```python
ln_Z_first_order = np.zeros(n)
ln_Z_first_order[ref_idx] = ln_Z_ref

# Integrate forward from reference point (ref_idx to end)
for i in range(ref_idx + 1, n):
    d_ln_omega = ln_omega[i] - ln_omega[i - 1]
    phi_avg = 0.5 * (phi[i] + phi[i - 1])
    ln_Z_first_order[i] = ln_Z_first_order[i - 1] + (2.0 / np.pi) * phi_avg * d_ln_omega

# Integrate backward from reference point (ref_idx to start)
for i in range(ref_idx - 1, -1, -1):
    d_ln_omega = ln_omega[i + 1] - ln_omega[i]
    phi_avg = 0.5 * (phi[i] + phi[i + 1])
    ln_Z_first_order[i] = ln_Z_first_order[i + 1] - (2.0 / np.pi) * phi_avg * d_ln_omega
```

#### Navrhovaná náhrada:

```python
from scipy.integrate import cumulative_trapezoid

# Cumulative integral from start
integral_full = cumulative_trapezoid((2.0 / np.pi) * phi, ln_omega, initial=0)

# Shift to reference point
ln_Z_first_order = ln_Z_ref + (integral_full - integral_full[ref_idx])
```

#### Výhody:
- Redukce z 10 na 3 řádky
- Eliminace manuálních smyček
- Konzistence s `np.gradient()` již použitým na řádku 236
- Lepší čitelnost a údržba

---

### 2. Manuální výpočet křivosti (STŘEDNÍ priorita)

**Soubor:** `eis_analysis/drt/gcv.py`
**Řádky:** 164-180
**Úspora:** ~8 řádků

#### Současný kód:

```python
curvature = np.zeros(n)

for i in range(1, n - 1):
    # První derivace (centrální)
    d_rho = (rho[i + 1] - rho[i - 1]) / 2
    d_eta = (eta[i + 1] - eta[i - 1]) / 2

    # Druhá derivace (centrální)
    dd_rho = rho[i + 1] - 2 * rho[i] + rho[i - 1]
    dd_eta = eta[i + 1] - 2 * eta[i] + eta[i - 1]

    # Křivost
    numerator = d_rho * dd_eta - dd_rho * d_eta
    denominator = (d_rho**2 + d_eta**2)**1.5

    if abs(denominator) > 1e-15:
        curvature[i] = numerator / denominator
    else:
        curvature[i] = 0.0
```

#### Navrhovaná náhrada:

```python
# První derivace pomocí np.gradient()
d_rho = np.gradient(rho)
d_eta = np.gradient(eta)

# Druhá derivace
dd_rho = np.gradient(d_rho)
dd_eta = np.gradient(d_eta)

# Křivost (vektorizovaně)
numerator = d_rho * dd_eta - dd_rho * d_eta
denominator = (d_rho**2 + d_eta**2)**1.5
curvature = np.where(np.abs(denominator) > 1e-15, numerator / denominator, 0.0)
```

#### Výhody:
- Plně vektorizované
- Využívá `np.gradient()` konzistentně
- Eliminace for cyklu
- `np.where()` pro podmíněné přiřazení

---

### 3. Konstrukce regularizační matice L (STŘEDNÍ priorita)

**Soubor:** `eis_analysis/drt/core.py`
**Řádky:** 337-341
**Úspora:** ~4 řádky

#### Současný kód:

```python
L = np.zeros((n_tau - 2, n_tau))
for i in range(n_tau - 2):
    L[i, i] = 1
    L[i, i + 1] = -2
    L[i, i + 2] = 1
```

#### Navrhovaná náhrada (varianta A - numpy):

```python
# Tridiagonální matice [1, -2, 1] pro druhou derivaci
L = (np.diag(np.ones(n_tau - 2), 0)
     - 2 * np.diag(np.ones(n_tau - 2), 1)[:n_tau-2, :n_tau]
     + np.diag(np.ones(n_tau - 2), 2)[:n_tau-2, :n_tau])
```

#### Navrhovaná náhrada (varianta B - scipy.sparse, doporučeno):

```python
from scipy.sparse import diags

# Efektivnější pro velké matice
L = diags([1, -2, 1], [0, 1, 2], shape=(n_tau - 2, n_tau)).toarray()
```

#### Výhody:
- Jednořádková konstrukce
- `scipy.sparse.diags()` je efektivnější pro velké matice
- Jasný záměr (tridiagonální struktura)

---

### 4. Vektorizace Voigt matice (NÍZKÁ priorita)

**Soubor:** `eis_analysis/fitting/voigt_chain/solvers.py`
**Řádky:** 185-187
**Úspora:** ~2 řádky

#### Současný kód:

```python
for i, tau_i in enumerate(tau):
    tau_omega_sq = (tau_i * omega) ** 2
    A[:, col_offset + i] = 1.0 / (1 + tau_omega_sq)
```

#### Navrhovaná náhrada:

```python
# Broadcasting: tau (n_tau,) x omega (n_freq,) -> (n_tau, n_freq)
tau_omega_sq = (tau[:, np.newaxis] * omega[np.newaxis, :]) ** 2
A[:, col_offset:col_offset + len(tau)] = (1.0 / (1 + tau_omega_sq)).T
```

#### Výhody:
- Eliminace for cyklu
- Plné využití NumPy broadcasting
- Marginální zlepšení výkonu (n_tau je typicky malé)

**Poznámka:** Priorita je nízká, protože `n_tau` je typicky <= 20 a současný kód je dostatečně rychlý.

---

## Příklady kódu, který je již optimální

### Detekce píků (auto_suggest.py:152-157)

```python
from scipy.signal import find_peaks

peaks, properties = find_peaks(
    gamma,
    height=np.max(gamma) * DRT_PEAK_HEIGHT_THRESHOLD,
    distance=min_distance,
    prominence=np.max(gamma) * DRT_PEAK_PROMINENCE_THRESHOLD
)
```

Toto je **vzorové použití** scipy.

### Derivace fáze (zhit.py:236)

```python
d_phi_d_ln_omega = np.gradient(phi, ln_omega)
```

Správné použití `np.gradient()` pro neekvidistantní data.

### NNLS s fallbacky (solvers.py:16-112)

```python
def robust_nnls(A, b, max_iter=None):
    # Try scipy.optimize.nnls first
    # Fallback to scipy.optimize.lsq_linear
    # Last resort: pseudoinverse + clipping
```

Robustní implementace s graceful degradation.

---

## Doporučení pro implementaci

### Vysoká priorita (doporučeno implementovat)

1. **zhit.py** - Kumulativní integrace
   - Jednoduchá změna
   - Významná redukce kódu
   - Konzistence s existujícím stylem

### Střední priorita (zvážit)

2. **gcv.py** - Výpočet křivosti
   - Zlepší čitelnost
   - Menší riziko off-by-one chyb

3. **core.py** - Regularizační matice
   - Elegantnější konstrukce
   - Zvážit scipy.sparse pro budoucí škálovatelnost

### Nízká priorita (kosmetické)

4. **solvers.py** - Vektorizace Voigt matice
   - Marginální přínos
   - Současný kód je dostatečný

---

## Závěr

Kód EIS Analysis Toolkit je **vysoce kvalitní** z pohledu využití vědeckých Python knihoven. Většina numerických operací již používá optimální scipy/numpy funkce.

Identifikované příležitosti jsou převážně **kosmetické zlepšení**, které zvýší čitelnost a konzistenci, ale nezmění zásadně výkon ani funkcionalitu.

**Celkové hodnocení využití scipy/numpy: A-**

Jediná významnější příležitost je náhrada manuální kumulativní integrace v Z-HIT validaci, která by redukovala 10 řádků kódu na 3.
