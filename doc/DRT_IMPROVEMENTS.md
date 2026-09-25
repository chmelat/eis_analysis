# Návrhy Vylepšení eis_analysis.py
## Inspirováno pyDRTtools, zachovává jednoduchost

---

## 1. Automatický Výběr λ pomocí GCV ⭐⭐⭐⭐⭐

**Užitečnost:** Kritická - eliminuje manuální tuning
**Složitost:** Střední (~50 řádků)
**Závislosti:** Žádné nové

### Princip
```python
GCV(λ) = (||Z - A·x(λ)||²) / (trace(I - K(λ)))²

kde K(λ) = A @ inv(A^T·A + λ·M) @ A^T
```

### Implementace
- Funkce `compute_gcv(lambda_val, A, b, M)` → vrací GCV score
- Funkce `find_optimal_lambda(A, b, M)` → minimalizuje GCV v log-prostoru
- Použití: `scipy.optimize.minimize_scalar` v bounds [1e-5, 1.0]

### Výhody
✅ Objektivní volba λ bez experimentování
✅ Robustní pro různé typy dat
✅ Rychlé (~1-2s navíc)

---

## 2. Fit R_∞ jako Parametr ⭐⭐⭐⭐

**Stav (0.41.0): přehodnoceno, jen jako diagnostika.** Numerická studie
`DRT_RINF_L_ANALYSIS_2026-09-25.md` ukázala, že bod níže situaci podceňuje:
- "Současný stav" už neplatí. R_inf je HF medián (až 5 bodů), nebo preset
  z `--ri-fit` (fit R-L-(R|Q) s kontrolou identifikovatelnosti, 0.40.0).
- Nejde o "minimální změnu". λ se musí vybírat na rozšířeném systému (jinak
  je bias -20 %), mez R_inf = 0 je třeba ošetřit a nevážený fit nechá malé
  R_inf zmizet.
- "Teoreticky přesnější" platí jen proti mediánu. Proti `--ri-fit` z 0.40.0
  prohrává ve všech určitelných případech (C -2.9 % proti -0.5 %) a na reálných
  datech padá na mez 0.
- Indukčnost, kterou studie označila za větší problém, je v modelu od 0.41.0
  (`calculate_drt(inductance=...)`).

Společný odhad má smysl nanejvýš jako diagnostika vedle externího R_inf
(D2 ve studii). Původní návrh:

**Užitečnost:** Vysoká - přesnější odhad R_inf
**Složitost:** Nízká (~10 řádků)
**Závislosti:** Žádné nové

### Současný stav
```python
R_inf = Z.real[high_freq_idx]  # Odhad z jednoho bodu
b = np.concatenate([Z.real - R_inf, Z.imag])
```

### Vylepšení
```python
# Rozšířit A o sloupec pro R_∞
A_re_extended = np.column_stack([np.ones(len(freq)), A_re])
A_im_extended = np.column_stack([np.zeros(len(freq)), A_im])
A_combined = np.vstack([A_re_extended, A_im_extended])

# Řešit pro x = [R_∞, γ₁, ..., γ_N]
x, residual = nnls(A_reg, b_reg)
R_inf = x[0]
gamma = x[1:]
```

### Výhody
✅ Simultánní fit - teoreticky přesnější
✅ Minimální změna kódu
✅ Žádná změna rychlosti

---

## 3. Re-Im Cross Validation ⭐⭐⭐

**Užitečnost:** Střední - robustnější λ pro šumová data
**Složitost:** Střední (~40 řádků)
**Závislosti:** Žádné nové

### Princip
1. Řeš DRT pouze z Re části → γ_re, R_∞_re
2. Řeš DRT pouze z Im části → γ_im
3. Křížová predikce:
   - Predikuj Im z Re parametrů: `Z_im_pred = A_im @ γ_re`
   - Predikuj Re z Im parametrů: `Z_re_pred = R_∞_re + A_re @ γ_im`
4. Score = `||Z_re - Z_re_pred||² + ||Z_im - Z_im_pred||²`

### Implementace
```python
def compute_re_im_cv(lambda_val, A_re, A_im, Z_re, Z_im, M):
    # Fit z Re dat
    x_re, _ = nnls(A_re_reg, Z_re)
    R_inf_re, gamma_re = x_re[0], x_re[1:]

    # Fit z Im dat
    x_im, _ = nnls(A_im_reg, Z_im)
    gamma_im = x_im[1:]  # bez R_inf

    # Křížová predikce
    Z_re_pred = R_inf_re + A_re @ gamma_im
    Z_im_pred = A_im @ gamma_re

    # Score
    return np.sum((Z_re - Z_re_pred)**2) + np.sum((Z_im - Z_im_pred)**2)
```

### Výhody
✅ Testuje konzistenci mezi Re a Im
✅ Robustnější k outlierům v jedné části
✅ Stále rychlé (~2-3s)

---

## 4. Diagnostické Metriky ⭐⭐⭐⭐

**Užitečnost:** Vysoká - objektivní hodnocení kvality
**Složitost:** Velmi nízká (~20 řádků)
**Závislosti:** Žádné nové

### Metriky k přidání

```python
def compute_diagnostics(Z_exp, Z_fit, gamma, tau, lambda_reg):
    """Vypočítá diagnostické metriky DRT fitu."""

    # 1. Relativní chyba Re a Im samostatně
    rel_error_re = np.mean(np.abs(Z_exp.real - Z_fit.real) / np.abs(Z_exp.real)) * 100
    rel_error_im = np.mean(np.abs(Z_exp.imag - Z_fit.imag) / np.abs(Z_exp.imag)) * 100

    # 2. Chi-squared (normalizovaná rezidua)
    chi_squared = np.sum(np.abs(Z_exp - Z_fit)**2 / np.abs(Z_exp)**2) / len(Z_exp)

    # 3. R-squared
    Z_mean = np.mean(np.abs(Z_exp))
    SS_res = np.sum(np.abs(Z_exp - Z_fit)**2)
    SS_tot = np.sum((np.abs(Z_exp) - Z_mean)**2)
    R_squared = 1 - SS_res / SS_tot

    # 4. Integrovaný polarizační odpor
    R_pol = np.trapz(gamma, np.log(tau))

    # 5. Efektivní počet píků (Shannon entropyová míra)
    gamma_norm = gamma / np.sum(gamma)
    gamma_norm = gamma_norm[gamma_norm > 0]
    entropy = -np.sum(gamma_norm * np.log(gamma_norm))
    n_eff_peaks = np.exp(entropy)

    # 6. Regularizační penalizace
    if lambda_reg > 0:
        # Druhá derivace
        d2_gamma = np.diff(gamma, n=2)
        smoothness = np.sum(d2_gamma**2)
    else:
        smoothness = 0.0

    return {
        'rel_error_re': rel_error_re,
        'rel_error_im': rel_error_im,
        'chi_squared': chi_squared,
        'R_squared': R_squared,
        'R_pol': R_pol,
        'n_eff_peaks': n_eff_peaks,
        'smoothness': smoothness
    }
```

### Výstup
```
DRT Diagnostika:
  Re chyba:        0.82 %  ✓
  Im chyba:        1.15 %  ⚠
  χ²:              0.0024  ✓
  R²:              0.9987  ✓
  R_pol:           5843 Ω
  Efektivní píky:  2.3
  Smoothness:      0.0012
```

### Výhody
✅ Objektivní hodnocení kvality
✅ Detekce problémů (vysoká Im chyba → možná induktivita)
✅ Srovnatelnost mezi různými měřeními

---

## 5. L-Curve Vizualizace ⭐⭐

**Užitečnost:** Střední - pomáhá pochopit volbu λ
**Složitost:** Nízká (~30 řádků)
**Závislosti:** Žádné nové

### Implementace
```python
def plot_l_curve(A, b, M, lambda_range=None):
    """Vykresli L-curve pro volbu regularizačního parametru."""

    if lambda_range is None:
        lambda_range = np.logspace(-5, 0, 30)

    residuals = []
    regularizations = []

    for lam in lambda_range:
        A_reg = np.vstack([A, np.sqrt(lam) * M])
        b_reg = np.concatenate([b, np.zeros(M.shape[0])])

        x, _ = nnls(A_reg, b_reg)

        residuals.append(np.linalg.norm(b - A @ x))
        regularizations.append(lam * np.linalg.norm(M @ x))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # L-curve
    ax1.loglog(residuals, regularizations, 'o-', linewidth=2)
    ax1.set_xlabel('||Z - A·γ|| (residual)')
    ax1.set_ylabel('λ||M·γ|| (smoothness)')
    ax1.set_title('L-Curve')
    ax1.grid(True, alpha=0.3)

    # Lambda vs metrics
    ax2.semilogx(lambda_range, residuals, 'o-', label='Residual', linewidth=2)
    ax2_twin = ax2.twinx()
    ax2_twin.semilogx(lambda_range, regularizations, 's-',
                      color='orange', label='Regularization', linewidth=2)
    ax2.set_xlabel('λ')
    ax2.set_ylabel('Residual', color='blue')
    ax2_twin.set_ylabel('Regularization', color='orange')
    ax2.set_title('Trade-off křivky')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
```

### Použití
```bash
python eis_analysis.py data.DTA --plot-lcurve
```

### Výhody
✅ Vizuální pochopení trade-off mezi fit a smoothness
✅ Pomáhá najít "loket" (optimal λ)
✅ Srovnání různých λ bez opakovaného spouštění

---

## 6. Toeplitz Optimalizace ⭐⭐

**Užitečnost:** Střední - 2-5× zrychlení
**Složitost:** Nízká (~15 řádků)
**Závislosti:** `scipy.linalg.toeplitz`

### Kdy použít
Pouze když:
- Frekvence jsou log-uniformní (std < 1% mean)
- `N_freqs == N_taus`

### Implementace
```python
from scipy.linalg import toeplitz

# Kontrola uniformity
d_ln_freq = np.diff(np.log(frequencies))
is_uniform = np.std(d_ln_freq) / np.mean(d_ln_freq) < 0.01

if is_uniform and len(frequencies) == n_tau:
    # Vypočti pouze první řádek a sloupec
    first_row = [...]  # A_re[0, :]
    first_col = [...]  # A_re[:, 0]
    A_re = toeplitz(first_col, first_row)
else:
    # Původní výpočet
    A_re = ...
```

### Výhody
✅ 2-5× rychlejší sestavení matice pro velká N
✅ Automatická detekce použitelnosti

---

## Doporučené Pořadí Implementace

### Fáze 1: Quick Wins (1-2 hodiny)
1. **Diagnostické metriky** - okamžitá hodnota
2. **Fit R_∞ jako parametr** - malá změna, velký přínos

### Fáze 2: Game Changer (3-4 hodiny)
3. **GCV automatický výběr λ** - eliminuje hlavní slabinu

### Fáze 3: Pokročilé (2-3 hodiny)
4. **Re-Im cross validation** - alternativa ke GCV
5. **L-curve vizualizace** - debugging tool

### Fáze 4: Optimalizace (1 hodina)
6. **Toeplitz optimalizace** - pouze pokud potřeba

---

## Příklad Vylepšené Workflow

```bash
# Před vylepšením
python eis_analysis.py data.DTA --lambda 0.1
# Kontrola grafu, ručně upravit λ...
python eis_analysis.py data.DTA --lambda 0.05
# Opět kontrola...

# Po vylepšení
python eis_analysis.py data.DTA --auto-lambda gcv
# ✓ Automaticky najde optimální λ = 0.073
# ✓ Zobrazí diagnostiku: Re 0.8%, Im 1.1%, R² 0.998
# ✓ R_inf fitován současně: 95.3 Ω
```

---

## Co NEIMPLEMENTOVAT (příliš složité)

❌ **RBF diskretizace** - Vyžaduje přepis celé matice A, komplikované integrály
❌ **CVXOPT QP solver** - NNLS je rychlejší a dostačující
❌ **Bayesian DRT** - Vyžaduje MCMC sampling, tisíce řádků
❌ **Hilbert transform validation** - Komplexní matematika
❌ **Nearest PD correction** - Zbytečné s NNLS

---

## Souhrn

| Vylepšení | Užitečnost | Složitost | Čas | Priorita |
|-----------|------------|-----------|-----|----------|
| GCV auto λ | ⭐⭐⭐⭐⭐ | Střední | 3h | 🔥 KRITICKÁ |
| Fit R_∞ | ⭐⭐⭐⭐ | Velmi nízká | 1h | 🔥 Vysoká |
| Diagnostika | ⭐⭐⭐⭐ | Velmi nízká | 1h | 🔥 Vysoká |
| Re-Im CV | ⭐⭐⭐ | Střední | 2h | Střední |
| L-curve plot | ⭐⭐ | Nízká | 1h | Nízká |
| Toeplitz | ⭐⭐ | Nízká | 1h | Nízká |

**Doporučení:** Implementuj minimálně první 3 (GCV, Fit R_∞, Diagnostika) = ~5 hodin práce, transformační přínos.
