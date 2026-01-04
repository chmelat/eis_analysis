# Technická Analýza Metody DRT v pyDRTtools

**Autor:** Analýza kódu pyDRTtools
**Datum:** 2025-12-10
**Verze:** pyDRTtools 0.2

---

## Exekutivní Souhrn

Tento dokument poskytuje podrobnou technickou analýzu metody pro výpočet distribuce relaxačních časů (Distribution of Relaxation Times, DRT) z dat elektrochemické impedanční spektroskopie (EIS) implementované v balíčku pyDRTtools. Metoda využívá radiální bázové funkce (RBF) pro diskretizaci a Tikhonov regularizaci (ridge regression) s automatickým výběrem regularizačního parametru pomocí různých cross-validation metod.

**Klíčové charakteristiky:**
- Bázový přístup: RBF diskretizace s 8 typy bázových funkcí
- Regularizace: Ridge regression s penalizací 1. nebo 2. derivace
- Automatizace: 6 metod pro výběr optimálního λ (GCV, mGCV, rGCV, LC, re-im, k-fold)
- Omezení: Non-negativity constraint na γ(τ)
- Řešič: CVXOPT quadratic programming

---

## 1. Matematický Základ

### 1.1 Model Impedance

Impedance elektrochemického systému je modelována jako:

```
Z(ω) = R_∞ + jωL_0 + ∫_{-∞}^{∞} γ(ln τ)/(1 + jωτ) d(ln τ)
```

**Komponenty:**
- **R_∞**: Vysokofrekvenční ohmický odpor (solution resistance)
- **L_0**: Indukčnost (kabely, elektrody)
- **γ(τ)**: Distribution of Relaxation Times - hledaná funkce
- **τ**: Relaxační čas
- **ω = 2πf**: Úhlová frekvence

### 1.2 Separace na Reálnou a Imaginární Část

**Reálná část:**
```
Z_re(ω) = R_∞ + ∫_{-∞}^{∞} γ(ln τ)/(1 + ω²τ²) d(ln τ)
```

**Imaginární část:**
```
Z_im(ω) = ωL_0 - ∫_{-∞}^{∞} γ(ln τ)·ωτ/(1 + ω²τ²) d(ln τ)
```

### 1.3 RBF Diskretizace

Distribuce γ(ln τ) se aproximuje jako lineární kombinace radiálních bázových funkcí:

```
γ(ln τ) ≈ Σ_{i=1}^{N} x_i · φ(ln τ - ln τ_i)
```

kde **φ(·)** je zvolená RBF funkce a **x_i** jsou koeficienty k určení.

---

## 2. Radiální Bázové Funkce (RBF)

### 2.1 Dostupné Typy RBF

| RBF Typ | Definice | Vlastnosti |
|---------|----------|------------|
| **Gaussian** | `exp(-(εx)²)` | Nekonečně hladká, kompaktní podpora v praxi |
| **C0 Matern** | `exp(-ε\|x\|)` | C⁰ spojitá (nespojitá derivace) |
| **C2 Matern** | `exp(-ε\|x\|)(1+ε\|x\|)` | C² spojitá |
| **C4 Matern** | `exp(-ε\|x\|)(3+3ε\|x\|+(εx)²)/3` | C⁴ spojitá |
| **C6 Matern** | `exp(-ε\|x\|)(15+15ε\|x\|+6(εx)²+(εx)³)/15` | C⁶ spojitá |
| **Inverse Quadratic** | `1/(1+(εx)²)` | Globální vliv |
| **Inverse Quadric** | `1/√(1+(εx)²)` | Globální vliv |
| **Cauchy** | `1/(1+ε\|x\|)` | Těžké ocasové rozdělení |
| **Piecewise Linear** | Lineární interpolace | Nejjednodušší, nehladké |

### 2.2 Shape Factor (ε)

**Automatický výpočet** (`compute_epsilon`, basics.py:96-137):

```python
# 1. Výpočet FWHM koeficientu (Full Width at Half Maximum)
FWHM_coeff = 2 * fsolve(rbf - 0.5, initial_guess)[0]

# 2. Průměrný logaritmický rozestup
delta = mean(diff(log(1/freq)))

# 3. Shape factor
epsilon = coeff * FWHM_coeff / delta
```

**Význam:**
- `coeff` (typicky 0.5): Uživatelem definovaný multiplikátor
- Menší ε → širší RBF → více vyhlazenější DRT
- Větší ε → užší RBF → více detailní DRT

**Default konfigurace:** `coeff=0.5`, `shape_control='FWHM Coefficient'`

---

## 3. Sestavení Diskretizačních Matic

### 3.1 Matice A_re (Reálná Impedance)

**Dimenze:** `[N_freqs × N_taus]`

**Element A_re[p, q]** (`g_i`, basics.py:31-61):

```python
alpha = 2*π*freq_p*tau_q

integrand = lambda x: 1/(1 + alpha² * exp(2*x)) * φ(x)

A_re[p,q] = ∫_{-50}^{50} integrand(x) dx
```

Integrál se počítá numericky pomocí `scipy.integrate.quad` s tolerancí `epsabs=1e-9, epsrel=1e-9`.

**Speciální případ - Piecewise Linear:**
```python
if q == 0:
    A_re[p,q] = 0.5/(1+(ω_p*τ_q)²) * log(τ_{q+1}/τ_q)
elif q == N_taus-1:
    A_re[p,q] = 0.5/(1+(ω_p*τ_q)²) * log(τ_q/τ_{q-1})
else:
    A_re[p,q] = 0.5/(1+(ω_p*τ_q)²) * log(τ_{q+1}/τ_{q-1})
```

### 3.2 Matice A_im (Imaginární Impedance)

**Element A_im[p, q]** (`g_ii`, basics.py:64-94):

```python
integrand = lambda x: alpha/(1/exp(x) + alpha² * exp(x)) * φ(x)

A_im[p,q] = -∫_{-50}^{50} integrand(x) dx
```

### 3.3 Rozšíření o R_∞ a L_0

**Finální struktura** (runs.py:176-185):

```python
# Bez indukčnosti (N_RL = 1)
A_re_final = [1, 0, ..., 0]  [DRT část]
             [1, 0, ..., 0]
             [1, 0, ..., 0]
             └─R_∞ sloupec

A_im_final = [0,   2πf_1, ..., 0]  [DRT část]
             [0,   2πf_2, ..., 0]
             [0,   2πf_3, ..., 0]

# S indukčností (N_RL = 2)
# Přidá se sloupec [0, 0, ..., 0]^T pro R_∞ v A_im
# a sloupec [2πf_1, 2πf_2, ...]^T pro L_0 v A_im
```

### 3.4 Toeplitz Optimalizace

**Podmínka použití** (basics.py:372):
```python
std_diff_freq = std(diff(log(1/freq)))
mean_diff_freq = mean(diff(log(1/freq)))
toeplitz_trick = (std_diff_freq / mean_diff_freq < 0.01) and (N_freqs == N_taus)
```

**Implementace:**
```python
# Místo výpočtu celé matice element po elementu:
R = [g_i(freq[0], tau_k, ε, rbf_type) for k in range(N)]
C = [g_i(freq_k, tau[0], ε, rbf_type) for k in range(N)]
A_re = toeplitz(C, R)
```

**Výhoda:** Redukce složitosti z O(N²) na O(N).

---

## 4. Regularizační Matice M

### 4.1 Účel Regularizace

Regularizace penalizuje "drsnost" řešení γ(τ) a zajišťuje:
- Numerickou stabilitu (ill-posed problém)
- Vyhlazenost DRT spektra
- Redukci přeučení (overfitting)

### 4.2 Matice První Derivace (M₁)

**Dimenze:** `[N_taus × N_taus]`

**Element M₁[n, m]** (`inner_prod_rbf_1`, basics.py:142-193):

```
M₁[n,m] = ⟨dφ_n/dx, dφ_m/dx⟩ = ∫ φ'_n(x) · φ'_m(x) dx
```

**Pro Gaussian RBF:**
```python
a = ε * log(freq_n / freq_m)
M₁[n,m] = -ε * (-1 + a²) * exp(-a²/2) * sqrt(π/2)
```

**Pro C2 Matern:**
```python
M₁[n,m] = (ε/6) * (3 + 3|a| - |a|³) * exp(-|a|)
```

### 4.3 Matice Druhé Derivace (M₂)

**Element M₂[n, m]** (`inner_prod_rbf_2`, basics.py:196-248):

```
M₂[n,m] = ⟨d²φ_n/dx², d²φ_m/dx²⟩
```

**Pro Gaussian RBF:**
```python
M₂[n,m] = ε³ * (3 - 6a² + a⁴) * exp(-a²/2) * sqrt(π/2)
```

### 4.4 Speciální Případ: Piecewise Linear

**První derivace** (basics.py:485-494):
```python
L = zeros(N_taus-1, N_taus)
for i in range(N_taus-1):
    delta = log(tau[i+1]/tau[i])
    L[i, i] = -1/delta
    L[i, i+1] = 1/delta

M₁ = L^T @ L
```

**Druhá derivace** (basics.py:539-555):
```python
L = zeros(N_taus-2, N_taus)
for i in range(N_taus-2):
    delta = log(tau[i+1]/tau[i])
    L[i, i] = 1/delta²
    L[i, i+1] = -2/delta²
    L[i, i+2] = 1/delta²

M₂ = L^T @ L
```

### 4.5 Rozšíření Matice M

**Finální struktura pro R_∞ a L_0** (runs.py:184-185):
```python
M_final = [0, 0, ..., 0]  # R_∞ není penalizován
          [0, 0, ..., 0]  # L_0 není penalizován
          [0, 0,        ]
          [0, 0,   M    ]  # Pouze DRT část je penalizována
```

---

## 5. Optimalizační Problém

### 5.1 Ridge Regression Formulace

**Primární problém:**
```
minimize:  ||Z_exp - A·x||² + λ||M·x||²
subject to: x ≥ 0
```

kde:
- **Z_exp** = [Z_re; Z_im] (stacked vektor)
- **A** = [A_re; A_im] (stacked matice)
- **x** = [R_∞, L_0, γ₁, γ₂, ..., γ_N]^T
- **λ** > 0 regularizační parametr

### 5.2 Kvadratický Program

**Transformace na QP** (`quad_format_combined`, basics.py:592-616):

```
minimize: (1/2) x^T H x + c^T x
subject to: -x ≤ 0  (tj. x ≥ 0)
```

kde:
```python
H = 2 * (A_re^T @ A_re + A_im^T @ A_im + λ*M)
H = (H + H^T) / 2  # Vynucení symetrie

c = -2 * (Z_re^T @ A_re + Z_im^T @ A_im)
```

### 5.3 Řešení pomocí CVXOPT

**Implementace** (`solve_gamma`, basics.py:619-644):

```python
# Omezení: G @ x ≤ h
G = -I  # Identita matice
h = 0   # Nulový vektor

# Řešení
sol = cvxopt.solvers.qp(
    P=matrix(H),      # Kvadratický term
    q=matrix(c),      # Lineární term
    G=matrix(G),      # Nerovnostní matice
    h=matrix(h)       # Pravá strana
)

x = array(sol['x']).flatten()
```

**Output:** Vektor `x = [R_∞, L_0, γ₁, ..., γ_N]`

---

## 6. Automatický Výběr Regularizačního Parametru

### 6.1 Přehled Metod

| Metoda | Referenční Funkce | Výpočetní Složitost | Doporučení |
|--------|-------------------|---------------------|------------|
| **GCV** | `compute_GCV` | O(N³) | Obecné použití |
| **mGCV** | `compute_mGCV` | O(N³) | Malé datasety (N<50) |
| **rGCV** | `compute_rGCV` | O(N³) | Robustní k outlierům |
| **LC** | `compute_LC` | O(N³) | Vizuální interpretace |
| **re-im** | `compute_re_im_cv` | O(N³) × 2 | Využívá Re/Im rozdíl |
| **kf** | `compute_kf_cv` | O(k·N³) | Nejvíce obecné CV |

### 6.2 GCV - Generalized Cross Validation

**Score funkce** (parameter_selection.py:30-75):

```python
K = A @ inv(A^T @ A + λ*M) @ A^T

GCV(λ) = (1/n * ||‌(I - K)·Z||²) / (1/n * trace(I - K))²
```

**Interpretace:**
- Čitatel: Průměrná chyba rekonstrukce
- Jmenovatel: Efektivní počet stupňů volnosti (squared)
- Minimum: Optimální trade-off mezi fit a hladkostí

**Implementace pomocí Cholesky:**
```python
A_in = A^T @ A + λ*M

# Kontrola pozitivní definitnosti
if not is_PD(A_in):
    A_in = nearest_PD(A_in)

# Cholesky faktorizace
L = cholesky(A_in)
inv_A_in = inv(L)^T @ inv(L)

K = A @ inv_A_in @ A^T
```

### 6.3 mGCV - Modified GCV

**Score funkce** (parameter_selection.py:78-128):

```python
ρ = 1.3  if n_cv < 50 else 2.0

mGCV(λ) = (1/n * ||‌(I - K)·Z||²) / (1/n * trace(I - ρ·K))²
```

**Stabilizační parametr ρ:**
- Snižuje bias u malých datasetů
- Penalizuje příliš velké λ

### 6.4 rGCV - Robust GCV

**Score funkce** (parameter_selection.py:131-187):

```python
μ₂ = (1/n) * trace(K^T @ K)
ξ = 0.2  if n_cv < 100 else 0.3

rGCV(λ) = (ξ + (1-ξ)·μ₂) * GCV(λ)
```

**Robustní korekce:**
- `μ₂` měří "sílu" predikce
- `ξ` váhový parametr pro robustnost

### 6.5 LC - L-Curve Method

**Princip:**
- Vykresli: log(||Z - Ax||²) vs. log(λ||Mx||²)
- Hledej maximum zakřivení ("loket")

**Score funkce** (parameter_selection.py:398-480):

```python
η = log(||Z - Ax||²)
θ = log(λ||Mx||²)

# Derivace
η' = dη/d(log λ)
θ' = dθ/d(log λ)

# Zakřivení
κ = (η'·θ'' - θ'·η'') / (η'² + θ'²)^(3/2)

LC_score = -κ  # Maximalizuj zakřivení
```

### 6.6 re-im Cross Validation

**Algoritmus** (parameter_selection.py:190-333):

```python
# 1. Řeš z reálné části
x_re = solve_QP(A_re, Z_re, M, λ)
R_∞, γ_re = x_re[0], x_re[1:]

# 2. Řeš z imaginární části
x_im = solve_QP(A_im, Z_im, M, λ)
L_0, γ_im = x_im[0], x_im[1:]

# 3. Křížová predikce
γ_re_cv = [R_∞, γ_im]  # Re parametry + Im DRT
γ_im_cv = [0, γ_re]    # Im parametry + Re DRT

# 4. Score
score = ||Z_re - A_re @ γ_re_cv||² + ||Z_im - A_im @ γ_im_cv||²
```

**Výhoda:** Testuje konzistenci mezi Re a Im částí.

### 6.7 k-fold Cross Validation

**Algoritmus** (parameter_selection.py:336-395):

```python
kf = KFold(n_splits=5, shuffle=True)
score = 0

for train_idx, test_idx in kf.split(Z):
    # Train na train_idx
    x_train = solve_QP(A[train_idx], Z[train_idx], M, λ)

    # Test na test_idx
    Z_pred = A[test_idx] @ x_train
    score += ||Z[test_idx] - Z_pred||² / N_test

return score / n_splits
```

### 6.8 Optimalizace λ

**Společný framework** (`optimal_lambda`, basics.py:647-700):

```python
# Definuj bounds
bounds = [(log(1e-7), log(1.0))]

# Minimalizuj score funkci
result = minimize(
    fun=score_function,
    x0=log(lambda_init),
    args=(A_re, A_im, Z_re, Z_im, M, ...),
    method='SLSQP',
    bounds=bounds,
    options={'maxiter': 2000}
)

λ_optimal = exp(result.x)
```

**Logaritmická parametrizace:**
- Optimalizuje v `log(λ)` prostoru
- Zajišťuje λ > 0
- Lepší numerická stabilita

---

## 7. Numerické Aspekty a Optimalizace

### 7.1 Cholesky Dekompozice

**Motivace:**
- Přímá inverze matice: O(N³), numericky nestabilní
- Cholesky: O(N³/3), stabilnější pro symetrické pozitivně definitní matice

**Použití v kódu:**
```python
# Místo: inv(A^T @ A + λ*M)
L = cholesky(A^T @ A + λ*M)
inv_L = inv(L)
inv_ATA = inv_L^T @ inv_L
```

### 7.2 Kontrola Pozitivní Definitnosti

**Implementace** (nearest_PD.py):
```python
def is_PD(B):
    try:
        cholesky(B)
        return True
    except LinAlgError:
        return False

def nearest_PD(A):
    # Higham algoritmus pro nejbližší PD matici
    B = (A + A^T) / 2
    _, s, V = svd(B)
    H = V^T @ diag(s) @ V
    A2 = (B + H) / 2
    A3 = (A2 + A2^T) / 2

    if is_PD(A3):
        return A3

    # Iterativní korekce...
```

**Kdy se používá:**
- Numerické chyby v sestavení M
- Extrémní hodnoty λ
- Špatně podmíněné matice A

### 7.3 Symetrizace Matic

```python
H = 2 * (A_re^T @ A_re + A_im^T @ A_im + λ*M)
H = (H + H^T) / 2  # Vynucení symetrie
```

**Důvod:** CVXOPT vyžaduje přesně symetrické matice.

### 7.4 Integrace s Vysokou Přesností

```python
integrate.quad(
    func=integrand,
    a=-50, b=50,
    epsabs=1e-9,  # Absolutní tolerance
    epsrel=1e-9   # Relativní tolerance
)
```

**Limity integrace [-50, 50]:**
- V log-prostoru pokrývají ~exp(±50) ≈ 10^±22 řádů
- Prakticky nekonečné intervaly

### 7.5 Výpočetní Složitost

| Operace | Složitost | Poznámka |
|---------|-----------|----------|
| Sestavení A_re/A_im (brute force) | O(N²_freqs × N_taus) | Každý element = integrál |
| Sestavení A_re/A_im (Toeplitz) | O(N) | Když je log-spaced |
| Sestavení M | O(N²_taus) | Analytické vzorce |
| Cholesky dekompozice | O(N³/3) | Jednou za iteraci λ |
| QP řešení (CVXOPT) | O(N³) | Interior point metoda |
| **GCV evaluace** | **O(N³)** | Bottleneck |
| **Celková optimalizace λ** | **O(iter × N³)** | iter ≈ 10-20 |

**Pro N=100:**
- Čas na evaluaci GCV: ~0.1s
- Celková optimalizace: ~1-2s
- S Toeplitz: ~0.5-1s

---

## 8. Implementační Workflow

### 8.1 High-Level Pipeline

```
EIS Data (freq, Z_re, Z_im)
    ↓
EIS_object.from_file()
    ↓
┌─────────────────────────────────┐
│ 1. Compute epsilon              │ → compute_epsilon()
│ 2. Assemble A_re, A_im          │ → assemble_A_re/im()
│ 3. Assemble M (1st or 2nd der.) │ → assemble_M_1/2()
│ 4. Extend for R_∞, L_0          │ → np.hstack([...])
└─────────────────────────────────┘
    ↓
┌─────────────────────────────────┐
│ 5. Optimize λ                   │ → optimal_lambda()
│    - Choose CV method           │   (GCV/mGCV/rGCV/LC/re-im/kf)
│    - SLSQP minimization         │
└─────────────────────────────────┘
    ↓
┌─────────────────────────────────┐
│ 6. Solve QP                     │ → solve_gamma()
│    - Format H, c                │   quad_format_combined()
│    - CVXOPT QP solver           │   cvxopt.solvers.qp()
│    - Extract x                  │
└─────────────────────────────────┘
    ↓
┌─────────────────────────────────┐
│ 7. Post-process                 │
│    - Split: R_∞, L_0, γ         │
│    - Map to fine grid           │ → x_to_gamma()
│    - Reconstruct Z              │
└─────────────────────────────────┘
    ↓
DRT spectrum γ(τ), fitted Z(ω)
```

### 8.2 Příklad Použití

```python
from pyDRTtools.runs import EIS_object, simple_run

# 1. Načtení dat
data = EIS_object.from_file('eis_data.txt')
# Formát: freq  Z_real  Z_imag

# 2. Spuštění DRT analýzy
result = simple_run(
    entry=data,
    rbf_type='Gaussian',           # RBF typ
    data_used='Combined Re-Im Data', # Použij Re i Im
    induct_used=1,                  # Zahrnout indukčnost
    der_used='2nd order',           # M₂ regularizace
    cv_type='GCV',                  # GCV pro výběr λ
    reg_param=1e-3,                 # Použito jen když cv_type='custom'
    shape_control='FWHM Coefficient',
    coeff=0.5
)

# 3. Přístup k výsledkům
tau = result.out_tau_vec     # Časové konstanty
gamma = result.gamma         # DRT spektrum
R_inf = result.R_inf        # Ohmický odpor
L_0 = result.L_0            # Indukčnost
lambda_opt = result.lambda_value  # Optimální λ
Z_fit = result.Z_DRT        # Rekonstruovaná impedance

# 4. Vizualizace
result.plot_DRT()
```

### 8.3 Klíčové Funkce a jejich Role

```
pyDRTtools/
│
├── basics.py              ← Core DRT výpočty
│   ├── g_i()             → Element A_re
│   ├── g_ii()            → Element A_im
│   ├── compute_epsilon() → Shape factor
│   ├── assemble_A_re()   → Sestaví A_re
│   ├── assemble_A_im()   → Sestaví A_im
│   ├── assemble_M_1/2()  → Regularizační matice
│   ├── quad_format_*()   → QP formulace
│   ├── solve_gamma()     → CVXOPT wrapper
│   └── optimal_lambda()  → λ optimalizace
│
├── parameter_selection.py ← CV metody
│   ├── compute_GCV()
│   ├── compute_mGCV()
│   ├── compute_rGCV()
│   ├── compute_LC()
│   ├── compute_re_im_cv()
│   └── compute_kf_cv()
│
├── runs.py               ← High-level API
│   ├── EIS_object        → Data container
│   ├── simple_run()      → Ridge regression
│   ├── Bayesian_run()    → Bayesian s HMC
│   └── BHT_run()         → Hilbert transform
│
└── nearest_PD.py         ← Numerická stabilita
    ├── is_PD()
    └── nearest_PD()
```

---

## 9. Parametry a Jejich Vliv

### 9.1 Parametry RBF Diskretizace

| Parametr | Rozsah | Default | Vliv |
|----------|--------|---------|------|
| `rbf_type` | 8 typů | 'Gaussian' | Tvar bázové funkce |
| `coeff` | 0.1-2.0 | 0.5 | Šířka RBF (menší = širší) |
| `epsilon` | Auto | - | Shape factor (odvozeno z coeff) |
| `N_taus` | 20-200 | N_freqs | Počet kolokačních bodů |

**Doporučení:**
- **Gaussian**: Univerzální, dobře funguje pro většinu systémů
- **C2 Matern**: Když jsou očekávány ostré peaky
- **Piecewise Linear**: Nejrychlejší, ale nejméně hladké

### 9.2 Regularizační Parametry

| Parametr | Volby | Default | Vliv |
|----------|-------|---------|------|
| `der_used` | '1st order', '2nd order' | '2nd order' | Stupeň vyhlazenosti |
| `cv_type` | GCV/mGCV/rGCV/LC/re-im/kf/custom | 'GCV' | Metoda výběru λ |
| `reg_param` | 1e-7 - 1.0 | Auto | Manuální λ (když cv_type='custom') |

**Doporučení:**
- **2nd order**: Silnější vyhlazení, doporučeno pro šumová data
- **1st order**: Méně vyhlazenosti, zachová ostřejší struktury
- **GCV**: Nejrobustnější pro obecné použití
- **re-im**: Dobrá volba když Re i Im mají podobnou kvalitu

### 9.3 Parametry Dat

| Parametr | Volby | Default | Použití |
|----------|-------|---------|---------|
| `data_used` | 'Combined Re-Im', 'Re Data', 'Im Data' | 'Combined Re-Im' | Která část impedance |
| `induct_used` | 0/1/2 | 1 | 0:bez L₀, 1:s L₀, 2:vyřaď induktivní data |

---

## 10. Validace a Testování

### 10.1 Kontrolní Body

**Numerická stabilita:**
```python
# 1. Pozitivní definitnost
assert is_PD(A^T @ A + λ*M), "Matice není PD!"

# 2. Symetrie
assert np.allclose(H, H.T), "H není symetrická!"

# 3. Kondice matice
cond_number = np.linalg.cond(H)
assert cond_number < 1e10, f"Špatně podmíněná: {cond_number}"
```

**Fyzikální konzistence:**
```python
# 1. Non-negativita γ
assert np.all(gamma >= 0), "Negativní DRT!"

# 2. Kladný R_∞
assert R_inf > 0, "Negativní odpor!"

# 3. Kvalita fitu
residual = np.linalg.norm(Z_exp - Z_fit) / np.linalg.norm(Z_exp)
assert residual < 0.1, f"Velká chyba fitu: {residual:.2%}"
```

### 10.2 Syntetická Data

**Generování testovacích dat:**
```python
# Dva Voigt elementy v sérii
def synthetic_impedance(freq, R_inf, R1, tau1, R2, tau2):
    omega = 2*np.pi*freq
    Z1 = R1 / (1 + 1j*omega*tau1)
    Z2 = R2 / (1 + 1j*omega*tau2)
    return R_inf + Z1 + Z2

# Přidání šumu
Z_noisy = Z_true + noise_level * (randn(N) + 1j*randn(N))
```

**Validační metriky:**
```python
# 1. Peak detection
from scipy.signal import find_peaks
peaks, _ = find_peaks(gamma, height=0.01*max(gamma))
assert len(peaks) == 2, "Nedetekoval oba peaky!"

# 2. Integrovaný odpor
R_pol = np.trapz(gamma, np.log(tau))
R_total = R_inf + R_pol
assert abs(R_total - R_expected) < 0.01*R_expected
```

---

## 11. Omezení a Problémy

### 11.1 Známá Omezení

**1. Non-negativity Constraint**
```
Problém: γ(τ) ≥ 0 je vynuceno
Důsledek: Nelze zachytit induktivní smyčky v DRT prostoru
Řešení: Pokročilejší metody (např. unconstrained fit + interpretace)
```

**2. Ill-posed Problém**
```
Problém: Malé změny v Z → velké změny v γ
Důsledek: Citlivost na šum, potřeba regularizace
Řešení: Automatický výběr λ pomocí CV metod
```

**3. Výběr N_taus**
```
Problém: Příliš málo → špatné rozlišení; příliš mnoho → overfitting
Doporučení: N_taus ≈ N_freqs
```

**4. Hranice τ**
```
Problém: DRT mimo rozsah měřených frekvencí je nespolehlivá
Řešení: Omezit interpretaci na τ ∈ [1/(2πf_max), 1/(2πf_min)]
```

### 11.2 Numerické Problémy

**Špatně podmíněné matice:**
```python
# Příznaky
cond(A^T @ A) > 1e10  # Velmi vysoké číslo kondice

# Příčiny
- Velmi úzké RBF (vysoké epsilon)
- Nerovnoměrný spacing frekvencí
- Extrémně malé/velké λ

# Řešení
- Použít nearest_PD() korekci
- Zvážit Tikhonov regularizaci s větším λ
- Změnit RBF typ nebo coeff
```

**Selhání CVXOPT:**
```python
# Symptom
sol['status'] != 'optimal'

# Možné příčiny
- Nesymetrická H matice
- Nesprávné bounds
- Numericky singulární problém

# Debug
print("H symetrie:", np.linalg.norm(H - H.T))
print("H eigenvalues:", np.linalg.eigvals(H)[:5])
```

### 11.3 Fyzikální Interpretace

**Pozor na:**
1. **Artefakty na okrajích:** DRT na τ << 1/(2πf_max) nebo τ >> 1/(2πf_min) může být nefyzikální
2. **Over-smoothing:** Příliš velké λ → sloučení separátních procesů
3. **Under-smoothing:** Příliš malé λ → rozpad jednoho procesu na více peaků
4. **Pseudo-peaky:** Numerické artefakty vypadající jako fyzikální procesy

---

## 12. Srovnání s Alternativními Metodami

### 12.1 pyDRTtools vs. Equivalent Circuit Fitting

| Aspekt | pyDRTtools (DRT) | EC Fitting |
|--------|------------------|------------|
| **Model** | Model-free | Vyžaduje předpoklad o obvodech |
| **Flexibilita** | Vysoká | Nízká (musí znát topologii) |
| **Interpretace** | Přímá (peaky = procesy) | Nepřímá (R, C, Q parametry) |
| **Robustnost** | Střední (závisí na λ) | Vysoká (když je model správný) |
| **Výpočetní čas** | Střední (O(N³)) | Rychlá (O(N) iterací) |
| **Overfitting** | Kontrolováno regularizací | Kontrolováno počtem elementů |

### 12.2 Ridge vs. LASSO vs. Bayesian

| Metoda | Regularizace | Sparse γ? | Uncertainty? | Čas |
|--------|--------------|-----------|--------------|-----|
| **Ridge (pyDRTtools)** | L2 (M norm) | Ne | Ne | O(N³) |
| **LASSO** | L1 | Ano | Ne | O(N³) iterace |
| **Bayesian (pyDRTtools)** | Hierarchical prior | Ne | Ano (credible intervals) | O(N³ × samples) |

**pyDRTtools implementace:**
- `simple_run()`: Ridge regression (tento report)
- `Bayesian_run()`: HMC sampling pro credible intervals
- `BHT_run()`: Bayesian Hilbert Transform pro validaci

---

## 13. Best Practices

### 13.1 Příprava Dat

```python
# 1. Kontrola kvality dat
plt.plot(Z_re, -Z_im, 'o')  # Nyquist plot
# Zkontroluj: uzavřené oblouky, žádné skoky

# 2. Odstranění induktivních dat (pokud nevhodné)
mask = Z_im < 0
freq, Z = freq[mask], Z[mask]

# 3. Log-spacing kontrola
delta_log_f = np.diff(np.log(freq))
print(f"Log spacing CV: {np.std(delta_log_f)/np.mean(delta_log_f):.3f}")
# < 0.01 → Toeplitz optimalizace aktivní
```

### 13.2 Volba Parametrů

```python
# 1. Start s defaulty
result = simple_run(data, cv_type='GCV', rbf_type='Gaussian')

# 2. Zkontroluj fit
plt.subplot(121)
plt.plot(freq, Z_exp.real, 'o', label='Exp')
plt.plot(freq, result.Z_DRT.real, '-', label='Fit')
plt.subplot(122)
plt.plot(freq, Z_exp.imag, 'o')
plt.plot(freq, result.Z_DRT.imag, '-')

# 3. Když je špatný fit → zkus jiné cv_type
for method in ['GCV', 'mGCV', 'rGCV', 're-im']:
    r = simple_run(data, cv_type=method)
    print(f"{method}: λ={r.lambda_value:.2e}, residual={...}")

# 4. Když je over/under-smoothed → manuální λ
result = simple_run(data, cv_type='custom', reg_param=1e-2)
```

### 13.3 Interpretace DRT

```python
# 1. Identifikace peaků
from scipy.signal import find_peaks
peaks_idx, properties = find_peaks(gamma, height=0.05*max(gamma))
tau_peaks = tau[peaks_idx]
gamma_peaks = gamma[peaks_idx]

# 2. Charakterizace procesů
for i, (t, g) in enumerate(zip(tau_peaks, gamma_peaks)):
    f_char = 1/(2*np.pi*t)
    R_process = np.trapz(gamma[near peak], np.log(tau[near peak]))
    print(f"Process {i}: f={f_char:.2e} Hz, R={R_process:.2f} Ω")

# 3. Celkový polarizační odpor
R_pol = np.trapz(gamma, np.log(tau))
print(f"R_pol = {R_pol:.2f} Ω")
print(f"R_total = R_inf + R_pol = {result.R_inf + R_pol:.2f} Ω")
```

### 13.4 Validace Výsledků

```python
# 1. Residuální analýza
residual = Z_exp - Z_fit
residual_norm = np.abs(residual) / np.abs(Z_exp)

plt.semilogx(freq, 100*residual_norm, 'o')
plt.ylabel('Relative Error (%)')
plt.axhline(1, color='r', linestyle='--', label='1%')

# 2. Kramers-Kronig check (pokud dostupné)
from impedance.validation import linKK
M, mu, Z_kk, res_real, res_imag = linKK(freq, Z_exp)
# res < 1% → data jsou konzistentní

# 3. Porovnání různých metod
results = {}
for rbf in ['Gaussian', 'C2 Matern', 'C4 Matern']:
    results[rbf] = simple_run(data, rbf_type=rbf)

# Plot overlay
for name, res in results.items():
    plt.semilogx(res.out_tau_vec, res.gamma, label=name)
# Podobné tvary → robustní závěr
```

---

## 14. Závěry

### 14.1 Silné Stránky Metody

✅ **Model-free přístup:** Nevyžaduje a priori znalost počtu relaxačních procesů

✅ **Automatizace:** Cross-validation metody automaticky vybírají optimální regularizaci

✅ **Numerická stabilita:** Cholesky dekompozice, PD korekce, Toeplitz optimalizace

✅ **Flexibilita RBF:** 8 typů bázových funkcí pro různé aplikace

✅ **Dobře zdokumentováno:** Kód navazuje na 10+ peer-reviewed publikací

✅ **Open-source:** MIT licence, aktivní vývoj

### 14.2 Oblasti pro Zlepšení

⚠️ **Performance:** O(N³) složitost limituje použití pro N > 200

⚠️ **Non-negativity:** Vynucení γ ≥ 0 může skrýt některé fyzikální jevy

⚠️ **Interpretace:** Vyžaduje expertní znalost EIS pro správnou interpretaci DRT

⚠️ **Dokumentace kódu:** Chybí docstringy u některých funkcí

⚠️ **Testy:** Absence unit testů, pouze tutorial notebooky

### 14.3 Použitelnost

**Doporučeno pro:**
- Exploratorní analýzu EIS dat
- Identifikaci počtu relaxačních procesů
- Přípravu pro equivalent circuit fitting
- Srovnání různých experimentálních podmínek

**Nedoporučeno pro:**
- Real-time analýzu (příliš pomalé)
- Data s extrémním šumem (bez robustní varianty)
- Systémy s dominantní induktivitou

### 14.4 Budoucí Směry

1. **GPU akcelerace** CVXOPT QP solveru pro N > 500
2. **Adaptive RBF** s automatickým výběrem epsilon pro různé τ oblasti
3. **Unconstrained DRT** pro zachycení induktivních smyček
4. **Uncertainty quantification** bez full Bayesian samplingingu
5. **Integration** s Kramers-Kronig validací v jednom workflow

---

## 15. Reference

### Klíčové Publikace (v kódu citované)

[1] T.H. Wan, M. Saccoccio, C. Chen, F. Ciucci, *Influence of the discretization methods on the distribution of relaxation times deconvolution: Implementing radial basis functions with DRTtools*, Electrochimica Acta **184** (2015) 483-499.

[2] M. Saccoccio, T.H. Wan, C. Chen, F. Ciucci, *Optimal regularization in distribution of relaxation times applied to electrochemical impedance spectroscopy: Ridge and lasso regression methods*, Electrochimica Acta **147** (2014) 470-482.

[3] J. Liu, T.H. Wan, F. Ciucci, *A Bayesian view on the Hilbert transform and the Kramers-Kronig transform of electrochemical impedance data*, Electrochimica Acta **357** (2020) 136864.

[4] A. Maradesa, B. Py, T.H. Wan, M.B. Effat, F. Ciucci, *Selecting the regularization parameter in the distribution of relaxation times*, Journal of the Electrochemical Society **170** (2023) 030502.

### Cross-Validation Reference

[5] G. Wahba, *A comparison of GCV and GML for choosing the smoothing parameter*, Annals of Statistics **13** (1985) 1378-1402.

[6] P.C. Hansen, D.P. O'Leary, *The use of the L-curve in the regularization of discrete ill-posed problems*, SIAM Journal on Scientific Computing **14** (1993) 1487-1503.

### Implementace

[7] CVXOPT: https://cvxopt.org/
[8] pyDRTtools GitHub: https://github.com/ciuccislab/pyDRTtools
[9] pyDRTtools Manual: pyDRTtools_manual.pdf

---

## Přílohy

### A. Kompletní Seznam Parametrů simple_run()

```python
simple_run(
    entry,                    # EIS_object instance
    rbf_type='Gaussian',     # 'Gaussian'|'C0 Matern'|'C2 Matern'|'C4 Matern'|
                             # 'C6 Matern'|'Inverse Quadratic'|'Inverse Quadric'|
                             # 'Cauchy'|'PWL'
    data_used='Combined Re-Im Data',  # 'Combined Re-Im Data'|'Re Data'|'Im Data'
    induct_used=1,           # 0: bez L_0 | 1: s L_0 | 2: discard inductive data
    der_used='2nd order',    # '1st order'|'2nd order'
    cv_type='GCV',           # 'GCV'|'mGCV'|'rGCV'|'LC'|'re-im'|'kf'|'custom'
    reg_param=1e-3,          # Manuální λ (jen když cv_type='custom')
    shape_control='FWHM Coefficient',  # 'FWHM Coefficient'|'Shape Factor'
    coeff=0.5                # FWHM koeficient nebo přímý shape factor
)
```

### B. Typické Hodnoty a Rozsahy

| Veličina | Typický Rozsah | Jednotky | Poznámka |
|----------|----------------|----------|----------|
| freq | 10⁻³ - 10⁶ | Hz | EIS měřicí rozsah |
| τ | 10⁻⁶ - 10³ | s | τ = 1/(2πf) |
| R_∞ | 0.1 - 100 | Ω | Závisí na systému |
| L_0 | 10⁻⁹ - 10⁻⁶ | H | Typicky nH-μH |
| γ(τ) | 0 - 100 | Ω | Nerestrikované |
| λ_opt | 10⁻⁵ - 10⁻¹ | - | Z CV optimalizace |
| epsilon | 1 - 100 | - | Auto-computed |

### C. Diagnostické Grafy

**1. Nyquist Plot s Fitem**
```python
plt.plot(Z_exp.real, -Z_exp.imag, 'o', label='Experimental')
plt.plot(Z_fit.real, -Z_fit.imag, '-', label='DRT Fit')
plt.xlabel('$Z_{re}$ / Ω')
plt.ylabel('$-Z_{im}$ / Ω')
plt.axis('equal')
```

**2. DRT Spektrum**
```python
plt.semilogx(tau, gamma, 'k-', linewidth=2)
plt.xlabel(r'$\tau$ / s')
plt.ylabel(r'$\gamma(\tau)$ / Ω')
plt.xlim([1/(2*π*freq.max()), 1/(2*π*freq.min())])
```

**3. Residuální Plot**
```python
res = np.abs(Z_exp - Z_fit) / np.abs(Z_exp) * 100
plt.semilogx(freq, res, 'o')
plt.ylabel('Relative Error (%)')
plt.axhline(1, color='r', linestyle='--')
```

**4. L-Curve (když cv_type='LC')**
```python
lambdas = np.logspace(-7, 0, 50)
residuals, regularizations = [], []
for lam in lambdas:
    x = solve_gamma(A_re, A_im, Z_re, Z_im, M, lam)
    residuals.append(np.linalg.norm(Z - A@x))
    regularizations.append(lam * np.linalg.norm(M@x))

plt.loglog(residuals, regularizations, 'o-')
plt.xlabel('||Z - Ax||')
plt.ylabel('λ||Mx||')
```

---

**Konec reportu**

*Vygenerováno z analýzy pyDRTtools source code, verze 0.2 (2024)*
