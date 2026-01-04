# Lineární fit Voigt chainu - matematický popis

## 1. Model impedance

### 1.1 Voigt element (K element)

Voigt element je paralelní kombinace rezistoru R a kapacitoru C:

```
Z_K(omega) = R / (1 + j*omega*tau)
```

kde `tau = R*C` je časová konstanta.

Rozklad na reálnou a imaginární část:

```
Re[Z_K] = R / (1 + (omega*tau)^2)
Im[Z_K] = -R*omega*tau / (1 + (omega*tau)^2)
```

### 1.2 Voigt chain (řetězec K elementů)

Kompletní model impedance:

```
Z(omega) = R_s + sum_{k=1}^{M} R_k / (1 + j*omega*tau_k) + j*omega*L
```

kde:
- `R_s` - sériový odpor
- `M` - počet Voigt elementů
- `R_k, tau_k` - parametry k-tého Voigt elementu
- `L` - sériová indukčnost (artefakt měření)

Rozklad:

```
Re[Z] = R_s + sum_{k=1}^{M} R_k / (1 + (omega*tau_k)^2)

Im[Z] = sum_{k=1}^{M} -R_k*omega*tau_k / (1 + (omega*tau_k)^2) + omega*L
```

## 2. Linearizace problému

### 2.1 Klíčová myšlenka

Pro **fixní** hodnoty `tau_k` je impedance **lineární** v parametrech `[R_s, R_1, ..., R_M, L]`.

Definujeme bázové funkce:

```
h_k(omega) = 1 / (1 + (omega*tau_k)^2)       pro Voigt elementy
g_k(omega) = -omega*tau_k / (1 + (omega*tau_k)^2)   pro imaginární část
```

Pak:

```
Re[Z] = R_s * 1 + sum_k R_k * h_k(omega)
Im[Z] = sum_k R_k * g_k(omega) + L * omega
```

### 2.2 Maticový zápis

Pro N frekvenčních bodů `omega_1, ..., omega_N` sestavíme matice:

**Matice pro reálnou část A_re** (N x (1 + M + 1)):

```
A_re = | 1   h_1(omega_1)  h_2(omega_1)  ...  h_M(omega_1)   0        |
       | 1   h_1(omega_2)  h_2(omega_2)  ...  h_M(omega_2)   0        |
       | :        :             :                  :         :        |
       | 1   h_1(omega_N)  h_2(omega_N)  ...  h_M(omega_N)   0        |
```

**Matice pro imaginární část A_im** (N x (1 + M + 1)):

```
A_im = | 0   g_1(omega_1)  g_2(omega_1)  ...  g_M(omega_1)   omega_1  |
       | 0   g_1(omega_2)  g_2(omega_2)  ...  g_M(omega_2)   omega_2  |
       | :        :             :                  :           :      |
       | 0   g_1(omega_N)  g_2(omega_N)  ...  g_M(omega_N)   omega_N  |
```

**Vektor parametrů:**

```
x = [R_s, R_1, R_2, ..., R_M, L]^T
```

**Soustava rovnic:**

```
A_re @ x = Re[Z_data]
A_im @ x = Im[Z_data]
```

## 3. Vážení bodů (Weighting)

### 3.1 Motivace

EIS data mají typicky velký dynamický rozsah (|Z| může být 10 Ohm až 10000 Ohm).
Bez vážení by velké impedance dominovaly fit.

### 3.2 Váhové schéma

Definujeme váhy `w_i` pro každý frekvenční bod:

| Schéma | Váha w_i | Popis |
|--------|----------|-------|
| uniform | 1 | Všechny body stejně |
| sqrt | 1/sqrt(\|Z_i\|) | Kompromis |
| proportional | 1/\|Z_i\| | Lin-KK standard |
| square | \|Z_i\|^2 | Důraz na vysoké Z |

**Normalizace:** Váhy jsou normalizovány tak, aby `mean(w) = 1`.

### 3.3 Aplikace vah

Vážené matice a vektory:

```
A_re_weighted[i, :] = w_i * A_re[i, :]
A_im_weighted[i, :] = w_i * A_im[i, :]

b_re_weighted[i] = w_i * Re[Z_data[i]]
b_im_weighted[i] = w_i * Im[Z_data[i]]
```

## 4. Řešení lineárního problému

### 4.1 Typy fitu

**fit_type = 'real':**
Fituje pouze reálnou část, L extrahuje z imaginárního residua.

```
min_x ||A_re_w @ x - b_re_w||^2
```

Pak L z residua:
```
L = pinv(omega * w) @ (Im[Z_data] - A_im[:, :-1] @ x[:-1]) * w
```

**fit_type = 'complex':**
Fituje obě části současně.

```
min_x ||A_re_w @ x - b_re_w||^2 + ||A_im_w @ x - b_im_w||^2
```

Ekvivalentně:
```
A_combined = [A_re_w]    b_combined = [b_re_w]
             [A_im_w]                 [b_im_w]

min_x ||A_combined @ x - b_combined||^2
```

### 4.2 Metody řešení

**NNLS (Non-Negative Least Squares):**
- Vynucuje `R_k >= 0` (fyzikální omezení)
- Používá scipy.optimize.nnls
- Default pro `allow_negative=False`

```
x = argmin_x ||A @ x - b||^2   s.t. x >= 0
```

**Pseudoinverze:**
- Povoluje záporné `R_k` (indikátor overfitu)
- Používá numpy.linalg.pinv
- Pro `allow_negative=True` (Lin-KK styl)

```
x = pinv(A) @ b = (A^T A)^{-1} A^T b
```

## 5. Generování mřížky tau

### 5.1 Logaritmické rozložení

Časové konstanty jsou rozloženy logaritmicky přes frekvenční rozsah:

```
tau_min = 1 / (2*pi*f_max)
tau_max = 1 / (2*pi*f_min)

tau_k = 10^(log10(tau_min) + (k-1)/(M-1) * log10(tau_max/tau_min))
```

### 5.2 Hustota bodů

Parametr `n_per_decade` určuje počet tau hodnot na dekádu:

```
M = n_per_decade * (log10(f_max) - log10(f_min)) + 1
```

Typické hodnoty: 2-3 body na dekádu.

### 5.3 Rozšíření rozsahu

Parametr `extend_decades` rozšiřuje tau směrem k nižším frekvencím:

```
f_min_extended = f_min / 10^(extend_decades)
```

## 6. Pruning (ořezávání malých elementů)

### 6.1 Kombinovaná strategie

Po fitu jsou malé R_k odstraněny pomocí kombinovaného thresholdu:

```
threshold_relative = prune_threshold * max(|R_k|)
threshold_absolute = 0.001 * sum(|R_k|)

threshold_effective = min(threshold_relative, threshold_absolute)

keep_mask = |R_k| >= threshold_effective
```

### 6.2 Motivace

- Relativní threshold: odstraní elementy zanedbatelné vůči největšímu
- Absolutní threshold: chrání malé ale významné elementy když jeden dominuje

## 7. Metrika mu pro optimalizaci M

### 7.1 Definice (Schönleber et al. 2014)

```
mu = 1 - sum(|R_k| for R_k < 0) / sum(|R_k| for R_k > 0)
```

### 7.2 Interpretace

- `mu = 1.0`: Všechny R_k >= 0, žádný overfit
- `mu < 1.0`: Některé R_k < 0, indikace overfitu
- `mu < 0.85`: Typický threshold pro "dostatečný" počet elementů

### 7.3 Algoritmus optimalizace M

```
M = 3
while mu > mu_threshold and M < max_M:
    M += 1
    tau = generate_tau_grid_fixed_M(freq, M)
    x = solve_linear(freq, Z, tau, allow_negative=True)
    mu = calc_mu(x)
return M
```

## 8. Výstup

### 8.1 Obvod

Výsledný obvod má strukturu:

```
R(R_s) - K(R_1, tau_1) - K(R_2, tau_2) - ... - K(R_M, tau_M) - L(L)
```

### 8.2 Parametry

Vektor parametrů pro nelineární optimalizaci:

```
params = [R_s, R_1, tau_1, R_2, tau_2, ..., R_M, tau_M, L]
```

## Reference

1. Schönleber, M. et al. "A Method for Improving the Robustness of linear
   Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)

2. Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance
   Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
