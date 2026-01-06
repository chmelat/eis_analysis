# Plan: Skew-Normal Mixture Model pro DRT Peak Detection

**Status:** Naplanovano, ceka na implementaci

## Cil

Implementovat alternativni peak detection metodu s **asymetrickymi piky** v log(tau) prostoru.

## Motivace

Soucasny GMM model predpoklada symetricke piky (Gaussiany) v log(tau) prostoru. Realne relaxacni procesy vsak casto vykazuji asymetricke distribuce casovych konstant - napr. difuzni procesy s dlouhym ocasem k vyssim tau.

## Matematicke pozadi

### Soucasny stav: GMM

```
GMM fituje: log10(tau) ~ Normal(mu, sigma^2)
Piky jsou symetricke v log(tau) zobrazeni
```

### Novy model: Skew-Normal Mixture

Skew-normal distribuce:
```
f(x; xi, omega, alpha) = (2/omega) * phi((x-xi)/omega) * Phi(alpha * (x-xi)/omega)
```

kde:
- xi = lokace (podobne jako mu)
- omega = skala (podobne jako sigma)
- alpha = parametr asymetrie (-inf, +inf)
  - alpha = 0: symetricky (Gauss)
  - alpha > 0: pravy ocas delsi (k vyssim tau)
  - alpha < 0: levy ocas delsi (k nizsim tau)
- phi = PDF normalni distribuce
- Phi = CDF normalni distribuce

### Proc skew-normal?

1. **Fyzikalni motivace**: Realne relaxacni procesy casto maji asymetrickou distribuci casovych konstant
2. **Zobecneni GMM**: alpha=0 redukuje na Gaussian, zpetna kompatibilita
3. **scipy.stats.skewnorm**: Existujici implementace, neni treba psat od nuly
4. **Interpretovatelnost**: alpha primo ukazuje smer a miru asymetrie

## Rozhodnuti

1. **Implementacni pristup**: scipy.optimize s L-BFGS-B
2. **Alpha**: Individualni pro kazdy pik (4n-1 parametru)
3. **Fallback**: Kdyz skew-normal fit selze -> fallback na GMM (alpha=0)

## Implementace

### Soubory k modifikaci

1. `eis_analysis/drt/peaks.py` - nova funkce `skew_normal_peak_detection()`
2. `eis_analysis/drt/core.py` - pridani `--peak-method skew` volby
3. `eis.py` - CLI parametr

### Nova funkce: skew_normal_peak_detection()

```python
def skew_normal_peak_detection(
    tau: NDArray[np.float64],
    gamma: NDArray[np.float64],
    n_components_range: Tuple[int, int] = (1, 6),
    bic_threshold: float = 10.0
) -> Tuple[List[Dict], Optional[object], List[float]]:
    """
    Detekuje asymetricke piky v DRT spektru pomoci Skew-Normal Mixture.

    Rozsiruje GMM o parametr asymetrie (alpha) pro kazdy pik.
    """
```

### Algoritmus fitu

```python
from scipy.stats import skewnorm
from scipy.optimize import minimize

def fit_skew_mixture(X, weights, n_components):
    """
    Fit skew-normal mixture pomoci scipy.optimize.

    Parametry pro kazdy pik: (xi, omega, alpha, weight)
    Celkem 4*n - 1 parametru (-1 protoze vahy se scitaji na 1)
    """

    def neg_log_likelihood(params):
        # Rozbal parametry
        components = unpack_params(params, n_components)

        # Mixture log-likelihood
        ll = 0
        for x, w in zip(X, weights):
            p = sum(comp_w * skewnorm.pdf(x, a, loc=xi, scale=omega)
                    for xi, omega, a, comp_w in components)
            ll += w * np.log(max(p, 1e-300))
        return -ll

    # Initial guess z GMM
    initial_guess = get_initial_from_gmm(X, weights, n_components)

    # Bounds
    bounds = get_bounds(n_components)

    result = minimize(neg_log_likelihood, initial_guess,
                      method='L-BFGS-B', bounds=bounds)
    return result
```

### BIC pro model selection

```python
# Skew-normal mixture ma 4*n - 1 parametru
k = 4 * n_components - 1
bic = -2 * log_likelihood + k * np.log(n_samples)
```

### Vystupy

Rozsireni GMM vystupu o asymetrii:
```python
peaks.append({
    'tau_center': tau_center,
    'tau_bounds': (tau_lower, tau_upper),  # Asymetricke!
    'log_tau_std': omega,
    'skewness': alpha,  # NOVE
    'weight': weight,
    'R_estimate': R_estimate,
    'f_center': f_center
})
```

### CLI integrace

```bash
# Nova volba
eis data.DTA --peak-method skew

# Volby peak-method:
# - scipy (default): rychle, bez modelu
# - gmm: symetricke piky, BIC
# - skew: asymetricke piky, BIC (NOVE)
```

## Odhad slozitosti

- peaks.py: nova funkce ~150 radku
- core.py: male zmeny ~20 radku
- eis.py: CLI parametr ~5 radku
- Dokumentace: aktualizace GMM_PEAK_DETECTION.md

## Zavislosti

- scipy.stats.skewnorm (uz je v scipy, zadna nova zavislost)
- numpy, scipy.optimize

## Testovani

1. Synteticka data s asymetrickymi piky
2. Porovnani s GMM na symetrickych datech (alpha ~ 0)
3. Realna EIS data s difuznimi procesy

## Reference

- Azzalini, A. (1985). "A class of distributions which includes the normal ones". Scandinavian Journal of Statistics.
- scipy.stats.skewnorm dokumentace
