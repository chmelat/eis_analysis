# Lin-KK implementace v eis_analysis

## Prehled

Modul `eis_analysis.validation.kramers_kronig` poskytuje nativni implementaci Lin-KK testu (Linear Kramers-Kronig) pro validaci kvality EIS dat. Implementace je zalozena na metode Schonleber et al. (2014) a nevyzaduje externi zavislosti.

**Reference:**
- Schonleber, M. et al. "A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
- Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
- Yrjana, V. and Bobacka, J. "Implementing Kramers-Kronig validity testing using pyimpspec." Electrochim. Acta 504, 144951 (2024)

---

## 1. Teoreticky zaklad

Lin-KK test fituje impedancni data pomoci rady Voigtovych elementu s pevnymi casovymi konstantami:

```
Z_fit(w) = R_s + sum_{k=1}^{M} R_k / (1 + jw*tau_k) + jw*L
```

Kde:
- `R_s` = seriovy odpor [Ohm]
- `R_k` = odpor k-teho Voigtova elementu [Ohm]
- `tau_k` = casova konstanta (fixni, logaritmicky rozlozena) [s]
- `L` = seriova induktance [H]
- `M` = pocet elementu (optimalizovan pomoci mu metriky)

**Klicove:** Casove konstanty tau_k jsou fixni. Fitovany jsou pouze odpory R_k a seriove komponenty (R_s, L).

---

## 2. Distribuce casovych konstant

### Zakladni rozsah

```python
tau_min = 1 / (2*pi*f_max)
tau_max = 1 / (2*pi*f_min)
```

Casove konstanty jsou logaritmicky rozlozeny v intervalu [tau_min, tau_max]:

```
tau_k = 10^[log10(tau_min) + (k-1)/(M-1) * log10(tau_max/tau_min)]
```

### Rozsireni rozsahu (extend_decades)

Parametr `extend_decades` umoznuje rozsirit rozsah tau o zadany pocet dekad:

```python
tau_min_extended = tau_min * 10^(-extend_decades)
tau_max_extended = tau_max * 10^(+extend_decades)
```

- `extend_decades = 0.0`: zadne rozsireni (standardni Lin-KK)
- `extend_decades = 0.4`: rozsireni o 0.4 dekady na kazde strane
- `extend_decades = 1.0`: rozsireni o 1 dekadu na kazde strane

---

## 3. Mu metrika a volba poctu elementu

### Algoritmus

1. Zacni s M=3 Voigtovymi elementy
2. Fituj pomoci pseudoinverze (povol zaporne R_k)
3. Vypocitej mu metriku
4. Pokud mu > threshold, zvys M a opakuj
5. Kdyz mu <= threshold, optimalni M nalezeno

### Mu metrika (Schonleber 2014)

```
mu = 1 - (sum|R_k| pro R_k < 0) / (sum|R_k| pro R_k >= 0)
```

**Interpretace:**
- mu -> 1: Vsechny R_k kladne, dobry fit
- mu -> 0.85: Prah pro zastaveni (vychozi hodnota)
- mu < 0.85: Prilis mnoho zapornych R_k, overfit

**Proc zaporne R_k znamenaji overfit?**
- Fyzikalne by mely byt vsechny R_k >= 0 (odpor nemuze byt zaporny)
- Zaporne R_k vznikaji, kdyz model ma prilis mnoho volnosti
- Model zacina fitovat sum misto skutecneho signalu

---

## 4. Metriky kvality

### Pseudo chi-squared (Boukamp 1995)

```
chi2_ps = sum_i w_i * [(Z'_exp - Z'_fit)^2 + (Z''_exp - Z''_fit)^2]
```

kde `w_i = 1/|Z_i|^2` je Boukampova vaha.

### Odhad sumu (Yrjana & Bobacka 2024)

```
noise_estimate = sqrt(chi2_ps * 5000 / N) [%]
```

kde N je pocet bodu.

### Rezidua

```
res_real = (Z'_exp - Z'_fit) / |Z_exp|
res_imag = (Z''_exp - Z''_fit) / |Z_exp|
```

**Interpretace:**
- |res| < 1%: Vynikajici kvalita dat
- |res| < 5%: Prijatelna kvalita
- |res| >= 5%: Mozne artefakty, nelinearita, nestabilita

### Vazeni pri fittingu

Pouzivame `proportional` vazeni (w = 1/|Z|) podle Schonleber (2014), nikoli `modulus` (w = 1/|Z|^2) podle Boukamp (1995):

- **1/|Z|** - vyrovnane relativni vazeni pres cele spektrum
- **1/|Z|^2** - silny duraz na vysoke frekvence (nizke |Z|), nizkofrekvencni oblast ma maly vliv

Pro typicka EIS data je 1/|Z| vhodnejsi, protoze nizkofrekvencni oblast casto obsahuje klicove elektrochemicke informace (prenos naboje, difuze) a bezne artefakty (drift, nestacionarita).

**Pozn.:** Pseudo chi^2 se pocita s vahou 1/|Z|^2 podle Boukampa - rozdil je pouze ve fittingu, ne ve vysledne metrice.

---

## 5. Automaticka optimalizace extend_decades

### Motivace

Optimalni rozsah casovych konstant zavisi na datech. Prilis uzky rozsah muze vest k vysokym reziduum na okrajich frekvencniho spektra, zatimco prilis siroky rozsah zvysuje pocet parametru.

### Implementace

Funkce `find_optimal_extend_decades()` pouziva grid search:

```python
def find_optimal_extend_decades(
    frequencies, Z, M,
    search_range=(-1.0, 1.0),
    n_evaluations=11
):
    """
    Hleda optimalni extend_decades minimalizaci pseudo chi^2.

    Algoritmus:
    1. Vytvor mrizku hodnot v search_range
    2. Pro kazdy bod: spocitej fit a chi^2
    3. Vrat hodnotu s minimalnim chi^2
    """
```

### CLI pouziti

```bash
eis data.dta --auto-extend
```

Vystup:
```
Optimizing extend_decades in range (-1.0, 1.0)...
Optimal extend_decades: 0.400
Lin-KK native: M=36, mu=0.8377
  extend_decades: 0.400
  Mean |res_real|: 0.03%
  Mean |res_imag|: 4.81%
  Pseudo chi^2: 3.79e-01
  Estimated noise: 4.40%
```

---

## 6. API reference

### KKResult dataclass

```python
@dataclass
class KKResult:
    M: int                           # Pocet Voigtovych elementu
    mu: float                        # Mu metrika
    Z_fit: NDArray[np.complex128]    # Fitovana impedance
    residuals_real: NDArray          # Rezidua realne casti
    residuals_imag: NDArray          # Rezidua imaginarni casti
    pseudo_chisqr: float             # Pseudo chi^2 (Boukamp 1995)
    noise_estimate: float            # Odhad sumu [%]
    extend_decades: float            # Pouzite rozsireni tau
    inductance: Optional[float]      # Induktance [H]
    figure: Optional[plt.Figure]     # Vizualizace

    # Convenience properties
    @property
    def mean_residual_real(self) -> float: ...  # Prumer |res_real| [%]
    @property
    def mean_residual_imag(self) -> float: ...  # Prumer |res_imag| [%]
    @property
    def is_valid(self) -> bool: ...             # True pokud rezidua < 5%
```

### kramers_kronig_validation()

```python
def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = False
) -> Optional[KKResult]:
    """
    Provede KK validaci EIS dat.

    Parameters
    ----------
    frequencies : array
        Merene frekvence [Hz]
    Z : array
        Komplexni impedance [Ohm]
    mu_threshold : float
        Prah pro mu metriku (default: 0.85)
    max_M : int
        Maximalni pocet Voigtovych elementu (default: 50)
    auto_extend_decades : bool
        Automaticka optimalizace extend_decades (default: False)

    Returns
    -------
    KKResult nebo None pri chybe
    """
```

### lin_kk_native()

```python
def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    fit_type: str = 'real',
    weighting: str = 'proportional',
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (-1.0, 1.0)
) -> Tuple[int, float, NDArray, NDArray, NDArray, Optional[float], float, float]:
    """
    Nativni Lin-KK implementace s plnou kontrolou parametru.

    Returns
    -------
    M, mu, Z_fit, res_real, res_imag, L_value, chi2_ps, extend_decades
    """
```

### Pomocne funkce

```python
def compute_pseudo_chisqr(Z_exp, Z_fit) -> float:
    """Vypocet pseudo chi^2 (Boukamp 1995)."""

def estimate_noise_percent(chi2_ps, n_points) -> float:
    """Odhad sumu z pseudo chi^2 (Yrjana & Bobacka 2024)."""

def reconstruct_impedance(frequencies, elements, tau, L_value, include_L) -> NDArray:
    """Rekonstrukce impedance z Voigtovych elementu."""

def find_optimal_extend_decades(frequencies, Z, M, ...) -> Tuple[...]:
    """Nalezeni optimalniho extend_decades grid search."""
```

---

## 7. Priklady pouziti

### Python API

```python
from eis_analysis.validation import kramers_kronig_validation

# Zakladni validace
result = kramers_kronig_validation(frequencies, Z)
print(f"M = {result.M}, mu = {result.mu:.3f}")
print(f"Pseudo chi^2: {result.pseudo_chisqr:.2e}")
print(f"Estimated noise: {result.noise_estimate:.1f}%")
print(f"Valid: {result.is_valid}")

# S automatickou optimalizaci extend_decades
result = kramers_kronig_validation(
    frequencies, Z,
    auto_extend_decades=True
)
print(f"Optimal extend_decades: {result.extend_decades:.3f}")
```

### CLI

```bash
# Zakladni KK validace
eis data.dta --no-drt --no-fit

# S automatickou optimalizaci extend_decades
eis data.dta --no-drt --no-fit --auto-extend

# Pouze KK validace (bez vizualizace)
eis data.dta --no-show --no-drt --no-fit --auto-extend
```

---

## 8. Porovnani metod fitu

### fit_type = 'real' (vychozi)

1. Fituje pouze Re(Z/|Z|) -> ziska [R_s, R_1, ..., R_M]
2. Vypocita Z_fit z techto parametru
3. Fituje residuum Im(Z - Z_fit) -> ziska L

**Vyhody pro validaci:**
- Realna cast je obvykle mene zasumena
- Imaginarni rezidua slouzi jako diagnostika KK compliance
- Pokud data splnuji KK relace, imaginarni cast by mela automaticky sedet

### fit_type = 'complex'

Fituje Re(Z) a Im(Z) soucasne. **Nedoporuceno pro validaci** - muze maskovat KK nekonzistence.

---

## 9. Interpretace vysledku

### Priklad: Dobra data

```
Lin-KK native: M=36, mu=0.8377
  extend_decades: 0.400
  Mean |res_real|: 0.03%
  Mean |res_imag|: 4.81%
  Pseudo chi^2: 3.79e-01
  Estimated noise: 4.40%
Data quality is good (residuals < 5%)
```

- Nizke rezidua v realne casti (0.03%) - vynikajici shoda
- Vyssi rezidua v imaginarni casti (4.81%) - prijatelne, pod prahem 5%
- Odhadovany sum 4.40% - odpovida ocekavani pro typicka EIS data

### Priklad: Problematicka data

```
Lin-KK native: M=32, mu=0.8436
  extend_decades: 0.400
  Mean |res_real|: 0.18%
  Mean |res_imag|: 13.34%
  Pseudo chi^2: 2.44e+00
  Estimated noise: 11.17%
! Data may contain artifacts (residuals >= 5%)
```

- Vysoke rezidua v imaginarni casti (13.34%) - nad prahem 5%
- Vysoky odhadovany sum (11.17%)
- Mozne priciny: nestacionarita, artefakty, nelinearita

---

## 10. Literatura

1. Schonleber, M., Klotz, D., and Ivers-Tiffee, E. (2014). A Method for Improving the Robustness of linear Kramers-Kronig Validity Tests. Electrochim. Acta, 131, 20-27.

2. Boukamp, B.A. (1995). A Linear Kronig-Kramers Transform Test for Immittance Data Validation. J. Electrochem. Soc., 142, 1885-1894.

3. Yrjana, V. and Bobacka, J. (2024). Implementing Kramers-Kronig validity testing using pyimpspec. Electrochim. Acta, 504, 144951.
