# Kriticke srovnani: Z-HIT implementace

**Projekt**: eis_analysis vs pyimpspec
**Datum**: 2026-01-09 (aktualizace)
**Autor**: Analyza provedena pomoci Claude Code

---

## Prehled

Srovnani Z-HIT (Z-Hilbert Impedance Transform) validace v:
- **eis_analysis** (beta verze, 337 radku)
- **pyimpspec** (produkcni verze, ~1500 radku)

---

## 1. Matematicky zaklad

### Spolecny princip

Obe implementace vychazi z Hilbertovy transformace:

```
ln|Z(w)| = ln|Z(w_ref)| + (2/pi) * integral[phi(w) d(ln w)] - (pi/6) * d(phi)/d(ln w)
```

Korekce druheho radu `-(pi/6) * d(phi)/d(ln w)` pochazi z Ehm et al. (2001).

**Poznamka k implementaci:** Kod pouziva `gamma = -pi/6` s plusem (`+ gamma * derivative`), coz je matematicky ekvivalentni.

### Srovnani implementace

| Aspekt | eis_analysis | pyimpspec |
|--------|--------------|-----------|
| Integrace | Lichobezniky v ln(w) | Spline integrator |
| Derivace | Centralni diference | Spline derivative |
| Korekce 2. radu | Ano (-pi/6) | Ano (-pi/6) |
| Komplexni rezidua | Ano (Re + Im) | Ano (Re + Im) |
| Pseudo chi^2 | Ano | Ano |
| Odhad sumu | Ano | Ano |
| Optimalizace offsetu | Ano (vahova funkce) | Ano (lmfit) |

Obe implementace pouzivaji stejnou korekci `-(pi/6) * d(phi)/d(ln w)`.

---

## 2. Smoothing fazovych dat

### eis_analysis

**Zadne smoothing** - pouziva surova fazova data.

Jedina moznost: `edge_padding` pro okrajove efekty:
- 'reflect' (default)
- 'constant'
- 'none'

### pyimpspec

**5 metod smoothingu:**

| Metoda | Popis | Parametry |
|--------|-------|-----------|
| none | Bez smoothingu | - |
| savgol | Savitzky-Golay filtr | num_points, polynomial_order |
| lowess | LOWESS | num_points, num_iterations |
| modsinc | Modifikovany sinc kernel (Schmid 2023) | degree (2,4,6,8,10), bandwidth |
| whithend | Whittaker-Henderson | polynomial_order, num_points |

**Auto mode:** Testuje vsechny metody a vybere nejlepsi (min chi^2).

### Hodnoceni

**Chybi v eis_analysis:**
- Smoothing fazovych dat je kriticke pro zasumena data
- Bez smoothingu mohou byt rekonstruovane hodnoty nestabilni
- **Doporuceni: Pridat alespon Savitzky-Golay filtr**

---

## 3. Interpolace faze

### eis_analysis

**Zadna interpolace** - pouziva diskretni body primo.

Integrace pomoci lichobezniku mezi sousednimi body.

### pyimpspec

**4 spline metody:**

| Metoda | Popis |
|--------|-------|
| akima | Akima spline (neoscilatni) |
| makima | Modifikovany Akima (default) |
| cubic | Kubicky spline |
| pchip | PCHIP (zachovava monotonii) |

Interpolator poskytuje:
- `integrate(a, b)` - numericka integrace
- `derivative()` - prvni derivace

**Auto mode:** Testuje vsechny metody.

### Hodnoceni

**Chybi v eis_analysis:**
- Spline interpolace zlepsuje presnost integrace
- Hladsi derivace pro korekci 2. radu
- **Doporuceni: Pridat Akima nebo PCHIP interpolaci**

---

## 4. Offset (referencni bod)

### eis_analysis

**Dva rezimy:**

1. **Pevny referencni bod** (default):
```python
ref_freq = sqrt(f_min * f_max)  # geometricky prumer
ref_idx = argmin(|frequencies - ref_freq|)
ln_Z_ref = ln(|Z[ref_idx]|)
```

2. **Optimalizovany offset** (`--zhit-optimize-offset`):
```python
# Gaussovo okno v log-frekvencnim prostoru
weights = exp(-0.5 * ((log10(f) - center) / sigma)^2)
offset = sum(weights * (ln_Z_exp - ln_Z_fit)) / sum(weights)
```

Parametry:
- `offset_center`: stred okna (default: 1.5 → ~31.6 Hz)
- `offset_width`: sirka okna v dekadach (default: 3.0)

### pyimpspec

**Optimalizovany offset:**
```python
# lmfit minimalizace
offset = minimize(sum(weights * (ln_Z_fit + offset - ln_Z_exp)^2))
```

Parametry:
- `center`: stred okna (default: log10(31.6 Hz))
- `width`: sirka okna (default: 3 dekady)
- `window`: okenna funkce (hann, hamming, boxcar, ...)
- `weights`: vlastni vahy

### Hodnoceni

Obe knihovny podporuji optimalizaci offsetu s vahami. eis_analysis pouziva analyticky vypocet, pyimpspec numerickou optimalizaci (lmfit).

---

## 5. Kvalitativni metriky

### eis_analysis

```python
# Amplitudova rezidua
residuals_mag = (Z_mag - Z_mag_reconstructed) / Z_mag * 100  # [%]

# Komplexni rezidua (normalizovane modulem)
residuals_real = (Z.real - Z_fit.real) / |Z|
residuals_imag = (Z.imag - Z_fit.imag) / |Z|

# Pseudo chi^2 (Boukamp 1995)
pseudo_chisqr = sum((Z_exp - Z_fit)^2 / |Z_exp|^2)

# Odhad sumu (Yrjana & Bobacka 2024)
noise_estimate = sqrt(pseudo_chisqr * 5000 / N)

# Kvalita (0-1 skala)
quality = max(0, 1 - mean(|residuals_mag|) / 5.0)
```

**Interpretace:**
- quality > 0.9: Vynikajici
- quality > 0.7: Dobre
- quality > 0.5: Prijatelne
- quality = 0: Rezidua >= 5%

**ZHITResult properties:**
- `mean_residual_real`, `mean_residual_imag`, `mean_residual_mag`
- `is_valid` (True pokud mean_residual_mag < 5%)

### pyimpspec

```python
pseudo_chisqr = sum((Z_exp - Z_fit)^2 / |Z_exp|^2)

residuals_real = (Z_real_exp - Z_real_fit) / |Z_exp|
residuals_imag = (Z_imag_exp - Z_imag_fit) / |Z_exp|
```

**Metriky:**
- Pseudo chi^2 (Boukamp 1995)
- Separatni rezidua pro realnou a imaginarni cast
- Rezidua normalizovana modulem (ne komponentou)

### Hodnoceni

Obe knihovny poskytuji stejne metriky:
- Pseudo chi^2 pro kvantitativni hodnoceni
- Komplexni rezidua (Re + Im) pro diagnostiku
- Odhad sumu z pseudo chi^2

**Poznamka:** Z-HIT odhad sumu je horni hranice (zahrnuje chybu numericke integrace).

---

## 6. Paralelizace a vykon

### eis_analysis

**Sekvencni zpracovani** - jednoducha implementace.

### pyimpspec

**Paralelni zpracovani:**
```python
multiprocessing.Pool.imap_unordered()
```

- Auto mode testuje kombinace paralelne
- `num_procs` parametr pro kontrolu

### Hodnoceni

Pro jednoduchou validaci neni paralelizace nutna. Relevantni pouze pro "auto" mode.

---

## 7. API a integrace

### eis_analysis

```python
def zhit_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    ref_freq: Optional[float] = None,
    quality_threshold: float = 5.0,
    optimize_offset: bool = False,
    offset_center: float = 1.5,
    offset_width: float = 3.0
) -> ZHITResult

@dataclass
class ZHITResult:
    Z_mag_reconstructed: NDArray[np.float64]  # Rekonstruovana amplituda
    Z_fit: NDArray[np.complex128]             # Rekonstruovana komplexni Z
    residuals_mag: NDArray[np.float64]        # Amplitudova rezidua [%]
    residuals_real: NDArray[np.float64]       # Realna rezidua (fraction)
    residuals_imag: NDArray[np.float64]       # Imaginarni rezidua (fraction)
    pseudo_chisqr: float                      # Pseudo chi^2
    noise_estimate: float                     # Odhad sumu [%]
    quality: float                            # Kvalita (0-1)
    ref_freq: float                           # Referencni frekvence [Hz]
    figure: Optional[plt.Figure]              # Vizualizace
    # Properties: mean_residual_real, mean_residual_imag, mean_residual_mag, is_valid
```

**CLI:**
- Z-HIT bezi automaticky spolu s Lin-KK (default)
- `--no-zhit`: vypnout Z-HIT
- `--zhit-optimize-offset`: pouzit optimalizaci offsetu

### pyimpspec

```python
def perform_zhit(
    data: DataSet,
    smoothing="modsinc",
    interpolation="makima",
    window="auto",
    num_points=3,
    polynomial_order=2,
    num_iterations=3,
    center=1.5,
    width=3.0,
    weights=None,
    admittance=False,
    num_procs=-1
) -> ZHITResult
```

**ZHITResult dataclass:**
- frequencies, impedances, residuals
- pseudo_chisqr
- smoothing, interpolation, window (pouzite metody)
- get_nyquist_data(), get_bode_data(), get_residuals_data()
- to_statistics_dataframe()

### Hodnoceni

Obe knihovny poskytuji strukturovany vysledek (dataclass) s komplexnimi metrikami.

**Rozdily:**
- pyimpspec: vice smoothing/interpolace metod, auto mode
- eis_analysis: jednodussi API, analyticky vypocet offsetu

---

## 8. Souhrnne hodnoceni

### Silne stranky eis_analysis

1. **Jednoduchost** - snadne pochopeni a pouziti
2. **Rychlost** - bez zbytecne komplexity
3. **Korekce 2. radu** - spravne implementovana
4. **CLI integrace** - Z-HIT bezi automaticky spolu s Lin-KK
5. **Kompletni metriky** - pseudo chi^2, komplexni rezidua, odhad sumu
6. **Strukturovany vysledek** - ZHITResult dataclass s convenience properties
7. **Optimalizace offsetu** - analyticky vypocet s Gaussovym oknem

### Slabe stranky eis_analysis (potencialni rozsireni)

| Priorita | Funkce | Popis |
|----------|--------|-------|
| Stredni | Smoothing faze | Savitzky-Golay pro zasumena data |
| Stredni | Spline interpolace | Akima/PCHIP pro lepsi presnost |
| Nizka | Auto mode | Automaticky vyber parametru |

### Silne stranky pyimpspec

1. **Robustnost** - 5 smoothing + 4 interpolace metod
2. **Auto mode** - automaticky vyber nejlepsi kombinace
3. **Flexibilita** - mnoho parametru
4. **Kompletni metriky** - pseudo chi^2, komplexni rezidua

### Slabe stranky pyimpspec

1. **Komplexita** - obtiznejsi pochopeni
2. **Zavislosti** - statsmodels, lmfit
3. **Overhead** - pomalejsi pro jednoduche pouziti

---

## 9. Doporuceni pro eis_analysis

### Implementovano (od v0.10.0)

Nasledujici funkce byly implementovany od puvodni verze tohoto reportu:

1. **ZHITResult dataclass** - strukturovany vysledek s convenience properties
2. **Pseudo chi^2** - `compute_pseudo_chisqr()` (Boukamp 1995)
3. **Komplexni rezidua** - `residuals_real`, `residuals_imag`
4. **Odhad sumu** - `estimate_noise_percent()` (Yrjana & Bobacka 2024)
5. **Optimalizace offsetu** - `--zhit-optimize-offset` s Gaussovym oknem
6. **Rekonstruovana komplexni Z** - `Z_fit = |Z_recon| * exp(j*phi)`

### Potencialni rozsireni (stredni priorita)

1. **Smoothing fazovych dat**
   ```python
   from scipy.signal import savgol_filter
   phi_smooth = savgol_filter(phi, window_length, polyorder)
   ```
   Pro zasumena data muze zlepsit stabilitu rekonstrukce.

2. **Spline interpolace**
   ```python
   from scipy.interpolate import Akima1DInterpolator
   interp = Akima1DInterpolator(ln_omega, phi)
   integral = interp.integrate(a, b)
   ```
   Zlepseni presnosti integrace a hladssi derivace.

### Nizka priorita

3. **Auto mode**
   - Automaticky vyber parametru testovanim kombinaci
   - Pro eis_analysis zatim neni potreba (jednodussi pristup funguje dobre)

---

## 10. Zaver

Od verze 0.10.0 je eis_analysis implementace Z-HIT kompletni a produkcne pripravena:

- **Korektni matematicky zaklad** - Hilbertova transformace s korekci 2. radu (Ehm et al. 2001)
- **Kompletni metriky** - pseudo chi^2, komplexni rezidua, odhad sumu
- **Strukturovany vysledek** - ZHITResult dataclass s convenience properties
- **Optimalizace offsetu** - analyticky vypocet s vahovym oknem

Hlavni rozdily oproti pyimpspec:

| Aspekt | eis_analysis | pyimpspec |
|--------|--------------|-----------|
| Smoothing | Zadny (raw data) | 5 metod + auto |
| Interpolace | Lichobezniky | 4 spline metody |
| Offset | Analyticky (Gauss) | lmfit optimalizace |
| Slozitost | Jednoduche | Komplexni |

Pro typicka EIS data je eis_analysis plne dostatecna. Smoothing a spline interpolace mohou byt pridany jako potencialni rozsireni pro zasumena data, ale nejsou kriticke pro bezne pouziti.

---

## Reference

1. Ehm, W., Kaus, R., Schiller, C.A., Strunz, W. (2001). The evaluation of electrochemical impedance spectra using a modified logarithmic Hilbert transform. Journal of Electroanalytical Chemistry 499, 216-225.

2. Boukamp, B.A. (1995). A Linear Kronig-Kramers Transform Test for Immittance Data Validation. J. Electrochem. Soc. 142, 1885-1894.

3. Schmid, M. (2023). Modified sinc kernel smoothing for impedance spectroscopy. ACS Measurement Science Au.
