# Generování syntetických dat s Warburgovým chováním

Dokumentace pro `generate_warburg_data.py` - nástroj pro simulaci EIS dat vykazujících difúzní impedanci.

## 📋 Přehled

Script umožňuje generování realistických EIS dat pro systémy s:
- **Warburgovým chováním** (difúzní impedance)
- **Volitelným šumem** (0-10%)
- **Vlastním rozlišením** (body na dekádu)
- **Dvěma typy obvodů** (jednoduchý Randles, složitý)
- **Dvěma typy Warburgu** (semi-infinite, finite)

Export do **Gamry DTA** a/nebo **CSV** formátu kompatibilního s `eis_v2.py`.

## 🚀 Rychlý start

### Základní použití

```bash
# Jednoduchý Randles obvod, 1% šum, 10 bodů/dekáda
python3 generate_warburg_data.py output.csv --format csv

# S vizualizací
python3 generate_warburg_data.py output.csv --format csv --plot
```

### Doporučený workflow

```bash
# 1. Generuj data
python3 generate_warburg_data.py warburg_data.csv --format csv --noise 1.5 --ppd 12 --verbose

# 2. Analyzuj pomocí eis_v2.py
# POZNÁMKA: Pro Warburg data NEPOUŽÍVEJ --auto-lambda (viz sekce DRT Analýza níže)
python3 eis_v2.py warburg_data.csv --lambda 1.0 --peak-method gmm -v
```

## 📖 Parametry

### Povinné

| Parametr | Popis | Příklad |
|----------|-------|---------|
| `output` | Výstupní soubor | `warburg_data.csv` |

### Volitelné

| Parametr | Zkratka | Default | Popis |
|----------|---------|---------|-------|
| `--format` | `-f` | `dta` | Formát výstupu: `dta`, `csv`, `both` |
| `--circuit` | `-c` | `randles` | Typ obvodu: `randles`, `complex` |
| `--warburg-type` | `-w` | `semi-infinite` | Typ Warburgu: `semi-infinite`, `finite` |
| `--noise` | `-n` | `1.0` | Úroveň šumu [%] |
| `--ppd` | | `10` | Počet bodů na dekádu |
| `--f-min` | | `0.01` | Minimální frekvence [Hz] |
| `--f-max` | | `1e5` | Maximální frekvence [Hz] |
| `--seed` | | `None` | Random seed (reprodukovatelnost) |
| `--plot` | `-p` | `False` | Zobrazit grafy |
| `--verbose` | `-v` | `False` | Verbose výstup |

## 🔬 Typy obvodů

### 1. Randles (`--circuit randles`)

**Obvod:** `R_s-p(R_ct,C_dl)-W`

**Název:** Pojmenováno po **J.E.B. Randlesovi** (1947), který tento model popsal. Je to **klasický, standardní ekvivalentní obvod** v elektrochemii pro systémy s difúzí. Není to náhodná volba (random), ale konkrétní, dobře definovaný obvod!

Typický pro elektrochemické systémy s jedním elektrodovým procesem a difúzí.

**Parametry:**
- R_s = 10 Ω (solution resistance)
- R_ct = 100 Ω (charge transfer resistance)
- C_dl = 20 µF (double layer capacitance)
- σ = 50 Ω·s^(-1/2) (Warburg coefficient)

**Použití:**
```bash
python3 generate_warburg_data.py randles_data.csv --format csv --circuit randles
```

**Charakteristika:**
- Jeden polokruh v Nyquist diagramu
- 45° přímka při nízkých frekvencích (Warburg)
- Fáze → -45° při f → 0

### 2. Complex (`--circuit complex`)

**Obvod:** `R_s-p(R_ct1,C_dl1)-p(R_ct2,C_dl2)-W`

**Název:** "Complex" je prostě **popisný název** pro složitější obvod s více procesy. Není to standardní elektrochemická terminologie (na rozdíl od "Randles"). Přesnější by bylo "dual-process" nebo "two-RC-Warburg", ale pro jednoduchost používáme "complex".

Systém s více elektrodovými procesy (různé časové konstanty) a difúzí.

**Parametry:**
- R_s = 15 Ω
- R_ct1 = 50 Ω, C_dl1 = 10 µF (rychlejší proces)
- R_ct2 = 150 Ω, C_dl2 = 50 µF (pomalejší proces)
- σ = 30 Ω·s^(-1/2)

**Použití:**
```bash
python3 generate_warburg_data.py complex_data.csv --format csv --circuit complex
```

**Charakteristika:**
- Dva překrývající se polokruhy
- 45° přímka při nízkých frekvencích
- Vhodný pro testování GMM peak detection

## ⚡ Typy Warburg elementu

### Semi-infinite (`--warburg-type semi-infinite`)

**Model:** Neomezená difúze (tlustá vrstva elektrolytu)

**Impedance:**
```
Z_W = σ/√ω - jσ/√ω = σ(1-j)/√ω
```

**Charakteristika:**
- 45° přímka v celém difúzním režimu
- Typický pro systémy s velkým objemem elektrolytu
- Default volba

**Použití:**
```bash
python3 generate_warburg_data.py data.csv --format csv --warburg-type semi-infinite
```

### Finite (`--warburg-type finite`)

**Model:** Omezená difúze (bounded, tenká vrstva)

**Impedance:**
```
Z_W = R_W · tanh(√(jωτ_W)) / √(jωτ_W)
```

**Charakteristika:**
- Přechod z 45° na kapacitní chování při nízkých frekvencích
- Typický pro tenké elektrolytické vrstvy, membrány
- Saturace při f → 0

**Použití:**
```bash
python3 generate_warburg_data.py data.csv --format csv --warburg-type finite
```

## 📊 Příklady použití

### 1. Základní simulace pro testování

```bash
# Čistá data bez šumu
python3 generate_warburg_data.py clean_data.csv --format csv --noise 0 --plot

# S realistickým šumem
python3 generate_warburg_data.py noisy_data.csv --format csv --noise 2.0 --plot
```

### 2. Vysoké rozlišení pro publikace

```bash
# 15 bodů/dekáda, široký frekvenční rozsah
python3 generate_warburg_data.py high_res.csv \
    --format csv \
    --ppd 15 \
    --f-min 0.001 \
    --f-max 1e6 \
    --noise 0.5 \
    --verbose
```

### 3. Reprodukovatelné simulace

```bash
# Použití seed pro identické výsledky
python3 generate_warburg_data.py reproducible.csv \
    --format csv \
    --noise 1.5 \
    --seed 42
```

### 4. Export do obou formátů

```bash
# Vytvoří data.DTA i data.csv
python3 generate_warburg_data.py data --format both --circuit complex --ppd 12
```

### 5. Testování různých úrovní šumu

```bash
# Generuj sérii s různým šumem
for noise in 0 0.5 1.0 2.0 5.0; do
    python3 generate_warburg_data.py "noise_${noise}.csv" \
        --format csv \
        --noise $noise \
        --seed 123 \
        --ppd 12
done

# Analyzuj všechny
for f in noise_*.csv; do
    python3 eis_v2.py "$f" --auto-lambda --save "${f%.csv}" --no-show -v
done
```

### 6. Validace DRT analýzy

```bash
# Generuj známý systém
python3 generate_warburg_data.py validation.csv \
    --format csv \
    --circuit randles \
    --noise 1.0 \
    --ppd 15 \
    --verbose \
    --plot

# Analyzuj s GMM
python3 eis_v2.py validation.csv \
    --auto-lambda \
    --peak-method gmm \
    --auto-circuit \
    -v
```

## 📈 Interpretace výstupů

### Konzolový výstup

```
✓ CSV data exportována do: test.csv
  Počet bodů: 70
  Frekvenční rozsah: 1.00e-02 - 1.00e+05 Hz
  |Z| rozsah: 10.05 - 365.23 Ω

Statistiky dat (--verbose):
  Z' rozsah: [10.05, 310.45] Ω
  -Z'' rozsah: [0.12, 195.67] Ω
  |Z| rozsah: [10.05, 365.23] Ω
  Phase rozsah: [-55.3, -0.8] °

  Warburg charakteristika (f < 1 Hz):
    Průměrná fáze: -20.5° (ideální -45° pro čistý Warburg)
```

### Co sledovat

| Parametr | Očekávaná hodnota | Význam |
|----------|-------------------|--------|
| Fáze při VF | ~ 0° | Čistě rezistivní (R_s + R_ct) |
| Fáze při NF (f < 1 Hz) | ~ -20° až -45° | Warburg + kapacitní |
| -Z'' maximum | Odpovídá R_ct/2 | Časová konstanta RC |
| Průměrná fáze < 1 Hz | -20° až -30° | Kombinace Warburg + C_dl |

**Poznámka:** Průměrná fáze není čistých -45° protože:
1. C_dl dává kapacitní příspěvek (-90°)
2. Warburg dává -45°
3. Výsledek je jejich kombinace

## ⚠️ DRT Analýza Warburg Dat - Důležité Upozornění

### Problém s --auto-lambda

**NEPOUŽÍVEJ `--auto-lambda` pro Warburg data!** GCV (Generalized Cross Validation) má zásadní problém s difúzními systémy.

### Proč GCV selhává?

```
Warburg:  Z_W(ω) = σ(1-j)/√ω  →  γ(τ) ∝ 1/√τ  (diverguje k ∞)
RC pík:   Z_RC ~ R/(1+jωτ)    →  γ(τ) = ostrý pík
```

**Co se stane:**
1. Warburg vytváří **širokou distribuci** až do nekonečna
2. Difúzní komponenta **dominuje reziduální chybě** v celém rozsahu frekvencí
3. GCV optimalizuje celkovou chybu → volí **příliš malé λ**
4. Warburg ocas zahlcuje spektrum, **RC píky mizí**

### Typické příznaky

Když použiješ `--auto-lambda` na Warburg data:
```
✗ DRT spektrum dominuje divergující ocas při dlouhých τ
✗ RC píky jsou potlačené nebo pod detekčním prahem
✗ Nalezené λ je velmi malé (< 1e-4)
✗ Peak detection nachází jen 0-1 píky (nebo artefakty)
✗ První (správný) RC pík je jen mírně větší než druhý (artefakt)
```

### ✅ Správný postup pro Warburg data

```bash
# 1. MÍSTO auto-lambda použij manuální λ = 1.0 až 2.0
python3 eis_v2.py warburg_data.DTA --lambda 1.0 --peak-method gmm -v

# 2. Zkus rozsah hodnot a vyber nejlepší vizuálně
for lambda in 0.5 1.0 1.5 2.0; do
    python3 eis_v2.py warburg_data.DTA --lambda $lambda --save "lambda_$lambda" --no-show
done

# 3. S GMM peak detection pro robustní detekci RC píků
python3 eis_v2.py warburg_data.DTA --lambda 1.5 --peak-method gmm -v --plot

# 4. Případně omeź frekvenční rozsah při generování
python3 generate_warburg_data.py data.DTA --f-min 0.1 --f-max 1e5 --ppd 15
```

### Doporučené hodnoty λ

| Typ systému | Doporučené λ | GCV použitelnost |
|-------------|--------------|------------------|
| Čisté RC (Voigt) | 0.01 - 0.1 | ✅ GCV funguje dobře |
| RC + Q | 0.1 - 0.5 | ✅ GCV použitelné |
| **RC + Warburg (Randles)** | **0.5 - 2.0** | ⚠️ **Manuální λ doporučeno** |
| Complex + Warburg | 1.0 - 3.0 | ⚠️ Manuální λ doporučeno |
| Pouze Warburg | > 1.0 | ❌ GCV nepoužívat |

### Proč je to tak?

**Toto není bug**, ale fundamentální omezení DRT metody:

- **DRT předpokládá:** Z(ω) = suma diskrétních RC relaxací
- **Warburg je:** spojitá distribuce (difúze), ne diskrétní relaxace
- **Výsledek:** DRT se snaží aproximovat Warburg jako nekonečně širokou distribuci

→ Pro čisté difúzní systémy použij DFRT (Distribution of Diffusion Times), ne DRT.

### Další tipy

1. **R+L fit pro R_inf:** `--rl-fit` dává robustnější odhad
2. **Vyšší ppd:** 12-15 bodů/dekáda zlepší rozlišení RC píků
3. **Nižší šum:** < 2% pro čistší spektra
4. **GMM peak detection:** Robustnější než scipy.find_peaks

## 🎯 Použití v rešerši

### 1. Testování algoritmů

```bash
# Generuj referenční data
python3 generate_warburg_data.py reference.csv \
    --format csv \
    --circuit complex \
    --noise 0 \
    --ppd 20 \
    --seed 100

# Test různých metod detekce píků
python3 eis_v2.py reference.csv --peak-method scipy --no-show -v
python3 eis_v2.py reference.csv --peak-method gmm --no-show -v
```

### 2. Benchmarking šumu

```bash
# Test robustnosti vůči šumu
for noise in 0.1 0.5 1.0 2.0 5.0 10.0; do
    python3 generate_warburg_data.py "benchmark_${noise}.csv" \
        --format csv \
        --noise $noise \
        --seed 42 \
        --ppd 15

    python3 eis_v2.py "benchmark_${noise}.csv" \
        --auto-lambda \
        --peak-method gmm \
        --save "results_${noise}" \
        --no-show
done
```

### 3. Demonstrace pro výuku

```bash
# Interaktivní demo s vizualizací
python3 generate_warburg_data.py demo.csv \
    --format both \
    --circuit randles \
    --noise 1.0 \
    --ppd 12 \
    --plot

# Následná analýza
python3 eis_v2.py demo.csv --auto-lambda --auto-circuit -v
```

## 🔧 Řešení problémů

### Problém: "ModuleNotFoundError: No module named 'matplotlib'"

**Řešení:**
```bash
sudo apt install python3-matplotlib
# nebo
pip install matplotlib --break-system-packages
```

### Problém: Data nelze načíst v eis_v2.py

**Řešení:** Použij CSV formát místo DTA:
```bash
# Generuj jako CSV
python3 generate_warburg_data.py data.csv --format csv

# Načti
python3 eis_v2.py data.csv --auto-lambda -v
```

### Problém: Warburg charakteristika není vidět

**Možné příčiny:**
1. Příliš úzký frekvenční rozsah → rozšiř `--f-min` a `--f-max`
2. Příliš vysoký šum → sniž `--noise`
3. Málo bodů → zvyš `--ppd`

**Řešení:**
```bash
python3 generate_warburg_data.py data.csv \
    --format csv \
    --f-min 0.001 \
    --f-max 1e6 \
    --ppd 15 \
    --noise 0.5 \
    --plot
```

### Problém: Píky v DRT nejsou správně detekovány

**Řešení:** Použij GMM metodu:
```bash
python3 eis_v2.py data.csv --peak-method gmm --auto-lambda -v
```

## 📚 Teorie

### Warburg impedance

Semi-infinite Warburg popisuje difúzi v neomezené vrstvě:

```
Z_W(ω) = σ/√ω - jσ/√ω
```

kde **σ** (Warburg koeficient) závisí na:
- Difúzním koeficientu **D** [cm²/s]
- Koncentraci **c** [mol/cm³]
- Počtu elektronů **n**

```
σ = RT / (n²F²A√2) · (1/(c·√D))
```

### Randles obvod

Nejjednodušší model elektrochemické buňky:

```
      R_s
       |
    ---|---
    |     |
   R_ct  C_dl
    |     |
    ---|---
       |
      Z_W
```

**Fyzikální význam:**
- **R_s**: Ohmický odpor roztoku
- **R_ct**: Aktivační odpor přenosu náboje
- **C_dl**: Kapacita dvojvrstvy
- **Z_W**: Difúzní impedance

### Charakteristické frekvence

| Proces | Časová konstanta | f_char | Pozice v spektru |
|--------|------------------|--------|------------------|
| Dvojvrstva | τ_dl = R_ct·C_dl | 1/(2πτ_dl) | VF polokruh |
| Difúze | τ_d = l²/D | 1/(2πτ_d) | NF přímka 45° |

## 📖 Reference

- Orazem, M.E., Tribollet, B.: *Electrochemical Impedance Spectroscopy* (2008), Chapter 8
- Bard, A.J., Faulkner, L.R.: *Electrochemical Methods* (2001), Chapter 10
- Warburg, E.: *Ann. Physik Chem.* 67 (1899) 493

## 🔗 Související nástroje

- **eis_v2.py** - Analýza EIS dat (DRT, KK, fitting)
- **eis_analysis/io/synthetic.py** - Generování Voigtových dat (bez Warburgu)

---

**Autor:** EIS Analysis Toolkit v2.0
**Datum:** 2025-12-13
