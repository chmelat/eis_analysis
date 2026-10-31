# Generování syntetických dat s Warburgovým chováním

Dokumentace pro `example/generate_warburg_data.py` - nástroj pro simulaci EIS dat vykazujících difúzní impedanci.

## 📋 Přehled

Script umožňuje generování realistických EIS dat pro systémy s:
- **Warburgovým chováním** (difúzní impedance)
- **Volitelným šumem** (0-10%)
- **Vlastním rozlišením** (body na dekádu)
- **Dvěma typy obvodů** (jednoduchý Randles, složitý)
- **Dvěma typy Warburgu** (semi-infinite, finite)

Export do **Gamry DTA** a/nebo **CSV** formátu kompatibilního s `eis.py`.

## 🚀 Rychlý start

### Základní použití

```bash
# Jednoduchý Randles obvod, 1% šum, 10 bodů/dekáda
python3 example/generate_warburg_data.py output.csv --format csv

# S vizualizací
python3 example/generate_warburg_data.py output.csv --format csv --plot
```

### Doporučený workflow

```bash
# 1. Generuj data
python3 example/generate_warburg_data.py warburg_data.csv --format csv --noise 1.5 --ppd 12 --verbose

# 2. Analyzuj pomocí eis.py
# λ se vybírá automaticky; ruční --lambda jen jako volba (viz sekce DRT analýza níže)
python3 eis.py warburg_data.csv -v
```

Hodnoty λ v tomto dokumentu jsou v bezrozměrné škále od v0.47 (viz
GCV_IMPLEMENTATION.md, „Škála λ“). Původní λ = 1.0 odpovídá na výchozí mřížce
(`n_tau = 100`) zhruba 3e-5; přepočet 3e-5 na jednotku staré škály byl změřen
na příkladových spektrech a platí jen přibližně.

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
python3 example/generate_warburg_data.py randles_data.csv --format csv --circuit randles
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
python3 example/generate_warburg_data.py complex_data.csv --format csv --circuit complex
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
python3 example/generate_warburg_data.py data.csv --format csv --warburg-type semi-infinite
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
python3 example/generate_warburg_data.py data.csv --format csv --warburg-type finite
```

## 📊 Příklady použití

### 1. Základní simulace pro testování

```bash
# Čistá data bez šumu
python3 example/generate_warburg_data.py clean_data.csv --format csv --noise 0 --plot

# S realistickým šumem
python3 example/generate_warburg_data.py noisy_data.csv --format csv --noise 2.0 --plot
```

### 2. Vysoké rozlišení pro publikace

```bash
# 15 bodů/dekáda, široký frekvenční rozsah
python3 example/generate_warburg_data.py high_res.csv \
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
python3 example/generate_warburg_data.py reproducible.csv \
    --format csv \
    --noise 1.5 \
    --seed 42
```

### 4. Export do obou formátů

```bash
# Vytvoří data.DTA i data.csv
python3 example/generate_warburg_data.py data --format both --circuit complex --ppd 12
```

### 5. Testování různých úrovní šumu

```bash
# Generuj sérii s různým šumem
for noise in 0 0.5 1.0 2.0 5.0; do
    python3 example/generate_warburg_data.py "noise_${noise}.csv" \
        --format csv \
        --noise $noise \
        --seed 123 \
        --ppd 12
done

# Analyzuj všechny
for f in noise_*.csv; do
    python3 eis.py "$f" --save "${f%.csv}" --no-show -v
done
```

### 6. Validace DRT analýzy

```bash
# Generuj známý systém
python3 example/generate_warburg_data.py validation.csv \
    --format csv \
    --circuit randles \
    --noise 1.0 \
    --ppd 15 \
    --verbose \
    --plot

# Analyzuj s GMM
python3 eis.py validation.csv \
    --peak-method gmm \
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

## ⚠️ DRT analýza Warburg dat

### Automatický výběr λ

Bez `--lambda` vybírá `eis.py` λ automaticky (větší z GCV a rohu
L-křivky) a pro Warburg data je to doporučené nastavení. Změřeno
2026-10-31 na spektrech z tohoto generátoru (výchozí obvody, 0.01 Hz až
100 kHz, 12 bodů/dekádu, šum 0 až 3 %, 5 seedů); „RC píky“ = pík do 0.3
dekády od skutečného τ = R·C:

| Obvod | λ (auto) | RC píky auto / ruční 3e-5 | Chyba rekonstrukce auto / ruční |
|-------|----------|---------------------------|---------------------------------|
| Randles, semi-infinite | ~1.3e-9 | 1/1 / 1/1 | 1.1-4.1 % / 2.8-4.9 % |
| Randles, finite | 2e-11 až 7e-6 | 1/1 / 1/1 | 0.05-3.9 % / 1.7-4.1 % |
| Complex, scipy | 1.6e-9 až 2e-8 | 2/2 / 2/2 | 0.6-3.8 % / 2.1-4.4 % |
| Complex, GMM | totéž | 1.6-2/2 / 1.0-1.2/2 | totéž |

- Na semi-infinite Warburgu je automatická λ malá (~1e-9), RC píky se
  přesto najdou. Ruční λ = 3e-5 je nenajde lépe: má vyšší chybu
  rekonstrukce, víc odporu v nejpomalejší dekádě τ (0.62 proti 0.48
  R_pol) a s GMM u dvou oblouků jeden z nich ztratí.
- S GMM dává automatická λ zhruba o jeden pík navíc, který neodpovídá
  žádnému RC (3 až 4.6 píků proti 2 až 3). Posoudit je pomůže
  `--lambda-probe` (pík stabilní / hraniční / artefakt).
- Starší verze tohoto dokumentu doporučovala ruční λ s tím, že GCV volí
  λ, při kterém RC píky mizí; s dnešním hybridním výběrem λ se to na
  těchto datech nepotvrdilo.

### Proč Warburg v DRT vypadá jako ocas

```
Warburg:  Z_W(ω) = σ(1-j)/√ω  →  γ(τ) ∝ 1/√τ  (roste k dlouhým τ)
RC pík:   Z_RC ~ R/(1+jωτ)    →  γ(τ) = ostrý pík
```

DRT předpokládá Z(ω) jako součet RC relaxací; Warburg je spojitá
distribuce (difúze), kterou DRT aproximuje širokou distribucí až k
okraji mřížky τ. To není chyba metody, ale její omezení: odpor ve
Warburgově ocasu není proces. Pro čistě difúzní systémy se hodí DFRT
(Distribution of Diffusion Times), ne DRT.

### Ruční λ (volitelně)

Pro citlivostní analýzu nebo srovnání s literaturou:

```bash
# Rozsah hodnot, porovnat výsledky
for lambda in 1e-7 1e-6 1e-5 3e-5; do
    python3 eis.py warburg_data.DTA --lambda $lambda --save "lambda_$lambda" --no-show
done

# Stabilita píků vůči λ (λ·10^±0.5 a λ·10^±1)
python3 eis.py warburg_data.DTA --lambda-probe -v

# Případně omez frekvenční rozsah při generování
python3 example/generate_warburg_data.py data.DTA --f-min 0.1 --f-max 1e5 --ppd 15
```

### Další tipy

1. **R_inf:** výchozí odhad je fit R-L-(R|Q) horní dekády spektra;
   `--rinf-plot` ho vykreslí.
2. **Vyšší ppd:** 12-15 bodů/dekáda zlepší rozlišení RC píků.
3. **Nižší šum:** < 2 % pro čistší spektra.
4. **Metoda píků:** na těchto datech scipy (výchozí) našel všechny RC
   píky, GMM s automatickou λ ve 2 obloucích 1.6 až 2 ze 2 a s falešným
   píkem navíc; GMM není automaticky lepší volba.

## 🎯 Použití v rešerši

### 1. Testování algoritmů

```bash
# Generuj referenční data
python3 example/generate_warburg_data.py reference.csv \
    --format csv \
    --circuit complex \
    --noise 0 \
    --ppd 20 \
    --seed 100

# Test různých metod detekce píků
python3 eis.py reference.csv --peak-method scipy --no-show -v
python3 eis.py reference.csv --peak-method gmm --no-show -v
```

### 2. Benchmarking šumu

```bash
# Test robustnosti vůči šumu
for noise in 0.1 0.5 1.0 2.0 5.0 10.0; do
    python3 example/generate_warburg_data.py "benchmark_${noise}.csv" \
        --format csv \
        --noise $noise \
        --seed 42 \
        --ppd 15

    python3 eis.py "benchmark_${noise}.csv" \
        --peak-method gmm \
        --save "results_${noise}" \
        --no-show
done
```

### 3. Demonstrace pro výuku

```bash
# Interaktivní demo s vizualizací
python3 example/generate_warburg_data.py demo.csv \
    --format both \
    --circuit randles \
    --noise 1.0 \
    --ppd 12 \
    --plot

# Následná analýza
python3 eis.py demo.csv -v
```

## 🔧 Řešení problémů

### Problém: "ModuleNotFoundError: No module named 'matplotlib'"

**Řešení:**
```bash
sudo apt install python3-matplotlib
# nebo
pip install matplotlib --break-system-packages
```

### Problém: Data nelze načíst v eis.py

**Řešení:** Použij CSV formát místo DTA:
```bash
# Generuj jako CSV
python3 example/generate_warburg_data.py data.csv --format csv

# Načti
python3 eis.py data.csv -v
```

### Problém: Warburg charakteristika není vidět

**Možné příčiny:**
1. Příliš úzký frekvenční rozsah → rozšiř `--f-min` a `--f-max`
2. Příliš vysoký šum → sniž `--noise`
3. Málo bodů → zvyš `--ppd`

**Řešení:**
```bash
python3 example/generate_warburg_data.py data.csv \
    --format csv \
    --f-min 0.001 \
    --f-max 1e6 \
    --ppd 15 \
    --noise 0.5 \
    --plot
```

### Problém: Píky v DRT nejsou správně detekovány

**Řešení:** Nejdřív posuď stabilitu píků vůči λ, pak porovnej obě metody
detekce (scipy je výchozí, GMM není automaticky lepší, viz sekce DRT
analýza Warburg dat). Odpor na pomalém okraji mřížky τ je u Warburg dat
ocas difúze, ne proces.
```bash
python3 eis.py data.csv --lambda-probe -v
python3 eis.py data.csv --peak-method gmm -v
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

- **eis.py** - Analýza EIS dat (DRT, KK, fitting)
- **eis_analysis/io/synthetic.py** - Generování Voigtových dat (bez Warburgu)

---

**Autor:** EIS Analysis Toolkit v2.0
**Datum:** 2025-12-13
