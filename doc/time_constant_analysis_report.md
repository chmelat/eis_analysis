# Analýza časových konstant: DRT vs Circuit Fitting

**Datum:** 2026-01-05
**Verze softwaru:** EIS Analysis v0.9.3
**Autor analýzy:** Claude Code

---

## 1. Cíl analýzy

Porovnat přesnost identifikace časových konstant relaxačních procesů pomocí různých metod:
1. **DRT spektrum + GMM detekce píků** (s různou regularizací λ)
2. **Přímý circuit fitting** impedančních dat

Zjistit vliv parametru regularizace λ na výsledky a porozumět rozdílům mezi metodami.

---

## 2. Syntetický model

**Obvod:** Rs - (R0||Q0) - (R1||Q1)

### Parametry modelu:
```
Rs = 10 Ω
R0 = 1.00e+05 Ω,  Q0 = 1.00e-06 S·s^n,  n0 = 0.60
R1 = 8.00e+05 Ω,  Q1 = 3.00e-05 S·s^n,  n1 = 0.43
```

### Teoretické časové konstanty:

Pro CPE element platí: **τ = (R·Q)^(1/n)**

**Proces R0||Q0 (rychlý):**
- τ₀ = (1.00e+05 · 1.00e-06)^(1/0.6) = (0.1)^(1.667)
- **τ₀ = 0.0215 s = 21.5 ms**
- **f₀ = 7.4 Hz**

**Proces R1||Q1 (pomalý):**
- τ₁ = (8.00e+05 · 3.00e-05)^(1/0.43) = (24)^(2.326)
- **τ₁ = 1624 s (27.1 minut)**
- **f₁ = 0.000098 Hz**

### Frekvenční rozsah měření:
- **f_min ≈ 0.01 Hz** (dolní mez)
- **f_max ≈ 100 kHz** (horní mez)
- **Rozsah: ~7 dekád**

**KRITICKÉ POZOROVÁNÍ:** Proces R1||Q1 (f₁ = 0.0001 Hz) je **mimo měřený rozsah** - vyžaduje f_min ≤ 0.0001 Hz.

---

## 3. Metodika

### 3.1 DRT analýza s GMM detekcí píků

**Testované konfigurace:**
1. **λ = 2.8e-05** (automatický výběr, hybridní GCV + L-curve)
2. **λ = 0.2** (silná regularizace, manuální zadání)

**Parametry GMM:**
- BIC threshold: 10.0 (výchozí)
- Rozsah komponent: 1-6
- Klasifikace typů: --classify-terms

### 3.2 Circuit fitting

**Metoda:**
- Diferenciální evoluce (randtobest1bin) + least squares refinement
- Struktura obvodu: R() - (R()|Q()) - (R()|Q())
- Fitování přímo impedančních dat Z(ω)

**DŮLEŽITÉ:** Circuit fitting je **nezávislý na parametru λ** - používá pouze měřená impedanční data.

---

## 4. Výsledky

### 4.1 Srovnávací tabulka časových konstant

| Metoda / Proces | **τ₀ (rychlý)** | **τ₁ (pomalý)** | Komentář |
|-----------------|-----------------|-----------------|----------|
| **TEORETICKÝ MODEL** | **21.5 ms** | **1624 s** | Referenční hodnoty |
| | (f = 7.4 Hz) | (f = 0.0001 Hz) | τ₁ mimo měřený rozsah! |
| | | | |
| **DRT + GMM (λ = 2.8e-05)** | | | |
| → Pík 1 | 7.2 ms | - | 3× kratší než teorie |
| → Pík 2 | - | 15.9 s | **Artifact** (R ≈ 0 Ω) |
| Relativní chyba | -67% | -99% | Pík 2 není reálný proces |
| | | | |
| **DRT + GMM (λ = 0.2)** | | | |
| → Pík 1 | 4.5 ms | - | 4.8× kratší než teorie |
| → Pík 2 | - | 12.7 s (R=72.8kΩ) | Částečný záchyt procesu |
| Relativní chyba | -79% | -99% | Hladší spektrum |
| | | | |
| **CIRCUIT FITTING (λ auto)** | | | |
| Fitted τ | 21.5 ms | 665 s | Náhodný běh 1 |
| Relativní chyba | **0%** ✅ | -59% | Nestabilní pro τ₁ |
| | | | |
| **CIRCUIT FITTING (λ = 0.2)** | | | |
| Fitted τ | 21.5 ms | 1794 s | Náhodný běh 2 |
| Relativní chyba | **0%** ✅ | **+10%** ✅ | Výborná shoda |

### 4.2 GMM detekce - detaily

#### λ = 2.8e-05 (automatická regularizace):
```
BIC optimum: n = 2 komponenty (BIC = 66.15)

Pík 1:
  τ = 7.21e-03 s,  f = 22.1 Hz
  Šířka (σ): 1.511 dekád (VELMI ŠIROKÝ)
  Hranice: [6.85e-06, 7.59] s
  R ≈ 140,474 Ω
  Klasifikace: CPE (high confidence)

Pík 2:
  τ = 15.9 s,  f = 0.01 Hz
  Šířka (σ): 0.001 dekád (EXTRÉMNĚ ÚZKÝ)
  Hranice: [15.8, 16.0] s
  R ≈ 0 Ω  ← BOUNDARY ARTIFACT!
  Klasifikace: VOIGT (high confidence) - chybná!
```

**Interpretace:**
- Pík 1: Kompozitní pík pokrývající měřitelnou část systému
- Pík 2: Artefakt na hranici měření (R ≈ 0 dokazuje, že není reálný)

#### λ = 0.2 (silná regularizace):
```
BIC optimum: n = 2 komponenty (BIC = 463.76)

Pík 1:
  τ = 4.46e-03 s,  f = 35.7 Hz
  Šířka (σ): 1.743 dekád (VELMI ŠIROKÝ)
  Hranice: [1.46e-06, 13.6] s
  R ≈ 211,586 Ω
  Klasifikace: CPE (high confidence)

Pík 2:
  τ = 12.7 s,  f = 0.0125 Hz
  Šířka (σ): 0.095 dekád (ÚZKÝ)
  Hranice: [8.21, 19.7] s
  R ≈ 72,760 Ω  ← MÁ FYZIKÁLNÍ VÝZNAM!
  Klasifikace: VOIGT (high confidence)
```

**Interpretace:**
- Pík 1: Stále kompozitní, ale hladší spektrum
- Pík 2: Již není čistý artefakt - zachycuje část pomalého procesu

### 4.3 Circuit fitting - detaily

#### Běh s automatickou regularizací:
```
Fitted parametry:
  R1 = 6.03e+05 Ω,  Q0 = 3.05e-05,  n0 = 0.448
  R2 = 1.01e+05 Ω,  Q1 = 1.00e-06,  n1 = 0.597

Časové konstanty:
  τ_fit1 = (6.03e5 · 3.05e-05)^(1/0.448) = 665 s
  τ_fit2 = (1.01e5 · 1.00e-06)^(1/0.597) = 21.5 ms

Fit error: 1.67% (excellent)
```

#### Běh se silnou regularizací (λ = 0.2):
```
Fitted parametry:
  R1 = 8.48e+05 Ω,  Q0 = 3.11e-05,  n0 = 0.437
  R2 = 1.01e+05 Ω,  Q1 = 9.95e-07,  n1 = 0.599

Časové konstanty:
  τ_fit1 = (8.48e5 · 3.11e-05)^(1/0.437) = 1794 s
  τ_fit2 = (1.01e5 · 9.95e-07)^(1/0.599) = 21.5 ms

Fit error: 1.31% (excellent)
```

**Pozorování:**
- τ_fit2 je **identický v obou bězích** (21.5 ms = teoretická hodnota)
- τ_fit1 se **liší**: 665 s vs. 1794 s (druhý běh blíže k teorii)
- Rozdíl způsoben **stochastickou náhodou** diferenciální evoluce, **NE** parametrem λ!

---

## 5. Klíčová zjištění

### 5.1 Vliv regularizace λ na DRT spektrum

**Silnější regularizace (λ = 0.2) vs. automatická (λ = 2.8e-05):**

✅ **Pozitiva:**
- Hladší DRT spektrum (méně artificiálních oscilací)
- Pík 2 získal fyzikální význam (R = 72.8 kΩ vs. R ≈ 0 Ω)
- Sloučení malých artificiálních píků
- Lepší reprezentace pomalého procesu (i když stále mimo rozsah)

❌ **Negativa:**
- Větší rozmazání píků (σ = 1.743 vs. 1.511 dekád)
- Posunutí τ dál od teoretické hodnoty (4.5 ms vs. 7.2 ms)
- Vyšší BIC hodnota (463.76 vs. 66.15) - horší statistické skóre

**Závěr:** Pro data s procesy mimo frekvenční rozsah je silnější regularizace výhodnější - vytváří realističtější spektrum bez artificiálních píků.

### 5.2 DRT + GMM vs. Circuit fitting

| Aspekt | DRT + GMM | Circuit fitting |
|--------|-----------|-----------------|
| **Závislost na λ** | Ano - silně ovlivněno | Ne - nezávislé |
| **Přesnost τ₀** | Slabá (3-5× chyba) | Perfektní (0% chyba) |
| **Přesnost τ₁** | Velmi špatná (>99% chyba) | Nestabilní (10-60% chyba) |
| **Stabilita** | Deterministická (pro dané λ) | Stochastická (DE náhoda) |
| **Extrapolace** | Nedokáže (mimo rozsah = artifact) | Částečně možná (používá model) |
| **Vizualizace** | Ano - γ(τ) spektrum | Ne - jen parametry |

**Doporučení:**
- **DRT + GMM:** Vizualizace procesů, identifikace počtu relaxací
- **Circuit fitting:** Kvantitativní odhad parametrů R, Q, n

### 5.3 Proces uvnitř vs. mimo frekvenční rozsah

**R0||Q0 (τ₀ = 21.5 ms, f = 7.4 Hz) - UVNITŘ rozsahu:**
- ✅ Circuit fitting: 100% přesnost (21.5 ms v obou bězích)
- ⚠️ DRT + GMM: 3-5× chyba, ale detekuje proces

**R1||Q1 (τ₁ = 1624 s, f = 0.0001 Hz) - MIMO rozsah:**
- ⚠️ Circuit fitting: 10-60% chyba, nestabilní (závisí na náhodě)
- ❌ DRT + GMM: Detekuje artefakt nebo částečný záchyt (>99% chyba)

**Fundamentální omezení:** Procesy mimo měřený frekvenční rozsah **nelze spolehlivě identifikovat** z dat.

### 5.4 Proč DRT nedává přesné časové konstanty

1. **Regularizace rozmazává píky:**
   - Nutná pro stabilitu inverzního problému
   - Větší λ → širší píky → posunuté τ

2. **GMM fituje distribuci, ne fyzikální proces:**
   - Hledá střed hmoty Gaussovské komponenty
   - CPE vytváří asymetrické distribuce → GMM centrum ≠ skutečné τ

3. **Kompozitní píky:**
   - Jeden GMM pík může reprezentovat více fyzikálních procesů
   - Pík 1 pokrývá 6 dekád → obsahuje více než jeden proces

4. **Boundary efekty:**
   - Na hranicích frekvenčního rozsahu vznikají artefakty
   - R ≈ 0 je spolehlivý indikátor artefaktu

---

## 6. Závěry a doporučení

### 6.1 Hlavní závěry

1. **λ parametr ovlivňuje POUZE DRT spektrum, NIKOLI circuit fitting**
   - DRT + GMM: silně závislé na λ
   - Circuit fitting: zcela nezávislé (fituje přímo Z(ω) data)

2. **Circuit fitting je přesnější pro kvantitativní odhady τ**
   - Pro procesy uvnitř rozsahu: perfektní přesnost
   - Pro procesy mimo rozsah: nestabilní, ale lepší než DRT

3. **DRT + GMM je lepší pro vizualizaci a počet procesů**
   - Poskytuje spektrum γ(τ) - přehledná vizualizace
   - Identifikuje počet relaxačních procesů
   - Ale časové konstanty jsou méně přesné

4. **Silnější regularizace (λ ~ 0.1-0.5) je vhodná pro:**
   - Data s procesy mimo frekvenční rozsah
   - Eliminaci artificiálních píků
   - Realističtější spektrum (i když méně "ostrý")

### 6.2 Doporučení pro praxi

#### Při analýze reálných EIS dat:

1. **Nejprve spusťte DRT s automatickou λ:**
   ```bash
   eis data.DTA --peak-method gmm --classify-terms
   ```
   - Zjistěte počet relaxačních procesů
   - Identifikujte typ (CPE vs. VOIGT)
   - Zkontrolujte boundary artefakty (R ≈ 0, σ < 0.01)

2. **Pokud vidíte mnoho malých píků nebo boundary artefakty:**
   ```bash
   eis data.DTA --peak-method gmm --classify-terms --lambda 0.2
   ```
   - Silnější regularizace sloučí artefakty
   - Realističtější spektrum

3. **Pro kvantitativní parametry použijte circuit fitting:**
   ```bash
   eis data.DTA --circuit "R() - (R()|Q()) - (R()|Q())"
   ```
   - Přesnější odhady R, Q, n
   - Ale vyžaduje znalost správné struktury obvodu

4. **Porovnejte metody:**
   ```bash
   eis data.DTA --peak-method scipy  # Najde všechny features
   ```
   - scipy: liberální (všechny píky)
   - GMM: konzervativní (jen významné)
   - Pravda je obvykle mezi nimi

#### Při tvorbě syntetických dat pro testování:

1. **Vyberte časové konstanty v měřitelném rozsahu:**
   - Pro f_min = 0.01 Hz: τ_max ~ 100 s
   - Pro f_max = 100 kHz: τ_min ~ 1 μs

2. **Doporučený testovací model:**
   ```python
   Rs = 10 Ω
   R0 = 1e5 Ω, Q0 = (1e-6, n=0.6)   # τ = 21.5 ms ✓
   R1 = 5e4 Ω, Q1 = (1e-4, n=0.8)   # τ = 0.84 s  ✓
   ```
   Oba procesy v rozsahu 0.01 Hz - 100 kHz.

3. **Pokud chcete testovat extrapolaci:**
   - Použijte aktuální model (τ₁ = 1624 s)
   - Ale očekávejte velké chyby pro proces mimo rozsah
   - Použijte silnější regularizace (λ ~ 0.2)

### 6.3 Budoucí vylepšení

**Potenciální rozšíření software:**

1. **Adaptivní regularizace:**
   - Automaticky zvýšit λ pokud detekováno mnoho píků s R ≈ 0
   - Nebo pokud píky na boundary (první/poslední dekáda)

2. **Upozornění na boundary artefakty:**
   - Warning pokud pík má R < 0.1 · R_pol
   - Warning pokud pík je v první/poslední dekádě τ

3. **Porovnání DRT τ vs. fitted τ:**
   - Automatický výpočet τ z fitted parametrů
   - Srovnávací tabulka (jako v tomto reportu)

4. **Guidance pro volbu λ:**
   - Doporučení na základě σ píků
   - Pokud average σ > 1.5 dekád → možná silnější λ

---

## 7. Přílohy

### 7.1 Numerické výpočty

**Teoretické τ pro CPE elementy:**

```
τ₀ = (R0 · Q0)^(1/n0)
   = (1.00e+05 · 1.00e-06)^(1/0.6)
   = (0.1)^(1.667)
   = 0.0215 s = 21.5 ms

τ₁ = (R1 · Q1)^(1/n1)
   = (8.00e+05 · 3.00e-05)^(1/0.43)
   = (24)^(2.326)
   = 1624 s
```

**Fitted τ z circuit fittingu (λ = 0.2 běh):**

```
τ_fit1 = (R1 · Q0)^(1/n0)
       = (8.48e+05 · 3.11e-05)^(1/0.437)
       = (26.4)^(2.288)
       = 1794 s

τ_fit2 = (R2 · Q1)^(1/n1)
       = (1.01e+05 · 9.95e-07)^(1/0.599)
       = (0.10045)^(1.669)
       = 0.0215 s = 21.5 ms
```

### 7.2 Spuštěné příkazy

```bash
# Automatická regularizace
eis --peak-method gmm --classify-terms > out_gmm.txt

# Silná regularizace
eis --peak-method gmm --classify-terms --lambda 0.2 > out_gmm_lambda02.txt

# Scipy pro srovnání
eis --peak-method scipy > out_scipy.txt
```

### 7.3 Výstupní soubory

- `out_gmm.txt` - GMM analýza s auto λ
- `out_gmm_lambda02.txt` - GMM analýza s λ = 0.2
- `out_gmm_classify.txt` - GMM s klasifikací typů
- `out_scipy.txt` - Scipy peak detection
- Grafy: PDF soubory s DRT spektry, Nyquist ploty, fits

---

**Konec reportu**

Tento report demonstruje důležitost pochopení rozdílu mezi DRT analýzou (vizualizace) a circuit fittingem (kvantifikace), a ukazuje omezení obou metod při práci s procesy mimo měřený frekvenční rozsah.
