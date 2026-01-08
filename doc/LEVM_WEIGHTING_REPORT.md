# Systém vážení experimentálních dat v LEVM

## Úvod

Program LEVM (Levenberg-Marquardt) používá vážený komplexní nelineární nejmenší čtverce (CNLS) pro fitování impedančních dat. Systém vážení umožňuje přizpůsobit důležitost jednotlivých datových bodů podle jejich očekávané přesnosti.

Reziduální funkce má tvar:

```
χ² = Σ [ (Y_měřené - Y_model)² / σ² ]
```

kde σ je standardní odchylka (váha) pro každý bod.

## Parametr IRCH - výběr schématu vážení

Schéma vážení se volí parametrem **IRCH** na řádku 3 vstupního souboru. Implementace je v subrutině `RWTS` v souboru `LV2.FOR` (řádky 384-570).

### Přehled schémat

| IRCH | Název | Typ dat | Popis |
|------|-------|---------|-------|
| 0 | Input Weights | C, R, I | Váhy načteny ze vstupního souboru |
| 1 | Unity | C, R, I | Jednotková váha pro všechny body |
| 2 | Power Law | C, R, I | Váhy podle velikosti jednotlivých složek |
| 3 | Modulus | pouze C | Váhy podle modulu komplexního čísla |
| 4 | Spinolo (1174) | pouze C | Pro Solatron 1174 FRA |
| 5 | Orazem-Agarwal | pouze C | Pro Solatron 1250/1286 FRA |
| 6 | Orazem Alt | pouze C | Pro Solatron 1250 & PAR 273 |

**Poznámka:** C = komplexní data, R = pouze reálná část, I = pouze imaginární část

---

## Detailní popis jednotlivých schémat

### IRCH = 0: Input Weights (vstupní váhy)

Váhy (standardní odchylky) jsou načteny přímo ze vstupního datového souboru.

- **Použití:** Když máte experimentálně určené nejistoty měření
- **Požadavky:** Parametry U (P31) a XI (P32) musí být fixovány
- **Výstupní zpráva:** `"WEIGHTS (STANDARD DEVIATIONS) ARE READ IN"`

### IRCH = 1: Unity Weighting (jednotkové vážení)

Všechny datové body mají stejnou váhu rovnou 1.

```
σ = 1 pro všechny body
```

- **Použití:** Exploratorní fit, srovnatelná přesnost všech měření
- **Výhoda:** Jednoduchost, žádné předpoklady o chybách
- **Nevýhoda:** Velké hodnoty dominují fit, malé hodnoty mohou být špatně fitovány
- **Výstupní zpráva:** `"UNIT WEIGHTING ASSIGNED TO EACH POINT"`

### IRCH = 2: Power Law Weighting (mocninné vážení)

Reálná a imaginární složka jsou váženy **odděleně** podle vlastní velikosti.

```
σ²(Re) = |Re|^(2·XI)
σ²(Im) = |Im|^(2·XI)
```

kde XI = P(32) je exponent.

- **Použití:** Obecné proporcionální vážení
- **Parametr:** P(32) = XI (typicky 0 až 2)
- **Výstupní zpráva:** `"WEIGHTS INVOLVE THE MAGNITUDES OF DATA OR FUNCTION VALUES RAISED TO THE POWER [XI]"`

### IRCH = 3: Modulus Weighting (vážení modulem)

Obě složky komplexní impedance jsou váženy **stejně** podle společného modulu.

```
|Z|² = Re² + Im²
σ²(Re) = σ²(Im) = |Z|^(2·XI)
```

- **Použití:** Komplexní impedanční data s konstantní relativní chybou
- **Výhoda:** Konzistentní zacházení s oběma složkami
- **Omezení:** Pouze pro komplexní data (DATTYP = 'C')
- **Výstupní zpráva:** `"MODULUS WEIGHTING: RESULTS RAISED TO THE POWER [XI]"`

### IRCH = 4: Spinolo Weighting (Solatron 1174 FRA)

Specializované vážení pro frekvenční analyzátor Solatron 1174.

```
σ²(Re) = [AK1S·(2·Re² + Im²) + AK2S]^XI
σ²(Im) = [AK1S·(2·Im² + Re²) + AK2S]^XI
```

kde AK1S = 1.0×10⁻⁶, AK2S = 1.0×10⁻⁸

- **Omezení:** Pouze komplexní data
- **Reference:** Spinolo et al.

### IRCH = 5, 6: Orazem-Agarwal Weighting

Pokročilé vážení pro Solatron 1250/1286 a PAR 273 s frekvenčně závislými odpory.

```
σ = OAA·|Im| + OAB·|Re| + OAG·|Z|²·RFI
```

kde:
- OAA = 8.12×10⁻⁴
- OAB = 9.33×10⁻⁴
- OAG = 2.31×10⁻⁴
- RFI = frekvenčně závislý odpor (4 oblasti)

**Frekvenční oblasti:**
| Oblast | Frekvence | Odpor |
|--------|-----------|-------|
| 1 | i ≤ P(34) | 1/P(33) |
| 2 | P(34) < i ≤ P(36) | 1/P(35) |
| 3 | P(36) < i ≤ P(38) | 1/P(37) |
| 4 | i > P(38) | 1/P(39) |

---

## Klíčové parametry vážení

### P(31) = U (aditivní člen)

Přidává konstantní člen k váze pro prevenci nulových vah:

```
σ² = U² + σ²_vypočtená
```

- **Účel:** Noise floor, stabilizace pro malé hodnoty
- **Typická hodnota:** 0 nebo malé kladné číslo

### P(32) = XI (exponent)

Určuje mocninu ve váhovacím vztahu:

```
σ² ∝ |Y|^(2·XI)
```

**Povolený rozsah:** -4 ≤ XI ≤ 4 (kontrolováno v kódu)

**Typické hodnoty:**

| XI | Váha w | Předpoklad o chybě | Označení |
|----|--------|-------------------|----------|
| 0 | w = 1 | Konstantní absolutní | Unity |
| 0.5 | w ∝ 1/\|Z\| | σ ∝ √\|Z\| | Modulus |
| **1.0** | w ∝ 1/\|Z\|² | σ ∝ \|Z\| (konst. relativní) | **Proportional** |
| 2 | w ∝ 1/\|Z\|⁴ | Silný důraz na malé hodnoty | - |

**Doporučení:** XI = 1 je nejběžnější volba pro impedanční spektroskopii (předpoklad konstantní relativní chyby měření).

---

## Power Law vs Modulus - srovnání

| Vlastnost | Power Law (IRCH=2) | Modulus (IRCH=3) |
|-----------|-------------------|------------------|
| Váha Re | Podle \|Re\| | Podle \|Z\| |
| Váha Im | Podle \|Im\| | Podle \|Z\| |
| Nezávislost složek | Ano | Ne (stejná váha) |
| Typ dat | C, R, I | Pouze C |

**Praktický příklad:**

Pro bod s Re = 1000 Ω, Im = 10 Ω, XI = 1:

- **Power Law:** σ(Re) ∝ 1000, σ(Im) ∝ 10
- **Modulus:** σ(Re) = σ(Im) ∝ √(1000² + 10²) ≈ 1000

V tomto případě Modulus dává imaginární složce menší relativní váhu než Power Law.

---

## Function vs Data Weighting

### Záporné IRCH

Pokud je IRCH záporné, program přepne na **function weighting**:

```
IRCH < 0 → IWT = 1, IRCH = |IRCH|
```

- **Data weighting (IWT=0):** Váhy počítány z měřených hodnot
- **Function weighting (IWT=1):** Váhy počítány z hodnot modelu

**Použití function weighting:**
- Iterativní zpřesňování vah během fitu
- Robustnější pro data s outliers
- Vyžaduje dobrý počáteční odhad parametrů

---

## Iterační fáze (JIT parametr)

Fitování probíhá ve fázích řízených parametrem JIT:

| JIT | Fáze | Popis |
|-----|------|-------|
| 1 | Hlavní fit | Unity nebo data-based weighting |
| 2 | Function weighting | Aktivní pokud IFP ≠ 0 |
| 3 | Residual weighting | Aktivní pokud IRE > 0 |
| 4 | Optimalizace | Pouze komplexní data, pokud IOP > 0 |

---

## Reference na zdrojové soubory

| Soubor | Obsah | Klíčové řádky |
|--------|-------|---------------|
| `LV2.FOR` | Subrutina RWTS | 384-570 |
| `LV2.FOR` | Výpočet reziduí v FCN | 114 |
| `LV0.FOR` | Vstup parametrů, výpis schématu | 157-390 |
| `LV1.FOR` | Iterační řízení (JIT) | 34-99 |

---

## Doporučení pro praxi

1. **Začněte s IRCH=1** (unity) pro počáteční exploratorní fit
2. **Přejděte na IRCH=3, XI=1** (proportional modulus) pro finální fit komplexních dat
3. **Použijte IRCH=2** pokud potřebujete nezávislé vážení Re a Im složek
4. **Pro přístrojově specifické vážení** použijte IRCH=4,5,6 s odpovídajícím FRA
5. **U (P31) > 0** může pomoci stabilizovat fit pro data s velmi malými hodnotami

---

*Report vygenerován na základě analýzy zdrojového kódu LEVM verze 8.06*
