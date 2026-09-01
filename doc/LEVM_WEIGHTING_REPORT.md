# Systém vážení experimentálních dat v LEVM

## Úvod

Program LEVM (Levenberg-Marquardt) používá vážený komplexní nelineární nejmenší čtverce (CNLS) pro fitování impedančních dat. Systém vážení umožňuje přizpůsobit důležitost jednotlivých datových bodů podle jejich očekávané přesnosti.

Reziduální funkce má tvar:

```
χ² = Σ [ (Y_měřené - Y_model)² / σ² ]
```

kde σ je standardní odchylka (váha) pro každý bod. Rezidua se počítají v `FCN`
(`LV2.FOR:114`) jako `FVEC(I) = (Y(I) - FN(I))/R(I)`, kde `R` drží σ.

> **Poznámka k zápisu.** Podprogram `RWTS` počítá interně **σ²** v poli `FT`
> a teprve na konci (`LV2.FOR:566-568`) vrací `FT = sqrt(FT)`, tedy σ. Vzorce
> níže jsou uvedené pro σ² tak, jak je kód počítá.

Ověřeno proti `LV0.FOR`, `LV1.FOR` a `LV2.FOR`.

## Parametr IRCH - výběr schématu vážení

Schéma vážení se volí parametrem **IRCH** na řádku 3 vstupního souboru. Implementace je v subrutině `RWTS` v souboru `LV2.FOR` (řádky 384-568).

### Přehled schémat

| IRCH | Název | Typ dat | Popis |
|------|-------|---------|-------|
| 0 | Input Weights | C, R, I | Váhy načteny ze vstupního souboru |
| 1 | Unity | C, R, I | Jednotková váha — **pouze ve fázi JIT = 1** |
| 2 | Power Law | C, R, I | Váhy podle velikosti jednotlivých složek |
| 3 | Modulus | pouze C | Váhy podle modulu komplexního čísla |
| 4 | Spinolo (1174) | pouze C | Pro Solatron 1174 FRA |
| 5 | Orazem-Agarwal | C; R/I jen při IRCH = −5 | Pro Solatron 1250/1286 FRA |
| 6 | Orazem Alt | C; R/I jen při IRCH = −6 | Pro Solatron 1250 & PAR 273 |

**Poznámka:** C = komplexní data, R = pouze reálná část, I = pouze imaginární část

---

## Detailní popis jednotlivých schémat

### IRCH = 0: Input Weights (vstupní váhy)

Váhy (standardní odchylky) jsou načteny přímo ze vstupního datového souboru.

- **Použití:** Když máte experimentálně určené nejistoty měření
- **Požadavky:** Parametry U (P31) a XI (P32) musí být fixovány — platí pro celé IRCH < 2
- **Výstupní zpráva:** `"WEIGHTS (STANDARD DEVIATIONS) ARE READ IN"`

Při IRCH = 0 se `RWTS` **vůbec nevolá** (`LV2.FOR:86`). Místo toho se v `FCN`
provede `FJ(I) = R(I)`, tedy načtené σ se použijí beze změny. Žádná z úprav
popsaných níže (U, normalizace, tvarování) se neuplatní.

### IRCH = 1: Unity Weighting (jednotkové vážení)

Všechny datové body mají stejnou váhu rovnou 1.

```
σ = 1 pro všechny body
```

- **Použití:** Exploratorní fit, srovnatelná přesnost všech měření
- **Výhoda:** Jednoduchost, žádné předpoklady o chybách
- **Nevýhoda:** Velké hodnoty dominují fit, malé hodnoty mohou být špatně fitovány
- **Výstupní zpráva:** `"UNIT WEIGHTING ASSIGNED TO EACH POINT"`

> **Unity vážení neplatí po celou dobu běhu.** První řádek větve je
> (`LV2.FOR:410-412`):
> ```fortran
> 100 CONTINUE
>     IF(JIT.GE.2) GOTO 200
> ```
> Ve všech fázích JIT ≥ 2 (function weighting, residual weighting) se tedy
> automaticky přepne na **power-law vážení** podle IRCH = 2. Jednotkové vážení
> platí jen v hlavním fitu.

Aditivní člen U (P(31)) se u této větve **neuplatňuje** — kód skáče rovnou na
konec, mimo blok, který U přičítá.

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
- **Omezení:** Pouze pro komplexní data (DATTYP = 'C'); jinak program skončí
  hláškou `DATTYP MUST BE C FOR THIS IRCH VALUE`
- **Výstupní zpráva:** `"MODULUS WEIGHTING: RESULTS RAISED TO THE POWER [XI]"`

### IRCH = 4: Spinolo Weighting (Solatron 1174 FRA)

Specializované vážení pro frekvenční analyzátor Solatron 1174.

```
σ²(Re) = [AK1S·(2·Re² + Im²) + AK2S]^XI
σ²(Im) = [AK1S·(2·Im² + Re²) + AK2S]^XI
```

kde AK1S = 1.0×10⁻⁶, AK2S = 1.0×10⁻⁸

- **Omezení:** Pouze komplexní data (stejná kontrola jako u IRCH = 3)
- **Reference:** Spinolo et al.

### IRCH = 5, 6: Orazem-Agarwal Weighting

Pokročilé vážení pro Solatron 1250/1286 a PAR 273 s frekvenčně závislými odpory.

```
σ = OAA·|Im| + OAB·|Re| + OAG·|Z|²·RFI
σ² = σ · σ
```

**Konstanty se pro IRCH = 5 a IRCH = 6 liší** (`LV2.FOR:394-395`):

| Konstanta | IRCH = 5 (1250 / 1286) | IRCH = 6 (1250 & PAR 273) |
|-----------|------------------------|---------------------------|
| OAA | 8.12×10⁻⁴ | 1.0×10⁻²⁵ |
| OAB | 9.33×10⁻⁴ | 6.9966×10⁻³ |
| OAG | 2.31×10⁻⁴ | 5.3903×10⁻⁶ |

RFI je převrácená hodnota rozsahového odporu, který se mění po **úsecích
indexu datového bodu** (ne podle frekvence):

| Oblast | Podmínka na index bodu i | RFI |
|--------|--------------------------|-----|
| 1 | i ≤ I1 | 1/P(33) |
| 2 | I1 < i ≤ I2 | 1/P(35) |
| 3 | I2 < i ≤ I3 | 1/P(37) |
| 4 | i > I3 | 1/P(39) |

kde `I1 = int(|P(34)|)`, `I2 = int(P(36))`, `I3 = int(P(38))`.
**Nulová hodnota znamená „až do konce dat" (MD)**, tedy vypnutí dalších oblastí.
Příslušné P(33), P(35), P(37), P(39) musí být kladné, jinak program skončí
hláškou `P(nn) MUST BE POSITIVE FOR THIS IRCH VALUE`.

**Omezení na typ dat:** kontrola je `IF(DATTYP.NE.'C'.AND.IWT.EQ.0)`. Pro reálná
či imaginární data je tedy nutné **function weighting**, tj. IRCH = −5 nebo −6.
Chybová hláška: `IRCH MUST BE -5 FOR DATTYP UNEQUAL TO C`.

**Neimplementováno:** volby PFIT (IPF ≠ 0) pro vážení v magnitudě a fázi.
Program v takovém případě skončí hláškou
`THIS INPUT OPTION IS NOT YET IMPLEMENTED`.

Aditivní člen U (P(31)) a normalizace se u IRCH = 5, 6 **neuplatňují**.

---

## Klíčové parametry vážení

### P(31) = U (aditivní člen)

Přidává konstantní člen k váze pro prevenci nulových vah:

```
σ² = U² + σ²_vypočtená
```

- **Účel:** Noise floor, stabilizace pro malé hodnoty
- **Typická hodnota:** 0 nebo malé kladné číslo
- **Platí pouze pro IRCH = 2, 3, 4.** Blok, který U přičítá (`LV2.FOR:527`),
  je dosažitelný jen z těchto tří větví. U IRCH = 0, 1, 5, 6 nemá P(31) žádný efekt.

### P(32) = XI (exponent)

Určuje mocninu ve váhovacím vztahu:

```
σ² ∝ |Y|^(2·XI)
```

**Povolený rozsah:** −4 ≤ XI ≤ 4. Kontroluje se v `FCN` (`LV2.FOR:43-51`),
ale **pouze když je XI volný parametr** (IXI = 1). Překročení nastaví
`IFLAG = -16` resp. `-15` a fit se ukončí.

**Typické hodnoty:**

| XI | Váha w | Předpoklad o chybě |
|----|--------|-------------------|
| 0 | w = 1 | Konstantní absolutní chyba |
| 0.5 | w ∝ 1/\|Z\| | σ ∝ √\|Z\| |
| **1.0** | w ∝ 1/\|Z\|² | σ ∝ \|Z\| (konstantní relativní chyba) |
| 2 | w ∝ 1/\|Z\|⁴ | Silný důraz na malé hodnoty |

**Doporučení:** XI = 1 je nejběžnější volba pro impedanční spektroskopii (předpoklad konstantní relativní chyby měření).

XI se stane volným parametrem nastavením `NFREE(32) > 0` — to je povoleno jen
pro IRCH > 1 (`LV0.FOR:257`).

### Normalizace vah — IXW

Když je XI volné **nebo** je zapnuté function weighting, nastaví se `IXW = 1`
(`LV0.FOR:278-279`) a váhy se navíc **normalizují geometrickým průměrem**
(`LV2.FOR:528-543`):

```
G = ( Π_i σ²_i )^(1/N)        (geometrický průměr přes všechny body)
σ²_i ← σ²_i / G
```

Bez toho by fit s volným XI mohl minimalizovat χ² prostým zvětšováním vah.
Normalizace je proto nutná, ale znamená, že **absolutní hodnota vah není při
IXW = 1 zachována** — porovnatelná jsou jen relativní čísla.

Platí opět jen pro větve IRCH = 2, 3, 4.

### Tvarování vah na okrajích — ISW

Nastavením **IGACC < 0** na řádku 3 se zapne tvarování vah (`LV0.FOR:614-645`).
Program se interaktivně zeptá na čtyři frekvence `FL1, FL2, FH1, FH2` a sestaví
násobící faktor SHW(i):

| Rozsah frekvence | SHW |
|------------------|-----|
| f ≤ FL1 | `(FL1·FL2)^0.25 / √f` |
| FL1 < f ≤ FL2 | `FL2^0.25 / f^0.25` |
| FL2 < f ≤ FH1 | `1` (nezměněno) |
| FH1 < f ≤ FH2 | `f^0.25 / FH1^0.25` |
| f > FH2 | `√f / √(FH1·FH2)` |

Na konci `RWTS` se pak provede `σ ← σ · SHW` (`LV2.FOR:559-563`). Účelem je
potlačit vliv krajních frekvenčních bodů, které bývají nejméně spolehlivé.

### Robustní ořez reziduí — ROE / RKE

Parametr **ROE** se čte z řádku 2 vstupního souboru. Je-li ROE > 0, proběhne
fit dvakrát: po prvním průchodu se nastaví `RKE = ROE · σ_fit` (`LV1.FOR:597`)
a fit se zopakuje s **ořezanými rezidui** (`LV2.FOR:116-127`):

```
|r| ≤ RKE  →  FVEC = r            (nejmenší čtverce)
|r| > RKE  →  FVEC = sign(r)·RKE  (absolutní hodnota)
```

Jde o robustní M-odhad Huberova typu, který omezuje vliv odlehlých bodů.
Ve výpisu se objeví hláška
`OBJECTIVE FUNCTION, S, INVOLVES LS/AB SWITCH VALUE, RKE=`.
Záporná hodnota ROE se ignoruje (`LV0.FOR:463` ji vynuluje).

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

Pokud je IRCH záporné, program přepne na **function weighting** (`LV0.FOR:157-168`):

```
IRCH < 0 → IWT = 1, IRCH = |IRCH|
```

- **Data weighting (IWT=0):** Váhy počítány z měřených hodnot — `CALL RWTS(...,Y,FJ)`
- **Function weighting (IWT=1):** Váhy počítány z hodnot modelu — `CALL RWTS(...,FN,FJ)`

**Použití function weighting:**
- Iterativní zpřesňování vah během fitu
- Robustnější pro data s outliers
- Vyžaduje dobrý počáteční odhad parametrů

Function weighting zapíná i `IXW = 1`, tedy normalizaci popsanou výše.
Ve výpisu se režim značí jako `UNIT` / `DATA` / `FUNC`.

**Optimalizace (JIT = 4) je při function weightingu zakázána** — program vypíše
`OPTIMIZATION DISABLED FOR FUNCTION WEIGHTING!: USE UWT OR PWT!` a fázi přeskočí.

---

## Iterační fáze (JIT parametr)

Fitování probíhá ve fázích řízených parametrem JIT (`SPFIT`, `LV1.FOR:570-700`):

| JIT | Fáze | Podmínka aktivace |
|-----|------|-------------------|
| 1 | Hlavní fit | vždy |
| 2 | Function weighting (FWT) | IFP ≠ 0 |
| 3 | Residual weighting (RWT) | IRE > 0 **a zároveň ROE = 0** |
| 4 | Optimalizace | IOP > 0, pouze DATTYP = 'C', a **ne** při function weightingu |

Fáze 2 a 3 iterují váhy, dokud se σ nestabilizuje: JIT = 2 končí při
relativní změně < 10⁻³, JIT = 3 při `|1 − σ| < 10⁻²`, nejvýše však po `IFP`
resp. `|IRE|` iteracích.

**Fáze 4 (optimalizace)** přeškáluje váhy reálné a imaginární části podle poměru
jejich směrodatných odchylek (`LV1.FOR:628-641`):

```
r = SDWR/SDWI
σ_Im ← σ_Im · √(2/(1 + r²))
σ_Re ← σ_Re · r · √(2/(1 + r²))
```

Opakuje se, dokud `|1 − r| < 10⁻²` nebo dokud počet iterací nedosáhne IOP.
Účelem je vyrovnat příspěvky obou složek k χ².

Pro JIT > 3 se `RWTS` už nevolá (`LV2.FOR:86`).

---

## Staged free/fixed weighting (volné XI)

`MAINCLC` (`LV1.FOR:33-129`) obsahuje samostatný mechanismus, který se aktivuje
při **CELCAP < 0** a nenulovém `NFREE(32)`. Střídavě optimalizuje dvě disjunktní
sady parametrů:

- **Sada W:** pouze váhovací parametry P(31) a P(32)
- **Sada F:** pouze modelové parametry P(1)–P(30), s P(31) a P(32) fixovanými

Cyklus běží nejvýše **75 iterací** a končí, jakmile `|XI_nové − XI_staré| < 10⁻⁵`.
Průběh se vypisuje pod hlavičkou `******** STAGED FREE/FIXED WEIGHTING ********`.

Vedlejší efekty, které je nutné znát:

- Při **IRCH = 1** se uvnitř tohoto režimu nastaví `IRCH = 2` a `IWT = 1`,
  tedy efektivně **IRCH = −2**. Zvolili jste unity vážení, ale poběží power-law
  function weighting.
- Při `CELCAP = -2` nebo `-4` se `NFREE(31)` vynutí na 2.
- Při `CELCAP ≤ -3` se v druhé polovině cyklu vynutí `IRCH = 2` s vypnutým
  IXI, IXW i IWT.

---

## Volné parametry a vážení — NFREE

| NFREE | Význam |
|-------|--------|
| 0 | Fixní |
| **1** | **Volný, omezený na kladné hodnoty** — při záporné hodnotě se aktivuje exponenciální penalizační funkce (`LV2.FOR:881` a `LV2.FOR:1004`) |
| **2** | **Volný bez omezení — smí jít do záporu** |
| 3 | Volný, a po fitu se jeho vliv odečte z dat do OUTIN (`LV0.FOR:931`) |

> **Pozor na obrácený význam 1 a 2.** Penalizace se aplikuje na NFREE = **1**,
> nikoli 2. Pokud fit hlásí singulární matici, nekonverguje, nebo se odhad
> parametru zmenšuje k nule, přepněte jeho NFREE z 1 na **2**, aby směl jít
> do záporu.

Pro **IRCH < 2** musí být U i XI fixní. Jinak program vypíše
`BAD CHOICE OF FREE/FIX FOR U AND XI` a obě NFREE vynutí na 0 (`LV0.FOR:250`).

---

## Reference na zdrojové soubory

| Soubor | Obsah | Klíčové řádky |
|--------|-------|---------------|
| `LV2.FOR` | Subrutina `RWTS` | 384-568 |
| `LV2.FOR` | Rozhodnutí, čím vážit (data vs. model) | 86-96 |
| `LV2.FOR` | Výpočet reziduí v `FCN`, ořez RKE | 108-127 |
| `LV2.FOR` | Kontrola rozsahu XI | 43-51 |
| `LV2.FOR` | Normalizace vah (IXW) | 528-543 |
| `LV2.FOR` | Tvarování vah (ISW) | 559-563 |
| `LV2.FOR` | Penalizace záporných parametrů (NFREE = 1) | 881-886, 1004-1014 |
| `LV0.FOR` | IRCH < 0 → IWT, nastavení IXI/IXW | 157-168, 250-281 |
| `LV0.FOR` | Výpis zvoleného schématu | 352-392 |
| `LV0.FOR` | Příprava tvarování vah (ISW/SHW) | 612-645 |
| `LV1.FOR` | Staged free/fixed weighting | 33-129 |
| `LV1.FOR` | Řízení fází JIT, RKE, optimalizace | 570-700 |

---

## Doporučení pro praxi

1. **Začněte s IRCH=1** (unity) pro počáteční exploratorní fit
2. **Přejděte na IRCH=3, XI=1** (proportional modulus) pro finální fit komplexních dat
3. **Použijte IRCH=2** pokud potřebujete nezávislé vážení Re a Im složek
4. **Pro přístrojově specifické vážení** použijte IRCH=4,5,6 s odpovídajícím FRA
5. **U (P31) > 0** může pomoci stabilizovat fit pro data s velmi malými hodnotami —
   ale jen při IRCH = 2, 3 nebo 4
6. **Nefunguje-li fit**, nejprve zkontrolujte NFREE: parametr, který potřebuje jít
   do záporu, musí mít NFREE = 2, ne 1
7. **Pro data s odlehlými body** zvažte ROE > 0 (robustní ořez) místo ručního
   odstraňování bodů
