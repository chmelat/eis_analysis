# O Circuit (OSUB) — podrobný popis

Rozšíření přehledové kapitoly v `CIRCUITS.md`. Ověřeno proti `LV6.FOR`
(`OSUB`, `GEDAE`, `GDAEFN1`, `CDFG`, `HNFG`), `LV0.FOR` (vstup, transformace,
závěrečný dielektrický výpis) a `LV5.FOR` (`DISTEL`). Jména parametrů odpovídají
manuálu, str. 5-29 až 5-32.

Notace topologie je stejná jako v `CIRCUITS.md`: **`-` sériové spojení**
(pomlčka, ne minus), **`|` paralelní spojení**.

---

## 1. K čemu obvod je

O-obvod není klasický náhradní obvod se soustředěnými prvky. Je to **rámec pro
fitování disperzní odezvy materiálu**, kde disperze může být dvojího původu:

| Zkratka | Typ | Fyzikální původ | Parametry |
|---------|-----|-----------------|-----------|
| **CSD** | conductive-system dispersion | pohyblivé nosiče náboje (ionty, vakance, elektrony, díry) v neuspořádaném materiálu | P(3)–P(10) |
| **DSD** | dielectric-system dispersion | relaxace dipólů, dielektrická odezva | P(13)–P(20) |

Obě části se dají zapnout nezávisle nebo současně, každá s jiným modelem
disperze. K tomu obvod obsahuje samostatnou **elektrodovou část** (P(11),
P(21)–P(30)) pro jevy na rozhraní vzorek–elektroda.

Autor ho v `README.txt` doporučuje jako první volbu pro pevné materiály:

> *It is recommended that one should use the O-circuit for fitting data for
> solid materials when practical.*

Manuál (str. 5-30) dodává, že je zvlášť vhodný pro dielektrická a
vysokoodporová data. V `FITTESTS` má 61 testovacích souborů — víc než všechny
ostatní obvody dohromady kromě H.

**Podmínka:** N (počet parametrů, řádek 3) musí být **≥ 40**.

---

## 2. Topologie

### Blokové schéma

```
                    ┌── R_INF - Z_CSD ──┐
                    │                   │
 Z = [  ────────────┼── Z_DSD ──────────┼──────── ] - Z_elektroda - [série]
                    │                   │
                    ├── G_D = P(12) ────┤
                    │                   │
                    └── C_INF = P(2) ───┘
```

### Přesná rekonstrukce z kódu

**Objemová část** (`LV6.FOR:393`, výpočet přes admitanci):

```
Y_objem = 1/( P(1) - Z_CSD )  +  Y_DSD  +  P(12)  +  jω·CELC·P(2)
Z_objem = 1 / Y_objem
```

topologicky:

```
Z_objem = ( R_INF - Z_CSD ) | Z_DSD | G_D | C_INF
```

kde `R_INF = P(1)`, `G_D = P(12)` (vodivost), `C_INF = P(2)`.

`Z_CSD` počítá blok zvolený parametrem P(35), `Y_DSD` blok zvolený P(40).
DSD blok vrací komplexní kapacitu nebo permitivitu, kterou kód převede na
admitanci (`LV6.FOR:384`):

```
Y_DSD = CELC · ω · ( −Im(ε) + j·Re(ε) )
```

**Elektrodová část** (`LV6.FOR:428-481`):

```
Z_elektroda = ( R2 - ( (R3 - DCE3) | C3 ) ) | C2 | G_E
```

kde `R2 = P(26)`, `C2 = P(27)`, `R3 = P(28)`, `C3 = P(29)`, `G_E = P(11)`
a DCE3 je distribuovaný element `DISTEL(P(21),P(22),P(23),P(24), NDE3=P(25))`.

**Celkem:**

```
Z = Z_objem - Z_elektroda - [sériový prvek P(30) podle P(37)]
```

| P(37) | Význam P(30) |
|-------|--------------|
| > 0 | sériový odpor R<sub>s</sub> |
| −4 (obecně < 0 a > −9) | sériová indukčnost L |
| −16 (obecně < −9 a > −19) | sériová kapacita C<sub>s</sub> |
| jinak | žádný sériový prvek |

Manuál (str. 5-31): *„although P(30) may be a negative quantity in the input
file, it is actually positive in the fitting circuit."*

### Poznámka k P(11) a P(12)

Komentář ve zdroji u P(11) říká `HERE, P(11) IS A PARALLEL RESISTANCE`
(`LV6.FOR:426`), ale kód ji přičítá k **admitanci**:

```fortran
YE = (1.D0/ZR) + IOMEGA*C2 + P11
```

Manuál to potvrzuje jednoznačně (str. 5-31): *„the parallel conductances (or
conductivities for specific data), GE and GD, are treated as conductances, not
resistances, in the fitting."* **Komentář ve zdroji je zavádějící — P(11) i
P(12) jsou vodivosti.**

---

## 3. Volba modelu: NCH1 = P(35), NCH2 = P(40)

Dvě nezávislé volby. P(35) vybírá model pro CSD blok, P(40) pro DSD blok.
Obě používají tutéž číselnou škálu (`LV6.FOR:49-53`):

| NCH | Model | Popis |
|----:|-------|-------|
| 0 | — | blok se přeskočí; P(1) a P(2) zůstávají účinné |
| 1 | **GED** | generalized exponential distribution of relaxation times, tři podtypy |
| 2 | **CD** | Cole-Davidson s cutoffem |
| 3 | **HN** | Havriliak-Negami s cutoffem (bez CD) |
| 4 | **DE** | libovolný distribuovaný element ze `DISTEL` |
| 5, 6 | DE + coupling | jako 4, navíc Ngai coupling model na vysokých frekvencích |

Modely 1–3 se počítají **přímou integrací přes distribuci relaxačních časů
s konečnými mezemi (cutoff)**. To je jejich hlavní přednost: zaručuje fyzikální
realizovatelnost, kterou analytické tvary bez cutoffu nemají, a umožňuje
počítat i přechodovou odezvu. Manuál pro ně používá souhrnnou zkratku **DWC**
(dispersion with cutoff).

### Dvojí význam P(10) a P(20)

**Podle zvoleného NCH mění dva parametry význam** (manuál, str. 5-29):

| Parametr | NCH = 1 (GED) | NCH = 4 (DE) |
|----------|---------------|--------------|
| P(10) | PHIC2 — druhý exponent φ₂ | **NDEC** — číslo elementu NELEM |
| P(20) | PHID2 — druhý exponent φ₂ | **NDED** — číslo elementu NELEM |

Pro CD a HN (NCH = 2, 3) **musí být P(10), resp. P(20), nastaveno na 6**
(`LV6.FOR:56`).

Podobně **P(8) a P(18) (UC2, UD2) existují jen v režimu GED.** Kód je při
NCH > 1 a NELEM ≠ 9 vynuluje (`LV6.FOR:194,205`):

```fortran
IF(NCH1.GT.1.AND.NDE1.NE.9) P(8) = 0
IF(NCH2.GT.1.AND.NDE2.NE.9) P(18) = 0
```

### Automatické vypnutí bloku

Kód blok sám vypne, pokud v něm nejsou data (`LV6.FOR:198,209`):

```fortran
IF(P(5).EQ.0.D0.AND.P(6).EQ.0.D0)  NCH1 = 0     ! CSD
IF(P(15).EQ.0.D0.AND.P(16).EQ.0.D0) NCH2 = 0    ! DSD
```

a zároveň doplní chybějící γ na 1: `IF(P(3).EQ.0) P(3) = 1`, obdobně P(4),
P(13), P(14).

### Kombinace použité v testech

Ze 61 souborů v `FITTESTS/CKT.O`:

| (P35, P40) | Počet | Význam |
|------------|------:|--------|
| (4, 0) | 41 | jen CSD, přes libovolný DCE — **nejběžnější případ** |
| (4, 4) | 8 | CSD i DSD, obojí přes DCE |
| (0, 4) | 5 | jen DSD |
| (3, 0) | 4 | jen CSD, HN s cutoffem |
| (0, 0) | 2 | jen R_INF a C_INF (šablony) |
| (2, 0), (1, 0), (0, 1) | po 1 | CD, GED |

Neplatná hodnota vede k `PAUSE 'WRONG VALUE OF P(35)'`, resp. `P(40)`.

---

## 4. Kompletní tabulka parametrů

Jména podle manuálu str. 5-29/5-30. Kde je dvojice, první platí pro DWC
(NCH 1–3), druhé pro DE (NCH = 4).

### Objemová část

| P | Jméno | Význam |
|---|-------|--------|
| 1 | RINF | vysokofrekvenční limita odporu, sériově před CSD blokem |
| 2 | CINF | vysokofrekvenční kapacita C<sub>g</sub> (nebo ε<sub>∞</sub>, je-li ATEMP < 0) |
| 12 | GD | paralelní **vodivost** přes celý objem (např. svodová) |

### CSD blok (identifikátor „C")

| P | Jméno | Význam |
|---|-------|--------|
| 3 | GAMC1 | γ₁ — tvarový exponent distribuce (jen GED, HN) |
| 4 | GAMC2 | γ₂ — totéž pro druhou větev |
| 5 | RD, RDE | (R₀ − R<sub>INF</sub>), měřítko disperze; často ρ₀ |
| 6 | TAUC, TDEC | τ₀ — charakteristický relaxační čas |
| 7 | UC1, UDEC | U₁ — mez cutoffu, znaménko volí podtyp (viz 5.1) |
| 8 | UC2 | U₂ — druhá mez, **jen GED** |
| 9 | PHIC1, PDEC | φ₁ — exponent |
| 10 | PHIC2 **/ NDEC** | φ₂ (GED) **nebo NELEM** (DE) |
| 35 | NCH1 | volba modelu CSD |

### DSD blok (identifikátor „D")

| P | Jméno | Význam |
|---|-------|--------|
| 13 | GAMD1 | γ₁ |
| 14 | GAMD2 | γ₂ |
| 15 | CD, CDE | (C₀ − C<sub>∞</sub>) — měřítko dielektrické disperze |
| 16 | TAUD, TDED | τ na dielektrické úrovni |
| 17 | UD1, UDED | U₁ |
| 18 | UD2 | U₂, **jen GED** |
| 19 | PHID1, PDED | φ₁ |
| 20 | PHID2 **/ NDED** | φ₂ (GED) **nebo NELEM** (DE) |
| 40 | NCH2 | volba modelu DSD |

### Elektrodová část

| P | Jméno | Význam |
|---|-------|--------|
| 11 | GE | paralelní **vodivost** v elektrodové části (nebo kapacita, je-li P(38) = 64) |
| 21–24 | RDE3, TDE3, UDE3, PDE3 | parametry elektrodového DCE |
| 25 | NDE3 | typ elektrodového DCE |
| 26 | R2 | |
| 27 | C2 | |
| 28 | R3 | |
| 29 | C3 | |
| 30 | R / L / C | sériový prvek, roli volí P(37) |

### Řídicí a transformační

| P | Význam |
|---|--------|
| 33 | > 0: objemový podíl pro Bruggemanovo efektivní medium (viz 8) |
| 36 | ≠ 0: odečtení σ₀ z dat (viz 10.1) |
| 37 | volba role P(30) |
| 38 | škálování frekvenční osy; hodnota 64 přepíná elektrodovou část (viz 7) |
| 39 | škálování vodivostních dat |
| 41–44 | Ngai coupling (NCH = 5, 6), nebo parametry DCE4 při P(38) = 64 |
| 45 | NELEM pro DCE4 při P(38) = 64 |

---

## 5. Modely disperze

### 5.1 GED — generalized exponential (NCH = 1)

Distribuce relaxačních časů se integruje v proměnné **y = ln(τ/τ₀)**.
Integrand je `GDAEFN1` (`LV6.FOR:814-842`):

```
G(y) ∝ exp( −sgn(φ)·|φ·y|^γ / γ )
```

a odezva ve frekvenční doméně:

```
I(ω) = ∫ G(y) / (1 + (ω·τ₀·e^y)²) dy        (reálná část, ICHG = 0)
       ∫ G(y)·e^y / (1 + (ω·τ₀·e^y)²) dy    (imaginární část, ICHG = 1)
```

Výsledek `ZI0C = RN·(Re − j·ω·τ₀·Im)`, škálovaný `RDAE = P(5)`.

**Exponent γ určuje tvar distribuce:**

| γ | Distribuce |
|---|-----------|
| 1 | čistě exponenciální — **EDAE** |
| 2 | **gaussovská** DRT |
| jiné | zobecněná exponenciální |

Manuál uvádí užitečný rozsah γ zhruba **0,5 až 4**.

**Tři podtypy podle znamének U₁ a U₂** (`LV6.FOR:69-76`, meze integrace
`LV6.FOR:743`):

| U₁ | U₂ | Meze integrace | Chování |
|----|----|----------------|---------|
| < 0 | = 0 | −\|U₁\| … 0 | **asymetrická**, jen γ₁; při γ₁ = 1 je to EDAE1 |
| > 0 | = 0 | −\|U₁\| … +\|U₁\| | **symetrická**, jen γ₁; γ₁ = 1 → EDAE2, γ₁ = 2 → gaussovská DRT |
| < 0 | > 0 | −\|U₁\| … 0 a 0 … U₂ | **nejobecnější** — různá U, φ i γ pro obě větve |

φ₁, γ₁, U₁ ovlivňují **vysokofrekvenční** oblast, φ₂, γ₂, U₂ **nízkofrekvenční**.
Manuál dodává praktické vodítko: −φ₁ je log-log sklon vysokofrekvenční části
(zvlášť je-li \|φ₁\| ≤ 0,5), φ₂ sklon nízkofrekvenční části.

**Normalizace** (`LV6.FOR:595-625`) používá Γ(1/γ) a uzavřený tvar
(1 − e^(−φ\|U\|))/φ; pro γ = 1 je exaktní i při konečném cutoffu.

### 5.2 CD — Cole-Davidson s cutoffem (NCH = 2)

Integrand `CDFG` (`LV6.FOR:844-880`), substituce `x = 1/(1 + y^(1/(1−φ)))`,
normalizační konstanta `CDN = sin(φπ)/(π(1−φ))`. Integruje se metodou
`QROMO` (otevřená, protože integrand má singularitu na okraji) v mezích
0 … (e^\|U₁\| − 1)^(1−φ).

Vyžaduje **P(10) = 6**.

### 5.3 HN — Havriliak-Negami s cutoffem (NCH = 3)

Integrand `HNFG` (`LV6.FOR:883-918`), y = ln(τ/τ₀):

```
θ(x) = atan2( sin(γπ), x^γ + cos(γπ) )

G(y) = (1/π) · x^(γφ) · sin(φ·θ) / ( x^(2γ) + 2x^γ·cos(γπ) + 1 )^(φ/2)
```

To je standardní HN distribuce relaxačních časů. Integruje se `QROMB`
v symetrických mezích −\|U₁\| … +\|U₁\|.

Vyžaduje **P(10) = 6**. Zahrnuje HN bez CD; čistý CD je NCH = 2.

### 5.4 DE — libovolný distribuovaný element (NCH = 4)

Nejpoužívanější volba. Zavolá se prostě

```fortran
ZT = DISTEL(P(5),P(6),P(7),P(9), NDE1, FREQ(I))     ! CSD, NDE1 = P(10)
YT = DISTEL(P(15),P(16),P(17),P(19), NDE2, FREQ(I)) ! DSD, NDE2 = P(20)
```

Pozor: **P(8) se nepoužije** (je vynulován) a **φ je P(9), ne P(10)** —
P(10) nese číslo elementu.

Elementy skutečně použité v `FITTESTS/CKT.O`:

| NELEM | Počet | Model |
|------:|------:|-------|
| 10 | 18 | KWW (Williams-Watts) bez cutoffu, rychlý |
| 37 | 7 | KWW DRT pro libovolné β |
| 6 | 7 | Havriliak-Negami model 1 |
| 12 | 6 | EDAE1 |
| 7 | 5 | HN model 2 (v současném kódu totožné s 6) |
| 33 | 4 | Arrheniovský teplotní fit |
| 35 | 4 | KWW s cutoffem, β₁ = 1/3 (model UN) |
| 36 | 3 | KWW s cutoffem, i pro velmi malá β |
| 9 | 3 | zobecněný konečný Warburg |
| 27 | 2 | PCPE (nearly constant loss) |
| 25, 32 | po 1 | efektivní medium s PCPE, KWW s cutoffem |

Manuál (str. 4-10) k volbě mezi 10 a 36: *„See NELEM = 36 for KWW response with
or without cutoff. Use NELEM = 10 when no cutoff and increased convergence speed
but somewhat reduced accuracy is desired."*

### 5.5 NCH = 5, 6 — Ngai coupling model

Jako NCH = 4, ale nad prahovou frekvencí ω ≥ 1/P(41) se odezva nahradí
Debyeovým členem (`LV6.FOR:294-297` (CSD), `LV6.FOR:342-345` (DSD)):

```fortran
IF(FREQ(I).GE.1/P(41)) ZT = P(42)/(1.D0 + IOMEGA*P(43))
```

P(41) je minimální τ, P(42) a P(43) jsou síla a čas Debyeova členu.
Manuál doporučuje **místo toho použít cutoff model** (NCH = 1–3).

---

## 6. MODE / MDE — CSD0 vs CSD1 a přechodová odezva

MODE je 6. pole na řádku 3 vstupního souboru; kód si ho ukládá jako
`MDE = MODE` (`LV0.FOR:118`).

| MODE | Význam |
|------|--------|
| 0, 1 | běžná **CSD0** (nebo DSD) odezva |
| < 0 | **CSD1** odezva (JRM/Moynihan) |
| −1 | CSD1a |
| −3 | CSD1b |
| −6 | CD a HN dávají CSD1a |
| \|MODE\| = 8 | **přechodová (transientní) odezva** místo frekvenční |

**Znaménko MODE nemá vliv na DSD odezvu typu GED.**

Pro \|MODE\| = 8 je nutné `DATTYP = R` (fit jen reálné části). V kódu se to
projeví přepnutím `IWT = 0` (`LV6.FOR:578-583`), což změní integrand
z Lorentzova jádra `1/(1+(ωτ)²)` na exponenciální `exp(−t/τ)`.

Manuál (str. 4-10) k použití s KWW:

> *With the O-circuit, use MODE = MDE = 1 for KD or K0 response, and MODE = −1
> for K1-model response. Use NELEM = 35 for the K1 model with beta1 fixed at
> 1/3: the quasi-universal UN model.*

V testech: MODE = −1 (26×), 1 (17×), −2 (10×), 0 (9×), −6 (1×).

### Dva CSD elementy v sérii

Manuál (str. 5-32) popisuje případ, kdy chcete v jednom modelu CSD1 i CSD0.
Při MODE < 0 jsou ve výchozím stavu **oba** elementy CSD1. Druhý (DCE3) lze
donutit k CSD0 charakteru:

| Element | Jak vynutit CSD0 |
|---------|------------------|
| HN | \|ATEMP\| > 900 |
| KWW, NELEM = 10 | P(23) < −900, fixované |
| libovolný DCE při NCH1 = 4 | ATEMP < −900 |

V kódu je to vidět na `LV6.FOR:436-440`:

```fortran
IF(P(23).LT.-9.0D2.OR.ATEMP.LT.-9.0D2) THEN
    MDE = 1        ! vynutí CSD0 pro elektrodový DCE
ELSE
    MDE = -1
ENDIF
```

Příklady: `OTSCOMP1`, `OTSCOMP2`.

---

## 7. Elektrodová část

Standardní tvar (P(38) ≠ 64):

```
Z_elektroda = ( R2 - ( (R3 - DCE3) | C3 ) ) | C2 | G_E
```

Je-li `NDE3 = 0` **a zároveň** `R3 = 0`, zbude jen C3 (nebo nic, je-li i C3 = 0).

**Nejčastější volba NDE3 v testech:**

| NDE3 | Počet | Element |
|-----:|------:|---------|
| 2 | 9 | CPE |
| **26** | 6 | **SCPE** — power law #7, elektrodové jevy, Z úroveň, sériový |
| 6 | 1 | Havriliak-Negami |
| 10 | 1 | KWW |

SCPE (NELEM = 26) je kanonická volba pro elektrodové jevy — `README.txt` na něj
odkazuje jako na „the SCPE electrode-effects contribution".

### Varianta P(38) = 64

Při `P(38) = 64` a `P(35) ≠ 0` se elektrodová část přepne (`LV6.FOR:409-426`):

```
DCE4 = DISTEL( P(41), P(42), P(43), P(44), NDE4 = P(45) )
Z_EL = DCE4 - C_P11             ! P(11) je zde KAPACITA, ne vodivost
```

a `Z_EL` se zapojí paralelně dovnitř elektrodové části (přes `YEL`), zatímco
`P11` se v hlavním vztahu vynuluje. **P(11) tedy v tomto režimu mění význam
z vodivosti na kapacitu** a P(41)–P(45) přestávají patřit Ngai coupling modelu.

---

## 8. Efektivní medium — P(33)

Je-li **P(33) > 0**, kód nepoužije obvyklé skládání, ale Bruggemanův vztah pro
efektivní medium (`LV6.FOR:352-372`). CSD admitance se převede na permitivitu
(`YZ ← YZ/(jω·CELCAP)`), pak:

```
ε_eff = 3·P(33)·(ε_D − ε_C) / ( ε_D + 2ε_C − P(33)·(ε_D − ε_C) )
ε_m   = P(2) + ε_C·(1 + ε_eff)
Z     = 1 / ( jω·CELCAP·ε_m )
```

P(33) je objemový podíl inkluzí. V `FITTESTS/CKT.O` tuto volbu nepoužívá
ani jeden soubor — je nevyzkoušená.

---

## 9. Závěrečný dielektrický výpis

Nad rámec obvyklých fitovacích výstupů vypíše O-obvod ještě blok
s **momenty distribuce a limitními dielektrickými hodnotami**
(`LV0.FOR:1461-1595`). Je to hlavní přidaná hodnota O-obvodu oproti
prostému fitu.

### Kdy se objeví

Rozhodující je příznak **INDE**, nastavený v `LV0.FOR:797-813`. Blok se vypíše
jen při `INDE = 0`, protože `LV0.FOR:1353` obsahuje `IF(INDE.GE.1) GOTO 669`:

```fortran
IF(FUN.EQ.OO.AND.(NELEM.EQ.7.OR.NELEM.EQ.10.OR.NELEM.EQ.32
*  .OR.NELEM.EQ.36)) THEN
        INDE = 0
  ELSEIF(IRE.LT.0) THEN
        INDE = 1
  ELSE
        INDE = 2
ENDIF
C
IF(NELEM.EQ.7.AND.P(7).NE.1.D0) THEN
        P(10) = 6.D0
        NELEM = 6
        INDE = 1
ENDIF
```

Z toho plynou dvě tvrdé podmínky:

1. **NELEM musí být 7, 10, 32 nebo 36.** Žádná jiná hodnota blok nevypíše —
   ani 35, ani 37, ani 6. NELEM se bere z P(10), je-li P(35) = 4, a z P(20),
   je-li P(40) = 4 (`LV0.FOR:787-795`).
2. **Je-li NELEM = 7, musí být P(7) = 1.** Jinak kód přepíše NELEM na 6
   a blok potlačí. (Doloženo dvojicí `OTSM1YC7` s P(7) = 1, která blok dostane,
   a `OTSM1Y07` s P(7) = 0,7996, která ne — jinak jsou konfigurace shodné.)

Dále musí platit (`LV0.FOR:1355-1370`):

3. **ne** současně P(5) = 0, P(6) = 0 a P(40) = 0
4. `|MODE| ≤ 3` nebo `|MODE| = 6`; pro `|MODE| = 8` se jde jinou větví
5. `IRE < 0` (`LV0.FOR:862`) — jinak se nevytvoří ani `OUTIN`

Pro DSD samotné platí navíc pravidlo z hlavičky `LV6.FOR:109`: **je-li P(5) ≠ 0,
žádný závěrečný výpis se neobjeví.** Chcete-li ho pro čistě DSD fit, nastavte
P(5) = 0.

> **Ani splnění všech podmínek výpis nezaručuje.** Za nimi následuje ještě
> spleť příznaků `INL`, `NELEM0` a `NCHK` (`LV0.FOR:1369-1400`), do níž zasahuje
> i elektrodový element. Změřeno na korpusu: **blok se objevil u 17 z 61 testů**
> v `FITTESTS/CKT.O`, a všech 17 má P(35) = 4 s NELEM ∈ {7, 10, 32, 36}.
>
> Názorný protipříklad jsou `OTSCOMP1` a `OTSCOMP2`. Mají shodné P(35) = 4,
> NELEM = 10, MODE = −1 a liší se **jen elektrodovým elementem**
> (NDE3 = 6 vs. 10) — a přesto blok vypíše pouze `OTSCOMP2`. Podobně
> `OTSTMOY.CD` splňuje NELEM = 7 i P(7) = 1, a blok stejně nedostane.
>
> Nevypíše-li se blok, zkontrolujte nejdřív NELEM. Pokud sedí, jde
> pravděpodobně o tuto nedokumentovanou závislost a spolehlivá cesta je
> vyjít z konfigurace některého ze 17 fungujících testů (viz část 12).

### Co blok obsahuje

```
        XXM1           XX1            XX2          XX3
        AVTAU          XRAT           RN           AIN

        RHOO           EDINF          EINF         E0
        TAO            ETAU           ECINF        EC0
```

| Veličina | Vzorec / význam |
|----------|-----------------|
| `XXM1` | ⟨x⁻¹⟩ — normalizovaný moment distribuce, x = τ/τ₀ |
| `XX1` | ⟨x⟩ |
| `XX2` | ⟨x²⟩ |
| `XX3` | ⟨x³⟩ |
| `AVTAU` | τ₀·⟨x⟩ — střední relaxační čas |
| `XRAT` | ε<sub>C0</sub>/ε<sub>C∞</sub> — poměr limitních hodnot |
| `RN`, `AIN` | normalizační ukazatele, viz níže |
| `RHOO` | ρ₀ (obvykle P(5)) |
| `TAO` | τ₀ (obvykle P(6)) |
| `ETAU` | ε<sub>τ</sub> = τ₀·REFI/CELCAP, kde REFI = ρ₀/(P(1)+ρ₀)² |
| `ECINF` | ε<sub>C∞</sub> = ε<sub>τ</sub>/⟨x⁻¹⟩ |
| `EC0` | ε<sub>C0</sub> = ε<sub>τ</sub>·⟨x⟩ |
| `EDINF` | ε<sub>D∞</sub> = P(2)/CELCAP (nebo přímo P(2), je-li ATEMP < 0) |
| `EINF` | ε<sub>∞</sub> = ε<sub>D∞</sub> + ε<sub>C∞</sub> |
| `E0` | ε<sub>0</sub> = ε<sub>D∞</sub> + ε<sub>C0</sub> |


### Ověřený příklad

Běh `FITTESTS/CKT.O/OTSM1M10` (P(35) = 4, NELEM = 10, MODE = −1) dá:

```
        XXM1           XX1            XX2          XX3
        AVTAU          XRAT           RN           AIN
     1.66667E-01    6.00000E+01    1.00800E+04    3.32641E+06
     2.40083E-04    1.00000E+01    0.00000E+00    0.00000E+00

        RHOO           EDINF          EINF         E0
        TAO            ETAU           ECINF        EC0
     5.42157E+07    2.70012E+01    3.20027E+01    7.70159E+01
     4.00138E-06    8.33576E-01    5.00146E+00    5.00146E+01
```

Kontrola vztahů z tabulky výše na těchto číslech:

| Vztah | Dosazení | Výsledek |
|-------|----------|----------|
| ε<sub>C∞</sub> = ε<sub>τ</sub>/⟨x⁻¹⟩ | 0,833576 / 0,166667 | 5,0015 ✓ |
| ε<sub>C0</sub> = ε<sub>τ</sub>·⟨x⟩ | 0,833576 · 60 | 50,015 ✓ |
| XRAT = ε<sub>C0</sub>/ε<sub>C∞</sub> | 50,015 / 5,0015 | 10,00 ✓ |
| ε<sub>∞</sub> = ε<sub>D∞</sub> + ε<sub>C∞</sub> | 27,0012 + 5,0015 | 32,0027 ✓ |
| ε<sub>0</sub> = ε<sub>D∞</sub> + ε<sub>C0</sub> | 27,0012 + 50,015 | 77,016 ✓ |
| AVTAU = τ₀·⟨x⟩ | 4,00138×10⁻⁶ · 60 | 2,4008×10⁻⁴ ✓ |

RN a AIN jsou zde nulové, protože při NCH = 4 se momenty berou z `DISTEL`,
ne z `GEDAE` — a nula podle manuálu znamená „nepočítáno".

Momenty počítá `GEDAE` pětinásobnou kvadraturou pro ICHG = −1…3
(`LV6.FOR:648-680`) — proto je běh s tímto výpisem znatelně pomalejší.

### RN a AIN jako diagnostika cutoffu

Tohle je nejužitečnější a nejméně zřejmá věc v celém výpisu. Manuál
(str. 5-32) a hlavička `LV6.FOR:118-141`:

- RN a AIN jsou poměry normalizace spočítané kvadraturou ku normalizaci
  spočítané uzavřeným vzorcem. **Ideál je 1.**
- Když se počítá jen CSD nebo jen DSD, jsou si rovné. Při obou je AIN pro CSD
  a RN pro DSD (počítá se jako druhá). Nula znamená „nepočítáno".
- **Pro GED (NCH = 1) s γ₁ = 1** je vzorec exaktní i při konečném cutoffu,
  takže RN a AIN vyjdou v podstatě 1.
- **Pro γ₁ ≠ 1** je odchylka od 1 přímou mírou vlivu konečného cutoffu.
  Zvětšujte \|U₁\|, dokud se nepřiblíží 1.
- **Pro CD a HN (NCH = 2, 3)** je normalizace exaktní jen v limitě bez cutoffu,
  takže RN > 1 vždy. Hodnota **1,01 už znamená dobrou aproximaci**. Očekávaný
  rozsah U₁ je 5 až 25; při U₁ > 30 je cutoff prakticky bez vlivu.

Výpočet pro NCH = 2, 3 je velmi pomalý pro U₁ > 15, hodně datových bodů
a IGACC > 3. Manuál: **IGACC = 2 nebo 3 dává pro HN a CD dostatečnou přesnost.**
Záporné U₁ u NCH = 2, 3 nemění výsledek, ale zapne výpis postupu na obrazovku
(`LV6.FOR:791`).

### Varovné hlášky bloku

| Hláška | Příčina |
|--------|---------|
| `WRONG INPUT: NO DEBYE HERE` | P(7) = 1 **a zároveň** P(9) = 1 (`LV0.FOR:1512`) |
| `***** SINCE P(1) NOT ZERO, FINAL EINF= EDINF ******` | P(1) ≠ 0 — ε<sub>∞</sub> nelze rozložit |

---

## 10. Speciální operace jen pro O-obvod

Všechny vyžadují **MAXFEV = 1 nebo 2** (`LV0.FOR:1050`) a probíhají při zápisu
souboru `OUTIN`. Po použití **parametr vždy vraťte na nulu**, jinak se
transformace bude opakovat.

### 10.1 Odečtení σ₀ (stejnosměrné vodivosti)

Nastavte `P(36) ≠ 0`, `DFIT = Y`, `MAXFEV = 1`. Kód provede (`LV0.FOR:1072-1076`):

```fortran
Y(I)  = Y(I) - 1.D0/P(5)                          ! odečte sigma0 = 1/rho0
Y(IM) = -Y(I)/(TUPI*FREQ(I)*EPV)                  ! prepocet na eps''
```

kde `EPV = 8.85418782D-14` F/cm. **Jedním během vznikne obojí** — reálný sloupec
je σ′ bez σ₀ a imaginární −ε″. Výsledná křivka ε″ je pro vhodný frekvenční
rozsah pěkně vrcholová.

> `README.txt` popisuje tento postup jako dvoukrokový (nejdřív DFIT = Y, pak
> nový běh s DFIT = E). Kód v současné verzi počítá oba sloupce naráz.
> Pořadí operací má jeden důsledek, který stojí za pozor: `Y(IM)` se počítá
> z **už zmenšeného** `Y(I)`, ne z původní hodnoty.

Postup předpokládá, že jste předtím fitem odstranili elektrodové jevy.

### 10.2 Škálování frekvenční osy — P(38)

| Podmínka | Efekt |
|----------|-------|
| `0 < P(38) ≤ 10¹⁷` | `FREQ ← P(38)·FREQ` |
| `P(38) > 10¹⁷` | `FREQ ← P(6)·FREQ` — tj. ω·τ₀; nový fit dá τ₀ = 1 |

Používá se po převodu dat na úhlovou frekvenci (volba Data → Hz/rad v LEVMW).
V nepřítomnosti elektrodových jevů a P(2) bude nová fitovaná hodnota P(6)
jednotková.

### 10.3 Škálování vodivosti — P(39)

Nejdřív převeďte data na komplexní vodivost (`DFIT = Y`, `MAXFEV = 1`, běh),
pak nastavte P(39):

| Podmínka | Efekt |
|----------|-------|
| `P(39) < 0` a `\|P(39)\| ≤ 10¹⁷` | `Y ← \|P(39)\|·Y` |
| `\|P(39)\| > 10¹⁷`, `DFIT = Y` | `Y ← P(5)·Y` — nový fit dá P(5) = 1 |

> **Nesoulad mezi `README.txt` a kódem.** README u obou transformací uvádí práh
> **10¹⁸**, kód testuje **10¹⁷** (`LV0.FOR:1054,1057,1065,1067`). Pro hodnoty
> mezi 10¹⁷ a 10¹⁸ se chování liší od dokumentace. Bezpečné je držet se
> hodnot buď výrazně pod 10¹⁷, nebo výrazně nad 10¹⁸.

### 10.4 Odečtení GD a CINF

Manuál, str. 5-32: nastavte `MAXFEV = 2`, `NFREE = 3` pro P(12), spusťte.
Vliv parametru se odečte z dat do `OUTIN`.

Důležité: **odečítání GD a CINF se musí dělat na úrovni Y**, odečítání RINF
na úrovni Z.

### 10.5 Postupné vypínání složek modelu

Diagnostický postup z `README.txt`, který ukáže, čím která část modelu
k odezvě přispívá:

1. Proveďte fit a prohlédněte grafy (Run → Launch LEVMVIEW).
2. Vraťte se na hlavní stránku a nastavte **MAXFEV = 0**.
3. Vynulujte jednu nebo více nebulkových složek:
   - `P(2) = 0` — odstraní C<sub>∞</sub>
   - `P(25) = 0` — odstraní elektrodový DCE (SCPE)
4. Spusťte a znovu otevřete LEVMVIEW.

Pro model CK1S tak uvidíte vedle sebe odezvu CK1S a odezvu K1S, CK1 nebo K1
podle toho, co jste vypnuli. (Písmena: **C** = s C<sub>∞</sub>, **K1** = základní
CSD1 model, **S** = se SCPE elektrodou.)

---

## 11. Praktický postup

### Nutná nastavení

| Volba | Hodnota | Proč |
|-------|---------|------|
| N (řádek 3) | **≥ 40** | O-obvod používá P(35), P(40) a dál |
| IRE | **< −9** (typicky −11) | vytvoří `OUTIN`, bez něj nejde test idempotence ani závěrečný dielektrický výpis |
| CELCAP | správná hodnota | **povinné při ATEMP < 0**, kdy jsou P(2) a P(15) na úrovni ε; obecně vždy pro reálná data |
| IOP | −4 až −5 pro experimentální data, −8 až −9 pro velmi přesná | volí toleranci konvergence |
| IGACC | 2 nebo 3 | přesnost kvadratur; víc je u HN/CD velmi pomalé |
| P(35), P(40) | konzistentní s nenulovými parametry | jinak *„a run error will probably occur"* |

### ATEMP a jednotky

| ATEMP | P(2), P(15) jsou |
|-------|------------------|
| ≥ 0 a < 1000 | kapacity ve faradech |
| < 0 | dielektrické konstanty (ε<sub>∞</sub>, ε<sub>0</sub>) — **CELCAP musí být správně** |

### Typický pracovní postup

1. **Vyberte model** — pro pevné elektrolyty začněte NCH1 = 4 s NELEM = 10
   (KWW, rychlé) nebo 36 (KWW s cutoffem, přesnější).
2. **Nastavte MODE** — 1 pro K0/KD odezvu, −1 pro K1.
3. **První fit bez elektrod** — NDE3 = 0, jen objemová část. Zkontrolujte
   `SD RC` a relativní směrodatné odchylky.
4. **Přidejte elektrodovou složku** — NDE3 = 26 (SCPE) je kanonická volba.
5. **Ověřte cutoff** — v závěrečném výpisu se podívejte na RN a AIN. Jsou-li
   výrazně nad 1, zvětšete \|U₁\|.
6. **Test idempotence** — spusťte `OUTIN` znovu a porovnejte (viz `TESTING.md`).
7. **Diagnostika složek** — postup 10.5.

### Časté chyby a jejich hlášky

| Hláška | Příčina |
|--------|---------|
| `WRONG VALUE OF P(35)` / `P(40)` | NCH mimo rozsah 0–6 |
| `WRONG INPUT: NO DEBYE HERE` | P(7) = 1 a P(9) = 1 zároveň |
| `!!!!!** SINGULAR MATRIX - BEWARE **!!!!!` | příliš mnoho volných parametrů, nebo parametr, který model nepoužívá |
| `#!!!!! BAD CORRELATION MATRIX !!!!!#` | volný parametr, který se ve fitovací funkci vůbec nevyskytuje — typicky nesoulad mezi NCH a nenulovými parametry |
| běh je extrémně pomalý | NCH = 2 nebo 3 s U₁ > 15, nebo IGACC > 3 |
| parametr klesá k nule, fit nekonverguje | parametr potřebuje jít do záporu — nastavte NFREE = **2** (ne 1) |

---

## 12. Reprezentativní testovací soubory

Z `FITTESTS/CKT.O`. Sloupec „Co ukazuje" vychází z hlavičky a konfigurace.

| Soubor | Konfigurace | Co ukazuje |
|--------|-------------|------------|
| `OTST1`, `OTST2` | P35 = 4, NELEM = 6, MODE = 0 | základní CSD fit s HN elementem |
| `OTSTHN` | P35 = 3, MODE = −1 | HN **s cutoffem** přímou integrací |
| `OTSTHN.FUL` | P35 = 3, MODE = −6 | totéž, CD/HN dávají CSD1a |
| `OTSTCD.A` | NELEM = 6, MODE = −1 | Cole-Davidson |
| `OTSTIME.CD`, `OTSTIME.HN` | MAXFEV = 10 / 4 | přechodová odezva |
| `OTSTKWWE.GED` | P35 = 1, γ = 0,586 | GED s neceločíselným γ |
| `OTSP1Y32/35/36` | NELEM = 32 / 35 / 36 | tři varianty KWW s cutoffem |
| `OTSZDRTA.K*` | NELEM = 37, DATTYP = R | přímý výpočet KWW DRT |
| `OTSCOMP1`, `OTSCOMP2` | NDE3 = 6 / 10, MODE = −1 | **dva CSD elementy v sérii** (viz 6) |
| `OTSTNOW1-4.422` | P35 = 4, P40 = 4 | současný CSD i DSD fit |
| `OTSM1E25`, `OTSM1E27` | P35 = 0, P40 = 4, NELEM = 25 / 27 | jen DSD, efektivní medium a PCPE |
| `OMFTMP1`, `OMFTMP2` | NELEM = 33, DATTYP = R, MAXFEV = 255 | **Arrheniovský fit** teplotní závislosti (ρ₀ a τ₀) |
| `le225e35`, `le225m36`, `225Z36el` | NELEM = 35 / 36, NDE3 = 26 | model UN s β₁ = 1/3 a SCPE elektrodou |
| `OTSP1Y35.0`, `OTSZDRT.K0` | MAXFEV = 0 | jen simulace, bez fitu |

---

## 13. Úskalí

- **P(10) a P(20) mají dvojí význam.** Při NCH = 4 je to NELEM, jinak φ₂.
  Záměna je nejčastější zdroj hlášky o špatné korelační matici.
- **P(8) a P(18) jsou v režimu DE tiše vynulovány.** Jejich hodnota ve vstupu
  se ignoruje (kromě NELEM = 9).
- **P(11) a P(12) jsou vodivosti**, přestože komentář ve zdroji u P(11) mluví
  o odporu. Manuál to potvrzuje.
- **P(38) = 64 mění význam P(11) na kapacitu** a zabírá P(41)–P(45).
  Nekombinujte s Ngai coupling modelem.
- **Práh 10¹⁷ vs 10¹⁸** u P(38) a P(39) — README a kód se rozcházejí (10.3).
- **Transformace P(36), P(38), P(39) se aplikují při každém běhu**, dokud je
  parametr nenulový. Po použití je vždy vraťte na nulu.
- **Bez `IRE < −9` nedostanete `OUTIN` ani závěrečný dielektrický výpis.**
- **Pro čistě DSD fit musí být P(5) = 0**, jinak se závěrečný výpis neobjeví.
- **P(33) (efektivní medium) nemá v `FITTESTS` jediný test.** Před použitím
  ověřte na simulovaných datech.
