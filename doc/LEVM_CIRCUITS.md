# LEVM Circuit Topologies

Přehled podporovaných obvodových modelů v programu LEVM v moderní textové notaci.
Ověřeno proti zdrojovému kódu `LV4.FOR` (A–K), `LV6.FOR` (O), `LV7.FOR` (R),
`LV8.FOR` (S), `LV9.FOR` (T) a `LV5.FOR` (DCE).

## Legenda

| Symbol | Význam |
|--------|--------|
| `-` | Sériové spojení (**pomlčka, ne minus**) |
| `\|` | Paralelní spojení |
| `DCEn` | Distribuovaný element (slot n), typ určen parametrem NDE |
| `(RA,CA)` | Augmentovaný DCE — DCE paralelně s článkem RA\|CA (viz `SDEA`) |
| `DAE` | Distribution of Activation Energies |
| `DRT` | Distribution of Relaxation Times |

Notace platí v blocích **Topologie**. Uvnitř vzorců pro impedanci (tam, kde se
vyskytuje `/`, `sqrt`, `coth`, mocniny) jde o běžnou algebru, kde `+` znamená
sčítání a `−` odčítání.

## Přehled obvodů

Obvod se volí znakem **FUN** na řádku 2 vstupního souboru. Dispatch je v `MODEL`
(`LV2.FOR:887`). Platné hodnoty jsou pouze tyto:

| FUN | Podprogram | Popis | Testy |
|-----|-----------|-------|-------|
| A | `ASUB` | 6 podobvodů v sérii | `FITTESTS/CKT.A` |
| B | `BSUB` | vnořené RC s DCE | `FITTESTS/CKT.B` |
| C | `CSUB` | 2 augmentované + 1 regulární DCE | `FITTESTS/CKT.C` |
| D | `DSUB` | vnořené RC/DCE s DAE | `FITTESTS/CKT.D` |
| E | `ESUB` | 5 DCE slotů | `FITTESTS/CKT.E` |
| F | `FSUB` | 4 sériové větve paralelně | `FITTESTS/CKT.F` |
| G | `GSUB` | obecná reálná fitovací funkce (není obvod) | `FITTESTS/CKT.G` |
| H | `HSUB` | blocking/conductive PNP model | `FITTESTS/CKT.H` |
| I | `ISUB` | přechodová EDAE odezva v čase | `FITTESTS/CKT.I` |
| J | `JSUB` | 1 augmentovaný + 3 regulární DCE | `FITTESTS/CKT.J` |
| K | `KSUB` | diskrétní DRT, až 11 párů | `FITTESTS/CKT.K` |
| O | `OSUB` | CSD/DSD disperzní odezva | `FITTESTS/CKT.O` |
| R | `RSUB` | diskrétní DRT, až 19 párů | `FITTESTS/CKT.R` |
| S | `SSUB` | elektrochemický superkondenzátor | `FITTESTS/CKT.S` |
| T | `TSUB` | **stub — program se ukončí** | `FITTESTS/CKT.T` |

Neplatná hodnota FUN vede k `*** PROGRAM TERMINATED: <x> IS INVALID MODEL`.

> **T-obvod není implementován.** `TSUB` (`LV9.FOR:9`) obsahuje pouze
> `PRINT *,'###### T CIRCUIT: CURRENTLY DUMMY! ####'` a `STOP`.

## Vnější struktura — pozor, není univerzální

Obvody A, B, C, D, F, J mají shodnou vnější obálku:

```
L - ( vnitřní_struktura | (RP|CP) )
```

| Parametr | Popis | Vypnutí |
|----------|-------|---------|
| P(28) = RP | Paralelní odpor (vnější) | RP = 0 |
| P(29) = CP | Paralelní kapacita (vnější) | CP = 0 |
| P(30) = L | Sériová indukčnost | L = 0 |

**U ostatních obvodů mají P(28)–P(30) jiný význam:**

| Obvod | P(28) | P(29) | P(30) |
|-------|-------|-------|-------|
| E | R3 (vnitřní odpor) | CP | L — **RP neexistuje** |
| G | XLG (dolní mez GDAE) | XUG (horní mez GDAE) | XTA (přepínač znaménka v exponentech) |
| H | RP | CP | ALINDUC (ICH=0) **nebo** PHIN (ICH<3 nebo ≥5) |
| K | **GP = vodivost**, ne odpor | CP | \|P(30)\| = L; **P(30) < 0 ⇒ sériový odpor** |
| O | R3 (elektrodová část) | C3 (elektrodová část) | sériový prvek, roli určuje P(37) |
| R | (součást DRT páru) | CP | roli určuje P(37): R / G / L / C |
| S | nepoužito | nepoužito | nepoužito |
| T | nepoužito | nepoužito | nepoužito |

U obvodů O a R určuje **P(37)** roli P(30):

| P(37) | Význam P(30) |
|-------|--------------|
| > 0 | sériový odpor R<sub>s</sub> (u R-obvodu použijte P(37) = 10) |
| = 0 | paralelní vodivost G (pouze R-obvod) |
| < 0 a > −9 | sériová indukčnost L (typicky P(37) = −4) |
| < −9 a > −19 | sériová kapacita C<sub>s</sub> (typicky P(37) = −16) |

---

## A Circuit (ASUB)

**Popis:** 6 subcircuits v sérii, paralelně s RP/CP, v sérii s L.

**Topologie:**
```
L - ((
    (R1|C1) - (R2|C2) - (R3|C3) -
    DCE1 -
    (RA|CA|DCE2) -
    ((R5 - DCE3 - C5) | (R4|C4))
) | (RP|CP))
```

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE3 |
| 2 | C1 | 17 | TDE3 |
| 3 | R2 | 18 | UDE3 |
| 4 | C2 | 19 | PDE3 |
| 5 | R3 | 20 | NDE3 |
| 6 | RDE1 | 21 | C3 |
| 7 | TDE1 | 22 | RA |
| 8 | UDE1 | 23 | CA |
| 9 | PDE1 | 24 | R4 |
| 10 | NDE1 | 25 | C4 |
| 11 | RDE2 | 26 | R5 |
| 12 | TDE2 | 27 | C5 |
| 13 | UDE2 | 28 | RP |
| 14 | PDE2 | 29 | CP |
| 15 | NDE2 | 30 | L |

Poznámka: C5 je v **sérii** s R5 a DCE3, ne paralelně.

---

## B Circuit (BSUB)

**Popis:** Vnořené RC s distribuovanými elementy.

**Topologie:**
```
L - ((
    R1 - (RA|CA|DCE1) - (
        (R2 - DCE2 - ((R3 - DCE3) | C3)) | (C2 - DCE4)
    )
) | (RP|CP))
```

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE3 |
| 2 | RA | 17 | TDE3 |
| 3 | CA | 18 | UDE3 |
| 4 | R2 | 19 | PDE3 |
| 5 | C2 | 20 | NDE3 |
| 6 | RDE1 | 21 | RDE4 |
| 7 | TDE1 | 22 | TDE4 |
| 8 | UDE1 | 23 | UDE4 |
| 9 | PDE1 | 24 | PDE4 |
| 10 | NDE1 | 25 | NDE4 |
| 11 | RDE2 | 26 | R3 |
| 12 | TDE2 | 27 | C3 |
| 13 | UDE2 | 28 | RP |
| 14 | PDE2 | 29 | CP |
| 15 | NDE2 | 30 | L |

### Korelovaná přeparametrizace při UDE3 < 0

`LV4.FOR:195-215`. Pokud je **P(18) = UDE3 záporné**, změní se význam několika parametrů:

- **UDE3 < −9:** P(4) a P(5) se stanou *normalizovanými* hodnotami:
  `R2 = P(1)·P(4)`, `C2 = P(29)·P(5)`.
- **NDE3 = 13 nebo 15** (obecná homogenní difúze): navíc
  `C3 = P(27)/R2`, `R3 = P(26)/C3`, a při UDE3 ≠ 2 také `RDE3 = P(16)·R2`.
  Pokud je při tom R2 = 0 nebo P(27) = 0, program skončí hláškou `WRONG INPUTS`.

---

## C Circuit (CSUB)

**Popis:** Dva augmentované DCE + jeden regulární DCE.

**Topologie:**
```
L - ((
    R1 -
    (((R3|C3) - (RA1|CA1|DCE1)) | (R2 - C2)) -
    ((R4 - DCE3 - (R5|C5)) | (RA2|CA2|DCE2))
) | (RP|CP))
```

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE3 |
| 2 | R2 | 17 | TDE3 |
| 3 | C2 | 18 | UDE3 |
| 4 | R3 | 19 | PDE3 |
| 5 | C3 | 20 | NDE3 |
| 6 | RDE1 | 21 | RA1 |
| 7 | TDE1 | 22 | CA1 |
| 8 | UDE1 | 23 | RA2 |
| 9 | PDE1 | 24 | CA2 |
| 10 | NDE1 | 25 | R4 |
| 11 | RDE2 | 26 | R5 |
| 12 | TDE2 | 27 | C5 |
| 13 | UDE2 | 28 | RP |
| 14 | PDE2 | 29 | CP |
| 15 | NDE2 | 30 | L |

---

## D Circuit (DSUB)

**Popis:** Vnořené R-C/DCE obvody s DAE (rozdělení aktivačních energií).

**Topologie:**
```
L - ((
    R1 -
    ((RA|CA) | DAE) -
    ((R2 - (R3|C2)) | (((R4 - DCE4)|C4) - ((R5 - DCE5)|C5)))
) | (RP|CP))
```

**DAE typy** (`LV4.FOR:508`):
- P(15) = DCH = 0: Exponenciální DAE (EDAE)
- P(15) ≠ 0: Gaussovské DAE (GDAE)

U EDAE mají hodnoty > 800 speciální význam „fixováno / odvozeno":
`PHI1 > 800 ⇒ PHI1 = PHI2`; `U2 > 800 ⇒ U2 = U1`; `PHI2 > 800 ⇒ PHI2 = −PHI1`.
Záporné TDAE znamená `TDAE = |TDAE|` s násobitelem `exp(−U1)`.
U GDAE platí obdobně `XL > 800 ⇒ XL = XU`.

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE4 |
| 2 | RA | 17 | TDE4 |
| 3 | CA | 18 | UDE4 |
| 4 | RDAE | 19 | PDE4 |
| 5 | TDAE | 20 | NDE4 |
| 6 | U1 (log) | 21 | RDE5 |
| 7 | U2 (log) | 22 | TDE5 |
| 8 | φ1 | 23 | UDE5 |
| 9 | φ2 | 24 | PDE5 |
| 10 | R2 | 25 | NDE5 |
| 11 | C2 | 26 | R5 |
| 12 | R3 | 27 | C5 |
| 13 | R4 | 28 | RP |
| 14 | C4 | 29 | CP |
| 15 | DCH | 30 | L |

Pozor: C2 (P(11)) je paralelně s **R3**, ne s R2 — větev je `R2 - (R3|C2)`.

---

## E Circuit (ESUB)

**Popis:** 5 DCE slotů — maximální flexibilita pro distribuované elementy.

**Topologie** (`LV4.FOR:812-846`):
```
L - ((
    DCE1 - R1 - (
        DCE4 | ( DCE2 - R2 - ( DCE5 | (R3 - DCE3) ) )
    )
) | CP)
```

**Poznámka:** Tento obvod nemá RP — P(28) je vnitřní odpor R3. Vnější paralelní
prvek je pouze CP.

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | RDE1 | 16 | RDE4 |
| 2 | TDE1 | 17 | TDE4 |
| 3 | UDE1 | 18 | UDE4 |
| 4 | PDE1 | 19 | PDE4 |
| 5 | NDE1 | 20 | NDE4 |
| 6 | RDE2 | 21 | RDE5 |
| 7 | TDE2 | 22 | TDE5 |
| 8 | UDE2 | 23 | UDE5 |
| 9 | PDE2 | 24 | PDE5 |
| 10 | NDE2 | 25 | NDE5 |
| 11 | RDE3 | 26 | R1 |
| 12 | TDE3 | 27 | R2 |
| 13 | UDE3 | 28 | R3 |
| 14 | PDE3 | 29 | CP |
| 15 | NDE3 | 30 | L |

---

## F Circuit (FSUB)

**Popis:** Čtyři sériové větve zapojené vzájemně paralelně.

**Topologie:**
```
L - (( Z3 | Z4 | Z5 | Z6 ) | (RP|CP))

Z3 = (RA|CA|DCE1) - R1 - C1
Z4 = ((R3 - C4) | C3) - R2 - C2
Z5 = (R6|C6) - (R7|C7) - C5
Z6 = (R9|C9) - DCE2 - R8
```

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RA |
| 2 | C1 | 17 | R3 |
| 3 | R2 | 18 | C3 |
| 4 | C2 | 19 | C4 |
| 5 | CA | 20 | C5 |
| 6 | RDE1 | 21 | R6 |
| 7 | TDE1 | 22 | C6 |
| 8 | UDE1 | 23 | R7 |
| 9 | PDE1 | 24 | C7 |
| 10 | NDE1 | 25 | R8 |
| 11 | RDE2 | 26 | R9 |
| 12 | TDE2 | 27 | C9 |
| 13 | UDE2 | 28 | RP |
| 14 | PDE2 | 29 | CP |
| 15 | NDE2 | 30 | L |

> **Pozor na chybu ve zdroji.** Testy na nulu pro větve R7\|C7 a R9\|C9 se ptají
> na nesprávný parametr (`LV4.FOR:912-913`: `RC70` testuje **R2** místo R7,
> `RC90` testuje **R3** místo R9). Důsledek: nastavíte-li R7 = 0 (nebo R9 = 0)
> při nenulovém R2 (resp. R3), příslušný kondenzátor se ze součtu ztratí místo
> aby zbyl samotný. Chcete-li čistý kondenzátor, použijte velmi malý nenulový odpor.

---

## G Circuit (GSUB)

**Popis:** Není obvod. Je to obecná fitovací funkce jedné reálné proměnné,
`Y = F(x)`, kde x je hodnota z frekvenčního sloupce (může to být frekvence, čas
nebo teplota).

**Vyžaduje DATTYP = R nebo I.** Při DATTYP = C skončí hláškou
`G CKT NOT APPLICABLE FOR COMPLEX FIT`.

Výsledek je součet osmi nezávislých členů; každý se vypne nastavením svého
prvního parametru na 0:

```
YT = YT1 + YT2 + YT3 + YT4 + YT5 + YT6 + YT7 + YT8   (+ volitelný člen dle P(40))
```

| Člen | Parametry | Tvar |
|------|-----------|------|
| YT1 | P(1)–P(5) = AY,BY,CY,DY,EY | EY≠0: `AY·x^BY·exp(CY(DY−x)/(DY(x−EY)))` (exponent ořezán na ±100); EY=0: `AY·exp(−x/DY)` |
| YT2 | P(6)–P(8) = FY,GY,HY | `FY(GY−x)/(GY(x−HY))` |
| YT3 | P(9),P(10) = PXY,QY | `PXY·x^QY` |
| YT4 | P(11)–P(13) = RY,SY,TY | `TY(RY + SY·x)` |
| YT5 | P(14)–P(16) = WY,YY,ZY | tvar dle P(40), viz níže |
| YT6 | P(17),P(18) = AX,BX | `AX·exp(−x/BX)` |
| YT7 | P(19)–P(21) = CX,DX,EX | `CX·exp(−(x/DX)^EX)` (stretched exponential) |
| YT8 | P(22)–P(24) = FX,GX,HX | `FX·GX·x^(GX−1)·HX^(−GX)·exp(−(x/HX)^GX)` (Weibull) |

**YT5 podle P(40):**

| P(40) | Tvar YT5 |
|-------|----------|
| 0 | `WY/(1 + (x/YY)^ZY)` |
| > 0 | `WY/(1 + x/YY)^ZY` |
| −1 | `WY/((1 + (x/YY)^ZY)^AX)` |
| −2 | `WY/((1 + x/BX + (x/YY)^ZY)^AX)` |

Alternativně: je-li WY ≠ 0, XTA = P(30) > 0 a YY = 0, počítá se **GDAE přechodová
odezva** kvadraturou přes meze −P(28)…P(29) s parametry P(25)=XI, P(26)=GAM,
P(27)=DEL.

**Přídavné fitovací funkce podle P(40)** (přičítají se k YT):

| P(40) | Přídavný člen |
|-------|---------------|
| −4 | polynom `P41 + x(P42 + x(P43 + x·P44))` |
| −8 | polynom v 1/x: `P41 + x⁻¹(P42 + x⁻¹(P43 + x⁻¹P44))` |
| −16 | `P41·sin(P42·x^P43)` |
| −32 | `P41·atan(P42·x^P43)` |
| −64 | `P41·tanh(P42·x^P43)` |
| −128 | `P46·x^P47·exp(−(P48·x)^P49)` (nahrazuje YT, nepřičítá se) |

**XTA = P(30) < 0** přepíná v členech YT1, YT6, YT7 exponenty z `exp(−x/D)` na
`exp(−x·D)`.

**Zbývající parametry:** P(25) = XI, P(26) = GAM, P(27) = DEL, P(28) = XLG,
P(29) = XUG, P(30) = XTA. Používají se P(41)–P(44) a P(46)–P(49).

---

## H Circuit (HSUB)

**Popis:** Blocking/conductive model — Poisson-Nernst-Planck (PNP) kontinuální
difúze a generačně-rekombinační modely. Nejsložitější obvod v LEVM.

**Význam parametrů P(11)–P(30) závisí na hodnotě `ICH = DINT(P(25))`.**
Není možné uvést jednu parametrovou tabulku. Viz manuál str. 4-9, 4-10
a 5-15 až 5-19.

**Pevné sloty (nezávislé na ICH):**

| P | Prvek |
|---|-------|
| 1–5 | DCE2: RDE2, TDE2, UDE2, PDE2, NDE2 |
| 6–10 | DCE3: RDE3, TDE3, UDE3, PDE3, NDE3 |
| 21 | R2 |
| 22 | C2 |
| 23 | R3 |
| 24 | C3 |
| 25 | **ICH** — volba PNP modelu |
| 28 | RP |
| 29 | CP |

**Volba ICH:**

| ICH | Model | IBC |
|-----|-------|-----|
| 0 | bez PNP jádra; P(30) = ALINDUC (sériová indukčnost) | — |
| 1 | blokující, jednodruhová mobilita | 1 |
| 2 | blokující, generace-rekombinace přes EKG/EKR | 2 |
| 3 | částečně blokující | 1 |
| 4 | částečně blokující, varianta | 2 |
| 5 | dvoumobilitní PNP s RH0/RHIN | — |

Při ICH ≠ 0 je **P(30) = PHIN**, ne indukčnost.

**Znaménko P(11)** volí vnitřní parametrizaci (`LV4.FOR:1271`):
`P(11) ≥ 0 ⇒ INN = 0`, délka vzorku `ELL = P(27)`;
`P(11) < 0 ⇒ INN = 1`, `ELL = P(12)`, `AN0 = −P(11)`.

**Příklady významů při ICH = 1 (IBC = 1):** P(17) = RA, P(18) = CA, P(19) = R1
nebo EMUN (mobilita) podle INN, P(20) = C1, P(26) = EEPS, P(27) = ELL.
**Při IBC = 2:** P(13) = EZN, P(14) = EZP, P(15) = EKG, P(16) = EKR,
P(17) = EMUN, P(18) = EMUP, P(19) = ELL, P(20) = EEPS.
**Při ICH = 5:** P(26), P(27) = RH0(1..2), P(17), P(18) = RHIN(1..2),
P(12) = ELL, P(13) = PIZ, P(14) = PIM, P(15) = AKRAT, P(16) = XI, P(29) = XIA(2).

**Topologie (větev P(40) < 4):**
```
ZA = C3 | (R3 - DCE3)
ZB = C2 | (R2 - ZA - DCE2)
ZC = Z_PNP | (RA|CA)            (při ICH = 0 jen RA|CA)
ZT = ALINDUC - ( (ZC - ZB) | (RP|CP) )
```

**Větev P(40) ≥ 4** používá jiné složení a na závěr přidává sériový prvek
řízený P(37), s hodnotou v **P(33)**:

| P(37) | Význam P(33) |
|-------|--------------|
| ≥ 0 | sériový odpor R<sub>s</sub> |
| < 0 a > −9 | sériová indukčnost (typicky P(37) = −4) |
| < −9 a > −19 | sériová kapacita (typicky P(37) = −16) |

**ATEMP** je teplota; při ATEMP < 0 se použije `TEMP = −0.1·ATEMP`.

> **Pozor — RP/CP u H-obvodu nefunguje.** Blok ve větvi P(40) < 4
> (`LV4.FOR:1784-1793`) je poškozený v obou směrech:
>
> ```fortran
> IF(RP0) THEN                          ! RP = 0 a CP <> 0
>   YC = CP*IOMEGA
> ELSE
>    IF(RP1) THEN                       ! RP = 0 (a CP = 0)
>       YC = (0,0)
>       YC = (CP*RP*IOMEGA + 1.D0)/RP   ! delení nulou
>    ENDIF                              ! chybí ELSE pro RP <> 0
> END IF
> ```
>
> Při RP ≠ 0 se neprovede nic a YC si ponese hodnotu z předchozí iterace;
> při RP = 0 a CP = 0 se dělí nulou. Komentář ve zdroji tuto část označuje
> jako `NO LONGER USED`. Nespoléhejte u H-obvodu na RP/CP — starší varianta
> `LV4n.FOR` tuto větev má správně.

---

## I Circuit (ISUB)

**Popis:** Není obvod. Počítá **přechodovou (transientní) EDAE odezvu v časové
doméně**. Sloupec „frekvence" obsahuje čas.

**Vyžaduje DATTYP = R.** Imaginární výstup je vždy 0.

Odezva se počítá přes neúplnou gama funkci ve `TRAGA`/`GINCO` z argumentů
`y = t/BY` a `AA = 1 − CY`, `1 − DY`:

```
F(t) = AY · TRSY(1)      kde TRSY = TRAGA(AA, t/BY, PP, QQ)
       PP = (1, EY),  QQ = (FY, 1)
```

**Parametry:**
| P | Prvek |
|---|-------|
| 1 | AY — amplituda |
| 2 | BY — časová konstanta (dělitel času) |
| 3 | CY — první exponent (AA₁ = 1 − CY) |
| 4 | DY — druhý exponent (AA₂ = 1 − DY) |
| 5 | EY — váha druhého členu |
| 6 | FY — váha prvního členu |

P(7)–P(30) se nepoužívají.

---

## J Circuit (JSUB)

**Popis:** Jeden augmentovaný a tři regulární distribuované elementy.
DCE1 je zapojen **paralelně** k celému vnitřnímu obvodu, spolu s RP a CP.

**Topologie:**
```
RL - L - ((
    R1 - DCE2 - ( (R4 - DCE4 - (R5|C5)) | (RA|CA|DCE3) )
) | (RP|CP) | DCE1)
```

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE3 |
| 2 | RA | 17 | TDE3 |
| 3 | CA | 18 | UDE3 |
| 4 | R4 | 19 | PDE3 |
| 5 | RL (sériový odpor) | 20 | NDE3 |
| 6 | RDE1 | 21 | RDE4 |
| 7 | TDE1 | 22 | TDE4 |
| 8 | UDE1 | 23 | UDE4 |
| 9 | PDE1 | 24 | PDE4 |
| 10 | NDE1 | 25 | NDE4 |
| 11 | RDE2 | 26 | R5 |
| 12 | TDE2 | 27 | C5 |
| 13 | UDE2 | 28 | RP |
| 14 | PDE2 | 29 | CP |
| 15 | NDE2 | 30 | L |

Na rozdíl od ostatních obvodů má J navíc **RL = P(5)**, čistý sériový odpor
přidávaný až úplně nakonec: `ZT = ZT + jωL + RL`.

> `JSUB` při každém posledním datovém bodě otevře, přepíše a zavře soubor
> `LNL` v pracovním adresáři (`LV4.FOR:2134-2139`). Jde o pozůstatek ladění;
> soubor lze ignorovat.

---

## K Circuit (KSUB)

**Popis:** Diskrétní **DRT** (distribuce relaxačních časů) — až 11 párů
(R, τ) v sérii, plus DCE, celé paralelně s GP a CP.

**Pro více než 11 párů použijte R-obvod.**

**Parametry:**
| P | Prvek |
|---|-------|
| 1–22 | až 11 párů: P(1),P(2) · P(3),P(4) · … · P(21),P(22) |
| 23–27 | DCE1: RDE1, TDE1, UDE1, PDE1, NDE1 |
| 28 | **GP — paralelní vodivost** (ne odpor!) |
| 29 | CP |
| 30 | \|P(30)\| = L; při P(30) < 0 je to kladný **sériový odpor** |
| 36 | ≠ 0: přesouvá CP paralelně k DCE (> 0) nebo k DRT části (< 0) |
| 38 | \|P(38)\| = 1 ⇒ **všechna τ musí být fixní**, jinak program skončí; jinak volba splinu v `DXSPS` |
| 39 | > 0: nastavuje minimální hodnotu τ |
| 40 | ≥ 0 konduktivní DRT, < 0 dielektrická DRT; \|P(40)\| volí interpretaci páru |

**Interpretace páru podle |P(40)|:**

| \|P(40)\| | Význam páru (P(k), P(k+1)) | Admitance členu |
|-------|---------------------------|------------------|
| 0, 1 | (R, C) — τ je součin R·C | `R/(1 + jωRC)` |
| 2, 3 | (R, τ) — τ přímo | `R/(1 + jωτ)` |
| 4 | (R, τ) — τ přímo, škálováno | `τ·R/(1 + jωτ)` |

Počet použitých pozic udává **ICF** (z `/CM12/`). Je-li `ATEMP < 0` a
`P(40) < 0`, jsou parametry DRT dielektrické konstanty, ne kapacity.

Je-li **P(21) < 0**, přidá se navíc paralelní RC větev z `|P(21)|` a `P(22)`.

---

## O Circuit (OSUB)

> Podrobný popis O-obvodu — všechny modely disperze, dielektrický výpis,
> speciální transformace a praktické postupy — je v samostatném dokumentu
> **`O_CIRCUIT.md`**. Zde je jen přehled.

**Popis:** Nejvýkonnější a autorem doporučovaný model pro pevné materiály.
Odděleně nebo společně počítá **CSD** (conductive-system dispersive) a **DSD**
(dielectric-system dispersive) odezvu, každou z pěti volitelných modelů.
Má nejvíc testovacích souborů (`FITTESTS/CKT.O`).

### Volba modelu

CSD volí **P(35)**, DSD volí **P(40)**. Obě používají stejné kódy (NCH):

| NCH | Model |
|-----|-------|
| 0 | žádná odezva — blok se přeskočí |
| 1 | GE — generalized exponential DAE (tři podtypy, viz níže) |
| 2 | CD — Cole-Davidson |
| 3 | HN — Havriliak-Negami (bez CD) |
| 4 | DC — libovolný distribuovaný element přes `DISTEL` |
| 5, 6 | jako 4, navíc Ngai coupling model (P(41)–P(43)) |

Např. `[P(35),P(40)] = [1,3]` znamená GE pro CSD a HN pro DSD; `[2,0]` je samotná
CD odezva CSD. Pro CD a HN **musí být P(10) resp. P(20) nastaveno na 6**.

**Tři podtypy GE** podle znamének U1 a U2:

| U1 | U2 | Chování |
|----|----|---------|
| < 0 | = 0 | asymetrické, jen GAM1; GAM1 = 1 ⇒ EDAE1 |
| > 0 | = 0 | symetrické, jen GAM1; GAM1 = 1 ⇒ EDAE2, GAM1 = 2 ⇒ gaussovská DRT |
| < 0 | > 0 | nejobecnější, různá U, φ i γ pro HF a LF oblast |

φ1, γ1, U1 ovlivňují vysokofrekvenční oblast, φ2, γ2, U2 nízkofrekvenční.

### Parametry

| P | Význam |
|---|--------|
| 1 | R<sub>∞</sub> — vysokofrekvenční limita odporu, přičítá se sériově k CSD |
| 2 | C<sub>∞</sub> / ε<sub>∞</sub> (epsilon, je-li ATEMP < 0) |
| 3, 4 | CSD: γ1, γ2 |
| 5–10 | CSD: R<sub>DAE</sub>(ρ₀), τ₀, U1, U2, φ1, φ2 |
| 10 | při NCH = 4: NELEM pro CSD (a P(9) je φ) |
| 11 | paralelní vodivost GE v elektrodové části (kapacita, je-li P(38) = 64) |
| 12 | paralelní vodivost přes objemovou část |
| 13, 14 | DSD: γ1, γ2 |
| 15–20 | DSD: ε<sub>DAE</sub>, τ, U1, U2, φ1, φ2 |
| 20 | při NCH = 4: NELEM pro DSD |
| 21–24 | elektrodový DCE3: R, T, U, φ |
| 25 | NDE3 — typ elektrodového DCE |
| 26 | R2 |
| 27 | C2 |
| 28 | R3 |
| 29 | C3 |
| 30 | sériový prvek, roli určuje P(37) (viz úvodní tabulka) |
| 33 | > 0: objemový podíl pro efektivní medium (Bruggeman) |
| 35 | volba CSD modelu (NCH) |
| 36, 38, 39 | transformační přepínače, viz `README.txt` |
| 40 | volba DSD modelu (NCH) |
| 41–43 | Ngai coupling: prahová ω = 1/P(41), pak odezva `P(42)/(1 + jω·P(43))` |
| 45 | NDE4 (při P(38) = 64) |

### Topologie

```
Z_objem     = ( P(1) - Z_CSD ) | Z_DSD | G12 | C_inf
Z_elektroda = ( R2 - ( (R3 - DCE3) | C3 ) ) | C2 | G11

Z = Z_objem - Z_elektroda - [sériový prvek P(30) dle P(37)]

G11, G12 = paralelní vodivosti P(11), P(12);  C_inf = P(2)
```

Při **P(33) > 0** se místo toho použije Bruggemanův vztah pro efektivní medium:
`ε_eff = 3·P(33)·(ε_DSD − ε_CSD)/(ε_DSD + 2ε_CSD − P(33)(ε_DSD − ε_CSD))`,
a výsledek se převede zpět na impedanci.

**MDE (= vstupní MODE):** při MDE ≥ 0 se počítá běžná CSD0 odezva, při MDE < 0
CSD1 (JRM/Moynihan). Znaménko MDE nemá vliv na DSD. Při **|MDE| = 8** se počítá
přechodová odezva místo frekvenční.

Speciální operace pro O-obvod (odečtení σ₀ přes P(36), škálování přes P(38)
a P(39)) jsou popsané v kořenovém `README.txt`.

---

## R Circuit (RSUB)

**Popis:** Diskrétní **DRT** — až **19 paralelních RC podobvodů v sérii**.
Rozšířená varianta K-obvodu.

**Parametry:**
| P | Prvek |
|---|-------|
| 1–28 | prvních 14 párů: P(1),P(2) · … · P(27),P(28) |
| 41–50 | dalších 5 párů (index KK > 28 se mapuje na P(KK+12)) |
| 29 | CP |
| 30 | role určena P(37): sériové R / paralelní G / sériové L / sériové C |
| 37 | volba role P(30) — viz úvodní tabulka |
| 38 | \|P(38)\| = 1 ⇒ všechna τ musí být fixní; jinak volba splinu |
| 39 | > 0: nastavuje minimální hodnotu τ |
| 40 | znaménko a velikost: viz níže |

**Interpretace páru podle |P(40)|** (`LV7.FOR:70-97`):

| \|P(40)\| | Význam | Člen |
|-------|--------|------|
| 0, 1 | (R, C), τ = R·C | `R/(1 + jωτ)` |
| 2, 3 | (R, τ) | `R/(1 + jωτ)` |
| 4 | (R, τ), škálováno | `τ·R/(1 + jωτ)` |
| 5 | přechodová odezva, „frekvence" je čas | `R·exp(−t/τ)` |

**P(40) < 0:** DRT reprezentuje komplexní dielektrickou konstantu (ATEMP < 0)
nebo komplexní kapacitu (ATEMP > 0). Jinak jde o kapacitu/impedanci.

Objeví-li se v páru τ = 0, výpočet se zkrátí (`GOTO 732`).

---

## S Circuit (SSUB)

**Popis:** Model elektrochemického superkondenzátoru — porézní elektroda jako
homogenní přenosové vedení.

**Topologie:**
```
Z_e = RE - LE                    (podélná impedance elektrody)
Z_1 = (RB - CF) | RSD | CDL      (příčná impedance póru)

Z = sqrt(Z_e · Z_1) · coth( sqrt(Z_e / Z_1) )
```

**Parametry:**
| P | Prvek |
|---|-------|
| 1 | RB — odpor faradaické větve v sérii s CF |
| 2 | RE — podélný odpor elektrody |
| 3 | RF — **načteno, ale nepoužito** |
| 4 | CF — faradaická kapacita |
| 5 | RSD — odpor samovybíjení |
| 6 | CDL — kapacita elektrické dvojvrstvy |
| 7 | LE — podélná indukčnost elektrody |

P(8)–P(30) se nepoužívají. Obvod nemá RP, CP ani L v obvyklém smyslu.

---

## T Circuit (TSUB)

**Neimplementováno.** `LV9.FOR:9` obsahuje pouze:

```fortran
PRINT *,'    ######  T CIRCUIT: CURRENTLY DUMMY!   ####'
PAUSE
STOP
```

Nastavení FUN = T tedy program ukončí. Adresář `FITTESTS/CKT.T` obsahuje jen
prázdnou šablonu `TTMP.T`.

---

## Co jsou DCE (Distributed Circuit Elements)?

**DCE = Distribuované obvodové prvky** — zobecnění ideálních R, C, L pro reálné systémy.

### Proč nestačí ideální R, C, L?

Ideální prvky předpokládají:
- Homogenní materiál
- Hladký povrch
- Jednu časovou konstantu

Reálné systémy mají:
- Drsné/porézní povrchy
- Nehomogenní strukturu
- Distribuci časových konstant
- Difúzní procesy

### Srovnání

| Vlastnost | Ideální C | DCE (např. CPE) |
|-----------|-----------|-----------------|
| Impedance | Z = 1/(jωC) | Z = 1/(T·(jω)^φ) |
| Exponent | −1 (fixní) | −φ (volný, 0 < φ ≤ 1) |
| Nyquist | Svislá přímka | Šikmá přímka |
| Fyzika | Deskový kondenzátor | Drsná elektroda |

### Kdy použít DCE?

| Systém | Typický DCE |
|--------|-------------|
| Koroze, pasivní vrstvy | CPE (NDE=2,3) |
| Difúze v elektrolytu | Warburg (NDE=9) |
| Porézní elektrody | CPE, Warburg |
| Dielektrická relaxace | Cole-Cole, H-N (NDE=6,7) |
| Polymery, skla | KWW (NDE=10, 32, 35–37) |
| Baterie, palivové články | Warburg, CPE |

### Fyzikální interpretace parametrů

**CPE jako příklad:**
```
Z_CPE = 1 / (Q · (jω)^φ)
```

| φ | Chování | Fyzikální význam |
|---|---------|------------------|
| 1.0 | Ideální C | Hladký povrch |
| 0.9 | Mírně neideální | Slabá drsnost |
| 0.8 | CPE | Porézní struktura |
| 0.5 | Warburg-like | Difúze |
| 0.0 | Ideální R | Čistě odporové |

### DCE vs složený obvod

**Důležité:** DCE **nelze** přesně nahradit konečným počtem R a C!

```
CPE ≠ R - C
CPE ≠ R | C
CPE ≈ nekonečný RC žebřík (teoreticky)
```

Proto LEVM implementuje DCE jako speciální matematické funkce, ne jako kombinace R/C.

---

## Parametry DCE

Každý DCE má 4 parametry + typ:

| Parametr | Obvyklý význam |
|----------|----------------|
| RDE | Odpor / škálování |
| TDE | Časová konstanta τ |
| UDE | Exponent nebo mez cutoffu U |
| PDE | Exponent φ (phi) |
| NDE | Typ elementu (0 = vypnuto, 1–37) |

> **Tento význam není univerzální.** Každý typ NELEM si čtyři vstupy interpretuje
> po svém. Např. u NELEM = 17 je T = C, U = R2, φ = L; u NELEM = 33 je R
> aktivační energie v eV a „frekvenční" sloupec obsahuje teplotu. Před použitím
> neobvyklého NELEM si vždy přečtěte příslušnou větev v `LV5.FOR`.

**Znaménko NDE:** větvení používá `IABS(NELEM)` (`LV5.FOR:174`), takže záporné
hodnoty vybírají stejný typ. Znaménko slouží jako příznak pro některé prvky.

**Vliv MODE (MDE):** u řady prvků (6, 29, 15, 16 a dalších) přepíná záporná
hodnota MODE výpočet na epsilon úroveň a mění význam vstupu R — např. u NELEM = 6
se R stane ε<sub>C∞</sub> (`LV5.FOR:274-278`), u NELEM = 29 při MDE = −1 se stane
ε<sub>C0</sub>.

### Přehled typů DCE (NDE)

Kompletní seznam podle větvení v `DISTEL` (`LV5.FOR:174-181`):

| NDE | Model | Poznámka |
|-----|-------|----------|
| 0 | Vypnuto | `DISTEL = 0` |
| 1 | R-C paralel | `Z = R/(1 + jωR·T)`, kde **T je τ**, ne kapacita |
| 2 | CPE model 1 | CPE paralelně s (R + C<sub>U</sub>); viz níže |
| 3 | CPE model 2 | `Z = R/(R·T·jω)^φ` |
| 4 | ZARC-Cole model 1 | `Z = R/(1 + R^U·T·(jω)^φ)` |
| 5 | ZARC-Cole model 2 | `Z = R/(1 + (jTω)^φ)` |
| 6 | Havriliak-Negami 1 | `Z = R·(1 + (jTω)^U)^(−φ)` |
| 7 | Havriliak-Negami 2 | **totožné s 6** — kód je `700 GOTO 600` |
| 8 | Havriliak-Negami 3 | `Z = R·sin(πU/2)·(1 + (jTω)^U)^(−φ)` |
| 9 | Generalized finite Warburg | Rovnice A.10–A.12 z práce JRM #252; chování řídí P8 a P18 |
| 10 | Williams-Watts (KWW) | Aproximace, čtyři úrovně přesnosti |
| 11 | Jonscher | Zobecněno o reálnou část |
| 12 | EDAE1 | Exponenciální DAE, obecný případ |
| 13 | EDAE2 | Exponenciální DAE, symetrický případ |
| 14 | GDAE | Gaussovská DAE, symetrický případ |
| 15 | Obecná homogenní difúze | Makroskopické parametry |
| 16 | Obecná homogenní difúze | Mikroskopické parametry |
| 17 | Ideální prvky R, C, R2, L | `(R\|C) \| (R2 - L)`; T = C, U = R2, φ = L |
| 18 | Dissado-Hill | Z úroveň |
| 19 | Dissado-Hill | epsilon úroveň |
| 20 | Power law #1 | `Y = R·ω^T + i·U·ω^φ` |
| 21 | Power law #2 | |
| 22 | Power law #3 | `Y = R·(iω)^T + U·i·ω^φ` |
| 23 | Ladder network #1 | Nekonečný žebřík se stejnými R a C |
| 24 | Ladder network #2 | |
| 25 | Efektivní medium s PCPE | |
| 26 | Power law #7 — **SCPE** | Elektrodové jevy, Z úroveň, sériové |
| 27 | Power law #8 — **PCPE** | Nearly constant loss, epsilon úroveň, paralelní |
| 28 | NCL / konstantní ztráta | Jen v OSUB, DSD část; U = 0, φ = 0; T = 0 ⇒ konstantní ztráta |
| 29 | Generalized Davidson-Cole | `Z = R(1+U)/(1 + U(1 + jTω)^φ)` |
| 30 | EDAE1 asymetrická, exaktní | Přesné log výrazy pro φ = 1 nebo 0; **φ musí být fixní** |
| 31 | GBEM EMA | Efektivní medium; φ musí být ±2, ±4 nebo ±8 |
| 32 | **KWW** s cutoffem | U = −ln(τs/τc); φ na vstupu je nepodstatné, ale musí být fixní |
| 33 | **Arrhenius teplotní fit** | **Není KWW.** ω = teplota v K, R = E v eV, `DATTYP` musí být R |
| 34 | **Generalized exponential DAE** | **Není KWW.** |
| 35 | KWW s cutoffem | U = −ln(τs/τc) |
| 36 | KWW s cutoffem | Varianta pro velmi malá β |
| 37 | KWW DRT pro libovolné β | `DATTYP = R`; frekvenční sloupec je τ, DRT je v reálném sloupci |

### Nejpoužívanější DCE podrobně

#### CPE (Constant Phase Element) - NDE = 2, 3

**NDE = 2** není čistý CPE. Kód (`LV5.FOR:200-216`) počítá:

```
Y = Y_R  +  T·(jω)^φ           kde  Y_R = jωU / (1 + R·jωU)
Z = 1 / Y
```

tedy CPE **paralelně** s větví R v sérii s kondenzátorem U. Zjednodušený tvar
`Z = 1/(T·(jω)^φ)` platí jen pro R = 0 a U = 0.

**NDE = 3:** `Z = R / (R·T·jω)^φ`

**Použití:**
- Neideální kapacitní chování elektrod
- Koroze a pasivní vrstvy
- Double-layer na drsném povrchu

**Nyquist:** Přímka pod úhlem φ·90° od reálné osy

**Efektivní kapacita:** `Ceff = (Q · R^(1−φ))^(1/φ)`

---

#### Warburg (difúze) - NDE = 9

Klasický tvar konečného Warburgu je `Z = R·tanh((jTω)^φ)/(jTω)^φ`, ale
implementace v LEVM je podstatně obecnější: realizuje rovnice A.10–A.12 z práce
JRM #252 včetně dvou rozšíření.

- **P18 = 4:** celková impedance zahrnuje dvě vážené rozhraní admitance
  paralelně — jednu s členem U²(jωτ)^φ, druhou s φ fixovaným na 1.
- **P18 = 8:** součet těchto dvou členů se použije pro výpočet impedance rozhraní.
- Váhování používá **P8 = AA** ve tvaru AM = (1 − AA) a AA.

Ostatní volby řízené P8 viz manuál str. 4-9 a 4-10.

**Varianty podle okrajových podmínek:**
- Nekonečná difúze: φ = 0.5, přímka 45°
- Konečná difúze (blokující): tanh → saturace
- Konečná difúze (transmisivní): coth varianta

**Použití:** difúze iontů v elektrolytu, lithiové baterie, palivové články,
elektrochromní materiály.

---

#### Havriliak-Negami - NDE = 6, 7, 8

**Vzorec (NDE = 6 i 7):** `Z = R · (1 + (jωτ)^U)^(−φ)`

Zde je **U = α** (šířka distribuce) a **φ = β** (asymetrie).

**NDE = 7 je v současném kódu totožné s NDE = 6** — větev 700 obsahuje pouze
`GOTO 600`. Volba mezi nimi tedy nemá žádný efekt.

**NDE = 8** přidává faktor `sin(πU/2)`: `Z = R·sin(πU/2)·(1 + (jωτ)^U)^(−φ)`.

**Speciální případy:**
| α (U) | β (φ) | Model |
|---|---|-------|
| 1 | 1 | Debye (ideální) |
| α | 1 | Cole-Cole |
| 1 | β | Davidson-Cole |
| α | β | Havriliak-Negami (obecný) |

**Použití:** dielektrická relaxace, polymery, skla, biologické tkáně.

Při MODE < 0 se výpočet přepne na epsilon úroveň a R znamená ε<sub>C∞</sub>.

---

#### KWW (Kohlrausch-Williams-Watts) - NDE = 10, 32, 35, 36, 37

**Časová doména:** `φ(t) = exp(−(t/τ)^β)`

| NDE | Varianta |
|-----|----------|
| 10 | aproximativní výpočet, čtyři úrovně přesnosti |
| 32 | s cutoffem, U = −ln(τs/τc); φ na vstupu nepodstatné, musí být fixní |
| 35 | s cutoffem, U = −ln(τs/τc) |
| 36 | s cutoffem, vhodné i pro velmi malá β |
| 37 | přímý výpočet DRT pro libovolné β; `DATTYP = R`, frekvenční sloupec obsahuje τ |

Pro NDE = 32, 35, 36 platí: je-li U ≤ −40, cutoff nemá v žádném praktickém
frekvenčním rozsahu efekt.

**Použití:** stretched exponential relaxace, amorfní materiály, iontová vodivost
ve sklech.

**Poznámka:** LEVM počítá KWW přes DRT (distribuci relaxačních časů).

**NDE = 33 a 34 nejsou KWW** — 33 je Arrheniovský teplotní fit, 34 je
generalized exponential DAE.

---

## Poznámky k použití

### 1. Vypnutí prvku

Nastavení parametru na 0 **neznamená vždy odstranění větve**:

- **DCE** se spolehlivě vypne jen nastavením **NDE = 0**.
- V `SDEA` (`LV4.FOR:719`) platí: `CA = 0 ⇒ ZA = RA` (čistý odpor),
  `RA = 0 a CA ≠ 0 ⇒ ZA = 1/(jωCA)` (čistý kondenzátor). Nulou tedy prvek
  neodstraníte, jen přepnete typ.
- Stejná logika je i v jednotlivých obvodech: např. v `ASUB` znamená
  `R1 = 0` při `C1 ≠ 0` čistý kondenzátor, ne zkrat.
- Pro odstranění vlivu parametru z **dat** slouží NFREE = 3, viz níže.

### 2. Výběr obvodu

Parametr **FUN** v INFL — pouze A–K, O, R, S, T (T je stub).

### 3. Volné parametry — NFREE

| NFREE | Význam |
|-------|--------|
| 0 | Fixní |
| **1** | **Volný, omezený na kladné hodnoty** — při záporné hodnotě se aktivuje exponenciální penalizační funkce (`LV2.FOR:881`, `LV2.FOR:1004`) |
| **2** | **Volný bez omezení — smí jít do záporu** |
| 3 | Volný, a po fitu se jeho vliv **odečte z dat** do souboru OUTIN (`LV0.FOR:931`) |

> **Toto je nejčastější zdroj problémů.** Pokud fit nekonverguje, hlásí
> singulární matici, nebo se odhad parametru zmenšuje k nule, je pravděpodobné,
> že parametr potřebuje jít do záporu — přepněte jeho NFREE z 1 na **2**.
> Manuál: *„If any of your free fit parameters either are or might need to be
> negative … set the NFREE index number for that parameter to 2 rather than 1."*

Parametry **P(33)–P(40) jsou vždy fixní** — výpis je uvádí pod hlavičkou
`THE FOLLOWING PARAMETERS ARE ALWAYS FIXED` (`LV0.FOR:425`).

Pro IRCH < 2 musí být U (P(31)) i XI (P(32)) fixní, jinak program vypíše
`BAD CHOICE OF FREE/FIX FOR U AND XI` a vynutí je na 0 (`LV0.FOR:250`).

### 4. Typ dat

**DATTYP** = `C` komplexní, `R` pouze reálná část, `I` pouze imaginární.

Obvody s povinným omezením:
- **G:** pouze R nebo I (při C se program ukončí)
- **I:** pouze R
- **NDE = 33 a 37:** pouze R
