# LEVM Circuit Topologies

Přehled podporovaných obvodových modelů v programu LEVM v moderní textové notaci.

## Legenda

| Symbol | Význam |
|--------|--------|
| `+` | Sériové spojení |
| `\|\|` | Paralelní spojení |
| `DCEn` | Distribuovaný element (slot n), typ určen parametrem NDE |
| `(RA,CA)` | Augmentovaný DCE - DCE paralelně s RC článkem |
| `DAE` | Distribution of Activation Energies |

## Společné prvky všech obvodů

Každý obvod má na vnější úrovni:

```
L + ( vnitřní_struktura ) || (RP||CP)
```

| Parametr | Popis | Vypnutí |
|----------|-------|---------|
| P(28) = RP | Paralelní odpor (vnější) | RP = 0 |
| P(29) = CP | Paralelní kapacita (vnější) | CP = 0 |
| P(30) = L | Sériová indukčnost | L = 0 |

---

## A Circuit (ASUB)

**Popis:** 6 subcircuits v sérii, paralelně s RP/CP, v sérii s L.

**Topologie:**
```
L + (
    (R1||C1) + (R2||C2) + (R3||C3) +
    DCE1 +
    (RA||CA||DCE2) +
    ((R5 + DCE3 + C5) || (R4||C4))
) || (RP||CP)
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

---

## B Circuit (BSUB)

**Popis:** Vnořené RC s distribuovanými elementy.

**Topologie:**
```
L + (
    R1 + (RA||CA||DCE1) + (
        (R2 + DCE2 + ((R3 + DCE3) || C3)) || (C2 + DCE4)
    )
) || (RP||CP)
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

---

## C Circuit (CSUB)

**Popis:** Dva augmentované DCE + jeden regulární DCE.

**Topologie:**
```
L + (
    R1 +
    (((R3||C3) + (RA1||CA1||DCE1)) || (R2 + C2)) +
    ((R4 + DCE3 + (R5||C5)) || (RA2||CA2||DCE2))
) || (RP||CP)
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
L + (
    R1 +
    ((RA||CA) || DAE) +
    ((R2 + (R3||C2)) || (((R4 + DCE4)||C4) + ((R5 + DCE5)||C5)))
) || (RP||CP)
```

**DAE typy:**
- P(15) = 0: Exponenciální DAE (EDAE)
- P(15) ≠ 0: Gaussovské DAE (GDAE)

**Parametry:**
| P | Prvek | P | Prvek |
|---|-------|---|-------|
| 1 | R1 | 16 | RDE4 |
| 2 | RA | 17 | TDE4 |
| 3 | CA | 18 | UDE4 |
| 4 | RDAE | 19 | PDE4 |
| 5 | TDAE | 20 | NDE4 |
| 6 | U1 | 21 | RDE5 |
| 7 | U2 | 22 | TDE5 |
| 8 | φ1 | 23 | UDE5 |
| 9 | φ2 | 24 | PDE5 |
| 10 | R2 | 25 | NDE5 |
| 11 | C2 | 26 | R5 |
| 12 | R3 | 27 | C5 |
| 13 | R4 | 28 | RP |
| 14 | C4 | 29 | CP |
| 15 | DCH | 30 | L |

---

## E Circuit (ESUB)

**Popis:** 5 DCE slotů - maximální flexibilita pro distribuované elementy.

**Topologie:**
```
L + (
    DCE1 + R1 +
    (DCE2 || ((R3 + DCE3) || DCE4)) +
    R2 + DCE5
) || CP
```

**Poznámka:** Tento obvod nemá RP, pouze CP.

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

## Co jsou DCE (Distributed Circuit Elements)?

**DCE = Distribuované obvodové prvky** - zobecnění ideálních R, C, L pro reálné systémy.

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

| Vlastnost | Ideální C | DCE (např. Q) |
|-----------|-----------|-----------------|
| Impedance | Z = 1/(jωC) | Z = 1/(T·(jω)^φ) |
| Exponent | -1 (fixní) | -φ (volný, 0 < φ ≤ 1) |
| Nyquist | Svislá přímka | Šikmá přímka |
| Fyzika | Deskový kondenzátor | Drsná elektroda |

### Kdy použít DCE?

| Systém | Typický DCE |
|--------|-------------|
| Koroze, pasivní vrstvy | Q (NDE=2,3) |
| Difúze v elektrolytu | Warburg (NDE=9) |
| Porézní elektrody | Q, Warburg |
| Dielektrická relaxace | Cole-Cole, H-N (NDE=6,7) |
| Polymery, skla | KWW (NDE=10, 32-37) |
| Baterie, palivové články | Warburg, Q |

### Fyzikální interpretace parametrů

**Q jako příklad:**
```
Z_Q = 1 / (Q · (jω)^φ)
```

| φ | Chování | Fyzikální význam |
|---|---------|------------------|
| 1.0 | Ideální C | Hladký povrch |
| 0.9 | Mírně neideální | Slabá drsnost |
| 0.8 | Q | Porézní struktura |
| 0.5 | Warburg-like | Difúze |
| 0.0 | Ideální R | Čistě odporové |

### DCE vs složený obvod

**Důležité:** DCE **nelze** přesně nahradit konečným počtem R a C!

```
Q ≠ R + C
Q ≠ R || C
Q ≈ nekonečný RC žebřík (teoreticky)
```

Proto LEVM implementuje DCE jako speciální matematické funkce, ne jako kombinace R/C.

---

## Parametry DCE

Každý DCE má 4 parametry + typ:

| Parametr | Význam |
|----------|--------|
| RDE | Odpor / škálování |
| TDE | Časová konstanta τ |
| UDE | Exponent U |
| PDE | Exponent φ (phi) |
| NDE | Typ elementu (1-37) |

### Přehled typů DCE (NDE)

| NDE | Model | Vzorec |
|-----|-------|--------|
| 0 | Vypnuto | Z = 0 |
| 1 | R-C paralel | Z = R/(1 + jωRC) |
| 2 | Q model 1 | Z = 1/(T·(jω)^φ) |
| 3 | Q model 2 | Z = 1/(jTω)^φ |
| 4 | Z-C model 1 | Z = R/(1 + R^U·T·(jω)^φ) |
| 5 | Z-C model 2 | Z = R/(1 + (jTω)^φ) |
| 6 | Havriliak-Negami 1 | Z = R/(1 + (jTω)^U)^φ |
| 7 | Havriliak-Negami 2 | Z = R/(1 + (jRTω)^U)^φ |
| 8 | Havriliak-Negami 3 | Z = R·sin(χ)/(1 + (jRTω)^U)^φ |
| 9 | Generalized Warburg | Z = R·tanh((jTω)^φ)/(jTω)^φ |
| 10 | Williams-Watts | KWW response |
| 11 | Jonscher | Generalized with real part |
| 12 | EDAE1 | Exponential DAE, asymmetric |
| 13 | EDAE2 | Exponential DAE, symmetric |
| 14 | GDAE | Gaussian DAE |
| 15 | Diffusion (macro) | General diffusion DCE |
| 16 | Diffusion (micro) | General diffusion DCE |
| 17 | RLC paralel | C \|\| R \|\| (R + L) |
| 18 | Dissado-Hill Z | Z level |
| 19 | Dissado-Hill ε | Epsilon level |
| 20-28 | PLE, Ladder, etc. | Various specialized models |
| 29 | Modified Davidson-Cole | Z = R(1+U)/(1+U(1+jTω)^φ) |
| 30 | Exact EDAE | For φ = 0 or 1 |
| 31 | GBEM EMA | Effective medium |
| 32-37 | KWW variants | Various β values |

### Nejpoužívanější DCE podrobně

#### Q (Constant Phase Element) - NDE = 2, 3

**Vzorec:** `Z = 1 / (Q · (jω)^φ)`

**Použití:**
- Neideální kapacitní chování elektrod
- Koroze a pasivní vrstvy
- Double-layer na drsném povrchu

**Nyquist:** Přímka pod úhlem φ·90° od reálné osy

**Efektivní kapacita:** `Ceff = (Q · R^(1-φ))^(1/φ)`

---

#### Warburg (difúze) - NDE = 9

**Vzorec:** `Z = R · tanh((jTω)^φ) / (jTω)^φ`

**Varianty podle okrajových podmínek:**
- Nekonečná difúze: φ = 0.5, přímka 45°
- Konečná difúze (blokující): tanh → saturace
- Konečná difúze (transmisivní): coth varianta

**Použití:**
- Difúze iontů v elektrolytu
- Lithiové baterie
- Palivové články
- Elektrochromní materiály

---

#### Havriliak-Negami (Cole-Cole zobecnění) - NDE = 6, 7

**Vzorec:** `Z = R / (1 + (jωτ)^α)^β`

**Speciální případy:**
| α | β | Model |
|---|---|-------|
| 1 | 1 | Debye (ideální) |
| α | 1 | Cole-Cole |
| 1 | β | Davidson-Cole |
| α | β | Havriliak-Negami (obecný) |

**Použití:**
- Dielektrická relaxace
- Polymery
- Skla
- Biologické tkáně

---

#### KWW (Kohlrausch-Williams-Watts) - NDE = 10, 32-37

**Časová doména:** `φ(t) = exp(-(t/τ)^β)`

**Použití:**
- Stretched exponential relaxace
- Amorfní materiály
- Iontová vodivost ve sklech

**Poznámka:** LEVM počítá KWW přes DRT (distribuci relaxačních časů)

---

## Poznámky k použití

1. **Vypnutí prvku:** Nastavit příslušný parametr = 0
2. **Výběr obvodu:** Parametr FUN v INFL (A, B, C, D, E, ...)
3. **Volné parametry:** NFREE array (0 = fixed, 1 = free, 2 = free + non-negative)
4. **Komplexní fit:** DATTYP = 'C', pouze reálná: 'R', pouze imaginární: 'I'
