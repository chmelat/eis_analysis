# Testování v projektu LEVM

Jak je v tomto repozitáři řešené ověřování správnosti kódu. Sestaveno ze
zdrojových souborů (`LV0`–`LV9.FOR`), obsahu `FITTESTS/`, dávkových souborů
v `LEVMMISC/` a manuálu; empirické údaje pocházejí z běhu popsaného v části 8.

---

## 1. Shrnutí

**V projektu není žádný testovací framework, žádné očekávané výstupy k porovnání
a žádné automatické vyhodnocení pass/fail.** Není tu build systém ani CI.

Místo toho existuje testovací přístup postavený na čtyřech pilířích:

| Pilíř | Kde |
|-------|-----|
| **Korpus referenčních vstupů** — 145 spustitelných fitů napříč 15 obvody | `FITTESTS/` |
| **Dávkové ovladače** pro opakované a sériové běhy | `LEVMMISC/LEVMISC/*.BAT`, `LEVMMISC/MKLIO/*.BAT` |
| **Diagnostika zabudovaná v programu** — kódy ukončení, detekce singulární matice, metriky kvality fitu | `LV1.FOR`, `LV2.FOR` |
| **Ověřovací postupy popsané v manuálu** — round-trip simulace, test idempotence, perturbace odhadů | `LEVMMANUAL.pdf`, `README.txt` |

Klíčová vlastnost, na které celý přístup stojí: **LEVM umí sám generovat přesná
data ze zadaného modelu** (`MAXFEV = 0`). Program je tak sám sobě generátorem
testovacích dat — může si vyrobit data z parametrů, které zná, a pak je zkusit
nafitovat zpět. Vyhodnocení shody ale vždy dělá člověk očima nad výpisem.

---

## 2. Korpus FITTESTS

169 souborů v 15 adresářích `CKT.<písmeno>`, jeden na fitovací obvod.

| Kategorie | Počet | Poznámka |
|-----------|------:|----------|
| Spustitelné vstupy (hlavička + data) | 145 | |
| Šablony (jen hlavička, bez dat) | 17 | `*TMP.T`, `AZC.T`, `HTopInfile`, `GPNPATOP` |
| Nevstupní soubory | 7 | surová data z přístrojů, seznamy frekvencí |

### Rozložení podle obvodů

| Obvod | Spustitelných | Šablon | | Obvod | Spustitelných | Šablon |
|-------|--------------:|-------:|-|-------|--------------:|-------:|
| A | 5 | 2 | | J | 1 | 1 |
| B | 7 | 1 | | K | 8 | 1 |
| C | 1 | 1 | | O | 61 | 2 |
| D | 3 | 1 | | R | 19 | 1 |
| E | 5 | 1 | | **S** | **0** | 1 |
| F | 1 | 1 | | **T** | **0** | 1 |
| G | 12 | 1 | | | | |
| H | 21 | 1 | | | | |

Pokrytí je silně nerovnoměrné. O-obvod (61 testů) a H-obvod (21) mají řádově
víc než C, F, I a J (po jednom). **S-obvod a T-obvod nemají ani jeden
spustitelný test** — u T je to konzistentní, protože `TSUB` je stub, který
program ukončí, ale S-obvod (model superkondenzátoru) je plně implementovaný
a neotestovaný.

### Konvence pojmenování

| Vzor | Význam |
|------|--------|
| `<X>TMP.T` | prázdná šablona hlavičky pro obvod X — všechny parametry nulové |
| `<X>.T` | hlavička konkrétního testu bez datové části (např. `AZC.T`) |
| `<X>TST`, `<X>TST.*` | kanonický test obvodu X |
| `*DAT.B`, `*.DAT` | jen datová část, ke složení se šablonou |
| `*.fra`, `*.U1x` | surový výstup z přístroje (Solartron FRA) |
| `LVFREQS` | seznam frekvencí pro generování simulací |

### Integritní vlastnost korpusu

Hodnota **MD** na řádku 3 (počet datových bodů) **souhlasí se skutečným počtem
datových řádků u všech 145 spustitelných souborů** — včetně 12 souborů, kde je
MD záporné (volba přesnějšího ε<sub>V</sub>, viz `README.txt`). To je jediná
strojově ověřitelná invarianta, kterou korpus má, a drží.

Datová část toleruje **prázdné pořadové číslo** v prvních pěti sloupcích
(formát `I5` přečte prázdno jako 0 a hodnota se nepoužívá). Např. `AZC` čísluje,
`BTST` ne — obojí je platné.

---

## 3. Čtyři druhy testovacích souborů

**a) Šablony.** Prázdná hlavička s vybraným obvodem, ze které se skládá nový
vstup. Nejsou to testy.

**b) Syntetická data se známým výsledkem.** Titulek na prvním řádku nese exaktní
hodnoty, kterými byla data vygenerována:

```
 AZC A CKT. EXACT PAR. VALUES: ZC(2D6,1D0,0.3),1D+6,1D-12
   BTST AZC B CKT ZC(2D6,1D0,0.3),1D-12,1D+6
KWW EXACT FREQ RESP. 1, 1E-4, 0.5
```

**Titulek je jediné místo, kde je uložena očekávaná hodnota.** Je to volný text
pro člověka, ne strojově čitelný formát. Takto anotovaných souborů je málo —
`AZC`, `AZC.T`, `BTST`, `BTS0Y03`, `CTST`, `ETST`, `FTST`,
`RTSTWC0.9`, `RTSTWC1.9`, `RTSTWD.19`, `RTSTT.9`, `OTSTHN`,
`RANDLESexactParamValues`, `DCDCJPNPAEXACTZN` a `AZCE`/`AZC1` v `LEVMMISC/LEVMISC`.

**c) Reálná měřená data.** Většina korpusu. Titulek identifikuje vzorek,
teplotu a datum měření, ne očekávaný výsledek:

```
 NOW NACL 23 degree C EPSY DATA 10/25/92
  ZIRCONIA 500C Z DATA: Boukamp/van Hassel 8/23/89
MOYNIHAN LAS 24 c EPSILON DATA 3/3/94. HN EXACT DATA: NO CUTOFF
```

U těchto souborů neexistuje objektivní kritérium správnosti. Slouží jako
regresní vzorky — očekává se, že fit doběhne a dá „rozumné" hodnoty.

**d) Kontrastní a křížové dvojice.** Nejzajímavější část korpusu:

| Dvojice | Co ověřuje |
|---------|-----------|
| `CKT.H/goodtRANDHfit` vs `CKT.H/poorRANDHfiterrdat` | jak vypadá dobrý a špatný fit stejných dat |
| `CKT.K/KTSTKWW1.11` vs `CKT.R/RTSTKWW1.19` | **křížová kontrola mezi obvody** — táž KWW data fitovaná K-obvodem s 11 páry a R-obvodem s 19; titulek KTSTKWW1.11 to říká výslovně: `SAME INPUT DATA AS THAT IN RTSTKWW1.19` |
| `CKT.B/RANDLESexactParamValues` vs `RANDLESwithRanErr` | týž model s exaktními daty a s přidaným náhodným šumem |
| `CKT.A/AZC` vs `LEVMMISC/LEVMISC/AZCE` | táž úloha s daty oříznutými na 2 desetinná místa a s plnou přesností |
| `CKT.O/OTSTNOW1.422` … `OTSTNOW4.422`, `OTSTNN1`…`OTSTNN4` | varianty jedné úlohy s postupně měněnou volbou modelu |

---

## 4. Zabudovaný mechanismus relativních chyb (IPAR < 0)

Program má jednu vestavěnou obdobu asserce. Je-li na řádku 2 vstupního souboru
**IPAR < 0** (sloupce 25–34), nastaví se `IORIG = 1` (`LV0.FOR:99-104`) a vstupní
hodnoty parametrů se uloží stranou jako `PEX` (`LV0.FOR:233-238`). Po fitu se pro
každý volný parametr spočítá (`LV2.FOR:641-650`):

```
PRELE(j) = ( fitovaná_hodnota − PEX(j) ) / PEX(j)
```

a do výstupu přibude čtvrtý řádek `RELATIVE ERRORS OF PARAMETERS` plus dvě
souhrnná čísla:

```
PRSDAV=  2.3702D-01      PRSDRMS=  3.3926D-01
```

`PRSDAV` je průměr absolutních relativních odchylek, `PRSDRMS` jejich RMS.

> **Přesně vzato to není oracle proti pravdě, ale míra posunu od vstupních
> hodnot.** Komentář ve zdroji říká `SET FOR "EXACT" INPUT PARAMETER REL SD
> CALCULATIONS` — záměr je vložit do vstupu *exaktní* hodnoty a nechat si
> spočítat, jak přesně je fit trefil. Ve všech čtyřech souborech, které
> to používají, jsou ale vstupní hodnoty záměrně rozhozené počáteční odhady,
> takže `PRSDAV` měří vzdálenost od startu, ne od pravdy.

Volbu používají jen čtyři soubory v celém repozitáři:
`FITTESTS/CKT.B/BTST`, `FITTESTS/CKT.B/BTS0Y03`,
`LEVMMISC/LEVMISC/AZCE`, `LEVMMISC/LEVMISC/AZC1`.

---

## 5. Testovací postupy

### 5.1 Round-trip: simulace → fit zpět

Nejsilnější postup, jaký projekt má. `MAXFEV = 0` vypne iterování a program
jen spočítá predikce modelu pro zadané parametry — vznikne exaktní datová sada.
Ta se pak fituje zpět a porovná se, zda vyjdou původní parametry.

Řídí to `LEVMMISC/LEVMISC/SIM.BAT`:

```bat
del SIM2
CALL MKFREQ.EXE            REM interaktivně vygeneruje seznam frekvencí -> LVFREQS
CALL MKINP SIM1 LVFREQS SIM2   REM hlavička SIM1 + frekvence -> vstup SIM2
CALL RNL SIM2              REM MAXFEV=0 -> simulace, výsledek do OUTIN
CALL RNL OUTIN             REM fit simulovaných dat
CTD                        REM grafické porovnání
```

Šablona `SIM1` má `MAXFEV = 0` a `IRE = -10`, aby se OUTIN automaticky vyrobil.
Manuál (str. 2-7) požaduje, aby alespoň jedna hodnota NFREE byla nenulová.

Obdoba pro sérii modelů je `LEVMMISC/MKLIO/MSIM.BAT`, která přes `LSIM.BAT`
projede sedm sad parametrů `SER2_O(1-7).INP` a vyrobí `DATA2-Z(1-7)`.

### 5.2 Test idempotence — dosažení lokálního minima

Manuál (str. 1-35) to označuje za základní test konvergence:

> *A useful test of reaching at least a local minimum is to always run with
> IRE < -9, so the file OUTIN is automatically produced. Then, when a run is
> finished, just do RNL OUTIN (or RNLO) and see if the same results as before
> are obtained.*

Automatizuje to `LEVMMISC/LEVMISC/RNLVOI.BAT`, který spustí fit třikrát za sebou
a tři výpisy slepí do jednoho souboru k ručnímu porovnání:

```bat
DEL OUTIN
CALL RNL %1
IF NOT EXIST OUTIN GOTO END
COPY PNTOUTL PNTLV1
CALL RNL OUTIN
COPY PNTOUTL PNTLV2
CALL RNL OUTIN
COPY PNTOUTL PNTLV3
COPY PNTLV1 + PNTLV2 + PNTLV3 PNTOUTL
```

Proto má naprostá většina souborů v `FITTESTS` nastavené `IRE = -11`.

### 5.3 Perturbace počátečních odhadů — test globálního minima

Manuál, str. 1-35:

> *To increase the chances that an absolute minimum has been found, change the
> input free parameter choices to values appreciably different from those
> originally used. Then do a fit and see if the same set of parameter estimates
> is obtained. If so, the chances are that the proper minimum has been found.
> If not, use the result with the lowest value of the standard deviation, SD.*

V korpusu je tento postup vidět na souborech jako `AZCE`, kde jsou počáteční
odhady záměrně posunuté o ~1 % (a P(6) dokonce o −45 %) proti exaktním hodnotám.

### 5.4 Sériový fit se statistikou

`LEVMMISC/MKLIO/MULTFIT.BAT` projede sedm datových sad týmž modelem a
`EXTRLVM.EXE` z každého běhu vytáhne vybrané parametry:

```bat
copy azc40 outin.
for %%f in (1 2 3 4 5 6 7) do call onefit2.bat data-z%%f.dat
```
```bat
REM onefit2.bat
MKLIO.EXE -aser1_i %1 %2 %3
levm
extrlvm -s1,6,7,9,29, -g par1.dat -aser1_o
```

Přepínač `-s1,6,7,9,29,` vybírá parametry P(1), P(6), P(7), P(9), P(29). Výstup
`PAR1.SUM` obsahuje podle `LEVMMISC/MKLIO/INDEX.TXT` tři řádky:

```
1. řádek: střední hodnoty parametrů
2. řádek: průměrné rms odchylky d = (x − x̄)
3. řádek: absolutní průměrné odchylky d/x̄
```

**To je nejblíž strojovému vyhodnocení, co v projektu je** — reprodukovatelnost
odhadu parametru napříč sadou měření vyjádřená čísly. Vyhodnocení, zda je
rozptyl přijatelný, ale opět dělá člověk.

### 5.5 Diagnostika při nefunkčním vstupu

Manuál, str. 1-36: pokud program vůbec nezačne iterovat, nastav `MAXFEV = 0`
a prohlédni predikce v `AUXPNTL` nebo `OUTIN`. Tím se otestuje **jen model**,
bez optimalizace. Když ani to nefunguje, je chyba ve volbách vstupu; když to
funguje, porovná se výstup s daty a hledají se oblasti špatné shody.

Doporučená druhá pomůcka: porovnat vstupní soubor řádek po řádku s `AZC` nebo
s jiným souborem z `FITTESTS`.

---

## 6. Ovladače a nástroje

### Dávkové soubory

| Soubor | Umístění | Co dělá |
|--------|----------|---------|
| `RNL.BAT` | `LEVMMISC/TESTRUN` | základní spouštěč: smaže staré výstupy, zkopíruje vstup na `INFL`, spustí `LEVM` |
| `RNLVOI.BAT` | `LEVMMISC/LEVMISC` | trojitý běh pro test idempotence (5.2) |
| `RPL.BAT` | `LEVMMISC/LEVMISC` | běh s měřením času, výstup rovnou na tiskárnu |
| `SIM.BAT` | `LEVMMISC/LEVMISC` | round-trip simulace (5.1) |
| `ONEMKE.BAT` | `LEVMMISC/MKLIO` | složí hlavičku + data do vstupního souboru |
| `CONLFL.BAT` | `LEVMMISC/MKLIO` | složí hlavičku z jednoho a data z druhého úplného souboru |
| `ONEFIT.BAT` / `ONEFIT2.BAT` | `LEVMMISC/MKLIO` | jeden fit + extrakce parametrů |
| `MULTFIT.BAT` | `LEVMMISC/MKLIO` | sériový fit se statistikou (5.4) |
| `LSIM.BAT` / `MSIM.BAT` | `LEVMMISC/MKLIO` | sériová simulace |

Klíčový detail `RNL.BAT`: **před každým během maže všechny výstupní soubory**
(`PNTOUTL`, `AUXPNTL`, `LINOUT`, `LVOUT`, …). LEVM je nezkracuje spolehlivě sám,
takže bez toho hrozí čtení zbytků z předchozího běhu. Kdo spouští LEVM ručně,
musí to udělat také.

### Nástroje

| Nástroj | Zdroj | Role v testování |
|---------|-------|------------------|
| `MKLIO.EXE` | **jen binárka** | skládá vstupní soubory z hlavičky a dat |
| `EXTRLVM.EXE` | **jen binárka** | vytahuje vybrané parametry z výsledku — nejblíž „aserci" |
| `MKFREQ.EXE` | `MKFREQ.FOR` | generuje seznam frekvencí pro simulaci |
| `MKLTP.EXE` | `MKLTP.FOR` | příprava teplotních dat |
| `DINPUT.EXE` | `DINPUT.FOR` | transformace dat, včetně **přidání náhodného šumu** do exaktních simulovaných dat |
| `GINPUT.EXE` | `GINPUT.FOR` | další transformace vstupů |
| `CTD.EXE` | **jen binárka** | grafické porovnání dat a fitu (DOS) |
| `LEVMVIEW.EXE` | **jen binárka** | totéž ve Windows verzi |

> Dvě nejdůležitější věci pro automatizaci — `MKLIO` (skládání vstupů) a
> `EXTRLVM` (extrakce výsledků) — jsou 16bitové DOS binárky **bez zdrojového
> kódu**. Na dnešním systému je nelze spustit ani přeložit. Kdo chce testování
> zautomatizovat, musí obě nahradit.

---

## 7. Kontroly zabudované v programu

### Kód ukončení optimalizátoru (INFO)

`LMDER` vrací `INFO`, které `SPFIT` přeloží na hlášku (`LV1.FOR:225-250`,
formáty 761–768):

| INFO | Hláška | Význam |
|-----:|--------|--------|
| 0 | `IMPROPER INPUT PARAMETERS` | špatný vstup — běh se ukončí přes `STOPX` |
| 1 | `...RELATIVE ERROR IN THE SUM OF SQUARES IS AT MOST FTOL` | konvergence podle FTOL |
| 2 | `...ERROR BETWEEN THE VECTOR OF FREE PARAMETERS AND THE DESIRED SOLUTION IS AT MOST XTOL` | konvergence podle XTOL |
| 3 | obojí | nejlepší případ |
| 4, 8 | `...FVEC IS ORTHOGONAL TO THE COLUMNS OF THE JACOBIAN` | text sám dodává: `THIS IS NOT A GOOD INDICATOR OF A CONVERGED MINIMUM` |
| 5 | `TERMINATION: NUMBER OF CALLS TO FCN ... EXCEEDS MAXFEV` | došly iterace; `LV0.FOR:1194` navíc hlásí `!!! INCOMPLETE CONVERGENCE: MAXFEV TOO SMALL !!!` |
| 6 | `TERMINATION: FTOL IS TOO SMALL` | dál už to nejde zlepšit |
| 7 | `TERMINATION: XTOL IS TOO SMALL` | dtto pro parametry |

Záporné `INFO` znamená chybu z `FCN`: **−15** = XI < −4, **−16** = XI > 4
(`LV2.FOR:43-51`, kontroluje se jen když je XI volný parametr).

### Tolerance

Nastavují se v `LV0.FOR:337-345` z parametru IOP (řádek 2):

| IOP | FTOL | XTOL |
|-----|------|------|
| 0 | 10⁻³⁰ | 10⁻⁴⁸ |
| < 0 | 10^(−\|IOP\|) | totéž |

`GTOL` je vždy 0. `IGACC` řídí přesnost kvadratur: `EPSG = 10^(−|IGACC|)`
(`LV0.FOR:649`). Manuál doporučuje IOP = −4 nebo −5 pro experimentální data,
−8 až −9 pro data s velmi malými chybami.

### Kontroly numerické korektnosti

| Kontrola | Kde | Hláška |
|----------|-----|--------|
| Nulový řádek v LU rozkladu | `LUDCMP`, `LV2.FOR:768-772` | `*** SINGULAR MATRIX ***`, nastaví `ISD = 1` |
| ISD propagované do výpisu | `SVDM`, `LV2.FOR:607-610` | `!!!!!** SINGULAR MATRIX - BEWARE **!!!!!` |
| Nekladný součin diagonál kovarianční matice | `SVDM`, `LV2.FOR:730-734` | počítá `IBAD`; při IBAD ≥ 1 vypíše `#!!!!! BAD CORRELATION MATRIX !!!!!#` a `PROBABLY AT LEAST ONE FREE PARAMETER IS NOT USED IN FITTING FUNCTION` |
| Velká relativní směrodatná odchylka | `SVDM`, `LV2.FOR:692-700` | **viz níže** |
| Divergence modelu | `FCN`, `LV2.FOR:102-107` | hlídá \|FN\| > 10⁷² |

> **Kontrola velkých relativních odchylek je od roku 1994 fakticky vypnutá.**
> Manuál (str. 1-36) ji popisuje jako standardní varování
> `NOTE: LARGE VALUE OF ONE OR MORE RELATIVE STANDARD DEVIATIONS`.
> Podmínka ve zdroji je ale `IF(PDAV.GT.1.D40)`, s komentářem
> `temporarily changed 1.d3 to 1.d40  8/17/94 !!!!!!!!!!!!!`.
> Práh 10⁴⁰ nemůže v praxi nastat, takže varování nikdy nevyskočí.
> V mém běhu celého korpusu (část 8) se neobjevilo ani jednou.

### Metriky kvality fitu ve výpisu

Každý běh zapisuje do `PNTOUTL`:

| Veličina | Význam |
|----------|--------|
| `PARAMETER ESTIMATES` | odhady volných parametrů |
| `ESTIMATED STD DEV OF PARAMETERS` | jejich směrodatné odchylky |
| `ESTIMATED RELATIVE STD DEV OF PARAMETERS` | relativní odchylky — hlavní ukazatel spolehlivosti |
| `PDAV`, `PDRMS` | průměr a RMS relativních směrodatných odchylek |
| `ESTIMATED PARAMETER CORRELATION MATRIX` | korelace mezi parametry; hodnoty blízké ±1 znamenají nerozlišitelné parametry |
| `SD WR`, `SD WI`, `SD WC` | vážené směrodatné odchylky reálné, imaginární a komplexní části |
| `SD RR`, `SD RI`, `SD RC` | totéž nevážené |
| `SSE` | součet čtverců |
| `NDF`, `FQF` | stupně volnosti a Akaikeho kritérium `FQF = K·ln(SUMSQ) + 2·NFREI` (`LV1.FOR:368`) |
| `PRSDAV`, `PRSDRMS` | jen při IPAR < 0, viz část 4 |

Manuál doporučuje jako primární kritérium porovnání modelů **SD RC** a **FQF**,
a jako varovný signál relativní směrodatnou odchylku parametru řádu jednotek
a výš.

---

## 8. Empirické ověření korpusu

Celý korpus jsem spustil na Linuxu s vlastním překladem (podrobnosti v kořenovém
`CLAUDE.md`: tři úpravy zdroje kvůli přenositelnosti a **povinné**
`-fno-automatic -finit-local-zero`, bez nichž běh tiše vrací samé `NaN`).

```
gfortran -std=legacy -fno-automatic -finit-local-zero \
  -o LEVMFORT LV0.FOR LV1.FOR LV2.FOR LV3.FOR LV4.FOR \
              LV5.FOR LV6.FOR LV7.FOR LV8.FOR LV9.FOR
```

Každý ze 145 spustitelných souborů zkopírován na `INFL` (s převodem CRLF → LF),
spuštěn s časovým limitem 60 s a prázdným stdin.

### Výsledky

| Ukazatel | Hodnota |
|----------|---------|
| Spuštěno | 145 |
| Nenulový návratový kód nebo timeout | **0** |
| Prázdný `PNTOUTL` | **0** |
| Vyroben `OUTIN` | 138 |
| `NaN` ve výstupu | 5 |
| Hlášení singulární matice | 6 |
| **Celkový čas** | **45,1 s** |
| Nejpomalejší | `CKT.E/TSTN36ft` (10,2 s) |

Další pomalé: `CKT.O/225Z36el` 6,3 s, `CKT.O/OTSP1Y36` 5,4 s,
`CKT.O/le225m36` 3,7 s. Všechny používají NELEM = 36 (KWW s cutoffem), který
počítá kvadraturou. Zbytek korpusu běží v desítkách milisekund.

### Rozložení hlášek o ukončení

| Hláška | Počet |
|--------|------:|
| `TERMINATION: XTOL IS TOO SMALL` | 118 |
| `TERMINATION: FTOL IS TOO SMALL` | 7 |
| `ALGORITHM TERMINATES ... RELATIVE ERROR IN THE SUM OF SQUARES` | 5 |
| `ALGORITHM TERMINATES ... FVEC IS ORTHOGONAL` | 5 |
| `TERMINATION: NUMBER OF CALLS TO FCN ... EXCEEDS MAXFEV` | 4 |
| `ALGORITHM TERMINATES ... ERROR BETWEEN THE VECTOR OF FREE PARAMETERS` | 3 |
| žádná | 3 |

> **`XTOL IS TOO SMALL` u 118 ze 145 testů není selhání.** Skoro všechny soubory
> v korpusu mají IOP = 0, takže XTOL = 10⁻⁴⁸ — mimo dosah dvojité přesnosti.
> Optimalizátor doiteruje, kam to v `double` jde, a ohlásí, že dál to nejde.
> Je to normální úspěšné ukončení v této konfiguraci. Kdo by chtěl testy
> vyhodnocovat podle návratového kódu, musel by nastavit IOP na −4 až −9.

### Problémové soubory

Všech 5 výskytů `NaN` a 5 ze 6 hlášení singulární matice připadá na **H-obvod**:

```
CKT.H/1MOBP25=1'BLOCKINGEPHI=1
CKT.H/1MOBP25=1'BLOCKINGZPHI=0.9
CKT.H/1MOBP25=1'PBLOCKING0.01ZPHI=0.9
CKT.H/1MOBP25=1'PBLOCKING0.01ZPHI=1
CKT.H/GEL25p53PNPorigdat64HP25=M1gpnpdeb3
```

Šesté hlášení je `CKT.H/HTST`. Tyto soubory zároveň nevyrobí `OUTIN`. Vzhledem
k poškozenému bloku RP/CP v `HSUB` (viz `CIRCUITS.md`, oddíl H) stojí za
prověření, jestli jde o vlastnost dat, nebo o důsledek té chyby.

Bez `OUTIN` skončily ještě `CKT.C/CTST` a `CKT.F/FTST`, ale z jiného důvodu:
mají na řádku 2 `IRE = 0`, kdežto `OUTIN` se zapisuje až při `IRE < -9`. Jsou to
jediné dva soubory v korpusu, které tuto volbu nemají nastavenou, takže se na
nich nedá provést test idempotence z části 5.2.

### Ověřený referenční příklad

`FITTESTS/CKT.A/AZC` proti hodnotám z vlastního titulku
(`EXACT PAR. VALUES: ZC(2D6,1D0,0.3),1D+6,1D-12`):

| P | Exaktní | Nafitováno | Rel. odchylka |
|---|---------|------------|---------------|
| 1 | 1,0×10⁶ | 9,9822×10⁵ | −1,8×10⁻³ |
| 6 | 2,0×10⁶ | 1,9916×10⁶ | −4,2×10⁻³ |
| 7 | 1,0 | 9,8491×10⁻¹ | −1,5×10⁻² |
| 9 | 0,3 | 2,9827×10⁻¹ | −5,8×10⁻³ |
| 29 | 1,0×10⁻¹² | 9,9996×10⁻¹³ | −4,0×10⁻⁵ |

`PDRMS = 7,93×10⁻³`. Odchylky jsou v rámci 1–2 odhadnutých směrodatných
odchylek. Manuál (str. 2-3) vysvětluje proč: **data v `AZC` jsou oříznutá na dvě
platná místa**.

Táž úloha s neoříznutými daty (`LEVMMISC/LEVMISC/AZCE`) dá:

| P | Exaktní | Nafitováno |
|---|---------|------------|
| 1 | 1,0×10⁶ | 1,0000×10⁶ |
| 6 | 2,0×10⁶ | 2,0000×10⁶ |
| 7 | 1,0 | 1,0000 |
| 9 | 0,3 | 3,0000×10⁻¹ |
| 29 | 1,0×10⁻¹² | 1,0000×10⁻¹² |

`PDRMS = 7,99×10⁻¹⁴` — shoda na 14 platných číslic, a to ze startu posunutého
o 1 % (u P(6) o −45 %). **`AZCE` je tedy nejpřísnější numerický test, jaký
projekt obsahuje**, a přitom neleží v `FITTESTS`, ale v `LEVMMISC/LEVMISC`.

---

## 9. Co v projektu chybí

| Mezera | Důsledek |
|--------|----------|
| **Žádné referenční výstupy** | není s čím porovnat `PNTOUTL`; regrese se pozná jen okem |
| **Žádné strojově čitelné očekávané hodnoty** | exaktní parametry jsou volný text v titulku |
| **Žádné pass/fail** | návratový kód je 0 i při `NaN` a singulární matici |
| **S a T bez testů** | `SSUB` je implementovaný a nikdy neběžel; `TSUB` je stub |
| **Vypnutá kontrola `PDAV`** | varování popsané v manuálu nikdy nevyskočí (část 7) |
| **`MKLIO` a `EXTRLVM` bez zdrojů** | skládání vstupů i extrakci výsledků nelze na dnešním systému spustit |
| **Interaktivní `PAUSE`** | řada chybových cest čeká na `go` od uživatele; blokuje dávkové běhy |
| **Tvarování vah vyžaduje vstup z klávesnice** | `IGACC < 0` se ptá na čtyři frekvence (`LV0.FOR:614-621`), takže testovat ho lze jen ručně |
| **Pevná jména souborů v CWD** | dva běhy nelze pustit paralelně ve stejném adresáři |
| **Nerovnoměrné pokrytí** | 61 testů pro O-obvod, po jednom pro C, F, I, J |

---

## 10. Jak z toho udělat regresní test

Není součástí projektu, ale z výše zjištěného to je přímočaré — korpus má
oracly i dostatečnou rychlost (45 s na 145 fitů).

**Nejlevnější varianta.** Vzít soubory s exaktními hodnotami v titulku
(`AZCE`, `AZC1`, `AZC`, `BTST`, `CTST`, `ETST`, `FTST`, `RTSTWC*`, `RTSTT.9`)
a porovnávat blok `PARAMETER ESTIMATES` proti hodnotám vypsaným v titulku,
s tolerancí odvozenou z `ESTIMATED RELATIVE STD DEV`. `AZCE` samo o sobě
odhalí prakticky každou regresi v jádru — je citlivé na 14 platných číslic.

**Systematičtější varianta.** Vygenerovat jednorázově referenční `PNTOUTL`
pro všech 145 souborů a porovnávat numericky s tolerancí (textově to nejde,
poslední číslice se mezi překladači liší). Pak stačí sledovat `PDRMS`, `SD RC`
a `FQF`.

**Co je nutné ošetřit v obou případech:**

- převést CRLF na LF (Fortran čte pevné sloupce, ale `RNL.BAT` to nedělal)
- před každým během smazat výstupní soubory, jak to dělá `RNL.BAT`
- pouštět s prázdným stdin a timeoutem kvůli `PAUSE`
- nevyhodnocovat podle návratového kódu — ten je vždy 0; číst `PNTOUTL`
- `TERMINATION: XTOL IS TOO SMALL` nepovažovat za chybu (část 8)
- hlídat `NaN` a řetězec `SINGULAR` jako tvrdé selhání

Pro S-obvod a nové testy dalších obvodů lze data vyrobit round-trip postupem
z části 5.1: nastavit `MAXFEV = 0`, zvolit parametry, vygenerovat data, pak je
fitovat zpět a doplnit exaktní hodnoty do titulku podle zavedené konvence.
