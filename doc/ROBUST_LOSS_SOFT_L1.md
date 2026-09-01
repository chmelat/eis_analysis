# Robustní ztrátová funkce (soft-L1): proč ji zavést

Datum: 2026-09-01
Stav kódu: v0.28.0
Souvisí s: `doc/ZSCOPE_COMPARISON.md` (položka "Robustní loss + sandwich
kovariance", priorita vysoká), `doc/WEIGHTING_AND_STATISTICS.md`

---

## 1. Shrnutí

Fitování obvodů dnes minimalizuje součet čtverců vážených reziduí. To znamená,
že vliv bodu roste s **druhou mocninou** jeho odchylky, takže jediný
kontaminovaný bod může přebít desítky dobrých. Na reálném měření čistého
rezistoru 10.9 MOhm to posunulo výslednou hodnotu o 19.6 %; se soft-L1 zbyla
odchylka 1.4 %, a to bez jakéhokoli ručního zásahu do dat.

Cena na čistých datech je přitom prakticky nulová: v kontrolovaném pokusu se
medián chyby parametru zhoršil z 0.10 % na 0.12 %.

Doporučení: zavést `loss='soft_l1'` jako volitelný režim fitování
(`--robust-loss`), s `f_scale` odvozeným z MAD reziduí prvního průchodu.

---

## 2. Problém: jeden bod přepíše výsledek

Reziduum se v `circuit.py:365` skládá takto:

```
r_i = (Z_i - Z_model,i) * w_i,      w_i = 1/|Z_i|  (modulus, výchozí)
```

a `scipy.optimize.least_squares` minimalizuje `0.5 * sum(r_i^2)`. Kvadrát je
tu podstatný: bod s reziduem 10 sigma přispěje do gradientu stokrát víc než
bod s 1 sigma. Estimátor nejmenších čtverců má bod zvratu (breakdown point)
roven nule, tedy jediné dostatečně vzdálené pozorování dokáže odsunout
výsledek libovolně daleko.

U modulus vážení je to ještě citlivější směrem dolů. Bod, jehož naměřené
`|Z|` je **menší** než modelové, má malé `|Z_i|`, tedy velkou váhu `w_i`, a
navíc velké reziduum. Obojí se násobí.

### Naměřený případ: rezistor 11 MOhm

Soubor `EISPOT-sada_rezistoru.DTA` (potenciostat REF620, 10 mV rms, 72 bodů).
Vzorek je obyčejný rezistor, takže správná odpověď je známá: medián `Re(Z)`
pod 10 kHz je **1.0925e7 Ohm**.

Jeden bod ve spektru leží na 50.22 Hz, tedy na síťové frekvenci, a má
`Re = 2.24e6` místo očekávaných 1.09e7 - odchylka 79.5 % při mediánové
odchylce všech ostatních bodů 1.0 %. Je to nejhorší bod z 51 v pásmu pod
10 kHz.

Fit `R()|C()` na pásmu pod 20 kHz:

| postup | R [Ohm] | odchylka | chyba fitu |
|---|---|---|---|
| kvadratická ztráta (dnešní stav) | 8.7849e6 | 19.6 % | 26.42 % |
| ruční smazání bodu na 50 Hz | 1.0829e7 | 0.9 % | 7.93 % |

Jeden bod ze 54 posunul výsledek o pětinu. A protože model tím sedí o 20 %
vedle, vykazuje pak **každý dobrý bod** reziduum kolem 20 %, což nafoukne
celkovou chybu z 8 % na 26 %. Diagnostika pak vypadá, jako by byla vadná
celá sada, ačkoli je vadný jeden bod.

---

## 3. Co je soft-L1

`least_squares` minimalizuje `0.5 * sum(rho(z_i))`, kde `z_i = (r_i/f_scale)^2`.
Výchozí `rho(z) = z` dává klasické nejmenší čtverce. Soft-L1 je

```
rho(z) = 2 * (sqrt(1 + z) - 1)
```

s chováním:

- pro `|r| << f_scale` platí `rho(z) ~ z`, tedy shodné s nejmenšími čtverci;
- pro `|r| >> f_scale` platí `rho(z) ~ 2*sqrt(z) = 2*|r|/f_scale`, tedy růst
  **lineární** místo kvadratického.

Je to hladká (všude diferencovatelná) aproximace L1 normy, odtud "soft".
Influence funkce `psi(r) = rho'(r) = r / sqrt(1 + (r/f_scale)^2)` je omezená:
ať je bod jakkoli daleko, jeho příspěvek do gradientu je shora ohraničený.
To je celá podstata - odlehlý bod se nemaže, jen přestane rozhodovat.

### Volba f_scale

`f_scale` je práh přechodu mezi kvadratickým a lineárním režimem, v jednotkách
rezidua. Rozumná automatická volba je robustní odhad rozptylu reziduí z
prvního (neroubustního) průchodu:

```
f_scale = 1.4826 * MAD(r),     MAD(r) = median(|r - median(r)|)
```

Konstanta 1.4826 dělá z MAD nestranný odhad směrodatné odchylky pro normální
rozdělení. Měření citlivosti na této volbě (rezistor 11 MOhm, pásmo pod
20 kHz):

| f_scale | R [Ohm] | odchylka | chyba fitu |
|---|---|---|---|
| 1x MAD | 1.0767e7 | 1.4 % | 13.98 % |
| 2x MAD | 1.0651e7 | 2.5 % | 14.42 % |
| 3x MAD | 1.0531e7 | 3.6 % | 14.95 % |
| 5x MAD | 1.0307e7 | 5.7 % | 16.14 % |
| 10x MAD | 9.8397e6 | 9.9 % | 19.08 % |

Chování je monotónní a bez ostrých zlomů: čím větší `f_scale`, tím blíž ke
kvadratické ztrátě, jak se očekává. `1x MAD` je zde nejlepší a dostane se na
1.4 % proti 19.6 % u dnešního stavu, tedy zachytí většinu toho, co dá ruční
smazání bodu (0.9 %) - ale automaticky a bez zásahu do dat.

---

## 4. Naměřený přínos a cena

Anekdota nestačí, proto kontrolovaný pokus se známou pravdou. Obvod
`R(10)-(R(1000)|C(1e-6))`, 50 bodů v rozsahu 0.1 Hz - 100 kHz, 1 % šumu,
40 nezávislých realizací. Kontaminace = náhodně vybrané body vynásobené
faktorem 1.5. Sleduje se medián absolutní chyby parametru `R_ct`:

| kontaminovaných bodů | kvadratická ztráta | soft-L1 (f_scale = 1x MAD) |
|---|---|---|
| 0 % | 0.10 % | 0.12 % |
| 2 % | 0.63 % | 0.12 % |
| 5 % | 0.91 % | 0.13 % |
| 10 % | 2.30 % | 0.19 % |

Dvě věci, které z toho plynou:

**Cena na čistých datech je zanedbatelná.** Rozdíl 0.10 % vs. 0.12 % je řádově
pod nejistotou, se kterou se v EIS pracuje. Robustní ztráta tedy neplatí
znatelnou daň za pojistku, kterou poskytuje. To je hlavní argument pro to,
aby mohla být zapnutá i rutinně.

**Přínos roste s kontaminací.** Při 10 % vadných bodů je chyba parametru
dvanáctkrát menší. Přitom 10 % není nerealistické číslo - stačí síťové
rušení plus pár bodů na okrajích rozsahu.

---

## 5. Co soft-L1 neřeší

Aby nevznikla falešná očekávání, tohle robustní ztráta **nespraví**:

**Drift a nestacionaritu.** Odchylky na nízkých frekvencích způsobené tím, že
se vzorek během skenu měnil, nejsou náhodné odlehlé body, ale systematický
trend. Robustní ztráta je potlačí stejně ochotně jako chybu měření, čímž
zamaskuje reálnou vlastnost experimentu. Na tohle slouží Z-HIT a KK
diagnostika, ne jiná ztrátová funkce.

**Špatný model.** Fit `R()-L()` na data, jejichž `|Z|` klesá, selže s jakoukoli
ztrátovou funkcí, protože model takový průběh neumí vytvořit. Robustní ztráta
neopraví strukturální nesoulad.

**Většinovou kontaminaci.** Soft-L1 má bod zvratu výrazně lepší než nula, ale
ne 50 %. Když je vadná většina bodů, "odlehlé" jsou najednou ty dobré. Na to
by byl potřeba LTS (least trimmed squares), což je pro rozsah tohoto projektu
zbytečné - většinová kontaminace znamená vadné měření, ne potřebu lepšího
estimátoru.

**Nedostatečnou citlivost přístroje.** V případě rezistoru 11 MOhm teče při
10 mV rms proud 0.92 nA, což je na hranici rozlišení. Robustní ztráta zlepší
odhad z těchto dat, ale nenahradí zvýšení amplitudy.

---

## 6. Návaznost: kovariance

Kovariance se dnes počítá jako `s^2 * (J^T J)^-1`
(`eis_analysis/fitting/covariance.py:68`). Ta formule platí za předpokladu, že
rezidua jsou nezávislá, stejně rozdělená a gaussovská s rozptylem, který
vážení předpokládá. Odlehlé body ten předpoklad porušují, takže intervaly
spolehlivosti vycházejí příliš úzké - u ZScope měřili pokrytí 89 % místo
nominálních 95 % při 5 % kontaminovaných bodů.

Sandwich estimátor ten předpoklad nepotřebuje:

```
Cov = A^-1 B A^-1

A = sum psi'(r_i) * J_i J_i^T     (bread - očekávaná křivost)
B = sum psi(r_i)^2 * J_i J_i^T    (meat  - empirický rozptyl skóre)
```

Pokud jsou rezidua opravdu i.i.d. gaussovská a ztráta lineární, pak
`psi(r) = r`, `psi' = 1`, `B = s^2 * A` a výraz se zkrátí zpět na
`s^2 (J^T J)^-1`. Sandwich je tedy striktní zobecnění, ne alternativa.

Obojí spolu souvisí, ale **nemusí se zavádět najednou**. Soft-L1 sám o sobě
zlepší bodové odhady parametrů; sandwich sám o sobě zlepší jejich nejistoty.
Doporučené pořadí je soft-L1 první, protože chybný bodový odhad je horší
problém než optimistický interval kolem něj.

---

## 7. Implementace

### Kde zasáhnout

Ztrátová funkce musí být konzistentní na všech místech, která hledají optimum,
jinak by jedna fáze táhla k jinému řešení než druhá:

| místo | co to je | zásah |
|---|---|---|
| `fitting/circuit.py:419` | hlavní `least_squares` | `loss=`, `f_scale=` |
| `fitting/diffevo.py:69` | cílová funkce DE (`sum r^2`) | aplikovat stejné `rho` |
| `fitting/diffevo.py:395` | doleštění po DE | `loss=`, `f_scale=` |
| `fitting/multistart.py` | volá `fit_equivalent_circuit` | zdědí automaticky |

Pozor na `diffevo.py:69`: kdyby DE minimalizovalo prostý součet čtverců a
robustní byla jen doleštovací fáze, DE by přistálo v pánvi vychýlené
odlehlými body a `least_squares` by z ní už nevystoupilo.

### Past ve scipy

`res.jac`, které `least_squares` vrací, **není** čistý Jacobián reziduí, když
je zapnutá robustní ztráta. Funkce `scale_for_robust_loss_function`
(scipy 1.17.1, `optimize/_lsq/common.py`) ho násobí `sqrt(rho' + 2*rho''*f^2)`
na místě, zatímco `res.fun` zůstává nescalované.

Důsledek pro `circuit.py:467`, kde se `opt_result.jac` předává do
`compute_covariance_matrix`: po přepnutí na soft-L1 by ta matice už byla
"bread" A ze sandwiche, zatímco "meat" B by se musela dopočítat zvlášť z
`res.fun` a nescalovaného Jacobiánu. Kdyby se to pustilo naivně jako dnes,
vyšlo by něco třetího - ani naivní kovariance, ani sandwich.

Nejjednodušší bezpečná varianta pro první krok: při zapnuté robustní ztrátě
počítat kovarianci z Jacobiánu vyčísleného zvlášť (`make_jacobian_function`
z `fitting/jacobian.py`) a hlásit, že intervaly jsou v tom režimu
orientační, dokud nepřibude sandwich.

### Navrhované CLI

```
--robust-loss              zapne soft_l1 s f_scale = 1.4826*MAD reziduí
                           prvního průchodu
--robust-f-scale NÁSOBEK   ruční přenásobení té hodnoty (výchozí 1.0)
```

Výchozí stav vypnutý, dokud se chování neověří na širší sadě měření. Výstup
fitu by měl uvést, že režim byl aktivní, a kolik bodů dostalo váhu pod
nějakou hranicí - to je zároveň diagnostika, která ukáže, které body ztráta
potlačila.

### Ověření

Prototyp použitý pro čísla v tomto dokumentu reprodukuje
`fit_equivalent_circuit` na čtyři platné číslice, když se `loss` nechá na
`'linear'` (stejná vážení, stejný analytický Jacobián, stejné meze, `x_scale`
i metrika chyby). To je minimální podmínka, aby byla čísla porovnatelná; při
implementaci by měl stejnou roli plnit regresní test.

---

## 8. Reference

- Huber, P. J. (1964). "Robust Estimation of a Location Parameter."
  *Annals of Mathematical Statistics* 35, 73-101.
- Rousseeuw, P. J., Leroy, A. M. (1987). *Robust Regression and Outlier
  Detection.* Wiley.
- scipy.optimize.least_squares, parametry `loss` a `f_scale`.
- `doc/ZSCOPE_COMPARISON.md` - položka, ze které tento návrh vychází.
- `doc/WEIGHTING_AND_STATISTICS.md` - definice vážení a metrik chyby.
