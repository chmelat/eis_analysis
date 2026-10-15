# Robustní ztrátová funkce (soft-L1): proč ji zavést

Datum: 2026-09-01, aktualizace 2026-10-10 a 2026-10-15
Stav kódu: v0.28.0, odkazy a měření ověřeny na v0.58.0; sekce 4.1-4.2
na 508c06a (v0.60.3 a nevydané změny kola 6)
Stav návrhu: **odloženo (2026-10-15)**, neimplementováno. Nálezy cold
agentů T13/T14 řeší přepínač `--exclude-flagged` (878b498, sekce 4.3);
soft-L1 se vrátí na pořad při podmínkách v 4.3
Souvisí s: `doc/ZSCOPE_COMPARISON.md` (položka "Robustní loss + sandwich
kovariance", priorita vysoká), `doc/WEIGHTING_AND_STATISTICS.md`,
`doc/STRESS_TEST.md` (pokrytí CI), `eis_analysis/validation/outliers.py`,
`doc/COLD_AGENT.md` (kolo 6)

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

Dva cold agenti (2026-10-15, sekce 4.1) ukázali, že nejde jen o jeden vadný
bod u rezistoru: u dvou blízkých oblouků posune 5-6 % vadných bodů rozdělení
polarizačního odporu o 25-50 %, a agent to neopraví, ani když ho knihovna
na vadné body upozorní. Soft-L1 v prototypu vrátí výsledek k pravdě (107/73
místo 142/38 při pravdě 100/80). Alternativa - detekovat vadné body
z reziduí fitu a fitovat bez nich (N3) - fungovala také, ale byla zamítnuta
(sekce 4.2).

**Odloženo 2026-10-15 (sekce 4.3):** jednodušší alternativa, přepínač
`--exclude-flagged`, který vyřadí body z kontroly po bodech, vyřešila T13
i se slabším agentem (Haiku: 105/75) a nechává statistiku fitu
standardní. Soft-L1 by navíc přinesl jen potlačení vadných bodů, které
kontrola nevypíše; to dnes neodůvodňuje zásah do jádra.

---

## 2. Problém: jeden bod přepíše výsledek

Reziduum se v `circuit.py:366-378` skládá takto:

```
r_i = (Z_i - Z_model,i) * w_i,      w_i ~ 1/|Z_i|  (modulus, výchozí)
```

Váhy jsou od v0.28.0 normované na `mean(w |Z|) = 1`
(`compute_residual_weights`, `fitting/diagnostics.py:57`), rezidua jsou tedy
v Ohm. Na argumentaci ani na volbu `f_scale` níže to vliv nemá: `f_scale` se
odvozuje ze stejných reziduí.

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

Ověřeno znovu na v0.58.0 (2026-10-10, soubor
`eis_analysis_data/R_11M/EISPOT-sada_rezistoru.DTA`, výchozí LM fit): obě
řádky tabulky vycházejí na všechny uvedené číslice stejně.

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

### 4.1 Blízké oblouky a řídké vadné body (cold agenti T13, T14)

Úlohy kola 6 (`doc/COLD_AGENT.md`), pravda a generátor v
`~/coldbox-private/gen_r6.py` a `gen_r6b.py`:

| úloha | data | vadné body |
|---|---|---|
| T13 | R_s 2 - ZARC(100, 1 ms, 0.85) - ZARC(80, 3.16 ms, 0.80), 81 bodů, šum 0.1 % | 5, Gauss sigma 2.5 % na Re i Im |
| T14 | R_s 3 - ZARC(120, 2 ms, 0.82) - ZARC(90, 7.96 ms, 0.78), 97 bodů, šum 0.15 % | 6, totéž |

Oblouky jsou 0.5, resp. 0.6 dekády od sebe, takže R1 a R2 jsou silně
antikorelované (r = -0.9996 u T14) a jejich součet je určený dobře, ale
rozdělení ne. Pár vadných bodů pak rozdělení posune daleko za uvedenou
nejistotu, přestože fit sedí "na šum" (poměr fit/šum 1.0, protože odhad
šumu Lin-KK obsahuje tytéž body):

| fit dvou ZARC | T13: R1 / R2 [Ohm] | T14: R1 / R2 [Ohm] |
|---|---|---|
| pravda | 100 / 80 | 120 / 90 |
| nejmenší čtverce (dnes) | 142 +- 8 / 38 +- 8 | 150 +- 11 / 60 +- 11 |
| soft-L1, f_scale = 1.4826 MAD | 107 / 73 | 128 / 82 |
| nejmenší čtverce bez vadných bodů | 105 / 75 | 122 / 89 |
| nejmenší čtverce, stejná data bez vadných bodů | - | 127 / 83 |

Soft-L1 dostane fit zpět na úroveň, kterou omezuje šum (poslední řádek:
čistá data T14 dávají 127/83, fit bez vadných bodů má nejistotu R +-5.5). Prototyp:
`scipy.optimize.least_squares(loss='soft_l1')` na reziduích s vahou 1/|Z|,
start z fitu nejmenších čtverců, `f_scale` z jeho reziduí; není to
implementace podle sekce 7 (bez DE, bez kovariance).

**Špatný model se nezamaskuje.** Fit jednoho ZARC (chybí oblouk) má
rezidua robustního fitu stejně strukturovaná jako u nejmenších čtverců:
rho1 0.68 / 0.80 (T13) a 0.87 / 0.72 (T14), správný model -0.02 / +0.02
a -0.04 / -0.02. Stávající diagnostika tvaru reziduí
(`fitting/residual_diagnostics.py`) tedy na robustním fitu čte model
správně. Na nejmenších čtvercích to nejde: správný model stažený vadnými
body nechá na ostatních bodech hladký nesoulad (T13: rho1 0.74 / 0.55 po
vynechání vadných bodů), stejně jako by ho nechal špatný model. Netestováno:
nesoulad omezený na pár bodů na kraji okna (chybějící L, sériové C), který
by soft-L1 mohl potlačit jako vadné body (sekce 5).

**Upozornění nestačí.** Oba agenti (Haiku) určili počet procesů a R_s
správně, R1/R2 ne. Agent T14 měl k dispozici varování, že odhad šumu táhne
pár bodů (0.72 % proti 0.29 % z mediánu), i seznam podezřelých bodů
z kontroly po bodech, obojí citoval, a body nevyřadil: "Podezřelé body
z kontroly po bodech jsem nevyřazoval." Uživatel, který má z výstupu
udělat ruční zásah, ho neudělá. Robustní režim musí být jedním
přepínačem, na který varování přímo ukážou (sekce 7).

### 4.2 Proč ne detekce vadných bodů a nový fit (N3, zamítnuto)

Alternativa zkoušená 2026-10-15 (větev `n3-isolated-refit`, smazána):
z reziduí fitu vybrat body nad 6x medián |dZ|/|Z|, jen pokud jsou rezidua
jinak ve tvaru šumu (|rho1| < 0.3), fit bez nich zopakovat a podle tvaru
jeho reziduí rozhodnout, zda šlo jen o vadné body, nebo i o špatný model.
Na T14 dala 122/89, tedy správný výsledek přímo v agentově výstupu.

Zamítnuta, protože:

- **Složitost:** tři kalibrované prahy (faktor 6, |rho1| 0.3, kritérium
  nového fitu), fit uvnitř diagnostiky, výpis v CLI se třemi větvemi, ~260
  řádků a nové veřejné API.
- **Verdikt chybuje v obou směrech:** správný model označil za správný
  v 88 % případů, špatný model zakrytý vadnými body odhalil v 86 %.
- **Kalibrace a ověření stály na jednom modelu vadných bodů:** 6 % bodů
  rovnoměrně, Gauss sigma 2.5 %, šum úměrný |Z| - na tomtéž modelu jako
  T13 a T14, takže úspěch na T14 nebyl nezávislý. Shluky (síťové rušení),
  jeden velký bod, chyba jen ve fázi, těžké chvosty a konstantní šum
  netestovány.
- **Soft-L1 je standardní estimátor** (Huber 1964), ne heuristika laděná
  na vlastní úloze. Nevyžaduje rozhodnutí "vadný bod, nebo model" a podle
  4.1 dělá totéž (107/73 proti 105/75 u T13).

Zůstává z toho: kalibrační skripty (stress spektra, tři špatné obvody,
s vadnými body i bez) se hodí k ověření soft-L1, a poznatek, že tvar
reziduí fitu nejmenších čtverců nerozliší vadné body od špatného modelu.

---

### 4.3 Odloženo: `--exclude-flagged` jako alternativa (2026-10-15)

Přepínač `--exclude-flagged {kk,both,any}` (878b498, README "Per-point
residual check") vyřadí body, které vypíše kontrola po bodech
(`validation/outliers.py`; výchozí `kk`: body od Lin-KK, protože Z-HIT
rozmaže vadný bod do sousedů), ze všeho za validací, stejnou maskou jako
`--f-min/--f-max`. Kontrola po bodech a varování u reziduí ho jmenují.

Výsledky (cold agenti, `doc/COLD_AGENT.md`, kolo 6):

| úloha | fit všech bodů | `--exclude-flagged kk` | soft-L1 (prototyp, 4.1) | pravda |
|---|---|---|---|---|
| T13 | 142 / 38 | 105 / 75 (5 bodů) | 107 / 73 | 100 / 80 |
| T14 | 150 / 60 | 123 / 87 (2 body) | 128 / 82 | 120 / 90 |

Haiku, který v krocích 1-2 vadné body ponechal (jen varování), s
přepínačem vyřešil T13 za 0.11 USD (krok 4); Opus si stejný zásah dřív
skriptoval sám.

Proč stačí (zatím):

| | `--exclude-flagged` | soft-L1 |
|---|---|---|
| body, které kontrola vypíše | vyřadí | potlačí |
| body, které kontrola nevypíše (T14: 4 ze 6, chyba ~10-20x šum) | zůstanou, táhnou dál | potlačí |
| drift v nejnižší dekádě | varuje, rozhodne uživatel | potichu potlačí, drift zamaskuje (sekce 5) |
| +- a CI, AIC/BIC, poměr fit/šum | zůstávají standardní | sandwich, robustní srovnání modelů, nová definice (sekce 6, 7) |
| zásah do kódu | hotovo (maska rozsahu analýzy) | šest míst v jádru vč. DE a archivu, past ve scipy (sekce 7) |
| průhlednost | výpis vyřazených bodů, záznam v JSON fitu | plynulé váhy |

Přínos soft-L1 navíc je tedy jen v bodech, které jsou na detekci malé, ale
fit ještě táhnou. Na T14 i tak stačilo vyřadit dva největší. Případ
rezistoru 11 MOhm (sekce 2: jeden bod o 80 %) pokryje přepínač také - takový
bod je nad absolutním prahem 5 %. Na měřených spektrech v repozitáři (ZrO2,
EISPOT-test1) převažuje nesoulad modelu, ne vadné body.

Vrátit se k soft-L1, když:

- se na měřených datech ukáže, že po `--exclude-flagged` zůstává podstatné
  vychýlení od bodů, které kontrola nevypíše (shluky síťového rušení, mnoho
  menších vadných bodů, chyba jen ve fázi);
- měření na kalibračních spektrech (fit všech bodů / `--exclude-flagged kk`
  / soft-L1, chyba parametrů; generátor vadných bodů, špatné obvody
  a prototyp soft-L1 z kola 6 v `~/coldbox-private/r6/calib/`, mimo
  repozitář) ukáže velký rozdíl mezi druhou a třetí variantou - zatím
  neprovedeno;
- bude potřeba sandwich kovariance tak jako tak (pokrytí CI 0.90 i bez
  vadných bodů, sekce 6) - pak odpadne velká část ceny.

Zbytek dokumentu (sekce 5-7) zůstává jako připravený návrh implementace.

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

**Lokální nesoulad modelu.** Chybí-li prvek, jehož projev je omezen na pár
bodů na kraji okna (sériová L na nejvyšších frekvencích, sériové C na
nejnižších), soft-L1 může tyto body potlačit stejně jako vadné a fit pak
vypadá správně. Výpis bodů se sníženou vahou (sekce 7) to musí ukázat:
souvislá řada bodů na kraji okna je model nebo drift, ne vadné body.
Ověřit (viz Ověření).

**Většinovou kontaminaci.** Soft-L1 má bod zvratu výrazně lepší než nula, ale
ne 50 %. Když je vadná většina bodů, "odlehlé" jsou najednou ty dobré. Na to
by byl potřeba LTS (least trimmed squares), což je pro rozsah tohoto projektu
zbytečné - většinová kontaminace znamená vadné měření, ne potřebu lepšího
estimátoru.

**Rozpor s reportem odlehlých bodů.** Ve stejné verzi (v0.28.0) vznikl
`validation/outliers.py` (`--max-residual`): z reziduí Lin-KK a Z-HIT
vypíše podezřelé body, ale záměrně je nemaže ani jim nesnižuje váhu, protože
odchylky na nízkých frekvencích bývají drift. Soft-L1 dělá právě to, co ten
modul odmítá: potichu snižuje váhu. Aby se rozpor neztratil, musí fit
v robustním režimu vypsat body, kterým ztráta výrazně snížila váhu (viz
sekce 7), stejným způsobem jako report odlehlých bodů, a uživatel musí
vidět, jestli leží izolovaně (vadný bod), nebo v trendu k nízkým frekvencím
(drift, robustní ztráta nevhodná).

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

Vlastní měření (stress test, invariant F3, `doc/STRESS_TEST.md`, známý
limit 9, 2026-10-09): linearizovaný 95% interval obsahuje pravdu u 0.90
parametrů, a to **na datech bez odlehlých bodů** (proporcionální šum,
dobře podmíněné fity), na všech úrovních šumu. Intervaly jsou tedy úzké
i bez kontaminace, nejspíš kvůli nelinearitě modelu nebo nesouladu vážení
a skutečného šumu. Sandwich zachytí to druhé, nelinearitu ne. Soft-L1 sama
pokrytí nezlepší vůbec, to je úloha kovariance. F3 je hotové měřidlo:
jakákoli změna kovariance se dá ověřit tím, jak posune 0.90 k 0.95.

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
| `fitting/circuit.py:416` | hlavní `least_squares` (přes `least_squares_normalized`, `fitting/optimizer.py`, který kwargs propustí) | `loss=`, `f_scale=` |
| `fitting/diffevo.py:98` | cílová funkce DE (`sum r^2`), ze které se plní i archiv kandidátů | aplikovat stejné `rho` |
| `fitting/diffevo.py:440` | doleštění po DE, jedno nastavení pro bod DE i kandidáty archivu | `loss=`, `f_scale=` |
| `fitting/diffevo.py:457, 478` | cena doleštěného fitu `sum(r.fun**2)` | robustní cena (`r.cost`, nebo `sum rho`) |
| `fitting/de_archive.py:168` | `choose()`: výběr mezi bodem DE, jeho doleštěním a archivem, `ADOPT_MARGIN`, nejednoznačné alternativy | porovnávat robustní ceny |
| `fitting/multistart.py` | volá `fit_equivalent_circuit` | zdědí automaticky |

Pozor na `diffevo.py:98`: kdyby DE minimalizovalo prostý součet čtverců a
robustní byla jen doleštovací fáze, DE by přistálo v pánvi vychýlené
odlehlými body a `least_squares` by z ní už nevystoupilo.

Stejná past je o krok dál ve výběru (`diffevo.py:457, 478`,
`de_archive.choose`): od v0.53.0 se po DE doleštuje víc startů a vítěz se
vybírá podle kvadratické ceny. S robustní ztrátou by tak mohl vyhrát
kandidát, kterého táhnou odlehlé body, přestože robustní doleštění našlo
lepší. Komentář v `diffevo.py:481-484` to říká obecně: výběr musí používat
právě optimalizovanou cenu. Kvadratická zůstává jen `cost_floor`
(zaokrouhlovací práh, nezávislý na ztrátě).

### Past ve scipy

`res.jac`, které `least_squares` vrací, **není** čistý Jacobián reziduí, když
je zapnutá robustní ztráta. Funkce `scale_for_robust_loss_function`
(scipy 1.17.1, `optimize/_lsq/common.py`; ověřeno 2026-10-10 i v 1.10.1
z apt a 1.18.1 ve vývojovém venv) ho násobí `sqrt(rho' + 2*rho''*f^2)` na místě, zatímco
`res.fun` zůstává nescalované. `least_squares_normalized` navíc `res.jac`
přepočítá zpět do původních proměnných (`jac / scale`), takže škálovaný
Jacobián by prošel dál a nic by na něj neupozornilo.

Důsledek pro `circuit.py:458`, kde se `opt_result.jac` předává do
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

Varování, která dnes na vadné body upozorňují, mají `--robust-loss` jmenovat
(sekce 4.1: agent upozornění četl a nic neudělal): varování Lin-KK, že odhad
šumu táhne pár bodů (`kramers_kronig.py`, `NOISE_DOMINATED_RATIO`),
seznam podezřelých bodů (`cli/handlers/validation.py`, `report_outliers`)
a varování "Residuals are not random" (`cli/handlers/fitting.py`).

Dál zasáhne robustní ztráta do míst, která dnes počítají se součtem
čtverců: srovnání kandidátů AIC/BIC (`model_comparison.py`) musí v robustním
režimu porovnávat robustní cenu, poměr fit/šum (`compare_fit_to_noise`)
a `fit_error_rel` odpovídají kvadratické ztrátě a je třeba rozhodnout, zda je
počítat ze všech bodů, nebo jen z bodů s plnou vahou.

### Ověření

Prototyp použitý pro čísla v tomto dokumentu reprodukuje
`fit_equivalent_circuit` na čtyři platné číslice, když se `loss` nechá na
`'linear'` (stejná vážení, stejný analytický Jacobián, stejné meze, `x_scale`
i metrika chyby). To je minimální podmínka, aby byla čísla porovnatelná; při
implementaci by měl stejnou roli plnit regresní test.

Ověření naslepo (poučení z 4.2): `f_scale` a hranici pro výpis bodů se
sníženou vahou určit na jednom typu kontaminace a ověřit na mechanismech,
na kterých se nic neladilo - shluky sousedních bodů (síťové rušení),
jeden bod o 10-30 %, chyba jen ve fázi nebo jen v |Z|, těžké chvosty
(Studentovo t), konstantní absolutní šum, drift spolu s vadnými body,
chybějící prvek na kraji okna. Kritéria stanovit předem, např.: na
čistých datech se parametry nezhorší nad 0.12 % mediánové chyby (sekce 4),
s vadnými body je robustní fit blíž pravdě než nejmenší čtverce aspoň
v 80 % případů, špatný model zůstane ve tvaru reziduí viditelný. Reálná
data se známými vadami: rezistor 11 MOhm (sekce 2), ZrO2 21 d z 16.7.2026
s body označenými THD. Nakonec cold agent na úloze, kterou nenavrhl autor
implementace.

---

## 8. Reference

- Huber, P. J. (1964). "Robust Estimation of a Location Parameter."
  *Annals of Mathematical Statistics* 35, 73-101.
- Rousseeuw, P. J., Leroy, A. M. (1987). *Robust Regression and Outlier
  Detection.* Wiley.
- scipy.optimize.least_squares, parametry `loss` a `f_scale`.
- `doc/ZSCOPE_COMPARISON.md` - položka, ze které tento návrh vychází.
- `doc/WEIGHTING_AND_STATISTICS.md` - definice vážení a metrik chyby.
