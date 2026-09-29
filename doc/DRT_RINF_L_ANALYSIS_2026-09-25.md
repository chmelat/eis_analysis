# R_inf a L jako proměnné DRT výpočtu - analýza

**Datum:** 2026-09-25 | **Verze:** 0.37.0 | **Typ:** numerická studie proveditelnosti (Claude Code)

**Otázka:** Je lepší do DRT výpočtu zahrnout R_inf (případně indukčnost L)
jako proměnné, místo externího odhadu R_inf?

**Souvislosti:** hodnotí návrh `DRT_IMPROVEMENTS.md` bod 2 ("Fit R_inf jako
parametr") a navazuje na `archive/DRT_MATH_AUDIT_2026-06-27.md` F8 (R_inf jako medián)
a F9 (rozsah tau-mřížky) a na `archive/AUDIT_ri_fit_2026-09-25.md` (externí odhad
`--ri-fit`). Produkční kód se nemění, všechna čísla pocházejí z prototypu
mimo repozitář (popis níže).

**Stav k verzi 0.41.0:** D1 je hotové (`calculate_drt(inductance='auto')`,
CLI `--drt-inductance`). Chyba rekonstrukce A2 / C / C2 s přesným R_inf klesla
na 1.3 / 1.3 / 1.2 %, flat a B se nezměnily. Měření před implementací viz
"Výsledky k D1 (0.41.0)" níže.

**Stav k verzi 0.40.0 (revize 2026-09-25):**
- **D1 platí beze změny.** Přeměřeno na 0.40.0 s přesným R_inf (5 seedů):
  chyba rekonstrukce A2 / C / C2 je 5.8 / 17.3 / 14.9 % při výchozím vážení
  `sqrt` a 5.3 / 18.1 / 15.3 % bez vážení. Vážení z 0.38 problém neřeší.
- **D4 je hotové** (0.40.0): `--ri-fit` je fit R-L-(R|Q) s kontrolou
  identifikovatelnosti. Sloupce `--ri-fit` v tabulkách níže popisují *starý*
  estimátor. Nový dává na stejných případech C -0.5 %, C2 -1.0 %, flat -0.7 %,
  B -0.5 % a `real_gamry_example.DTA` označí jako neurčitelný (medián
  s varováním). Společný odhad (D2) tím proti `--ri-fit` na určitelných
  případech prohrává všude.
- Opraveny dva závěry, které vlastní data nepodporovala: mez R_inf = 0 není
  detektor neurčitelnosti (výsledek 3) a kombinace "odhad R_inf + volné L"
  nebyla měřená (D4).
- Nová otevřená otázka: neurčitelné případy teď padají na HF medián, který je
  u otevřeného oblouku sám zkreslený (F8, viz výsledek 3).

---

## Závěr

1. **L do modelu DRT zařadit: ano.** Dnešní DRT neumí vytvořit Im(Z) > 0,
   takže indukční data fituje špatně i s *přesně správným* R_inf (chyba
   rekonstrukce 5-18 %). S volným L klesne na 0.2-3.8 % a L vyjde na 1-2 %
   přesně. Riziko: volné L pohltí i jiné chyby modelu, proto je nutné ho
   vracet a hlídat.
2. **R_inf jako proměnná: jen jako diagnostika.** Když se λ vybere znovu
   (GCV) na *rozšířeném* systému, je společný odhad lepší než medián u tří
   z pěti syntetických případů. Proti `--ri-fit` z 0.40.0 prohrává ve všech
   určitelných případech (C -2.9 % proti -0.5 %, C2 +15.7 % proti -1.0 %).
   Na reálných datech s nesouladem modelu padá na mez R_inf = 0. Není to
   náhrada externího odhadu. Je to druhý, nezávislý odhad a jeho rozdíl proti
   externímu může sloužit jako diagnostika.
3. **Identifikovatelnost se tím nevyřeší.** Oblouk nad f_max nebo silně
   neuzavřený CPE oblouk neurčí žádná varianta (+940 % až +1400 %).
4. Návrh v `DRT_IMPROVEMENTS.md` bod 2 ("~10 řádků, žádná změna") situaci
   podceňuje. Bez nového výběru λ je bias -20 % (případ C), chybí ošetření
   meze a nevážený fit dělá R_inf u dat s velkým dynamickým rozsahem
   prakticky neviditelným.

---

## Současný stav

Stav při studii (0.37.0), `drt/linear_system.py:52-90`, `drt/core.py`:

```
b = [Re Z - R_inf, Im Z]         R_inf pevné (medián / --ri-fit / preset)
A = [A_re; A_im]                 jen gamma, bez L
tau in [1/w_max, 1/w_min]        bez rozšíření
min ||A g - b||^2 + lambda ||D2 g||^2,  g >= 0,  bez vážení
```

Chyba rekonstrukce: `mean(|Z - Z_rec| / |Z|)` (dnes `core.py:364`, práh
varování 10 % na `:367`).

Změny od té doby (0.40.0):
- **Vážení:** výchozí `weighting='sqrt'` (1/sqrt|Z|, od 0.38). Studie počítala
  bez vážení, tedy dnešní `'uniform'`.
- **tau-mřížka:** volitelné rozšíření za *pomalý* konec (`tau_extend_decades`,
  od 0.39). Rychlý konec se záměrně nerozšiřuje, v souladu s výsledkem 4.
- **R_inf:** stále pevné. `calculate_drt` bere buď `r_inf_preset`, nebo HF
  medián (`rinf_estimation.hf_median`). Parametr `use_rl_fit` je od 0.40
  odstraněný. `--ri-fit` odhad počítá předem a předá ho jako preset.
- **L:** v DRT stále chybí.

---

## Prototyp

Rozšířený systém, x = [R_inf, L', gamma]:

```
A_re_ext = [1, 0,         A_re]
A_im_ext = [0, w/w_max,   A_im]          L = L' / w_max  (škálování kvůli podmíněnosti)
D_ext    = [0, 0,         D2]            R_inf a L se neregularizují
nnls([W*A_ext; sqrt(lambda)*D_ext], [W*b; 0])     R_inf >= 0, L >= 0
```

- λ: buď převzaté z dnešního běhu s přesným R_inf ("λ fix"), nebo znovu
  vybrané `find_optimal_lambda_hybrid` na rozšířeném systému ("GCV").
- W: jednotkové (jako dnes) nebo 1/|Z| normalizované na průměr 1.
- Varianta "správné R + L": R_inf pevně na skutečné hodnotě, L volné. Odděluje
  vliv L od vlivu R_inf.
- Ověřeno: bez sloupce L a s pevným R_inf dává prototyp přesně stejnou chybu
  rekonstrukce jako `calculate_drt` (4.46 / 17.58 / 14.65 / 0.07 %).

Syntetická spektra (Rs = 10 Ohm, 10 bodů/dekádu). Značení navazuje na
`archive/AUDIT_ri_fit_2026-09-25.md`, ale **flat a A1 jsou tu jiná spektra** než
v auditu (flat má navíc ZARC a končí na 100 kHz, A1 má dva oblouky):

| Případ | Model | Rozsah |
|---|---|---|
| flat | Rs + RC(100, 10 ms) + ZARC(50, 100 us, 0.9) | 0.1 Hz - 100 kHz |
| A2 | Rs + jw 1 uH + RC(100, 100 us) | 0.1 Hz - 1 MHz |
| C | Rs + jw 10 uH + RC(100, 10 us) | 0.1 Hz - 1 MHz |
| C2 | jako C, ZARC n = 0.7 | 0.1 Hz - 1 MHz |
| B | Rs + ZARC(100, 100 us, 0.8) | 0.1 Hz - 100 kHz |
| D | Rs = 0.05 + ZARC(100, 1 ms, 0.6) | 0.1 Hz - 100 kHz |
| A1 | Rs + jw 0.1 uH + RC(100, 10 ns) + RC(100, 1 ms), první oblouk nad f_max | 0.1 Hz - 1 MHz |

Šum: 1 % komplexní gaussovský, 20 seedů (`default_rng(0..19)`).

---

## Výsledky

### 1. Vliv L (R_inf pevně na správné hodnotě)

Chyba rekonstrukce v %, šum 1 %, průměr přes 20 seedů:

| Případ | dnes (bez L) | + volné L | nalezené L | skutečné L |
|---|---|---|---|---|
| flat | 0.89 | 0.88 | 62 nH | 0 |
| A2 | **5.00** | 1.16 | 1006 nH | 1000 nH |
| C | **17.89** | 3.76 | 10168 nH | 10000 nH |
| C2 | **15.06** | 0.90 | 10119 nH | 10000 nH |
| B | 1.07 | 1.04 | 284 nH | 0 |

Dnešní DRT má na indukčních datech chybu rekonstrukce 5-18 % a u C a C2
tím překračuje vlastní varovný práh 10 % (`core.py:271`), a to s dokonalým
R_inf. NNLS s gamma >= 0 vytváří jen Im < 0. Indukční část dat proto nemá
kam dát a deformuje gamma.

**Fiktivní L:** na datech bez indukčnosti vyjde volné L 42-330 nH (šum,
chyba okraje tau-mřížky). U `real_gamry_example.DTA` vyšlo společně s R_inf = 0
L = 8.8 uH (viz výsledek 3). L tedy kromě indukčnosti pohlcuje i chyby modelu.

### 2. R_inf jako proměnná, syntetická data

Chyba R_inf v % (průměr / směrodatná odchylka), šum 1 %, 20 seedů:

| Případ | medián (dnes) | `--ri-fit` | společně, λ fix | **společně, GCV** | společně, GCV, 1/\|Z\| |
|---|---|---|---|---|---|
| flat | +2.8 | +0.1 | +0.2 / 0.5 | **+0.2 / 0.5** | +0.1 / 0.7 |
| A2 | -0.2 | -0.2 | -2.2 / 0.5 | **-1.5 / 0.4** | -0.4 / 0.6 |
| C | +0.9 | +10.2 | **-20.1** / 0.5 | **-2.9 / 0.7** | -0.7 / 0.7 |
| C2 | +40.6 | +88.2 | +18.1 / 3.7 | **+15.7 / 4.1** | +12.9 / 3.1 |
| B | +19.0 | +0.8 | +0.4 / 1.2 | **-0.3 / 1.4** | +4.5 / 0.6 |

(Sloupce medián a `--ri-fit` jsou z jednoho seedu, ostatní průměr přes 20.
Rozdíly v jednotkách procent mezi nimi proto nejsou průkazné. Sloupec
`--ri-fit` je estimátor z 0.37. Estimátor z 0.40.0 dává na
spektrech auditu (5 seedů) flat -0.7, A2 -0.1, C -0.5, C2 -1.0 a B -0.5 %.)

- **Výběr λ je rozhodující.** S λ převzatým z běhu s pevným R_inf vychází C
  na -20 %: regularizace, která penalizuje gamma a R_inf ne, přesune váhu
  z rozmazaného píku do R_inf. S λ vybraným na rozšířeném systému -2.9 %.
- Proti mediánu je společný odhad lepší u flat, B a C2 a mírně horší u A2 a C.
  Proti `--ri-fit` z 0.37 byl lepší u C (-2.9 proti +10.2) a C2 (+16 proti +88).
  Proti `--ri-fit` z 0.40.0 (C -0.5, C2 -1.0) je horší.
- Vážení 1/|Z| pomáhá u indukčních případů (C -0.7 %), ale zhoršuje B (+4.5 %).
- Oproti modelu R-L-ZARC z `AUDIT_ri_fit` (C2: 0.0 %) je společný odhad
  horší. Rozmazání gamma regularizací se projeví i v R_inf.

### 3. R_inf jako proměnná, neurčitelná a reálná data

| Data | medián | `--ri-fit` | společně (GCV) | společně (GCV, 1/\|Z\|) | chyba rekonstrukce: dnes / společně / vážené |
|---|---|---|---|---|---|
| D (Rs = 0.05) | +3289 % | +27 % | +1032 % | +1018 % | - |
| A1 (oblouk nad f_max) | +998 % | +996 % | +942 % | +942 % | - |
| example_eis_data.csv (Rs = 10) | +3262 % | +258 % | +1356 % | +1384 % | - / 1.28 / 1.38 |
| EISPOT-test1.DTA | 1.446 | 1.180 | **0 (mez)** | 1.231 | 46.9 / 51.1 / **8.6** |
| real_gamry_example.DTA | 1402 | 825.9 | **0 (mez)**, L = 8.8 uH | 0 (mez) | 7.8 / 4.2 / 3.6 |

- **EISPOT-test1.DTA** má blokující nízkofrekvenční chování (Im = -1e6 Ohm
  při 3 mHz). To DRT bez sériové kapacity nepopíše a dnešní chyba rekonstrukce
  je 47 %. Při neváženém fitu je R_inf (~1 Ohm) proti datům řádu 1e6 Ohm
  v nejmenších čtvercích neviditelné a NNLS ho nastaví na nulu. S vážením
  1/|Z| vychází 1.231 Ohm a chyba rekonstrukce klesne na 8.6 %.
- **real_gamry_example.DTA** nemá uzavřený oblouk ani na f_max (fáze -59 deg).
  `AUDIT_ri_fit` ho označil za neurčitelný. Společný fit dává R_inf = 0
  a lepší rekonstrukci (4.2 % proti 7.8 %). Dnešní medián 1402 Ohm přitom leží
  **nad** Re(Z) na f_max (826 Ohm), protože bere pět bodů (2087, 1740, 1402,
  1100, 826 Ohm), z nichž nižší frekvence mají vyšší Re. Jde o konkrétní
  příklad k F8. Od 0.40.0 na tento medián padá i `--ri-fit`, protože soubor
  správně označí jako neurčitelný. Zkreslený fallback je tak otevřená otázka.
- Uvíznutí R_inf na mezi 0 je příznak, že R_inf z dat neurčíme nebo že model
  DRT data nepopisuje, a nemá se vracet jako výsledek. **Není to ale detektor
  neurčitelnosti:** skutečně neurčitelné případy D, A1 a CSV na mez nespadnou
  a vracejí +942 % až +1384 %. Mez je postačující příznak, ne nutný.
  Detektor neurčitelnosti externího odhadu je stderr R_s v `--ri-fit` (0.40.0).

### 4. Rozšíření tau-mřížky (F9) a volné R_inf se nesnášejí

Prvek s w tau << 1 je čistě reálný, a tedy k nerozeznání od R_inf. Při
rozšíření mřížky o dekádu nad 1/w_max (λ fix) vycházely nestabilní výsledky:
B se šumem -97 %, C2 -30 % místo +19 %, A1 +559 % místo +943 %. Pokud bude
R_inf volné, horní (krátký) konec mřížky by se rozšiřovat neměl. Rozšíření
dolního konce (F9, pomalé procesy) tím dotčeno není.

---

## Proč to tak vychází

- **Kramers-Kronig:** Im(Z) R_inf neobsahuje. Gamma se určuje hlavně z Im
  a R_inf je to, co zbude v Re. Když celé rozdělení leží v tau-mřížce (flat,
  B), funguje to dobře. Když jeho část leží mimo mřížku na krátkém konci (A1,
  CSV, D), připadne R_inf a žádná metoda to nepozná.
- **Regularizace je nesymetrická:** penalizuje gamma, ne R_inf ani L. Při
  větším λ se "levnější" vysvětlení přesouvá do nepenalizovaných parametrů.
  Proto se λ musí vybírat na rozšířeném systému.
- **Vážení:** nevážený fit vidí jen body s velkým |Z|. R_inf je malé
  a ovlivňuje právě body s malým |Z|, takže při dynamickém rozsahu 1e6
  zmizí.

---

## Návrhy

### Výsledky k D1 (0.41.0)

Měřeno před implementací, sloupec L přidaný k matici z `_build_drt_matrices`
(vážení `sqrt`, auto-λ vybrané na rozšířeném systému), šum 1 %, 3 seedy:

- **Fiktivní L neškodí.** Na datech bez indukčnosti (flat, B, D) vychází
  L = 0-600 nH, ale R_pol se mění o méně než 0.1 %, píky se posunou nejvýš
  o 0.06 dekády a chyba rekonstrukce se nezhorší. Na reálných souborech
  (`EISPOT-test1.DTA`, `real_gamry_example.DTA`, `example_eis_data.csv`)
  vychází L = 0.
- **Kritérium pro `'auto'`:** všechny indukční případy (A2, C, C2) mají
  Im(Z) > 0 v horní dekádě, neindukční žádný. `EISPOT-test1.DTA`, kde
  `--ri-fit` najde L = 326 nH maskované obloukem (Im(f_max) < 0), dává
  i s vynuceným L hodnotu L = 0. Zvoleno "aspoň jeden bod horní dekády
  s Im > 0" (`DRT_INDUCTANCE_DECADES = 1`).
- **Tvrzení D4 změřeno:** `--ri-fit` (0.40.0) + volné L dává na A2, C a C2
  chybu rekonstrukce 1.1-1.4 % a R_pol do 0.8 %, tedy stejně dobře jako
  s přesným R_inf.

**D1 - L do DRT [HOTOVO v 0.41.0].** Nepenalizovaný sloupec jw/w_max, L >= 0.
Stejný mechanismus už používá Lin-KK (`estimate_R_linear(include_L=True)`).
Podle paměti autora mají DRTtools (Wan, Saccoccio, Chen, Ciucci 2015)
a pyDRTtools u indukčnosti volby "bez L / s L / zahodit indukční body".
Před zavedením ověřit v literatuře.
- **Před zavedením změřit dopad fiktivního L** na gamma a R_pol u dat bez
  indukčnosti (B: L = 284 nH). Studie ho vyčíslila, ale jeho vliv neměřila.
  Podle toho rozhodnout, zda L zapínat vždy, nebo jen při Im > 0 na HF konci.
- L vracet v `DRTResult` a vypisovat.
- Varovat, když L > 0 a data na HF konci nemají Im > 0 (fiktivní L). Jinou
  možností je zapínat L automaticky jen tehdy, když indukční body existují.
- Chybu rekonstrukce počítat včetně jwL.

**D2 - Společný odhad R_inf jen jako diagnostika.** Jako zdroj R_inf nemá
smysl, protože `--ri-fit` z 0.40.0 je na určitelných případech přesnější.
Pokud se zavede:
- λ vybírat na rozšířeném systému (GCV/L-křivka na [1, 0, A]),
- R_inf na mezi 0 (aktivní omezení NNLS) hlásit jako varování, že model DRT
  data nepopisuje nebo R_inf nelze určit. Nevracet ho a použít externí odhad.
  Nespoléhat na mez jako na detektor neurčitelnosti (výsledek 3),
- při externím R_inf počítat společný odhad paralelně a varovat, když se liší
  víc než o X % (práh určit na reálných datech, synteticky se dobré případy
  liší o <= 3 %),
- horní konec tau-mřížky nerozšiřovat (výsledek 4).

**D3 - Vážení DRT fitu (zčásti vyřešeno v 0.38).** 1/|Z| pomohlo C
(-2.9 -> -0.7 %) i EISPOT (rekonstrukce 47 -> 8.6 %), ale zhoršilo B
(-0.3 -> +4.5 %). Studie vážení proběhla v 0.38 a za výchozí zvolila `sqrt`
(1/sqrt|Z|), ne 1/|Z|. Jak se se `sqrt` chová argument o "neviditelném" R_inf
u EISPOT, se neměřilo. Je to relevantní jen pro D2.

**D4 - Externí R_inf zlepšit podle `archive/AUDIT_ri_fit_2026-09-25.md` N1/N2.
[HOTOVO v 0.40.0]** R-L-(R|Q) s kontrolou identifikovatelnosti přes stderr R_s
je hlavní cestou pro R_inf, D2 je nanejvýš diagnostika. Původní tvrzení, že
kombinace "dobré externí R_inf + volné L" byla nejlepší ve všech určitelných
případech, měřením podložené nebylo. Tabulka 1 počítá s *přesným* R_inf,
ne s odhadem. Kombinaci "`--ri-fit` + volné L" je třeba změřit v rámci D1.

**D5 - Aktualizovat `DRT_IMPROVEMENTS.md` bod 2 [HOTOVO v 0.41.0]** odkazem na tento
dokument, aby nesliboval "minimální změnu kódu" a "teoreticky přesnější" bez
výhrad. Jeho "současný stav" (odhad R_inf z jednoho bodu) neplatí ani pro
HF medián.

**D6 - Fallback pro neurčitelné R_inf [HOTOVO v 0.41.0: Re(Z) na f_max].**
Změřeno (šum 1 %, 20 seedů, jen seedy s příznakem): D +3295 -> +2483 %,
A1 beze změny (+997 %), `real_gamry` 1402 -> 826 Ohm. Na plochém, čistě R-L
a Warburgově konci, kde fit příznak dostane jen občas, je cena -0.6 až -0.7 %
místo ~0 %. Varianta min(medián, Re(f_max)) vyšla stejně, zvolena jednodušší.
Výchozí medián DRT (F8) se nemění. Původní text: HF medián z pěti
bodů je u otevřeného oblouku zkreslený nahoru (`real_gamry_example.DTA`:
1402 Ohm proti 826 Ohm na f_max). Rozhodnout, zda fallback nahradit, například
hodnotou Re(Z) na f_max, která je u kapacitního konce horní mezí R_inf.

---

## Omezení studie

- Syntetická data jsou z modelů RC/ZARC. Reálné spektrum se skutečně známým
  R_inf a indukčností v sadě chybí.
- λ pro "společně, GCV" vybírá `find_optimal_lambda_hybrid` s dnešním
  rozsahem (1e-5, 1). U vážené varianty vybírá GCV λ pro jinak škálovaný
  systém. Robustnost výběru λ napříč daty ověřena nebyla.
- Vliv na tvar gamma(tau) a na detekci píků (scipy/GMM) se neměřil, jen R_inf,
  L, R_pol a chyba rekonstrukce.
- Sloupce medián a `--ri-fit` v tabulce 2 jsou z jednoho seedu, ostatní
  z dvaceti.
- Prototyp není v repozitáři. Popis systému stačí k jeho rekonstrukci, čísla
  ale nejde přímo znovu spustit. Revize k 0.40.0 přeměřila jen tabulku 1
  ("dnes") a hodnoty `--ri-fit`, a to přes veřejné API.
