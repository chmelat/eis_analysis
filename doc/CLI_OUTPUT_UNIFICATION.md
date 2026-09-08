# Sjednocení CLI výstupu — plán refaktoringu

Datum: 2026-09-08
Rozsah: pět modulů skupiny C (viz níže) + `format_voigt_report()`.
Rozhodnuto: **zpětná kompatibilita není omezením, prioritu má čistota designu.**

Navazuje na analýzu z 2026-09-01. Ta pojmenovala problém, tento dokument
říká, co se s ním udělá a v jakém pořadí.

---

## 1. Cíl

Jedna konvence pro veškerý uživatelský výstup:

> Modul počítá a vrací dataclass. Tiskne výhradně vrstva CLI.

Deklarovaná je už dnes v docstringu `drt/results.py` a dodržují ji DRT
(`DRTDiagnostics`), Kramers-Kronig (`KKResult`) a Z-HIT (`ZHITResult`).
Vzor příjemce je `cli/handlers/drt.py` — přečte výsledek, vytiskne sekci,
projde `warnings` a vypíše je (řádky 59, 187, 244).

Konvence B (`format_voigt_report()` vrací string) se do A složí: funkce je
čistá prezentace a patří do `cli/handlers/`.

## 2. Co je vlastně dluh

Původní počty logovacích volání míchaly tři různé věci. Rozpad podle úrovně:

| Modul | info | warning | debug |
|---|---:|---:|---:|
| `analysis/oxide.py` | 43 | 24 | 2 (+2 error) |
| `fitting/voigt_chain/fitting.py` | 33 | 2 | 6 |
| `io/data_loading.py` | 20 | 10 | 6 |
| `fitting/auto_suggest.py` | 14 | 6 | 0 |
| `fitting/voigt_chain/mu_optimization.py` | 14 | 6 | 0 |
| **celkem** | **124** | **48** | **14** |

- **`logger.info` (124) = vlastní dluh.** Je to UI: hlavičky sekcí, tabulky
  výsledků, průběh výpočtu po krocích. Tohle se celé stěhuje do CLI.
- **`logger.debug` (14) zůstává na místě.** Trasování pro vývojáře, ne
  uživatelský výstup; `-v` je pro něj správný vypínač a `data_loading`
  ho používá korektně (delimiter, indexy sloupců, přeskočené řádky).
- **`logger.warning` (48) se dělí podle kritéria R1.** Výhrady k datům nebo
  k výsledku jdou do `warnings: List[str]` na výsledku — tak to dělá
  `DRTDiagnostics` a tak to CLI už umí vytisknout. Provozní události
  (chyba I/O, nepovedený parse) zůstávají v logu. Viz R1 níže.

## 3. Pravidla cílové konvence

1. Modul nemá `logger.info`. Tečka.
2. Diagnostika, kterou dnes modul tiskne, je pole výsledné dataclass —
   včetně údajů o *postupu* (zvolená větev, počty, prahy), ne jen finálních
   čísel. Bez toho se sekce v CLI nedá složit.
3. Varování, které kvalifikuje nějaký výsledek, jde do jeho
   `warnings: List[str]`. Varování o operaci, která žádný výsledek
   nevyprodukovala, zůstává `logger.warning`. Kritérium a hranice: R1.
4. `logger.debug` je povolen kdekoliv.
5. Funkce vrací dataclass, ne tuple ani dict.
6. Oddělovač: jednotně `log_separator()` (50 znaků). Inline `"=" * 60`
   zmizí spolu s `logger.info`, takže se to vyřeší samo — jen se hlídá,
   aby nové řádky v handlerech psaly `log_separator()`.

## 4. Pořadí

Každá etapa je samostatně vydatelná (zelené testy, vlastní commit, vlastní
řádek v CHANGELOG). Pořadí je podle poměru přínos/cena, ne podle velikosti.

| # | Modul | info | Hlavní cena | Odhad |
|---|---|---:|---|---|
| 1 | `io/data_loading.py` | 20 | `log_metadata()` je printer ve veřejném API | S |
| 2 | `fitting/voigt_chain/mu_optimization.py` | 14 | volá ho i KK — knihovní volající | S |
| 3 | `fitting/auto_suggest.py` | 14 | dict -> dataclass, přesun `format_voigt_report` | M |
| 4 | `fitting/voigt_chain/fitting.py` | 33 | tuple -> dataclass, 4 kroky postupu | M |
| 5 | `analysis/oxide.py` | 43 | 96 caplog assertů v `tests/test_oxide.py` | L |

### Etapa 1 — `io/data_loading.py` — HOTOVO

Tři odchylky od původního návrhu, zjištěné až nad kódem:

- Dataclass se jmenuje `LoadResult`, ne `LoadedData` — ten název už patří
  kontejneru CLI v `cli/utils.py` (drží navíc `title` a `ocv_data` a pokrývá
  i syntetická data, kde žádný soubor není). Nese `metadata`, protože
  `load_data()` je stejně parsuje kvůli kontrole useknutého sweepu; CLI tím
  přestalo číst soubor podruhé.
- Převeden i `read_gamry_native()`, který plán nezmiňoval. Bez toho by
  výhrada "hlavička nepojmenovává sloupce, beru standardní pořadí" neměla
  kam jít a kritérium R1 by v tomto modulu neplatilo.
- `print_metadata()` si v CLI ponechává `"=" * 60`, aby seděl se sousedním
  blokem syntetických dat. Sjednocení šířky oddělovače (pravidlo 6) je
  samostatný průchod na konci; udělat ho teď by výstup rozhodilo uprostřed
  refaktoringu.

Jediná změna výstupu: výhrada o useknutém sweepu je jeden řádek místo dvou
(stejná informace). Ověřeno diffem výstupu CLI proti HEAD na `.DTA`,
useknutém `.DTA` a `.csv` s nerozpoznanými sloupci. Syntetická data nejsou
seedovaná, takže se u nich porovnávala jen struktura výstupu.

Původní zadání: `load_data()` a `load_csv_data()` tisknou "Loaded N points" a rozsah
frekvencí; `log_metadata(metadata)` je čistý printer, který žije v `io/`
a je exportován ve veřejném API (`eis_analysis/__init__.py:29`).

- `load_data()` / `load_csv_data()` -> vrací `LoadedData(frequencies, Z,
  filename, warnings)` místo tuple.
- `log_metadata()` se přesouvá do `cli/data_handling.py` jako
  `print_metadata()` a mizí z `io/__init__.py` i z kořenového `__all__`.
- Validační varování (malý rozsah frekvencí, duplicity, useknutý sweep)
  jdou do `warnings`. Modul je zároveň jediný, kde žijí obě kategorie R1:
  `:123` (vadná ZCURVE hlavička), `:196` a `:384` ("Error reading file"
  v `except` větvi, která vrací `None`) zůstávají `logger.warning`.
- Volající: `cli/data_handling.py:107-109`.

Verifikace: `python3 -m pytest tests/test_io_data_loading.py` (15 caplog
míst přepsat na `result.warnings`), `python3 eis.py example/*.DTA` dává
stejný výstup jako před změnou.

### Etapa 2 — `fitting/voigt_chain/mu_optimization.py`

Nejsilnější argument z celého seznamu: `find_optimal_M_mu()` volá
`validation/kramers_kronig.py:411`, tedy modul, který sám konvenci A
dodržuje. KK vrací čistý `KKResult` a přitom pod ním teče 14 řádků
cizího `info` výstupu.

- Návrat 6-tuple `(M, mu, tau, elements, L_value, C_value)` -> `MuOptimization`
  s týmiž poli plus `iterations: List[MuIteration]` (M, mu na iteraci — dnes
  řádek 208) a `warnings`.
- Volající: `voigt_chain/fitting.py:403` a `kramers_kronig.py:411`.
- Tisk iterační tabulky přebírá `cli/handlers/fitting.py`; KK ji netiskne
  (dnes ji tiskne, aniž by o tom KK handler věděl — to je ta chyba).

Verifikace: `python3 -m pytest tests/test_kramers_kronig.py tests/test_voigt_chain.py`;
`python3 eis.py example/*.DTA --validate` nesmí vypsat žádný mu řádek mimo
sekci fittingu.

### Etapa 3 — `fitting/auto_suggest.py`

Jediný modul, který dělá B i C zároveň.

- `analyze_voigt_elements()` -> vrací `VoigtSuggestion` dataclass místo
  `dict` (pole: `R_inf`, `R_pol`, `peak_method`, `elements: List[VoigtElement]`,
  `quality`, `warnings`).
- `format_voigt_report()` -> `cli/handlers/drt.py` jako `_print_voigt_report()`,
  mizí z `fitting/__init__.py:87,143` i z kořenového API.
- Sekce "Automatic circuit suggestion from DRT" (řádky 108-110) se přesune
  do handleru; dnes se vytiskne dřív, než handler stihne svou hlavičku.

Verifikace: `python3 -m pytest tests/test_cli_integration.py -k voigt`.

### Etapa 4 — `fitting/voigt_chain/fitting.py`

33 `info` řádků je narace postupu "Step 1..4" (mu optimalizace, tau grid,
regrese, prořezání, stavba obvodu) plus závěrečné shrnutí.

- `fit_voigt_chain_linear()` -> `VoigtChainFit(circuit, initial_params,
  diagnostics, warnings)`.
- `VoigtChainDiagnostics` nese to, co dnes tečou řádky 390-568: zvolená
  větev (mu vs. pevná mřížka), parametry mřížky, metoda a váhování regrese,
  `R_s`, rozsah `R_i`, reziduum, `L`, prahy prořezání a počty
  před/po, finální počet parametrů.
- Volající: `cli/handlers/fitting.py:416`, `tests/test_voigt_chain.py` (6 míst),
  `tests/test_cli_integration.py:470`.

Verifikace: `python3 -m pytest tests/test_voigt_chain.py tests/test_cli_integration.py`;
diff výstupu `python3 eis.py example/*.DTA --voigt` proti uložené referenci.

### Etapa 5 — `analysis/oxide.py`

Nejdražší, dělat naposledy. `cli/handlers/oxide.py:69,76` dnes návratovou
hodnotu **vůbec nezachytí** — `OxideAnalysisResult` existuje, ale veškerý
uživatelský výstup vzniká uvnitř modulu. Je to nejčistší případ konvence C
v repu.

- `OxideAnalysisResult` musí povyrůst o to, co se dnes jen tiskne:
  `candidates: List[dict]` (výpis nalezených kapacitních prvků, řádky
  485-500), `selection_reason: str` (proč zvítězil dominantní prvek, řádky
  323-367), `mode` ('circuit' | 'hf_estimate', řádek 660), `cc_regime`
  (režim C*(omega) z `_cc_capacitance_regime`) a `warnings`.
- `_log_cc_capacitance_choice()` -> vrací text místo logování.
- Handler `cli/handlers/oxide.py` (dnes 82 řádků) vyroste o tisk sekce,
  tj. zhruba to, co dnes dělají řádky 601-646 a 785-908 modulu.
- **96 caplog assertů v `tests/test_oxide.py`** se přepisuje na kontrolu
  polí výsledku. Většina je mechanická; testy typu
  `test_estimate_permittivity_does_not_log_thickness` (řádek 60) se převedou
  na `assert result.thickness_nm is None`, což je i významově správnější —
  dnes testují formátování, ne chování.

Verifikace: `python3 -m pytest tests/test_oxide.py`; ruční diff sekce
"Oxide layer analysis" před/po.

## 5. Dopad na testy

146 odkazů na `caplog` v šesti souborech. Dotčené etapami:

| Soubor | caplog | Etapa |
|---|---:|---|
| `tests/test_oxide.py` | 96 | 5 |
| `tests/test_cli_integration.py` | 21 | 3, 4 |
| `tests/test_io_data_loading.py` | 15 | 1 |
| `tests/test_kramers_kronig.py` | 4 | 2 |
| `tests/test_outliers.py` | 7 | — |
| `tests/test_residual_diagnostics.py` | 3 | — |

Po refaktoringu má `caplog` zůstat jen tam, kde se testuje CLI vrstva nebo
`logger.warning`/`debug`. Assert na `logger.info` v testu knihovního modulu
je po dokončení etapy chyba.

Vedlejší přínos: testy přestanou být závislé na formátování řetězců.

## 6. Rozhodnutí k potvrzení

**R1 — kam s `logger.warning` (48 volání).** Rozhoduje se, jestli je varování
*událost v běhu programu*, nebo *součást výsledku*. Obojí je obhajitelné:
knihovna, která emituje log record a nechá aplikaci rozhodnout, co s ním, je
standardní pythonovská praxe; jenže hlášky v tomto repu tak nevypadají.

```
oxide.py:479  "No capacitive element (C, Q, K, CC) found in circuit"
oxide.py:480  "Falling back to high-frequency estimate..."
oxide.py:706  "Positive imaginary impedance (inductive) - result may be invalid"
```

To nejsou události, to jsou vědecké výhrady k vrácenému číslu. Když
`analyze_oxide_layer()` vrátí `thickness_nm = 42.0`, údaj "vzniklo to
z nouzového odhadu a může být neplatné" je část odpovědi, ne poznámka
na okraj.

**Kritérium: existuje výsledek, který to varování kvalifikuje?**

- **Ano -> do `warnings` na tom výsledku.** Prakticky celý `oxide.py`
  (24 z 48), malý rozsah frekvencí, duplicitní frekvence, useknutý sweep,
  záporná `R_i`. Tyto výhrady se tím stanou daty: dají se serializovat do
  protokolu, testovat bez `caplog` a CLI je vytiskne uvnitř své sekce,
  nikoliv v okamžiku vzniku.
- **Ne -> zůstává `logger.warning`.** Chyba I/O, nepovedený parse, funkce
  vracející `None` (`data_loading.py:123, :196, :384`). Není kam to dát;
  žádný výsledek nevzniká.

Co tím padá: dnes se varování z `data_loading` vysypou na výstup v okamžiku
vzniku, tedy klidně dřív, než CLI vytiskne hlavičku sekce, do které patří —
táž chyba pořadí, jaká se řešila v v0.28.1. Kdyby se stěhovalo jen `info`,
u varování by zůstala neopravená.

Zbývající riziko: varování, které se stalo daty, tiše zmizí, když je volající
nevytiskne. To je horší druh regrese než hlučný `logger.warning`, protože se
neprojeví. Proto se v každé etapě kontroluje výstup CLI proti stavu před
změnou, ne jen zelené testy.

**R2 — jak daleko s veřejným API.** Etapy 1 a 3 ruší tři exportované názvy
(`log_metadata`, `format_voigt_report` a tuple návraty). Zpětná
kompatibilita není omezením, takže se ruší bez deprecation shimu; verze
skočí na 0.32.0 a CHANGELOG dostane sekci **Breaking changes**.

**R3 — kdy vydat.** Buď pět patch verzí (0.31.4 ... 0.32.0), nebo jedna
0.32.0 na konci. Doporučení: **jedna 0.32.0**, protože etapy 1 a 3 samy
o sobě rozbíjejí API a vydávat rozbité API pětkrát po sobě nemá komu
posloužit. Commity zůstávají po etapách.
