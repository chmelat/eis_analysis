# Plan: meze parametru podle plochy elektrody

Stav: **zamitnuto (2026-10-08), neimplementovano**. Navrh z 2026-10-07,
vychozi verze eis_analysis v0.57.0. Duvody v sekci Rozhodnuti; zbytek
dokumentu je puvodni navrh.

## Rozhodnuti

1. **Plocha je nespolehlivy vstup.** Operatori ji casto nezaznamenavaji,
   takze v metadatech zustava vychozi hodnota bez vztahu ke skutecne plose.
   Z DTA nelze poznat, zda plochu nekdo nastavil (Gamry `AREA` zapisuje
   vzdy). Spatne vyplnena plocha by meze tise posunula o tolik dekad, o kolik
   se lisi, a ze spektra se to poznat neda. Zbylo by jen explicitni
   `--area`, coz je pro jediny hypoteticky pripad (vysoke R na male plose)
   prilis.
2. **Stress test plan nepodporuje.** Invariant Babs testuje zmenu jednotek
   (Z x k), ne geometrii. Plocha zavislost na jednotkach neresi; navrzeny
   invariant "Z x k spolu s A / k" by prosel jen konstrukci.
3. **Alternativa, pokud bude potreba:** az se najde realne spektrum, jehoz
   fit konci na horni mezi R (1e10 Ohm), rozsirit `RESISTANCE_RANGE` na mez
   pouziteho potenciostatu, bez plochy. Paralelni R je identifikovatelne jen
   pri namerenem |Z| ~ R, takze mez vadi jen spektrum s |Z| > 1e10. Cena:
   DE prohledava R v log prostoru o 2 dekady sirsi (14 -> 16).

## Kontext

Meze fitu (`PARAMETER_BOUNDS`, `fitting/bounds.py`) jsou absolutni a
fyzikalni: R 1e-4..1e10 Ohm, C 1e-15..0.1 F atd. Popisuji ale cely vzorek,
ne material. Stejny material na elektrode 1 dm2 a 1 mm2 (pomer ploch 1e4)
ma R 1e4x vetsi a C, Q 1e4x mensi; n, tau a plosne normovane hodnoty
(R*A v Ohm*cm2, C/A v F/cm2) jsou stejne. Pevne meze v Ohm tak nemohou
pokryt vsechny velikosti vzorku: kvalitni oxid s 1e9 Ohm*cm2 ma na 1 mm2
(0.01 cm2) R = 1e11 Ohm, nad horni mezi 1e10. Proto se meze v historii
opakovane rozsirovaly (komentare u sigma a R_W v `bounds.py`).

Stress test (`doc/STRESS_TEST.md`, invariant Babs) ukazuje druhy dusledek:
spektrum x k (fyzikalne vzorek s plochou A/k) dava jiny fit, protoze je
jinak daleko od mezi - i kdyz na ne nenarazi (trust region scipy skaluje
krok podle vzdalenosti k mezim).

Meze odvozene ze spektra byly zamitnuty (2026-10-07): v jednom spektru muze
byt R_ser blizko nule i R oxidu blizko nekonecnu a mez z |Z| by mohla
oriznout fyzikalne smysluplne optimum. Plocha je jina vec: neodvozuje se ze
spektra, je to znama geometrie vzorku.

**Cil:** fyzikalni meze vztazene na jednotku plochy, prepoctene plochou
elektrody na meze konkretniho vzorku. Bez zadane plochy beze zmeny.

## Navrh

Dnesni `PARAMETER_BOUNDS` se cte jako meze pro referencni plochu
A_ref = 1 cm2. Pro plochu A se mez parametru nasobi (A / A_ref)^p, kde p je
mocnina plochy daneho typu parametru:

| p | Parametry | Proc |
|---|---|---|
| -1 | R, R_W, sigma (W), sigma_GE, A_DQ | impedance, Z ~ 1/A |
| +1 | C, Q, G, C_inf, dC, C_YG | admitance / kapacita, ~ A |
| 0 | n, tau, tau_W, tau_GE, tau_CC, alpha_CC, n_DQ, tau_DQ, U_DQ, p_YG, tau_YG | bezrozmerne nebo cas, na plose nezavisi |
| 0 | L | indukcnost kabelu a privodu, ne elektrody |

Pri A = 1 cm2 (vychozi) jsou meze presne dnesni, tedy zadna zmena chovani.
Nove meze zustavaji fyzikalni (zadna informace ze spektra), jen se posouvaji
se znamou geometrii.

Mocnina se vede u kazdeho typu primo v `PARAMETER_BOUNDS` (nebo v tabulce
vedle), aby novy prvek nemohl zapomenout ji uvest; test overi, ze ji ma
kazdy klic.

### API

- `generate_simple_bounds(param_labels, area_cm2=1.0)`.
- `area_cm2: float = 1.0` v `fit_equivalent_circuit`, `fit_circuit_diffevo`
  a `fit_circuit_multistart`; predava se do `generate_simple_bounds` (DE
  a multistart si meze generuji samy, `bounds=` z v0.57.0 to nepokryva).
  Explicitni `bounds=` u `fit_equivalent_circuit` ma prednost.
- Validace: `area_cm2` konecne a > 0, jinak ValueError.

### CLI

`--area` dnes patri do skupiny oxidu a pouziva ho jen `--analyze-oxide`
(vychozi: metadata DTA, jinak 1.0). Navrh: presunout do obecne skupiny,
stejnou hodnotu pouzit pro meze fitu i pro analyzu oxidu, ve vypisu fitu
uvest pouzitou plochu a jeji zdroj.

## Otevrene otazky (rozhodnout pred implementaci)

1. **Plocha z metadat DTA:** pouzit ji pro meze automaticky, nebo jen
   explicitni `--area`? Automaticky = spravne meze bez prace, ale fit
   souboru s plochou != 1 cm2 se tise zmeni oproti v0.57.0.
2. **A_ref = 1 cm2:** odpovidaji dnesni meze opravdu vzorku ~1 cm2? Hodnoty
   v komentarich `bounds.py` (mOhm baterie, Zr oxid EISPOT-M136113-4) byly
   ladeny na realnych vzorcich - zjistit jejich plochy; kdyz se lisi,
   prepocitat meze na 1 cm2, nebo A_ref zvolit jinak.
3. **Odpor elektrolytu R_s** se neskaluje presne s 1/A (u mikroelektrod
   rozlozeni proudu, R ~ 1/r, ne 1/r2). Meze jsou siroke (14 dekad), takze
   to nejspis nevadi; overit na malem vzorku.
4. **Plosne normovana data (Ohm*cm2):** odpovidaji plose 1 cm2 - zadna
   zmena; zminit v dokumentaci.

## Postup

1. Rozhodnout otazky 1-4.
2. Mocniny plochy v `bounds.py`, `generate_simple_bounds(..., area_cm2)`,
   unit test: kazdy klic ma mocninu; pro A = 1 meze beze zmeny; pro
   A = 0.01 R x 100, C / 100, n a tau beze zmeny.
3. `area_cm2` ve trech fit funkcich, regresni test: oxid 1e9 Ohm*cm2 na
   0.01 cm2 (R = 1e11 Ohm) uz nekonci na mezi.
4. CLI `--area` (podle otazky 1), vypis plochy u fitu, test CLI.
5. Stress test: invariant Babs nahradit variantou s plochou skalovanou
   spolu s daty (Z x k, A / k), ktera ma platit presne i s fyzikalnimi
   mezemi; pridat DE, az bude invariant G.
6. README (CLI), `doc/PYTHON_API.md`, docstring `bounds.py`, CHANGELOG,
   bump (nova funkce; pri automatickem pouziti metadat zmena chovani).

## Overeni

- `python3 -m pytest tests/` a `-m slow` (ZScope) beze zmeny pri A = 1.
- Novy invariant v `tests/stress.py` bez selhani mimo tridy tesny/slaby.
- Realny maly vzorek (pokud je k dispozici): fit s `--area` vs bez.
