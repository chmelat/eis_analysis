# Matematický audit modulu DRT

**Datum:** 2026-06-27
**Auditovaná verze:** 0.13.17
**Rozsah:** `eis_analysis/drt/` — `core.py`, `gcv.py`, `peaks.py`,
`term_classification.py` (+ dotčená závislost `rinf_estimation`)
**Typ:** Code-level review matematické korektnosti + empirické ověření
(Claude Code)
**Zaměření:** Správnost použitých řešení (diskretizace DRT, výběr λ, GCV/
L-curve, GMM detekce píků, klasifikace termů). **Není** to audit stylu,
typů ani výkonu — ten řeší `AUDIT_2026-06-23.md`.

---

## 0. Hlavní sdělení (verdikt)

**Jádro DRT je matematicky správné.** Diskretizace Fredholmovy rovnice,
konstrukce jádra `A`, znaménková konvence, Tikhonovova regularizace 2.
derivací a rekonstrukce `Z` jsou implementovány korektně a empiricky ověřeny
(jednopíkový Voigt: poloha τ s chybou 6 %, R_pol 100.04 vs 100.0, rekonstrukční
chyba 0.18 %).

**Slabší místa nejsou v jádře, ale v navazujících heuristikách** — výběru λ
pro NNLS a zejména v **post-processingu píků** (GMM přes replikaci bodů,
double-counting odporů, klasifikace termů podle šířky píku). Tyto části dávají
správný řád, ale jejich kvantitativní výstupy (počet píků, šířky, n u CPE,
dílčí R_i) jsou **závislé na volbě λ a na arbitrárních škálách** a nelze je
brát jako objektivní měření.

**Kritická interakce (empiricky doložena):** auto-λ optimalizuje fit dat, což
pro nízkošumová data tlačí λ na spodní mez rozsahu → DRT degeneruje na úzké
„spiky" (v testu jen **4 nenulové body γ**). Veškerý analýza tvaru píků
(GMM, šířky, klasifikace) pak pracuje s téměř prázdným vstupem. Pro
zašuměná reálná data je λ větší a problém mírnější, ale hranice není nikde
hlídána.

---

## 1. Přehled nálezů dle závažnosti

| # | Nález | Oblast | Závažnost | Dopad |
|---|---|---|---|---|
| F1 | GMM přes celočíselnou replikaci bodů — arbitrární škála BIC + background floor | `peaks.py` | **Vysoká** (metodická) | Počet píků, `bic_threshold` nemá stabilní význam |
| F2 | Klasifikace termů podle šířky píku je závislá na λ | `term_classification.py` | **Vysoká** (metodická) | Voigt vs CPE rozhodnutí nestabilní |
| F3 | Interakce auto-λ → spiky DRT → degenerovaný post-processing | `core.py`+`gcv.py` | **Střední** | GMM/klasifikace na nízkošumových datech bez smyslu |
| F4 | Double-counting odporů u překrývajících se píků — ✅ OPRAVENO (v0.13.17) | `peaks.py`, `core.py` | **Střední** | Σ R_i ≠ R_pol, nadhodnocené dílčí odpory |
| F5 | GCV používá NNLS reziduum, ale lineární trace(K) | `gcv.py` | **Střední** (dokum.) | Efektivní DoF nezohledňuje aktivní omezení |
| F6 | L-curve roh přes argmax(křivost) — znaménko/orientace neověřeny | `gcv.py` | **Střední** | Riziko výběru nesprávného „rohu" |
| F7 | λ na hranici rozsahu se tiše akceptuje (bez varování) | `gcv.py` | Nízká–stř. | Pinning na mez není reportován |
| F8 | R_inf jako medián 5 nejvyšších frekvencí — nadhodnocení při HF obloucích | `core.py` | Nízká–stř. (záměr) | Posun celého γ přes `b = Z'−R_inf` |
| F9 | τ-mřížka omezena jen na měřený rozsah (bez auto-extend) | `core.py` | Nízká | Píky u krajů τ uříznuty/zkresleny |
| F10 | R_pol integrováno lichoběžníkem, jádro používá obdélník | `core.py` | Nízká | Zanedbatelná nekonzistence (≤ 1 interval) |
| F11 | DRT nereprezentuje induktivní data (γ≥0, A_im≤0) | `core.py` | Nízká (záměr) | Body Z''>0 zvyšují reziduum; jen varování |
| F12 | Větev `uncertain` v klasifikaci je mrtvý kód | `term_classification.py` | Kosmetika | Dosažitelná jen přesnou rovností |
| F13 | Testy λ-výběru jsou jen smoke; matematika neověřena; matice duplikovány | `tests/` | Střední (proces) | Drift test↔produkce neviditelný |

---

## 2. Co je matematicky správné (ověřeno)

Aby audit nebyl jen seznam výhrad — následující je korektní a empiricky
potvrzené:

### 2.1 Diskretizace a jádro `A`

Model `Z(ω) = R_inf + ∫ γ(τ)/(1+iωτ) d ln τ` je diskretizován obdélníkovým
pravidlem na log-rovnoměrné τ-mřížce:

```
A_re[k,m] = Δlnτ / (1 + (ω_k τ_m)²)
A_im[k,m] = -ω_k τ_m · Δlnτ / (1 + (ω_k τ_m)²)
```

(`core.py:329-330`). To přesně odpovídá reálné a imaginární části jádra.
Δlnτ je konstanta (logspace), takže obdélníkové pravidlo je konzistentní.
Pravá strana `b = [Z'−R_inf, Z'']` (`core.py:335`) sedí na model. ✓

### 2.2 Znaménková konvence

`A_im ≤ 0` a `γ ≥ 0` ⟹ modelové `Z'' ≤ 0`, což odpovídá kapacitnímu
(EIS) chování s `Z.imag < 0`. Vnitřně konzistentní. ✓

### 2.3 Rekonstrukce

`Z_recon = R_inf + (A_re + i·A_im) @ γ` (`core.py:769`) dává zpět
`Z' + i·Z''`. Ověřeno: rekonstrukční chyba 0.18 % pro ideální Voigt. ✓

### 2.4 Regularizace

`L` je operátor 2. derivace `[1, -2, 1]` tvaru `(n_tau-2, n_tau)`
(`core.py:339-342`) — standardní hladkostní Tikhonov. Regularizovaná
soustava `[A; √λ·L]` řešená přes `scipy.optimize.nnls` je korektní
non-negative Tikhonov. ✓

### 2.5 GCV trace — lineární část

`K = A·(AᵀA + λ·LᵀL)⁻¹·Aᵀ`, `trace(I−K) = n − trace(K)`
(`gcv.py:71-79`). Pro **lineární** Tikhonov je to korektní hat-matice a
efektivní počet stupňů volnosti. Ověřeno: `trace(K)` kladná, monotónně
klesá s λ, `trace(I−K)` zůstává kladná (112…127 pro λ 1e-5…1, n=140). ✓
(Omezení viz F5.)

---

## 3. Detailní nálezy

### F1 — GMM přes celočíselnou replikaci bodů je metodicky nekorektní (Vysoká)

`peaks.py:68-94`:

```python
weights_normalized = weights / np.mean(weights)
weights_int = np.maximum(1, np.round(weights_normalized).astype(int))
X_weighted = np.repeat(X, weights_int, axis=0)
...
gmm.fit(X_weighted)
bic = gmm.bic(X_weighted)
```

DRT γ(τ) je **vážená hustota**, ne vzorek. GMM se sem aplikuje tak, že se
každý τ-bod replikuje celočíselně úměrně γ. Tři problémy:

1. **BIC nemá stabilní škálu.** `BIC = -2·lnL + k·ln(N)`, kde `N` = počet
   replikovaných bodů. `N` je umělé číslo dané normalizací (`/mean`) a
   zaokrouhlením. Penalizace složitosti `k·ln(N)` proto nemá fyzikální ani
   statistický základ a **`bic_threshold=10` není srovnatelný napříč
   datasety**. Replikace navíc nepřináší nezávislou informaci, takže
   `lnL` je „přepočítaná" — likelihood je nafouknutá počtem kopií.

2. **`max(1, …)` floor injektuje pozadí.** Každý bin s γ>0, byť nepatrným,
   dostane ≥1 vzorek. To přidá uniformní pozadí, které **táhne odhady μ a
   σ ke středu τ-mřížky a nafukuje zdánlivé šířky** píků (přímý vstup do F2).

3. **Kvantizace ztrácí informaci** — γ < mean/2 se zaokrouhlí na 1 nebo 0.

**Dopad:** počet detekovaných píků a jejich šířky jsou částečně artefakt
metody, ne dat. **Doporučení:** použít vážený EM (GMM s `sample_weight`,
pokud verze sklearn podporuje — novější ano), nebo fitovat hustotu jiným
nástrojem (vážený KDE / NNLS dekonvoluce na gaussovský slovník). Pro
model-selection použít kritérium, jehož `N` je počet skutečných měření
(frekvencí), ne replik.

---

### F2 — Klasifikace termů podle šířky píku je závislá na regularizaci (Vysoká)

`term_classification.py:75, 108-137`:

```python
peak_width_decades = np.log10(tau_upper) - np.log10(tau_lower)  # = 4σ
...
if peak_width_decades < VOIGT_WIDTH_THRESHOLD:   # 1.5
    type = 'voigt'
elif peak_width_decades > VOIGT_WIDTH_THRESHOLD:
    type = 'cpe'
```

**Fyzikální problém:** ideální R‖C (Voigt) má v DRT **Diracovu δ-funkci**
(nulová šířka). Veškerá pozorovaná šířka píku pochází z **regularizace λ**
a diskretizace — nikoli z fyziky. Šířka tedy **není vlastnost elementu**, ale
vlastnost zvoleného λ:

- malé λ → úzké spiky → i pravý CPE vypadá jako Voigt;
- velké λ → široké píky → i pravý Voigt vypadá jako CPE.

Absolutní prahy (1.5 / 2.0 dekády) proto zaměňují „hodně regularizace" za
„CPE". Komentář v kódu sám uvádí očekávanou šířku ~1.14 dekády pro single RC
— jenže to platí jen pro jednu konkrétní úroveň vyhlazení.

**Doporučení:** klasifikaci buď (a) vázat na λ-normalizovanou šířku (porovnat
s šířkou, kterou by stejné λ dalo na syntetické δ), nebo (b) přestat tvrdit
typ z DRT a ponechat určení n na circuit fittingu (jak už kód dělá pro
`R_estimate`). Minimálně dokumentovat, že výstup je orientační.

---

### F3 — auto-λ → spiky DRT → degenerovaný post-processing (Střední)

**Empiricky doloženo.** Pro čistý jednopíkový Voigt vrátil hybridní výběr
`λ ≈ 1.5e-5` (spodní mez rozšířeného rozsahu), což dalo DRT s **pouhými
4 nenulovými body γ**. GMM s `n_components` až 6 pak fituje téměř prázdný
vstup; klasifikace šířky ztrácí smysl.

Příčina: výběr λ (GCV i L-curve) optimalizuje **shodu s daty**. Pro
nízkošumová data je optimální λ → 0 (data jdou fitovat „na spike"), protože
δ-řešení je pro ideální Voigt korektní. Post-processing píků ale implicitně
předpokládá **hladké** γ.

**Doporučení:** (a) hlídat „spikiness" výsledku (např. počet nenulových
binů / efektivní šířka) a varovat, když je DRT příliš řídký pro analýzu
tvaru; (b) pro účely detekce píků zvážit mírně vyšší λ než pro rekonstrukci
(rozpojit „λ pro fit" a „λ pro tvarovou analýzu"); (c) propojit s F7
(detekce λ na hranici).

---

### F4 — Double-counting odporů u překrývajících se píků (Střední) — ✅ OPRAVENO (v0.13.17)

Dvě nezávislé cesty integrují **celkové** γ přes okno každého píku:

- `peaks.py:189-199` (GMM): `R_estimate = trapz(γ[idx_lower:idx_upper], lnτ)`
  přes ±2σ okno daného píku.
- `core.py:178-209` (`_estimate_peak_resistance`, scipy): integrace γ od
  prahu 10 % výšky vlevo/vpravo.

Pokud se píky překrývají, **oblast překryvu se započítá vícekrát** a
`Σ R_i > R_pol`. U GMM je správné integrovat **příspěvek dané komponenty**
(`weight_i · R_pol`, protože Σ weight = 1), nikoli celkové γ v okně.

**Doporučení:** pro GMM použít `R_i = weight_i · R_pol_from_gamma`
(rozklad jednotky), což zaručí `Σ R_i = R_pol`. Pro scipy cestu rozdělit
osu v minimech mezi píky.

**Oprava (v0.13.17):** obě dokumentované cesty opraveny.
- GMM (`peaks.py`): `R_estimate = weight_i · R_pol` (Σ weight = 1) →
  `Σ R_i = R_pol` přesně.
- scipy (`core.py:_estimate_peak_resistance`): τ-osa rozdělena v údolích
  (minimech γ) mezi sousedními píky; díky aditivitě `trapz` přes sdílený
  uzel je `Σ R_i` rovno R_pol pokrytého rozsahu bez překryvu. Odstraněn
  nepoužívaný parametr `tolerance`. Regresní testy
  `tests/test_peak_resistance.py` (4) zamykají `Σ R_i = R_pol`.

**POZN. mimo rozsah F4:** stejný okénkový vzor (integrace γ od prahu 10 %
výšky) je i v `fitting/auto_suggest.py:269-289`, kde nastavuje **počáteční
R pro circuit fitting**. Tato cesta nebyla v F4 pojmenována a mění chování
fittingu — neopraveno, k samostatnému rozhodnutí.

---

### F5 — GCV pro NNLS: nekonzistentní reziduum vs. trace (Střední, dokumentováno)

`gcv.py:52-86`: reziduum `||b − A·x||` se počítá z **NNLS** řešení
(s omezením x≥0), ale `trace(K)` z **lineární** hat-matice (bez omezení).
Efektivní DoF lineárního Tikhonova nadhodnocuje skutečné DoF NNLS, protože
počítá i sloupce, které NNLS vynuloval (aktivní omezení). Kód to **poctivě
komentuje** (`gcv.py:65-67`) a hybrid to mitiguje L-curve metodou.

Matematika trace je přitom korektní pro lineární případ (viz 2.5). Jde tedy
o **vědomou aproximaci**, ne chybu. **Doporučení:** pokud má GCV hrát větší
roli, odvodit efektivní DoF z aktivní množiny NNLS (počet nenulových /
volných parametrů, příp. projekce na aktivní sloupce). Jinak ponechat a
spoléhat na L-curve (současný stav je rozumný).

---

### F6 — L-curve roh přes argmax(křivost): znaménko a orientace neověřeny (Střední)

`gcv.py:170-176, 197-208`:

```python
numerator = d_rho * dd_eta - dd_rho * d_eta
denominator = (d_rho**2 + d_eta**2)**1.5
curvature = numerator / denominator
...
corner_idx_inner = np.argmax(curvature_inner)   # maximum KLADNÉ křivosti
```

Roh L-křivky = bod maximální křivosti. Použité `argmax` ale hledá **maximum
znaménkové křivosti**, jejíž **znaménko závisí na pořadí (ρ, η) a směru
průchodu** (rostoucí λ). Při opačné orientaci je roh v **minimu** (nejvíce
záporná křivost), a `argmax` by trefil opačný (plochý) konec. Standardní
(Hansenova) definice řeší orientaci explicitně.

Není to nutně chyba — pro zdejší konvenci (ρ=reziduum první, rostoucí λ)
to *může* být správně — ale **není to nikde ověřeno testem** na křivce se
známým rohem. Vzhledem k tomu, že hybrid často preferuje právě L-curve λ,
je chybný roh přímou chybou výběru λ.

**Doporučení:** přidat unit test s analytickou/synteticky L-křivkou se
známým rohem a ověřit, že `find_lcurve_corner` jej trefí; případně
převzít robustní variantu (např. maximalizace |κ| jen na konvexní větvi,
nebo Hansenova trojúhelníková metoda).

---

### F7 — λ na hranici rozsahu se tiše akceptuje (Nízká–střední)

`gcv.py:254-263` rozšiřuje jemné hledání, pokud minimum padne na okraj
(`×/÷10`), ale výsledné λ může pořád ležet na nové hranici a **vrátí se bez
varování**. (Kontrast: `peaks.py:158-163` na hraniční optimum upozorňuje.)
Empiricky λ skončilo na 1.5e-5 (mez), což je signál „dej míň regularizace,
než dovolím" — uživatel to nevidí.

**Doporučení:** detekovat, že vybrané λ ≈ mez rozsahu, a přidat varování do
diagnostiky (analogicky GMM). Souvisí s F3.

---

### F8 — R_inf jako medián 5 nejvyšších frekvencí (Nízká–střední, záměr)

`core.py:225-227`: `R_inf = median(Z.real[top-5 f])`. Při existenci
vysokofrekvenčních oblouků (rychlé procesy) to **nadhodnocuje** pravé R_inf,
a protože `b = Z' − R_inf`, posune se **celé γ**. K dispozici jsou robustnější
volby (`use_rl_fit`, `use_voigt_fit`). Jde o vědomou defaultní heuristiku;
stačí dokumentovat citlivost.

---

### F9 — τ-mřížka omezená jen na měřený rozsah (Nízká)

`core.py:317-319`: `τ ∈ [1/(2πf_max), 1/(2πf_min)]`. Píky s τ u krajů
(zejména pomalé procesy u nejnižších f) jsou uříznuty / zkresleny okrajovým
efektem 2. derivace. Pro KK validaci byl auto-extend zaveden (v0.13.10), pro
DRT nikoli. Pro DRT je užší mřížka konzervativní volba, ale stojí za zmínku
do dokumentace.

---

### F10 — R_pol: lichoběžník vs. obdélník jádra (Nízká)

`core.py:750`: `R_pol_from_gamma = trapz(γ, lnτ)`, zatímco jádro a tím i
modelové `Z'(0)−R_inf = Σγ_m·Δlnτ` používají obdélník. Rozdíl je ≤ jeden
krajní interval (ověřeno: 100.04 vs 100.0). Pro konzistenci by `R_pol`
mohlo být `Σγ_m·Δlnτ` (vlastní R_pol modelu). Zanedbatelné.

---

### F11 — DRT nereprezentuje induktivní data (Nízká, záměr)

`A_im ≤ 0` a `γ ≥ 0` ⟹ body s `Z'' > 0` (induktivní) jsou
nereprezentovatelné; NNLS je „zploští" a zvýší reziduum. Kód na induktivní
podíl varuje (`core.py:420-424`), ale stejně fituje. Inherentní limitace
γ≥0 DRT. OK, jen dokumentovat.

---

### F12 — Mrtvá větev `uncertain` (Kosmetika)

`term_classification.py:108-147`: podmínky `< THRESHOLD` a `> THRESHOLD`
pokrývají vše kromě přesné rovnosti, takže `else: 'uncertain'` je prakticky
nedosažitelná. Buď zavést pásmo nejistoty kolem prahu (`THRESHOLD ± δ`),
nebo větev odstranit.

---

### F13 — Testy λ-výběru jsou jen smoke; matice duplikovány (Střední, proces)

`tests/test_hybrid_lambda.py` ověřuje pouze, že λ je kladné, konečné a
„v rozumném rozsahu". **Neexistuje** test, který by:
- ověřil rekonstrukci známého γ / poloh píků z `calculate_drt`;
- testoval `compute_gcv_score`, `compute_lcurve_curvature`,
  `find_lcurve_corner` (F6!);
- pokryl `peaks.py` (GMM) a `term_classification.py` (0 testů).

Navíc test **re-implementuje `build_drt_matrices`** místo importu
`core._build_drt_matrices`, takže odchylka mezi testovou a produkční
konstrukcí matic by zůstala neviditelná.

**Doporučení:** přidat korektnostní testy (recovery známých spekter),
test rohu L-křivky a importovat produkční konstrukci matic.

---

## 4. Priority k řešení (přínos/úsilí)

1. **F6 + F13** — přidat test rohu L-křivky a korektnostní testy DRT
   (recovery γ). Nízké úsilí, vysoká hodnota: zamkne to současné chování
   a odhalí případný špatný roh. Importovat `_build_drt_matrices` do testů.
2. ✅ **F4** — opravit double-counting: u GMM `R_i = weight_i · R_pol`,
   u scipy rozdělení v údolích. **Hotovo (v0.13.17)**; zvážit ještě
   stejnou opravu v `auto_suggest.py` (mimo rozsah F4).
3. **F1 + F3 + F7** — ozdravit GMM/λ pipeline: vážený EM místo replikace,
   detekce a varování na λ u meze a na „příliš řídký" DRT.
4. **F2** — přehodnotit klasifikaci šířkou (λ-normalizace nebo přesun
   rozhodnutí na circuit fitting); minimálně dokumentovat omezení.
5. **F5, F8–F12** — dokumentovat jako vědomá omezení / drobnosti; opravit
   F12 a F10 při příležitosti.

---

## 5. Závěr

Numerické jádro DRT (sestavení a řešení regularizované soustavy) je
**správné a spolehlivé**. Rizika leží v **interpretační nadstavbě**: GMM
detekce píků a klasifikace termů produkují čísla, která vypadají objektivně,
ale jsou závislá na λ a na arbitrárních škálách (replikace, absolutní prahy
šířky). Nejde o početní chyby, ale o **metodické předpoklady, které nejsou
splněny** pro řídká/nízkošumová spektra a nejsou hlídány. Nejvýnosnější
kroky jsou levné: testy korektnosti (F13/F6) a oprava double-countingu (F4);
zbytek je o ozdravení heuristik a poctivé dokumentaci omezení.
