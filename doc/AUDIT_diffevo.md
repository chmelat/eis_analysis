# Audit modulu `diffevo` — matematická korektnost

**Datum:** 2026-06-30
**Rozsah:** `eis_analysis/fitting/diffevo.py` a navazující moduly
(`jacobian.py`, `covariance.py`, `diagnostics.py`, `bounds.py`,
`circuit_builder.py`)
**Metoda:** statická analýza + numerické ověření analytického jakobiánu
proti centrálním diferencím + spuštění `tests/test_diffevo.py` (18/18 prošlo)

---

## Co je ověřeno jako správné

### Analytický jakobián (`jacobian.py`)
Numericky ověřen proti centrálním konečným diferencím, shoda 1e-6 až 1e-9:

| Případ                | max. rel. chyba |
|-----------------------|-----------------|
| R-C, R-L              | ~1e-9           |
| Randles, R\|Q, R-K, R-G | ~1e-6         |
| R-Wo, R-W             | ~1e-8           |
| vnořený R\|(R-C)      | ~3e-9 ✱         |

✱ Hrubě reportovaná chyba 1.6e+02 byl falešný poplach — týká se prvku, kde
analytická hodnota je ~1.6e-10 (fakticky nula) a FD dá přesnou 0; absolutní
i sloupcové chyby jsou ~1e-9.

Konkrétně potvrzeno:
- Derivace všech prvků (R, C, L, Q, W, Wo, K, G) — včetně
  `dZ/dn = -Z·ln(jω)` s `ln(jω) = ln(ω) + jπ/2`, řetězového pravidla u Wo a
  Gerischera.
- Řetězové pravidlo pro Series (`dZ/dp = dZi/dp`) i Parallel
  (`dZ/dp = (Z/Zi)²·dZi/dp`), rekurzivně pro vnořené obvody.
- **Pořadí parametrů** je konzistentní mezi `impedance()` a `circuit_jacobian()`
  (obě iterují `self.elements` a krájí podle `len(elem.get_all_params())`) —
  sloupce jakobiánu odpovídají parametrům.
- **Znaménko reziduí**: residual = (Z_data − Z_fit)·w ⟹ ∂/∂p = −∂Z_fit/∂p·w.
- **Konzistence objektivů DE↔LS**: DE minimalizuje `Σ w²[(ΔRe)²+(ΔIm)²]`,
  `least_squares` minimalizuje `½Σ(Δ·w)²` — stejná účelová funkce (až na ½).

---

## Nálezy

### 1. [STŘEDNÍ] Výběrová metrika ≠ optimalizovaná účelová funkce
**Místo:** `diffevo.py:337`, `diffevo.py:340-352`

Rozhodnutí mezi DE a LS výsledkem (i hodnota `improvement`) používá
`compute_fit_metrics`, což je **vážená střední relativní chyba** (L1-typu), ne
vážené SSR (L2), které se reálně optimalizuje. Protože `trf` je sestupná metoda
z DE bodu, LS téměř vždy sníží SSR, ale tato jiná metrika se může nepatrně
zhoršit → kód zahodí lepší LS výsledek a zaloguje „Refinement worsened fit".

**Oprava:** výběr provádět na téže veličině, která se minimalizuje (vážené SSR).

### 2. [STŘEDNÍ] Kovariance počítaná v nesprávném bodě (větev „refinement worsened")
**Místo:** `diffevo.py:346-351`, `diffevo.py:357-366`

Když se zvolí DE výsledek, je `params_opt_free = de_result.x` a
`final_residuals` se vyhodnotí v DE bodě, **ale** `jacobian=ls_result.jac` je
jakobián v *LS* řešení. Tím se `s² = RSS/dof` (z DE reziduí) kombinuje
s `(JᵀJ)⁻¹` z jiného parametrického vektoru → nevalidní kovariance a stderr.

**Oprava:** přepočítat analytický jakobián v `params_opt_free`. (V běžné větvi,
kdy vyhraje LS, je vše konzistentní.)

### 3. [NÍZKÁ–STŘEDNÍ] Sémantika `condition_number`
**Místo:** `covariance.py:131`, docstring `CovarianceResult`

Počítá se `cond(J) = S[0]/S[-1]`, ale docstring/pole tvrdí „condition number of
JᵀJ". Kovariance závisí na JᵀJ, jehož číslo podmíněnosti je `cond(J)²`. Práh
„well-conditioned" `< 1e10` se tak aplikuje na `cond(J)`, čímž připustí
`cond(JᵀJ)` až ~1e20.

**Oprava:** buď umocnit, nebo opravit dokumentaci — interpretace je zavádějící.

### 4. [NÍZKÁ] Regularizace malých singulárních hodnot
**Místo:** `covariance.py:153-158`

Pro `s ≤ threshold` se nastaví `S_inv_sq = 1/threshold²` místo Moore-Penroseovy
0. Nafukuje to rozptyl neidentifikovatelných směrů o arbitrární hodnotu závislou
na `S[0]`. Jako „signalizace neidentifikovatelnosti přes velký stderr" je to
obhajitelné (bezpečnější než pseudoinverze, jež by falešně tvrdila nulovou
nejistotu), ale magnituda je arbitrární a není dokumentována jako záměrný odklon
od pseudoinverze.

### 5. [NÍZKÁ] Nekonzistentní stupně volnosti mezi moduly
**Místo:** `covariance.py:118` vs `covariance.py:249`

`compute_covariance_matrix` používá `dof = 2·n_freq − p` (správně pro rozdělená
komplexní rezidua), ale `compute_confidence_interval` používá `dof = n_freq − p`.
Z `diffevo` se CI nevolá, ale při kombinaci jinde dostane t-multiplikátor
nekonzistentní dof.

### 6. [NÍZKÁ] Reportovací metrika dvojnásobně počítá 1/|Z|
**Místo:** `diagnostics.py:155-159`

`compute_fit_metrics` s `modulus` váhováním tvoří
`Σ (1/|Z|)·|ΔZ|/|Z| / Σ(1/|Z|)`, tj. efektivně 1/|Z|² důraz na |ΔZ|. Pro
„relativní chybu" je to možná silnější než zamýšleno — vhodné potvrdit záměr.

---

## Shrnutí

Jádro matematiky (analytické derivace, řetězová pravidla, zarovnání parametrů,
konzistence DE/LS objektivu) je **korektní a numericky ověřené**. Skutečné vady
jsou dvě logické nekonzistence: výběr výsledku podle jiné metriky než se
optimalizuje (#1) a kovariance v nesprávném bodě v okrajové větvi (#2). Zbytek
jsou dokumentační/sémantické nepřesnosti.

| # | Závažnost     | Stav               |
|---|---------------|--------------------|
| 1 | Střední       | Opraveno (v0.16.2) |
| 2 | Střední       | Opraveno (v0.16.2) |
| 3 | Nízká–střední | Opraveno (v0.16.3) |
| 4 | Nízká         | Opraveno (v0.16.4) |
| 5 | Nízká         | Opraveno (v0.16.5) |
| 6 | Nízká         | Otevřeno           |
