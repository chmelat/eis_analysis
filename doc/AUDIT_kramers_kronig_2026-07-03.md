# Audit modulu Kramers-Kronig validace

**Datum:** 2026-07-03 | **Verze:** 0.16.18 | **Typ:** kritické code-level review (Claude Code)

## Stav oprav (aktualizováno 2026-09-02)

| Nález | Stav | Verze | Commit |
|-------|------|-------|--------|
| K1    | Opraveno | 0.16.19 | 2354c2d |
| K2    | Opraveno | 0.16.20 | 67ee6b9 |
| K10   | Opraveno | 0.21.6 | 9f9991b |
| K3-K9, K11 | Otevřené | — | — |

Stav K3-K9 a K11 ověřen k 2026-09-02 (v0.30.0). Čísla řádků v nálezech
níže jsou z data auditu a od té doby se posunula — K8 je dnes na ř. 532.

**Rozsah:** `validation/kramers_kronig.py` a podpůrné jádro ve
`fitting/voigt_chain/` (mu_optimization, tau_grid, fitting.estimate_R_linear,
validation), konzumace v `cli/handlers/validation.py`, testy
(`tests/test_kramers_kronig.py`, 24 testů). Modul zhit.py mimo rozsah.

## Verdikt

Matematické jádro je věrné literatuře: chi²_ps dle Boukamp (1995), tau-grid,
mu metrika a iterace M dle Schönleber (2014), extrakce L z imaginárního
rezidua korektní vážené LS (numericky ověřeno). Výchozí KK cesta
(fit_type='real', modulus) počítá správně. Rizika: numericky rozbitá
rekonstrukce R_s ve fit_type='imag' (dosažitelná z CLI), matoucí sémantika
mu a diagnostika, která zůstává jen v logu.

## Nálezy

### K1 — R_s rekonstrukce ve fit_type='imag' je numericky rozbitá (Vysoká)

`voigt_chain/fitting.py` ř. 227-231. Numericky ověřeno na KK-konzistentních
syntetických datech (R_s = 100 + Voigt(5000, 5e-3)):

| weighting | vrácené R_s | správně |
|-----------|-------------|---------|
| modulus   | −4513       | ~100    |
| sqrt      | −15731      | ~100    |
| uniform   | −36596      | ~100    |

Dvě složené chyby: (a) `z_re_fit` se počítá z *vážené* matice A_real
a „odvažuje" se násobením |Z| — platné jen pro w = 1/|Z| **bez** normalizace
vah na průměr 1 (ř. 130), která vnáší faktor mean(1/|Z|); (b) pro
uniform/sqrt/proportional je násobení |Z| špatně úplně. Reálná část Z_fit
je pak nesmysl. Výchozí KK cesta ('real') postižená není; chyba je
dosažitelná přes `--voigt-chain --voigt-fit-type imag` a přímé API.
**Náprava:** z_re_fit stavět z nevážené design matice (jako v 'real' větvi
pro L) a |Z| násobení odstranit; regresní test se syntetickými daty.

**OPRAVENO (0.16.19, commit 2354c2d):** predikce se staví z nevážené
matice (váhy vyděleny zpět); po opravě R_s = 99.85-100.00 pro všechna
čtyři vážení (M=15). Poznatek z implementace: zbytková odchylka ~30 %
u uniform při M=10 je diskretizační chyba tau-gridu, ne formule —
od M=15 pod 0.3 %. Regresní testy v tests/test_voigt_chain.py
(všechna vážení, mΩ případ, end-to-end imag fit).

### K2 — Matoucí sémantika mu (Střední)

Vrácené mu je z konstrukce algoritmu **vždy < mu_threshold** při úspěšném
běhu — je to hodnota při zastavení (první M, kde mu klesla pod práh), ne
kvalita dat. Docstringy (KKResult, LinKKResult) ale říkají „mu (1.0 = no
overfit, <0.85 = overfit)" a CLI hodnotu tiskne do titulku grafu — uživatel
čte „mu=0.76" jako indikaci špatných dat. **Náprava:** opravit docstringy
a popis v grafu/logu (mu je stop hodnota Lin-KK iterace).

**OPRAVENO (0.16.20, commit 67ee6b9):** sjednoceno na sémantiku stop
hodnoty v docstringech, README, PYTHON_API i CLI help textech; CLI log
nově „mu=0.73 (Lin-KK stop, threshold 0.85)", graf „stop mu=". Nad rámec
auditu nalezeno a opraveno obrácené tvrzení ve find_optimal_M_mu
(„lower values -> fewer elements" — správně: nižší práh zastaví iteraci
později, tedy více elementů) a chybný popis v PYTHON_API („>0.85 =
good"). Regresní test v tests/test_kramers_kronig.py.

### K3 — Diagnostika se nepropaguje do dat (Střední)

- „Reached max_M, model may still be overfit" existuje jen v logu
  `find_optimal_M_mu`; do LinKKResult/KKResult.warnings se nedostane —
  v rozporu s deklarací „all diagnostics returned as data" v headeru.
- Obráceně: KKResult.warnings („Data may contain artifacts") CLI nikdy
  nečte a kvalitu si počítá vlastní logikou (duplicitní kritérium).

### K4 — Práh 5 % triplikovaný; is_valid (průměr) nesedí s grafem (per bod) (Střední)

`is_valid` testuje *průměr* |reziduí| < 5 %; graf kreslí pásmo ±5 % *per
bod*. Spektrum s 10% rezidui na okraji a 1 % jinde projde jako validní,
i když body leží mimo pásmo. Práh není pojmenovaná konstanta (dataclass
properties 2x + plot); properties mean_residual_*/is_valid duplikované
v KKResult i LinKKResult.

### K5 — Wrapper vždy staví matplotlib figuru (Střední)

`kramers_kronig_validation` nemá `plot=False` — porušuje princip oddělení
vizualizace od algoritmů (lin_kk_native je čisté). V dávkovém API
zpracování se hromadí otevřené figury.

### K6 — Zmatený komentář v tau_grid (Nízká)

`tau_grid.py` ř. 149-150: „Lin-KK standard: tau = 1/omega, not
1/(2*pi*f)" — jenže 1/omega ≡ 1/(2*pi*f), identické veličiny. Kód je
správně (přesně Schönleber), komentář tvrdí neexistující rozpor.

### K7 — Nedokumentované konstanty (Nízká)

- 5000 v `estimate_noise_percent` (= 100²/2, tj. sigma = 100*sqrt(chi²/2N);
  odvození chybí v kódu),
- tolerance 0.001 a n_evaluations=11 v `find_optimal_extend_decades`.

### K8 — Mrtvá větev (Nízká)

`kramers_kronig.py` ř. 486: `if lkk.Z_fit is None` — lin_kk_native nikdy
None nevrátí (buď raise, nebo pole).

### K9 — Tiché polykání výjimek (Nízká)

Catch-all `except Exception` → `KKResult(error=...)` jen
s `logger.debug` — i programátorská chyba se promění ve „validace
selhala" a debug log není běžně vidět. Minimálně logger.warning.

### K10 — Duplikace výpočtu vah (Nízká)

`estimate_R_linear` má vlastní kopii `diagnostics.compute_weights`
(identické větve i normalizace na průměr 1).

**OPRAVENO (0.21.6, commit 9f9991b):** kopie odstraněna,
`voigt_chain/fitting.py` importuje `compute_weights` z `..diagnostics`.

### K11 — CLI import privátního symbolu (Nízká)

`cli/handlers/validation.py` importuje `_quality_label` ze zhit.

## V pořádku (ověřeno)

chi²_ps přesně dle Boukampa (váha 1/|Z|²); noise formule konzistentní
se sigma = 100*sqrt(chi²/2N); tau-grid tau_min = 1/omega_max,
tau_max = 1/omega_min, log rozestup, přesně M bodů (Schönleber); extenze
jen k nižším frekvencím, záporné hodnoty clipnuté; mu metrika i iterace
M=3..max_M dle Lin-KK; extrakce L z imaginárního rezidua při real fitu
(vážené LS) numericky ověřena; reconstruct_impedance indexuje správně
[R_s, R_1..R_M, L]; vstupní validace hlídá empty/NaN/nekladné frekvence;
auto-extend s tie-breakem k nulové extenzi; 24 testů včetně auto-extend
regresí.

## Priority

1. ~~K1~~ — HOTOVO (0.16.19).
2. K3 — propagace varování do dat (~~K2~~ HOTOVO, 0.16.20).
3. K4 + K5 — konzistence kritéria a plot=False.
4. K6-K9, K11 — drobnosti průběžně (~~K10~~ HOTOVO, 0.21.6).
