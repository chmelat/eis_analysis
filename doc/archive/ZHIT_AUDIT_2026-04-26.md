# Z-HIT Validation Audit (2026-04-26)

Audit modulu `eis_analysis/validation/zhit.py` na základě externího code review.
Hodnotí 5 bodů z review + 1 bod nalezený navíc.

**Status: 5/6 opraveno ve verzi 0.13.6 (2026-04-26). Zbývá bod 2.**

## Souhrn

Kolega má pravdu prakticky ve všech 5 bodech. Body 1, 3 a 5 jsou chyby/díry,
které je potřeba opravit; body 2 a 4 jsou návrhové slabiny.

| # | Bod                                              | Závažnost | Stav         |
|---|--------------------------------------------------|-----------|--------------|
| 1 | Setřídění narušuje API contract                  | Vysoká    | ✅ Opraveno  |
| 2 | `np.gradient` zesiluje šum ve fázi               | Střední   | ⏳ Otevřeno   |
| 3 | Chybí `np.unwrap` na fázi                        | Nízká*    | ✅ Opraveno  |
| 4 | Magic number `center=1.5` v offset windowu       | Střední   | ✅ Opraveno  |
| 5 | Práh kvality 5 % je příliš shovívavý             | Střední   | ✅ Opraveno  |
| 6 | `is_valid` ignoruje `quality_threshold` (extra)  | Nízká     | ✅ Opraveno  |

\* Nízká pravděpodobnost, ale vysoký dopad pokud nastane.

---

## 1. Setřídění narušuje API contract — ✅ OPRAVENO (v0.13.6)

`zhit_validation` (řádky 327-329) data setřídí vzestupně, ale výstupní pole
v `ZHITResult` jsou ve setříděném pořadí:

```python
sort_idx = np.argsort(frequencies)
frequencies = frequencies[sort_idx]
Z = Z[sort_idx]
```

Uživatel, který předá data ve standardním sestupném pořadí (Gamry, BioLogic),
dostane zpět:

- `Z_mag_reconstructed`, `residuals_mag`, `residuals_real`, `residuals_imag`,
  `Z_fit` — všechny ve vzestupném pořadí

Pokud chce vykreslit rezidua proti svým původním frekvencím, vyrobí si
nesmysly. **Reálný footgun.**

**Oprava** je triviální — buď invertovat na konci přes
`inv_idx = np.argsort(sort_idx)`, nebo (čistěji) vyžadovat setříděný vstup
v docstringu a třídění nedělat tady, ale v CLI vrstvě.

## 2. `np.gradient` zesiluje šum ve fázi — ⏳ OTEVŘENO

Řádek 227:

```python
d_phi_d_ln_omega = np.gradient(phi, ln_omega)
```

Korekce druhého řádu se násobí `gamma = -π/6 ≈ -0.524` (řádek 233), což není
zanedbatelné. Pokud má fáze šum (typicky vysokofrekvenční konec kvůli
kabeláži/CPE), derivace ho spike-uje a propíše se přímo do magnitudy.

**Návrh** s volitelným Savitzky-Golay smoothingem před derivací je rozumný
kompromis. Případně argument typu `phase_smooth_window: int = 0` (0 = bez
smoothingu, zachová stávající chování), aby to bylo backward-compatible.

## 3. Chybí `np.unwrap` na fázi — ✅ OPRAVENO (v0.13.6)

Řádek 333:

```python
phi = np.arctan2(Z.imag, Z.real)
```

U typických pasivních systémů je fáze v intervalu [-π/2, 0], takže k přeskoku
nedojde. Ale:

- Induktivní smyčky (vysokofrekvenční kabely, baterie) mohou jít mimo
  [-π/2, π/2]
- Šum těsně u ±π hranice způsobí umělé skoky

`np.unwrap` je 2 znaky kódu navíc s nulovým rizikem regrese pro normální data:

```python
phi = np.unwrap(np.arctan2(Z.imag, Z.real))
```

Udělat bez výjimky.

## 4. Magic number `center=1.5` — ✅ OPRAVENO (v0.13.6)

Řádky 113, 191, 284: default `center=1.5` (log10(31.6 Hz)).

Argument JE vystaven jako parametr, takže uživatel ho může přepsat. Ale default
je pevný a může být mimo měřený rozsah (např. spektrum 1 mHz - 100 mHz pro
koroze povlaků). Pak Gaussova váha vyhasne a offset bude nepředvídatelný.

**Lepší default**: medián `log10(frequencies)` uvnitř
`_calculate_offset_weighted`, pokud není explicitně předán. Behavior pro
existující testy s plným rozsahem bude prakticky stejný (medián log-f vyjde
typicky okolo 1-2).

## 5. Práh kvality 5 % je příliš shovívavý — ✅ OPRAVENO (v0.13.6)

Řádky 388, 398-401, 105-106:

```python
quality = max(0.0, 1.0 - mean_abs_residual_mag / quality_threshold)  # default 5.0
...
if mean_abs_residual_mag < 5.0:
    logger.info("Data quality is good (residuals < 5%)")
...
def is_valid(self) -> bool:
    return self.mean_residual_mag < 5.0
```

Pro KK-konzistentní data jsou rezidua typicky **0.1–0.5 %**. 5 % už je signál,
že něco není v pořádku (drift, nelineární odezva, šum). Plošná zpráva "good"
je matoucí.

**Návrh** — odstupňovat:

| Rezidua    | Hodnocení                                    |
|------------|----------------------------------------------|
| < 0.5 %    | excellent                                    |
| < 1.0 %    | good                                         |
| < 2.5 %    | acceptable                                   |
| < 5.0 %    | marginal — check for drift/nonlinearity      |
| ≥ 5.0 %    | poor                                         |

A `is_valid` vázat na realističtější hranici (např. 1-2 %), nebo to udělat
parametrické.

## 6. Interní nekonzistence v `is_valid` (extra) — ✅ OPRAVENO (v0.13.6)

Řádky 105-106:

```python
@property
def is_valid(self) -> bool:
    """Check if data passes validation (magnitude residuals < 5%)."""
    return self.mean_residual_mag < 5.0
```

Hardcoded `< 5.0` ignoruje uživatelem nastavený `quality_threshold` v
`zhit_validation`. Buď uložit threshold do `ZHITResult` jako field a v
`is_valid` ho použít, nebo dokumentovat, že `is_valid` má fixní hranici a
přidat property `is_valid_at(threshold)`.

---

## Doporučené pořadí oprav

1. ✅ **Bod 3** (`np.unwrap`) — 1 řádek, nulové riziko *(v0.13.6)*
2. ✅ **Bod 1** (un-sort výstupů) — pár řádků, opravuje API contract *(v0.13.6)*
3. ✅ **Bod 5 + 6** (kvalitativní stupně, propojit `is_valid` s `quality_threshold`) *(v0.13.6)*
4. ✅ **Bod 4** (default `center` = medián log-f) *(v0.13.6)*
5. ⏳ **Bod 2** (volitelný smoothing fáze) — přidává API surface, nech jako poslední

Po každém kroku spustit `python3 -m pytest tests/ -k zhit` a ověřit, že
existující testy pořád procházejí (nebo upravit testy, kde se mění očekávaný
output).
