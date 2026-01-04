# Parser ekvivalentnich obvodu

Tento dokument popisuje implementaci parseru pro definici ekvivalentnich obvodu
v EIS Analysis Toolkit.

## Architektura

Parser pouziva **operator overloading** pristup inspirovany knihovnou
EISAnalysis.jl (Julia). Misto textoveho parseru se obvod definuje primo
v Pythonu pomoci operatoru.

### Klicove moduly

```
eis_analysis/fitting/
├── circuit_elements.py   # Zakladni elementy (R, C, Q, L, W, Wo, K)
├── circuit_builder.py    # Kombinatory (Series, Parallel)
└── circuit.py            # Fitting funkce
```

### Tridy

**CircuitElement** (abstraktni bazova trida):
- Reprezentuje jeden obvodovy prvek
- Definuje operatory `-` (serie) a `|` (paralel)
- Uchovava parametry a jejich fixed/free status

**Series** a **Parallel** (kompozitni obvody):
- Reprezentuji seriove/paralelni spojeni
- Rekurzivne obsahuji dalsi elementy nebo kompozity
- Implementuji `impedance()`, `get_all_params()`, `update_params()`


## Parsing pomoci eval()

Retezec obvodu se parsuje v `eis.py` funkci `parse_circuit_expression()`:

```python
def parse_circuit_expression(expr: str):
    safe_namespace = {
        'R': R, 'C': C, 'Q': Q,
        'L': L, 'W': W, 'Wo': Wo, 'K': K
    }
    circuit = eval(expr, {"__builtins__": {}}, safe_namespace)
    return circuit
```

Pouziti `eval()` s omezenym namespace zajistuje:
- Pristup pouze k definovanym elementum
- Zadne built-in funkce (bezpecnost)
- Plna syntaxe Pythonu pro vyrazy


## Priorita operatoru

**Kriticka informace:** Protoze se pouziva Pythonovy `eval()`, plati
**Pythonova priorita operatoru**, ne intuitivni priorita pro elektricke obvody.

### Pythonova priorita (od nejvyssi po nejnizsi)

| Priorita | Operator | Vyznam v obvodech |
|----------|----------|-------------------|
| 13       | `*`, `/` | Skalovani parametru |
| 12       | `-`, `+` | **Seriove spojeni** |
| 8        | `\|`     | **Paralelni spojeni** |

**Operator `-` ma VYSSI prioritu nez `|`!**

### Dusledky

Bez zavorek se vyraz parsuje takto:

```python
# Uzivatel napise:
R(1) - R(2)|C(3) - R(4)|C(5)

# Python interpretuje (kvuli priorite - > |):
R(1) - (R(2) | (C(3) - (R(4) | C(5))))

# Vysledna struktura:
#   R(1)
#     │
#   Serie
#     │
#   R(2) ══╗
#          ║ Parallel
#   C(3)   ║
#     │    ║
#   Serie  ║
#     │    ║
#   R(4) ══╬══╗
#          ║  ║ Parallel
#   C(5) ══╝══╝
```

Toto NENI zamysleny obvod `R - (R||C) - (R||C)`!

### Spravny zapis

Vzdy pouzivejte explicitni zavorky kolem paralelnich kombinaci:

```python
# Spravne:
R(1) - (R(2)|C(3)) - (R(4)|C(5))

# Vysledna struktura:
# R(1) - [R(2)||C(3)] - [R(4)||C(5)]
```

### Doporuceni

1. **Vzdy** pouzivejte zavorky kolem `(R|C)` kombinaci
2. Nezapomente na zavorky i u slozitejsich paralelnich kombinaci
3. Pri nejistote zkontrolujte strukturu pomoci `print(circuit)`


## Podporovane elementy

### R - Rezistor

```python
Z_R = R
```

| Parametr | Jednotka | Default | Popis |
|----------|----------|---------|-------|
| R        | Ohm      | 100     | Odpor |

```python
R(100)      # 100 Ohm, volny parametr
R("100")    # 100 Ohm, fixovany parametr
R()         # default 100 Ohm
```

### C - Kapacitor

```python
Z_C = 1 / (j * omega * C)
```

| Parametr | Jednotka | Default | Popis |
|----------|----------|---------|-------|
| C        | F        | 1e-6    | Kapacita |

```python
C(1e-6)     # 1 uF
C("1e-6")   # fixovany
```

### L - Induktor

```python
Z_L = j * omega * L
```

| Parametr | Jednotka | Default | Popis |
|----------|----------|---------|-------|
| L        | H        | 1e-6    | Induktance |

```python
L(1e-6)     # 1 uH
```

### Q - Constant Phase Element (CPE)

```python
Z_Q = 1 / (Q * (j * omega)^n)
```

| Parametr | Jednotka     | Default | Popis |
|----------|--------------|---------|-------|
| Q        | F * s^(n-1)  | 1e-4    | Q koeficient |
| n        | -            | 0.8     | Q exponent (0-1) |

Specialni pripady:
- n = 1: idealny kapacitor
- n = 0.5: Warburgova difuze
- n = 0: idealny rezistor

```python
Q(1e-4, 0.8)      # typicke Q
Q(1e-4, "0.9")    # Q volny, n fixovany
Q("1e-4", "0.9")  # oba fixovane
```

### W - Warburg (semi-infinite)

```python
Z_W = sigma / sqrt(omega) * (1 - j)
```

| Parametr | Jednotka     | Default | Popis |
|----------|--------------|---------|-------|
| sigma    | Ohm*s^(-1/2) | 50      | Warburg koeficient |

```python
W(50)       # sigma = 50
```

### Wo - Warburg Open (bounded)

```python
Z_Wo = R_W * tanh(sqrt(j*omega*tau)) / sqrt(j*omega*tau)
```

| Parametr | Jednotka | Default | Popis |
|----------|----------|---------|-------|
| R_W      | Ohm      | 100     | Warburg odpor |
| tau_W    | s        | 1.0     | Difuzni casova konstanta |

```python
Wo(100, 1.0)    # R_W=100, tau=1s
```

### K - Voigt element (R||C s tau parametrizaci)

```python
Z_K = R / (1 + j*omega*tau)
```

Ekvivalent (R || C) kde C = tau/R.

| Parametr | Jednotka | Default | Popis |
|----------|----------|---------|-------|
| R        | Ohm      | 1000    | Odpor |
| tau      | s        | 1e-4    | Casova konstanta |

Vyhody tau parametrizace:
- tau primo urcuje charakteristickou frekvenci: f = 1/(2*pi*tau)
- R a tau jsou nezavislejsi (lepsi numericka podminienost)
- Konzistentni s DRT notaci

```python
K(1000, 1e-4)   # R=1k, tau=100us, f=1.59kHz, C=100nF
```


## Fixovane parametry

Parametry predane jako **string** jsou automaticky fixovane pri fittingu:

```python
# Volne parametry (budou fitovany):
R(100)
Q(1e-4, 0.8)

# Fixovane parametry (nebudou fitovany):
R("0.86")           # R_inf z predchoziho mereni
Q("1e-4", 0.8)    # Q fixovane, n volne
Q("1e-4", "0.9")  # oba fixovane
```

Pouziti:
- Fixovani R_inf z vysokofrekvencniho mereni
- Fixovani geometrickych parametru
- Postupny fitting (nejdriv cast, pak celek)


## Operatory

### Serie (`-`)

```python
Z_total = Z1 + Z2 + ... + Zn
```

```python
R(100) - C(1e-6)                    # R a C v serii
R(10) - (R(100)|C(1e-6))            # R seriove s Voigtem
```

### Parallel (`|`)

```python
1/Z_total = 1/Z1 + 1/Z2 + ... + 1/Zn
```

```python
R(1000) | C(1e-6)                   # Voigt element
R(100) | Q(1e-4, 0.8)             # R||Q
```

### Skalovani (`*`)

```python
2 * R(100)      # = R(200)
0.5 * C(1e-6)   # = C(5e-7)
```

### Exponent (`**`) - pouze Q

```python
Q(1e-4) ** 0.9    # zmeni n na 0.9
```


## Priklady obvodu

### Jednoduchy Voigt element

```python
R(100) - (R(5000) | C(1e-6))
```

Struktura: R_s - (R_p || C_p)

### Randles circuit

```python
R(10) - ((R(100) - W(50)) | Q(1e-4, 0.8))
```

Struktura: R_s - ((R_ct - W) || Q_dl)

### Dva Voigt elementy (oxidova vrstva)

```python
R("1.69") - (R(4000) | Q(5e-8, 0.98)) - (R(1200) | Q(3e-8, 0.96))
```

Struktura: R_s(fix) - (R1 || Q1) - (R2 || Q2)

### K element (alternativa k Voigt)

```python
R(1) - K(1000, 1e-4) - K(500, 1e-3)
```


## Interni reprezentace

Po parsovani vznikne stromova struktura:

```python
circuit = R(100) - (R(5000) | C(1e-6))
print(type(circuit))  # <class 'Series'>
print(circuit)        # R(100) - (R(5000) | C(1e-6))
```

Interni struktura:
```
Series
├── R(100)
└── Parallel
    ├── R(5000)
    └── C(1e-6)
```

### Metody

```python
# Vsechny parametry (pro fitting)
circuit.get_all_params()        # [100, 5000, 1e-6]

# Ktere jsou fixovane
circuit.get_all_fixed_params()  # [False, False, False]

# Nazvy parametru
circuit.get_param_labels()      # ['R', 'R', 'C']

# Vypocet impedance
Z = circuit.impedance(frequencies, params)

# Update po fittingu
circuit.update_params(fitted_params)
```


## Omezeni soucasne implementace

### 1. Priorita operatoru

Jak je popsano vyse, `-` ma vyssi prioritu nez `|`. Reseni: pouzivat zavorky.

### 2. Zadne bounds v syntaxi

Bounds pro parametry se nedaji specifikovat v circuit stringu.
Pouzivaji se defaultni bounds podle typu elementu.

### 3. Zadna validace struktury

Parser nevaliduje, zda je obvod fyzikalne smysluplny.
Napriklad `C(1e-6) | C(1e-6)` je syntakticky spravne, ale nesmyslne.

### 4. Zalezitost na poradi

Poradi parametru v `get_all_params()` zavisi na poradi elementu
ve vyrazu. Pri zmene struktury se meni i indexy parametru.


## Bezpecnost

Parser pouziva `eval()` s omezenym namespace:

```python
eval(expr, {"__builtins__": {}}, safe_namespace)
```

- `__builtins__: {}` - zadne built-in funkce
- `safe_namespace` - pouze obvodove elementy

Toto zabranuje:
- Spusteni libovolneho kodu
- Pristupu k souborovemu systemu
- Import modulu

Presto se doporucuje pouzivat pouze duveryhodne vstupy.


## Mozna budouci rozsireni

### Vlastni parser

Implementace rekurzivniho sestupoveho parseru by umoznila:
- Spravnou prioritu operatoru (| > -)
- Validaci struktury
- Lepsí chybove hlasky
- Bounds v syntaxi: `R(100, bounds=(10, 1000))`

Gramatika:
```
expr   := term ('-' term)*
term   := factor ('|' factor)*
factor := element | '(' expr ')'
element := R(...) | C(...) | Q(...) | ...
```

Odhadovana slozitost: ~100 radku kodu.

### Alternativni operatory

Zmena operatoru na ty s opacnou prioritou v Pythonu:

| Novy | Stary | Priorita |
|------|-------|----------|
| `*`  | `\|`  | 13 (vyssi) |
| `+`  | `-`   | 12 (nizsi) |

Breaking change pro uzivatele.


## Reference

- Schonleber, M. et al. "A Method for Improving the Robustness of linear
  Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
- Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance
  Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
