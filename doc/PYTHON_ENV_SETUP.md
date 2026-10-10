# Migrace na aktualni Python a knihovny (venv pres uv)

Provedeno 2026-10-10 na Orange Pi (Debian 12, aarch64): system Python 3.11.2
s numpy 1.24.2 / scipy 1.10.1 z apt -> venv s Pythonem 3.15.0, numpy 2.5.3,
scipy 1.18.1, matplotlib 3.11.2. Systemove balicky zustavaji, nic se nemaze.

## Postup (Linux)

```bash
# 1. uv (stahuje hotove buildy Pythonu, bez kompilace) do ~/.local/bin
curl -LsSf https://astral.sh/uv/install.sh | sh
uv --version

# 2. venv s nejnovejsim Pythonem
uv venv ~/venvs/eis --python 3.15

# 3. projekt (editovatelne) + dev nastroje + scikit-learn (GMM)
cd ~/cesta/k/eis_analysis
uv pip install --python ~/venvs/eis/bin/python -e ".[dev]" scikit-learn

# 4. prikaz `eis` z venv (stary skript schovat, pokud existuje)
[ -e ~/.local/bin/eis ] && mv ~/.local/bin/eis ~/.local/bin/eis.system-python
ln -s ~/venvs/eis/bin/eis ~/.local/bin/eis
hash -r
```

`~/.local/bin` musi byt v `PATH` (na Debianu/Ubuntu je, pokud adresar
existoval pri prihlaseni; jinak se odhlasit a prihlasit).

## Overeni

```bash
eis --version                                    # verze projektu, z venv
~/venvs/eis/bin/python -c "import tkinter"       # okna s grafy (Tk)
source ~/venvs/eis/bin/activate
python3 -m pytest tests/                         # vse passed
python3 -m pytest tests/ -m stress               # smoke proti baseline
```

## Pouzivani

- `eis data.DTA` funguje odkudkoli bez aktivace (skript ma natvrdo Python
  z venv). Zmeny v kodu se projevi hned (editovatelna instalace).
- Pro `pytest`, `mypy`, `ruff`, `python3 eis.py` nejdriv
  `source ~/venvs/eis/bin/activate` (opustit: `deactivate`).
- Dalsi balicek: `uv pip install --python ~/venvs/eis/bin/python <balicek>`.

## Pozor

- **Systemove `python3-numpy` / `python3-scipy` z apt nemazat**: zavisi na nich
  matplotlib, sympy a dalsi systemove balicky. Venv je zakryva, nevadi.
- **Ne `pip install --user --upgrade numpy`**: numpy 2.x v user site by rozbila
  matplotlib z apt (kompilovany proti numpy 1.x).
- Prubeh DE se s novou scipy lisi (jiny pocet iteraci), vysledne fity stejne.
  Stress baseline (`tests/stress_baseline.json`) overit `--check`.

## Navrat

```bash
rm ~/.local/bin/eis && mv ~/.local/bin/eis.system-python ~/.local/bin/eis
rm -rf ~/venvs/eis
```

## Windows (PowerShell)

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
uv venv $HOME\venvs\eis --python 3.15
cd C:\cesta\k\eis_analysis
uv pip install --python $HOME\venvs\eis\Scripts\python.exe -e ".[dev]" scikit-learn
# prikaz eis: $HOME\venvs\eis\Scripts\eis.exe (pridat Scripts do PATH,
# nebo aktivovat: $HOME\venvs\eis\Scripts\Activate.ps1)
```
