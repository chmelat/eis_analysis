#!/usr/bin/env python3
"""Reference spectra with known physics: driftlet's Randles electrode.

The spectra in `driftlet_randles.json` are not made from an equivalent
circuit. driftlet (a drift-diffusion-Poisson solver in JavaScript,
https://github.com/markblundeberg/driftlet, 0BSD) simulates Pt in 1 mM Fe3+/Fe2+ with KCl at rest,
and its Randles page computes each circuit element from the physics alone:
R_s from the conductivity, C_dl as Stern in series with Gouy-Chapman, R_ct
from k0, and a finite-length (transmissive) Warburg per iron ion. A fit that
recovers those values therefore is not checked against its own model, unlike
the synthetic data elsewhere in the suite.

The circuit is an approximation of the simulation, not its exact solution:
driftlet states it within 1 % of |Z| across its sliders (migration, the
iron's share in charging the double layer), measured here as 1e-5 to 9e-3.
The two Warburgs (tau 1.19x apart) are fitted as one Ws, since two are not
separable (DE splits them arbitrarily); its R_W is compared with their sum
and its tau must lie between theirs.

Data: six parameter sets, the corners and middle checked by driftlet's
test/demos.test.js plus the page's default, 10 mHz to 1 MHz at 10 points per
decade, Ohm m^2 scaled to 1 cm^2. Regenerate with `tests/driftlet_randles.mjs`.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting.diffevo import fit_circuit_diffevo
from eis_analysis.validation import kramers_kronig_validation

DATA = json.loads(Path(__file__).with_name('driftlet_randles.json').read_text())
FREQ = np.array(DATA['freq'])
# Electrode area that turns driftlet's Ohm m^2 into Ohm: 1 cm^2, so that the
# values sit in the range of a measured cell (R_s 0.05-5 Ohm, C_dl ~15 uF).
AREA = 1e-4

# Max relative parameter error [%] against the physics. The 1-Ws fits deviate
# by at most 0.72 % (R_s at 10 mM KCl, where migration is largest), which is
# the circuit's own error against the simulation; 1.5 % keeps 2x headroom.
TOLERANCE_PCT = 1.5
# Relative fit error [%] the circuit must reach: measured 3e-4 to 0.11 %, the
# simulation's departure from the Randles circuit (driftlet's own estimate of
# its Z error, relError, was at most 3e-4 of |Z| when the data were made).
MAX_FIT_ERROR_PCT = 0.2
# Below this R_W / R_ct the Warburg is a 0.1 %-scale tail behind the charge
# transfer resistance and not identifiable: the fit trades R_ct for a short-tau
# Ws (R_ct -> 0, same DC sum). Such a set checks only R_ct + R_W. The six sets
# have R_W / R_ct of 1e-3 (unresolved) and 0.096 to 960 (all resolved).
WARBURG_RESOLVED_RATIO = 0.01

CASES = DATA['spectra']
IDS = [f"k0={s['params']['k0']:g},salt={s['params']['salt']:g},L={s['params']['L']:.3g}"
       for s in CASES]


def _truth(case):
    ph = case['physics']
    return {'Rs': ph['Rs'] / AREA, 'Cdl': ph['Cdl'] * AREA, 'Rct': ph['Rct'] / AREA,
            'RW': ph['RW'] / AREA}


def _simulated(case):
    return (np.array(case['Z_re']) + 1j * np.array(case['Z_im'])) / AREA


def _fit(case):
    """DE fit of R-(C|(R-Ws)) to one simulated spectrum, started from nothing
    but the bounds (the expression's values only name the elements)."""
    circuit = parse_circuit_expression('R(1)-(C(1e-5)|(R(100)-Ws(100,1)))')
    result, _ = fit_circuit_diffevo(FREQ, _simulated(case), circuit, seed=0)
    return result.best_result


@pytest.mark.parametrize('case', CASES, ids=IDS)
def test_our_circuit_is_the_pages_circuit(case):
    """With the physics' values, our R-(C|(R-Ws-Ws)) is the page's Randles
    impedance: same Ws (tanh) form, units and sign of Im Z."""
    t = _truth(case)
    w3, w2 = case['physics']['Ws_Fe3'], case['physics']['Ws_Fe2']
    circuit = parse_circuit_expression(
        f"R({t['Rs']})-(C({t['Cdl']})|(R({t['Rct']})"
        f"-Ws({w3['R_W'] / AREA},{w3['tau_W']})-Ws({w2['R_W'] / AREA},{w2['tau_W']})))")
    Z = circuit.impedance(FREQ, circuit.get_all_params())
    page = (np.array(case['circuit_re']) + 1j * np.array(case['circuit_im'])) / AREA
    assert np.max(np.abs(Z - page) / np.abs(page)) < 1e-12


@pytest.mark.parametrize('case', CASES, ids=IDS)
def test_simulated_spectrum_passes_kramers_kronig(case):
    """A simulated linear, causal, stable system is KK-consistent."""
    result = kramers_kronig_validation(FREQ, _simulated(case))
    assert result.error is None and result.is_valid


@pytest.mark.slow
@pytest.mark.parametrize('case', CASES, ids=IDS)
def test_fit_recovers_the_physics(case):
    """DE on the simulated spectrum returns the elements the physics gives."""
    t = _truth(case)
    fit = _fit(case)
    Rs, Cdl, Rct, RW, tau = np.asarray(fit.params_opt, dtype=float)
    assert fit.fit_error_rel < MAX_FIT_ERROR_PCT, f'fit error {fit.fit_error_rel:.3g} %'

    errors = {'Rs': (Rs, t['Rs']), 'Cdl': (Cdl, t['Cdl'])}
    if t['RW'] / t['Rct'] >= WARBURG_RESOLVED_RATIO:
        errors.update(Rct=(Rct, t['Rct']), RW=(RW, t['RW']))
        taus = sorted(case['physics'][k]['tau_W'] for k in ('Ws_Fe3', 'Ws_Fe2'))
        assert taus[0] <= tau <= taus[1], f'tau {tau:.4g} outside {taus}'
    else:
        errors['Rct+RW'] = (Rct + RW, t['Rct'] + t['RW'])
    pct = {name: 100 * abs(got / want - 1) for name, (got, want) in errors.items()}
    assert max(pct.values()) < TOLERANCE_PCT, ', '.join(f'{n} {e:.3f} %' for n, e in pct.items())
