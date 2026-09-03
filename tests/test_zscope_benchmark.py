#!/usr/bin/env python3
"""Regression benchmark: parameter recovery on the four reference circuits
used by the ZScope benchmark suite.

The circuits, their ground-truth parameters and the frequency grid
(80 points, 1e-2 to 1e5 Hz, log-spaced) were taken from the published ZScope
benchmark (v2.2.0, `benchmarks/zview_test_data/`); see
doc/ZSCOPE_COMPARISON.md. The ground truth was verified by fitting their
noise-free CSVs: our circuits reproduce them to max|dZ|/|Z| ~ 5e-7, which is
the rounding of their 7-digit export.

Data is regenerated here rather than vendored. Their repository ships only a
CC BY 3.0 template license (for the HTML5UP website), so the license of the
CSVs themselves is unclear, and a test pointing at an untracked scratch
directory would skip everywhere except one machine. The circuits are what
carries the value, and they are public in the README.

Noise model: proportional Gaussian, sigma = level * |Z| on the real and the
imaginary part independently. This is the assumption behind `modulus`
weighting, which is the fitter's default. It is not a bit-exact clone of the
ZScope generator (theirs is somewhat narrower and not documented), so the
tolerances below are calibrated from our own runs, never copied from theirs.

Everything is deterministic: fixed seed for the noise, fixed seed for DE.
Runtime is roughly 20 s, dominated by DE converging on the noise-free cases.
"""

from functools import lru_cache

import numpy as np
import pytest

# Suppress matplotlib GUI (fit_circuit_diffevo builds a figure)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from eis_analysis.cli.utils import parse_circuit_expression
from eis_analysis.fitting.diffevo import fit_circuit_diffevo

# Excluded from the default run; see [tool.pytest.ini_options] in pyproject.
pytestmark = pytest.mark.slow

# Same grid as the ZScope benchmark data
FREQ = np.logspace(-2, 5, 80)

# name -> (circuit expression with ground-truth values, parameter labels)
CIRCUITS = {
    'Randles': (
        'R(10)-(R(200)|C(2e-5))',
        ['Rs', 'Rct', 'Cdl'],
    ),
    'Randles_Warburg': (
        'R(10)-((R(150)-W(80))|C(1.5e-5))',
        ['Rs', 'Rct', 'sigma_W', 'Cdl'],
    ),
    'CPE_Randles': (
        'R(10)-(R(200)|Q(2e-5,0.85))',
        ['Rs', 'Rct', 'Q', 'alpha'],
    ),
    'Two_TimeConstants': (
        'R(10)-(R(200)|C(1e-5))-(R(100)|C(1e-3))',
        ['Rs', 'R1', 'C1', 'R2', 'C2'],
    ),
}

# Max relative parameter error [%] tolerated at each noise level.
# Calibrated from the seed-0 run with roughly 2x headroom; the CPE case is
# looser because Q is strongly correlated with alpha, so noise moves it far
# more than the other parameters (Q alone reaches 4.7% at 2% and 10.9% at 5%,
# while Rs, Rct and alpha all stay under 2%). ZScope reports the same
# Q-alpha correlation as a physical interdependence, not a solver defect.
TOLERANCE_PCT = {
    ('Randles', 0.02): 2.0,
    ('Randles', 0.05): 4.0,
    ('Randles_Warburg', 0.02): 2.0,
    ('Randles_Warburg', 0.05): 4.0,
    ('CPE_Randles', 0.02): 8.0,
    ('CPE_Randles', 0.05): 18.0,
    ('Two_TimeConstants', 0.02): 2.0,
    ('Two_TimeConstants', 0.05): 4.0,
}

NOISE_LEVELS = [0.02, 0.05]


def _canonical(name, params):
    """Order the two RC branches of Two_TimeConstants by time constant.

    A series chain of Voigt elements is invariant under permutation of the
    branches, so the fitter is free to return them in either order (and does:
    on the ZScope 2% data it swapped R1<->R2, C1<->C2 at an otherwise correct
    RMSE of 2.0%). Sorting by tau = R*C makes the comparison well-defined
    instead of testing an arbitrary labelling.
    """
    p = list(params)
    if name != 'Two_TimeConstants':
        return p
    branches = sorted([(p[1] * p[2], p[1], p[2]), (p[3] * p[4], p[3], p[4])])
    return [p[0], branches[0][1], branches[0][2],
            branches[1][1], branches[1][2]]


def _truth(name):
    circuit = parse_circuit_expression(CIRCUITS[name][0])
    return _canonical(name, circuit.get_all_params())


def _clean_impedance(name):
    circuit = parse_circuit_expression(CIRCUITS[name][0])
    return circuit.impedance(FREQ, circuit.get_all_params())


@lru_cache(maxsize=None)
def _fit(name, noise):
    """Fit `name` at the given noise level. Cached so the whole module needs
    only one DE run per (circuit, noise) combination."""
    Z = _clean_impedance(name)
    if noise:
        rng = np.random.default_rng(0)
        sigma = noise * np.abs(Z)
        Z = Z + rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)

    circuit = parse_circuit_expression(CIRCUITS[name][0])
    result, _, _ = fit_circuit_diffevo(circuit, FREQ, Z, seed=0)
    plt.close('all')
    return result.best_result


def _errors_pct(name, fitted):
    truth = np.asarray(_truth(name), dtype=float)
    params = np.asarray(_canonical(name, fitted.params_opt), dtype=float)
    return 100.0 * np.abs(params - truth) / np.abs(truth)


def _report(name, errors):
    return ', '.join(f'{label}={err:.2f}%'
                     for label, err in zip(CIRCUITS[name][1], errors))


@pytest.mark.parametrize('name', list(CIRCUITS))
def test_exact_recovery_from_noise_free_data(name):
    """On noise-free data the optimizer must return the generating parameters.

    This is the strongest regression signal in the file: it has no
    statistical slack, so any change that degrades DE, the analytic Jacobian
    or an element's impedance shows up here immediately.
    """
    errors = _errors_pct(name, _fit(name, 0.0))
    assert np.max(errors) < 0.1, f'{name}: {_report(name, errors)}'


@pytest.mark.parametrize('noise', NOISE_LEVELS)
@pytest.mark.parametrize('name', list(CIRCUITS))
def test_parameter_recovery_with_noise(name, noise):
    """Parameters stay within a calibrated band under proportional noise."""
    errors = _errors_pct(name, _fit(name, noise))
    tolerance = TOLERANCE_PCT[(name, noise)]
    assert np.max(errors) < tolerance, (
        f'{name} at {noise:.0%} noise: {_report(name, errors)} '
        f'(tolerance {tolerance}%)')


@pytest.mark.parametrize('noise', NOISE_LEVELS)
@pytest.mark.parametrize('name', list(CIRCUITS))
def test_fit_error_reaches_noise_floor(name, noise):
    """The fit must converge to the noise floor - no better, no worse.

    Below the floor the model is following the noise; well above it the
    optimizer stopped early or found a local minimum. Measured values sit at
    1.19-1.20x the noise level for every circuit, which is the expected ratio
    for this noise model: the relative fit error is computed against |Z|
    while the noise enters the real and imaginary part independently.
    """
    fit_error = _fit(name, noise).fit_error_rel
    assert 100.0 * noise <= fit_error < 100.0 * noise * 1.4, (
        f'{name} at {noise:.0%} noise: fit error {fit_error:.2f}%')
