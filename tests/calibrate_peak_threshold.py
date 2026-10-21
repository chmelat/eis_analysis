#!/usr/bin/env python3
"""Calibration of DRT_PEAK_DCHI2_MIN (drt.significance) on synthetic spectra.

    python3 tests/calibrate_peak_threshold.py               # 248 runs, 4 processes
    python3 tests/calibrate_peak_threshold.py --seeds 2     # quick look

Eight families x (noise 0 once + 0.1 / 1 / 3 % x 10 seeds) = 248 runs of
calculate_drt with its defaults (auto lambda, sqrt weighting, noise = the
DRT's own residual). The DRT is computed once per run; the shoulder test is
then repeated for every threshold, since the threshold decides which smaller
maxima are rejected and with them the runs of the others (not only a cut on
one Delta chi^2 per candidate).

A significant peak within MATCH_DECADES of a true tau (one to one, nearest
first) finds that process; the rest are false peaks. The margin is the
largest finite Delta chi^2 of a candidate matched to no process.

Not a pytest module (no test_ prefix).
"""

import os

# One BLAS thread per process (see tests/stress.py)
if __name__ == '__main__':
    for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[_var] = '1'

import argparse  # noqa: E402
import logging  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from collections import defaultdict  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import eis_analysis.drt.core as drt_core  # noqa: E402
from eis_analysis.drt.significance import assess_peak_significance  # noqa: E402

FREQUENCIES = np.logspace(5, -2, 71)
OMEGA = 2 * np.pi * FREQUENCIES
NOISE_LEVELS = (0.0, 1e-3, 1e-2, 3e-2)
THRESHOLDS = (4.0, 9.0, 16.0, 25.0, 36.0)
# Half the spacing of the 2xZARC 1 decade pair: a single maximum anywhere
# between two processes finds one of them (one to one), and is not counted
# as false; a resolved peak shifts by < 0.1 decade
MATCH_DECADES = 0.5
# At most 4 CPU-heavy processes for long runs
MAX_WORKERS = 4


def _zarc(R, tau, n=1.0):
    return R / (1 + (1j * OMEGA * tau) ** n)


# name -> (series R, [(R, tau, n)]); the RC taus lie inside the window
# (1.6 us .. 16 s), 3RC is the spectrum of finding F1 (doc/CONTRACT_REVIEW.md)
FAMILIES = {
    '1RC': (10, [(100, 1e-3, 1.0)]),
    '1ZARC n0.8': (10, [(100, 1e-3, 0.8)]),
    '2RC 1:10': (10, [(100, 1e-4, 1.0), (1000, 1e-1, 1.0)]),
    '2RC 1:100': (10, [(100, 1e-4, 1.0), (1e4, 1e-1, 1.0)]),
    '2RC 1:1000': (10, [(100, 1e-4, 1.0), (1e5, 1e-1, 1.0)]),
    '2ZARC 0.5dec': (10, [(100, 1e-3, 0.9), (100, 10 ** -2.5, 0.9)]),
    '2ZARC 1dec': (10, [(100, 1e-3, 0.9), (100, 1e-2, 0.9)]),
    '3RC': (10, [(100, 1e-4, 1.0), (5, 2e-3, 1.0), (3, 0.3, 1.0)]),
}


def _match(true_taus, taus):
    """Indices into taus matched one to one to true_taus, nearest first."""
    pairs = sorted((abs(np.log10(t / s)), i, k)
                   for i, t in enumerate(true_taus) for k, s in enumerate(taus))
    used_true, used = set(), set()
    for dist, i, k in pairs:
        if dist <= MATCH_DECADES and i not in used_true and k not in used:
            used_true.add(i)
            used.add(k)
    return used


def run_case(job):
    family, noise, seed = job
    R_s, elements = FAMILIES[family]
    Z = R_s + sum(_zarc(R, tau, n) for R, tau, n in elements)
    e = np.random.default_rng(seed).standard_normal((2, len(Z)))
    Z = Z + noise * np.abs(Z) * (e[0] + 1j * e[1])

    captured = {}

    def capture(*args, **kwargs):
        captured['args'] = args
        return assess_peak_significance(*args, **kwargs)

    drt_core.assess_peak_significance = capture
    logging.disable(logging.WARNING)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        drt_core.calculate_drt(FREQUENCIES, Z)
    matrices = captured['args'][0]
    true_taus = [tau for _, tau, _ in elements]

    rows = {}
    for threshold in THRESHOLDS:
        sig = assess_peak_significance(*captured['args'], threshold=threshold)
        taus = matrices.tau[sig.candidates]
        matched = _match(true_taus, taus)
        significant = set(sig.significant.tolist())
        found = _match(true_taus, [matrices.tau[c] for c in sig.significant])
        spurious = [d for k, d in enumerate(sig.delta_chi2)
                    if np.isfinite(d) and k not in matched]
        rows[threshold] = {
            'missed': len(true_taus) - len(found),
            'false': len(significant) - len(found),
            'max_spurious': max(spurious, default=-np.inf),
        }
    return family, noise, seed, rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--seeds', type=int, default=10, help='seeds per noise level > 0')
    args = parser.parse_args(argv)

    jobs = [(f, noise, seed) for f in FAMILIES for noise in NOISE_LEVELS
            for seed in (range(args.seeds) if noise > 0 else [0])]
    with ProcessPoolExecutor(MAX_WORKERS) as pool:
        results = list(pool.map(run_case, jobs))
    print(f'{len(results)} runs')

    print('\nthreshold  missed  false  max spurious dchi2 (case)')
    for threshold in THRESHOLDS:
        family, noise, seed, rows = max(results, key=lambda r: r[3][threshold]['max_spurious'])
        print(f'{threshold:9g}  {sum(r[threshold]["missed"] for *_, r in results):6d}  '
              f'{sum(r[threshold]["false"] for *_, r in results):5d}  '
              f'{rows[threshold]["max_spurious"]:18.3g} ({family}, {noise:.1%}, seed {seed})')

    print('\nmissed / false per family and noise (threshold columns)')
    cells = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for family, noise, _, rows in results:
        for threshold, r in rows.items():
            cells[(family, noise)][threshold][0] += r['missed']
            cells[(family, noise)][threshold][1] += r['false']
    print(f'{"family":14s} {"noise":>6s}  ' + '  '.join(f'{t:>7g}' for t in THRESHOLDS))
    for (family, noise), by_t in cells.items():
        if any(m or f for m, f in by_t.values()):
            print(f'{family:14s} {noise:6.1%}  '
                  + '  '.join(f'{m:>3d}/{f:<3d}' for m, f in (by_t[t] for t in THRESHOLDS)))


if __name__ == '__main__':
    main()
