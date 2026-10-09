#!/usr/bin/env python3
"""Stress test: random, physically faithful spectra checked against invariants.

    python3 tests/stress.py                          # full run, all families
    python3 tests/stress.py --n 20 --family rc       # smaller run
    python3 tests/stress.py --family oxide --index 37 -v   # one case
    python3 tests/stress.py --json out.json          # also write every check
    python3 tests/stress.py --update-baseline        # full run -> stress_baseline.json
    python3 tests/stress.py --n 3 --check            # fail on a failure not in the baseline

Not a pytest module (no test_ prefix). Design and invariants:
doc/STRESS_TEST_PLAN.md.
"""

import os

# One BLAS thread per process, before numpy is imported: 4 processes x a
# multithreaded OpenBLAS oversubscribe the cores, and threaded reductions are
# not bit-identical, which would break invariant K and --index replay.
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_var] = '1'

import argparse  # noqa: E402
import json  # noqa: E402
import logging  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tests.stress_cases import FAMILIES, generate_case  # noqa: E402
from tests.stress_consistency import (CONSISTENCY_INVARIANTS, RATIO_INVARIANTS,  # noqa: E402
                                      consistency_rows)
from tests.stress_invariants import ANALYSES, INVARIANTS, check_case  # noqa: E402

# Cases per family of a full run; sized so that the run stays under ~1 h on
# 4 processes (measured ~4-5 s per case, see doc/STRESS_TEST_PLAN.md)
DEFAULT_N = 250

# At most 4 CPU-heavy processes for long runs
MAX_WORKERS = 4

# Known failures (family/index analysis:invariant) and aggregate rates of
# the last full run; --check compares against it (doc/STRESS_TEST_PLAN.md)
BASELINE = Path(__file__).with_name('stress_baseline.json')
# A rate fails --check when it moves the wrong way by more than this many
# binomial standard errors of the baseline rate (over the cases): a smaller
# move is within what redrawing the cases would do, so a change of code that
# shifts a few cases does not count as a regression.
RATE_SIGMAS = 3.0
# Rates where a rise is the regression (local minima); for the rest a drop is
RISING_RATES = ('F2',)

# Invariants with a fail/checked table; F3, Ifit and M are aggregate rates ('stat')
TABLE_INVARIANTS = INVARIANTS + tuple(i for i in CONSISTENCY_INVARIANTS if i not in ('F3', 'Ifit', 'M'))
FIT_INVARIANTS = INVARIANTS + ('F1', 'F2', 'F4', 'G')

# Deviations that are not failures, shown in brackets: absolute bounds
# (Babs) and fits ending in a local minimum (F2), a rate by class
SOFT = ('meze', 'lokmin')


def _cell(n: Counter) -> str:
    soft = sum(n[s] for s in SOFT)
    return f'{n["fail"]}/{n["pass"] + n["fail"] + soft}' + (f' [{soft}]' if soft else '')


# Quantiles of a consistency invariant's measured value, for calibration
QUANTILES = (0.5, 0.9, 0.99, 1.0)

# Identifiability classes (doc/STRESS_TEST_PLAN.md): neighbouring arcs closer
# than 1 decade are tight, over 2 decades loose; an arc under 5 % of the
# summed arc R is weak (lost in the noise of Re Z at 3 %)
TIGHT_SEP_DEC = 1.0
LOOSE_SEP_DEC = 2.0
WEAK_ARC_FRAC = 0.05


def _quiet():
    """Library warnings and log lines are UI; the invariants judge the results."""
    logging.disable(logging.CRITICAL)
    warnings.simplefilter('ignore')


def run_case(job):
    family, index = job
    _quiet()
    start = time.perf_counter()
    case = generate_case(family, index)
    rows, raw = check_case(case)
    rows += consistency_rows(case, raw)
    return {
        'family': family, 'index': index, 'expression': case.expression,
        'n_arcs': case.n_arcs, 'min_sep': case.min_sep, 'min_frac': case.min_frac,
        'noise': f'{case.noise_level:g} {case.noise_kind}',
        'n_points': len(case.frequencies), 'n_rejected': case.n_rejected,
        'seconds': time.perf_counter() - start,
        'checks': [dict(zip(('analysis', 'invariant', 'status', 'detail', 'value'), r))
                   for r in rows],
    }


def _row_key(result):
    return f'{result["family"]}/{result["n_arcs"]}'


def print_tables(results):
    """One table per invariant: rows family/n_arcs, columns analyses, cells
    failures/total (meze, lokmin in brackets)."""
    rows = sorted({_row_key(r) for r in results},
                  key=lambda k: (list(FAMILIES).index(k.split('/')[0]), k))
    for inv in TABLE_INVARIANTS:
        counts = defaultdict(Counter)
        for r in results:
            for c in r['checks']:
                if c['invariant'] == inv:
                    counts[(_row_key(r), c['analysis'])][c['status']] += 1
        print(f'\nInvariant {inv}  (fail/checked, [meze/lokmin], skip not counted)')
        print(f'{"":14}' + ''.join(f'{a:>14}' for a in ANALYSES))
        for row in rows:
            cells = []
            for a in ANALYSES:
                cells.append(f'{_cell(counts[(row, a)]):>14}')
            print(f'{row:14}' + ''.join(cells))


def identifiability_classes(result):
    """(separation class, smallest-arc class) of a case, doc/STRESS_TEST_PLAN.md."""
    if result['min_frac'] is None:
        return 'n/a', 'n/a'
    sep = result['min_sep']
    sep_class = ('jeden' if sep is None else 'tesny' if sep < TIGHT_SEP_DEC
                 else 'stredni' if sep <= LOOSE_SEP_DEC else 'volny')
    return sep_class, 'slaby' if result['min_frac'] < WEAK_ARC_FRAC else 'normalni'


def print_fit_classes(results):
    """Fit invariants by identifiability class: a failure in a tight or weak
    case is expected ill-posedness, in the other classes it is a bug."""
    columns = ('jeden', 'tesny', 'stredni', 'volny', 'slaby', 'normalni', 'n/a')
    counts = defaultdict(Counter)
    for r in results:
        classes = set(identifiability_classes(r))
        for c in r['checks']:
            if c['analysis'] == 'fit':
                for cls in classes:
                    counts[(c['invariant'], cls)][c['status']] += 1
    print('\nFit by identifiability class  (fail/checked, [meze/lokmin])')
    print(f'{"":8}' + ''.join(f'{c:>12}' for c in columns))
    for inv in FIT_INVARIANTS:
        cells = []
        for cls in columns:
            cells.append(f'{_cell(counts[(inv, cls)]):>12}')
        print(f'{inv:8}' + ''.join(cells))


def aggregate_rates(results):
    """'invariant noise-level' -> [hits, total, cases] of the rate invariants
    (F3, Ifit, M) and of F2's local minima. total counts parameters (F3) or
    points (M), which within a case are correlated; cases counts the cases."""
    rates = defaultdict(lambda: [0, 0, 0])
    for r in results:
        level = r['noise'].split()[0]
        for c in r['checks']:
            if c['status'] == 'stat':
                hit, total = map(int, c['detail'].split()[0].split('/'))
            elif c['invariant'] == 'F2' and c['status'] in ('pass', 'lokmin'):
                # not 'fail': an exception says nothing about local minima
                hit, total = int(c['status'] == 'lokmin'), 1
            else:
                continue
            rate = rates[f'{c["invariant"]} {level}']
            rate[0] += hit
            rate[1] += total
            rate[2] += 1
    return dict(sorted(rates.items()))


def print_calibration(results):
    """Distribution of measured / allowed of each threshold check by noise
    level (what its threshold is calibrated against), and the aggregate
    rates, with M's worst |dn| per case."""
    values = defaultdict(list)
    for r in results:
        level = r['noise'].split()[0]
        for c in r['checks']:
            if c.get('value') is None:
                continue
            if c['invariant'] in RATIO_INVARIANTS:
                values[(c['analysis'], c['invariant'], level)].append(c['value'])
            elif c['invariant'] == 'M':
                values[('n(f)', 'M |dn|', level)].append(c['value'])
    print('\nCalibration: measured / allowed by noise level, M as max |dn| (quantiles '
          + ', '.join(f'{q:g}' for q in QUANTILES) + ')')
    for (analysis, inv, level), v in sorted(values.items()):
        q = np.quantile(v, QUANTILES)
        print(f'  {analysis:6} {inv:8} noise {level:6} n={len(v):5}  '
              + '  '.join(f'{x:10.3g}' for x in q))
    print('\nAggregate rates (F2: local minima; F3: truth inside the 95 % CI;'
          ' Ifit: R_inf fitted on a closed HF end; M: |dn| <= 2 x uncertainty)')
    for key, (hit, total, _) in aggregate_rates(results).items():
        inv, level = key.split()
        print(f'  {inv:4} noise {level:6} {hit}/{total} = {hit / max(total, 1):.3f}')


def failure_keys(results):
    return {f'{r["family"]}/{r["index"]} {c["analysis"]}:{c["invariant"]}'
            for r in results for c in r['checks'] if c['status'] == 'fail'}


def write_baseline(results):
    BASELINE.write_text(json.dumps({
        'failures': sorted(failure_keys(results)),
        'rates': aggregate_rates(results)}, indent=1) + '\n')
    print(f'\nBaseline written to {BASELINE}')


def check_baseline(results, full_run):
    """False on a failure the baseline does not list, or, on a full run,
    a rate that moved the wrong way by more than RATE_SIGMAS."""
    base = json.loads(BASELINE.read_text())
    ran = {f'{r["family"]}/{r["index"]}' for r in results}
    expected = {k for k in base['failures'] if k.split()[0] in ran}
    got = failure_keys(results)
    new, gone = sorted(got - expected), sorted(expected - got)
    print(f'\nBaseline check: {len(new)} new failures, {len(gone)} gone')
    for k in new:
        print(f'  NEW  {k}')
    for k in gone:
        print(f'  gone {k}  (candidate for --update-baseline)')
    ok = not new
    if full_run:
        for key, (hit, total, cases) in aggregate_rates(results).items():
            if key not in base['rates'] or total == 0:
                continue
            b_hit, b_total, _ = base['rates'][key]
            p = b_hit / b_total
            # The standard error over cases, not points: the points of one
            # case move together. At least one case of slack, where the
            # baseline rate is 0 or 1.
            allowed = max(RATE_SIGMAS * np.sqrt(p * (1 - p) / cases), 1 / cases)
            worse = (hit / total - p) * (1 if key.split()[0] in RISING_RATES else -1)
            if worse > allowed:
                ok = False
                print(f'  RATE {key}: {hit / total:.3f} against {p:.3f} (allowed {allowed:.3f})')
    return ok


def print_failures(results):
    failures = [(r['family'], r['index'], c) for r in results for c in r['checks']
                if c['status'] == 'fail']
    print(f'\nFailures: {len(failures)}')
    for family, index, c in failures:
        print(f'  {family}/{index} {c["analysis"]}:{c["invariant"]}  {c["detail"]}')


def print_case(result):
    print(f'{result["family"]}/{result["index"]}: {result["expression"]}')
    print(f'  noise {result["noise"]}, {result["n_points"]} points, '
          f'n_arcs {result["n_arcs"]}, min_sep {result["min_sep"]}, '
          f'min_frac {result["min_frac"]}, rejected draws {result["n_rejected"]}')
    for c in result['checks']:
        print(f'  {c["analysis"]:6} {c["invariant"]:5} {c["status"]:5} {c["detail"]}')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--n', type=int, default=DEFAULT_N, help='cases per family')
    parser.add_argument('--family', action='append', choices=list(FAMILIES),
                        help='family to run (repeatable; default all)')
    parser.add_argument('--index', type=int, help='run only this case (with --family)')
    parser.add_argument('-v', '--verbose', action='store_true', help='print every check')
    parser.add_argument('--workers', type=int, default=MAX_WORKERS,
                        help=f'processes (default and maximum {MAX_WORKERS})')
    parser.add_argument('--json', type=Path, help='write every check to this file')
    # Exclusive: written first, the baseline would pass its own check
    baseline = parser.add_mutually_exclusive_group()
    baseline.add_argument('--update-baseline', action='store_true',
                          help=f'write {BASELINE.name} (full run only)')
    baseline.add_argument('--check', action='store_true',
                          help=f'exit 1 on a failure not in {BASELINE.name}')
    args = parser.parse_args(argv)

    families = list(dict.fromkeys(args.family or FAMILIES))   # a repeated --family runs once
    if args.index is not None:
        if len(families) != 1:
            parser.error('--index needs exactly one --family')
        jobs = [(families[0], args.index)]
    else:
        jobs = [(family, i) for family in families for i in range(args.n)]
    full_run = args.index is None and args.n == DEFAULT_N and set(families) == set(FAMILIES)
    if args.update_baseline and not full_run:
        parser.error('--update-baseline needs a full run (no --n, --family, --index)')
    if args.check and not BASELINE.exists():
        parser.error(f'no {BASELINE.name} yet: write it with a full run --update-baseline')

    start = time.perf_counter()
    if len(jobs) == 1:
        results = [run_case(jobs[0])]
    else:
        with ProcessPoolExecutor(max_workers=min(args.workers, MAX_WORKERS)) as pool:
            results = list(pool.map(run_case, jobs))
    elapsed = time.perf_counter() - start

    if args.verbose or len(results) == 1:
        for r in results:
            print_case(r)
    print_tables(results)
    print_fit_classes(results)
    print_calibration(results)
    print_failures(results)

    seconds = [r['seconds'] for r in results]
    print(f'\n{len(results)} cases in {elapsed:.0f} s '
          f'({sum(seconds) / len(seconds):.1f} s/case CPU, max {max(seconds):.1f} s); '
          f'rejected draws {sum(r["n_rejected"] for r in results)}')
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))
    if args.update_baseline:
        write_baseline(results)
    if args.check and not check_baseline(results, full_run):
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
