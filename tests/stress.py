#!/usr/bin/env python3
"""Stress test: random, physically faithful spectra checked against invariants.

    python3 tests/stress.py                          # full run, all families
    python3 tests/stress.py --n 20 --family rc       # smaller run
    python3 tests/stress.py --family oxide --index 37 -v   # one case
    python3 tests/stress.py --json out.json          # also write every check

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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tests.stress_cases import FAMILIES, generate_case  # noqa: E402
from tests.stress_invariants import ANALYSES, INVARIANTS, check_case  # noqa: E402

# Cases per family of a full run; sized so that the run stays under ~1 h on
# 4 processes (measured ~4-5 s per case, see doc/STRESS_TEST_PLAN.md)
DEFAULT_N = 250

# At most 4 CPU-heavy processes for long runs
MAX_WORKERS = 4

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
    rows = check_case(case)
    return {
        'family': family, 'index': index, 'expression': case.expression,
        'n_arcs': case.n_arcs, 'min_sep': case.min_sep, 'min_frac': case.min_frac,
        'noise': f'{case.noise_level:g} {case.noise_kind}',
        'n_points': len(case.frequencies), 'n_rejected': case.n_rejected,
        'seconds': time.perf_counter() - start,
        'checks': [dict(zip(('analysis', 'invariant', 'status', 'detail'), r)) for r in rows],
    }


def _row_key(result):
    return f'{result["family"]}/{result["n_arcs"]}'


def print_tables(results):
    """One table per invariant: rows family/n_arcs, columns analyses, cells
    failures/total (meze in brackets)."""
    rows = sorted({_row_key(r) for r in results},
                  key=lambda k: (list(FAMILIES).index(k.split('/')[0]), k))
    for inv in INVARIANTS:
        counts = defaultdict(Counter)
        for r in results:
            for c in r['checks']:
                if c['invariant'] == inv:
                    counts[(_row_key(r), c['analysis'])][c['status']] += 1
        print(f'\nInvariant {inv}  (fail/checked, [meze], skip not counted)')
        print(f'{"":14}' + ''.join(f'{a:>14}' for a in ANALYSES))
        for row in rows:
            cells = []
            for a in ANALYSES:
                n = counts[(row, a)]
                checked = n['pass'] + n['fail'] + n['meze']
                cell = f'{n["fail"]}/{checked}' + (f' [{n["meze"]}]' if n['meze'] else '')
                cells.append(f'{cell:>14}')
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
    print('\nFit by identifiability class  (fail/checked, [meze])')
    print(f'{"":8}' + ''.join(f'{c:>12}' for c in columns))
    for inv in INVARIANTS:
        cells = []
        for cls in columns:
            n = counts[(inv, cls)]
            cell = f'{n["fail"]}/{n["pass"] + n["fail"] + n["meze"]}' + (f' [{n["meze"]}]' if n['meze'] else '')
            cells.append(f'{cell:>12}')
        print(f'{inv:8}' + ''.join(cells))


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
    args = parser.parse_args(argv)

    families = args.family or list(FAMILIES)
    if args.index is not None:
        if len(families) != 1:
            parser.error('--index needs exactly one --family')
        jobs = [(families[0], args.index)]
    else:
        jobs = [(family, i) for family in families for i in range(args.n)]

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
    print_failures(results)

    seconds = [r['seconds'] for r in results]
    print(f'\n{len(results)} cases in {elapsed:.0f} s '
          f'({sum(seconds) / len(seconds):.1f} s/case CPU, max {max(seconds):.1f} s); '
          f'rejected draws {sum(r["n_rejected"] for r in results)}')
    if args.json:
        args.json.write_text(json.dumps(results, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
