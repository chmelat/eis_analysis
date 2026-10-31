#!/usr/bin/env python3
"""Fingerprint of the analyses' results, for a change meant to change nothing.

    python3 tests/fingerprint.py record ../before.pkl         # before the refactor
    python3 tests/fingerprint.py record ../after.pkl          # after it
    python3 tests/fingerprint.py compare ../before.pkl ../after.pkl

    # a committed refactor: the library of each side from a worktree
    git worktree add ../old <commit>^ && git worktree add ../new <commit>
    python3 tests/fingerprint.py record ../old.pkl --root ../old
    python3 tests/fingerprint.py record ../new.pkl --root ../new

Runs the stress test's analyses (tests/stress_invariants.py: DRT, Lin-KK,
Z-HIT, R_inf, LM fit, n(f)) once on the base spectrum of its first --n cases
per family, and the DE fit of invariant G (with its 95 % CI) on every
DE_EVERY-th case. It stores the values of their signatures, which carry the
stress test's tolerances, and the flag and text fields of each result
(is_valid, warnings, methods). Not covered: the LM fit's warnings and CI,
and whatever the CLI adds.

`--root` takes only `eis_analysis` from another checkout; the case generator
and the signatures are this checkout's code (as driftlet's
tools/fingerprint.js runs an older commit with the current tool). The
generator calls the library (circuits, bounds, impedance, fit starts), so
each case stores a hash of its input, and a case built differently is
reported as `input`, not as changed analyses. An old checkout that lacks a
library name the stress modules import cannot be recorded (record says so).

`compare` sorts each case and analysis into identical (every value equal,
NaN = NaN, -0.0 = 0.0), within the tolerances (with the largest relative
change) or beyond them, and exits 0 only when all are identical. Identity
holds only on the same Python, numpy, scipy, BLAS and machine; compare warns
when the recorded ones differ.

Records are pickles (multi-MB): keep them outside the repository, and load
only your own (pickle runs code at load).

Not a pytest module (no test_ prefix); its compare logic is checked by
tests/test_fingerprint.py. Idea: tools/fingerprint.js of driftlet,
https://github.com/markblundeberg/driftlet
"""

import argparse
import contextlib
import dataclasses
import functools
import hashlib
import inspect
import logging
import multiprocessing
import os
import pickle
import platform
import subprocess
import sys
import tempfile
import time
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
# A directory holding only an `eis_analysis` link to --root's package, put
# first on the workers' sys.path (they read it from the environment)
LIB_VAR = 'EIS_FINGERPRINT_LIB'

# Cases per family of a default record: 300 cases and 60 DE fits, under
# 3 min on 4 processes (RK3588); --n 250 is the stress test's full set
DEFAULT_N = 50
# Differences listed by compare after the table
SHOW = 20
# Hex digits kept of an input's SHA-1: 64 bits, so a chance collision among a
# record's at most 1500 cases is ~1e-13
HASH_DIGITS = 16
# A DE CI end may move by this share of the interval's width, as the stderr
# it is built from may move by STDERR_STDERR_SHARE of itself
CI_WIDTH_SHARE = 1e-2

CLASSES = ('identical', 'within', 'beyond', 'exception', 'missing', 'input')

# Fields of tests.stress_invariants.Entry, stored by name so a file does not
# depend on the Entry class (or its field order) of the checkout that wrote it
FIELDS = ('value', 'power', 'kind', 'finite', 'atol')


def _setup_path():
    """This checkout for `tests`, --root's package (if any) first for
    `eis_analysis`; called again after imports that put the repo first
    (tests/stress.py does)."""
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    lib = os.environ.get(LIB_VAR)
    if lib:
        if lib in sys.path:
            sys.path.remove(lib)
        sys.path.insert(0, lib)


def _entry(value, kind='exact', atol=None):
    return {'value': value, 'power': 0, 'kind': kind, 'finite': False, 'atol': atol}


def _flags(result):
    """Public flag and text fields and properties of a result object: the
    numbers are in the signatures, with their tolerances; these are compared
    with ==."""
    names = [f.name for f in dataclasses.fields(result)] if dataclasses.is_dataclass(result) else []
    names += [n for n, _ in inspect.getmembers(type(result), lambda v: isinstance(v, property))]
    out = {}
    for name in dict.fromkeys(n for n in names if not n.startswith('_')):
        try:
            value = getattr(result, name)
        except Exception as e:
            value = f'raises {type(e).__name__}'
        if isinstance(value, (np.bool_, np.integer)):
            value = value.item()
        if isinstance(value, (bool, int, str, type(None))) or (
                isinstance(value, (list, tuple)) and all(isinstance(x, (bool, str)) for x in value)):
            out[f'raw.{name}'] = _entry(value)
    return out


def _plain(sig, result=None):
    """A signature as dicts of the Entry fields, plus the flags of its 'raw'
    result (or of `result`)."""
    if isinstance(sig, Exception):
        return ('exception', type(sig).__name__, str(sig))
    out = {name: {f: getattr(e, f) for f in FIELDS} for name, e in sig.items() if e.kind != 'raw'}
    raw = sig['raw'].value if 'raw' in sig else result
    if raw is not None:
        out.update(_flags(raw))
    return out


def record_case(job):
    """(job, {analysis: plain signature}, the package that ran it)."""
    _setup_path()
    import eis_analysis   # before the tests modules, as in record
    logging.disable(logging.CRITICAL)   # library output is UI, the values are judged
    warnings.simplefilter('ignore')
    from tests.stress_cases import capacitive_end, de_start, fit_start, generate_case, k_powers
    from tests.stress_consistency import DE_EVERY, run_de
    from tests.stress_invariants import ANALYSES, fit_signature, run_analysis
    package = str(Path(eis_analysis.__file__).resolve().parent)
    try:
        case = generate_case(*job)
        # What the analyses are given, the library-built parts included
        h = hashlib.sha1(f'{case.expression} {capacitive_end(case)}'.encode())
        for array in (case.truth, case.frequencies, case.Z, fit_start(case), de_start(case)):
            h.update(np.ascontiguousarray(array).tobytes())
    except Exception as e:   # the rest of the record is kept
        return job, {'case': _plain(e)}, package
    out = {'_input': h.hexdigest()[:HASH_DIGITS]}
    out.update({name: _plain(run_analysis(name, case, case.frequencies, case.Z))
                for name in ANALYSES})
    if job[1] % DE_EVERY == 0:
        try:
            r = run_de(case)[0].best_result
            out['de'] = _plain(fit_signature(r, k_powers(case.labels)), r)
            low, high = r.params_ci_95
            atol = CI_WIDTH_SHARE * np.where(np.isfinite(high - low), high - low, 0.0)
            out['de'].update({'ci95_low': _entry(low, 'rel', atol),
                              'ci95_high': _entry(high, 'rel', atol)})
        except Exception as e:
            out['de'] = _plain(e)
    return job, out, package


def _git(root, *args):
    try:
        return subprocess.run(['git', '-C', root, *args], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def record(out: Path, n: int, root: Path) -> int:
    # Checked before the run: the file is written only after minutes of it
    if out.exists():
        # An old file passed off as a new record would compare on its own data
        print(f'record: {out} exists; give a new file')
        return 2
    if not out.parent.is_dir():
        print(f'record: no directory {out.parent}')
        return 2
    if not (root / 'eis_analysis').is_dir():
        print(f'record: no eis_analysis package in {root}')
        return 2
    with tempfile.TemporaryDirectory() if root != REPO else contextlib.nullcontext() as lib:
        if lib:
            (Path(lib) / 'eis_analysis').symlink_to(root / 'eis_analysis')
            os.environ[LIB_VAR] = lib
        _setup_path()
        # First: tests/stress.py puts the repo first on sys.path and imports the library
        import eis_analysis
        package = Path(eis_analysis.__file__).resolve().parent
        try:
            from tests.stress import BLAS_ENV, MAX_WORKERS
            from tests.stress_cases import FAMILIES
        except ImportError as e:
            print(f'record: the stress modules of this checkout need a name {package} lacks: {e}')
            return 2
        if not package.is_relative_to(root):
            # e.g. an editable install of another checkout found first
            print(f'record: eis_analysis imported from {package}, not from {root}')
            return 2
        from eis_analysis.version import __version__
        import scipy

        # The stress test's BLAS settings (one thread, one kernel) and one hash
        # seed; spawned workers read them at start, and this process computes
        # nothing
        os.environ.update(BLAS_ENV, PYTHONHASHSEED='0')
        start = time.perf_counter()
        with ProcessPoolExecutor(max_workers=MAX_WORKERS,
                                 mp_context=multiprocessing.get_context('spawn')) as pool:
            done = list(pool.map(record_case, [(f, i) for f in FAMILIES for i in range(n)]))
    if {pkg for _, _, pkg in done} != {str(package)}:
        print(f'record: workers ran {sorted({pkg for _, _, pkg in done})}, not {package}; '
              'nothing written')
        return 2
    cases = {job: result for job, result, _ in done}
    try:
        blas = np.show_config(mode='dicts')['Build Dependencies']['blas']
        blas = f"{blas.get('name')} {blas.get('version')}"
    except Exception:   # numpy < 1.26 or another layout
        blas = None
    meta = {'commit': _git(str(root), 'describe', '--always', '--dirty'), 'version': __version__,
            'python': platform.python_version(), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'machine': platform.machine(), 'blas': blas,
            'n': n, 'package': str(package)}
    out.write_bytes(pickle.dumps({'meta': meta, 'cases': cases}))
    print(f'{len(cases)} cases, {sum("de" in c for c in cases.values())} DE fits in '
          f'{time.perf_counter() - start:.0f} s from {package} ({meta["commit"]}, '
          f'{__version__}) -> {out}')
    return 0


@functools.cache
def _invariants():
    """Entry and compare of this checkout (compare needs no --root)."""
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    from tests import stress_invariants
    return stress_invariants


def _max_rel(ref, got):
    """(entry, largest relative change) over the numeric entries: element by
    element for 'rel' (parameters spanning decades), against the entry's
    max |value| otherwise (arrays with legitimate zeros)."""
    worst = (None, 0.0)
    for name, entry in ref.items():
        if entry.kind == 'exact' or entry.value is None or got[name].value is None:
            continue
        a, b = (np.atleast_1d(np.asarray(e.value, dtype=complex)) for e in (entry, got[name]))
        if a.shape != b.shape:
            continue
        both = np.isfinite(a) & np.isfinite(b)
        a, b = a[both], b[both]
        if a.size == 0:
            continue
        scale = np.abs(a) if entry.kind == 'rel' else np.full(a.shape, np.max(np.abs(a)))
        change = np.abs(b - a)
        rel = np.where(scale > 0, change / np.where(scale > 0, scale, 1.0), change)
        if np.max(rel) > worst[1]:
            worst = (name, float(np.max(rel)))
    return worst


def classify(ref, got):
    """(class, detail) of one case and analysis; ref and got are plain
    signatures, exceptions or None (not recorded)."""
    if ref is None or got is None:
        return 'missing', 'only in ' + ('the second' if ref is None else 'the first')
    if isinstance(ref, tuple) or isinstance(got, tuple):
        if ref == got:
            return 'identical', ''
        return 'exception', f'{ref[:2] if isinstance(ref, tuple) else "result"} -> ' \
                            f'{got[:2] if isinstance(got, tuple) else "result"}'
    if ref.keys() != got.keys():
        return 'beyond', f'entries {sorted(ref.keys() ^ got.keys())} on one side only'
    inv = _invariants()
    a, b = ({k: inv.Entry(**v) for k, v in s.items()} for s in (ref, got))
    if not inv.compare(a, b, bitwise=True):
        return 'identical', ''
    beyond = inv.compare(a, b)
    if beyond:
        return 'beyond', '; '.join(beyond)
    name, rel = _max_rel(a, b)
    return 'within', f'{name}: rel {rel:.1e}'


def compare_files(path_a: Path, path_b: Path) -> int:
    A, B = (pickle.loads(p.read_bytes()) for p in (path_a, path_b))
    for side, m in (('first', A['meta']), ('second', B['meta'])):
        print(f'{side}: {m["package"]} ({m["commit"]}, {m["version"]}), n {m["n"]}, '
              f'Python {m["python"]}, numpy {m["numpy"]}, scipy {m["scipy"]}, '
              f'{m.get("machine")}, {m.get("blas")}')
    other = [k for k in ('python', 'numpy', 'scipy', 'machine', 'blas')
             if A['meta'].get(k) != B['meta'].get(k)]
    if other:
        print(f'warning: other {", ".join(other)}: identical values are not expected, '
              'only the tolerances')
    if A['meta']['n'] != B['meta']['n']:
        print('note: other --n: the cases beyond the smaller one are missing on one side')

    counts = {}
    shown = []
    for key in sorted(A['cases'].keys() | B['cases'].keys()):
        a, b = A['cases'].get(key, {}), B['cases'].get(key, {})
        if '_input' in a and '_input' in b and a['_input'] != b['_input']:
            # Another spectrum, circuit or start: its analyses cannot be compared
            counts.setdefault('(case)', Counter())['input'] += 1
            shown.append(f'  {key[0]}/{key[1]} input differs: the library built another '
                         'case (circuit, bounds, impedance, fit start)')
            continue
        for analysis in sorted((a.keys() | b.keys()) - {'_input'}):
            cls, detail = classify(a.get(analysis), b.get(analysis))
            counts.setdefault(analysis, Counter())[cls] += 1
            if cls != 'identical':
                shown.append(f'  {key[0]}/{key[1]} {analysis} {cls}: {detail}')

    print(f'\n{"analysis":10s}' + ''.join(f'{c:>11s}' for c in CLASSES))
    for analysis, c in counts.items():
        print(f'{analysis:10s}' + ''.join(f'{c[k]:11d}' for k in CLASSES))
    if shown:
        print(f'\n{len(shown)} differences' + (f', first {SHOW}:' if len(shown) > SHOW else ':'))
        print('\n'.join(shown[:SHOW]))
        return 1
    print('\nall identical')
    return 0


def _positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f'needs at least 1 case per family, got {value}')
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    rec = sub.add_parser('record', help='run the analyses and store their values')
    rec.add_argument('out', type=Path)
    rec.add_argument('--n', type=_positive_int, default=DEFAULT_N, help='cases per family')
    rec.add_argument('--root', type=Path, default=REPO,
                     help='checkout whose eis_analysis is recorded (default: this one)')
    cmp_ = sub.add_parser('compare', help='compare two records')
    cmp_.add_argument('a', type=Path)
    cmp_.add_argument('b', type=Path)
    args = parser.parse_args(argv)
    if args.command == 'record':
        return record(args.out.resolve(), args.n, args.root.resolve())
    return compare_files(args.a, args.b)


if __name__ == '__main__':
    sys.exit(main())
