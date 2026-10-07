"""Analyses and invariants of the stress test (tests/stress.py).

Each analysis reduces its result to a *signature*: name -> Entry, where an
Entry says how the value behaves when Z -> k*Z and when the points are
reordered. The invariants then compare signatures, so one comparison serves
DRT, Lin-KK, Z-HIT, R_inf and the fit alike. See doc/STRESS_TEST_PLAN.md.
"""

from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from eis_analysis.drt.core import calculate_drt
from eis_analysis.fitting import fit_equivalent_circuit
from eis_analysis.fitting.bounds import generate_simple_bounds
from eis_analysis.rinf_estimation.estimate import estimate_rinf
from eis_analysis.validation.kramers_kronig import kramers_kronig_validation
from eis_analysis.validation.zhit import zhit_validation

from tests.stress_cases import ORDER, Case, fit_start, k_powers

# Exact invariants (B, C) allow this relative difference: an iterative solver
# on a rescaled or reordered problem rounds differently, but nothing a user
# could see survives at 1e-6.
RTOL = 1e-6

# Absolute floors where a relative comparison is meaningless because the
# value is rounding noise (noise-free spectra fit to ~1e-11):
RES_ATOL = 1e-9         # relative residuals: 1e-9 of |Z|
CHI2_ATOL = 1e-12       # pseudo chi^2: residuals moving by RES_ATOL shift it less
PERCENT_ATOL = 1e-6     # fit error and noise estimate, in percent
# A fitted parameter may move by this fraction of its own stderr: LM stops
# on a relative cost change (ftol 1e-8), which leaves a poorly determined
# parameter that much play, and a printed CI cannot show it. Its stderr
# gets ten times more, since it moves with every parameter. On noise-free
# data the stderr itself is rounding noise, so the floor is the change that
# moves the fit by RES_ATOL instead (stderr scaled from the fit's residual
# level to RES_ATOL).
PARAM_STDERR_SHARE = 1e-3
STDERR_STDERR_SHARE = 1e-2


@dataclass
class Entry:
    """One value of a signature.

    power : how it scales under Z -> k*Z (value * k**power); per parameter
        for the fit, so it may be an array
    kind : 'exact' compared with ==; 'num' allclose with atol relative to
        max|value| (arrays with legitimate zeros, like gamma); 'rel' elementwise
        rtol only (parameters spanning decades); 'pf' like 'num' but one value
        per frequency, so it follows a reordering of the points
    finite : must be finite when not None (invariant A)
    atol : absolute tolerance in the units of `value` (scales with it);
        None = 1e-9 max|value| for 'num' and 'pf', 0 for 'rel'
    """
    value: Any
    power: Any = 0
    kind: str = 'num'
    finite: bool = True
    atol: Any = None


Signature = Dict[str, Entry]


def _inductance_atol(f, Z) -> float:
    """An inductance whose reactance at f_max is below 1e-9 |Z| is zero."""
    return 1e-9 * float(np.max(np.abs(Z))) / (2 * np.pi * float(np.max(f)))


# --- analyses ---------------------------------------------------------------

def _drt(case: Case, f, Z, k: float) -> Signature:
    r = calculate_drt(f, Z, auto_lambda=True)   # CLI defaults otherwise
    peaks = sorted(r.peaks or [], key=lambda p: p['tau'])
    return {
        'success': Entry(r.success, kind='exact'),
        'lambda': Entry(r.lambda_used),
        'tau': Entry(r.tau),
        'gamma': Entry(r.gamma, 1),
        'R_inf': Entry(r.R_inf, 1),
        'R_pol': Entry(r.R_pol, 1),
        'L_series': Entry(r.L_series, 1, atol=_inductance_atol(f, Z)),
        'rec_error': Entry(r.reconstruction_error),
        'Z_rec': Entry(r.Z_reconstructed, 1, 'pf'),
        'peak_tau': Entry(np.array([p['tau'] for p in peaks]), 0, 'rel'),
        'peak_R': Entry(np.array([p['R_estimate'] for p in peaks]), 1, 'rel'),
        'peak_flags': Entry([(p.get('boundary_sensitive'), p.get('outside_window'))
                             for p in peaks], kind='exact'),
    }


def _linkk(case: Case, f, Z, k: float) -> Signature:
    r = kramers_kronig_validation(f, Z, include_C=case.family == 'blocking')
    return {
        'error': Entry(r.error, kind='exact', finite=False),
        'M': Entry(r.M, kind='exact'),
        'M_lower': Entry(r.M_lower, kind='exact'),
        'mu': Entry(r.mu),
        'chi2': Entry(r.pseudo_chisqr, atol=CHI2_ATOL),
        'noise': Entry(r.noise_estimate, atol=PERCENT_ATOL),
        'extend': Entry(r.extend_decades),
        'res_re': Entry(r.residuals_real, 0, 'pf', atol=RES_ATOL),
        'res_im': Entry(r.residuals_imag, 0, 'pf', atol=RES_ATOL),
        'Z_fit': Entry(r.Z_fit, 1, 'pf'),
        'L': Entry(r.inductance, 1, atol=_inductance_atol(f, Z)),
        'C': Entry(r.capacitance, -1),
        'elements': Entry(r.elements, 1),
    }


def _zhit(case: Case, f, Z, k: float) -> Signature:
    r = zhit_validation(f, Z)
    return {
        'quality': Entry(r.quality),
        'chi2': Entry(r.pseudo_chisqr, atol=CHI2_ATOL),
        'noise': Entry(r.noise_estimate, atol=PERCENT_ATOL),
        'res_re': Entry(r.residuals_real, 0, 'pf', atol=RES_ATOL),
        'res_im': Entry(r.residuals_imag, 0, 'pf', atol=RES_ATOL),
        'Z_fit': Entry(r.Z_fit, 1, 'pf'),
    }


def _rinf(case: Case, f, Z, k: float) -> Signature:
    r = estimate_rinf(f, Z)
    # A fitted R_inf gets the play of a fitted parameter (PARAM_STDERR_SHARE)
    stderr = r.R_inf_stderr if r.method == 'rlq_fit' else None
    return {
        'R_inf': Entry(r.R_inf, 1, atol=PARAM_STDERR_SHARE * stderr
                       if stderr is not None and np.isfinite(stderr) else None),
        'R_inf_hf': Entry(r.R_inf_hf, 1),
        'f_hf': Entry(r.f_hf),
        'method': Entry(r.method, kind='exact'),
        'n_window': Entry(len(r.f_window), kind='exact'),
    }


def _fit(case: Case, f, Z, k: float, absolute_bounds: bool = False) -> Signature:
    """One LM fit from truth x U(0.3, 3); the start scales with Z.

    The bounds scale with Z as well (PARAMETER_BOUNDS x k^power), so the
    optimizer alone is under test: absolute bounds steer scipy's trust
    region by their distance even far from them, and are measured apart
    (`absolute_bounds`, invariant Babs).
    """
    powers = k_powers(case.labels)
    lower, upper = (np.array(b) for b in generate_simple_bounds(case.labels))
    scale = 1.0 if absolute_bounds else k ** powers
    r, _ = fit_equivalent_circuit(f, Z, case.circuit(fit_start(case) * k ** powers),
                                  bounds=(lower * scale, upper * scale))
    stderr = np.where(np.isfinite(r.params_stderr), r.params_stderr, 0.0)
    resolution = stderr * RES_ATOL / max(r.fit_error_rel / 100, 1e-300)
    return {
        'params': Entry(r.params_opt, powers, 'rel',
                        atol=np.maximum(PARAM_STDERR_SHARE * stderr, resolution)),
        'stderr': Entry(r.params_stderr, powers, 'rel', finite=False,
                        atol=RTOL * np.abs(r.params_opt)
                        + np.maximum(STDERR_STDERR_SHARE * stderr, resolution)),
        'fit_error': Entry(r.fit_error_rel, atol=PERCENT_ATOL),
        'well_conditioned': Entry(r.is_well_conditioned, kind='exact'),
        'bound_status': Entry(list(r.bound_status or []), kind='exact'),
    }


ANALYSES: Dict[str, Callable[..., Signature]] = {
    'drt': _drt, 'linkk': _linkk, 'zhit': _zhit, 'rinf': _rinf, 'fit': _fit,
}


def run_analysis(name: str, case: Case, f, Z, k: float = 1.0, **kwargs):
    """Signature, or the exception it raised."""
    try:
        return ANALYSES[name](case, f, Z, k, **kwargs)
    except Exception as e:  # invariant A reports it
        return e


# --- comparison -------------------------------------------------------------

def _as_array(value) -> Optional[np.ndarray]:
    if value is None:
        return None
    return np.atleast_1d(np.asarray(value))


def _compare(name: str, ref: Entry, got: Entry, k: float,
             perm: Optional[np.ndarray], bitwise: bool) -> Optional[str]:
    """Mismatch description, or None. `got` was computed on k*Z[perm]."""
    if ref.kind == 'exact':
        return None if ref.value == got.value else f'{name}: {ref.value!r} -> {got.value!r}'
    a, b = _as_array(ref.value), _as_array(got.value)
    if a is None or b is None:
        return None if a is None and b is None else f'{name}: None mismatch'
    if a.shape != b.shape:
        return f'{name}: shape {a.shape} -> {b.shape}'
    if ref.kind == 'pf' and perm is not None:
        a = a[perm]
    k_scale = np.asarray(k, dtype=float) ** np.asarray(ref.power)
    expected = a * k_scale
    if bitwise:
        same = np.array_equal(a, b, equal_nan=True)
        return None if same else f'{name}: not bit-identical'
    finite = np.isfinite(expected)
    if not np.array_equal(finite, np.isfinite(b)) or not np.array_equal(
            expected[~finite], b[~finite], equal_nan=True):
        return f'{name}: non-finite pattern differs'
    if ref.atol is not None:
        atol = np.broadcast_to(np.asarray(ref.atol) * k_scale, finite.shape)[finite]
    expected, b = expected[finite], b[finite]
    if expected.size == 0:
        return None
    if ref.atol is None:
        atol = 0.0 if ref.kind == 'rel' else 1e-9 * np.max(np.abs(expected))
    excess = np.abs(b - expected) - (atol + RTOL * np.abs(expected))
    if np.all(excess <= 0):
        return None
    worst = int(np.argmax(excess))
    return (f'{name}[{worst}]: expected {expected.flat[worst]:.6g}, '
            f'got {b.flat[worst]:.6g} (rel {abs(b.flat[worst] / expected.flat[worst] - 1):.1e})')


def compare(ref: Signature, got: Signature, k: float = 1.0,
            perm: Optional[np.ndarray] = None, bitwise: bool = False) -> List[str]:
    return [m for name, entry in ref.items()
            if (m := _compare(name, entry, got[name], k, perm, bitwise))]


# --- invariants -------------------------------------------------------------
# Each returns (status, detail): status 'pass' | 'fail' | 'meze' | 'skip'.

Check = Tuple[str, str]


def invariant_A(sig) -> Check:
    """Robustness: no exception, no NaN/Inf, a non-empty result."""
    if isinstance(sig, Exception):
        return 'fail', f'{type(sig).__name__}: {sig}'
    if 'success' in sig and not sig['success'].value:
        return 'fail', 'empty result'
    if 'error' in sig and sig['error'].value is not None:
        return 'fail', f'error: {sig["error"].value}'
    for name, entry in sig.items():
        if entry.kind == 'exact' or not entry.finite or entry.value is None:
            continue
        a = _as_array(entry.value)
        if not np.all(np.isfinite(a)):
            return 'fail', f'{name} not finite'
    return 'pass', ''


def _mismatch(errors: List[str]) -> Check:
    return ('pass', '') if not errors else ('fail', '; '.join(errors[:3]))


def _scaled_check(name: str, case: Case, ref: Signature, k: float, **kwargs) -> Check:
    got = run_analysis(name, case, case.frequencies, case.Z * k, k, **kwargs)
    if isinstance(got, Exception):
        return 'fail', f'k={k:g}: {type(got).__name__}: {got}'
    status, detail = _mismatch(compare(ref, got, k))
    return status, f'k={k:g}: {detail}' if detail else ''


def invariant_B(name: str, case: Case, ref: Signature) -> List[Tuple[str, Check]]:
    """Units: Z -> k*Z scales every result by its power of k.

    For the fit also Babs: the same with the absolute PARAMETER_BOUNDS the
    library uses. A difference there is the known limit of absolute bounds,
    so it is counted as 'meze', never as a failure.
    """
    checks = []
    for k in (1e-3, 1e3):
        sign = '-' if k < 1 else '+'
        checks.append((f'B{sign}', _scaled_check(name, case, ref, k)))
        if name == 'fit':
            ref_abs = run_analysis(name, case, case.frequencies, case.Z, absolute_bounds=True)
            status, detail = (_scaled_check(name, case, ref_abs, k, absolute_bounds=True)
                              if not isinstance(ref_abs, Exception) else ('fail', str(ref_abs)))
            checks.append((f'Babs{sign}', ('meze' if status == 'fail' else status, detail)))
    return checks


def invariant_C(name: str, case: Case, ref: Signature) -> List[Tuple[str, Check]]:
    """Order: reversed and shuffled points give the same result."""
    n = len(case.frequencies)
    checks = []
    for label, perm in (('Crev', np.arange(n)[::-1]),
                        ('Cmix', case.rng(ORDER).permutation(n))):
        got = run_analysis(name, case, case.frequencies[perm], case.Z[perm])
        if isinstance(got, Exception):
            checks.append((label, ('fail', f'{type(got).__name__}: {got}')))
        else:
            checks.append((label, _mismatch(compare(ref, got, perm=perm))))
    return checks


def invariant_K(name: str, case: Case, ref: Signature) -> List[Tuple[str, Check]]:
    """Determinism: a second run with the same input and seeds is bit-identical."""
    got = run_analysis(name, case, case.frequencies, case.Z)
    if isinstance(got, Exception):
        return [('K', ('fail', f'{type(got).__name__}: {got}'))]
    return [('K', _mismatch(compare(ref, got, bitwise=True)))]


INVARIANTS = ('A', 'B-', 'B+', 'Babs-', 'Babs+', 'Crev', 'Cmix', 'K')


def check_case(case: Case, analyses=tuple(ANALYSES)) -> List[Tuple[str, str, str, str]]:
    """All invariants of one case: [(analysis, invariant, status, detail)]."""
    rows = []
    for name in analyses:
        ref = run_analysis(name, case, case.frequencies, case.Z)
        status, detail = invariant_A(ref)
        rows.append((name, 'A', status, detail))
        if status != 'pass':
            rows += [(name, inv, 'skip', 'A failed') for inv in INVARIANTS[1:]
                     if name == 'fit' or not inv.startswith('Babs')]
            continue
        for checker in (invariant_B, invariant_C, invariant_K):
            rows += [(name, inv, s, d) for inv, (s, d) in checker(name, case, ref)]
    return rows
