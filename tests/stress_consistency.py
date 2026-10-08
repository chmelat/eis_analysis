"""Consistency invariants of the stress test (tests/stress.py): D, E, F, I, M.

Unlike A, B, C, K (tests/stress_invariants.py), which compare an analysis
with itself, these compare it with the truth the case was generated from.
Each check returns (analysis, invariant, status, detail, value). For a
threshold check `value` is measured / allowed, so it passes at <= 1 and the
runner's quantiles show how much margin a threshold has (plan:
doc/STRESS_TEST_PLAN.md, calibrated with ~2x margin over the first run).
The allowance is always factor x sigma + floor x |Z_true|: the noise the
case was given, plus the method's own error on exact data. Absolute, not
relative to the measured |Z|, which at constant noise can be pure noise
near zero. Status 'stat' marks a per-case contribution to an aggregate rate
(F3, M), and 'lokmin' a fit that ended in a local minimum (F2): both are
rates, never failures (plan: F2 is a rate by identifiability class).

DRT, Lin-KK, Z-HIT and R_inf are not run again: check_case hands over the
results it computed (`raw`). A check whose analysis failed A is skipped;
A reports it. gamma's finiteness is A's too; plan item "every peak inside
the window or flagged" is not checked: the flags are computed from that
very distance (drt/estimation.py, _flag_boundary_peaks), so it cannot fail.
"""

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from eis_analysis.analysis.local_exponent import local_exponent
from eis_analysis.fitting import fit_circuit_multistart, fit_equivalent_circuit
from eis_analysis.fitting.bounds import generate_simple_bounds
from eis_analysis.fitting.diagnostics import compute_weights

from tests.stress_cases import MULTISTART, Case, fit_start
from tests.stress_invariants import RES_ATOL, RTOL, drt_peaks

# --- thresholds (PROVISIONAL until calibrated from the first run) ------------
# Pointwise residual checks: |dZ_i| <= NOISE_FACTOR * sigma_i + FLOOR * |Z_i|.
# NOISE_FACTOR covers the largest of ~2N Gaussian draws (~3.5 sigma at
# N = 150); FLOOR is the method's own error on exact data.
LINKK_NOISE_FACTOR, LINKK_FLOOR = 5.0, 1e-3
ZHIT_NOISE_FACTOR, ZHIT_FLOOR = 5.0, 1e-2
# DRT reconstruction, as an RMS over the points of |dZ_i| / allowance_i
DRT_REC_NOISE_FACTOR, DRT_REC_FLOOR = 1.5, 1e-3
# Single-value checks read one noisy point: 3 sigma of it on top of the tolerance
POINT_SIGMAS = 3.0
# A spectrum counts as closed at an end when its true phase there is above
# this (plan: -5 deg), so the remaining arc is a small part of Re Z.
CLOSED_PHASE_DEG = -5.0
DRT_DC_TOL = 0.05            # |R_inf + R_pol - Re Z(f_min)|, relative to Re Z(f_min)
DRT_PEAK_TOL_DEC = 0.15      # true tau to the nearest DRT peak
RINF_HF_TOL = 0.01           # R_inf above Re Z_true(f_max), relative
RINF_CLOSED_TOL = 0.05       # |R_inf - Rs| / Rs on a closed high-frequency end
F2_RTOL = 1e-6               # F2: cost(fit) <= cost(truth) * (1 + F2_RTOL) + floor


Row = Tuple[str, str, str, str, Optional[float]]


def _row(analysis: str, inv: str, ok: bool, detail: str, value) -> Row:
    return analysis, inv, 'pass' if ok else 'fail', detail if not ok else '', value


def _ratio_row(analysis: str, inv: str, measured, allowed, detail: str) -> Row:
    """Pass when measured <= allowed everywhere; value: the largest ratio."""
    ratio = np.atleast_1d(np.asarray(measured) / np.maximum(allowed, 1e-300))
    worst = int(np.argmax(ratio))
    return _row(analysis, inv, bool(ratio[worst] <= 1),
                f'{detail} (point {worst}, {ratio[worst]:.3g}x allowed)', float(ratio[worst]))


def _allowance(case: Case, factor: float, floor: float, idx=slice(None)):
    return factor * case.sigma[idx] + floor * np.abs(case.Z_clean[idx])


def _phase_deg(Z: complex) -> float:
    return float(np.degrees(np.angle(Z)))


# --- D: Kramers-Kronig -------------------------------------------------------

def check_D(case: Case, raw: Dict[str, Any], rows: List[Row]) -> None:
    """Every generated spectrum is KK-consistent, so Lin-KK (with the series
    C where the truth's end is capacitive) and Z-HIT fit it to its noise plus their own
    error on exact data."""
    Z_mag = np.abs(case.Z)      # both normalize their residuals by the measured |Z|
    if 'linkk' in raw:
        kk = raw['linkk']
        kk_abs = np.maximum(np.abs(kk.residuals_real), np.abs(kk.residuals_imag)) * Z_mag
        rows.append(_ratio_row('linkk', 'D', kk_abs,
                               _allowance(case, LINKK_NOISE_FACTOR, LINKK_FLOOR), f'M {kk.M}'))
    if 'zhit' in raw:
        zh = raw['zhit']
        rows.append(_ratio_row('zhit', 'D', np.abs(zh.residuals_mag) / 100 * Z_mag,
                               _allowance(case, ZHIT_NOISE_FACTOR, ZHIT_FLOOR),
                               f'quality {zh.quality:.3g}'))


# --- E: DRT ------------------------------------------------------------------

def _true_rc_taus(case: Case) -> np.ndarray:
    """tau = R C of each (R|C) arc of an `rc` case: every C follows its R."""
    return np.array([case.truth[i - 1] * case.truth[i]
                     for i, label in enumerate(case.labels) if label == 'C'])


def check_E(case: Case, raw: Dict[str, Any], rows: List[Row]) -> None:
    if 'drt' not in raw:
        return
    r = raw['drt']
    gamma_min = float(np.min(r.gamma) / max(np.max(r.gamma), 1e-300))
    rows.append(_row('drt', 'Eneg', gamma_min >= 0, f'min gamma / max = {gamma_min:.2e}', None))

    i_max = int(np.argmax(case.frequencies))
    # Only where the high-frequency end is closed: on an open one the DRT's
    # median R_inf overestimates and gamma >= 0 cannot make up for it - the
    # known limit behind the deferred --ri-fit default (rc/149: 3469 Ohm
    # against R_s = 1383), not what this check is for.
    if case.family == 'rc' and abs(_phase_deg(case.Z_clean[i_max])) < -CLOSED_PHASE_DEG:
        # sqrt(2) sigma: |dZ| has noise on Re and Im
        e = np.abs(case.Z - r.Z_reconstructed) / _allowance(
            case, DRT_REC_NOISE_FACTOR * np.sqrt(2), DRT_REC_FLOOR)
        rms = float(np.sqrt(np.mean(e ** 2)))
        rows.append(_row('drt', 'Erec', rms <= 1,
                         f'error {r.reconstruction_error:.3g} %, RMS {rms:.3g}x allowed', rms))

    i_min = int(np.argmin(case.frequencies))
    if _phase_deg(case.Z_clean[i_min]) > CLOSED_PHASE_DEG:
        # R_pol is extrapolated to DC, so anything between Re Z(f_min) and
        # the true DC resistance is right
        lo = float(case.Z_clean[i_min].real)
        hi = float(case.circuit().impedance(np.array([1e-6 * case.frequencies[i_min]]),
                                            list(case.truth))[0].real)
        total = r.R_inf + r.R_pol
        allowed = DRT_DC_TOL * lo + POINT_SIGMAS * case.sigma[i_min]
        rows.append(_ratio_row('drt', 'Edc', max(lo - total, total - hi, 0.0), allowed,
                               f'R_inf + R_pol = {total:.4g}, Re Z from {lo:.4g} (f_min) to {hi:.4g} (DC)'))

    if (case.family == 'rc' and case.noise_level == 0 and case.min_frac >= 0.05
            and (case.min_sep is None or case.min_sep >= 1.0)):
        peaks = np.log10([p.get('tau', p.get('tau_center')) for p in drt_peaks(r)])
        dist = [float(np.min(np.abs(peaks - lt))) if peaks.size else np.inf
                for lt in np.log10(_true_rc_taus(case))]
        rows.append(_ratio_row('drt', 'Epeak', dist, DRT_PEAK_TOL_DEC,
                               f'a true tau is {max(dist):.2f} dec from the nearest of {peaks.size} peaks'))


# --- F: circuit fit ----------------------------------------------------------

def _cost(case: Case, Z_model: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(np.abs((case.Z - Z_model) * weights) ** 2))


def check_F(case: Case, raw: Dict[str, Any], rows: List[Row]) -> None:
    # F1: noise-free, started at the truth -> stays there
    if case.noise_level == 0:
        r1, _ = fit_equivalent_circuit(case.frequencies, case.Z, case.circuit())
        rows.append(_ratio_row('fit', 'F1', np.abs(r1.params_opt / case.truth - 1), RTOL,
                               'parameter moved from the truth'))

    # F2: multistart from truth x U(0.3, 3) is no worse than the truth
    ms, Z_fit = fit_circuit_multistart(case.circuit(fit_start(case)), case.frequencies,
                                       case.Z, rng=case.rng(MULTISTART))
    best = ms.best_result
    w = compute_weights(case.Z, 'modulus')     # the multistart default
    cost_fit, cost_true = _cost(case, Z_fit, w), _cost(case, case.Z_clean, w)
    floor = RES_ATOL ** 2 * float(np.sum(np.abs(case.Z * w) ** 2))
    row = _ratio_row('fit', 'F2', cost_fit, cost_true * (1 + F2_RTOL) + floor,
                     f'cost {cost_fit / max(cost_true, floor):.4g}x truth, '
                     f'error {best.fit_error_rel:.3g} %')
    f2_ok = row[2] == 'pass'
    rows.append(row if f2_ok else ('fit', 'F2', 'lokmin', row[3], row[4]))

    # F3: coverage of the CI the library shows, only where it is meant to hold
    if (f2_ok and case.noise_level > 0 and case.noise_kind == 'proportional'
            and best.is_well_conditioned):
        lo, hi = best.params_ci_95
        covered = int(np.sum((lo <= case.truth) & (case.truth <= hi)))
        rows.append(('fit', 'F3', 'stat', f'{covered}/{len(case.truth)}',
                     covered / len(case.truth)))

    # F4: a parameter pinned at a bound (not a lower bound of 0, which the
    # library deliberately never reports) is reported. The multistart fits
    # with the default bounds.
    lower, upper = (np.asarray(b) for b in generate_simple_bounds(case.labels))
    p = best.params_opt
    at = (((lower > 0) & (np.abs(p - lower) <= RTOL * np.abs(lower)))
          | (np.abs(p - upper) <= RTOL * np.abs(upper)))
    if at.any():
        warned = bool(best.diagnostics and best.diagnostics.bounds_warnings)
        rows.append(_row('fit', 'F4', warned,
                         f'parameters {np.flatnonzero(at).tolist()} at a bound, no warning', None))


# --- I: R_inf ----------------------------------------------------------------

def check_I(case: Case, raw: Dict[str, Any], rows: List[Row]) -> None:
    """R_inf_range holds the true Rs (the n(f) map relies on it), R_inf is
    at most Re Z_true(f_max) (passivity), and Rs where the high-frequency
    end is closed. Cases with L are excluded: L lifts Re Z
    nowhere, but its reactance hides how close f_max is to closing."""
    if case.has_L or 'rinf' not in raw:
        return
    r = raw['rinf']
    i = int(np.argmax(case.frequencies))
    re_hf = float(case.Z_clean[i].real)
    noise = POINT_SIGMAS * case.sigma[i]
    # Irange: R_inf_range, which the n(f) map relies on, holds the true Rs
    Rs = case.truth[case.labels.index('R')]
    lo, hi = r.R_inf_range
    rows += [_ratio_row('rinf', 'Irange', max(lo - Rs, Rs - hi, 0.0), RINF_HF_TOL * Rs + noise,
                        f'Rs {Rs:.4g} outside R_inf_range {lo:.4g}..{hi:.4g} ({r.method})'),
             _ratio_row('rinf', 'Ihf', max(r.R_inf - re_hf, 0.0), RINF_HF_TOL * re_hf + noise,
                        f'R_inf {r.R_inf:.4g} above Re Z_true(f_max) {re_hf:.4g} ({r.method})')]
    if _phase_deg(case.Z_clean[i]) > CLOSED_PHASE_DEG:
        rows.append(_ratio_row('rinf', 'Iclosed', abs(r.R_inf - Rs), RINF_CLOSED_TOL * Rs + noise,
                               f'R_inf {r.R_inf:.4g} vs Rs {Rs:.4g} ({r.method})'))
        # Ifit: how often the window fit determines R_s on a closed end,
        # where an "R_inf not determined" warning is a false alarm - the
        # remaining precondition for making --ri-fit the default
        fitted = r.method == 'rlq_fit'
        rows.append(('rinf', 'Ifit', 'stat', f'{int(fitted)}/1 {r.method}', float(fitted)))


# --- M: n(f) map through the CLI path -----------------------------------------

def check_M(case: Case, raw: Dict[str, Any], rows: List[Row]) -> None:
    """Points the noisy map calls determined agree with the noise-free map
    within 2 x their uncertainty. Aggregate rate, plus the worst deviation."""
    if case.noise_level == 0 or 'rinf' not in raw:
        return
    est = raw['rinf']
    # As cli/handlers/local_exponent.py
    r = local_exponent(case.frequencies, case.Z, est.R_inf, est.L, R_inf_range=est.R_inf_range)
    Rs = case.truth[case.labels.index('R')]
    L_true = case.truth[case.labels.index('L')] if 'L' in case.labels else 0.0
    ref = local_exponent(case.frequencies, case.Z_clean, Rs, L_true)
    # Both determined: the reference near f_max is rounding noise otherwise
    use = r.valid & ref.valid
    if not use.any():
        return
    dev = np.abs(r.n - ref.n)[use]
    within = int(np.sum(dev <= 2 * r.n_uncertainty[use]))
    rows.append(('n(f)', 'M', 'stat', f'{within}/{use.sum()} max |dn| {np.max(dev):.3f}',
                 float(np.max(dev))))


# Each check appends to `rows`, so what it measured before an exception stays.
# On an exception the row goes to the check's main analysis and invariant.
CONSISTENCY: Tuple[Tuple[Callable[..., None], str, str], ...] = (
    (check_D, 'linkk', 'D'), (check_E, 'drt', 'Eneg'), (check_F, 'fit', 'F2'),
    (check_I, 'rinf', 'Ihf'), (check_M, 'n(f)', 'M'),
)
CONSISTENCY_INVARIANTS = ('D', 'Eneg', 'Erec', 'Edc', 'Epeak', 'F1', 'F2', 'F3', 'F4',
                          'Irange', 'Ihf', 'Iclosed', 'Ifit', 'M')
# Their value is measured / allowed (pass at <= 1), shown in the calibration table
RATIO_INVARIANTS = ('D', 'Erec', 'Edc', 'Epeak', 'F1', 'F2', 'Irange', 'Ihf', 'Iclosed')


def consistency_rows(case: Case, raw: Dict[str, Any]) -> List[Row]:
    """All consistency checks of one case; an exception is a failed check."""
    rows: List[Row] = []
    for check, analysis, inv in CONSISTENCY:
        try:
            check(case, raw, rows)
        except Exception as e:
            rows.append((analysis, inv, 'fail', f'{type(e).__name__}: {e}', None))
    return rows
