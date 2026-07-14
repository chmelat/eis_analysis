"""
Lambda-probe peak stability diagnostics for DRT analysis.

Re-solves the regularized NNLS problem at several lambda values around the
selected lambda* and tracks which peaks of the main solution persist. A peak
that appears only in a narrow lambda window is likely a regularization
artifact, while a peak stable across a decade of lambda reflects real
relaxation structure in the data.

Peak detection inside the probe always uses scipy.signal.find_peaks (fast and
deterministic), even when the main analysis uses GMM; GMM reference peaks are
matched against the scipy peaks of each probe solution by proximity in
log10(tau).
"""

import numpy as np
import logging
from typing import List, Tuple
from numpy.typing import NDArray
from scipy.signal import find_peaks

from .results import DRTMatrices, LambdaProbePoint, PeakStability, StabilityDiagnostics
from .estimation import _estimate_peak_resistance
from .linear_system import _solve_nnls
from ..fitting.config import DRT_PEAK_HEIGHT_THRESHOLD

logger = logging.getLogger(__name__)

# Probe lambdas relative to lambda*: half a decade and a full decade on each
# side. One decade spans the range where a physically meaningful peak should
# survive; the half-decade points detect drift before it disappears.
PROBE_EXPONENTS = (-1.0, -0.5, 0.5, 1.0)

# Absolute lambda bounds for the probe. Matches the extremes used elsewhere in
# the package: 1e-6 is below the GCV search floor (1e-5), 1.0 is its ceiling.
PROBE_LAMBDA_MIN = 1e-6
PROBE_LAMBDA_MAX = 1.0

# Peak matching tolerance in log10(tau) decades: half the distance to the
# nearest neighbouring reference peak, clamped to [FLOOR, CAP]. The floor
# absorbs ordinary grid-level drift (a few tau bins at n_tau=100 over a
# typical 7-decade range); the cap keeps an isolated peak from "matching"
# unrelated structure far away.
MATCH_TOLERANCE_FLOOR = 0.3
MATCH_TOLERANCE_CAP = 0.5

# Verdict thresholds. A stable peak must survive every successful probe with
# a position drift well inside the matching tolerance and an area (R) change
# small compared to typical DRT peak-area uncertainty (~25%). A peak missing
# from more than half of the probes is classified as an artifact.
STABLE_MAX_DRIFT_DECADES = 0.2
STABLE_MAX_R_VARIATION = 0.25


def _match_tolerances(ref_log_taus: NDArray) -> NDArray:
    """
    Per-peak matching tolerance in decades of tau.

    Half the log10 distance to the nearest neighbouring reference peak,
    clamped to [MATCH_TOLERANCE_FLOOR, MATCH_TOLERANCE_CAP]. A single peak
    gets the cap.
    """
    n = len(ref_log_taus)
    if n == 1:
        return np.array([MATCH_TOLERANCE_CAP])

    tolerances = np.empty(n)
    for i in range(n):
        gaps = np.abs(ref_log_taus - ref_log_taus[i])
        nearest = np.min(gaps[gaps > 0]) if np.any(gaps > 0) else np.inf
        tolerances[i] = min(max(nearest / 2, MATCH_TOLERANCE_FLOOR),
                            MATCH_TOLERANCE_CAP)
    return tolerances


def _match_peaks(ref_log_taus: NDArray, tolerances: NDArray,
                 probe_peaks: List[dict]) -> List[Tuple[int, dict]]:
    """
    Greedily match probe peaks to reference peaks by log10(tau) distance.

    Each probe peak is assigned to at most one reference peak and vice versa;
    pairs are taken in order of increasing distance, skipping pairs beyond the
    reference peak's tolerance.

    Returns list of (reference_index, probe_peak) pairs.
    """
    if len(probe_peaks) == 0 or len(ref_log_taus) == 0:
        return []

    probe_log_taus = np.log10([p['tau'] for p in probe_peaks])
    # All candidate pairs sorted by distance
    pairs = sorted(
        ((abs(probe_log_taus[j] - ref_log_taus[i]), i, j)
         for i in range(len(ref_log_taus))
         for j in range(len(probe_peaks))),
        key=lambda item: item[0]
    )

    matched: List[Tuple[int, dict]] = []
    used_refs = set()
    used_probes = set()
    for distance, i, j in pairs:
        if i in used_refs or j in used_probes:
            continue
        if distance >= tolerances[i]:
            continue
        matched.append((i, probe_peaks[j]))
        used_refs.add(i)
        used_probes.add(j)
    return matched


def _detect_probe_peaks(tau: NDArray, gamma: NDArray) -> List[dict]:
    """Scipy peak detection with the package-wide height threshold."""
    gamma_max = float(np.max(gamma)) if len(gamma) else 0.0
    if gamma_max <= 0:
        return []
    peaks_idx, _ = find_peaks(gamma, height=gamma_max * DRT_PEAK_HEIGHT_THRESHOLD)
    resistances = _estimate_peak_resistance(tau, gamma, peaks_idx)
    return [
        {'tau': float(tau[idx]), 'R_estimate': float(resistances[i])}
        for i, idx in enumerate(peaks_idx)
    ]


def _assess_peaks(reference_peaks: List[Tuple[float, float]],
                  probe_points: List[LambdaProbePoint]) -> List[PeakStability]:
    """
    Aggregate per-peak stability over the successful probe points.
    """
    successful = [p for p in probe_points if p.success]
    n_probes = len(successful)
    if n_probes == 0 or len(reference_peaks) == 0:
        return []

    ref_log_taus = np.log10([tau for tau, _ in reference_peaks])
    tolerances = _match_tolerances(ref_log_taus)

    persistence = np.zeros(len(reference_peaks), dtype=int)
    max_drift = np.zeros(len(reference_peaks))
    max_r_variation = np.zeros(len(reference_peaks))

    for point in successful:
        for i, probe_peak in _match_peaks(ref_log_taus, tolerances, point.peaks):
            persistence[i] += 1
            drift = abs(np.log10(probe_peak['tau']) - ref_log_taus[i])
            max_drift[i] = max(max_drift[i], drift)
            R_ref = reference_peaks[i][1]
            if R_ref > 0:
                variation = abs(probe_peak['R_estimate'] - R_ref) / R_ref
                max_r_variation[i] = max(max_r_variation[i], variation)

    results = []
    for i, (tau_ref, R_ref) in enumerate(reference_peaks):
        if persistence[i] <= n_probes // 2:
            verdict = 'artifact'
        elif (persistence[i] == n_probes
              and max_drift[i] < STABLE_MAX_DRIFT_DECADES
              and max_r_variation[i] < STABLE_MAX_R_VARIATION):
            verdict = 'stable'
        else:
            verdict = 'marginal'
        results.append(PeakStability(
            tau_ref=float(tau_ref),
            f_ref=float(1 / (2 * np.pi * tau_ref)),
            R_ref=float(R_ref),
            persistence=int(persistence[i]),
            n_probes=n_probes,
            tau_drift_decades=float(max_drift[i]),
            R_variation_rel=float(max_r_variation[i]),
            verdict=verdict
        ))
    return results


def probe_lambda_stability(matrices: DRTMatrices, lambda_star: float,
                           reference_peaks: List[Tuple[float, float]],
                           Z: NDArray, R_inf: float,
                           n_tau: int) -> StabilityDiagnostics:
    """
    Assess peak stability by re-solving the DRT at lambdas around lambda*.

    Parameters
    ----------
    matrices : DRTMatrices
        System matrices from the main DRT run (reused, not rebuilt).
    lambda_star : float
        Regularization parameter of the main solution.
    reference_peaks : list of (tau, R_estimate)
        Peaks of the main solution to track across the probe.
    Z : ndarray
        Complex impedance [Ohm] (for reconstruction error).
    R_inf : float
        High-frequency resistance used in the main run [Ohm].
    n_tau : int
        Number of tau grid points.

    Returns
    -------
    StabilityDiagnostics
        Probe solutions, per-peak stability, and summary warnings.
    """
    probe_lambdas = np.clip(lambda_star * 10.0 ** np.array(PROBE_EXPONENTS),
                            PROBE_LAMBDA_MIN, PROBE_LAMBDA_MAX)
    # Deduplicate after clipping and drop values ~equal to lambda* (the main
    # solution already covers lambda*).
    unique_lambdas: List[float] = []
    for lam in probe_lambdas:
        if np.isclose(lam, lambda_star, rtol=1e-3):
            continue
        if any(np.isclose(lam, existing, rtol=1e-3) for existing in unique_lambdas):
            continue
        unique_lambdas.append(float(lam))

    probe_points: List[LambdaProbePoint] = []
    for lam in unique_lambdas:
        solution = _solve_nnls(matrices.A, matrices.b, matrices.L, lam, n_tau, Z)
        if not solution.success or solution.gamma is None:
            message = '; '.join(solution.warnings) or 'solver failed'
            probe_points.append(LambdaProbePoint(
                lambda_value=lam, success=False, error=message))
            logger.debug(f"Lambda probe at {lam:.2e} failed: {message}")
            continue

        gamma = solution.gamma
        Z_reconstructed = R_inf + (matrices.A_re + 1j * matrices.A_im) @ gamma
        rel_error = float(np.mean(np.abs(Z - Z_reconstructed) / np.abs(Z)) * 100)
        probe_points.append(LambdaProbePoint(
            lambda_value=lam,
            success=True,
            gamma=gamma,
            gamma_max=float(np.max(gamma)),
            reconstruction_error_rel=rel_error,
            peaks=_detect_probe_peaks(matrices.tau, gamma)
        ))

    peak_stability = _assess_peaks(reference_peaks, probe_points)

    warnings = []
    n_failed = sum(1 for p in probe_points if not p.success)
    if n_failed:
        warnings.append(f"{n_failed} of {len(probe_points)} lambda probe "
                        f"solutions failed")
    n_unstable = sum(1 for p in peak_stability if p.verdict != 'stable')
    if n_unstable:
        warnings.append(
            f"{n_unstable} of {len(peak_stability)} peaks not stable across "
            f"lambda probe - interpret with caution"
        )

    return StabilityDiagnostics(
        lambda_star=lambda_star,
        probe_points=probe_points,
        peak_stability=peak_stability,
        warnings=warnings
    )
