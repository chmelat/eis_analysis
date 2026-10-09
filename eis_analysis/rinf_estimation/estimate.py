"""
R_inf estimation: R-L-(R|Q) fit over the top frequency decades.

Model: Z(omega) = R_s + j*omega*L + R_k / (1 + R_k*Q*(j*omega)^n)

One nonlinear fit covers inductive, capacitive and mixed high-frequency ends
alike; R_inf = R_s. The standard error of R_s decides whether the window
determines R_inf at all. When it does not, or when the fit cannot run, the
result falls back to Re(Z) at the highest frequency with Im(Z) <= 0, an upper
bound (every passive term adds Re >= 0 to R_s, and at Im = 0 a series L
contributes nothing), and says why in `warnings`. An inductive top is skipped:
above Im = 0 lead artifacts can pull Re(Z) below R_s, even below zero.

Clean design: No logging in core functions, all diagnostics returned as data.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting import FitResult, L, Q, R, fit_equivalent_circuit
from ..fitting.bounds import PARAMETER_BOUNDS

# Fit window: f >= f_max / 10**RINF_FIT_DECADES. Two decades is the compromise
# measured in doc/archive/AUDIT_ri_fit_2026-09-25.md: one decade resolves an arc just
# above f_max better (-6.8 % vs -20.4 %), but scatters ~10x more under 1 %
# noise on a flat end (+-2.1 % vs +-0.2 %).
RINF_FIT_DECADES = 2

# The model has 5 free parameters; 5 complex points give 10 real residuals,
# i.e. 5 degrees of freedom, the least for which the covariance (and hence
# the stderr the identifiability test relies on) means anything.
RINF_FIT_MIN_POINTS = 5

# Relative stderr of R_s above which R_inf counts as not determined by the
# window. On the audit's synthetic set (1 % noise) determinable cases stay at
# 0.2-3.6 %, non-determinable ones (arc above f_max, strongly open CPE arc,
# model mismatch) start at 14 %: 5 % leaves a ~2x margin on both sides.
RINF_REL_STDERR_MAX = 0.05

# HF median: up to 5 highest points, but no more than 10 % of the spectrum,
# so short spectra do not reach down into the arc.
HF_MEDIAN_MAX_POINTS = 5

# Range of the window fit, relative to the window's impedance: every term from
# min|Z|/RINF_BOUND_RANGE to max|Z|*RINF_BOUND_RANGE over the window. Above
# 1e6 times the largest |Z| an arc resistance is an open arc, and Q the same
# for a short, to better than 1e-6 - under any instrument's resolution (~1e-4
# of |Z|) - so these are the upper bounds. R and L are bounded below by their
# physical limit 0 instead: a relative floor cut off what the data still
# determine (noise-free 0.05 Ohm R_s in front of a GOhm film, |Z| 3e3..2.5e5
# Ohm: an L floor at 1e-6 |Z| biased R_s to 0.064 Ohm, and R_inf fell back
# to 243 Ohm). The lower end of the range only places the starts, so that
# none is 0. Absolute PARAMETER_BOUNDS made R_inf
# depend on the units: the open arc's R_k ran into R <= 1e10 Ohm (or Q into
# 0.1, L into 1e-4 H), and Z -> 1000 Z shifted R_inf on 51 % of the stress
# test's random spectra, by over 10 % on 7 % of them.
RINF_BOUND_RANGE = 1e6

# R_inf_upper: the HF bound widened by this many standard deviations of the
# noise at f_hf, so that a point shifted by noise still bounds R_s. The
# noise is read off the window fit's residuals near f_hf (model mismatch
# included, which only widens the bound). 3 sigma: a 0.3 % chance per
# spectrum that the noise at f_hf exceeds it. Also the width of
# R_inf_range around a fitted R_s, in units of its stderr.
RINF_BOUND_NOISE_SIGMAS = 3.0

# The noise at f_hf is read from the residuals of the points down to this
# many decades below it: near f_hf, where |Z| is close to |Z(f_hf)|, so
# proportional and constant noise read alike (over the whole window, the
# absolute residuals of a capacitive end are dominated by its bottom, up to
# 100x |Z(f_hf)|). Half a decade holds >= 3 points at >= 5 points/decade.
RINF_NOISE_LOCAL_DECADES = 0.5

# Initial n of the CPE: midway in the typical 0.6-1.0 range of real arcs.
_N_GUESS = 0.8

# Initial L when the top point is capacitive (Im <= 0): reactance 0.1 % of
# |Z| at f_max, so the fit starts near "no inductance".
_L_GUESS_REACTANCE_SHARE = 1e-3


@dataclass
class RinfResult:
    """R_inf estimate with the fit and the fallback it was chosen from.

    `R_inf` is the value to use: the fitted R_s when the window determines
    it (`method == 'rlq_fit'`), otherwise the HF upper bound `R_inf_hf`
    (`method == 'hf_bound'`, clipped at 0) with the reason in `warnings`.
    """
    R_inf: float  # [Ohm]
    method: str  # 'rlq_fit' | 'hf_bound'
    R_inf_hf: float  # HF upper bound of R_inf, see _hf_bound [Ohm]
    # R_inf_hf + RINF_BOUND_NOISE_SIGMAS x the noise at f_hf, clipped at 0:
    # what R_s cannot exceed even if f_hf is noisy. With no fit to read the
    # noise from: R_inf_hf, or when it is <= 0 and bounds nothing, the larger
    # of max Re(Z) over the window and |Z(f_hf)| [Ohm]
    R_inf_upper: float
    f_hf: float  # frequency of R_inf_hf; below f_max if the top is inductive [Hz]
    f_window: NDArray[np.float64]  # frequencies of the fit window [Hz]
    Z_window: NDArray[np.complexfloating]  # impedance of the fit window [Ohm]
    fit: Optional[FitResult] = None  # None if the fit did not run
    warnings: List[str] = field(default_factory=list)

    @property
    def R_inf_fit(self) -> Optional[float]:
        """Fitted R_s [Ohm], also when it was not used."""
        return float(self.fit.params_opt[0]) if self.fit is not None else None

    @property
    def R_inf_stderr(self) -> Optional[float]:
        """Standard error of R_s [Ohm]: an identifiability flag, not a CI.

        Under model mismatch it understates the real error.
        """
        return float(self.fit.params_stderr[0]) if self.fit is not None else None

    @property
    def R_inf_range(self) -> Tuple[float, float]:
        """Interval meant to hold the true R_s [Ohm]: (0, R_inf_upper) for a
        bound; around a fitted R_s, RINF_BOUND_NOISE_SIGMAS x its stderr,
        but at least +-RINF_REL_STDERR_MAX, as the stderr understates the
        error under model mismatch.

        In the stress test the bound held R_s every time, a fitted range
        not always: where the arc is still open at f_max (phase down to
        -82 deg) and the noise is low (<= 1 %), the window model's error
        can exceed the +-5 %. R_s was outside in 39 of 1417 spectra, the
        fit 6-34 % off, mostly below R_s (doc/STRESS_TEST.md, known
        limit 7)."""
        if self.method == 'rlq_fit' and self.fit is not None:
            shift = max(RINF_REL_STDERR_MAX * abs(self.R_inf),
                        RINF_BOUND_NOISE_SIGMAS * float(self.fit.params_stderr[0]))
            return self.R_inf - shift, self.R_inf + shift
        return 0.0, self.R_inf_upper

    @property
    def L(self) -> float:
        """Series inductance to subtract together with R_inf [H]: the window
        fit's when it determined R_inf, else 0 - a rejected fit's L beside
        the HF bound would be two inconsistent corrections."""
        if self.method == 'rlq_fit' and self.fit is not None:
            return float(self.fit.params_opt[1])
        return 0.0


def hf_median(frequencies: NDArray, Z: NDArray) -> Tuple[float, int]:
    """Median of Re(Z) over the highest frequencies.

    Returns ``(R_inf, n_points)``; see HF_MEDIAN_MAX_POINTS.
    """
    n = min(HF_MEDIAN_MAX_POINTS, max(1, len(frequencies) // 10))
    idx = np.argsort(frequencies)[-n:]
    return float(np.median(Z.real[idx])), n


def _hf_bound(frequencies: NDArray, Z: NDArray) -> Tuple[float, float, int]:
    """Re(Z) at the highest frequency with Im(Z) <= 0, f_max if there is none.

    Returns ``(R, f, index)``. Every point of the R-L-(R|Q) model is an upper bound
    of R_s; this one is the first below an inductive top, where lead
    artifacts can pull Re(Z) below R_s (redoxED flow cell: 0.168 Ohm at the
    Im = 0 crossing, 0.003 and -0.064 Ohm above it). On a capacitive top it
    is f_max, the tightest bound.
    """
    order = np.argsort(frequencies)[::-1]
    capacitive = np.flatnonzero(Z.imag[order] <= 0)
    i = order[capacitive[0]] if capacitive.size else order[0]
    return float(Z.real[i]), float(frequencies[i]), int(i)


def _window_range(f: NDArray, Z: NDArray) -> Dict[str, Tuple[float, float]]:
    """Range per parameter type, each term's impedance within
    min|Z|/RANGE..max|Z|*RANGE.

    R: the value itself. L: its reactance omega*L across the window. Q: its
    impedance 1/(Q omega^n) across the window and the whole n range. All scale
    with Z, so R_inf does not depend on the units.
    """
    Z_abs = np.abs(Z)
    lo, hi = float(Z_abs.min()) / RINF_BOUND_RANGE, float(Z_abs.max()) * RINF_BOUND_RANGE
    w_lo, w_hi = 2 * np.pi * float(f.min()), 2 * np.pi * float(f.max())
    n_lo, n_hi = PARAMETER_BOUNDS['n']
    w_n = [w**n for w in (w_lo, w_hi) for n in (n_lo, n_hi)]
    return {
        'R': (lo, hi),
        'L': (lo / w_hi, hi / w_lo),
        'Q': (1 / (hi * max(w_n)), 1 / (lo * min(w_n))),
        'n': (n_lo, n_hi),
    }


def _initial_circuit(f: NDArray, Z: NDArray, bounds: Dict[str, Tuple[float, float]]):
    """R - L - (R|Q) with a data-driven start, clipped into the range `bounds`.

    R_s: the lowest Re(Z) in the window (Re >= R_s for every term of the model).
    R_k: the spread of Re(Z). tau: 1/omega at the -Im(Z) maximum, which the
    (R|Q) arc peaks at. L: Im/omega at f_max if the top point is inductive.
    """
    def clip(value: float, label: str) -> float:
        return float(np.clip(value, *bounds[label]))

    omega = 2 * np.pi * f
    i_top = int(np.argmax(f))
    R_s = clip(Z.real.min(), 'R')
    R_k = clip(np.ptp(Z.real), 'R')
    tau = 1.0 / omega[int(np.argmax(-Z.imag))]
    if Z.imag[i_top] > 0:
        L_val = Z.imag[i_top] / omega[i_top]
    else:
        L_val = _L_GUESS_REACTANCE_SHARE * abs(Z[i_top]) / omega[i_top]
    return (R(R_s) - L(clip(L_val, 'L'))
            - (R(R_k) | Q(clip(tau**_N_GUESS / R_k, 'Q'), _N_GUESS)))


def estimate_rinf(frequencies: NDArray, Z: NDArray) -> RinfResult:
    """
    Estimate R_inf by an R-L-(R|Q) fit over the top RINF_FIT_DECADES decades.

    Parameters
    ----------
    frequencies : array_like of float
        Frequencies [Hz]
    Z : array_like of complex
        Complex impedance [Ohm]

    Returns
    -------
    RinfResult
        `R_inf` is the fitted R_s if its relative stderr is at most
        RINF_REL_STDERR_MAX, otherwise the HF upper bound of `_hf_bound`
        (reason in `warnings`).

    Raises
    ------
    ValueError
        If the arrays differ in length or hold no finite point.

    Notes
    -----
    No data-only method can tell an arc lying entirely above f_max from a
    flat high-frequency end; such spectra end with the HF upper bound.
    Measured against the 5-point HF median it replaces (audit cases,
    1 % noise): open CPE arc +2483 % instead of +3295 %,
    real_gamry_example.DTA 826 instead of 1402 Ohm; flat ends where the fit
    is flagged -0.7 % instead of ~0 %.
    """
    frequencies = np.asarray(frequencies, dtype=float)
    Z = np.asarray(Z, dtype=np.complex128)
    if frequencies.shape != Z.shape:
        raise ValueError(f"frequencies and Z differ in shape: "
                         f"{frequencies.shape} vs {Z.shape}")
    finite = np.isfinite(frequencies) & np.isfinite(Z) & (frequencies > 0)
    if not finite.any():
        raise ValueError("No finite data point with frequency > 0")
    frequencies, Z = frequencies[finite], Z[finite]

    window = frequencies >= frequencies.max() / 10**RINF_FIT_DECADES
    f_win, Z_win = frequencies[window], Z[window]

    R_hf, f_hf, i_hf = _hf_bound(frequencies, Z)
    # R_s >= 0, so a negative bound is no bound: the point is dominated by
    # noise (Re Z < 0 where |noise| > Re Z) or an artifact. 0 is then the
    # tightest bound the data support; R_inf_hf keeps the measured value.
    # Until a fit reads the noise: a bound <= 0 bounds nothing; every other
    # Re(Z) of the window is a bound too (the loosest the largest), and |Z|
    # at f_hf the scale of its noise - the larger of the two
    upper = R_hf if R_hf > 0 else max(float(np.max(Z_win.real)), float(abs(Z[i_hf])))
    result = RinfResult(R_inf=max(R_hf, 0.0), method='hf_bound', R_inf_hf=R_hf,
                        R_inf_upper=upper, f_hf=f_hf, f_window=f_win, Z_window=Z_win)
    if not finite.all():
        result.warnings.append(f"Ignored {int((~finite).sum())} non-finite point(s)")
    fallback = f"using the HF upper bound Re(Z) = {R_hf:.4g} Ohm at {f_hf:.3g} Hz"
    if f_hf < frequencies.max():
        fallback += " (inductive top above Im(Z) = 0 skipped)"
    if R_hf < 0:
        fallback += ("; it is negative (noise or an artifact exceeds R_s there), "
                     "so R_inf = 0")
    if len(f_win) < RINF_FIT_MIN_POINTS:
        result.warnings.append(
            f"Only {len(f_win)} point(s) in the top {RINF_FIT_DECADES} decades "
            f"(need >= {RINF_FIT_MIN_POINTS}); {fallback}")
        return result

    if not np.min(np.abs(Z_win)) > 0:
        result.warnings.append(f"Zero impedance in the top {RINF_FIT_DECADES} decades; {fallback}")
        return result

    window_range = _window_range(f_win, Z_win)
    circuit = _initial_circuit(f_win, Z_win, window_range)
    labels = circuit.get_param_labels()
    lower = [0.0 if label in ('R', 'L') else window_range[label][0] for label in labels]
    upper = [window_range[label][1] for label in labels]
    try:
        fit, Z_fit = fit_equivalent_circuit(f_win, Z_win, circuit, bounds=(lower, upper))
    except RuntimeError as e:
        result.warnings.append(f"R-L-(R|Q) fit failed ({e}); {fallback}")
        return result
    result.fit = fit
    # Noise per component at f_hf: RMS of the residuals near it (see
    # RINF_NOISE_LOCAL_DECADES; at least the 3 nearest points), scaled up for
    # the 2N residuals that the fit's parameters absorb
    local = np.abs(np.log10(f_win / f_hf)) <= RINF_NOISE_LOCAL_DECADES
    if local.sum() < 3:
        local = np.argsort(np.abs(np.log10(f_win / f_hf)))[:3]
    res = (Z_win - Z_fit)[local]
    dof = 2 * len(f_win) / max(2 * len(f_win) - len(fit.params_opt), 1)
    noise = float(np.sqrt(np.mean(np.abs(res) ** 2) / 2 * dof))
    result.R_inf_upper = max(R_hf + RINF_BOUND_NOISE_SIGMAS * noise, 0.0)

    R_fit, stderr = float(fit.params_opt[0]), float(fit.params_stderr[0])
    rel = stderr / R_fit
    # `not <=` also catches a NaN stderr (covariance could not be computed).
    if not rel <= RINF_REL_STDERR_MAX:
        result.warnings.append(
            f"R_inf is not determined by the top {RINF_FIT_DECADES} decades: "
            f"fit gives {R_fit:.4g} +- {stderr:.2g} Ohm "
            f"({100 * rel:.3g} % > {100 * RINF_REL_STDERR_MAX:.0f} %); "
            f"{fallback}")
        return result

    result.R_inf, result.method = R_fit, 'rlq_fit'
    return result
