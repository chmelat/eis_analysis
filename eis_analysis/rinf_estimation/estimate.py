"""
R_inf estimation: R-L-(R|Q) fit over the top frequency decades.

Model: Z(omega) = R_s + j*omega*L + R_k / (1 + R_k*Q*(j*omega)^n)

One nonlinear fit covers inductive, capacitive and mixed high-frequency ends
alike; R_inf = R_s. The standard error of R_s decides whether the window
determines R_inf at all. When it does not, or when the fit cannot run, the
result falls back to the HF median and says why in `warnings`.

Clean design: No logging in core functions, all diagnostics returned as data.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from ..fitting import FitResult, L, Q, R, fit_equivalent_circuit
from ..fitting.bounds import PARAMETER_BOUNDS

# Fit window: f >= f_max / 10**RINF_FIT_DECADES. Two decades is the compromise
# measured in doc/AUDIT_ri_fit_2026-09-25.md: one decade resolves an arc just
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

# Initial n of the CPE: midway in the typical 0.6-1.0 range of real arcs.
_N_GUESS = 0.8

# Initial L when the top point is capacitive (Im <= 0): a fraction of the
# typical cable inductance, so the fit starts near "no inductance".
_L_GUESS_CAPACITIVE = 1e-9  # [H]


@dataclass
class RinfResult:
    """R_inf estimate with the fit and the HF median it was chosen from.

    `R_inf` is the value to use: the fitted R_s when the window determines
    it (`method == 'rlq_fit'`), otherwise the HF median
    (`method == 'hf_median'`) with the reason in `warnings`.
    """
    R_inf: float  # [Ohm]
    method: str  # 'rlq_fit' | 'hf_median'
    R_inf_median: float  # [Ohm]
    n_median_points: int
    f_window: NDArray[np.float64]  # frequencies of the fit window [Hz]
    Z_window: NDArray[np.complex128]  # impedance of the fit window [Ohm]
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


def hf_median(frequencies: NDArray, Z: NDArray) -> Tuple[float, int]:
    """Median of Re(Z) over the highest frequencies.

    Returns ``(R_inf, n_points)``; see HF_MEDIAN_MAX_POINTS.
    """
    n = min(HF_MEDIAN_MAX_POINTS, max(1, len(frequencies) // 10))
    idx = np.argsort(frequencies)[-n:]
    return float(np.median(Z.real[idx])), n


def _clip(value: float, label: str) -> float:
    lo, hi = PARAMETER_BOUNDS[label]
    return float(np.clip(value, lo, hi))


def _initial_circuit(f: NDArray, Z: NDArray):
    """R - L - (R|Q) with a data-driven start.

    R_s: the lowest Re(Z) in the window (Re >= R_s for every term of the model).
    R_k: the spread of Re(Z). tau: 1/omega at the -Im(Z) maximum, which the
    (R|Q) arc peaks at. L: Im/omega at f_max if the top point is inductive.
    """
    omega = 2 * np.pi * f
    i_top = int(np.argmax(f))
    R_s = _clip(Z.real.min(), 'R')
    R_k = _clip(np.ptp(Z.real), 'R')
    tau = 1.0 / omega[int(np.argmax(-Z.imag))]
    L_val = Z.imag[i_top] / omega[i_top] if Z.imag[i_top] > 0 else _L_GUESS_CAPACITIVE
    return (R(R_s) - L(_clip(L_val, 'L'))
            - (R(R_k) | Q(_clip(tau**_N_GUESS / R_k, 'Q'), _N_GUESS)))


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
        RINF_REL_STDERR_MAX, otherwise the HF median (reason in `warnings`).

    Raises
    ------
    ValueError
        If the arrays differ in length or hold no finite point.

    Notes
    -----
    No data-only method can tell an arc lying entirely above f_max from a
    flat high-frequency end; such spectra end with the HF median.
    """
    frequencies = np.asarray(frequencies, dtype=float)
    Z = np.asarray(Z, dtype=complex)
    if frequencies.shape != Z.shape:
        raise ValueError(f"frequencies and Z differ in shape: "
                         f"{frequencies.shape} vs {Z.shape}")
    finite = np.isfinite(frequencies) & np.isfinite(Z) & (frequencies > 0)
    if not finite.any():
        raise ValueError("No finite data point with frequency > 0")
    frequencies, Z = frequencies[finite], Z[finite]

    window = frequencies >= frequencies.max() / 10**RINF_FIT_DECADES
    f_win, Z_win = frequencies[window], Z[window]

    R_median, n_median = hf_median(frequencies, Z)
    result = RinfResult(R_inf=R_median, method='hf_median', R_inf_median=R_median,
                        n_median_points=n_median, f_window=f_win, Z_window=Z_win)
    if not finite.all():
        result.warnings.append(f"Ignored {int((~finite).sum())} non-finite point(s)")
    if len(f_win) < RINF_FIT_MIN_POINTS:
        result.warnings.append(
            f"Only {len(f_win)} point(s) in the top {RINF_FIT_DECADES} decades "
            f"(need >= {RINF_FIT_MIN_POINTS}); using HF median")
        return result

    try:
        fit, _, _ = fit_equivalent_circuit(f_win, Z_win, _initial_circuit(f_win, Z_win),
                                           plot=False)
    except RuntimeError as e:
        result.warnings.append(f"R-L-(R|Q) fit failed ({e}); using HF median")
        return result
    result.fit = fit

    R_fit, stderr = float(fit.params_opt[0]), float(fit.params_stderr[0])
    rel = stderr / R_fit
    # `not <=` also catches a NaN stderr (covariance could not be computed).
    if not rel <= RINF_REL_STDERR_MAX:
        result.warnings.append(
            f"R_inf is not determined by the top {RINF_FIT_DECADES} decades: "
            f"fit gives {R_fit:.4g} +- {stderr:.2g} Ohm "
            f"({100 * rel:.3g} % > {100 * RINF_REL_STDERR_MAX:.0f} %); using HF median")
        return result

    result.R_inf, result.method = R_fit, 'rlq_fit'
    return result
