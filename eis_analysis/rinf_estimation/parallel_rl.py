"""
Parallel R|L at an inductive high-frequency end (Re(Z) rising toward f_max).

A series L adds only Im(Z) and every RC element's Re(Z) falls with
frequency, so R - L - (R|Q) cannot follow a Re(Z) that rises toward f_max,
as skin effect and eddy currents in a wound cell or long leads make it.
`check_parallel_rl` fits R - (R|L) - (R|Q) to the R_inf window and keeps it
when it is decisively better. estimate_rinf calls it after an accepted
R-L-(R|Q) fit; its contract (estimate.py) describes the result.

Clean design: No logging in core functions, all diagnostics returned as data.
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ..fitting import FitResult, L, Q, R
from ..fitting.circuit import _fit_local
from ..fitting.diagnostics import compute_information_criteria
from ..utils.impedance import _top_down
from .window import clip_to, window_bounds, window_range

# Checked when the top RINF_RL_TOP_DECADES of the window hold at least
# RINF_RL_MIN_INDUCTIVE points with Im(Z) > 0: the decade where an R|L
# corner above f_max shows (NMC811: 5-7 points there), and 3 points are the
# fewest that show a rise rather than one noisy point. R - (R|L) - (R|Q) is
# kept when it beats R - L - (R|Q) by RINF_RL_DBIC_MIN in BIC (decisive) and
# its R|L adds at least RINF_RL_MIN_RE_SHARE of Re(Z) at f_max. Measured
# 2026-11-02: 54 NMC811 Li-ion cells (Zenodo 15422339) dBIC 103-135, share
# ~8 %; 285 stress spectra with an inductive top decade (series L only)
# dBIC < 6 except one noise-free case (rc/83, 549), where a poorly
# converged R-L-(R|Q) lost to an R|L with R_L = 18 MOhm, i.e. a series L
# adding 0.15 % - the share keeps such a fit out.
RINF_RL_TOP_DECADES = 1
RINF_RL_MIN_INDUCTIVE = 3
RINF_RL_DBIC_MIN = 10.0
RINF_RL_MIN_RE_SHARE = 0.01

# R_L starts: |Z(f_max)| (a corner near f_max, as on the NMC cells) and
# RINF_RL_FAR_START x max|Z| of the window (near a series L, where the
# R-L-(R|Q) fit starts from); the better fit counts
RINF_RL_FAR_START = 10.0


@dataclass
class ParallelRL:
    """R - (R|L) - (R|Q) fit of the R_inf window that beat R - L - (R|Q)."""
    R_L: float  # [Ohm]
    L: float  # [H]
    R_s: float  # series R of this fit, the one that goes with the R|L [Ohm]
    dbic: float  # BIC of R - L - (R|Q) minus this fit's (modulus weighting)
    re_share: float  # Re of R|L at f_max / Re(Z(f_max))


def check_parallel_rl(f_win: NDArray, Z_win: NDArray, fit: FitResult,
                      Z_fit: NDArray) -> Optional[ParallelRL]:
    """
    R - (R|L) - (R|Q) on the window, kept when it beats the R-L-(R|Q) `fit`.

    Parameters
    ----------
    f_win, Z_win : ndarray
        The window both fits use [Hz], [Ohm] (at least 5 points)
    fit, Z_fit : FitResult, ndarray
        The R-L-(R|Q) fit of that window and its impedance, under the
        bounds of window.py

    Returns
    -------
    ParallelRL or None
        None when the top RINF_RL_TOP_DECADES of the window hold fewer than
        RINF_RL_MIN_INDUCTIVE points with Im(Z) > 0, when both fits raise
        RuntimeError, or when the better one (by BIC, modulus weighting)
        does not beat `fit` by RINF_RL_DBIC_MIN or its R|L adds less than
        RINF_RL_MIN_RE_SHARE of Re(Z) at f_max.

    Notes
    -----
    Starts: `fit`'s values, R_L at |Z(f_max)| (of several points at f_max
    the one with the lowest Re(Z), as hf_bound, so the order of the points
    does not matter) and at RINF_RL_FAR_START x max|Z|, clipped into
    the range of window.py; R and L bounded below by 0, as for `fit`. The series resistance check
    of fit_equivalent_circuit is skipped. All values scale with Z: R_L, L
    and R_s by the factor, dbic and re_share not.
    """
    top = f_win >= f_win.max() / 10**RINF_RL_TOP_DECADES
    if (Z_win.imag[top] > 0).sum() < RINF_RL_MIN_INDUCTIVE:
        return None
    R_s, L_val, R_k, Q_val, n_val = (float(v) for v in fit.params_opt)
    # Of several points at f_max the lowest Re(Z), as hf_bound: the result
    # must not depend on the order of the points
    i_top = int(_top_down(f_win, Z_win)[0])
    ranges = window_range(f_win, Z_win)

    best = None
    for R_L in (abs(Z_win[i_top]), RINF_RL_FAR_START * float(np.max(np.abs(Z_win)))):
        circuit = (R(R_s) - (R(clip_to(R_L, 'R', ranges)) | L(clip_to(L_val, 'L', ranges)))
                   - (R(R_k) | Q(Q_val, n_val)))
        labels = circuit.get_param_labels()
        try:
            # Without the series resistance check of fit_equivalent_circuit:
            # its refits per block would be discarded here
            fit_rl, Z_rl = _fit_local(f_win, Z_win, circuit,
                                      bounds=window_bounds(labels, ranges))
        except RuntimeError:
            continue
        bic_rl = compute_information_criteria(Z_win, Z_rl, 'modulus', len(labels))[2]
        if best is None or bic_rl < best[0]:
            best = (bic_rl, fit_rl)
    if best is None:
        return None
    bic_rl, fit_rl = best
    bic = compute_information_criteria(Z_win, Z_fit, 'modulus', len(fit.params_opt))[2]
    R_s_rl, R_L, L_rl = (float(v) for v in fit_rl.params_opt[:3])
    # Re of R_L | j omega L: R_L (omega L)^2 / (R_L^2 + (omega L)^2)
    x = 2 * np.pi * float(f_win[i_top]) * L_rl
    re_rl = R_L * x**2 / (R_L**2 + x**2) if R_L > 0 else 0.0
    re_top = float(Z_win.real[i_top])
    share = re_rl / re_top if re_top > 0 else 0.0
    dbic = float(bic - bic_rl)
    if not (dbic >= RINF_RL_DBIC_MIN and share >= RINF_RL_MIN_RE_SHARE):
        return None
    return ParallelRL(R_L=R_L, L=L_rl, R_s=R_s_rl, dbic=dbic, re_share=float(share))
