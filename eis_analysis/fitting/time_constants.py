"""
Time constants of the (R|C) and (R|Q) arcs of a fitted circuit, G in place
of R included.

The fit reports R, C, Q and n; the time constant of an arc - what a reader
of an impedance spectrum asks for first - follows from them, and so does its
uncertainty, through the covariance the fit already has.
"""

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from .circuit import FitResult
from .circuit_builder import Parallel, Series
from .circuit_elements import R, G, C, Q


@dataclass
class ArcTimeConstant:
    """
    Time constant of one parallel arc, (R|C) or (R|Q).

    Attributes
    ----------
    elements : str
        The arc by its parameter labels, e.g. 'R1|C0' or 'R2|Q0'
    tau : float
        Time constant [s]: R*C, or (R*Q)^(1/n) for a CPE - the peak of the
        ZARC, equal to R*C_eff with the Hsu-Mansfeld C_eff
    tau_stderr : float or None
        Standard error of tau [s], propagated from the fit covariance; None
        without a usable covariance or with one of the parameters at a bound
    stderr_note : str or None
        Why tau_stderr is None, e.g. 'no covariance' or 'n0 at its upper
        bound'; None when there is a stderr
    frequency : float
        Characteristic frequency 1/(2*pi*tau) [Hz], where -Z'' peaks
    """
    elements: str
    tau: float
    tau_stderr: Optional[float]
    stderr_note: Optional[str]
    frequency: float


def arc_time_constants(result: FitResult) -> List[ArcTimeConstant]:
    """
    Time constants of the (R|C) and (R|Q) arcs of a fitted circuit.

    An arc is a parallel combination of exactly one R (or G = 1/R) and one
    C or Q, in either order. A C or Q beside a longer branch is not one: in the Randles
    `C_dl | (R_ct - Ws)` the resistance beside C_dl at DC includes R_W, and
    R*C would not be the charge-transfer time constant. A K carries its tau
    as a parameter already.

    The uncertainty is the delta method in log space, with the full
    covariance (correlated R and C, fixed parameters at zero):
    ln tau = ln R + ln C, or (ln R + ln Q)/n, with ln R = -ln G.

    Parameters
    ----------
    result : FitResult
        A fit from fit_equivalent_circuit, or .best_result of a multistart
        or DE fit

    Returns
    -------
    list of ArcTimeConstant
        In the order of the circuit expression; empty without such an arc
        or without parameter labels
    """
    params = np.asarray(result.params_opt, dtype=float)
    labels = result.param_labels
    if labels is None:
        return []

    arcs = []

    def walk(node, offset: int) -> int:
        """Collect the arcs under node; return the offset past its parameters."""
        if isinstance(node, (Series, Parallel)):
            starts = []
            for elem in node.elements:
                starts.append(offset)
                offset = walk(elem, offset)
            if isinstance(node, Parallel) and len(node.elements) == 2:
                arc = _arc(node.elements, starts, params, labels, result)
                if arc is not None:
                    arcs.append(arc)
            return offset
        return offset + len(node.get_all_params())

    walk(result.circuit, 0)
    return arcs


def _arc(elements, starts, params, labels, result: FitResult) -> Optional[ArcTimeConstant]:
    """The arc of a two-element Parallel, or None if it is not R or G with C or Q."""
    kinds = [type(e) for e in elements]
    resistive = [k for k in kinds if k in (R, G)]
    if len(resistive) != 1 or not (C in kinds or Q in kinds):
        return None
    i_R = starts[kinds.index(resistive[0])]
    # G = 1/R: ln R = -ln G, so only the sign of its gradient changes
    sign = -1.0 if resistive[0] is G else 1.0
    R_val = params[i_R] ** sign
    if C in kinds:
        i_C = starts[kinds.index(C)]
        idx = [i_R, i_C]
        C_val = params[i_C]
        if not (params[i_R] > 0 and C_val > 0):
            return None
        tau = R_val * C_val
        grad = np.array([sign / params[i_R], 1 / C_val])
    else:
        i_C = starts[kinds.index(Q)]
        idx = [i_R, i_C, i_C + 1]
        Q_val, n = params[i_C], params[i_C + 1]
        if not (params[i_R] > 0 and Q_val > 0 and n > 0):
            return None
        log_RQ = np.log(R_val) + np.log(Q_val)
        tau = np.exp(log_RQ / n)
        grad = np.array([sign / (n * params[i_R]), 1 / (n * Q_val), -log_RQ / n**2])
    # A tiny n or an extreme R*Q under- or overflows exp(): no usable tau
    if not (np.isfinite(tau) and 0 < tau < np.inf and np.isfinite(1 / tau)):
        return None
    stderr, note = _tau_stderr(float(tau), grad, idx, labels, result)
    return ArcTimeConstant(elements=f"{labels[i_R]}|{labels[i_C]}", tau=float(tau),
                           tau_stderr=stderr, stderr_note=note,
                           frequency=float(1 / (2 * np.pi * tau)))


def _tau_stderr(tau: float, grad, idx: List[int], labels,
                result: FitResult) -> Tuple[Optional[float], Optional[str]]:
    """Delta-method stderr of tau and None, or None and the reason why not."""
    if result.cov is None:
        return None, 'no covariance'
    status = result.bound_status
    # A parameter at a bound has no Jacobian CI (the JSON export drops its
    # stderr for the same reason), so neither has anything built on it. A
    # fixed one is exact: its zero row in the covariance is the right answer.
    if status is not None:
        for i in idx:
            if status[i] not in ('', 'fixed'):
                return None, f"{labels[i]} at its {status[i]} bound"
    cov = np.asarray(result.cov)[np.ix_(idx, idx)]
    var_log = float(grad @ cov @ grad)
    if not np.isfinite(var_log) or var_log < 0:
        return None, 'covariance not usable'
    return float(tau * np.sqrt(var_log)), None
