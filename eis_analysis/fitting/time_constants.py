"""
Time constants of the (R|C) and (R|Q) arcs of a fitted circuit, G in place
of R included, and the sum of the resistances of arcs in series.

The fit reports R, C, Q and n; the time constant of an arc - what a reader
of an impedance spectrum asks for first - follows from them, and so does its
uncertainty, through the covariance the fit already has. The same covariance
gives the uncertainty of R1 + R2: for close arcs R1 and R2 are correlated
near -1, the split is poorly determined and the sum is not.

The same time constant puts interchangeable series blocks of a fit into one
order, see block_order.py.

CONTRACT

Callers. cli/handlers/fit_output.py and io/export.py call arc_time_constants
and arc_resistance_sums with a finished FitResult; the public API exports
both. block_order.py calls _arc_resistor and _arc_tau on one block.
From the FitResult only these are read, all aligned by parameter index:
  circuit       the circuit tree (Series, Parallel, elements); its
                parameters, depth-first in expression order, each element
                contributing len(get_all_params()), are params_opt in order
  params_opt    fitted values, fixed ones included
  param_labels  a label per parameter, or None
  cov           (n_params, n_params) over all parameters, or None; positive
                semi-definite (covariance.py builds it so); a fixed parameter
                always has a zero row and column, a non-identifiable free one
                (or every free one, with too few residuals) an inf row and
                column
  bound_status  one entry per parameter, '', 'lower', 'upper' or 'fixed',
                or None (no bound information: nothing counts as at a bound)
Nothing else of the fit (stderr, CI, series_split) is used, and nothing
passed in is modified.

What an arc is: a Parallel of exactly two elements, one R or G and one C or
Q, in either order, that is not itself a direct branch of another Parallel
(there its R and C sit beside the other branches). Anything else is not an
arc: a K (it carries tau as a parameter), (R|C|Q), (R|R), C beside a branch
(C_dl | (R_ct - W)). An arc inside a Series that is a branch of a Parallel
is still an arc; its tau is its own, not where the whole spectrum peaks.
R of a G arc is 1/G.

arc_time_constants(result) -> list of ArcTimeConstant, one per arc that has
a time constant, in expression order; empty when param_labels is None.
  tau          R*C, or (R*Q)^(1/n) for Q; exists only when R (or G), C (or Q)
               and n are finite and > 0 and tau and 1/tau are finite and > 0
  frequency    1/(2 pi tau)
  elements     '<label of R or G>|<label of C or Q>', the resistive one first
               whatever the order in the expression
  tau_stderr   delta method in log space with the full covariance block of
               the arc's parameters (R or G, C or Q, n):
               tau * sqrt(g^T Sigma g), g = gradient of ln tau
               (ln R + ln C, or (ln R + ln Q)/n; ln R = -ln G)
  stderr_note  None exactly when tau_stderr is not None

arc_resistance_sums(result) -> list of ArcResistanceSum, one per series
chain with at least two arcs, in expression order; empty when param_labels
is None. A chain is the set of direct members of one Series, a Series
directly inside a Series belonging to the same chain; a Series inside a
Parallel branch is a chain of its own. An arc counts in its chain whether or
not it has a time constant. Plain series elements (Rs, L) and K are not
part of the sum.
  R            sum of the arc resistances (R, or 1/G) as fitted; a chain
               with a non-finite term parameter or a non-finite sum (G = 0)
               gives no entry
  elements     the terms by label, ' + ' joined, a G as '1/<label>'
  stderr       delta method with the full covariance block of the terms:
               sqrt(sum of all entries of D Sigma D), D = diag(gradient),
               d/dR = 1, d(1/G)/dG = -1/G^2
  stderr_parts the stderr of each term alone, in the order of elements
  correlation  for exactly two terms with both parts > 0: (D Sigma D)_12
               scaled by the parts, the correlation of the two terms (R or
               1/G, so a G term flips the sign of cov(R, G)); else None
  stderr_note  None exactly when stderr is not None
  split_less_certain  stderr < min(stderr_parts); False without a stderr

What must hold:
  - No covariance (cov None): every stderr None, note 'no covariance'; tau
    and R are still reported.
  - A parameter of the arc (or of the sum) with bound_status 'lower' or
    'upper': that stderr None, the note names the parameter and the bound.
    A parameter at a bound elsewhere in the circuit does not matter.
  - Fixed parameters are exact: they contribute zero variance. An arc with
    every parameter fixed has stderr 0, not None. 'fixed' takes precedence
    over a bound.
  - An inf or NaN entry in the covariance block used (a non-identifiable
    free parameter, too few residuals): stderr None with a note, never inf,
    NaN or a finite number built from it. An inf row of a parameter outside
    the block does not matter.
  - A variance below zero only by rounding (terms cancelling at correlation
    -1) gives stderr 0, for tau as for the sum.
  - Independent R and C (zero off-diagonal): relative tau_stderr is
    sqrt((s_R/R)^2 + (s_C/C)^2); a sum of uncorrelated terms has
    stderr^2 = sum of stderr_parts^2. With correlation -1 and equal parts the
    sum's stderr is 0.
  - Describing the same arc with G = 1/R instead of R, and the covariance
    transformed accordingly, gives the same tau, tau_stderr and sum to
    rounding. Swapping the two elements of an arc changes nothing.
  - Non-finite or extreme values (NaN or inf parameter, n -> 0, R*Q under-
    or overflowing, 1/G overflowing): no tau, no sum or no stderr as above;
    never an exception and no numpy RuntimeWarning.

_arc_resistor(elements, starts) -> (index of the R or G among all
parameters, +1 for R or -1 for G), or None when the two elements are not an
arc pair; the caller passes the elements of a two-element Parallel.
_arc_tau(elements, starts, resistor, params) -> (tau, indices of
[R or G, C or Q(, n)], gradient of ln tau over them) under the same
existence rule as tau above, else None. Both may emit numpy warnings on
extreme values; callers run them under np.errstate.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from .circuit_builder import Parallel, Series
from .circuit_elements import R, G, C, Q

if TYPE_CHECKING:
    # Annotations only: result.py imports this module via series_resistance.py
    # and block_order.py
    from .result import FitResult



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
        Why tau_stderr is None, e.g. 'no covariance' or 'n0 near its upper
        bound'; None when there is a stderr
    frequency : float
        Characteristic frequency 1/(2*pi*tau) [Hz], where -Z'' peaks
    """
    elements: str
    tau: float
    tau_stderr: Optional[float]
    stderr_note: Optional[str]
    frequency: float


@dataclass
class ArcResistanceSum:
    """
    Sum of the resistances of the arcs in one series chain.

    Attributes
    ----------
    elements : str
        The resistances by their labels, e.g. 'R1 + R2' (a G arc as '1/G0')
    R : float
        The sum [Ohm]
    stderr : float or None
        Its standard error [Ohm], propagated from the full covariance; None
        without a usable covariance or with one of the parameters at a bound
    stderr_note : str or None
        Why stderr is None; None when there is a stderr
    stderr_parts : list of float or None
        Standard error of each resistance alone [Ohm], in the order of
        `elements`; with stderr below the smallest of them the split between
        the arcs is less certain than their sum
    correlation : float or None
        Correlation of the two resistances, for exactly two arcs; else None
    """
    elements: str
    R: float
    stderr: Optional[float]
    stderr_note: Optional[str]
    stderr_parts: Optional[List[float]]
    correlation: Optional[float]

    @property
    def split_less_certain(self) -> bool:
        """The sum is better determined than any of its parts (close arcs)."""
        if self.stderr is None or not self.stderr_parts:
            return False
        return self.stderr < min(self.stderr_parts)


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
    return [arc for arc, _, _ in _collect_arcs(result) if arc is not None]


def arc_resistance_sums(result: FitResult) -> List[ArcResistanceSum]:
    """
    Sum of the arc resistances, per series chain with at least two arcs.

    Every (R|C), (R|Q) pair (G for R) counts, grouped by the series chain it
    is a direct member of: only resistances in series with each other add
    up at DC. Also a pair without a time constant: its R is still in series.
    G = 0 makes the sum infinite (the chain blocks DC), and it is then not
    reported. A Series directly inside a Series is the
    same chain (the builder does not always flatten them); one inside a
    parallel branch is a chain of its own.

    The uncertainty is the delta method with the full covariance: d/dR = 1,
    d(1/G)/dG = -1/G^2. With R1 and R2 correlated near -1 (close arcs) it is
    far below the stderr of either.

    Parameters
    ----------
    result : FitResult
        A fit from fit_equivalent_circuit, or .best_result of a multistart
        or DE fit

    Returns
    -------
    list of ArcResistanceSum
        In the order of the circuit expression; empty without two arcs in
        series or without parameter labels
    """
    labels = result.param_labels
    chains: Dict[int, List[Tuple[int, float]]] = {}
    for _, resistor, chain in _collect_arcs(result):
        if chain is not None:
            chains.setdefault(chain, []).append(resistor)
    params = np.asarray(result.params_opt, dtype=float)
    sums = []
    # An extreme G overflows 1/G or its gradient: no sum, or no stderr,
    # rather than numpy warnings in the fit output
    with np.errstate(over='ignore', under='ignore', divide='ignore', invalid='ignore'):
        for resistors in chains.values():
            if len(resistors) >= 2:
                s = _resistance_sum(resistors, params, labels, result)
                if s is not None:
                    sums.append(s)
    return sums


def _collect_arcs(result: FitResult) -> List[Tuple[Optional[ArcTimeConstant],
                                                   Tuple[int, float], Optional[int]]]:
    """
    Every arc with its time constant, resistance (index, sign) and chain.

    One walk for both public functions. The time constant is None where it
    does not exist (C, Q or n out of range, overflow); the resistance still
    counts in the sum. The chain is an id shared by the direct members of
    one Series (and of a Series directly inside it); None for an arc that is
    not a direct member of a Series.
    """
    params = np.asarray(result.params_opt, dtype=float)
    labels = result.param_labels
    if labels is None:
        return []

    found = []
    chain_ids = itertools.count()

    def walk(node, offset: int, in_parallel: bool = False,
             chain: Optional[int] = None) -> int:
        """Collect the arcs under node; return the offset past its parameters.

        `chain` is the chain of node's parent when that is a Series, else None.
        """
        if isinstance(node, (Series, Parallel)):
            inner = None
            if isinstance(node, Series):
                inner = chain if chain is not None else next(chain_ids)
            starts = []
            for elem in node.elements:
                starts.append(offset)
                offset = walk(elem, offset, isinstance(node, Parallel), inner)
            # A Parallel directly inside a Parallel (the builder does not
            # always flatten them) has its R and C beside the other branches
            # too: R*C is not its time constant
            if isinstance(node, Parallel) and len(node.elements) == 2 and not in_parallel:
                resistor = _arc_resistor(node.elements, starts)
                if resistor is not None:
                    arc = _arc(node.elements, starts, resistor, params, labels, result)
                    found.append((arc, resistor, chain))
            return offset
        return offset + len(node.get_all_params())

    # _arc drops a non-finite tau; numpy's overflow and divide warnings on
    # the way there would only clutter the fit output
    with np.errstate(over='ignore', under='ignore', divide='ignore', invalid='ignore'):
        walk(result.circuit, 0)
    return found


def _resistance_sum(resistors, params, labels, result: FitResult) -> Optional[ArcResistanceSum]:
    """The sum of the (index, sign) resistances and its delta-method stderr.

    None when a term parameter or the sum is not finite (1/G of an extreme G).
    """
    idx = [i for i, _ in resistors]
    total = float(sum(params[i] ** sign for i, sign in resistors))
    # G = inf is a finite 1/G = 0 whose variance the gradient would erase
    if not (np.all(np.isfinite(params[idx])) and np.isfinite(total)):
        return None
    # d(R)/dR = 1; d(1/G)/dG = -1/G^2
    grad = np.array([1.0 if sign > 0 else -1.0 / params[i] ** 2 for i, sign in resistors])
    names = [labels[i] if sign > 0 else f"1/{labels[i]}" for i, sign in resistors]
    stderr = parts = correlation = None
    cov, note = _cov_block(idx, labels, result)
    if cov is not None:
        contrib, var = _delta_variance(cov, grad)
        if (np.all(np.isfinite(contrib)) and var >= 0
                and np.all(np.diag(contrib) >= 0)):
            stderr = float(np.sqrt(var))
            parts = [float(s) for s in np.sqrt(np.diag(contrib))]
            if len(idx) == 2 and parts[0] > 0 and parts[1] > 0:
                # + 0.0: an uncorrelated G arc would print as -0.00
                correlation = float(contrib[0, 1] / (parts[0] * parts[1])) + 0.0
        else:
            note = 'covariance not usable'
    return ArcResistanceSum(elements=' + '.join(names), R=total,
                            stderr=stderr, stderr_note=None if stderr is not None else note,
                            stderr_parts=parts, correlation=correlation)


def _arc_resistor(elements, starts) -> Optional[Tuple[int, float]]:
    """Index and sign (+1 R, -1 G) of the resistance of a two-element arc, or None."""
    kinds = [type(e) for e in elements]
    resistive = [k for k in kinds if k in (R, G)]
    if len(resistive) != 1 or not (C in kinds or Q in kinds):
        return None
    # G = 1/R: ln R = -ln G, so only the sign of its gradient changes
    return starts[kinds.index(resistive[0])], (-1.0 if resistive[0] is G else 1.0)


def _arc(elements, starts, resistor: Tuple[int, float], params, labels,
         result: FitResult) -> Optional[ArcTimeConstant]:
    """The time constant of an arc found by _arc_resistor, or None if it has none."""
    found = _arc_tau(elements, starts, resistor, params)
    if found is None:
        return None
    tau, idx, grad = found
    i_R, i_C = idx[0], idx[1]
    stderr, note = _tau_stderr(tau, grad, idx, labels, result)
    return ArcTimeConstant(elements=f"{labels[i_R]}|{labels[i_C]}", tau=tau,
                           tau_stderr=stderr, stderr_note=note,
                           frequency=float(1 / (2 * np.pi * tau)))


def _arc_tau(elements, starts, resistor: Tuple[int, float],
             params) -> Optional[Tuple[float, List[int], NDArray[np.float64]]]:
    """Tau of an arc found by _arc_resistor, its parameter indices [R, C or Q
    (, n)] and the gradient of ln tau over them; None if it has no tau."""
    i_R, sign = resistor
    kinds = [type(e) for e in elements]
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
    # A tiny n or an extreme R*Q under- or overflows exp(): no usable tau.
    # An inf n would pass as tau = 1 s with a zero gradient
    if not (np.all(np.isfinite([params[i] for i in idx]))
            and np.isfinite(tau) and 0 < tau < np.inf and np.isfinite(1 / tau)):
        return None
    return float(tau), idx, grad


def _cov_block(idx: List[int], labels,
               result: FitResult) -> Tuple[Optional[np.ndarray], Optional[str]]:
    """The covariance of the idx parameters, or None and the reason why not."""
    if result.cov is None:
        return None, 'no covariance'
    status = result.bound_status
    # A parameter at a bound has no Jacobian CI (the JSON export drops its
    # stderr for the same reason), so neither has anything built on it. A
    # fixed one is exact: its zero row in the covariance is the right answer.
    if status is not None:
        for i in idx:
            if i < len(status) and status[i] not in ('', 'fixed'):
                return None, f"{labels[i]} near its {status[i]} bound"
    return np.asarray(result.cov)[np.ix_(idx, idx)], None


def _delta_variance(cov, grad) -> Tuple[np.ndarray, float]:
    """The terms cov_ij g_i g_j and their sum, the delta-method variance.

    With terms correlated near -1 they cancel and rounding can leave the sum
    a hair below 0; double precision errs at ~1e-16 of the terms, so a sum
    above -1e-9 of the trace is a rounding zero.
    """
    contrib = cov * np.outer(grad, grad)
    var = float(contrib.sum())
    if -1e-9 * float(np.trace(contrib)) < var < 0:
        var = 0.0
    return contrib, var


def _tau_stderr(tau: float, grad, idx: List[int], labels,
                result: FitResult) -> Tuple[Optional[float], Optional[str]]:
    """Delta-method stderr of tau and None, or None and the reason why not."""
    cov, note = _cov_block(idx, labels, result)
    if cov is None:
        return None, note
    _, var_log = _delta_variance(cov, grad)
    if not np.isfinite(var_log) or var_log < 0:
        return None, 'covariance not usable'
    return float(tau * np.sqrt(var_log)), None
