"""
Time constants of the (R|C) and (R|Q) arcs of a fitted circuit, G in place
of R included, and the sum of the resistances of arcs in series.

The fit reports R, C, Q and n; the time constant of an arc - what a reader
of an impedance spectrum asks for first - follows from them, and so does its
uncertainty, through the covariance the fit already has. The same covariance
gives the uncertainty of R1 + R2: for close arcs R1 and R2 are correlated
near -1, the split is poorly determined and the sum is not.

The same time constant puts interchangeable series blocks of a fit into one
order, see canonical_block_order.
"""
from __future__ import annotations

import itertools
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .circuit_builder import Parallel, Series
from .circuit_elements import R, G, C, Q, K

if TYPE_CHECKING:
    # Annotations only: circuit.py imports canonical_block_order from here
    from .circuit import FitResult

logger = logging.getLogger(__name__)


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


def canonical_block_order(
    circuit,
    params: ArrayLike,
    fixed: ArrayLike,
    lower: ArrayLike,
    upper: ArrayLike,
) -> Tuple[NDArray[np.intp], NDArray[np.intp]]:
    """
    Permutation that puts interchangeable series blocks in ascending tau.

    Blocks in one series chain with the same structure, the same fixed
    parameters (mask and values) and the same bounds are interchangeable:
    swapping their values leaves the impedance and the model unchanged, so a
    fit returns them in any order - DE at random, least_squares in the order
    of its start. Sorted by tau, the fastest process (highest frequency)
    comes first, as auto_suggest orders the DRT peaks. This is the order in
    the spectrum, not the identity of a process: where the time constants of
    two processes cross in a series of spectra, the names swap there.

    Recognizing the symmetry is structural; the key is the physical tau:
    R*C or (R*Q)^(1/n) of an (R|C) or (R|Q) arc (G for R), a K's own tau.
    A group in which some block has no tau stays as it is.

    Parameters
    ----------
    circuit : Circuit
        The circuit; only its structure is used
    params : array_like of float
        Full parameter vector (fixed ones included)
    fixed : array_like of bool
        Fixed parameter mask
    lower, upper : array_like of float
        Bounds of the full parameter vector

    Returns
    -------
    perm : ndarray of int
        params[perm] is the vector in canonical order
    free_perm : ndarray of int
        The same over the free parameters only, e.g. for the columns of a
        Jacobian; the fixed masks of swapped blocks are equal, so free
        parameters only move onto free positions
    """
    values = np.asarray(params, dtype=float)
    fixed_arr = np.asarray(fixed, dtype=bool)
    lower_arr = np.asarray(lower, dtype=float)
    upper_arr = np.asarray(upper, dtype=float)
    perm = np.arange(len(values))

    def sort_chain(members: List[Tuple[int, object]]) -> None:
        groups: Dict[tuple, List[Tuple[int, int, object]]] = {}
        for start, node in members:
            n = len(node.get_all_params())  # type: ignore[attr-defined]
            span = slice(start, start + n)
            mask = fixed_arr[span]
            key = (_structure(node), tuple(mask), tuple(values[perm[span]][mask]),
                   tuple(lower_arr[span]), tuple(upper_arr[span]))
            groups.setdefault(key, []).append((start, n, node))
        for blocks in groups.values():
            if len(blocks) < 2:
                continue
            found = [_block_tau(node, values[perm[s:s + n]]) for s, n, node in blocks]
            taus = [t for t in found if t is not None]
            if len(taus) < len(blocks):
                continue
            # ponytail: stable sort on the exact tau, no relative tolerance.
            # Only an exact tie keeps the written order; nearly equal taus
            # are a coin flip in any order, a tolerance would not fix that.
            order = sorted(range(len(blocks)), key=lambda j: taus[j])
            old = perm.copy()
            for (s, n, _), j in zip(blocks, order):
                src = blocks[j][0]
                perm[s:s + n] = old[src:src + n]

    def walk(node, offset: int, chain: Optional[List[Tuple[int, object]]]) -> int:
        """Sort the chains under node; return the offset past its parameters.

        `chain` collects (offset, node) of the members of the series chain
        node is a member of; a Series directly inside a Series is the same
        chain (the builder does not always flatten them). None outside one.
        """
        if isinstance(node, Series):
            members: List[Tuple[int, object]] = [] if chain is None else chain
            for elem in node.elements:
                offset = walk(elem, offset, members)
            if chain is None:
                sort_chain(members)
            return offset
        if chain is not None:
            chain.append((offset, node))
        if isinstance(node, Parallel):
            # Inner chains first: a member's values are then already in order
            for elem in node.elements:
                offset = walk(elem, offset, None)
            return offset
        return offset + len(node.get_all_params())

    # _arc_tau drops a non-finite tau; its numpy warnings on the way there
    # would only clutter the fit output
    with np.errstate(over='ignore', under='ignore', divide='ignore', invalid='ignore'):
        walk(circuit, 0, None)

    free = np.flatnonzero(~fixed_arr)
    free_perm = np.searchsorted(free, perm[free])
    if not np.array_equal(perm, np.arange(len(perm))):
        logger.debug(f"Interchangeable blocks put in ascending tau: permutation {perm.tolist()}")
    return perm, free_perm


def _structure(node) -> object:
    """The tree of element types under node, without values."""
    if isinstance(node, (Series, Parallel)):
        return type(node).__name__, tuple(_structure(e) for e in node.elements)
    return type(node).__name__


def _block_tau(node, params: NDArray[np.float64]) -> Optional[float]:
    """Tau of an (R|C), (R|Q) (G for R) or K block from its own parameters, or None."""
    if isinstance(node, K):
        tau = params[1]
        return float(tau) if np.isfinite(tau) and tau > 0 else None
    if isinstance(node, Parallel) and len(node.elements) == 2:
        starts = [0, len(node.elements[0].get_all_params())]
        resistor = _arc_resistor(node.elements, starts)
        if resistor is not None:
            found = _arc_tau(node.elements, starts, resistor, params)
            return found[0] if found is not None else None
    return None


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

    None when the sum itself is not finite (1/G of an extreme G).
    """
    idx = [i for i, _ in resistors]
    total = float(sum(params[i] ** sign for i, sign in resistors))
    if not np.isfinite(total):
        return None
    # d(R)/dR = 1; d(1/G)/dG = -1/G^2
    grad = np.array([1.0 if sign > 0 else -1.0 / params[i] ** 2 for i, sign in resistors])
    names = [labels[i] if sign > 0 else f"1/{labels[i]}" for i, sign in resistors]
    stderr = parts = correlation = None
    cov, note = _cov_block(idx, labels, result)
    if cov is not None:
        contrib = cov * np.outer(grad, grad)
        var = float(contrib.sum())
        # With R1, R2 correlated near -1 the terms cancel and rounding can
        # leave var a hair below 0; double precision errs at ~1e-16 of the
        # terms, a non-PSD block by the size of the terms themselves
        if -1e-9 * float(np.trace(contrib)) < var < 0:
            var = 0.0
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
    # A tiny n or an extreme R*Q under- or overflows exp(): no usable tau
    if not (np.isfinite(tau) and 0 < tau < np.inf and np.isfinite(1 / tau)):
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


def _tau_stderr(tau: float, grad, idx: List[int], labels,
                result: FitResult) -> Tuple[Optional[float], Optional[str]]:
    """Delta-method stderr of tau and None, or None and the reason why not."""
    cov, note = _cov_block(idx, labels, result)
    if cov is None:
        return None, note
    var_log = float(grad @ cov @ grad)
    if not np.isfinite(var_log) or var_log < 0:
        return None, 'covariance not usable'
    return float(tau * np.sqrt(var_log)), None
