"""
Canonical order of interchangeable series blocks, and the series chain walk
and arc recognition it shares with the series resistance check.

Equal (R|C), (R|Q) (G for R) or K blocks in series swap values without
changing the impedance; canonical_block_order puts them in ascending tau.
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .circuit_builder import Parallel, Series
from .circuit_elements import K
from .time_constants import _arc_resistor, _arc_tau

logger = logging.getLogger(__name__)


def canonical_block_order(
    circuit,
    params: ArrayLike,
    fixed: ArrayLike,
    lower: ArrayLike,
    upper: ArrayLike,
) -> Tuple[NDArray[np.intp], NDArray[np.intp], bool]:
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

    The order of elements inside a block does not matter: (C|R) is the same
    block as (R|C). Recognizing the symmetry is structural; the key is the
    physical tau:
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
    grouped : bool
        True if the circuit has a group of interchangeable blocks, also
        when they already were in order
    """
    values = np.asarray(params, dtype=float)
    fixed_arr = np.asarray(fixed, dtype=bool)
    lower_arr = np.asarray(lower, dtype=float)
    upper_arr = np.asarray(upper, dtype=float)
    perm = np.arange(len(values))
    grouped = False

    def sort_chain(members: List[Tuple[int, Tuple[int, ...], Any]]) -> None:
        nonlocal grouped
        groups: Dict[tuple, List[Tuple[int, NDArray[np.intp], Any]]] = {}
        for start, _, node in members:
            structure, idx = _structure(node)
            # The block's parameters in the order of its structure, so that
            # the n-th slot of (R|C) and of (C|R) is the same parameter
            slots = start + np.asarray(idx, dtype=np.intp)
            mask = fixed_arr[slots]
            key = (structure, tuple(mask), tuple(values[perm[slots]][mask]),
                   tuple(lower_arr[slots]), tuple(upper_arr[slots]))
            groups.setdefault(key, []).append((start, slots, node))
        for blocks in groups.values():
            if len(blocks) < 2:
                continue
            arcs = [_block_arc(node, values[perm[s:s + len(slots)]]) for s, slots, node in blocks]
            taus = [a[2] for a in arcs if a is not None and a[2] is not None]
            if len(taus) < len(blocks):
                continue
            grouped = True
            # ponytail: stable sort on the exact tau, no relative tolerance.
            # Only an exact tie keeps the written order; nearly equal taus
            # are a coin flip in any order, a tolerance would not fix that.
            order = sorted(range(len(blocks)), key=lambda j: taus[j])
            old = perm.copy()
            for (_, slots, _), j in zip(blocks, order):
                perm[slots] = old[blocks[j][1]]

    # Inner chains first: a member's values are then already in order
    for members in _series_chains(circuit):
        sort_chain(members)

    free = np.flatnonzero(~fixed_arr)
    free_perm = np.searchsorted(free, perm[free])
    if not np.array_equal(perm, np.arange(len(perm))):
        logger.debug(f"Interchangeable blocks put in ascending tau: permutation {perm.tolist()}")
    return perm, free_perm, grouped


def _structure(node) -> Tuple[tuple, List[int]]:
    """The tree of element types under node, without values, and the offsets
    of its parameters in the same order.

    Series and parallel are commutative, so the children are sorted by their
    own structure: (R|C) and (C|R) give one tree, and the offsets map the
    parameters of either onto the same slots. Equal children keep their order.
    """
    if not isinstance(node, (Series, Parallel)):
        return (type(node).__name__,), list(range(len(node.get_all_params())))
    children = []
    offset = 0
    for elem in node.elements:
        structure, idx = _structure(elem)
        children.append((structure, [offset + i for i in idx]))
        offset += len(elem.get_all_params())
    children.sort(key=lambda child: child[0])
    return ((type(node).__name__, tuple(s for s, _ in children)),
            [i for _, idx in children for i in idx])


def _series_chains(circuit) -> List[List[Tuple[int, Tuple[int, ...], Any]]]:
    """Members (offset, path, node) of every series chain under circuit,
    inner chains first, the circuit's own chain last.

    A Series directly inside a Series is the same chain (the builder does not
    always flatten them). `path` indexes .elements down to the node.
    """
    chains = []

    def walk(node, offset: int, path: Tuple[int, ...], chain) -> int:
        if isinstance(node, Series):
            members = [] if chain is None else chain
            for i, elem in enumerate(node.elements):
                offset = walk(elem, offset, path + (i,), members)
            if chain is None:
                chains.append(members)
            return offset
        if chain is not None:
            chain.append((offset, path, node))
        if isinstance(node, Parallel):
            for i, elem in enumerate(node.elements):
                offset = walk(elem, offset, path + (i,), None)
            return offset
        return offset + len(node.get_all_params())

    walk(circuit, 0, (), None)
    return chains


def _block_arc(node, params: NDArray[np.float64]
               ) -> Optional[Tuple[int, Optional[float], Optional[float]]]:
    """(index of R in the node, R, tau) of an (R|C), (R|Q) (G for R, R = 1/G)
    or K block from its own parameters; R or tau None where not finite and
    positive. None for any other block."""
    if isinstance(node, K):
        return 0, _positive(params[0]), _positive(params[1])
    if not (isinstance(node, Parallel) and len(node.elements) == 2):
        return None
    starts = [0, len(node.elements[0].get_all_params())]
    resistor = _arc_resistor(node.elements, starts)
    if resistor is None:
        return None
    i, sign = resistor
    # _arc_tau drops a non-finite tau; numpy's warnings on the way there
    # would only clutter the fit output
    with np.errstate(over='ignore', under='ignore', divide='ignore', invalid='ignore'):
        found = _arc_tau(node.elements, starts, resistor, params)
        return i, _positive(params[i] ** sign), found[0] if found is not None else None


def _positive(value) -> Optional[float]:
    """value as float if finite and > 0, else None."""
    return float(value) if np.isfinite(value) and value > 0 else None


__all__ = ['canonical_block_order']
