"""
Series resistance that the data determine only together with an arc.

An arc - (R|C), (R|Q) (G for R), K - whose time constant lies far above the
measured frequencies adds its R as a constant: the data fix only the sum
Rs + R_k, and every split of it fits equally well. The fit can land
elsewhere, with the arc as a small feature fitting the noise; there the
covariance gives Rs a tight interval, because the branch on which the arc
degenerates to a plain resistance (C or Q -> 0) is not reachable along a
straight line from the fit. On 47 such AutoECM test spectra with 1 % noise
the 95 % CI of Rs held the truth in 22, with this check in 42.

The check refits with each arc merged into Rs. Worse by less than the
chi2(1) quantile, the profile of Rs is flat from its lower bound up to the
merged value (on the branch the arc's R takes any share of the sum), and
FitResult widens the CI of Rs to that range.
"""
from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.stats import chi2

from .bounds import generate_simple_bounds
from .circuit_builder import Parallel, Series
from .circuit_elements import K, R
from .diagnostics import compute_information_criteria
from .time_constants import _arc_resistor

if TYPE_CHECKING:
    # Annotations only: circuit.py imports this module
    from .circuit import FitResult

logger = logging.getLogger(__name__)

# A split is kept up to the widest CI FitResult offers (params_ci_99); the
# CI of each level then widens where dchi2 is below its own chi2(1) quantile
SPLIT_KEEP_LEVEL = 0.99
# The warning goes with the reported interval, params_ci_95
SPLIT_WARN_LEVEL = 0.95


@dataclass
class SeriesResistanceSplit:
    """
    An arc that can be merged into the series resistance at little cost.

    Attributes
    ----------
    index : int
        Index of the series resistance in params_opt
    arc_resistance : str
        Label of the arc's resistance (R1, G0, the R of a K)
    dchi2 : float
        (RSS_merged - RSS) / s^2 with s^2 = RSS / dof of the fit; the profile
        of the series resistance on the branch where the arc is a plain
        resistance. Negative when the merged fit is the better one.
    lower : float
        Lower bound of the series resistance, where that branch ends
    merged : FitResult
        Fit with the arc merged into the series resistance
    merged_index : int
        Index of the series resistance in merged.params_opt
    arc_index : int or None
        Index of the arc's R in params_opt, which takes the rest of the sum
        on that branch; None for a G
    arc_lower : float
        Lower bound of the arc's R
    """
    index: int
    arc_resistance: str
    dchi2: float
    lower: float
    merged: 'FitResult'
    merged_index: int
    arc_index: Optional[int] = None
    arc_lower: float = 0.0

    def widens(self, confidence_level: float) -> bool:
        """Whether the CI at this level spans the branch (dchi2 below chi2(1))."""
        return self.dchi2 < chi2.ppf(confidence_level, 1)


def series_resistance_split(
    result: 'FitResult',
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    weighting: str = 'modulus',
    bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
    use_analytic_jacobian: bool = True,
) -> Optional[SeriesResistanceSplit]:
    """
    The arc whose merging into the series resistance costs least, if cheap.

    Each free arc in the series chain of the fitted circuit is removed and
    its R added to the start of the series resistance, and the rest is
    refitted locally from the fit.

    Parameters
    ----------
    result : FitResult
        The fit; its circuit is not changed
    frequencies, Z : ndarray
        The data the fit ran on
    weighting : str
        Weighting of the fit
    bounds : (lower, upper), optional
        Bounds of the fit, default generate_simple_bounds
    use_analytic_jacobian : bool
        As in the fit

    Returns
    -------
    SeriesResistanceSplit or None
        None without exactly one free series R, or when every merge costs
        more than the chi2(1) quantile at SPLIT_KEEP_LEVEL
    """
    from .circuit import _fit_local  # circuit.py imports this module

    circuit = result.circuit
    members = _series_chain(circuit)
    if not members or result._dof <= 0:
        return None
    params = np.asarray(result.params_opt, dtype=float)
    fixed = circuit.get_all_fixed_params()
    series_r = [offset for offset, _, node in members if type(node) is R and not fixed[offset]]
    # ponytail: exactly one free series R. Two are degenerate with each other
    # already (stderr inf); G or a fixed R in series are not handled.
    if len(series_r) != 1:
        return None
    rs = series_r[0]
    labels = result.param_labels or circuit.get_param_labels()
    if bounds is None:
        bounds = generate_simple_bounds(circuit.get_param_labels())
    lower, upper = np.asarray(bounds[0], dtype=float), np.asarray(bounds[1], dtype=float)

    def rss(c, p) -> float:
        # k = 1: only the RSS is used, not the criteria
        return compute_information_criteria(Z, c.impedance(frequencies, list(p)), weighting, 1)[0]

    rss_fit = rss(circuit, params)
    s2 = rss_fit / result._dof
    if not s2 > 0:
        return None

    raw_labels = circuit.get_param_labels()
    best: Optional[SeriesResistanceSplit] = None
    for offset, path, node in members:
        n = len(node.get_all_params())
        arc = _arc_r(node, params[offset:offset + n])
        if arc is None or any(fixed[offset:offset + n]):
            continue
        i_R, R_k = arc
        # ponytail: the merged fit drops the arc (C or Q = 0) without checking
        # that its bounds reach that branch; a lower bound on C or Q so high
        # that the arc is not constant at f_max overstates the widened CI.
        # Upgrade: refit with C/Q fixed at its lower bound instead of removed.
        merged_circuit = deepcopy(circuit)
        parent: Any = merged_circuit
        for i in path[:-1]:
            parent = parent.elements[i]
        del parent.elements[path[-1]]
        keep = np.r_[0:offset, offset + n:len(params)]
        start = params.copy()
        # The sum can pass the upper bound of Rs; an initial_guess outside
        # the bounds would fail the fit instead of being clipped
        start[rs] = min(start[rs] + R_k, upper[rs])
        try:
            merged, _ = _fit_local(frequencies, Z, merged_circuit, weighting=weighting,
                                   initial_guess=list(start[keep]),
                                   use_analytic_jacobian=use_analytic_jacobian,
                                   bounds=(list(lower[keep]), list(upper[keep])))
        except (RuntimeError, ValueError) as exc:
            logger.debug(f"Merging {labels[offset + i_R]} into {labels[rs]} failed: {exc}")
            continue
        dchi2 = (rss(merged.circuit, merged.params_opt) - rss_fit) / s2
        if best is None or dchi2 < best.dchi2:
            i_arc = offset + i_R
            best = SeriesResistanceSplit(index=rs, arc_resistance=labels[i_arc],
                                         dchi2=float(dchi2), lower=float(lower[rs]),
                                         merged=merged,
                                         merged_index=rs if rs < offset else rs - n,
                                         arc_index=i_arc if raw_labels[i_arc] == 'R' else None,
                                         arc_lower=float(lower[i_arc]))
    if best is None or not best.widens(SPLIT_KEEP_LEVEL):
        return None
    return best


def check_series_resistance(
    result: 'FitResult',
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    weighting: str = 'modulus',
    bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None,
    use_analytic_jacobian: bool = True,
) -> None:
    """Set result.series_split and add its warning (series_resistance_split)."""
    split = series_resistance_split(result, frequencies, Z, weighting, bounds,
                                    use_analytic_jacobian)
    result.series_split = split
    if split is None or not split.widens(SPLIT_WARN_LEVEL) or result.diagnostics is None:
        return
    labels = result.param_labels or result.circuit.get_param_labels()
    name, arc = labels[split.index], split.arc_resistance
    low, high = result.params_ci_95
    if split.dchi2 < 0:
        cost = (f"fits better (dchi2 = {split.dchi2:.2f}), so this fit is not the "
                f"least-squares minimum")
    else:
        cost = f"is worse by dchi2 = {split.dchi2:.2f} (< {chi2.ppf(SPLIT_WARN_LEVEL, 1):.2f})"
    widened = f"{name} and {arc}" if split.arc_index is not None else name
    result.diagnostics.warnings.append(
        f"{name} is determined only as {name} + {arc}: the fit with that arc as a "
        f"plain resistance {cost}; the 95 % CI of {widened} is widened, "
        f"{name} to [{low[split.index]:.3g}, {high[split.index]:.4g}]")


def _series_chain(circuit) -> List[Tuple[int, Tuple[int, ...], Any]]:
    """(offset, path, node) of the members of the top series chain; a Series
    inside it is the same chain. `path` indexes .elements down to the node."""
    if not isinstance(circuit, Series):
        return []
    out = []

    def walk(node: Series, offset: int, path: Tuple[int, ...]) -> int:
        for i, elem in enumerate(node.elements):
            if isinstance(elem, Series):
                offset = walk(elem, offset, path + (i,))
                continue
            out.append((offset, path + (i,), elem))
            offset += len(elem.get_all_params())
        return offset

    walk(circuit, 0, ())
    return out


def _arc_r(node, params: NDArray[np.float64]) -> Optional[Tuple[int, float]]:
    """(index in the node, resistance) of an (R|C), (R|Q) (G for R) or K arc, or None."""
    if isinstance(node, K):
        i, value = 0, params[0]
    elif isinstance(node, Parallel) and len(node.elements) == 2:
        starts = [0, len(node.elements[0].get_all_params())]
        found = _arc_resistor(node.elements, starts)
        if found is None:
            return None
        i, sign = found
        value = params[i] ** sign
    else:
        return None
    return (i, float(value)) if np.isfinite(value) and value > 0 else None


__all__ = ['SeriesResistanceSplit', 'series_resistance_split', 'check_series_resistance']
