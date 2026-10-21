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

The merged fit is as good also when another arc takes the role of the one
removed (two arcs close in tau, a K beside an (R|C) of the same tau). The
profile of Rs is flat there too: the removed arc fits again as a fast plain
resistance, its tau at the lower bound. The warning then points at an arc
the data do not need.

CONTRACT

Callers. fit_equivalent_circuit (circuit.py) after _fit_local with the
fit's bounds (None = generate_simple_bounds), fitting/diffevo.py with the
full bounds of the DE fit, fitting/multistart.py on the chosen start only
with the default bounds (it has no others); all through
check_series_resistance, with the data, weighting and analytic-Jacobian
choice the fit ran on. The public API exports both functions.
From the FitResult it reads circuit (structure and fixed flags),
params_opt, param_labels (else the circuit's own labels), _dof and
diagnostics. bounds, when given, are (lower, upper) over all parameters.

series_resistance_split(result, frequencies, Z, weighting, bounds,
use_analytic_jacobian) -> SeriesResistanceSplit or None.
  Where it looks: only the circuit's own series chain (the direct members
  of the top Series, a Series directly inside it included), not chains
  inside a parallel branch. It needs exactly one free plain R among those
  members (the series resistance Rs); with none, two or more, a free G
  among those members, or the circuit not a Series it returns None.
  Candidates: every member that is an (R|C), (R|Q) (G for R) or K with no
  fixed parameter and a finite R > 0 (1/G for a G); its tau does not
  matter.
  Merge: the candidate is removed from a deep copy of the circuit, Rs starts
  at Rs + R_k (clipped to its upper bound), every other parameter at its
  fitted value, fixed flags kept, and the rest is refitted locally (no DE,
  no multistart, no series check of its own) within the fit's bounds of the
  parameters kept. A refit that raises RuntimeError or ValueError, or whose
  dchi2 is not finite, drops that candidate only.
  dchi2 = (RSS_merged - RSS_fit) / s^2, s^2 = RSS_fit / dof, both RSS with
  the fit's weighting (the objective the fit minimized), dof = result._dof.
  dof <= 0 or s^2 not > 0 (a perfect or non-finite fit): None.
  Result: the candidate with the smallest dchi2 (the first in expression
  order on a tie), returned only if widens(SPLIT_KEEP_LEVEL); else None.
  Two arcs that can take each other's role give merged fits equal up to
  the refit's convergence, so either may win; both are valid.
  Fields: index of Rs in params_opt; arc_resistance the label of the arc's
  R, G or K's R; dchi2 (negative when the merged fit is better: the fit was
  not the minimum); lower the lower bound of Rs; merged the merged
  FitResult; merged_index the index of Rs in merged.params_opt; arc_index
  the index of the arc's R (of an R arc or a K), None for a G; arc_lower
  the lower bound of that R.

SeriesResistanceSplit.widens(level) is dchi2 < chi2.ppf(level, 1): the
profile of Rs stays within the level's quantile along the merged branch.
FitResult._ci (result.py) uses it: where it widens, the CI of Rs runs from
lower to max(its own high, the merged fit's high of Rs at that level), and
so does the CI of the arc's R from arc_lower.

check_series_resistance(...) -> None: sets result.series_split to what
series_resistance_split returns (None included, overwriting an earlier
value) and, when that widens at SPLIT_WARN_LEVEL (0.95) and
result.diagnostics exists, appends exactly one warning naming Rs, the arc's
resistance and the widened 95 % CI of Rs. Within SPLIT_TIE_DCHI2 of 0 the
merged fit fits as well; below -SPLIT_TIE_DCHI2 the warning says the fit is
not the least-squares minimum.

What must hold:
  - The fit is not changed: params_opt, params_stderr, cov, its circuit
    (structure and values) and every other field but series_split and the
    one warning stay as they were.
  - A data set or circuit the check cannot handle gives None, never an
    exception; a failing refit is logged at debug level only.
  - The order of the candidates does not change the result: the smallest
    finite dchi2 wins whatever comes before or after it (ties as above).
  - dchi2 does not depend on the scale of the weights: multiplying every
    weighted residual by a constant leaves it unchanged.
  - widens is monotone in level: widening at 0.95 implies widening at 0.99.
"""
from __future__ import annotations

import logging
from copy import deepcopy
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.stats import chi2

from .bounds import generate_simple_bounds
from .circuit_builder import Series
from .circuit_elements import G, R
from .diagnostics import compute_information_criteria
from .block_order import _block_arc, _series_chains

if TYPE_CHECKING:
    # Annotations only: result.py imports this module
    from .result import FitResult

logger = logging.getLogger(__name__)

# A split is kept up to the widest CI FitResult offers (params_ci_99); the
# CI of each level then widens where dchi2 is below its own chi2(1) quantile
SPLIT_KEEP_LEVEL = 0.99
# The warning goes with the reported interval, params_ci_95
SPLIT_WARN_LEVEL = 0.95
# A merged fit this close to the fit (in dchi2) fits as well: in the
# degenerate case the refits land within +-1e-5 of each other, with either
# sign; far below the quantiles that decide the widening (3.84 and 6.63 at
# SPLIT_WARN_LEVEL and SPLIT_KEEP_LEVEL)
SPLIT_TIE_DCHI2 = 0.01


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
    # The circuit's own chain, the last of them
    members = _series_chains(circuit)[-1] if isinstance(circuit, Series) else []
    if not members or result._dof <= 0:
        return None
    params = np.asarray(result.params_opt, dtype=float)
    fixed = circuit.get_all_fixed_params()
    series_r = [offset for offset, _, node in members if type(node) is R and not fixed[offset]]
    # ponytail: exactly one free series R. Two are degenerate with each other
    # already (stderr inf); a free G in series is not handled, a fixed R or G
    # is a constant.
    if len(series_r) != 1 or any(type(node) is G and not fixed[offset]
                                 for offset, _, node in members):
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
        arc = _block_arc(node, params[offset:offset + n])
        if arc is None or arc[1] is None or any(fixed[offset:offset + n]):
            continue
        i_R, R_k, _ = arc
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
        # A non-finite refit would warn from the impedance; it is dropped below
        with np.errstate(all='ignore'):
            dchi2 = (rss(merged.circuit, merged.params_opt) - rss_fit) / s2
        if not np.isfinite(dchi2):
            logger.debug(f"Merging {labels[offset + i_R]} into {labels[rs]}: dchi2 {dchi2}")
            continue
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
    if split.dchi2 < -SPLIT_TIE_DCHI2:
        cost = (f"fits better (dchi2 = {split.dchi2:.3g}), so this fit is not the "
                f"least-squares minimum")
    elif split.dchi2 <= SPLIT_TIE_DCHI2:
        cost = f"fits as well (|dchi2| <= {SPLIT_TIE_DCHI2:g})"
    else:
        cost = f"is worse by only dchi2 = {split.dchi2:.3g} (< {chi2.ppf(SPLIT_WARN_LEVEL, 1):.2f})"
    widened = f"{name} and {arc}" if split.arc_index is not None else name
    result.diagnostics.warnings.append(
        f"{name} is determined only as {name} + {arc}: the fit without the arc of "
        f"{arc} {cost}. Either that arc acts as a plain resistance in series, or "
        f"another arc can take its place and the circuit has an arc too many. "
        f"The 95 % CI of {widened} is widened, "
        f"{name} to [{low[split.index]:.3g}, {high[split.index]:.4g}]")


__all__ = ['SeriesResistanceSplit', 'series_resistance_split', 'check_series_resistance']
