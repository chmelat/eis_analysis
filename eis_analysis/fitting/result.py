"""
Result of a circuit fit: FitResult and its FitDiagnostics.

The confidence intervals are computed here from the stderr, widened where
the series resistance check found a flat branch (series_resistance.py).
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from .circuit_builder import Parallel, Series
from .circuit_elements import CircuitElement
from .covariance import compute_confidence_interval
from .series_resistance import SeriesResistanceSplit

# Type alias for circuit
Circuit = Union[CircuitElement, Series, Parallel]


@dataclass
class FitDiagnostics:
    """Diagnostics from circuit fitting."""
    # Optimization info
    optimizer_status: int
    optimizer_message: str
    optimizer_success: bool
    n_function_evals: int
    jacobian_type: str  # 'analytic' or 'numeric'

    # Covariance info
    condition_number: float
    covariance_rank: int
    covariance_warning: Optional[str] = None

    # Bounds info. Full-parameter-space indices (fixed params included),
    # derived from FitResult.bound_status (classify_bound_status criterion) —
    # the same criterion that suppresses CIs for at-bound parameters.
    params_at_bounds: List[int] = field(default_factory=list)
    bounds_warnings: List[str] = field(default_factory=list)

    # Parameter diagnostics
    param_warnings: List[str] = field(default_factory=list)

    # General warnings
    warnings: List[str] = field(default_factory=list)


@dataclass
class FitResult:
    """
    Result from circuit fitting.

    Attributes
    ----------
    circuit : Circuit
        Original circuit object
    params_opt : ndarray of float
        Optimized parameters
    params_stderr : ndarray of float
        Standard errors of parameters (from covariance matrix)
    params_significance : ndarray of float or None
        Per-parameter sensitivity max|dln|Z|/dlnP|, aligned with params_opt.
        Answers a different question than params_stderr: not how precisely a
        parameter is determined, but whether it matters in this frequency
        window at all. ~1 means it dominates somewhere, below
        SIGNIFICANCE_NEGLIGIBLE the element may be dropped. None if the
        circuit has an element with no analytic derivative.
    fit_error_rel : float
        Relative fit error [%]
    fit_error_abs : float
        Absolute fit error [Ohm]
    quality : str
        Fit quality: 'excellent', 'good', 'acceptable', or 'poor'
    condition_number : float
        Condition number of the column-scaled normal matrix
        J_s^T J_s = cond(J_s)^2, where J_s has unit-norm columns
        (high = ill-conditioned). Scale-invariant: reflects genuine
        parameter correlation, not unit disparity between parameters.
    is_well_conditioned : bool
        True if cond(J_s^T J_s) < 1e10 (covariance reliable)
    cov : ndarray of float or None
        Covariance matrix of parameters (None if computation failed)
    diagnostics : FitDiagnostics or None
        Detailed diagnostics
    block_order : list of int or None
        Permutation that put interchangeable series blocks in ascending tau
        (canonical_block_order): params_opt is the optimizer's vector indexed
        by it, the identity if they already were. None when the circuit has
        no such blocks.
    series_split : SeriesResistanceSplit or None
        An arc the data cannot tell apart from part of the series resistance
        (series_resistance.py); params_ci_95 and params_ci_99 widen the CIs of
        the series resistance and the arc's R by it. None when no merge is that cheap.
    n_free_params : int
        Number of freely optimized parameters, with fixed parameters
        excluded. Needed by compute_information_criteria; kept explicit
        rather than derived from _dof so callers outside this module do not
        have to assume the residual count is 2*len(frequencies).
    _dof : int
        Residual degrees of freedom (internal, for CI computation).
        Matches CovarianceResult.dof = n_residuals - n_free_params.
    _ci_log_scale : list of bool or None
        Per-parameter mask for log-space CI (internal); True for positive
        scale parameters. None = no bound info, all CIs linear/symmetric.
    """
    circuit: Circuit
    params_opt: NDArray[np.float64]
    params_stderr: NDArray[np.float64]
    fit_error_rel: float
    fit_error_abs: float = 0.0
    quality: str = "unknown"
    condition_number: float = 1.0
    is_well_conditioned: bool = True
    cov: Optional[NDArray[np.float64]] = None
    diagnostics: Optional[FitDiagnostics] = None
    param_labels: Optional[List[str]] = None
    n_free_params: int = 0
    # Per-parameter bound status: '' (interior), 'lower', 'upper', or 'fixed'.
    # Aligned with params_opt; None means caller did not supply bound info.
    bound_status: Optional[List[str]] = None
    _dof: int = 0
    # Per-parameter mask: True for positive scale parameters whose CI is
    # computed in log space (see compute_confidence_interval). None means
    # no bound info was available; all CIs stay linear/symmetric.
    _ci_log_scale: Optional[List[bool]] = None
    # Per-parameter sensitivity max|dln|Z|/dlnP|, aligned with params_opt.
    # None when the circuit has an element with no analytic derivative.
    # See compute_significance in diagnostics.py.
    params_significance: Optional[NDArray[np.float64]] = None
    block_order: Optional[List[int]] = None
    series_split: Optional[SeriesResistanceSplit] = None

    def _ci(self, confidence_level: float) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        # Non-finite stderr yields (-inf, +inf) per parameter inside
        # compute_confidence_interval; other parameters keep finite CIs.
        low, high = compute_confidence_interval(
            self.params_opt, self.params_stderr, self._dof, confidence_level,
            log_scale=self._ci_log_scale
        )
        # The profile of the series R is flat from its bound to the merged
        # value, which the covariance at the fit cannot see
        split = self.series_split
        if split is not None and split.widens(confidence_level):
            top = split.merged._ci(confidence_level)[1][split.merged_index]
            # On that branch the arc's R takes the rest of the sum
            for i, bound in ((split.index, split.lower), (split.arc_index, split.arc_lower)):
                if i is not None:
                    low[i] = bound
                    high[i] = max(high[i], top)
        return low, high

    @property
    def params_ci_95(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """95% confidence intervals for parameters."""
        return self._ci(0.95)

    @property
    def params_ci_99(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """99% confidence intervals for parameters."""
        return self._ci(0.99)

    @property
    def all_warnings(self) -> List[str]:
        """Collect all warnings from diagnostics."""
        if self.diagnostics is None:
            return []
        warnings = []
        warnings.extend(self.diagnostics.warnings)
        warnings.extend(self.diagnostics.bounds_warnings)
        warnings.extend(self.diagnostics.param_warnings)
        if self.diagnostics.covariance_warning:
            warnings.append(self.diagnostics.covariance_warning)
        return warnings

    def __repr__(self) -> str:
        lines = []
        lines.append("Fit Result:")
        lines.append(f"  Circuit: {self.circuit}")
        lines.append(f"  Parameters: {self.params_opt}")
        lines.append(f"  Std errors: {self.params_stderr}")

        ci_low, ci_high = self.params_ci_95
        if not np.all(np.isinf(ci_low)):
            lines.append(f"  95% CI low: {ci_low}")
            lines.append(f"  95% CI high: {ci_high}")

        lines.append(f"  Fit error: {self.fit_error_rel:.2f}% (rel), {self.fit_error_abs:.2f} Ohm (abs)")
        lines.append(f"  Quality: {self.quality}")
        return '\n'.join(lines)
