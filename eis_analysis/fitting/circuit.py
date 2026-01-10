"""
Equivalent circuit fitting module for EIS analysis using operator overloading.

Clean design: No logging in core functions, all diagnostics returned as data.
CLI layer is responsible for user output.

Usage:
    from eis_analysis.fitting.circuit_elements import R, C, Q, W
    from eis_analysis.fitting import fit_equivalent_circuit

    # Build circuit using operators
    circuit = R(100) - (R(5000) | C(1e-6))

    # Fit to data
    result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
import warnings
from typing import Tuple, Union, List, Optional
from numpy.typing import NDArray
from scipy.optimize import least_squares, OptimizeWarning
from dataclasses import dataclass, field

from .circuit_elements import CircuitElement
from .circuit_builder import Series, Parallel
from .covariance import compute_covariance_matrix, compute_confidence_interval
from .bounds import generate_simple_bounds, check_bounds_proximity
from .diagnostics import compute_weights, check_parameter_diagnostics, compute_fit_metrics
from .jacobian import make_jacobian_function

logger = logging.getLogger(__name__)

# Type alias for circuit
Circuit = Union[CircuitElement, Series, Parallel]

# Valid weighting types
VALID_WEIGHTINGS = ['uniform', 'sqrt', 'proportional', 'modulus']


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

    # Bounds info
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
    fit_error_rel : float
        Relative fit error [%]
    fit_error_abs : float
        Absolute fit error [Ohm]
    quality : str
        Fit quality: 'excellent', 'good', 'acceptable', or 'poor'
    condition_number : float
        Condition number of Jacobian (high = ill-conditioned)
    is_well_conditioned : bool
        True if condition number < 1e10 (covariance reliable)
    cov : ndarray of float or None
        Covariance matrix of parameters (None if computation failed)
    diagnostics : FitDiagnostics or None
        Detailed diagnostics
    _n_data : int
        Number of data points (internal, for CI computation)
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
    _n_data: int = 0

    @property
    def params_ci_95(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """95% confidence intervals for parameters."""
        if np.any(np.isinf(self.params_stderr)) or np.any(np.isnan(self.params_stderr)):
            return (
                np.full_like(self.params_opt, -np.inf),
                np.full_like(self.params_opt, np.inf)
            )
        return compute_confidence_interval(
            self.params_opt, self.params_stderr, self._n_data, 0.95
        )

    @property
    def params_ci_99(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """99% confidence intervals for parameters."""
        if np.any(np.isinf(self.params_stderr)) or np.any(np.isnan(self.params_stderr)):
            return (
                np.full_like(self.params_opt, -np.inf),
                np.full_like(self.params_opt, np.inf)
            )
        return compute_confidence_interval(
            self.params_opt, self.params_stderr, self._n_data, 0.99
        )

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


@dataclass
class OptimizationSetup:
    """Internal data structure for optimization setup."""
    initial_guess: List[float]
    lower_bounds: Optional[List[float]]
    upper_bounds: Optional[List[float]]
    fixed_params: Optional[List[bool]]
    param_labels: Optional[List[str]]
    param_labels_indexed: Optional[List[str]]
    clipped_params: List[int] = field(default_factory=list)


def _prepare_optimization(circuit: Circuit, weighting: str) -> OptimizationSetup:
    """
    Prepare optimization setup: extract parameters, labels, bounds.

    Returns OptimizationSetup with all necessary configuration.
    """
    initial_guess = list(circuit.get_all_params())

    # Get parameter labels
    param_labels_raw = None
    param_labels_indexed = None
    if hasattr(circuit, 'get_param_labels'):
        param_labels_raw = circuit.get_param_labels()
        label_counts = {}
        param_labels_indexed = []
        for label in param_labels_raw:
            if label in label_counts:
                label_counts[label] += 1
            else:
                label_counts[label] = 0
            param_labels_indexed.append(f"{label}{label_counts[label]}")

    # Get fixed parameters
    fixed_params = None
    if hasattr(circuit, 'get_all_fixed_params'):
        fixed_params = circuit.get_all_fixed_params()

    # Generate bounds
    lower_bounds, upper_bounds = None, None
    clipped_params = []
    if param_labels_raw is not None:
        lower_bounds, upper_bounds = generate_simple_bounds(param_labels_raw)

        # Clip initial guess to bounds
        initial_guess_clipped = []
        for i, (ig, lb, ub) in enumerate(zip(initial_guess, lower_bounds, upper_bounds)):
            if ig < lb:
                initial_guess_clipped.append(lb)
                clipped_params.append(i)
            elif ig > ub:
                initial_guess_clipped.append(ub)
                clipped_params.append(i)
            else:
                initial_guess_clipped.append(ig)

        if clipped_params:
            initial_guess = initial_guess_clipped

    return OptimizationSetup(
        initial_guess=initial_guess,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        fixed_params=fixed_params,
        param_labels=param_labels_raw,
        param_labels_indexed=param_labels_indexed,
        clipped_params=clipped_params
    )


def fit_equivalent_circuit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    circuit: Circuit,
    weighting: str = 'modulus',
    initial_guess: Optional[List[float]] = None,
    plot: bool = True,
    verbose: bool = True,  # Kept for backward compatibility, ignored
    use_analytic_jacobian: bool = True
) -> Tuple[FitResult, NDArray[np.complex128], Optional[plt.Figure]]:
    """
    Fit equivalent circuit to impedance data.

    Parameters
    ----------
    frequencies : ndarray of float
        Measurement frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    circuit : Circuit
        Circuit built using operator overloading
    weighting : str, optional
        Weighting type: 'uniform', 'sqrt', 'modulus' (default), 'proportional'
    initial_guess : list of float, optional
        Override initial guess for parameters
    plot : bool, optional
        Create visualization plot (default: True)
    verbose : bool, optional
        Ignored (kept for backward compatibility)
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian (default: True)

    Returns
    -------
    result : FitResult
        Fitting results with all diagnostics
    Z_fit : ndarray of complex
        Predicted impedance from fit
    fig : matplotlib.figure.Figure or None
        Nyquist plot (None if plot=False)
    """
    # Validate weighting parameter
    if weighting not in VALID_WEIGHTINGS:
        raise ValueError(f"weighting must be one of {VALID_WEIGHTINGS}, got '{weighting}'")

    # Step 1: Prepare optimization
    setup = _prepare_optimization(circuit, weighting)

    # Save original circuit values for fixed parameters
    circuit_values = list(circuit.get_all_params())

    # Override initial guess if provided
    if initial_guess is not None:
        if len(initial_guess) != len(setup.initial_guess):
            raise ValueError(
                f"initial_guess has {len(initial_guess)} elements, "
                f"but circuit has {len(setup.initial_guess)} parameters"
            )
        # Merge: use circuit values for fixed params, initial_guess for free params
        merged_guess = []
        fixed = setup.fixed_params or [False] * len(initial_guess)
        for i, (ig, cv, is_fixed) in enumerate(zip(initial_guess, circuit_values, fixed)):
            merged_guess.append(cv if is_fixed else ig)
        setup = OptimizationSetup(
            initial_guess=merged_guess,
            lower_bounds=setup.lower_bounds,
            upper_bounds=setup.upper_bounds,
            fixed_params=setup.fixed_params,
            param_labels=setup.param_labels,
            param_labels_indexed=setup.param_labels_indexed,
            clipped_params=setup.clipped_params
        )

    initial_guess_list = setup.initial_guess
    lower_bounds = setup.lower_bounds
    upper_bounds = setup.upper_bounds
    fixed_params = setup.fixed_params
    param_labels = setup.param_labels_indexed

    # Step 2: Create residual function
    initial_guess_for_opt = initial_guess_list
    bounds_for_opt = lower_bounds, upper_bounds
    initial_guess_full = list(initial_guess_list)

    # Precompute weights
    weights = compute_weights(Z, weighting)

    if fixed_params is not None and any(fixed_params):
        initial_guess_for_opt = [v for v, f in zip(initial_guess_list, fixed_params) if not f]
        if lower_bounds is not None:
            bounds_for_opt = (
                [lb for lb, f in zip(lower_bounds, fixed_params) if not f],
                [ub for ub, f in zip(upper_bounds, fixed_params) if not f]
            )

        def reconstruct_params(free_params):
            full, idx = [], 0
            for is_fixed in fixed_params:
                full.append(initial_guess_full[len(full)] if is_fixed else free_params[idx])
                if not is_fixed:
                    idx += 1
            return full

        def residual(free_params):
            full_params = reconstruct_params(free_params)
            Z_pred = circuit.impedance(frequencies, full_params)
            return np.concatenate([
                (Z.real - Z_pred.real) * weights,
                (Z.imag - Z_pred.imag) * weights
            ])
    else:
        def residual(params):
            Z_pred = circuit.impedance(frequencies, list(params))
            return np.concatenate([
                (Z.real - Z_pred.real) * weights,
                (Z.imag - Z_pred.imag) * weights
            ])

    # Create Jacobian function
    jacobian_type = 'analytic'
    if use_analytic_jacobian:
        try:
            jac_func = make_jacobian_function(
                circuit, frequencies, weights,
                fixed_params=fixed_params,
                full_initial_guess=initial_guess_full
            )
        except NotImplementedError:
            jac_func = '2-point'
            jacobian_type = 'numeric'
    else:
        jac_func = '2-point'
        jacobian_type = 'numeric'

    # Step 3: Run optimization
    diag_warnings = []

    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", OptimizeWarning)

            bounds_arg = bounds_for_opt if bounds_for_opt[0] is not None else (-np.inf, np.inf)
            x_scale = np.maximum(np.abs(initial_guess_for_opt), 1e-10)

            opt_result = least_squares(
                residual,
                x0=initial_guess_for_opt,
                jac=jac_func,
                bounds=bounds_arg,
                max_nfev=10000,
                x_scale=x_scale
            )

            params_opt_free = opt_result.x

            # Reconstruct full parameters
            if fixed_params is not None and any(fixed_params):
                params_opt = np.array(reconstruct_params(params_opt_free))
            else:
                params_opt = params_opt_free

            if not opt_result.success:
                diag_warnings.append(f"Optimizer did not converge: {opt_result.message}")

            # Check optimizer warnings
            for warning in w:
                if issubclass(warning.category, OptimizeWarning):
                    diag_warnings.append(f"Optimizer warning: {warning.message}")

        # Step 4: Check bounds proximity
        bounds_warnings = []
        params_at_bounds = []
        if bounds_for_opt[0] is not None:
            # Check which parameters are at bounds (silent check)
            for i, (p, lb, ub) in enumerate(zip(params_opt_free, bounds_for_opt[0], bounds_for_opt[1])):
                if abs(p - lb) / max(abs(lb), 1e-10) < 0.01:
                    params_at_bounds.append(i)
                    bounds_warnings.append(f"Parameter {i} at lower bound")
                elif abs(p - ub) / max(abs(ub), 1e-10) < 0.01:
                    params_at_bounds.append(i)
                    bounds_warnings.append(f"Parameter {i} at upper bound")

        # Step 5: Compute covariance matrix
        n_params_total = len(fixed_params) if fixed_params is not None else len(params_opt)
        cov_result = compute_covariance_matrix(
            jacobian=opt_result.jac,
            residuals=opt_result.fun,
            n_params=n_params_total,
            fixed_params=fixed_params
        )
        params_stderr = cov_result.stderr

        # Step 6: Parameter diagnostics (silent check)
        param_warnings = []
        for i, (p, s) in enumerate(zip(params_opt, params_stderr)):
            if np.isinf(s) or np.isnan(s):
                param_warnings.append(f"Parameter {i}: undefined uncertainty")
            elif s > abs(p) * 2:
                param_warnings.append(f"Parameter {i}: very high relative uncertainty ({s/abs(p)*100:.0f}%)")

        # Step 7: Compute fit metrics
        Z_fit = circuit.impedance(frequencies, list(params_opt))
        fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

        # Step 8: Update circuit with fitted parameters
        if hasattr(circuit, 'update_params'):
            circuit.update_params(list(params_opt))

        # Build diagnostics
        diagnostics = FitDiagnostics(
            optimizer_status=opt_result.status,
            optimizer_message=opt_result.message,
            optimizer_success=opt_result.success,
            n_function_evals=opt_result.nfev,
            jacobian_type=jacobian_type,
            condition_number=cov_result.condition_number,
            covariance_rank=cov_result.rank,
            covariance_warning=cov_result.warning_message,
            params_at_bounds=params_at_bounds,
            bounds_warnings=bounds_warnings,
            param_warnings=param_warnings,
            warnings=diag_warnings
        )

        # Step 9: Create result object
        n_data = len(frequencies)
        result = FitResult(
            circuit=circuit,
            params_opt=params_opt,
            params_stderr=params_stderr,
            fit_error_rel=fit_error_rel,
            fit_error_abs=fit_error_abs,
            quality=quality,
            condition_number=cov_result.condition_number,
            is_well_conditioned=cov_result.is_well_conditioned,
            cov=cov_result.cov,
            diagnostics=diagnostics,
            param_labels=param_labels,
            _n_data=n_data
        )

    except Exception as e:
        logger.error(f"Fit failed: {type(e).__name__}: {e}")
        raise RuntimeError(f"Circuit fitting failed: {e}") from e

    # Step 10: Create visualization (only if requested)
    fig = None
    if plot:
        from ..visualization.plots import plot_circuit_fit

        f_min, f_max = frequencies.min(), frequencies.max()
        freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
        Z_fit_plot = circuit.impedance(freq_plot, list(params_opt))
        fig = plot_circuit_fit(frequencies, Z, Z_fit_plot, circuit, Z_fit_at_data=Z_fit)

    return result, Z_fit, fig


__all__ = [
    'fit_equivalent_circuit',
    'FitResult',
    'FitDiagnostics',
    'Circuit',
]
