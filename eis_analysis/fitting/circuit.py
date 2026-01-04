"""
Equivalent circuit fitting module for EIS analysis using operator overloading.

This module provides the main fitting function for equivalent circuits.
Circuit building uses operator overloading inspired by EISAnalysis.jl (Julia).

Usage:
    from eis_analysis.fitting.circuit_elements import R, C, Q, W
    from eis_analysis.fitting import fit_equivalent_circuit

    # Build circuit using operators
    circuit = R(100) - (R(5000) | C(1e-6))

    # Fit to data
    result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

Author: EIS Analysis Toolkit
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
import warnings
from typing import Tuple, Union, List, Optional
from numpy.typing import NDArray
from scipy.optimize import least_squares, OptimizeWarning
from dataclasses import dataclass

from .circuit_elements import CircuitElement
from .circuit_builder import Series, Parallel
from .covariance import compute_covariance_matrix, compute_confidence_interval
from .bounds import generate_simple_bounds, check_bounds_proximity
from .diagnostics import compute_weights, check_parameter_diagnostics, compute_fit_metrics, log_fit_results
from .jacobian import make_jacobian_function

logger = logging.getLogger(__name__)

# Type alias for circuit
Circuit = Union[CircuitElement, Series, Parallel]

# Valid weighting types
VALID_WEIGHTINGS = ['uniform', 'sqrt', 'proportional', 'square']


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
    _n_data : int
        Number of data points (internal, for CI computation)
    """
    circuit: Circuit
    params_opt: NDArray[np.float64]
    params_stderr: NDArray[np.float64]
    fit_error_rel: float
    fit_error_abs: float
    quality: str
    condition_number: float = 1.0
    is_well_conditioned: bool = True
    cov: Optional[NDArray[np.float64]] = None
    _n_data: int = 0

    @property
    def params_ci_95(self) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """
        95% confidence intervals for parameters.

        Returns
        -------
        ci_low, ci_high : tuple of ndarray
            Lower and upper bounds of 95% CI
        """
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
        """
        99% confidence intervals for parameters.

        Returns
        -------
        ci_low, ci_high : tuple of ndarray
            Lower and upper bounds of 99% CI
        """
        if np.any(np.isinf(self.params_stderr)) or np.any(np.isnan(self.params_stderr)):
            return (
                np.full_like(self.params_opt, -np.inf),
                np.full_like(self.params_opt, np.inf)
            )
        return compute_confidence_interval(
            self.params_opt, self.params_stderr, self._n_data, 0.99
        )

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


def _prepare_optimization(circuit: Circuit, weighting: str, verbose: bool = True) -> OptimizationSetup:
    """
    Prepare optimization setup: extract parameters, labels, bounds.

    Parameters
    ----------
    circuit : Circuit
        Circuit object with initial guess values
    weighting : str
        Weighting type for logging
    verbose : bool
        Whether to log detailed progress

    Returns
    -------
    setup : OptimizationSetup
        Prepared optimization configuration
    """
    initial_guess = circuit.get_all_params()

    if verbose:
        logger.info("=" * 50)
        logger.info("Equivalent circuit fitting")
        logger.info("=" * 50)
        logger.info(f"Circuit: {circuit}")

        weighting_labels = {
            'uniform': 'uniform (w = 1)',
            'sqrt': 'square root (w = 1/sqrt|Z|)',
            'proportional': 'proportional (w = 1/|Z|)',
            'square': 'square (w = |Z|^2)'
        }
        logger.info(f"Weighting: {weighting_labels[weighting]}")

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
        if verbose:
            logger.info(f"Parameters: {param_labels_indexed}")

    logger.debug(f"Initial guess: {initial_guess}")
    logger.debug(f"Number of parameters: {len(initial_guess)}")

    # Get fixed parameters
    fixed_params = None
    if hasattr(circuit, 'get_all_fixed_params'):
        fixed_params = circuit.get_all_fixed_params()
        n_fixed = sum(fixed_params)
        if n_fixed > 0 and verbose:
            logger.info(f"Fixed parameters: {n_fixed} of {len(fixed_params)}")
            for i, (is_fixed, val) in enumerate(zip(fixed_params, initial_guess)):
                if is_fixed:
                    label = param_labels_indexed[i] if param_labels_indexed else f"p{i}"
                    logger.info(f"  {label} = {val:.4g} (fixed)")

    # Generate bounds
    lower_bounds, upper_bounds = None, None
    if param_labels_raw is not None:
        lower_bounds, upper_bounds = generate_simple_bounds(param_labels_raw)
        if verbose:
            logger.info("Bounds: physically constrained (R: 0.1mOhm-10GOhm, C: 1fF-100mF, n: 0.3-1.0)")

        # Clip initial guess to bounds
        initial_guess_clipped = []
        clipped_any = False
        for i, (ig, lb, ub) in enumerate(zip(initial_guess, lower_bounds, upper_bounds)):
            if ig < lb:
                initial_guess_clipped.append(lb)
                clipped_any = True
                logger.warning(f"  Parameter {i} ({param_labels_raw[i]}): initial guess {ig:.2e} < lower bound {lb:.2e}")
            elif ig > ub:
                initial_guess_clipped.append(ub)
                clipped_any = True
                logger.warning(f"  Parameter {i} ({param_labels_raw[i]}): initial guess {ig:.2e} > upper bound {ub:.2e}")
            else:
                initial_guess_clipped.append(ig)

        if clipped_any:
            logger.warning("Initial guess was adjusted to fit within bounds")
            initial_guess = initial_guess_clipped

    return OptimizationSetup(
        initial_guess=initial_guess,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        fixed_params=fixed_params,
        param_labels=param_labels_raw,
        param_labels_indexed=param_labels_indexed
    )


def fit_equivalent_circuit(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    circuit: Circuit,
    weighting: str = 'uniform',
    initial_guess: Optional[List[float]] = None,
    plot: bool = True,
    verbose: bool = True,
    use_analytic_jacobian: bool = True
) -> Tuple[FitResult, NDArray[np.complex128], Optional[plt.Figure]]:
    """
    Fit equivalent circuit to impedance data.

    Parameters
    ----------
    frequencies : ndarray of float
        Measurement frequencies [Hz] (N points)
    Z : ndarray of complex
        Complex impedance [Ohm] (N points)
    circuit : Circuit
        Circuit built using operator overloading.
        Values in circuit serve as INITIAL GUESS for fitting.
    weighting : str, optional
        Weighting type for optimization. Options:
        - 'uniform' (default): all points equally important (w = 1)
        - 'sqrt': compromise weighting (w = 1/sqrt|Z|)
        - 'proportional': inverse weighting (w = 1/|Z|)
        Weighting affects which frequencies have more influence on fit.
        Uniform favors low frequencies (large |Z|),
        proportional equalizes importance across all frequencies.
    initial_guess : list of float, optional
        Override initial guess for parameters. If None (default), uses
        values from circuit object. Useful for multi-start optimization.
    plot : bool, optional
        Whether to create visualization plot (default: True).
        Set to False for multi-start optimization to avoid creating
        many unnecessary figures.
    verbose : bool, optional
        Whether to log detailed progress and results (default: True).
        Set to False for multi-start optimization to reduce output.
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian for least_squares optimization (default: True).
        More accurate and faster than numeric approximation.

    Returns
    -------
    result : FitResult
        Fitting results (optimal parameters, errors, quality)
    Z_fit : ndarray of complex
        Predicted impedance from fit [Ohm]
    fig : matplotlib.figure.Figure or None
        Nyquist plot with data and fit (None if plot=False)

    Examples
    --------
    >>> from eis_analysis.fitting.circuit_elements import R, C
    >>> # Voigt element
    >>> circuit = R(100) - (R(5000) | C(1e-6))
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)
    >>> print(result.params_opt)
    [98.5, 4823.2, 8.7e-7]

    >>> # Randles circuit
    >>> from eis_analysis.fitting.circuit_elements import Q, W
    >>> circuit = R(10) - (R(100) - W(50)) | Q(1e-4, 0.8)
    >>> result, Z_fit, fig = fit_equivalent_circuit(freq, Z, circuit)

    Notes
    -----
    Operators:
    - `-` : series connection (Z = Z1 + Z2)
    - `|` : parallel connection (1/Z = 1/Z1 + 1/Z2)
    - `*` : parameter scaling (2*R(100) = R(200))
    - `**`: exponent for Q/CPE (Q()**0.9)

    See Also
    --------
    circuit_elements : Definition of all circuit elements
    circuit_builder : Series and Parallel combinators
    """
    # Validate weighting parameter
    if weighting not in VALID_WEIGHTINGS:
        raise ValueError(f"weighting must be one of {VALID_WEIGHTINGS}, got '{weighting}'")

    # Step 1: Prepare optimization
    setup = _prepare_optimization(circuit, weighting, verbose)

    # Save original circuit values for fixed parameters
    circuit_values = list(circuit.get_all_params())

    # Override initial guess if provided (only for FREE parameters)
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
            param_labels_indexed=setup.param_labels_indexed
        )
    initial_guess = setup.initial_guess
    lower_bounds = setup.lower_bounds
    upper_bounds = setup.upper_bounds
    fixed_params = setup.fixed_params
    param_labels = setup.param_labels_indexed

    # Step 2: Create residual function
    initial_guess_for_opt = initial_guess
    bounds_for_opt = lower_bounds, upper_bounds
    initial_guess_full = list(initial_guess)

    # Precompute weights (used by both residual and Jacobian)
    weights = compute_weights(Z, weighting)

    if fixed_params is not None and any(fixed_params):
        initial_guess_for_opt = [v for v, f in zip(initial_guess, fixed_params) if not f]
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

    # Create Jacobian function if using analytic
    if use_analytic_jacobian:
        try:
            jac_func = make_jacobian_function(
                circuit, frequencies, weights,
                fixed_params=fixed_params,
                full_initial_guess=initial_guess_full
            )
        except NotImplementedError as e:
            if verbose:
                logger.warning(f"Analytic Jacobian failed: {e}, using numeric")
            jac_func = '2-point'
    else:
        jac_func = '2-point'

    # Step 3: Run optimization
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

            # Log optimizer status
            if verbose:
                logger.info(f"Optimizer status: {opt_result.status} - {opt_result.message}")
            if not opt_result.success:
                logger.warning("WARNING: Optimizer did not converge!")

            # Step 4: Check bounds proximity
            if bounds_for_opt[0] is not None:
                check_bounds_proximity(
                    params_opt, params_opt_free,
                    bounds_for_opt[0], bounds_for_opt[1],
                    fixed_params
                )

            # Step 5: Compute covariance matrix
            n_params_total = len(fixed_params) if fixed_params is not None else len(params_opt)
            cov_result = compute_covariance_matrix(
                jacobian=opt_result.jac,
                residuals=opt_result.fun,
                n_params=n_params_total,
                fixed_params=fixed_params
            )
            params_stderr = cov_result.stderr

            if verbose:
                logger.info(f"Covariance: cond={cov_result.condition_number:.2e}, rank={cov_result.rank}/{len(params_opt_free)}")
            if cov_result.warning_message:
                logger.warning(f"WARNING: {cov_result.warning_message}")

            # Check optimizer warnings
            for warning in w:
                if issubclass(warning.category, OptimizeWarning):
                    logger.warning(f"Optimizer warning: {warning.message}")

        # Step 6: Run diagnostics (warnings are always shown)
        if verbose:
            check_parameter_diagnostics(params_opt, params_stderr, cov_result.cov, param_labels)

        # Step 7: Compute fit metrics
        Z_fit = circuit.impedance(frequencies, list(params_opt))
        fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

        # Step 8: Update circuit with fitted parameters
        # This ensures circuit.params contains fitted values (not initial guesses)
        if hasattr(circuit, 'update_params'):
            circuit.update_params(list(params_opt))

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
            _n_data=n_data
        )

        # Step 10: Log results
        if verbose:
            ci_low, ci_high = result.params_ci_95
            log_fit_results(
                params_opt, params_stderr, ci_low, ci_high,
                fit_error_rel, fit_error_abs, quality, param_labels
            )

    except Exception as e:
        logger.error(f"Fit failed: {type(e).__name__}: {e}")
        raise RuntimeError(f"Circuit fitting failed: {e}") from e

    # Step 11: Create visualization (only if requested)
    fig = None
    if plot:
        from ..visualization.plots import plot_circuit_fit

        f_min, f_max = frequencies.min(), frequencies.max()
        freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
        Z_fit_plot = circuit.impedance(freq_plot, list(params_opt))
        # Pass Z_fit (at original frequencies) for residuals calculation
        fig = plot_circuit_fit(frequencies, Z, Z_fit_plot, circuit, Z_fit_at_data=Z_fit)

    return result, Z_fit, fig


__all__ = [
    'fit_equivalent_circuit',
    'FitResult',
    'Circuit',
]
