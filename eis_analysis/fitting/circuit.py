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
    result, Z_fit = fit_equivalent_circuit(freq, Z, circuit)

    # Plot (optional)
    from eis_analysis.visualization import plot_circuit_fit
    fig = plot_circuit_fit(freq, Z, result)
"""

import numpy as np
import logging
import warnings
from typing import Tuple, List, Optional, Sequence
from numpy.typing import NDArray
from scipy.optimize import OptimizeWarning
from dataclasses import dataclass, field

from .covariance import compute_covariance_matrix
from .bounds import (generate_simple_bounds, build_bound_status, log_scale_ci_mask,
                     validate_fixed_params)
from .diagnostics import compute_residual_weights, compute_fit_metrics, compute_significance
from .optimizer import least_squares_normalized
from .jacobian import make_jacobian_function
from .block_order import canonical_block_order
from .series_resistance import check_series_resistance
from .result import Circuit, FitDiagnostics, FitResult

logger = logging.getLogger(__name__)

# Valid weighting types
VALID_WEIGHTINGS = ['uniform', 'sqrt', 'proportional', 'modulus']


@dataclass
class OptimizationSetup:
    """Internal data structure for optimization setup."""
    initial_guess: List[float]
    lower_bounds: List[float]
    upper_bounds: List[float]
    fixed_params: List[bool]
    param_labels: List[str]
    param_labels_indexed: List[str]
    clipped_params: List[int] = field(default_factory=list)


def _prepare_optimization(
    circuit: Circuit,
    bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None
) -> OptimizationSetup:
    """
    Prepare optimization setup: extract parameters, labels, bounds.

    `bounds` overrides the PARAMETER_BOUNDS of generate_simple_bounds().
    Returns OptimizationSetup with all necessary configuration.
    """
    initial_guess = list(circuit.get_all_params())

    param_labels_raw = circuit.get_param_labels()
    # ponytail: O(n^2) in the parameter count, fine below ~1000 parameters;
    # a running Counter if a circuit ever gets that large.
    param_labels_indexed = [f"{label}{param_labels_raw[:i].count(label)}"
                            for i, label in enumerate(param_labels_raw)]

    fixed_params = circuit.get_all_fixed_params()
    if bounds is None:
        lower_bounds, upper_bounds = generate_simple_bounds(param_labels_raw)
    else:
        lower_bounds, upper_bounds = [float(b) for b in bounds[0]], [float(b) for b in bounds[1]]
        if not len(lower_bounds) == len(upper_bounds) == len(initial_guess):
            raise ValueError(
                f"bounds have {len(lower_bounds)} and {len(upper_bounds)} elements, "
                f"but circuit has {len(initial_guess)} parameters"
            )
        for label, lb, ub in zip(param_labels_indexed, lower_bounds, upper_bounds):
            if not lb < ub:  # also catches NaN
                raise ValueError(f"bounds of {label}: lower {lb} must be below upper {ub}")

    # Clip initial guess to bounds. Fixed parameters are exempt: they are
    # never optimized, so the bounds do not apply to them and clipping
    # would silently fit a value the caller did not ask for (DE keeps the
    # fixed value as given; validate_fixed_params() warns instead).
    clipped_params = []
    for i, (ig, lb, ub) in enumerate(zip(initial_guess, lower_bounds, upper_bounds)):
        if not fixed_params[i] and (ig < lb or ig > ub):
            initial_guess[i] = min(max(ig, lb), ub)
            clipped_params.append(i)

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
    Z: NDArray[np.complexfloating],
    circuit: Circuit,
    weighting: str = 'modulus',
    initial_guess: Optional[List[float]] = None,
    use_analytic_jacobian: bool = True,
    bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None
) -> Tuple[FitResult, NDArray[np.complexfloating]]:
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
    use_analytic_jacobian : bool, optional
        Use analytic Jacobian (default: True)
    bounds : (lower, upper), optional
        Per-parameter bounds in circuit order. Default: PARAMETER_BOUNDS by
        parameter type, which are absolute - pass bounds derived from the
        data where the result must not depend on the units.

    Returns
    -------
    result : FitResult
        Fitting results with all diagnostics
    Z_fit : ndarray of complex
        Predicted impedance from fit

    Notes
    -----
    After the fit each arc in series is merged into the series resistance
    and refitted (one local fit per arc, series_resistance.py); where that
    costs little, result.series_split widens the CI of the series R.
    """
    result, Z_fit = _fit_local(frequencies, Z, circuit, weighting, initial_guess,
                               use_analytic_jacobian, bounds)
    check_series_resistance(result, frequencies, Z, weighting, bounds, use_analytic_jacobian)
    return result, Z_fit


def _fit_local(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    circuit: Circuit,
    weighting: str = 'modulus',
    initial_guess: Optional[List[float]] = None,
    use_analytic_jacobian: bool = True,
    bounds: Optional[Tuple[Sequence[float], Sequence[float]]] = None
) -> Tuple[FitResult, NDArray[np.complexfloating]]:
    """fit_equivalent_circuit without the series resistance check, for the
    fits that check runs itself and for the starts of a multistart."""
    # Validate weighting parameter
    if weighting not in VALID_WEIGHTINGS:
        raise ValueError(f"weighting must be one of {VALID_WEIGHTINGS}, got '{weighting}'")

    # Step 1: Prepare optimization
    setup = _prepare_optimization(circuit, bounds)

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
        for i, (ig, cv, is_fixed) in enumerate(zip(initial_guess, circuit_values, setup.fixed_params)):
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

    # Raises when every parameter is fixed (empty optimization vector)
    fixed_warnings = validate_fixed_params(
        initial_guess_list, lower_bounds, upper_bounds, fixed_params, param_labels
    )

    # Step 2: Create residual function
    initial_guess_for_opt = initial_guess_list
    bounds_for_opt = lower_bounds, upper_bounds
    initial_guess_full = list(initial_guess_list)

    # Precompute weights
    weights = compute_residual_weights(Z, weighting)

    if any(fixed_params):
        initial_guess_for_opt = [v for v, f in zip(initial_guess_list, fixed_params) if not f]
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
    diag_warnings = list(fixed_warnings)

    # Report initial-guess values clipped into bounds by _prepare_optimization
    # (audit D3). With an explicit initial_guess override the clipped values
    # were discarded by the merge above, so there is nothing to report.
    if initial_guess is None:
        for i in setup.clipped_params:
            msg = (
                f"Initial guess for parameter {param_labels[i]} = {circuit_values[i]:.3e} "
                f"outside bounds, clipped to {initial_guess_list[i]:.3e}"
            )
            diag_warnings.append(msg)
            logger.warning(msg)

    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", OptimizeWarning)

            opt_result = least_squares_normalized(
                residual, initial_guess_for_opt, jac_func, bounds_for_opt,
                max_nfev=10000
            )

            params_opt_free = opt_result.x

            # Reconstruct full parameters
            if any(fixed_params):
                params_opt = np.array(reconstruct_params(params_opt_free))
            else:
                params_opt = params_opt_free

            # Interchangeable blocks in ascending tau, before anything is
            # derived from the parameters: the covariance, bound status and
            # warnings below then come out in that order on their own
            perm, free_perm, grouped = canonical_block_order(
                circuit, params_opt, fixed_params, lower_bounds, upper_bounds)
            params_opt = np.asarray(params_opt)[perm]
            opt_result.jac = opt_result.jac[:, free_perm]

            if not opt_result.success:
                diag_warnings.append(f"Optimizer did not converge: {opt_result.message}")

            # Check optimizer warnings
            for warning in w:
                if issubclass(warning.category, OptimizeWarning):
                    diag_warnings.append(f"Optimizer warning: {warning.message}")

        # Step 4: Check bounds proximity. Single criterion for the whole
        # result: bound_status (classify_bound_status) drives both the
        # params_at_bounds/bounds_warnings diagnostics and the per-parameter
        # tags in FitResult, so the two channels cannot disagree. Indices and
        # bounds are in full parameter space (fixed params included).
        bound_status = build_bound_status(
            params_opt, lower_bounds, upper_bounds, fixed_params
        )
        bounds_warnings = []
        params_at_bounds = []
        for i, status in enumerate(bound_status):
            if status not in ('lower', 'upper'):
                continue
            params_at_bounds.append(i)
            bound_val = lower_bounds[i] if status == 'lower' else upper_bounds[i]
            bounds_warnings.append(
                f"Parameter {param_labels[i]} = {params_opt[i]:.3e} near {status} "
                f"bound {bound_val:.1e}"
            )

        # Step 5: Compute covariance matrix
        cov_result = compute_covariance_matrix(
            jacobian=opt_result.jac,
            residuals=opt_result.fun,
            n_params=len(fixed_params),
            fixed_params=fixed_params
        )
        params_stderr = cov_result.stderr

        # Step 6: Parameter diagnostics (silent check)
        param_warnings = []
        for i, (p, s) in enumerate(zip(params_opt, params_stderr)):
            if np.isinf(s) or np.isnan(s):
                param_warnings.append(f"Parameter {param_labels[i]}: undefined uncertainty")
            elif s > abs(p) * 2:
                # p = 0 (a parameter on a zero bound) has no relative uncertainty;
                # p just above it overflows the ratio to inf, which is the message
                with np.errstate(over='ignore'):
                    rel = f"{s/abs(p)*100:.0f}%" if p != 0 else f"stderr {s:.3g} at 0"
                param_warnings.append(f"Parameter {param_labels[i]}: very high relative uncertainty ({rel})")

        # Step 7: Compute fit metrics
        Z_fit = circuit.impedance(frequencies, list(params_opt))
        fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, weighting)

        # Step 8: Update circuit with fitted parameters
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

        # Step 9: Create result object (bound_status computed in Step 4)
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
            n_free_params=len(initial_guess_for_opt),
            bound_status=bound_status,
            _dof=cov_result.dof,
            _ci_log_scale=log_scale_ci_mask(lower_bounds, upper_bounds),
            params_significance=compute_significance(circuit, frequencies, params_opt),
            block_order=perm.tolist() if grouped else None
        )

    except Exception as e:
        # Debug only: the caller gets the message with the exception, and the
        # series resistance check skips a merge that fails
        logger.debug(f"Fit failed: {type(e).__name__}: {e}")
        raise RuntimeError(f"Circuit fitting failed: {e}") from e

    return result, Z_fit


__all__ = [
    'fit_equivalent_circuit',
    'FitResult',
    'FitDiagnostics',
    'Circuit',
]
