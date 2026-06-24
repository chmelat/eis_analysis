"""
Equivalent circuit fitting handlers for the EIS CLI.

- run_circuit_fitting: dispatch to Voigt-chain or standard circuit fitting
- _fit_voigt_chain / _fit_standard_circuit: fitting backends
- _log_* helpers: optimizer and fit-result diagnostics logging
"""

import argparse
import logging
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import EISAnalysisError, save_figure, parse_circuit_expression
from ...fitting import (
    fit_equivalent_circuit,
    fit_circuit_multistart,
    fit_circuit_diffevo,
    fit_voigt_chain_linear,
    FitResult,
    MultistartResult,
    DiffEvoResult,
)
from ...fitting.diagnostics import compute_fit_metrics

logger = logging.getLogger(__name__)


# =============================================================================
# Fitting Diagnostics Logging
# =============================================================================

def _log_diffevo_diagnostics(diffevo_result: DiffEvoResult) -> None:
    """Log Differential Evolution optimization diagnostics."""
    diag = diffevo_result.diagnostics
    if diag is None:
        return

    log_separator()
    logger.info("Differential Evolution optimization")
    log_separator()

    # Settings - map strategy name to option number
    strategy_to_option = {'randtobest1bin': 1, 'best1bin': 2, 'rand1bin': 3}
    option_num = strategy_to_option.get(diag.strategy, 1)
    logger.info(f"  Strategy: {diag.strategy} (option {option_num})")
    logger.info(f"  Population: {diag.popsize} * n_params")
    logger.info(f"  Max iterations: {diag.maxiter}")
    logger.info(f"  Tolerance: {diag.tol}")
    logger.info(f"  Workers: {diag.workers}")
    logger.info(f"  Weighting: {diag.weighting} (w=1/|Z|)")

    # Get params info from best_result
    result = diffevo_result.best_result
    n_free = len(result.params_opt) - diag.n_fixed_params
    logger.info(f"  Parameters: {n_free} (free)")

    # Initial guess captured before DE ran (circuit.update_params overwrote
    # circuit.get_all_params() with the final fit, so we read the snapshot
    # stored in diagnostics instead).
    logger.info(f"  Initial guess: {np.array(diag.initial_guess)}")
    logger.info("")

    # DE progress
    logger.info("Running differential evolution...")
    logger.info(f"  DE converged: {diag.de_converged}")
    logger.info(f"  DE iterations: {diag.de_iterations}")
    logger.info(f"  DE evaluations: {diag.de_evaluations}")
    logger.info(f"  DE error: {diag.de_error:.3f}%")
    logger.info("")

    # Refinement
    jac_type = "analytic" if diag.jacobian_type == "analytic" else "numeric"
    logger.info(f"Refining with least_squares ({jac_type} Jacobian)...")
    logger.info(f"  Refined error: {diag.refined_error:.3f}%")
    improvement = (diag.de_error - diag.refined_error) / diag.de_error * 100 if diag.de_error > 0 else 0
    logger.info(f"  Improvement: {improvement:+.1f}%")
    logger.info("")

    # Summary
    log_separator()
    logger.info("Differential Evolution results")
    log_separator()
    logger.info(f"  Strategy: {diag.strategy}")
    logger.info(f"  Total evaluations: {diag.total_evaluations}")
    logger.info(f"  DE error: {diag.de_error:.3f}% -> Refined: {diag.refined_error:.3f}%")
    logger.info("")


def _log_multistart_diagnostics(multistart_result: MultistartResult) -> None:
    """Log Multi-start optimization diagnostics."""
    diag = multistart_result.diagnostics
    if diag is None:
        return

    log_separator()
    logger.info("Multi-start optimization")
    log_separator()

    # Settings
    logger.info(f"  Restarts: {diag.n_restarts}")
    logger.info(f"  Perturbation scale: {diag.scale} sigma")
    logger.info(f"  Perturbation method: {diag.perturbation_method}")
    logger.info(f"  Weighting: {diag.weighting}")
    logger.info(f"  Jacobian: {diag.jacobian_type}")
    logger.info(f"  Parallel: {diag.parallel}")
    logger.info("")

    # Results
    logger.info(f"Running {diag.n_restarts} optimization starts...")
    logger.info(f"  Successful: {diag.n_successful}/{diag.n_restarts}")
    logger.info(f"  Initial error: {diag.initial_error:.3f}%")
    logger.info(f"  Best error: {diag.best_error:.3f}%")
    logger.info(f"  Best start: #{diag.best_start_index}")
    logger.info(f"  Improvement: {multistart_result.improvement:+.1f}%")
    logger.info("")

    # Warnings
    for warning in diag.warnings:
        logger.warning(f"  {warning}")


def _log_fit_result(result: FitResult) -> None:
    """Log fit result with parameters and confidence intervals."""
    logger.info("")
    logger.info("Fit results:")
    logger.info("  Parameters:")

    # Get param labels
    labels = result.param_labels
    if labels is None:
        labels = [f"p{i}" for i in range(len(result.params_opt))]

    # Get 95% confidence intervals
    ci_low, ci_high = result.params_ci_95

    # Print each parameter. If a parameter sits at a bound or is fixed, the
    # Jacobian-based CI is not meaningful (locally non-quadratic surface), so
    # we suppress it and tag the line instead.
    bound_status = result.bound_status or [''] * len(result.params_opt)
    for i, (label, val, stderr) in enumerate(zip(labels, result.params_opt, result.params_stderr)):
        status = bound_status[i] if i < len(bound_status) else ''
        if status == 'lower' or status == 'upper':
            logger.info(f"    {label:5s} = {val:.2e}  [at {status} bound — CI not meaningful]")
        elif status == 'fixed':
            logger.info(f"    {label:5s} = {val:.2e}  [fixed]")
        elif np.isinf(stderr) or np.isnan(stderr):
            logger.info(f"    {label:5s} = {val:.2e} +/- inf")
        else:
            low, high = ci_low[i], ci_high[i]
            logger.info(f"    {label:5s} = {val:.2e} +/- {stderr:.2e}  [95% CI: {low:.2e}, {high:.2e}]")

    # Fit quality
    logger.info(f"  Fit error: {result.fit_error_rel:.2f}% (rel), {result.fit_error_abs:.2f} Ohm (abs)")
    logger.info(f"  Quality: {result.quality.capitalize()} (<10.0%)")

    # Warnings
    for warning in result.all_warnings:
        logger.warning(f"  {warning}")

    log_separator()


# =============================================================================
# Circuit Fitting
# =============================================================================

def run_circuit_fitting(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Run equivalent circuit fitting.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_fit, input, circuit, voigt_chain, weighting,
                       optimizer, multistart, de_*, save, format)

    Returns
    -------
    result : FitResult or None
        Fitting result
    fig : Figure or None
        Fit figure
    """
    # Check if fitting should be skipped
    if args.no_fit:
        return None, None

    # Auto-skip for custom data without circuit specification
    if args.input is not None and args.circuit is None and not args.voigt_chain:
        logger.warning("Fit skipped (custom data without --circuit or --voigt-chain)")
        logger.info("For fitting, use: --circuit 'R(100) - (R(1000) | C(1e-6))' "
                    "or --voigt-chain")
        return None, None

    # Default circuit for synthetic data: Rs - (R||Q) - (R||Q)
    circuit_expr = args.circuit
    if args.input is None and circuit_expr is None:
        circuit_expr = 'R() - (R() | Q()) - (R() | Q())'

    # Voigt chain fitting
    if args.voigt_chain:
        return _fit_voigt_chain(frequencies, Z, args)
    else:
        return _fit_standard_circuit(frequencies, Z, args, circuit_expr)


def _fit_voigt_chain(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Fit using Voigt chain linear method.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments for Voigt chain options

    Returns
    -------
    result : FitResult or None
        Fitting result
    fig : Figure or None
        Fit figure
    """
    log_separator()
    logger.info("Using Voigt chain linear fit")
    log_separator()

    try:
        circuit, initial_params = fit_voigt_chain_linear(
            frequencies, Z,
            n_per_decade=args.voigt_n_per_decade,
            extend_decades=args.voigt_extend_decades,
            include_L=not args.voigt_no_inductance,
            fit_type=args.voigt_fit_type,
            prune_threshold=args.voigt_prune_threshold,
            allow_negative=args.voigt_allow_negative,
            auto_optimize_M=args.voigt_auto_M,
            mu_threshold=args.voigt_mu_threshold,
            max_M=args.voigt_max_M,
            weighting=args.weighting
        )
        logger.info(f"Created circuit: {circuit}")
        logger.info(f"Number of parameters: {len(initial_params)}")
    except Exception as e:
        raise EISAnalysisError(f"Voigt chain linear fit error: {e}") from e

    log_separator()
    logger.info("Using linear result (no nonlinear optimization)")
    log_separator()

    # Compute Z_fit
    Z_fit = circuit.impedance(frequencies, initial_params)

    # Compute fit metrics using shared utility function
    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(
        Z, Z_fit, args.weighting
    )

    logger.info("Fit error:")
    logger.info(f"  Relative:  {fit_error_rel:.2f}%")
    logger.info(f"  Absolute:  {fit_error_abs:.2e} Ohm")
    logger.info(f"  Quality:   {quality}")

    # Update circuit with fitted parameters
    if hasattr(circuit, 'update_params'):
        circuit.update_params(initial_params)

    # Create FitResult
    result = FitResult(
        circuit=circuit,
        params_opt=np.array(initial_params),
        params_stderr=np.zeros(len(initial_params)),
        fit_error_rel=fit_error_rel,
        fit_error_abs=fit_error_abs,
        quality=quality,
        _n_data=len(frequencies)
    )

    # Generate interpolated frequencies for smooth curve
    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = circuit.impedance(freq_plot, initial_params)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(Z.real, -Z.imag, 'o', label='Data', markersize=6, alpha=0.7)
    ax.plot(Z_fit_plot.real, -Z_fit_plot.imag, '-', label='Linear fit', linewidth=2)
    ax.set_xlabel("Z' [Ohm]")
    ax.set_ylabel("-Z'' [Ohm]")
    ax.set_title(f"Voigt chain - Linear fit (error: {fit_error_rel:.2f}%)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    plt.tight_layout()

    save_figure(fig, args.save, 'fit', args.format)

    return result, fig


def _fit_standard_circuit(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    circuit_expr: str
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Fit using standard circuit with nonlinear optimization.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments for fitting options
    circuit_expr : str
        Circuit expression string

    Returns
    -------
    result : FitResult or None
        Fitting result
    fig : Figure or None
        Fit figure
    """
    try:
        circuit = parse_circuit_expression(circuit_expr)
    except ValueError as e:
        raise EISAnalysisError(
            f"Circuit parsing error: {e}\n"
            "Use syntax like: --circuit 'R(100) - (R(5000) | C(1e-6))'"
        ) from e

    # Log circuit expression
    logger.info("=" * 60)
    logger.info("Equivalent circuit")
    logger.info("=" * 60)
    logger.info(f"Circuit: {circuit_expr}")

    try:
        if args.optimizer == 'de':
            # Differential Evolution global optimization
            diffevo_result, Z_fit, fig = fit_circuit_diffevo(
                circuit,
                frequencies, Z,
                strategy=args.de_strategy,
                popsize=args.de_popsize,
                maxiter=args.de_maxiter,
                tol=args.de_tol,
                workers=args.de_workers,
                weighting=args.weighting,
                verbose=True,
                use_analytic_jacobian=not args.numeric_jacobian
            )
            result = diffevo_result.best_result

            # Log DE diagnostics
            _log_diffevo_diagnostics(diffevo_result)

        elif args.optimizer == 'multistart':
            # Multi-start optimization
            n_restarts = args.multistart if args.multistart > 0 else 16
            multistart_result, Z_fit, fig = fit_circuit_multistart(
                circuit,
                frequencies, Z,
                n_restarts=n_restarts,
                scale=args.multistart_scale,
                weighting=args.weighting,
                verbose=True,
                use_analytic_jacobian=not args.numeric_jacobian
            )
            result = multistart_result.best_result

            # Log multistart diagnostics
            _log_multistart_diagnostics(multistart_result)

        else:
            # Single fit (args.optimizer == 'single')
            result, Z_fit, fig = fit_equivalent_circuit(
                frequencies, Z,
                circuit,
                weighting=args.weighting,
                use_analytic_jacobian=not args.numeric_jacobian
            )

        # Log fit results
        _log_fit_result(result)

        save_figure(fig, args.save, 'fit', args.format)
        return result, fig

    except Exception as e:
        logger.error(f"Fitting error: {e}")
        logger.error("Try adjusting --circuit expression")
        logger.debug("Traceback:", exc_info=True)
        return None, None
