"""
Analysis workflow handlers for EIS CLI.

Each handler function corresponds to a step in the analysis pipeline:
- run_kk_validation: Kramers-Kronig validation
- run_zhit_validation: Z-HIT validation
- run_rinf_estimation: R_inf estimation
- run_drt_analysis: DRT analysis
- run_voigt_analysis: Voigt element analysis from DRT
- run_circuit_fitting: Equivalent circuit fitting
- run_oxide_analysis: Oxide layer analysis
"""

import argparse
import logging
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from .logging import log_separator
from .utils import EISAnalysisError, CLIDRTResult, save_figure, parse_circuit_expression
from ..validation import kramers_kronig_validation, zhit_validation
from ..drt import calculate_drt
from ..fitting import (
    fit_equivalent_circuit,
    fit_circuit_multistart,
    fit_circuit_diffevo,
    fit_voigt_chain_linear,
    analyze_voigt_elements,
    format_voigt_report,
    FitResult,
)
from ..fitting.diagnostics import compute_fit_metrics
from ..analysis import analyze_oxide_layer
from ..rinf_estimation import estimate_rinf_with_inductance

logger = logging.getLogger(__name__)


# =============================================================================
# Kramers-Kronig Validation
# =============================================================================

def run_kk_validation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[plt.Figure]:
    """
    Run Kramers-Kronig validation.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_kk, mu_threshold, auto_extend, extend_decades_max, save, format)

    Returns
    -------
    fig : Figure or None
        KK validation figure
    """
    if args.no_kk:
        return None

    result = kramers_kronig_validation(
        frequencies, Z,
        mu_threshold=args.mu_threshold,
        auto_extend_decades=args.auto_extend,
        extend_decades_range=(0.0, args.extend_decades_max)
    )
    if result is None:
        return None

    save_figure(result.figure, args.save, 'kk', args.format)
    return result.figure


# =============================================================================
# Z-HIT Validation
# =============================================================================

def run_zhit_validation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[plt.Figure]:
    """
    Run Z-HIT validation.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_zhit, zhit_optimize_offset, save, format)

    Returns
    -------
    fig : Figure or None
        Z-HIT validation figure
    """
    if args.no_zhit:
        return None

    result = zhit_validation(
        frequencies, Z,
        optimize_offset=args.zhit_optimize_offset
    )

    save_figure(result.figure, args.save, 'zhit', args.format)
    return result.figure


# =============================================================================
# R_inf Estimation
# =============================================================================

def run_rinf_estimation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[float], Optional[plt.Figure]]:
    """
    Run R_inf estimation if --ri-fit is specified.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: ri_fit, save, format)

    Returns
    -------
    R_inf : float or None
        Estimated R_inf value
    fig : Figure or None
        R_inf fit figure
    """
    if not args.ri_fit:
        return None, None

    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    try:
        R_inf, L_fit, circuit_rl, diag_rl, fig = estimate_rinf_with_inductance(
            frequencies, Z, verbose=False, plot=True
        )

        # Print results
        method = diag_rl.get('method', 'unknown')
        logger.info(f"R_inf = {R_inf:.3f} Ohm ({diag_rl['n_points_used']} HF points)")

        if 'r_squared' in diag_rl and diag_rl['r_squared'] > 0:
            logger.info(f"  Quality: R^2 = {diag_rl['r_squared']:.4f}")
        if 'L_nH' in diag_rl and diag_rl['L_nH'] > 0:
            logger.info(f"  Inductance: L = {diag_rl['L_nH']:.2f} nH")

        # For Voigt fit: show R_ct and C
        if 'voigt' in method and 'R_ct' in diag_rl and 'C_nF' in diag_rl:
            logger.info(f"  Voigt params: R_ct = {diag_rl['R_ct']:.2f} Ohm, "
                        f"C = {diag_rl['C_nF']:.2f} nF")
            if 'f_characteristic' in diag_rl:
                logger.info(f"  Characteristic freq: f_char = "
                            f"{diag_rl['f_characteristic']/1e6:.3f} MHz")

        # Comparison with median
        n_avg = min(5, max(1, len(frequencies) // 10))
        high_freq_indices = np.argsort(frequencies)[-n_avg:]
        R_inf_median = np.median(Z.real[high_freq_indices])
        diff_abs = R_inf - R_inf_median
        diff_pct = (diff_abs / R_inf_median * 100) if R_inf_median != 0 else 0
        logger.info(f"  For comparison: median = {R_inf_median:.3f} Ohm "
                    f"(diff: {diff_abs:+.3f} Ohm, {diff_pct:+.1f}%)")

        # Warnings
        if diag_rl.get('warnings'):
            for warning in diag_rl['warnings']:
                logger.warning(f"  {warning}")

        log_separator()
        save_figure(fig, args.save, 'ri_fit', args.format)

        return R_inf, fig

    except Exception as e:
        logger.error(f"R_inf estimation failed: {e}")
        logger.debug("Traceback:", exc_info=True)
        return None, None


# =============================================================================
# DRT Analysis
# =============================================================================

def run_drt_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    R_inf_computed: Optional[float],
    peak_method: str
) -> CLIDRTResult:
    """
    Run DRT analysis.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_drt, lambda_reg, n_tau, normalize_rpol, ri_fit,
                       gmm_bic_threshold, save, format)
    R_inf_computed : float or None
        Pre-computed R_inf from --ri-fit
    peak_method : str
        Peak detection method ('scipy' or 'gmm')

    Returns
    -------
    CLIDRTResult
        Container with tau, gamma, peaks, and figures
    """
    if args.no_drt:
        return CLIDRTResult(
            tau=None, gamma=None, peaks_gmm=None, fig_drt=None, fig_rinf=None
        )

    use_auto_lambda = args.lambda_reg is None

    tau, gamma, fig_drt, peaks_gmm, fig_rinf = calculate_drt(
        frequencies, Z,
        n_tau=args.n_tau,
        lambda_reg=args.lambda_reg,
        auto_lambda=use_auto_lambda,
        normalize_rpol=args.normalize_rpol,
        peak_method=peak_method,
        use_rl_fit=False,
        use_voigt_fit=args.ri_fit,
        r_inf_preset=R_inf_computed,
        gmm_bic_threshold=args.gmm_bic_threshold
    )

    save_figure(fig_drt, args.save, 'drt', args.format)

    # Save R_inf figure from DRT only if not already saved via --ri-fit
    if not args.ri_fit:
        save_figure(fig_rinf, args.save, 'ri_fit', args.format)

    return CLIDRTResult(
        tau=tau,
        gamma=gamma,
        peaks_gmm=peaks_gmm,
        fig_drt=fig_drt,
        fig_rinf=fig_rinf
    )


# =============================================================================
# Voigt Element Analysis
# =============================================================================

def run_voigt_analysis(
    drt_result: CLIDRTResult,
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> None:
    """
    Run Voigt element analysis from DRT results.

    Parameters
    ----------
    drt_result : CLIDRTResult
        DRT analysis results
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_drt, no_voigt_info, classify_terms)
    """
    if args.no_drt or args.no_voigt_info:
        return
    if drt_result.tau is None or drt_result.gamma is None:
        return

    try:
        voigt_info = analyze_voigt_elements(
            drt_result.tau, drt_result.gamma, frequencies, Z,
            peaks_gmm=drt_result.peaks_gmm,
            classify_terms=args.classify_terms
        )
        report = format_voigt_report(voigt_info)
        logger.info(report)

    except Exception as e:
        logger.warning(f"Voigt element analysis failed: {e}")
        logger.debug(f"Traceback: {e}", exc_info=True)


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
        else:
            # Single fit (args.optimizer == 'single')
            result, Z_fit, fig = fit_equivalent_circuit(
                frequencies, Z,
                circuit,
                weighting=args.weighting,
                use_analytic_jacobian=not args.numeric_jacobian
            )

        save_figure(fig, args.save, 'fit', args.format)
        return result, fig

    except Exception as e:
        logger.error(f"Fitting error: {e}")
        logger.error("Try adjusting --circuit expression")
        logger.debug("Traceback:", exc_info=True)
        return None, None


# =============================================================================
# Oxide Layer Analysis
# =============================================================================

def run_oxide_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    fitted_result: Optional[FitResult],
    metadata: Optional[dict]
) -> None:
    """
    Run oxide layer analysis.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: analyze_oxide, epsilon_r, area)
    fitted_result : FitResult or None
        Circuit fitting result
    metadata : dict or None
        DTA file metadata
    """
    if not args.analyze_oxide:
        return

    # Use area from metadata if available and not explicitly specified
    area_to_use = args.area
    if metadata is not None and metadata.get('area') is not None:
        if args.area == 1.0:  # Default value was not changed
            area_to_use = metadata['area']
            logger.info(f"Using area from DTA metadata: {area_to_use:.4f} cm^2")
        else:
            logger.info(f"Using explicitly specified area: {area_to_use:.4f} cm^2 "
                        f"(metadata: {metadata['area']:.4f} cm^2)")

    analyze_oxide_layer(
        frequencies, Z,
        epsilon_r=args.epsilon_r,
        area_cm2=area_to_use,
        fit_result=fitted_result
    )
