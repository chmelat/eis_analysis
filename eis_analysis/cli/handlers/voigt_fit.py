"""
Voigt-chain fitting handler for the EIS CLI (--voigt-chain).

- _fit_voigt_chain: linear Voigt-chain fit, its report and output
- _log_voigt_chain / _log_mu_optimization: how the chain was built
"""

import argparse
import logging
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import EISAnalysisError
from ...fitting import (
    fit_voigt_chain_linear,
    VoigtChainDiagnostics,
    FitResult,
)
from ...fitting.voigt_chain import MuOptimization
from ...fitting.diagnostics import compute_fit_metrics, compute_significance
from .fit_output import (
    _WEIGHTING_LABELS,
    _log_residual_diagnostics,
    _residual_diagnostics,
    _save_fit_data,
    _save_fit_figure,
)

logger = logging.getLogger(__name__)


def _log_mu_optimization(mu_opt: MuOptimization,
                         diag: VoigtChainDiagnostics) -> None:
    """Report the mu search behind the element count."""
    logger.info("Step 1: Auto-optimizing M using mu metric (Lin-KK)")
    logger.info(f"  mu threshold: {diag.mu_threshold}")
    logger.info(f"  Max M: {diag.max_M}")
    logger.info(f"  Fit type: {diag.fit_type}")
    logger.info(f"  Weighting: {_WEIGHTING_LABELS.get(diag.weighting, diag.weighting)}")
    logger.info(f"  Include L: {diag.include_L}")
    logger.info("  Allow negative R_i: True")
    logger.info("")

    for it in mu_opt.iterations:
        logger.info(f"  Iter {it.iteration:2d}: M={it.M:2d}, mu={it.mu:.4f}, "
                    f"residual={it.residual:.3e}, "
                    f"negative R_i={it.n_negative}/{it.n_R}")

    logger.info("")
    if not mu_opt.reached_max_M:
        logger.info(f"Optimal M found: M = {mu_opt.M}")
        logger.info(f"  mu = {mu_opt.mu:.4f} <= {diag.mu_threshold}")
    for warning in mu_opt.warnings:
        logger.warning(warning)

    if mu_opt.n_negative > 0:
        logger.info(f"  Negative R_i: {mu_opt.n_negative}/{mu_opt.n_R} "
                    f"({mu_opt.n_negative / mu_opt.n_R * 100:.1f}%)")

    log_separator()
    logger.info(f"  Optimal M: {mu_opt.M}, mu: {mu_opt.mu:.4f}")


def _log_voigt_chain(diag: VoigtChainDiagnostics) -> None:
    """Report the four steps behind the Voigt chain initial guess."""
    log_separator()
    logger.info("Voigt chain initial guess estimation (Lin-KK compatible)")
    log_separator()

    # Step 1 and 2: how the tau grid was chosen, and the regression on it
    if diag.mu_optimization is not None:
        _log_mu_optimization(diag.mu_optimization, diag)
        logger.info("Step 2: Refit with NNLS (R_i >= 0) for physical circuit")
    else:
        logger.info(f"Step 1: Generating tau grid ({diag.n_per_decade} per decade, "
                    f"+{diag.extend_decades} dec extension)")
        logger.info(f"  Generated {diag.n_tau_generated} time constants")
        for warning in diag.warnings:
            logger.warning(warning)

        method_str = ("pseudoinverse (allows R_i < 0)" if diag.allow_negative
                      else "NNLS (R_i >= 0)")
        weighting_labels = {
            'uniform': 'uniform (w=1)',
            'sqrt': 'sqrt (w=1/sqrt|Z|)',
            'modulus': 'modulus (w=1/|Z|, Lin-KK standard)',
            'proportional': 'proportional (w=1/|Z|^2)'
        }
        logger.info(f"Step 2: Linear regression (method: {method_str}, "
                    f"fit_type: {diag.fit_type})")
        logger.info(f"  Weighting: "
                    f"{weighting_labels.get(diag.weighting, diag.weighting)}")

    if diag.include_Rs:
        logger.info(f"  R_s (series): {diag.R_s:.3e} Ohm")
    else:
        logger.info("  R_s not included (set to 0)")
    if diag.R_i_min is not None:
        logger.info(f"  R_i range: [{diag.R_i_min:.3e}, {diag.R_i_max:.3e}] Ohm")
    logger.info(f"  Residual: {diag.residual:.3e}")
    if diag.include_L and diag.L_value is not None:
        logger.info(f"  L (inductance): {diag.L_value:.3e} H")

    # Step 3: pruning
    if diag.pruning_enabled:
        n_pruned = diag.n_before_prune - diag.n_after_prune
        logger.info("Step 3: Pruning small R_i")
        logger.info(f"  Relative threshold: {diag.prune_threshold * 100:.1f}% "
                    f"of max = {diag.threshold_relative:.3e} Ohm")
        logger.info(f"  Absolute threshold: 0.1% of R_pol = "
                    f"{diag.threshold_absolute:.3e} Ohm")
        logger.info(f"  Effective threshold: {diag.threshold_effective:.3e} Ohm")
        if n_pruned > 0:
            logger.info(f"  Removed {n_pruned} elements "
                        f"(kept {diag.n_after_prune}/{diag.n_before_prune})")
        else:
            logger.info("  No elements pruned (all above threshold)")
    else:
        logger.info("Step 3: Pruning disabled")

    # Step 4: the circuit
    logger.info(f"Step 4: Building circuit with {diag.n_elements} K elements")
    if diag.include_L and diag.L_value is not None:
        logger.info(f"  Added inductance: L = {diag.L_value:.3e} H")

    logger.info("")
    logger.info("Initial guess summary:")
    if diag.include_Rs:
        logger.info(f"  R_s = {diag.R_s:.3e} Ohm")
    if diag.n_elements > 0:
        logger.info(f"  {diag.n_elements} K elements:")
        logger.info(f"    R_i in [{diag.R_i_min_kept:.3e}, "
                    f"{diag.R_i_max_kept:.3e}] Ohm")
        logger.info(f"    tau_i in [{diag.tau_min:.3e}, {diag.tau_max:.3e}] s")
    if diag.include_L and diag.L_value is not None:
        logger.info(f"  Inductance: L = {diag.L_value:.3e} H")
    logger.info(f"  Total parameters: {diag.n_params}")
    log_separator()


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
        chain = fit_voigt_chain_linear(
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
        circuit, initial_params = chain.circuit, chain.initial_params
        _log_voigt_chain(chain.diagnostics)
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
    circuit.update_params(initial_params)

    # Create FitResult. The linear fit provides no uncertainty estimate, so
    # stderr is inf ("unknown", covariance.py convention), not 0 ("exact");
    # params_ci_95 then honestly returns (-inf, inf) instead of +/- 0 (audit D4).
    result = FitResult(
        circuit=circuit,
        params_opt=np.array(initial_params),
        params_stderr=np.full(len(initial_params), np.inf),
        fit_error_rel=fit_error_rel,
        fit_error_abs=fit_error_abs,
        quality=quality,
        # No Jacobian, so no conditioning: NaN, not the 1.0 default, and
        # not well conditioned - there is no covariance to trust
        condition_number=np.nan,
        is_well_conditioned=False,
        _dof=max(2 * len(frequencies) - len(initial_params), 1),
        # The linear fit gives no uncertainty, but significance needs none -
        # it is a property of the model at its parameters, not of the fit.
        # A pruned Voigt chain is exactly where "which of these 16 elements
        # earns its place?" is worth asking.
        params_significance=compute_significance(circuit, frequencies, initial_params)
    )

    # The Voigt chain is where "too few elements" is directly actionable:
    # --voigt-auto-M or a larger --voigt-n-per-decade is the fix.
    _log_residual_diagnostics(
        _residual_diagnostics(frequencies, Z, Z_fit, args.weighting),
        getattr(args, 'exclude_flagged', None) is None)

    _save_fit_data(frequencies, Z, Z_fit, result, args, 'fit')

    # The full circuit would make the title unreadable: a chain holds up to
    # max_M K elements. What sets this fit apart is that it is linear only.
    return result, _save_fit_figure(
        frequencies, Z, result, args, 'fit',
        title=f"Voigt chain - Linear fit (error: {fit_error_rel:.2f}%)")
