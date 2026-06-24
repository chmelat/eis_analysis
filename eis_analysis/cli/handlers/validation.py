"""
Data validation handlers for the EIS CLI.

- run_kk_validation: Kramers-Kronig validation
- run_zhit_validation: Z-HIT validation
"""

import argparse
import logging
from typing import Optional

import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..utils import save_figure
from ...validation import kramers_kronig_validation, zhit_validation
from ...validation.zhit import _quality_label

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

    logger.info("=" * 60)
    logger.info("Kramers-Kronig validation")
    logger.info("=" * 60)

    result = kramers_kronig_validation(
        frequencies, Z,
        mu_threshold=args.mu_threshold,
        auto_extend_decades=args.auto_extend,
        extend_decades_range=(0.0, args.extend_decades_max)
    )
    if not result.success:
        logger.warning(f"KK validation failed: {result.error}")
        return None

    # Summary (format consistent with Z-HIT validation)
    logger.info(f"KK: M={result.M}, mu={result.mu:.4f}, "
                f"extend_decades={result.extend_decades:.2f}")
    logger.info(f"  Mean |res_real|: {result.mean_residual_real:.2f}%")
    logger.info(f"  Mean |res_imag|: {result.mean_residual_imag:.2f}%")
    logger.info(f"  Pseudo chi^2: {result.pseudo_chisqr:.2e}")
    logger.info(f"  Estimated noise (upper bound): {result.noise_estimate:.2f}%")

    mean_abs_residual = max(result.mean_residual_real, result.mean_residual_imag)
    log_fn = logger.info if result.is_valid else logger.warning
    log_fn(f"Data quality: {_quality_label(mean_abs_residual)} "
           f"(max mean |res|={mean_abs_residual:.2f}%, threshold=5.0%)")

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
