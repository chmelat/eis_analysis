"""
R_inf estimation handler for the EIS CLI.

- run_rinf_estimation: high-frequency resistance estimation (--ri-fit)
"""

import argparse
import logging
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import save_figure
from ...rinf_estimation import estimate_rinf_with_inductance

logger = logging.getLogger(__name__)


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
        fit, fig = estimate_rinf_with_inductance(frequencies, Z, plot=True)

        # Print results
        logger.info(f"R_inf = {fit.R_inf:.3f} Ohm ({fit.n_points_used} HF points)")

        if fit.R_squared > 0:
            logger.info(f"  Quality: R^2 = {fit.R_squared:.4f}")
        if fit.L_nH > 0:
            logger.info(f"  Inductance: L = {fit.L_nH:.2f} nH")

        # Comparison with median
        n_avg = min(5, max(1, len(frequencies) // 10))
        high_freq_indices = np.argsort(frequencies)[-n_avg:]
        R_inf_median = np.median(Z.real[high_freq_indices])
        diff_abs = fit.R_inf - R_inf_median
        diff_pct = (diff_abs / R_inf_median * 100) if R_inf_median != 0 else 0
        logger.info(f"  For comparison: median = {R_inf_median:.3f} Ohm "
                    f"(diff: {diff_abs:+.3f} Ohm, {diff_pct:+.1f}%)")

        for warning in fit.warnings:
            logger.warning(f"  {warning}")

        log_separator()
        save_figure(fig, args.save, 'ri_fit', args.format)

        return fit.R_inf, fig

    except Exception as e:
        logger.error(f"R_inf estimation failed: {e}")
        logger.debug("Traceback:", exc_info=True)
        return None, None
