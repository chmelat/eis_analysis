"""
R_inf estimation handler for the EIS CLI.

- run_rinf_estimation: high-frequency resistance estimation (--ri-fit)
"""

import argparse
import logging
from typing import Optional, Tuple

import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import save_figure
from ...rinf_estimation import estimate_rinf
from ...visualization import plot_rinf_fit

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
        R_inf to hand to the DRT: the fitted R_s, or the HF upper bound
        when the fit does not determine it. None if --ri-fit is
        off or the data are unusable.
    fig : Figure or None
        R_inf fit figure
    """
    if not args.ri_fit:
        return None, None

    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    try:
        est = estimate_rinf(frequencies, Z)
    except ValueError as e:
        logger.error(f"R_inf estimation failed: {e}")
        return None, None

    fit = est.fit
    if fit is not None:
        f_win = est.f_window
        R_fit, stderr = fit.params_opt[0], fit.params_stderr[0]
        logger.info(f"R-L-(R|Q) fit, {f_win.min():.3g}-{f_win.max():.3g} Hz "
                    f"({len(f_win)} points): R_inf = {R_fit:.4g} +- {stderr:.2g} Ohm "
                    f"({100 * stderr / R_fit:.2g} %)")
        logger.info(f"  L = {fit.params_opt[1] * 1e9:.3g} nH, "
                    f"fit error {fit.fit_error_rel:.2g} %")
    logger.info(f"HF upper bound: Re(Z) = {est.R_inf_hf:.4g} Ohm at {est.f_hf:.3g} Hz")
    used = 'fit' if est.method == 'rlq_fit' else 'upper bound'
    logger.info(f"Using R_inf = {est.R_inf:.4g} Ohm ({used})")

    for warning in est.warnings:
        logger.warning(f"  {warning}")

    fig = plot_rinf_fit(est)
    save_figure(fig, args.save, 'ri_fit', args.format)

    return est.R_inf, fig
