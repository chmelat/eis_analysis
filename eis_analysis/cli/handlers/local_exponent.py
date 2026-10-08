"""
Local CPE exponent handler for the EIS CLI.

- run_local_exponent: map of n(f) from the real part of the admittance
  (--local-exponent)
"""

import argparse
import logging
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import draw_figure
from ...analysis import local_exponent, LocalExponentResult
from ...rinf_estimation import estimate_rinf
from ...visualization import plot_local_exponent

logger = logging.getLogger(__name__)


def run_local_exponent(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[LocalExponentResult]:
    """
    Map the local CPE exponent n(f) if --local-exponent is specified.

    R_inf and L come from `estimate_rinf`, independently of --ri-fit, so the
    map needs no circuit fit.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: local_exponent, save, format)

    Returns
    -------
    LocalExponentResult or None
        None if the option is off or the data are unusable
    """
    if not args.local_exponent:
        return None

    log_separator()
    logger.info("Local CPE exponent n(f) = d ln Re(Y) / d ln(omega)")
    log_separator()

    try:
        est = estimate_rinf(frequencies, Z)
        # L only from a fit the estimator accepted: a rejected fit's L paired
        # with the HF-bound R_inf would be two inconsistent corrections
        L = float(est.fit.params_opt[1]) if est.method == 'rlq_fit' and est.fit is not None else 0.0
        result = local_exponent(frequencies, Z, est.R_inf, L)
    except ValueError as e:
        logger.error(f"Local exponent failed: {e}")
        return None

    used = 'fit' if est.method == 'rlq_fit' else 'HF upper bound'
    logger.info(f"Subtracted: R_inf = {result.R_inf:.4g} Ohm ({used}), "
                f"L = {result.L * 1e9:.3g} nH")
    if est.method != 'rlq_fit':
        logger.warning("  R_inf is not determined, so n near f_max is unreliable")
    for warning in est.warnings:
        logger.warning(f"  {warning}")

    f, n, valid = result.frequencies, result.n, result.valid
    logger.info(f"Determined (uncertainty <= {result.uncertainty_max:g}) at {int(valid.sum())} of {len(f)} "
                f"points, window {result.window_decades:g} decade")
    if valid.any():
        logger.info(f"  n min = {result.n_min:.3f} at {result.f_n_min:.3g} Hz, "
                    f"max = {result.n_max:.3f} at {result.f_n_max:.3g} Hz, "
                    f"span = {result.span:.3f}")
        # One row per decade: the determined point nearest to each power of ten
        log_f = np.log10(f)
        shown = set()
        for decade in range(int(np.floor(log_f[valid].min())),
                            int(np.ceil(log_f[valid].max())) + 1):
            dist = np.where(valid, np.abs(log_f - decade), np.inf)
            i = int(np.argmin(dist))
            if dist[i] <= 0.5 and i not in shown:
                shown.add(i)
                logger.info(f"    {f[i]:10.3g} Hz   n = {n[i]:.3f} "
                            f"+- {result.n_uncertainty[i]:.3f}")

    for warning in result.warnings:
        logger.warning(f"  {warning}")

    draw_figure(lambda: plot_local_exponent(result), 'Local exponent',
                args.save, 'local_exponent', args.format)

    return result
