"""
Fit-to-noise comparison for the EIS CLI: the circuit fit against the Lin-KK
noise of the same points (compare_fit_to_noise).

- noise_adequacy: the comparison, or None without a noise estimate or on failure
- log_noise_adequacy: the 'Fit vs. noise' line of a fit block
- format_noise_ratio: the fit/KK cell of the comparison table

Shared by the single-fit block (fitting.py) and the comparison table
(model_comparison.py), so both read the same verdict.
"""

import logging
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from ...fitting.config import FIT_NOISE_RATIO_WARN, NOISE_FLOOR_PERCENT
from ...fitting.diagnostics import compare_fit_to_noise, NoiseAdequacyResult

logger = logging.getLogger(__name__)


def noise_adequacy(
    Z_measured: Optional[NDArray[np.complexfloating]],
    Z_fit: NDArray[np.complexfloating],
    noise_percent: Optional[float]
) -> Optional[NoiseAdequacyResult]:
    """Fit compared with the Lin-KK noise, or None without KK or on failure.

    Same contract as the residual diagnostics: a diagnostic must never cost
    the user a converged fit or the comparison table.
    """
    if Z_measured is None or noise_percent is None:
        return None
    try:
        return compare_fit_to_noise(Z_measured, Z_fit, noise_percent)
    except Exception as exc:
        logger.debug(f"Fit-to-noise comparison skipped: {exc}", exc_info=True)
        return None


def log_noise_adequacy(a: Optional[NoiseAdequacyResult]) -> None:
    """Log the fit-to-noise line. Pure formatter over NoiseAdequacyResult."""
    if a is None:
        return
    if np.isfinite(a.ratio):
        floor = (f", noise floor {NOISE_FLOOR_PERCENT:g}%"
                 if a.noise_percent < NOISE_FLOOR_PERCENT else "")
        logger.info(f"  Fit vs. noise: {a.fit_noise_percent:.2f}% vs "
                    f"{a.noise_percent:.2f}% (Lin-KK, same points{floor}), "
                    f"ratio {a.ratio:.1f} (limit {FIT_NOISE_RATIO_WARN:g})")
    for warning in a.warnings:
        logger.warning(f"  {warning}")


def format_noise_ratio(a: Optional[NoiseAdequacyResult]) -> str:
    """Six-character fit/KK cell: the ratio, '*' when inadequate, '-' without one.

    Ratios reach 1e4 on clean synthetic data; from 1000 on they switch to one
    significant digit so the cell keeps its width and the columns stay aligned.
    """
    if a is None or not np.isfinite(a.ratio):
        return f"{'-':>6}"
    number = f"{a.ratio:5.1f}" if a.ratio < 1000 else f"{a.ratio:5.0e}"
    return number + (' ' if a.adequate else '*')
