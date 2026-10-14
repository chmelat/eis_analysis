"""
Fit-to-noise comparison for the EIS CLI: the circuit fit against the Lin-KK
noise of the same points (compare_fit_to_noise).

- kk_noise_for_fit: the Lin-KK noise of the fitted points, or None
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
from ...validation import KKResult

logger = logging.getLogger(__name__)


def kk_noise_for_fit(
    kk_result: Optional[KKResult],
    mask: NDArray[np.bool_],
    fit_on: str
) -> Optional[float]:
    """
    Lin-KK noise [%] of the points the circuit is fitted to, or None.

    None without KK, and under --fit-on zhit/all: the fit then runs on the
    Z-HIT reconstruction, whose own error (up to ~3 % at relaxations) the
    Lin-KK noise of the measurement does not describe. Measured on a
    synthetic two-arc spectrum, the correct circuit fitted to the
    reconstruction lands at 2.9-5.0 against either the measured data or the
    reconstruction on drift-free data, and at 0.15-0.20 against the
    reconstruction under 30 % drift - no reference gives a usable ratio.

    On a capacitive low-frequency end without the series C
    (kk_result.capacitive_end) the noise is that of the rerun with it
    (kk_result.with_series_c): without C it is the Voigt chain's misfit (35 %
    on a reflective Warburg with 1 % noise), and a fit would look adequate
    against it. None when that rerun is missing.

    Never raises: a diagnostic must not end the pipeline.
    """
    if kk_result is not None and kk_result.capacitive_end:
        kk_result = kk_result.with_series_c
    if kk_result is None or fit_on != 'original':
        return None
    try:
        return kk_result.noise_estimate_on(mask)
    except ValueError as exc:
        logger.debug(f"Lin-KK noise of the fitted points unavailable: {exc}")
        return None


def noise_adequacy(
    Z: NDArray[np.complexfloating],
    Z_fit: NDArray[np.complexfloating],
    noise_percent: Optional[float]
) -> Optional[NoiseAdequacyResult]:
    """Fit compared with the Lin-KK noise, or None without it or on failure.

    Same contract as the residual diagnostics: a diagnostic must never cost
    the user a converged fit or the comparison table.
    """
    if noise_percent is None:
        return None
    try:
        return compare_fit_to_noise(Z, Z_fit, noise_percent)
    except Exception as exc:
        logger.debug(f"Fit-to-noise comparison skipped: {exc}", exc_info=True)
        return None


def log_noise_adequacy(a: Optional[NoiseAdequacyResult], fit_on: str = 'original',
                       note: Optional[str] = None, source: Optional[str] = None) -> None:
    """Log the fit-to-noise line. Pure formatter over NoiseAdequacyResult.

    Under --fit-on zhit/all, or with `note` (why the Lin-KK noise was not
    usable, see kk_noise_for_fit), says why there is no comparison - an
    absent line would read as an oversight. `source` names the Lin-KK run
    the noise came from, if not the plain one.
    """
    if fit_on != 'original':
        logger.info(f"  Fit vs. noise: not compared (--fit-on {fit_on}: the "
                    "Lin-KK noise does not describe the reconstruction's error)")
        return
    if a is None:
        if note:
            logger.info(f"  Fit vs. noise: not compared ({note})")
        return
    if np.isfinite(a.ratio):
        floor = (f", noise floor {NOISE_FLOOR_PERCENT:g}%"
                 if a.noise_percent < NOISE_FLOOR_PERCENT else "")
        logger.info(f"  Fit vs. noise: {a.fit_noise_percent:.2f}% vs "
                    f"{a.noise_percent:.2f}% ({source or 'Lin-KK'}, same points{floor}), "
                    f"ratio {a.ratio:.1f} (limit {FIT_NOISE_RATIO_WARN:g})")
    for warning in a.warnings:
        logger.warning(f"  {warning}")


def format_noise_ratio(a: Optional[NoiseAdequacyResult]) -> str:
    """Six-character fit/KK cell: the ratio, '*' when inadequate, '-' without one.

    Ratios reach 1e4 on clean synthetic data; from 1000 on (after rounding to
    one decimal, so 999.97 does not become the seven-character '1000.0') they
    switch to one significant digit so the columns stay aligned.
    """
    if a is None or not np.isfinite(a.ratio):
        return f"{'-':>6}"
    number = f"{a.ratio:5.1f}" if round(a.ratio, 1) < 1000 else f"{a.ratio:5.0e}"
    return number + (' ' if a.adequate else '*')
