"""
Kramers-Kronig validation for EIS data quality assessment.

Clean design: No logging in core functions, all diagnostics returned as data.

Two entry points, both returning KKResult:
1. lin_kk_native() - Native implementation using Voigt chain (no external
   dependencies); full control over the fit, raises on failure
2. kramers_kronig_validation() - High-level wrapper with the validation
   defaults; returns a failure as KKResult.error (plot it with
   visualization.plot_kk_validation)
"""

import numpy as np
import logging
from typing import Tuple
from numpy.typing import NDArray

from .noise import compute_pseudo_chisqr, estimate_noise_percent
from .kk_result import KKResult
from .kk_grid import (EXTEND_N_EVALUATIONS, _chi2_lower_M,
                      find_optimal_extend_decades, reconstruct_impedance)

logger = logging.getLogger(__name__)

# noise_estimate / noise_estimate_median above which a few points dominate the
# noise estimate (KKResult warning). Gaussian residuals give ~1.0. Calibration
# (stress-test spectra, 6 families x 40, own noise; share of spectra above):
#
#   noise   clean   + 6 % bad points, sigma 1 / 2.5 / 5 %
#   0.1 %   0.07         0.38 / 0.81 / 0.82
#   0.3 %   0.00         0.01 / 0.15 / 0.65
#   1 %     0.00         0.00 / 0.00 / 0.02
#
# 2.5 halves the clean 0.1 % share (0.03) and loses a third of the detections
# (0.57, 0.34 at 0.3 %/5 %). The clean cases above 2 are all a Lin-KK misfit
# of the low-frequency end (diffusion, oxide, anomalous), where "a few points
# dominate" is true as well. Example files: 0.97-1.24 with noise, T13 3.4,
# EISPOT-M136113-4 2.9, real_gamry_example (drift) 6.3.
NOISE_DOMINATED_RATIO = 2.0


def low_frequency_slope(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating]
) -> float:
    """
    Slope of log10(-Z'') vs log10(f) over the lowest measured decade.

    Negative when -Z'' still grows as the frequency drops: a capacitive end
    (-1 for a capacitor, -1/2 for a semi-infinite Warburg, around 0 for an
    arc whose maximum sits at the lowest frequency). Positive when the arc
    closes and the phase returns toward zero. This is the case for the
    series capacitance (include_C): it fixes a capacitive end, and on a
    closing end it can only absorb drift.

    Returns NaN with fewer than 3 points in the decade or when -Z'' <= 0 at
    any of them (an end at or past zero phase, inductive, or a noisy end
    crossing zero, as on drifting data).
    """
    order = np.argsort(frequencies)
    f, Z = frequencies[order], Z[order]
    in_decade = f <= 10 * f[0]
    minus_im = -Z.imag[in_decade]
    if in_decade.sum() < 3 or np.any(minus_im <= 0):
        return float('nan')
    return float(np.polyfit(np.log10(f[in_decade]), np.log10(minus_im), 1)[0])


def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'real',
    weighting: str = 'modulus',
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0)
) -> KKResult:
    """
    Native Lin-KK implementation using Voigt chain fitting.

    Implements the linear Kramers-Kronig test from Schönleber et al. (2014)
    without external dependencies.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Measured impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85)
    max_M : int, optional
        Maximum number of Voigt elements to try (default: 50). Capped so the
        fitted part keeps a degree of freedom: N - 2 for 'real' (R_s + M
        resistances; at M + 1 = N it interpolates Z' and the predicted Z''
        is unconstrained) and 'complex' (2N equations, so this is
        conservative), N - 2 - [L] - [C] for 'imag' (M resistances plus the
        L and C terms; R_s has no imaginary part).
    include_L : bool, optional
        Include series inductance L (default: True)
    include_C : bool, optional
        Include series capacitance C (default: False). Captures blocking
        (capacitive) low-frequency behavior, e.g. two-electrode cells
        (Schonleber Lin-KK 'add_cap'). A series C is KK-compliant but has
        zero real part, so the Voigt chain cannot represent it.
    fit_type : str, optional
        Fit type: 'real', 'imag', or 'complex'
    weighting : str, optional
        Weighting scheme: 'uniform', 'sqrt', 'proportional', 'modulus'
    auto_extend_decades : bool, optional
        Automatically optimize extend_decades (default: False)
    extend_decades_range : tuple, optional
        Search range for extend_decades optimization

    Returns
    -------
    KKResult
        M, mu, Z_fit, residuals, pseudo_chisqr, noise_estimate,
        extend_decades, inductance, capacitance, elements, tau (never
        with `error` set: failures raise).
        Note: the returned mu is the Lin-KK stop value and is expected
        to be below mu_threshold on normal termination — judge data
        quality by the residuals, not by mu. An extended tau grid is
        accepted only if its own mu is not below it.

    Raises
    ------
    ValueError
        With fewer than 5 points (7 for 'imag' with L and C): no M >= 3
        satisfies the cap on max_M.

    Notes
    -----
    The mu search starts at M_lower, where pseudo chi^2 levels off, not at
    M = 3 (see CHI2_PLATEAU_DECADES).

    References
    ----------
    Schönleber, M. et al. "A Method for Improving the Robustness of linear
    Kramers-Kronig Validity Tests." Electrochimica Acta 131, 20-27 (2014)
    """
    from ..fitting.voigt_chain import find_optimal_M_mu

    n_tail = int(include_L) + int(include_C) if fit_type == 'imag' else 0
    max_M = min(max_M, len(frequencies) - 2 - n_tail)
    if max_M < 3:
        raise ValueError(f"Lin-KK needs at least {5 + n_tail} points, got {len(frequencies)}")
    M_lower = _chi2_lower_M(frequencies, Z, max_M, include_L, include_C, fit_type, weighting)

    # Find optimal M using mu metric. Its progress stays on the result:
    # the KK section reports M and mu in its own summary line.
    mu_opt = find_optimal_M_mu(
        frequencies, Z,
        mu_threshold=mu_threshold,
        max_M=max_M,
        min_M=M_lower,
        extend_decades=0.0,
        include_Rs=True,
        include_L=include_L,
        include_C=include_C,
        fit_type=fit_type,
        allow_negative=True,
        weighting=weighting
    )
    M, mu, tau, elements = mu_opt.M, mu_opt.mu, mu_opt.tau, mu_opt.elements
    L_value, C_value = mu_opt.L_value, mu_opt.C_value

    extend_decades = 0.0
    warnings = list(mu_opt.warnings)

    # Optionally optimize extend_decades, never below the stop mu
    if auto_extend_decades:
        best = find_optimal_extend_decades(
            frequencies, Z, M,
            search_range=extend_decades_range,
            n_evaluations=EXTEND_N_EVALUATIONS,
            include_L=include_L,
            include_C=include_C,
            fit_type=fit_type,
            weighting=weighting,
            min_mu=mu
        )
        if best is None:
            warnings.append(f"Every extend_decades in {extend_decades_range} "
                            f"lowers mu below the stop value {mu:.4f}; "
                            f"kept the unextended tau grid")
        else:
            extend_decades, _, tau, elements, L_value, C_value = best

    # Reconstruct Z_fit from fitted parameters
    Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)

    # Calculate residuals normalized by |Z|
    Z_mag = np.abs(Z)
    Z_mag_safe = np.maximum(Z_mag, 1e-15)

    res_real = (Z.real - Z_fit.real) / Z_mag_safe
    res_imag = (Z.imag - Z_fit.imag) / Z_mag_safe

    # Compute pseudo chi-squared (Boukamp 1995)
    chi2_ps = compute_pseudo_chisqr(Z, Z_fit)
    noise_est = estimate_noise_percent(chi2_ps, len(Z))

    result = KKResult(
        M=M,
        mu=mu,
        Z_fit=Z_fit,
        residuals_real=res_real,
        residuals_imag=res_imag,
        pseudo_chisqr=chi2_ps,
        noise_estimate=noise_est,
        extend_decades=extend_decades,
        inductance=L_value,
        elements=elements,
        tau=tau,
        weighting=weighting,
        capacitance=C_value,
        M_lower=M_lower,
        warnings=warnings
    )

    # A capacitive low-frequency end (-Z'' still rising) is what the series C
    # fixes: blocking/2-electrode cells, Warburg diffusion, an arc continuing
    # past f_min. Without it the residuals and the noise measure the model,
    # not the data (a reflective Warburg: 80 of 81 points above 5 %, "noise"
    # 35 %, against 1.2 % with include_C). The test is the shape of the data,
    # not of the residuals: a Warburg fails on a few edge points with a 1.5 %
    # mean, and on a closing end the series C would only absorb drift
    # (doc/KK_INTUITION.md). Only a failed test: one that passes has a few
    # edge points at most, and its noise stays near the data's.
    # The slope alone also fires on drift (a phase drifting below 1 Hz bends
    # -Z'' upward on a closing arc), so the test is rerun with the series C
    # and kept: drift of |Z| or phase still fails it (14-16 % residuals at
    # the end of synthetic blocking cells), and its residuals then show where.
    # The flag stays either way - on a blocking cell with drift the noise
    # without C is ~50 %, the model's, not the data's.
    if not include_C and not result.is_valid:
        slope = low_frequency_slope(frequencies, Z)
        if slope < 0:
            result.capacitive_end = True
            try:
                result.with_series_c = lin_kk_native(
                    frequencies, Z, mu_threshold=mu_threshold, max_M=max_M,
                    include_L=include_L, include_C=True, fit_type=fit_type,
                    weighting=weighting, auto_extend_decades=auto_extend_decades,
                    extend_decades_range=extend_decades_range)
            except (ValueError, np.linalg.LinAlgError) as exc:
                # One more parameter than the run that succeeded: on very few
                # points it may not fit. The flag alone still holds.
                logger.debug(f"Lin-KK rerun with the series C failed: {exc}")
            with_C = result.with_series_c
            if with_C is None:
                advice = "Rerun with the series capacitance (include_C)."
            elif with_C.is_valid:
                advice = (f"With the series capacitance (include_C) the test "
                          f"passes, noise {with_C.noise_estimate:.1f} %.")
            else:
                advice = (f"The test fails with the series capacitance (include_C) "
                          f"as well, noise {with_C.noise_estimate:.1f} %: the data "
                          f"violate KK beyond the shape of the end - drift or "
                          f"another violation.")
            result.warnings.append(
                f"The low-frequency end is capacitive (-Z'' still rising, slope "
                f"{slope:.2f} per decade), which a Voigt chain cannot follow: the "
                f"residuals and the noise estimate ({noise_est:.1f} %) describe "
                f"the model, not the data. {advice}")
            if with_C is not None:
                result.warnings += [f"With the series capacitance: {w}"
                                    for w in with_C.warnings if w not in result.warnings]

    # A few large residuals (isolated bad points, local drift, or the chain
    # missing the very end of the spectrum) inflate the mean-square noise but
    # not the median. The fit-to-noise comparison divides by the former and
    # then passes a fit those points pulled (cold agent T13: 0.70 % against
    # 0.10 % without the 5 bad points, ratio 1.0 instead of 4.5). Not on a
    # capacitive end: its warning already says the noise is the model's. Not
    # on a sparse fit (KKResult.sparse_fit): its median says nothing.
    from ..fitting.config import NOISE_FLOOR_PERCENT
    noise_median = result.noise_estimate_median
    if (not result.capacitive_end and not result.sparse_fit and np.isfinite(noise_est)
            and noise_median >= NOISE_FLOOR_PERCENT
            and noise_est > NOISE_DOMINATED_RATIO * noise_median):
        result.warnings.append(
            f"The noise estimate ({noise_est:.2f} %) is dominated by a few points: "
            f"from the median residual it is {noise_median:.2f} %. Isolated bad "
            f"points, local drift or a local misfit of the Voigt chain (usually at "
            f"an end of the spectrum); the per-point residual check lists the "
            f"points that stand out individually. Compared with {noise_est:.2f} %, "
            f"a fit those points pull can still look adequate.")
    return result


def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = True,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0),
    include_C: bool = False
) -> KKResult:
    """
    Perform Kramers-Kronig validation test on EIS data.

    Uses native Lin-KK implementation (Schönleber et al. 2014) with the
    validation settings fixed: fit_type='real' (the imaginary residuals are
    the KK prediction), weighting='modulus' and a series L. Call
    lin_kk_native directly to change them.

    Parameters
    ----------
    frequencies : ndarray of float
        Measured frequencies [Hz]
    Z : ndarray of complex
        Complex impedance [Ohm]
    mu_threshold : float, optional
        Threshold for mu metric (default: 0.85)
    max_M : int, optional
        Maximum number of Voigt elements (default: 50)
    auto_extend_decades : bool, optional
        Automatically optimize extend_decades (default: True). Extends the
        Voigt time-constant grid below the lowest measured frequency only,
        which reduces spurious imaginary-part residuals when a relaxation
        continues past it (capacitive low-frequency tail). The high-frequency
        inductive tail is covered by the series L, not by the extension. A
        purely capacitive (blocking) tail needs include_C.
    extend_decades_range : tuple of float, optional
        Search range for extend_decades optimization
    include_C : bool, optional
        Include a series capacitance term 1/(j*omega*C) in the Lin-KK model
        (default: False; Schonleber 'add_cap'). Use for blocking (capacitive)
        low-frequency behavior, e.g. two-electrode cells, where the standard
        Voigt chain produces spurious imaginary residuals at low frequencies.

    Returns
    -------
    KKResult
        Result object with all diagnostics. Check result.success to verify
        validation completed successfully. Data quality is judged by the
        residuals (result.is_valid); result.mu is only the Lin-KK stop
        value (expected below mu_threshold on normal termination). A failed
        test on a capacitive low-frequency end without include_C sets
        result.capacitive_end, carries a warning saying so, and reruns the
        test with include_C into result.with_series_c (result.series_c_passes
        is its verdict): use its residuals and noise there.
    """
    try:
        return lin_kk_native(
            frequencies, Z,
            mu_threshold=mu_threshold,
            max_M=max_M,
            include_L=True,
            include_C=include_C,
            fit_type='real',
            weighting='modulus',
            auto_extend_decades=auto_extend_decades,
            extend_decades_range=extend_decades_range
        )
    # Only what bad data can cause: invalid input or too few points
    # (ValueError) and a failed SVD in the linear fit. A TypeError or
    # AttributeError is a caller's bug and must not pass as "KK failed".
    except (ValueError, np.linalg.LinAlgError) as e:
        logger.debug(f"KK validation error: {e}")
        return KKResult(error=str(e))
