"""
Kramers-Kronig validation for EIS data quality assessment.

Clean design: No logging in core functions, all diagnostics returned as data.

Provides two implementations:
1. lin_kk_native() - Native implementation using Voigt chain (no external dependencies)
2. kramers_kronig_validation() - High-level wrapper with visualization
"""

import numpy as np
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass, field
from typing import Tuple, Optional, List
from numpy.typing import NDArray

logger = logging.getLogger(__name__)

# The mu search starts at the first M whose log10(pseudo chi^2) is within
# CHI2_PLATEAU_DECADES of the lowest value over the next CHI2_PLATEAU_WINDOW
# M. Below that point the tau grid is too coarse and mu dips from
# discretization, not from fitting noise (doc/KRAMERS_KRONIG_REVIEW.md 1.1).
# 0.3 = a factor of 2 in chi^2, RMS residuals within sqrt(2) of what more
# elements reach; 1.0 was measured to stop too early on 1 % noise.
# The window is local on purpose: on drifting (non-KK) data chi^2 keeps
# creeping down as negative R_k absorb the drift, and the minimum over all M
# pushed the start to M ~ 45-49 in 2 of 6 synthetic drift spectra, hiding a
# 10 % drift. It must still span the chi^2 bumps of a coarse grid (a few M):
# 5 stopped early on 1 % noise, 8 matched the all-M minimum on every
# KK-compliant spectrum tested.
CHI2_PLATEAU_DECADES = 0.3
CHI2_PLATEAU_WINDOW = 8

# is_valid counts the points whose larger residual component, max(|res_real|,
# |res_imag|), exceeds KK_RESIDUAL_THRESHOLD [%] - the +-lines of the plot -
# and passes the spectrum when at most KK_MAX_FRACTION_ABOVE of them do.
# A mean over all points (the criterion before) hid local violations:
# example/real_gamry_example.DTA violates KK over 0.03-4 Hz (residual hump
# to 20 % with every fit_type) and passed with a 3.8 % mean, 19 of 72 points
# above 5 % (doc/KRAMERS_KRONIG_REVIEW.md 1.4).
# Both values are empirical, a loose "clearly broken" bound, not a noise-level
# test. The allowance tolerates a few edge points, where Lin-KK residuals grow.
# Any threshold between ~3 and ~15 % separates the example spectra, so
# neither value is a measured optimum.
KK_RESIDUAL_THRESHOLD = 5.0
KK_MAX_FRACTION_ABOVE = 0.05


def count_points_above_threshold(
    residuals_real: NDArray[np.float64],
    residuals_imag: NDArray[np.float64]
) -> int:
    """Count points whose larger |residual| component exceeds KK_RESIDUAL_THRESHOLD [%].

    NaN residuals count as above: a point that could not be fitted is not a pass.
    """
    per_point = 100 * np.maximum(np.abs(residuals_real), np.abs(residuals_imag))
    return int(np.sum(~(per_point <= KK_RESIDUAL_THRESHOLD)))


@dataclass
class KKResult:
    """Result of Kramers-Kronig validation.

    Attributes
    ----------
    M : int
        Number of Voigt elements used (0 if failed)
    M_lower : int
        First M the mu search tried: where pseudo chi^2 levels off (0 if failed)
    mu : float
        Lin-KK stop value: mu at the first M >= M_lower where it dropped
        below mu_threshold, so it is expected to be below the threshold on
        normal termination. Not a data-quality metric (judge quality by
        the residuals); mu > threshold only when max_M was reached.
        With extend_decades > 0 the returned model's own mu is at least
        this value: extensions that would lower it are rejected.
    Z_fit : NDArray[np.complex128] or None
        Fitted impedance
    residuals_real : NDArray[np.float64] or None
        Real part residuals (normalized by |Z|)
    residuals_imag : NDArray[np.float64] or None
        Imaginary part residuals (normalized by |Z|)
    pseudo_chisqr : float
        Pseudo chi-squared (Boukamp 1995)
    noise_estimate : float
        Estimated noise in percent (Yrjana & Bobacka 2024)
    extend_decades : float
        Tau range extension in decades (0.0 = no extension)
    inductance : Optional[float]
        Fitted series inductance [H]
    capacitance : Optional[float]
        Fitted series capacitance [F] (None unless include_C was requested)
    figure : Optional[plt.Figure]
        Visualization figure
    warnings : List[str]
        Warning messages
    error : Optional[str]
        Error message if validation failed
    """
    M: int = 0
    M_lower: int = 0
    mu: float = 0.0
    Z_fit: Optional[NDArray[np.complex128]] = None
    residuals_real: Optional[NDArray[np.float64]] = None
    residuals_imag: Optional[NDArray[np.float64]] = None
    pseudo_chisqr: float = 0.0
    noise_estimate: float = 0.0
    extend_decades: float = 0.0
    inductance: Optional[float] = None
    capacitance: Optional[float] = None
    figure: Optional[plt.Figure] = None
    warnings: List[str] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        """Check if KK validation completed successfully."""
        return self.Z_fit is not None and self.residuals_real is not None

    @property
    def mean_residual_real(self) -> float:
        """Mean absolute real residual in percent."""
        if self.residuals_real is None:
            return float('inf')
        return float(np.mean(np.abs(self.residuals_real)) * 100)

    @property
    def mean_residual_imag(self) -> float:
        """Mean absolute imaginary residual in percent."""
        if self.residuals_imag is None:
            return float('inf')
        return float(np.mean(np.abs(self.residuals_imag)) * 100)

    @property
    def n_above_threshold(self) -> int:
        """Points with max(|res_real|, |res_imag|) above KK_RESIDUAL_THRESHOLD."""
        if self.residuals_real is None or self.residuals_imag is None:
            return 0
        return count_points_above_threshold(self.residuals_real, self.residuals_imag)

    @property
    def is_valid(self) -> bool:
        """Check if data passes KK validation (at most KK_MAX_FRACTION_ABOVE
        of the points above KK_RESIDUAL_THRESHOLD %)."""
        if not self.success or self.residuals_real is None:
            return False
        return self.n_above_threshold <= KK_MAX_FRACTION_ABOVE * len(self.residuals_real)


@dataclass
class LinKKResult:
    """Result of Lin-KK native fitting.

    Attributes
    ----------
    M : int
        Number of Voigt elements used
    mu : float
        Lin-KK stop value: mu at the first M >= M_lower where it dropped
        below mu_threshold, so it is expected to be below the threshold on
        normal termination. Not a data-quality metric (judge quality by
        the residuals); mu > threshold only when max_M was reached.
        With extend_decades > 0 the returned model's own mu is at least
        this value: extensions that would lower it are rejected.
    Z_fit : NDArray[np.complex128]
        Fitted impedance
    residuals_real : NDArray[np.float64]
        Real part residuals (normalized by |Z|)
    residuals_imag : NDArray[np.float64]
        Imaginary part residuals (normalized by |Z|)
    pseudo_chisqr : float
        Pseudo chi-squared (Boukamp 1995)
    noise_estimate : float
        Estimated noise in percent
    extend_decades : float
        Tau range extension in decades
    inductance : Optional[float]
        Fitted series inductance [H]
    elements : NDArray[np.float64]
        Fitted elements [R_s, R_1, ..., R_M]
    tau : NDArray[np.float64]
        Time constants [s]
    M_lower : int
        First M the mu search tried: where pseudo chi^2 levels off
    weighting : str
        Weighting scheme used
    capacitance : Optional[float]
        Fitted series capacitance [F] (None unless include_C was requested;
        not part of the elements array)
    warnings : List[str]
        Caveats about the fit (e.g. max_M reached, extension rejected)
    """
    M: int
    mu: float
    Z_fit: NDArray[np.complex128]
    residuals_real: NDArray[np.float64]
    residuals_imag: NDArray[np.float64]
    pseudo_chisqr: float
    noise_estimate: float
    extend_decades: float
    inductance: Optional[float]
    elements: NDArray[np.float64]
    tau: NDArray[np.float64]
    M_lower: int
    weighting: str = 'modulus'
    capacitance: Optional[float] = None
    warnings: List[str] = field(default_factory=list)

    @property
    def mean_residual_real(self) -> float:
        """Mean absolute real residual in percent."""
        return float(np.mean(np.abs(self.residuals_real)) * 100)

    @property
    def mean_residual_imag(self) -> float:
        """Mean absolute imaginary residual in percent."""
        return float(np.mean(np.abs(self.residuals_imag)) * 100)

    @property
    def n_above_threshold(self) -> int:
        """Points with max(|res_real|, |res_imag|) above KK_RESIDUAL_THRESHOLD."""
        return count_points_above_threshold(self.residuals_real, self.residuals_imag)

    @property
    def is_valid(self) -> bool:
        """Check if data passes KK validation (at most KK_MAX_FRACTION_ABOVE
        of the points above KK_RESIDUAL_THRESHOLD %)."""
        return self.n_above_threshold <= KK_MAX_FRACTION_ABOVE * len(self.residuals_real)


def compute_pseudo_chisqr(
    Z_exp: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128]
) -> float:
    """
    Compute pseudo chi-squared (Boukamp 1995).

    Parameters
    ----------
    Z_exp : array
        Experimental impedance
    Z_fit : array
        Fitted impedance

    Returns
    -------
    float
        Pseudo chi-squared value

    References
    ----------
    Boukamp, B.A. "A Linear Kronig-Kramers Transform Test for Immittance
    Data Validation." J. Electrochem. Soc. 142, 1885-1894 (1995)
    """
    weight = 1.0 / (Z_exp.real**2 + Z_exp.imag**2)
    return float(np.sum(weight * (
        (Z_exp.real - Z_fit.real)**2 +
        (Z_exp.imag - Z_fit.imag)**2
    )))


def estimate_noise_percent(chi2_ps: float, n_points: int) -> float:
    """
    Estimate noise standard deviation from pseudo chi-squared.

    Based on Yrjana & Bobacka (2024).

    Parameters
    ----------
    chi2_ps : float
        Pseudo chi-squared value
    n_points : int
        Number of data points

    Returns
    -------
    float
        Estimated noise in percent

    References
    ----------
    Yrjana, V. and Bobacka, J. "Implementing Kramers-Kronig validity testing
    using pyimpspec." Electrochim. Acta 504, 144951 (2024)
    """
    return float(np.sqrt(chi2_ps * 5000 / n_points))


def reconstruct_impedance(
    frequencies: NDArray[np.float64],
    elements: NDArray[np.float64],
    tau: NDArray[np.float64],
    L_value: Optional[float],
    include_L: bool = True,
    C_value: Optional[float] = None
) -> NDArray[np.complex128]:
    """
    Reconstruct impedance from fitted Voigt elements.

    Parameters
    ----------
    frequencies : array
        Frequencies [Hz]
    elements : array
        Fitted elements [R_s, R_1, ..., R_M, (L)]
    tau : array
        Time constants [s]
    L_value : float or None
        Inductance value [H]
    include_L : bool
        Whether inductance is included in elements array
    C_value : float or None, optional
        Series capacitance [F]; adds 1/(j*omega*C). Unlike L, it is never
        part of the elements array (default: None = no series C term).

    Returns
    -------
    Z_fit : array
        Reconstructed complex impedance
    """
    omega = 2 * np.pi * frequencies
    R_s = elements[0]
    R_i_end = -1 if include_L else len(elements)
    R_i = elements[1:R_i_end]

    Z_fit = np.full_like(frequencies, R_s, dtype=complex)
    for r, t in zip(R_i, tau):
        Z_fit += r / (1 + 1j * omega * t)

    if include_L and L_value is not None:
        Z_fit += 1j * omega * L_value

    if C_value is not None:
        Z_fit += 1.0 / (1j * omega * C_value)

    return Z_fit


def find_optimal_extend_decades(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    M: int,
    search_range: Tuple[float, float] = (0.0, 1.0),
    n_evaluations: int = 11,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'real',
    weighting: str = 'modulus',
    min_mu: float = -np.inf
) -> Optional[Tuple[float, float, NDArray[np.float64], NDArray[np.float64], Optional[float], Optional[float]]]:
    """
    Find optimal extend_decades that minimizes pseudo chi-squared.

    Uses grid search over the specified range.

    Parameters
    ----------
    frequencies : array
        Measured frequencies [Hz]
    Z : array
        Measured impedance
    M : int
        Number of Voigt elements
    search_range : tuple
        Min and max extend_decades to search
    n_evaluations : int
        Number of grid points
    include_L : bool
        Include series inductance
    include_C : bool
        Include series capacitance (blocking low-frequency behavior)
    fit_type : str
        Fit type ('real', 'imag', 'complex')
    weighting : str
        Weighting scheme
    min_mu : float
        Reject candidates whose own mu falls below this value. Lin-KK passes
        its stop mu: at fixed M, a wider tau grid can fit with oscillating
        negative R_i, the overfit the mu criterion exists to prevent.
        -inf (default) keeps every candidate.

    Returns
    -------
    None
        Only with min_mu, when every candidate falls below it.
    optimal_extend_decades : float
        Value that minimizes chi^2
    min_chi2 : float
        Minimum chi^2 achieved
    tau : array
        Time constants for optimal extend_decades
    elements : array
        Fitted elements for optimal extend_decades
    L_value : float or None
        Inductance for optimal extend_decades
    C_value : float or None
        Series capacitance for optimal extend_decades
    """
    from ..fitting.voigt_chain import generate_tau_grid_fixed_M, estimate_R_linear, calc_mu

    candidates = np.linspace(search_range[0], search_range[1], n_evaluations)
    results = []

    for ext_dec in candidates:
        tau = generate_tau_grid_fixed_M(frequencies, M, extend_decades=ext_dec)
        elements, residual, L_value, C_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=True, include_L=include_L, include_C=include_C,
            fit_type=fit_type, allow_negative=True,
            weighting=weighting
        )

        # R_1..R_M sit right after R_s whether or not L follows them
        if calc_mu(elements[1:1 + len(tau)]) < min_mu:
            continue

        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)
        chi2 = compute_pseudo_chisqr(Z, Z_fit)
        results.append((ext_dec, chi2, tau, elements, L_value, C_value))

    if not results:
        return None

    # Find minimum chi^2
    min_chi2 = min(r[1] for r in results)
    tolerance = 0.001 * min_chi2
    near_optimal = [r for r in results if r[1] <= min_chi2 + tolerance]
    best = min(near_optimal, key=lambda x: abs(x[0]))
    return best[0], best[1], best[2], best[3], best[4], best[5]


def _chi2_lower_M(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    max_M: int,
    include_L: bool,
    include_C: bool,
    fit_type: str,
    weighting: str
) -> int:
    """First M whose pseudo chi^2 is within CHI2_PLATEAU_DECADES of the
    lowest over the next CHI2_PLATEAU_WINDOW M (unextended tau grid, as the
    mu search)."""
    from ..fitting.voigt_chain import generate_tau_grid_fixed_M, estimate_R_linear

    Ms = np.arange(3, max_M + 1)
    log_chi2 = np.empty(len(Ms))
    for i, M in enumerate(Ms):
        tau = generate_tau_grid_fixed_M(frequencies, M)
        elements, _, L_value, C_value = estimate_R_linear(
            frequencies, Z, tau,
            include_Rs=True, include_L=include_L, include_C=include_C,
            fit_type=fit_type, allow_negative=True, weighting=weighting
        )
        Z_fit = reconstruct_impedance(frequencies, elements, tau, L_value, include_L, C_value=C_value)
        # Exact synthetic data can fit to chi^2 = 0
        log_chi2[i] = np.log10(np.maximum(compute_pseudo_chisqr(Z, Z_fit), np.finfo(float).tiny))

    # A point with Z = 0 makes every chi^2 infinite and every comparison NaN:
    # no plateau, so start where the original Lin-KK does
    return int(next((M for i, M in enumerate(Ms)
                     if log_chi2[i] - log_chi2[i:i + CHI2_PLATEAU_WINDOW + 1].min() <= CHI2_PLATEAU_DECADES),
                    Ms[0]))


def lin_kk_native(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    include_L: bool = True,
    include_C: bool = False,
    fit_type: str = 'real',
    weighting: str = 'modulus',
    auto_extend_decades: bool = False,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0)
) -> LinKKResult:
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
    LinKKResult
        Dataclass containing M, mu, Z_fit, residuals, pseudo_chisqr,
        noise_estimate, extend_decades, inductance, elements, tau.
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
            n_evaluations=11,
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

    return LinKKResult(
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


def kramers_kronig_validation(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    mu_threshold: float = 0.85,
    max_M: int = 50,
    auto_extend_decades: bool = True,
    extend_decades_range: Tuple[float, float] = (0.0, 1.0),
    include_C: bool = False
) -> KKResult:
    """
    Perform Kramers-Kronig validation test on EIS data.

    Uses native Lin-KK implementation (Schönleber et al. 2014).

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
        value (expected below mu_threshold on normal termination).
    """
    try:
        lkk = lin_kk_native(
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
    except Exception as e:
        logger.debug(f"KK validation error: {e}")
        return KKResult(error=str(e))

    if lkk.Z_fit is None:
        return KKResult(error="KK fitting failed - could not fit Voigt chain")

    # Generate interpolated frequencies for smooth curve
    f_min, f_max = frequencies.min(), frequencies.max()
    freq_plot = np.logspace(np.log10(f_min), np.log10(f_max), 300)
    Z_fit_plot = reconstruct_impedance(freq_plot, lkk.elements, lkk.tau, lkk.inductance, include_L=True,
                                       C_value=lkk.capacitance)

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Fit comparison (Nyquist plot)
    ax1 = axes[0]
    ax1.plot(Z.real, -Z.imag, 'o', label='Data', markersize=4)
    ax1.plot(Z_fit_plot.real, -Z_fit_plot.imag, '-', label='KK fit', linewidth=2)
    ax1.set_xlabel("Z' [Ohm]")
    ax1.set_ylabel("-Z'' [Ohm]")
    ax1.set_title(f"Kramers-Kronig fit in real domain (M={lkk.M})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal', adjustable='datalim')

    # Residuals plot
    ax2 = axes[1]
    ax2.semilogx(frequencies, lkk.residuals_real * 100, 'o', label='Real', markersize=4)
    ax2.semilogx(frequencies, lkk.residuals_imag * 100, 's', label='Imaginary', markersize=4)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax2.axhline(y=KK_RESIDUAL_THRESHOLD, color='r', linestyle=':', alpha=0.5)
    ax2.axhline(y=-KK_RESIDUAL_THRESHOLD, color='r', linestyle=':', alpha=0.5)
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Residuals [%]")
    ax2.set_title(f"KK residuals (stop mu={lkk.mu:.3f}, chi^2={lkk.pseudo_chisqr:.2e}, noise~{lkk.noise_estimate:.1f}%)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return KKResult(
        M=lkk.M,
        M_lower=lkk.M_lower,
        mu=lkk.mu,
        Z_fit=lkk.Z_fit,
        residuals_real=lkk.residuals_real,
        residuals_imag=lkk.residuals_imag,
        pseudo_chisqr=lkk.pseudo_chisqr,
        noise_estimate=lkk.noise_estimate,
        extend_decades=lkk.extend_decades,
        inductance=lkk.inductance,
        capacitance=lkk.capacitance,
        figure=fig,
        warnings=lkk.warnings
    )
