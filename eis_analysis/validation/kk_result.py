"""
Result of a Lin-KK test (KKResult) and the pass/fail rule it applies.
"""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from numpy.typing import NDArray

from .noise import estimate_noise_percent, estimate_noise_percent_median


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


@dataclass
class KKResult:
    """Result of a Lin-KK fit: `lin_kk_native` and `kramers_kronig_validation`.

    The two differ only in how they fail: `lin_kk_native` raises, while
    `kramers_kronig_validation` returns a result with `error` set, empty
    arrays (None) and `success == False`.

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
    Z_fit : NDArray[np.complexfloating] or None
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
    elements : Optional[NDArray[np.float64]]
        Fitted elements [R_s, R_1, ..., R_M, L] (L only with include_L;
        the series C is never part of it)
    tau : Optional[NDArray[np.float64]]
        Time constants [s]
    weighting : str
        Weighting scheme of the fit
    warnings : List[str]
        Caveats about the fit (e.g. max_M reached, extension rejected)
    capacitive_end : bool
        The test failed without include_C on a capacitive low-frequency end
        (-Z'' still rising): the residuals and noise_estimate then describe
        the Voigt chain's failure, not the data. A warning says so.
    with_series_c : Optional[KKResult]
        With capacitive_end: the same test rerun with include_C. Its
        residuals and noise describe the data, so the per-point check and
        the fit-to-noise comparison use them in place of these. None when
        not run, or when the rerun raised.
    error : Optional[str]
        Error message if validation failed
    """
    M: int = 0
    M_lower: int = 0
    mu: float = 0.0
    Z_fit: Optional[NDArray[np.complexfloating]] = None
    residuals_real: Optional[NDArray[np.float64]] = None
    residuals_imag: Optional[NDArray[np.float64]] = None
    pseudo_chisqr: float = 0.0
    noise_estimate: float = 0.0
    extend_decades: float = 0.0
    inductance: Optional[float] = None
    capacitance: Optional[float] = None
    elements: Optional[NDArray[np.float64]] = None
    tau: Optional[NDArray[np.float64]] = None
    weighting: str = 'modulus'
    warnings: List[str] = field(default_factory=list)
    capacitive_end: bool = False
    with_series_c: Optional['KKResult'] = None
    error: Optional[str] = None

    @property
    def series_c_passes(self) -> Optional[bool]:
        """Does the test pass with the series C? False: drift or another
        violation beyond the shape of the end. None without the rerun."""
        return None if self.with_series_c is None else self.with_series_c.is_valid

    @property
    def success(self) -> bool:
        """Check if KK validation completed successfully."""
        return self.Z_fit is not None and self.residuals_real is not None

    @property
    def sparse_fit(self) -> bool:
        """More Voigt elements than half the points (M > N/2): the chain
        passes almost through most points, so the median residual collapses
        and says nothing about the noise. On a clean two-arc spectrum with
        10-12 points (M = N - 2) a residual over 8x the median then appears
        on most runs; on the stress-test spectra 9 % are sparse."""
        return self.residuals_real is not None and self.M > len(self.residuals_real) / 2

    @property
    def noise_estimate_median(self) -> float:
        """Noise [%] from the median residual (estimate_noise_percent_median):
        unlike noise_estimate not moved by a few large residuals. nan if
        validation failed."""
        if self.residuals_real is None or self.residuals_imag is None:
            return float('nan')
        return estimate_noise_percent_median(self.residuals_real, self.residuals_imag)

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
        per_point = 100 * np.maximum(np.abs(self.residuals_real), np.abs(self.residuals_imag))
        # NaN counts as above: a point that could not be fitted is not a pass
        return int(np.sum(~(per_point <= KK_RESIDUAL_THRESHOLD)))

    @property
    def is_valid(self) -> bool:
        """Check if data passes KK validation (at most KK_MAX_FRACTION_ABOVE
        of the points above KK_RESIDUAL_THRESHOLD %)."""
        if not self.success or self.residuals_real is None:
            return False
        return self.n_above_threshold <= KK_MAX_FRACTION_ABOVE * len(self.residuals_real)

    def noise_estimate_on(self, mask: NDArray[np.bool_]) -> float:
        """
        noise_estimate [%] over a subset of the points, e.g. a fit's frequency range.

        The residuals are normalized by |Z|, so their sum of squares over the
        subset is that subset's pseudo chi-squared. With every point selected
        this equals noise_estimate.

        Raises
        ------
        ValueError
            If validation failed, the mask does not match the residuals, or
            it selects no point
        """
        if self.residuals_real is None or self.residuals_imag is None:
            raise ValueError("noise_estimate_on: KK validation has no residuals")
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != self.residuals_real.shape:
            raise ValueError(f"noise_estimate_on: mask shape {mask.shape} does not "
                             f"match {self.residuals_real.shape} residuals")
        chi2 = float(np.sum(self.residuals_real[mask] ** 2
                            + self.residuals_imag[mask] ** 2))
        return estimate_noise_percent(chi2, int(np.count_nonzero(mask)))


# Kept for callers of the former separate class: it duplicated every field and
# property of KKResult (doc/KRAMERS_KRONIG_REVIEW.md 2.3).
LinKKResult = KKResult
