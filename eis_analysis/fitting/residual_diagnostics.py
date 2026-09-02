"""
Shape tests for fit residuals: are they noise, or is the model missing something?

A circuit fit reports how *large* its residuals are (`fit_error_rel`) but says
nothing about their *shape*, and the two questions have different answers. A
model missing a relaxation can still fit to 2% while its residuals march
smoothly across the spectrum; a correct model fitting to 2% scatters them at
random. Only the second is noise, and only the second justifies the
independence assumption AIC and BIC are built on.

Three statistics, computed separately for the real and imaginary parts:

- Lag-1 autocorrelation. ~0 for independent residuals, approaching 1 when
  consecutive points share a systematic offset.
- Runs test (Wald-Wolfowitz). The number of sign changes against the count
  expected from independent residuals, as a p-value.
- Dominant period from a Lomb-Scargle periodogram over log10(f), reported
  relative to the width of the measured window.

The first two answer "is this systematic?"; only the third says *what kind* -
a short period is a ripple (right elements, too few of them), a period near
the window width is a bump or a drift (an element missing or of the wrong
type). doc/WEIGHTING_AND_STATISTICS.md carries the reading for users.

Nothing here changes a fit or a criterion. In particular `n_eff` is reported
but deliberately not fed into AIC/BIC - see its note in `_series_diagnostics`.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.signal import lombscargle
from scipy.stats import norm

from .diagnostics import compute_weights

logger = logging.getLogger(__name__)

# A residual series is called systematic when the runs test rejects
# independence at this level. 0.01 rather than the customary 0.05: with 80-160
# points the test is sensitive enough that 5% flags spectra whose residuals are
# visibly random, and this diagnostic is worth acting on only when it is not
# borderline.
RUNS_P_THRESHOLD = 0.01

# The periodogram always returns a best period, so one is needed to tell a
# genuine ripple from the longest wavelength the window can hold. A bump or a
# monotone drift maximizes at the window width itself (measured: 5.98 and 7.00
# decades in a 7-decade window), whereas a real 1.5-decade ripple was recovered
# at 1.50. Half the window separates them with room to spare, and a period
# above it is reported as a trend rather than named as a period.
MAX_PERIOD_FRACTION = 0.5

# Below this normalized Lomb-Scargle power the peak is not far enough above
# what white noise produces. Over 2000 white-noise series of 80 points the peak
# power has median 0.10, p95 0.16 and max 0.30, so 0.2 still lets 0.9% through
# on its own; every systematic case tested landed at 0.33-0.98. The residual
# gate is not this threshold alone - the caller also requires the runs test to
# have rejected independence, which white noise does not.
MIN_PERIODOGRAM_POWER = 0.2

# Shortest period the periodogram is asked about, in decades, before the
# Nyquist floor below raises it. Anything shorter approaches point-to-point
# alternation, which is noise rather than structure, and the runs test already
# covers it.
MIN_PERIOD_DECADES = 0.4

# Periods shorter than this multiple of the sample spacing cannot be resolved,
# only aliased. Without the floor a 10-point spectrum over 6 decades (spacing
# 0.67) reports a pure monotone drift as a 0.62-decade "ripple" at power 0.98,
# for every seed tried - the alias sits at the sampling interval, and the CLI
# then prescribes adding a branch when the opposite fix is needed.
NYQUIST_PERIOD_FACTOR = 2.0


@dataclass
class SeriesDiagnostics:
    """
    Shape statistics for one residual series (the real or the imaginary part).

    Attributes
    ----------
    lag1_autocorr : float
        Lag-1 autocorrelation. ~0 for independent residuals; approaches 1 when
        the model carries a systematic offset across neighbouring frequencies.
    n_eff : float
        Effective sample size n*(1-rho)/(1+rho) for an AR(1) process.
        Informational only - see the note in `_series_diagnostics`.
    runs : int
        Observed number of sign runs.
    runs_z, runs_p : float
        Normal-approximation z statistic and two-sided p-value. A large
        negative z means far fewer sign changes than chance allows.
    period_decades : float or None
        Dominant period in decades of frequency, or None when the peak is a
        trend rather than a ripple (see MAX_PERIOD_FRACTION) or too weak to
        distinguish from noise (see MIN_PERIODOGRAM_POWER).
    period_over_window : float
        The same period divided by the width of the measured window, reported
        even when `period_decades` is None so the caller can see why.
    power : float
        Normalized Lomb-Scargle power at the peak, in [0, 1].
    is_systematic : bool
        Whether the runs test rejects independence at RUNS_P_THRESHOLD.
    """
    lag1_autocorr: float
    n_eff: float
    runs: int
    runs_z: float
    runs_p: float
    period_decades: Optional[float]
    period_over_window: float
    power: float
    is_systematic: bool


@dataclass
class ResidualDiagnostics:
    """
    Shape tests for a fit's residuals.

    Attributes
    ----------
    real, imag : SeriesDiagnostics or None
        Per-part statistics. None when the series was too short or too
        degenerate to test (see `warnings`).
    window_decades : float
        Width of the measured frequency window, log10(f_max/f_min).
    warnings : list of str
        Why a part could not be tested, if it could not.
    """
    real: Optional[SeriesDiagnostics] = None
    imag: Optional[SeriesDiagnostics] = None
    window_decades: float = 0.0
    warnings: List[str] = field(default_factory=list)

    @property
    def is_systematic(self) -> bool:
        """True when either part fails the runs test."""
        return any(s is not None and s.is_systematic
                   for s in (self.real, self.imag))


def lag1_autocorrelation(residuals: NDArray[np.float64]) -> float:
    """
    Lag-1 autocorrelation of a residual series.

    Parameters
    ----------
    residuals : ndarray of float
        Residuals ordered by frequency.

    Returns
    -------
    rho1 : float
        Sum(r_i * r_i+1) / Sum(r_i^2) over the centred series; 0.0 when the
        series is constant (no variance to correlate).
    """
    r = np.asarray(residuals, dtype=float)
    r = r - r.mean()
    denominator = float(np.sum(r ** 2))
    if denominator <= 0:
        return 0.0
    return float(np.sum(r[:-1] * r[1:]) / denominator)


def runs_test(residuals: NDArray[np.float64]) -> Tuple[int, float, float, float]:
    """
    Wald-Wolfowitz runs test on the signs of a residual series.

    Counts maximal blocks on one side of the median and compares them with
    the count expected from independent signs. Systematic residuals cross far
    less often than chance, giving a large negative z.

    The median rather than zero is the textbook reference, and here it is also
    the more useful one: a weighted fit does not force its residuals to zero
    mean, so a bump riding on an offset would keep one sign throughout and
    escape a zero-referenced test. On measured fits the two references agree
    closely (p = 1e-16 vs 8e-17 on an under-parametrized circuit, 0.89 vs 0.37
    on the correct one); they part company exactly on that bump, which the
    median flags at p = 1e-11 and zero misses at p = 0.15.

    What the median cannot see is a uniform bias - shifting every residual by
    a constant leaves the crossings where they were. Neither can
    `lag1_autocorrelation`, which centres the series first. That is by design:
    a constant offset has no shape, and its size is already what
    `fit_error_rel` reports. These three statistics exist for the opposite
    case, where the error is small but structured.

    scipy has no runs test (statsmodels does, and is not a dependency here),
    so the normal approximation is computed directly. It is accurate for the
    tens-to-hundreds of points an EIS spectrum has.

    Parameters
    ----------
    residuals : ndarray of float
        Residuals ordered by frequency. Values exactly at the median are
        dropped: they belong to neither side, and assigning them one would
        invent a run boundary.

    Returns
    -------
    runs, expected, z, p : int, float, float, float
        Observed runs, expected runs, z statistic and two-sided p-value.
        z and p are NaN when the test is not defined - fewer than three
        non-zero residuals, or all of one sign, which leaves the variance
        zero or negative.
    """
    r = np.asarray(residuals, dtype=float)
    signs = np.sign(r - np.median(r))
    signs = signs[signs != 0]

    n_pos = int(np.sum(signs > 0))
    n_neg = int(np.sum(signs < 0))
    n = n_pos + n_neg

    runs = 1 + int(np.sum(signs[1:] != signs[:-1])) if n > 0 else 0

    if n < 3 or n_pos == 0 or n_neg == 0:
        return runs, float('nan'), float('nan'), float('nan')

    expected = 2.0 * n_pos * n_neg / n + 1.0
    variance = (2.0 * n_pos * n_neg * (2.0 * n_pos * n_neg - n)
                / (n ** 2 * (n - 1.0)))
    if variance <= 0:
        return runs, expected, float('nan'), float('nan')

    z = (runs - expected) / np.sqrt(variance)
    p = 2.0 * norm.cdf(-abs(z))
    return runs, float(expected), float(z), float(p)


def dominant_period(
    log_freq: NDArray[np.float64],
    residuals: NDArray[np.float64]
) -> Tuple[float, float]:
    """
    Dominant period of a residual series in decades of frequency.

    A Lomb-Scargle periodogram rather than a fitted sinusoid: the residual of
    a misspecified circuit is generally not a sinusoid, so the fit would
    return a period equal to the window width together with a confidence
    interval that looks tight. The periodogram reports the same peak without
    dressing it as an estimate, and the caller decides whether the peak is a
    ripple or a trend from `period / window`.

    Parameters
    ----------
    log_freq : ndarray of float
        log10 of the frequencies, ordered.
    residuals : ndarray of float
        Residuals aligned with `log_freq`.

    Returns
    -------
    period_decades, power : float, float
        Period at the peak and its normalized power in [0, 1]. Both NaN when
        the spectrum is too sparse or too narrow to resolve any period, or
        when the series has no variance.
    """
    r = np.asarray(residuals, dtype=float)
    r = r - r.mean()
    window = float(log_freq[-1] - log_freq[0])

    # Shortest resolvable period: the Nyquist limit of the actual sampling,
    # never below MIN_PERIOD_DECADES. Median spacing rather than mean - a
    # spectrum with one large gap is still sampled at its usual rate.
    spacing = float(np.median(np.diff(log_freq)))
    shortest = max(MIN_PERIOD_DECADES, NYQUIST_PERIOD_FACTOR * spacing)

    if window <= shortest or np.sum(r ** 2) <= 0:
        return float('nan'), float('nan')

    # One candidate period per ~2% of the window: fine enough that the peak
    # is located to better than the reporting precision, cheap at this size.
    periods = np.linspace(shortest, window, 500)
    power = lombscargle(log_freq, r, 2 * np.pi / periods, normalize=True)
    peak = int(np.argmax(power))
    return float(periods[peak]), float(power[peak])


def _series_diagnostics(
    log_freq: NDArray[np.float64],
    residuals: NDArray[np.float64],
    window: float
) -> SeriesDiagnostics:
    """Run all three statistics on one residual series."""
    rho = lag1_autocorrelation(residuals)
    runs, expected, z, p = runs_test(residuals)
    period, power = dominant_period(log_freq, residuals)

    n = len(residuals)
    # Effective sample size for an AR(1) process. Reported because it puts a
    # number on how much independent information the spectrum really carries,
    # but deliberately NOT substituted for n in AIC/BIC. Three reasons: it is
    # computed from the residuals, so it differs per candidate and their BICs
    # would sit on different scales; it does nothing when the model is
    # adequate (rho ~ 0), which is when model selection happens; and where it
    # does bite it degenerates - a measured rho of 0.987 gives n_eff = 1.0,
    # i.e. ln(n_eff) = 0 and no complexity penalty at all. Correlated EIS
    # residuals are model error, not AR(1) noise, and this diagnostic exists
    # to expose that rather than absorb it.
    n_eff = n * (1.0 - rho) / (1.0 + rho) if rho > -1.0 else float(n)

    fraction = period / window if window > 0 and np.isfinite(period) else float('nan')
    is_ripple = (np.isfinite(fraction)
                 and fraction < MAX_PERIOD_FRACTION
                 and power >= MIN_PERIODOGRAM_POWER)

    return SeriesDiagnostics(
        lag1_autocorr=rho,
        n_eff=float(n_eff),
        runs=runs,
        runs_z=z,
        runs_p=p,
        period_decades=period if is_ripple else None,
        period_over_window=fraction,
        power=power,
        is_systematic=bool(np.isfinite(p) and p < RUNS_P_THRESHOLD),
    )


def analyze_residuals(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complex128],
    Z_fit: NDArray[np.complex128],
    weighting: str = 'modulus'
) -> ResidualDiagnostics:
    """
    Test whether a fit's residuals look like noise.

    Parameters
    ----------
    frequencies : ndarray of float
        Frequencies [Hz]. Sorted internally - the statistics read consecutive
        points as neighbours, so an unsorted spectrum would be tested in
        whatever order the instrument wrote it.
    Z : ndarray of complex
        Measured impedance [Ohm].
    Z_fit : ndarray of complex
        Fitted impedance at the same frequencies.
    weighting : str, optional
        Weighting used for the fit (default 'modulus'), applied so the tested
        residuals are the ones the optimizer actually minimized.

    Returns
    -------
    diagnostics : ResidualDiagnostics
        Per-part statistics, or empty parts with a warning when the spectrum
        is too short to test.

    Notes
    -----
    The real and imaginary parts are tested separately. The optimizer keeps
    its residuals as one concatenated vector of length 2n
    (`fitting/circuit.py`), which must not be used here: the join between the
    two halves is not a frequency step, and the two parts carry different
    systematic structure anyway. This matches how residuals are held
    elsewhere - `KKResult` keeps `residuals_real` and `residuals_imag` apart,
    and so does the residual panel in `visualization/plots.py`.

    Examples
    --------
    >>> import numpy as np
    >>> f = np.logspace(-2, 5, 80)
    >>> Z = np.full(80, 100 + 0j)
    >>> ripple = 3 * np.sin(2 * np.pi * np.log10(f) / 1.5)
    >>> d = analyze_residuals(f, Z, Z - ripple, weighting='uniform')
    >>> d.is_systematic
    True
    >>> round(d.real.period_decades, 1)
    1.5

    A constant offset is deliberately *not* flagged - it has no shape, and its
    size is what `fit_error_rel` already reports:

    >>> analyze_residuals(f, Z, Z * 1.01).is_systematic
    False
    """
    frequencies = np.asarray(frequencies, dtype=float)
    Z = np.asarray(Z, dtype=complex)
    Z_fit = np.asarray(Z_fit, dtype=complex)

    if len(frequencies) < 4:
        return ResidualDiagnostics(warnings=[
            f"Residual diagnostics need at least 4 points, got {len(frequencies)}"
        ])
    if not np.all(frequencies > 0):
        return ResidualDiagnostics(warnings=[
            "Residual diagnostics need positive frequencies (log10 scale)"
        ])

    # sort_by_frequency takes a pair; the fit has to follow the same
    # permutation, so the order is applied to all three directly.
    order = np.argsort(frequencies)
    frequencies, Z, Z_fit = frequencies[order], Z[order], Z_fit[order]

    weights = compute_weights(Z, weighting)
    residuals_real = weights * (Z.real - Z_fit.real)
    residuals_imag = weights * (Z.imag - Z_fit.imag)

    log_freq = np.log10(frequencies)
    window = float(log_freq[-1] - log_freq[0])

    return ResidualDiagnostics(
        real=_series_diagnostics(log_freq, residuals_real, window),
        imag=_series_diagnostics(log_freq, residuals_imag, window),
        window_decades=window,
    )


__all__ = [
    'ResidualDiagnostics',
    'SeriesDiagnostics',
    'analyze_residuals',
    'lag1_autocorrelation',
    'runs_test',
    'dominant_period',
    'RUNS_P_THRESHOLD',
    'MAX_PERIOD_FRACTION',
    'MIN_PERIODOGRAM_POWER',
    'NYQUIST_PERIOD_FACTOR',
]
