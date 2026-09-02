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
- Shape decomposition: the slope of a least-squares line through the
  residuals, and the period and amplitude of the strongest wave left once
  that line is removed.

The first two answer "is this systematic?"; the third says *what kind*, and
says it as two numbers rather than as a choice between them. A trend means an
element is missing or of the wrong type; a wave means the elements are right
but too few. Residuals routinely carry both at once, so both are always
reported with their own significance - the slope with a p-value, the wave
with its periodogram power - and the reader compares them.
doc/WEIGHTING_AND_STATISTICS.md carries the reading for users.

Nothing here changes a fit or a criterion. In particular `n_eff` is reported
but deliberately not fed into AIC/BIC - see its note in `_series_diagnostics`.
"""

import logging
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.signal import lombscargle
from scipy.stats import linregress, norm

from .diagnostics import compute_weights

logger = logging.getLogger(__name__)

# A residual series is called systematic when the runs test rejects
# independence at this level. 0.01 rather than the customary 0.05: with 80-160
# points the test is sensitive enough that 5% flags spectra whose residuals are
# visibly random, and this diagnostic is worth acting on only when it is not
# borderline.
RUNS_P_THRESHOLD = 0.01

# Below this normalized Lomb-Scargle power the peak is not far enough above
# what white noise produces. Over 2000 white-noise series of 80 points the peak
# power has median 0.10, p95 0.16 and max 0.30, so 0.2 still lets 0.9% through
# on its own; every systematic case tested landed at 0.33-0.98. Nothing is
# suppressed by it - the power is always reported - it is the level the reader
# compares the printed power against before acting on a wave.
MIN_PERIODOGRAM_POWER = 0.2

# Shortest period the periodogram is asked about, in decades, before the
# Nyquist floor below raises it. Anything shorter approaches point-to-point
# alternation, which is noise rather than structure, and the runs test already
# covers it.
MIN_PERIOD_DECADES = 0.4

# Periods shorter than this multiple of the sample spacing cannot be resolved,
# only aliased - and the alias lands at the sampling interval, where it reads
# as a short ripple and inverts the advice the CLI gives
# (test_sparse_spectrum_does_not_alias_a_drift_into_a_ripple).
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
    slope : float
        Slope of a least-squares line through the residuals, per decade of
        frequency, in the units of the weighted residuals. Times the window it
        gives the trend's total swing, which is what the reader compares with
        `amplitude` - same units, the trend's full swing against the wave's
        half-swing - to see which of the two shapes is the larger part of this
        residual. That is the question a single verdict used to answer by
        discarding one of them.
    slope_p : float
        p-value of the slope. NaN when the series has no variance. Optimistic
        for correlated residuals - it assumes independent points, which is the
        very thing the runs test is there to check - so read it as a ranking
        of the two parts, not as an exact probability.
    period_decades : float
        Period of the strongest wave left after the line is removed, in
        decades of frequency. NaN only when no period is resolvable at all
        (window narrower than the shortest testable period, or no variance).
        A weak wave is still reported; `power` says how weak.
    amplitude : float
        Amplitude of that wave in the units of the weighted residuals, implied
        by the share of the variance it explains. An estimate, not a fitted
        parameter: whatever the residual carries at other periods counts into
        that variance, so a wave riding on curvature reads a little low
        (measured: 0.84 for a wave of 1.0 on top of a drift).
    power : float
        Normalized Lomb-Scargle power at the peak, in [0, 1] - the wave's
        significance, to be read against MIN_PERIODOGRAM_POWER.
    is_systematic : bool
        Whether the runs test rejects independence at RUNS_P_THRESHOLD.
    """
    lag1_autocorr: float
    n_eff: float
    runs: int
    runs_z: float
    runs_p: float
    slope: float
    slope_p: float
    period_decades: float
    amplitude: float
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
) -> Tuple[float, float, float]:
    """
    Dominant period, power and amplitude of a residual series.

    A Lomb-Scargle periodogram locates the period rather than a fitted
    sinusoid: the residual of a misspecified circuit is generally not a
    sinusoid, so a direct fit would return whatever period the window allows
    together with a confidence interval that looks tight. The amplitude at that
    peak follows from its power, and is the number the reader compares
    against the trend's span.

    Pass the residuals with any linear trend already removed, otherwise the
    peak is the trend.

    Parameters
    ----------
    log_freq : ndarray of float
        log10 of the frequencies, ordered.
    residuals : ndarray of float
        Residuals aligned with `log_freq`.

    Returns
    -------
    period_decades, power, amplitude : float, float, float
        Period at the peak, its normalized power in [0, 1], and the amplitude
        implied by that power. All NaN when the spectrum is too sparse or too
        narrow to resolve any period, or when the series has no variance.
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
        return float('nan'), float('nan'), float('nan')

    # One candidate period per ~2% of the window: fine enough that the peak
    # is located to better than the reporting precision, cheap at this size.
    periods = np.linspace(shortest, window, 500)
    power = lombscargle(log_freq, r, 2 * np.pi / periods, normalize=True)
    peak = int(np.argmax(power))
    period = float(periods[peak])

    # Amplitude from the normalized power rather than a second least-squares
    # fit at the peak: normalized power is the fraction of the variance the
    # sinusoid explains, and a sinusoid of amplitude A has variance A^2/2.
    # Fitting sin and cos directly gives the same answer where the period is
    # comfortably resolved, but the basis is ill-conditioned near the Nyquist
    # floor, where it returned amplitudes several times the residual range
    # (measured: 11.3 on residuals spanning 6). This form cannot exceed
    # sqrt(2 * variance) by construction.
    amplitude = np.sqrt(2.0 * power[peak] * np.mean(r ** 2))
    return period, float(power[peak]), float(amplitude)


def _series_diagnostics(
    log_freq: NDArray[np.float64],
    residuals: NDArray[np.float64]
) -> SeriesDiagnostics:
    """Run all the statistics on one residual series."""
    rho = lag1_autocorrelation(residuals)
    runs, expected, z, p = runs_test(residuals)

    # The line is fitted and removed before the periodogram runs, which is
    # what lets both shapes be reported at once: otherwise a drift and a wave
    # compete for the same peak, the drift wins it at a period near the window
    # width, and any genuine wave underneath it disappears.
    trend = linregress(log_freq, residuals)
    period, power, amplitude = dominant_period(
        log_freq, residuals - (trend.slope * log_freq + trend.intercept))

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

    return SeriesDiagnostics(
        lag1_autocorr=rho,
        n_eff=float(n_eff),
        runs=runs,
        runs_z=z,
        runs_p=p,
        slope=float(trend.slope),
        slope_p=float(trend.pvalue),
        period_decades=period,
        amplitude=amplitude,
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
    >>> round(d.real.period_decades, 1), round(d.real.amplitude, 1)
    (1.5, 3.0)

    >>> analyze_residuals(f, Z, Z * 1.01).is_systematic   # offset has no shape
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
        real=_series_diagnostics(log_freq, residuals_real),
        imag=_series_diagnostics(log_freq, residuals_imag),
        window_decades=window,
    )


__all__ = [
    'ResidualDiagnostics',
    'SeriesDiagnostics',
    'analyze_residuals',
    'lag1_autocorrelation',
    'runs_test',
    'dominant_period',
    'MIN_PERIODOGRAM_POWER',
]
