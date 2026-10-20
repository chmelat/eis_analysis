"""
Console and file output of the circuit fit for the EIS CLI.

- _log_fit_result: parameters, CIs, significance, fit-to-noise ratio, residuals
- _log_diffevo_diagnostics / _log_multistart_diagnostics: optimizer reports
- _save_fit_data / _save_fit_figure: fit JSON/CSV and figure under --save
"""

import argparse
import logging
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import draw_figure
from ...fitting import (
    FitResult,
    MultistartResult,
    DiffEvoResult,
    arc_time_constants,
    arc_resistance_sums,
)
from ...fitting.diagnostics import NoiseAdequacyResult
from ...fitting.residual_diagnostics import (
    MIN_PERIODOGRAM_POWER,
    RUNS_P_THRESHOLD,
    ResidualDiagnostics,
    analyze_residuals,
)
from ...fitting.config import (FIT_QUALITY_EXCELLENT_ERROR, FIT_QUALITY_GOOD_ERROR,
                               FIT_QUALITY_ACCEPTABLE_ERROR, SIGNIFICANCE_NEGLIGIBLE)
from ...io import fit_result_record, save_fit_result
from ...visualization import plot_circuit_fit
from .fit_noise import log_noise_adequacy

logger = logging.getLogger(__name__)


# The weight each --weighting puts on a point, for the log lines
_WEIGHTING_LABELS = {
    'uniform': 'uniform (w=1)',
    'sqrt': 'sqrt (w=1/sqrt|Z|)',
    'modulus': 'modulus (w=1/|Z|)',
    'proportional': 'proportional (w=1/|Z|^2)',
}


def _mark_zhit_fit(fig: plt.Figure, args: argparse.Namespace) -> None:
    """
    Note on the fit figure that the fit ran on the Z-HIT reconstruction.

    The "Data" curve is whatever impedance array the fit was handed, so under
    --fit-on it would otherwise pass the reconstruction off as the measurement.
    A footer rather than a suptitle: the panel titles sit at the top.
    """
    # getattr, not args.fit_on: run_circuit_fitting is exported, and a caller
    # building the namespace by hand predates --fit-on. This call sits inside
    # draw_figure's broad except, so an AttributeError there would silently
    # cost the figure of a fit that had already succeeded.
    if getattr(args, 'fit_on', 'original') == 'original':
        return
    fig.text(0.5, 0.005,
             'The "Data" curve is the Z-HIT reconstruction of |Z| from the '
             'phase, not the raw measurement',
             ha='center', va='bottom', fontsize=8, style='italic', color='0.35')


def _save_fit_figure(
    frequencies: NDArray,
    Z: NDArray,
    result: FitResult,
    args: argparse.Namespace,
    save_suffix: str,
    title: Optional[str] = None
) -> Optional[plt.Figure]:
    """Draw and save the fit figure; a figure error costs the figure, not the fit."""
    def draw() -> plt.Figure:
        fig = plot_circuit_fit(frequencies, Z, result, title=title)
        _mark_zhit_fit(fig, args)
        return fig

    return draw_figure(draw, 'Fit', args.save, save_suffix, args.format)


def _save_fit_data(
    frequencies: NDArray,
    Z: NDArray,
    Z_fit: NDArray,
    result: FitResult,
    args: argparse.Namespace,
    save_suffix: str
) -> None:
    """Write the fit as JSON and CSV under --save; an export error costs the files, not the fit."""
    if args.save is None:
        return
    try:
        # getattr: a caller building the namespace by hand need not have
        # these options, and a missing one must not cost the fit
        context = {
            'input': args.input,
            'optimizer': 'linear' if args.voigt_chain else args.optimizer,
            'fit_on': getattr(args, 'fit_on', 'original'),
            'f_min': getattr(args, 'f_min', None),
            'f_max': getattr(args, 'f_max', None),
            'exclude_flagged': getattr(args, 'exclude_flagged', None),
            'excluded_frequencies': getattr(args, 'excluded_frequencies', None),
        }
        record = fit_result_record(frequencies, Z, Z_fit, result, args.weighting, context)
        paths = save_fit_result(f"{args.save}_{save_suffix}", record,
                                frequencies, Z, Z_fit)
    except Exception as e:
        logger.warning(f"Fit result export failed: {e}")
        logger.debug("Traceback:", exc_info=True)
        return
    for path in paths:
        logger.info(f"Saved: {path}")


# =============================================================================
# Fitting Diagnostics Logging
# =============================================================================

def _log_diffevo_diagnostics(diffevo_result: DiffEvoResult) -> None:
    """Log Differential Evolution optimization diagnostics."""
    diag = diffevo_result.diagnostics
    if diag is None:
        return

    log_separator()
    logger.info("Differential Evolution optimization")
    log_separator()

    # Settings - map strategy name to option number
    strategy_to_option = {'randtobest1bin': 1, 'best1bin': 2, 'rand1bin': 3}
    option_num = strategy_to_option.get(diag.strategy, 1)
    logger.info(f"  Strategy: {diag.strategy} (option {option_num})")
    logger.info(f"  Population: {diag.popsize} * n_params")
    logger.info(f"  Max iterations: {diag.maxiter}")
    logger.info(f"  Tolerance: {diag.tol}")
    logger.info(f"  Workers: {diag.workers}")
    logger.info(f"  Weighting: {_WEIGHTING_LABELS.get(diag.weighting, diag.weighting)}")

    # Get params info from best_result
    result = diffevo_result.best_result
    n_free = len(result.params_opt) - diag.n_fixed_params
    logger.info(f"  Parameters: {n_free} (free)")
    if diag.log_search_params:
        logger.info(f"  Search space: log10 for {', '.join(diag.log_search_params)} "
                    "(bounds span many decades)")

    # Initial guess captured before DE ran (circuit.update_params overwrote
    # circuit.get_all_params() with the final fit, so we read the snapshot
    # stored in diagnostics instead).
    logger.info(f"  Initial guess: {np.array(diag.initial_guess)}")
    logger.info("")

    # DE progress
    logger.info("Running differential evolution...")
    logger.info(f"  DE converged: {diag.de_converged}")
    logger.info(f"  DE iterations: {diag.de_iterations}")
    logger.info(f"  DE evaluations: {diag.de_evaluations}")
    logger.info(f"  DE error: {diag.de_error:.3f}%")
    logger.info("")

    # Refinement
    jac_type = "analytic" if diag.jacobian_type == "analytic" else "numeric"
    logger.info(f"Refining with least_squares ({jac_type} Jacobian)...")
    logger.info(f"  Refined error: {diag.refined_error:.3f}%")
    # Improvement is measured on the optimized objective (weighted SSR), which
    # is what the DE/refinement selection uses. The error above is a weighted
    # relative error and can move the other way; reporting it as "improvement"
    # would contradict the choice the optimizer actually made.
    logger.info(f"  Improvement (SSR): {diffevo_result.improvement:+.3f}%")
    if not diag.archive_checked:
        logger.info("  Archive check: off (--no-archive-check)")
    else:
        outcome = "the result is one of them" if diag.archive_used else "none fitted better"
        logger.info(f"  Archive check: {diag.archive_candidates} candidate(s) from early "
                    f"generations refined, {outcome}")
    logger.info("")

    # Summary
    log_separator()
    logger.info("Differential Evolution results")
    log_separator()
    logger.info(f"  Strategy: {diag.strategy}")
    logger.info(f"  Total evaluations: {diag.total_evaluations}")
    logger.info(f"  DE error: {diag.de_error:.3f}% -> Refined: {diag.refined_error:.3f}%")
    logger.info("")


def _log_multistart_diagnostics(multistart_result: MultistartResult) -> None:
    """Log Multi-start optimization diagnostics."""
    diag = multistart_result.diagnostics
    if diag is None:
        return

    log_separator()
    logger.info("Multi-start optimization")
    log_separator()

    # Settings
    logger.info(f"  Restarts: {diag.n_restarts}")
    logger.info(f"  Perturbation scale: {diag.scale} sigma")
    logger.info(f"  Perturbation method: {diag.perturbation_method}")
    logger.info(f"  Weighting: {diag.weighting}")
    logger.info(f"  Jacobian: {diag.jacobian_type}")
    logger.info(f"  Parallel: {diag.parallel}")
    logger.info("")

    # Results
    logger.info(f"Running {diag.n_restarts} optimization starts...")
    logger.info(f"  Successful: {diag.n_successful}/{diag.n_restarts}")
    logger.info(f"  Initial error: {diag.initial_error:.3f}%")
    logger.info(f"  Best error: {diag.best_error:.3f}%")
    logger.info(f"  Best start: #{diag.best_start_index}")
    # On the weighted SSR the best start is selected on, as for DE; the errors
    # above are weighted relative errors and can move the other way.
    logger.info(f"  Improvement (SSR): {multistart_result.improvement:+.1f}%")
    logger.info("")

    # Warnings
    for warning in diag.warnings:
        logger.warning(f"  {warning}")


def _fmt_pair(real: float, imag: float, fmt: str) -> str:
    """Format one statistic for both parts as 'Re / Im'."""
    return " / ".join("n/a" if not np.isfinite(v) else format(v, fmt)
                      for v in (real, imag))


def _residual_diagnostics(
    frequencies: NDArray[np.float64],
    Z: NDArray[np.complexfloating],
    Z_fit: NDArray[np.complexfloating],
    weighting: str
) -> Optional[ResidualDiagnostics]:
    """Residual shape tests, or None if they fail.

    Catches its own errors: in _fit_standard_circuit an exception would land
    in the handler's try/except, which reports a fitting error and returns
    no result, and _fit_voigt_chain calls it outside any try. A diagnostic
    must never cost the user a converged fit.
    """
    try:
        return analyze_residuals(frequencies, Z, Z_fit, weighting)
    except Exception as exc:
        logger.debug(f"Residual diagnostics skipped: {exc}", exc_info=True)
        return None


def _log_residual_diagnostics(d: Optional[ResidualDiagnostics],
                              suggest_exclude: bool = True) -> None:
    """Log the residual shape tests. Pure formatter over ResidualDiagnostics.

    `suggest_exclude` False when --exclude-flagged already ran: naming it
    again would send the reader in a loop instead of to the model.
    """
    if d is None:
        return
    if d.real is None or d.imag is None:
        for warning in d.warnings:
            logger.warning(f"  {warning}")
        return

    logger.info(
        f"  Residuals: rho1 = "
        f"{_fmt_pair(d.real.lag1_autocorr, d.imag.lag1_autocorr, '+.2f')} (Re/Im), "
        f"runs p = {_fmt_pair(d.real.runs_p, d.imag.runs_p, '.1e')}")

    if not d.is_systematic:
        return

    # Both shapes, always, with their own significance - never a choice
    # between them. Residuals commonly carry a trend and leftover structure at
    # once, and naming only the stronger one prescribes half the repair: the
    # trend says an element is missing or of the wrong type, the structure
    # says the elements are right but too few. The reader compares the trend's
    # span with the structure's amplitude (same units) to see which dominates,
    # and the p-value and power to see whether each is there at all.
    logger.warning(f"  Residuals are not random "
                   f"(Re/Im, {d.window_decades:.1f} decade window):")
    # The runs test is what raised the warning; trend and structure below only
    # describe the shape. Cold agent T14: with both insignificant and no
    # trigger named, the warning read as contradicting itself.
    logger.warning(
        f"    runs test: p = {_fmt_pair(d.real.runs_p, d.imag.runs_p, '.1e')} "
        f"(fails below {RUNS_P_THRESHOLD})")
    span_re, span_im = (abs(s.slope) * d.window_decades for s in (d.real, d.imag))
    logger.warning(
        f"    trend:     {_fmt_pair(d.real.slope, d.imag.slope, '+.2g')} per decade, "
        f"span {_fmt_pair(span_re, span_im, '.2g')} "
        f"(p = {_fmt_pair(d.real.slope_p, d.imag.slope_p, '.1e')})")
    # No period on this line - see `residual_structure` for why.
    logger.warning(
        f"    structure: amplitude {_fmt_pair(d.real.amplitude, d.imag.amplitude, '.2g')} "
        f"(power {_fmt_pair(d.real.power, d.imag.power, '.2f')}, "
        f"noise < {MIN_PERIODOGRAM_POWER})")
    logger.warning("    A trend means an element is missing or of the wrong type. "
                   "Structure left after")
    logger.warning("    it means the right elements, too few of them - the residual "
                   "plot shows where.")
    # Cold agent T13: 5 bad points among 81 failed the runs test of the right
    # circuit, and the agent added an element. The per-point check lists them.
    if suggest_exclude:
        logger.warning("    Isolated bad points fail it too: if the per-point residual "
                       "check lists points,")
        logger.warning("    refit without them (--exclude-flagged kk) before changing "
                       "the model.")


def _significance_suffix(significance: Optional[NDArray[np.float64]], i: int) -> str:
    """Format one parameter's significance for the parameter line, '' if absent."""
    if significance is None or i >= len(significance):
        return ""
    if significance[i] < SIGNIFICANCE_NEGLIGIBLE:
        return f"  S={significance[i]:.3f} [negligible - element may be omitted]"
    return f"  S={significance[i]:.2f}"


def _log_fit_result(result: FitResult,
                    residuals: Optional[ResidualDiagnostics] = None,
                    noise: Optional[NoiseAdequacyResult] = None,
                    fit_on: str = 'original',
                    noise_note: Optional[str] = None,
                    noise_source: Optional[str] = None,
                    suggest_exclude: bool = True,
                    n_excluded: int = 0,
                    n_flagged_kept: int = 0) -> None:
    """Log fit result with parameters, CIs, fit-to-noise ratio and residual shape.

    `n_excluded` is the number of points --exclude-flagged left out,
    `n_flagged_kept` the number of KK-flagged points a fit without it keeps;
    with either, a note says which fit to report.
    """
    logger.info("")
    logger.info("Fit results:")
    # The labels below count each parameter type from the left (K's R is an
    # R too); the fitted expression maps them to the elements directly
    logger.info(f"  Fitted: {result.circuit}")
    if result.block_order is not None:
        # Without this line a result in another order than the one written
        # would look like a swap
        logger.info("  Equal blocks are reported in ascending time constant (fastest first)")
    logger.info("  Parameters:")

    # Get param labels
    labels = result.param_labels
    if labels is None:
        labels = [f"p{i}" for i in range(len(result.params_opt))]

    # Get 95% confidence intervals
    ci_low, ci_high = result.params_ci_95

    # Print each parameter. If a parameter sits at a bound or is fixed, the
    # Jacobian-based CI is not meaningful (locally non-quadratic surface), so
    # we suppress it and tag the line instead.
    bound_status = result.bound_status or [''] * len(result.params_opt)
    significance = result.params_significance
    split = result.series_split
    for i, (label, val, stderr) in enumerate(zip(labels, result.params_opt, result.params_stderr)):
        status = bound_status[i] if i < len(bound_status) else ''
        # Significance is a separate question from the CI - it stays on the
        # line even where the CI is suppressed, since "at a bound" and
        # "irrelevant to the fit" are exactly the pair worth telling apart.
        sig = _significance_suffix(significance, i)
        if status == 'lower' or status == 'upper':
            logger.info(f"    {label:5s} = {val:.2e}  [near {status} bound — CI not meaningful]{sig}")
        elif status == 'fixed':
            logger.info(f"    {label:5s} = {val:.2e}  [fixed]{sig}")
        elif np.isinf(stderr) or np.isnan(stderr):
            logger.info(f"    {label:5s} = {val:.2e} +/- inf{sig}")
        else:
            low, high = ci_low[i], ci_high[i]
            # The stderr is the local one, the CI spans the flat branch
            widened = (" widened, see warnings" if split is not None and split.widens(0.95)
                       and i in (split.index, split.arc_index) else "")
            logger.info(f"    {label:5s} = {val:.2e} +/- {stderr:.2e}  "
                        f"[95% CI{widened}: {low:.2e}, {high:.2e}]{sig}")

    # The time constant is what a reader of an arc asks for first, and its
    # error needs the covariance (R and C are correlated)
    arcs = arc_time_constants(result)
    if arcs:
        logger.info("  Time constants:")
        for arc in arcs:
            err = (f" +/- {arc.tau_stderr:.2e}" if arc.tau_stderr is not None
                   else f" (no stderr: {arc.stderr_note})")
            logger.info(f"    {arc.elements:11s} tau = {arc.tau:.2e} s{err}  "
                        f"(f = {arc.frequency:.2e} Hz)")
    # The sum of the arc resistances in series: with close arcs R1 and R2 are
    # correlated near -1, each is poorly determined and the sum is not - the
    # number cold agents (T13, T14) computed by hand from the covariance
    sums = arc_resistance_sums(result)
    if sums:
        logger.info("  Arc resistances in series:")
    for s in sums:
        if s.stderr is None:
            logger.info(f"    {s.elements:11s} R = {s.R:.4e} Ohm (no stderr: {s.stderr_note})")
            continue
        notes = ([f"corr {s.correlation:+.2f}"] if s.correlation is not None else []) + \
            (["the split is less certain than the sum"] if s.split_less_certain else [])
        note = f"  ({': '.join(notes)})" if notes else ""
        logger.info(f"    {s.elements:11s} R = {s.R:.4e} Ohm +/- {s.stderr:.2e}{note}")

    # Fit quality
    logger.info(f"  Fit error: {result.fit_error_rel:.2f}% (rel), {result.fit_error_abs:.2f} Ohm (abs)")
    # Threshold bounds match the tiers in compute_fit_metrics.
    quality_thresholds = {
        'excellent': f"<{FIT_QUALITY_EXCELLENT_ERROR:.1f}%",
        'good': f"<{FIT_QUALITY_GOOD_ERROR:.1f}%",
        'acceptable': f"<{FIT_QUALITY_ACCEPTABLE_ERROR:.1f}%",
        'poor': f">={FIT_QUALITY_ACCEPTABLE_ERROR:.1f}%",
    }
    threshold = quality_thresholds.get(result.quality)
    logger.info(f"  Quality: {result.quality.capitalize()}"
                + (f" ({threshold})" if threshold else ""))
    log_noise_adequacy(noise, fit_on, noise_note, noise_source)

    if residuals is not None:
        _log_residual_diagnostics(residuals, suggest_exclude)
        # Here, not in _log_residual_diagnostics: --voigt-chain shares that
        # block but prints no uncertainties. s^2 (J^T J)^-1 assumes independent
        # residuals from a correct model; a missing element breaks both.
        if residuals.is_systematic:
            logger.warning("    The +/- and CIs above assume random residuals; with "
                           "these they are too narrow.")

    # Warnings
    for warning in result.all_warnings:
        logger.warning(f"  {warning}")

    # Cold agents T13, T14: with a fit with and one without the points, both
    # quoted the shift between them as the uncertainty of the split. The README
    # says otherwise, but only text at the point of decision was acted on.
    if n_excluded:
        logger.info(f"  Without the {n_excluded} point(s) --exclude-flagged left out: report "
                    "this fit with its +/-")
        logger.info("  and those frequencies. The shift against a fit with all points is "
                    "their pull,")
        logger.info("  not an uncertainty to quote as a range; a shift beyond the +/- means "
                    "points the")
        logger.info("  check did not list may still pull.")
    # Cold-agent A/B (T13, T14, 12 runs): all computed the fit without the
    # points, but 7 reported this one or a range of both - the note above
    # sits only under the other fit
    if n_flagged_kept:
        logger.warning(f"  This fit includes the {n_flagged_kept} point(s) Lin-KK flagged "
                       "above the lowest decade.")
        logger.warning("  Report the fit without them (--exclude-flagged kk) with its +/-, "
                       "not this one")
        logger.warning("  and not a range of both.")

    log_separator()
