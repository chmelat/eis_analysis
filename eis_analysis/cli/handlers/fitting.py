"""
Equivalent circuit fitting handlers for the EIS CLI.

- run_circuit_fitting: dispatch to Voigt-chain or standard circuit fitting
- _fit_standard_circuit: DE, multi-start or single fit of a --circuit
- _compare_circuits: several --circuit candidates ranked by AIC/BIC

Output (console, JSON/CSV, figure) is in fit_output, the Voigt chain in
voigt_fit.
"""

import argparse
import logging
from typing import Optional, Tuple

import matplotlib.pyplot as plt
from numpy.typing import NDArray

from ..logging import log_separator
from ..utils import EISAnalysisError, parse_circuit_expression
from ...fitting import (
    fit_equivalent_circuit,
    fit_circuit_multistart,
    fit_circuit_diffevo,
    FitResult,
)
from .model_comparison import score_candidates, log_comparison
from .fit_noise import noise_adequacy
from .fit_output import (
    _log_diffevo_diagnostics,
    _log_fit_result,
    _log_multistart_diagnostics,
    _residual_diagnostics,
    _save_fit_data,
    _save_fit_figure,
)
from .voigt_fit import _fit_voigt_chain

logger = logging.getLogger(__name__)


# =============================================================================
# Circuit Fitting
# =============================================================================

def run_circuit_fitting(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    noise_percent: Optional[float] = None,
    noise_note: Optional[str] = None
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Run equivalent circuit fitting.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments (uses: no_fit, input, circuit, voigt_chain, weighting,
                       optimizer, multistart, de_*, fit_on, save, format;
                       optional: f_min, f_max, exclude_flagged, and
                       excluded_frequencies - the points select_excluded
                       chose, set on args by the pipeline (eis.py) for the
                       fit JSON)
    noise_percent : float, optional
        Lin-KK noise of these points [%] (kk_noise_for_fit). Without it the
        fit is not compared with the noise. Standard circuits only: a Voigt
        chain is a Lin-KK model itself, its ratio is ~1 by construction.
    noise_note : str, optional
        With noise_percent None: why (e.g. Lin-KK failed on a capacitive end),
        printed in place of the comparison. With noise_percent: the Lin-KK
        run it came from, if not the plain one (e.g. 'Lin-KK with the series C').

    Returns
    -------
    result : FitResult or None
        Fitting result
    fig : Figure or None
        Fit figure
    """
    # Check if fitting should be skipped
    if args.no_fit:
        return None, None

    # Auto-skip for custom data without circuit specification
    if args.input is not None and args.circuit is None and not args.voigt_chain:
        logger.warning("Fit skipped (custom data without --circuit or --voigt-chain)")
        logger.info("For fitting, use: --circuit 'R(100) - (R(1000) | C(1e-6))' "
                    "or --voigt-chain")
        return None, None

    # Voigt chain stays a separate branch: its K(R, tau) elements sit on a
    # fixed tau grid that is not optimized, so its parameter count is not
    # comparable with a nonlinearly fitted circuit's (see the plan note in
    # doc/ZSCOPE_COMPARISON.md).
    if args.voigt_chain:
        return _fit_voigt_chain(frequencies, Z, args)

    # --circuit is action='append', so argparse always yields a list. A caller
    # building the namespace by hand may still pass a bare string; iterating
    # that would silently slice it into single characters, so normalize.
    circuit_exprs = args.circuit
    if isinstance(circuit_exprs, str):
        circuit_exprs = [circuit_exprs]

    # Default circuit for synthetic data: Rs - (R||Q) - (R||Q)
    if args.input is None and circuit_exprs is None:
        circuit_exprs = ['R() - (R() | Q()) - (R() | Q())']

    # A single candidate keeps the original path: one fit, one figure named
    # exactly as before, no comparison table.
    if len(circuit_exprs) == 1:
        return _fit_standard_circuit(frequencies, Z, args, circuit_exprs[0],
                                     noise_percent=noise_percent, noise_note=noise_note)

    return _compare_circuits(frequencies, Z, args, circuit_exprs,
                             noise_percent=noise_percent, noise_note=noise_note)


def _compare_circuits(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    circuit_exprs: list,
    noise_percent: Optional[float] = None,
    noise_note: Optional[str] = None
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Fit every candidate circuit and return the one BIC selects.

    A failing candidate is reported in the table and skipped rather than
    taking the whole run with it.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments
    circuit_exprs : list of str
        Candidate circuit expressions, in the order given on the command line
    noise_percent : float, optional
        See run_circuit_fitting

    Returns
    -------
    result : FitResult or None
        Fit of the BIC-selected circuit
    fig : Figure or None
        Figure of the BIC-selected circuit
    """
    if not circuit_exprs:
        logger.warning("No circuit candidates to compare")
        return None, None

    candidates = []
    figures = []
    for index, expr in enumerate(circuit_exprs, start=1):
        try:
            # Figures are saved per candidate, otherwise they overwrite one file
            result, fig = _fit_standard_circuit(
                frequencies, Z, args, expr, save_suffix=f'fit_{index}',
                noise_percent=noise_percent, noise_note=noise_note
            )
        except EISAnalysisError as e:
            # An unparsable expression must cost its own candidate, not the
            # whole run - the other fits are still valid and worth reporting.
            logger.error(f"Candidate {index} ({expr}): {e}")
            result, fig = None, None
        candidates.append((expr, result))
        figures.append(fig)

    scores = score_candidates(frequencies, Z, candidates, args.weighting,
                              noise_percent=noise_percent)
    log_comparison(scores, len(frequencies), args.weighting)

    if scores[0].result is None:
        logger.error("All candidate circuits failed to fit")
        return None, None

    return scores[0].result, figures[scores[0].candidate - 1]


def _fit_standard_circuit(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    circuit_expr: str,
    save_suffix: str = 'fit',
    noise_percent: Optional[float] = None,
    noise_note: Optional[str] = None
) -> Tuple[Optional[FitResult], Optional[plt.Figure]]:
    """
    Fit using standard circuit with nonlinear optimization.

    Parameters
    ----------
    frequencies : ndarray
        Frequency array [Hz]
    Z : ndarray
        Complex impedance [Ohm]
    args : argparse.Namespace
        CLI arguments for fitting options
    circuit_expr : str
        Circuit expression string
    save_suffix : str, optional
        Filename suffix for --save (default: 'fit'). Comparison runs pass a
        per-candidate suffix so the figures do not overwrite one file.
    noise_percent : float, optional
        See run_circuit_fitting

    Returns
    -------
    result : FitResult or None
        Fitting result
    fig : Figure or None
        Fit figure
    """
    try:
        circuit = parse_circuit_expression(circuit_expr)
    except ValueError as e:
        raise EISAnalysisError(
            f"Circuit parsing error: {e}\n"
            "Use syntax like: --circuit 'R(100) - (R(5000) | C(1e-6))'"
        ) from e

    # Log circuit expression
    log_separator()
    logger.info("Equivalent circuit")
    log_separator()
    logger.info(f"Circuit: {circuit_expr}")

    try:
        if args.optimizer == 'de':
            # Differential Evolution global optimization
            diffevo_result, Z_fit = fit_circuit_diffevo(
                frequencies, Z,
                circuit,
                strategy=args.de_strategy,
                popsize=args.de_popsize,
                maxiter=args.de_maxiter,
                tol=args.de_tol,
                workers=args.de_workers,
                weighting=args.weighting,
                archive_check=not args.no_archive_check,
                use_analytic_jacobian=not args.numeric_jacobian,
                seed=args.seed
            )
            result = diffevo_result.best_result

            # Log DE diagnostics
            _log_diffevo_diagnostics(diffevo_result)

        elif args.optimizer == 'multistart':
            # Multi-start optimization; parser guarantees multistart is None
            # (use default) or a positive int
            n_restarts = args.multistart if args.multistart else 16
            multistart_result, Z_fit = fit_circuit_multistart(
                frequencies, Z,
                circuit,
                n_restarts=n_restarts,
                scale=args.multistart_scale,
                weighting=args.weighting,
                use_analytic_jacobian=not args.numeric_jacobian,
                rng=args.seed
            )
            result = multistart_result.best_result

            # Log multistart diagnostics
            _log_multistart_diagnostics(multistart_result)

        else:
            # Single fit (args.optimizer == 'single')
            result, Z_fit = fit_equivalent_circuit(
                frequencies, Z,
                circuit,
                weighting=args.weighting,
                use_analytic_jacobian=not args.numeric_jacobian
            )

        # Log fit results. Z_fit is bound by every optimizer branch above.
        _log_fit_result(result,
                        _residual_diagnostics(frequencies, Z, Z_fit, args.weighting),
                        noise_adequacy(Z, Z_fit, noise_percent),
                        getattr(args, 'fit_on', 'original'),
                        # noise_note is a reason without a noise, a source with one
                        noise_note if noise_percent is None else None,
                        noise_note if noise_percent is not None else None,
                        getattr(args, 'exclude_flagged', None) is None,
                        len(getattr(args, 'excluded_frequencies', None) or []))

    except Exception as e:
        logger.error(f"Fitting error: {e}")
        logger.error("Try adjusting --circuit expression")
        logger.debug("Traceback:", exc_info=True)
        return None, None

    _save_fit_data(frequencies, Z, Z_fit, result, args, save_suffix)
    return result, _save_fit_figure(frequencies, Z, result, args, save_suffix)
