#!/usr/bin/env python3
"""
EIS Analysis with DRT
=====================

CLI tool for electrochemical impedance spectroscopy analysis.

Version: Imported from eis_analysis.version (single source of truth)

Features:
- Kramers-Kronig validation (Lin-KK method)
- DRT analysis with automatic lambda selection (GCV + L-curve)
- GMM peak detection in DRT spectrum (--peak-method gmm)
- Equivalent circuit fitting with operator syntax
- Oxide layer thickness estimation
- Robust R_inf estimation (--ri-fit)

Usage:
    eis                             # synthetic data demo
    eis data.DTA                    # Gamry file (auto-lambda default)
    eis data.DTA --lambda 0.5       # manual regularization
    eis data.DTA --peak-method gmm  # GMM peak detection
    eis data.DTA --ri-fit -v        # robust R_inf estimation
    eis data.DTA --voigt-chain --analyze-oxide  # oxide analysis

    eis --help                      # help
"""

import argparse
import logging
import os
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from numpy.typing import NDArray

# Import refactored modules
from eis_analysis import (
    # Version info
    get_version_string,
    # I/O
    load_data,
    load_csv_data,
    parse_dta_metadata,
    parse_ocv_curve,
    log_metadata,
    generate_synthetic_data,
    # Validation
    kramers_kronig_validation,
    # DRT
    calculate_drt,
    # Fitting (new operator overloading approach v2.2.0)
    fit_equivalent_circuit,
    fit_circuit_multistart,  # Multi-start optimization
    fit_circuit_diffevo,  # Differential evolution optimization
    R, C, Q, L, W, Wo, K,  # Circuit elements
    fit_voigt_chain_linear,  # Voigt chain linear fitting
    # Analysis
    analyze_oxide_layer,
    # Visualization
    visualize_data,
    visualize_ocv,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Custom Logging
# =============================================================================

class InfoFormatter(logging.Formatter):
    """Formatter for INFO level - no prefix, clean output."""
    def format(self, record):
        return record.getMessage()


class WarningFormatter(logging.Formatter):
    """Formatter for WARNING level - prefix with '!'."""
    def format(self, record):
        return f"! {record.getMessage()}"


class ErrorFormatter(logging.Formatter):
    """Formatter for ERROR/CRITICAL level - prefix with '!!'."""
    def format(self, record):
        return f"!! {record.getMessage()}"


class DebugFormatter(logging.Formatter):
    """Formatter for DEBUG level - prefix with '[DEBUG]'."""
    def format(self, record):
        return f"[DEBUG] {record.getMessage()}"


class LevelFilter(logging.Filter):
    """Filter that accepts only specific log levels."""
    def __init__(self, levels):
        super().__init__()
        self.levels = levels if isinstance(levels, (list, tuple)) else [levels]

    def filter(self, record):
        return record.levelno in self.levels


# =============================================================================
# Logging Helpers
# =============================================================================

def log_separator(length: int = 50, char: str = "=") -> None:
    """Log a separator line for visual clarity."""
    logger.info(char * length)


# =============================================================================
# Exceptions
# =============================================================================

class EISAnalysisError(Exception):
    """Base exception for EIS analysis errors."""
    pass


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class LoadedData:
    """Container for loaded EIS data."""
    frequencies: NDArray[np.float64]
    Z: NDArray[np.complex128]
    title: str
    metadata: Optional[dict]
    ocv_data: Optional[dict] = None  # OCV data from DTA file


@dataclass
class DRTResult:
    """Container for DRT analysis results."""
    tau: Optional[NDArray]
    gamma: Optional[NDArray]
    peaks_gmm: Optional[list]
    fig_drt: Optional[plt.Figure]
    fig_rinf: Optional[plt.Figure]


# =============================================================================
# Helper Functions
# =============================================================================

def save_figure(
    fig: Optional[plt.Figure],
    prefix: Optional[str],
    suffix: str,
    fmt: str = 'png'
) -> None:
    """
    Save figure to file if fig and prefix are provided.

    Parameters
    ----------
    fig : Figure or None
        Matplotlib figure to save
    prefix : str or None
        File prefix (from --save argument)
    suffix : str
        File suffix (e.g., 'nyquist_bode', 'kk', 'drt')
    fmt : str
        Output format: 'png', 'pdf', 'svg', 'eps' (default: 'png')
    """
    if fig is None or prefix is None:
        return

    filepath = f"{prefix}_{suffix}.{fmt}"
    try:
        # For bitmap formats use dpi, for vector formats use bbox_inches
        if fmt == 'png':
            fig.savefig(filepath, dpi=150, bbox_inches='tight')
        else:
            # Vector formats (pdf, svg, eps)
            fig.savefig(filepath, bbox_inches='tight')
        logger.info(f"Saved: {filepath}")
    except Exception as e:
        logger.error(f"Error saving figure: {e}")


def parse_circuit_expression(expr: str):
    """
    Parse circuit expression string into circuit object.

    Parameters
    ----------
    expr : str
        Circuit expression, e.g. "R(100) - (R(5000) | C(1e-6))"

    Returns
    -------
    circuit : Circuit object
        Parsed circuit

    Examples
    --------
    >>> circuit = parse_circuit_expression("R(100) - (R(5000) | C(1e-6))")
    >>> print(circuit.get_all_params())
    [100, 5000, 1e-6]
    """
    # Safe namespace for eval - only circuit elements
    safe_namespace = {
        'R': R,
        'C': C,
        'Q': Q,
        'L': L,
        'W': W,
        'Wo': Wo,
        'K': K,
    }

    try:
        circuit = eval(expr, {"__builtins__": {}}, safe_namespace)
        return circuit
    except Exception as e:
        raise ValueError(f"Cannot parse circuit expression '{expr}': {e}")


# =============================================================================
# Setup Functions
# =============================================================================

def setup_logging(args: argparse.Namespace) -> None:
    """
    Configure logging based on command line arguments.

    Output behavior:
    - Default: INFO + WARNING on stdout, ERROR on stderr
    - Quiet (-q): WARNING on stdout, ERROR on stderr (no INFO)
    - Verbose (-v): DEBUG on stderr + default behavior

    Prefixes:
    - INFO: no prefix (clean output)
    - WARNING: "! " prefix
    - ERROR: "!! " prefix
    - DEBUG: "[DEBUG] " prefix
    """
    import sys

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)  # Allow all levels, filter per handler

    # Clear any existing handlers
    root_logger.handlers.clear()

    # Determine what to show
    show_info = not args.quiet
    show_debug = args.verbose >= 1

    # Handler for INFO -> stdout (no prefix)
    if show_info:
        info_handler = logging.StreamHandler(sys.stdout)
        info_handler.setLevel(logging.INFO)
        info_handler.addFilter(LevelFilter(logging.INFO))
        info_handler.setFormatter(InfoFormatter())
        root_logger.addHandler(info_handler)

    # Handler for WARNING -> stdout (prefix "!")
    warning_handler = logging.StreamHandler(sys.stdout)
    warning_handler.setLevel(logging.WARNING)
    warning_handler.addFilter(LevelFilter(logging.WARNING))
    warning_handler.setFormatter(WarningFormatter())
    root_logger.addHandler(warning_handler)

    # Handler for ERROR/CRITICAL -> stderr (prefix "!!")
    error_handler = logging.StreamHandler(sys.stderr)
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(ErrorFormatter())
    root_logger.addHandler(error_handler)

    # Handler for DEBUG -> stderr (prefix "[DEBUG]")
    if show_debug:
        debug_handler = logging.StreamHandler(sys.stderr)
        debug_handler.setLevel(logging.DEBUG)
        debug_handler.addFilter(LevelFilter(logging.DEBUG))
        debug_handler.setFormatter(DebugFormatter())
        root_logger.addHandler(debug_handler)


class OnePerLineHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter that puts each option on a separate line in usage."""

    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = 'usage: '

        # If custom usage is provided, use it
        if usage is not None:
            usage = usage % dict(prog=self._prog)
            return f'{prefix}{usage}\n\n'

        # Otherwise build usage with one option per line
        prog = '%(prog)s' % dict(prog=self._prog)
        lines = [f'{prefix}{prog}']

        for action in actions:
            if action.option_strings:
                # Optional argument
                option = action.option_strings[0]
                if action.nargs == 0:
                    # Flag (store_true/store_false)
                    lines.append(f'              [{option}]')
                elif action.metavar:
                    lines.append(f'              [{option} {action.metavar}]')
                elif action.dest:
                    lines.append(f'              [{option} {action.dest.upper()}]')
                else:
                    lines.append(f'              [{option}]')
            elif not action.option_strings and action.dest != 'help':
                # Positional argument
                if action.nargs == '?':
                    lines.append(f'              [{action.dest}]')
                else:
                    lines.append(f'              {action.dest}')

        return '\n'.join(lines) + '\n\n'


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='EIS analysis with DRT',
        usage='eis [input] [options]',
        formatter_class=OnePerLineHelpFormatter,
        epilog="""
Examples:
  eis                                Synthetic data demo
  eis --ri-fit data.csv              Analyze CSV data file
  eis data.DTA                       Analyze Gamry DTA file
  eis data.DTA --circuit 'R()-(R()|Q())-(R()|Q())'
                                     Fit equivalent circuit
        """
    )

    # Input
    parser.add_argument('input', nargs='?', default=None,
                        help='Input file (.DTA for Gamry, .csv for CSV). '
                             'Without argument, synthetic data is used.')
    # Frequency filtering
    parser.add_argument('--f-min', type=float, default=None,
                        help='Minimum frequency [Hz] - data below will be cut off')
    parser.add_argument('--f-max', type=float, default=None,
                        help='Maximum frequency [Hz] - data above will be cut off')

    # Circuit fitting
    parser.add_argument('--circuit', '-c', type=str, default=None,
                        help='Equivalent circuit for fitting. '
                             'Syntax: R(100) - (R(5000) | C(1e-6))  [- = series, | = parallel].')
    parser.add_argument('--weighting', type=str, default='sqrt',
                        choices=['uniform', 'sqrt', 'proportional', 'square'],
                        help='Weighting type for fitting (default: sqrt)')

    # Multi-start
    parser.add_argument('--multistart', type=int, default=0, metavar='N',
                        help='Multi-start optimization with N restarts (default: 0 = disabled)')
    parser.add_argument('--multistart-scale', type=float, default=2.0,
                        help='Perturbation scaling in sigma units (default: 2.0)')

    # Differential Evolution options
    parser.add_argument('--optimizer', type=str, default='de',
                        choices=['single', 'multistart', 'de'],
                        help='Optimizer: de (differential evolution, default), multistart, or single (one local fit)')
    parser.add_argument('--de-strategy', type=int, default=1, choices=[1, 2, 3],
                        help='DE strategy: 1=randtobest1bin (default), 2=best1bin, 3=rand1bin')
    parser.add_argument('--de-popsize', type=int, default=15,
                        help='DE population size multiplier (default: 15)')
    parser.add_argument('--de-maxiter', type=int, default=1000,
                        help='DE maximum generations (default: 1000)')
    parser.add_argument('--de-tol', type=float, default=0.01,
                        help='DE convergence tolerance (default: 0.01)')
    parser.add_argument('--de-workers', type=int, default=1,
                        help='DE parallel workers (default: 1, use -1 for all CPUs)')
    parser.add_argument('--numeric-jacobian', action='store_true',
                        help='Use numeric Jacobian instead of analytic (fallback for custom elements)')

    # DRT options
    parser.add_argument('--lambda', '-l', dest='lambda_reg', type=float,
                        default=None,
                        help='Manual regularization parameter for DRT. '
                             'Without this, automatic selection (GCV + L-curve) is used.')
    parser.add_argument('--normalize-rpol', action='store_true',
                        help='Normalize γ(τ) by R_pol so that ∫γ(τ)d(ln τ) = 1.')
    parser.add_argument('--n-tau', '-n', type=int, default=100,
                        help='Number of points on tau axis for DRT (default: 100)')
    parser.add_argument('--no-voigt-info', action='store_true',
                        help='Do not display Voigt element analysis from DRT')
    parser.add_argument('--classify-terms', action='store_true',
                        help='Classify term types from DRT peaks. Enables GMM detection.')

    # Peak detection
    parser.add_argument('--peak-method', type=str, default='scipy',
                        choices=['scipy', 'gmm'],
                        help='Peak detection method: scipy (fast) or gmm (robust). Default: scipy')
    
    # KK options
    parser.add_argument('--mu-threshold', type=float, default=0.85,
                        help='μ metric threshold for Lin-KK test (default: 0.85)')

    # Skip options
    parser.add_argument('--no-kk', action='store_true',
                        help='Skip Kramers-Kronig validation')
    parser.add_argument('--no-drt', action='store_true',
                        help='Skip DRT analysis')
    parser.add_argument('--no-fit', action='store_true',
                        help='Skip equivalent circuit fitting')
    parser.add_argument('--ocv', action='store_true',
                        help='Enable OCV (Open Circuit Voltage) curve visualization')

    # Output options
    parser.add_argument('--save', '-s', type=str, default=None,
                        help='Save plots to files with this prefix')
    parser.add_argument('--format', '-f', type=str, default='png',
                        choices=['png', 'pdf', 'svg', 'eps'],
                        help='Output format for saved plots (default: png). '
                             'Use pdf/svg/eps for vector graphics.')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display plots (useful with --save)')
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Show debug messages on stderr')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Quiet mode - hide INFO messages, show only warnings and errors')

    # Oxide analysis
    parser.add_argument('--analyze-oxide', action='store_true',
                        help='Perform oxide layer analysis')
    parser.add_argument('--epsilon-r', type=float, default=22.0,
                        help='Relative permittivity of oxide (default: 22 for ZrO2)')
    parser.add_argument('--area', type=float, default=1.0,
                        help='Electrode area in cm² (default: 1.0)')

    # R_inf estimation
    parser.add_argument('--ri-fit', action='store_true',
                        help='Perform robust R_inf estimation before DRT analysis')

    # Voigt chain options
    parser.add_argument('--voigt-chain', action='store_true',
                        help='Use automatic Voigt chain fitting via linear regression.')
    parser.add_argument('--voigt-n-per-decade', type=int, default=3,
                        help='Time constants per decade for --voigt-chain (default: 3)')
    parser.add_argument('--voigt-extend-decades', type=float, default=0.0,
                        help='Extend τ range by N decades (default: 0.0)')
    parser.add_argument('--voigt-prune-threshold', type=float, default=0.01,
                        help='Threshold for removing small R_i (default: 0.01)')
    parser.add_argument('--voigt-allow-negative', action='store_true',
                        help='Allow negative R_i values (Lin-KK style)')
    parser.add_argument('--voigt-no-inductance', action='store_true',
                        help='Do not include series inductance L')
    parser.add_argument('--voigt-fit-type', type=str, default='complex',
                        choices=['real', 'imag', 'complex'],
                        help='Fit type: complex (default), real, or imag')
    parser.add_argument('--voigt-auto-M', action='store_true',
                        help='Auto-optimize M elements using μ metric')
    parser.add_argument('--voigt-mu-threshold', type=float, default=0.85,
                        help='μ threshold for --voigt-auto-M (default: 0.85)')
    parser.add_argument('--voigt-max-M', type=int, default=50,
                        help='Maximum M elements for --voigt-auto-M (default: 50)')

    return parser.parse_args()


# =============================================================================
# Data Loading and Filtering
# =============================================================================

def load_eis_data(args: argparse.Namespace) -> LoadedData:
    """
    Load EIS data from file or generate synthetic data.

    Parameters
    ----------
    args : Namespace
        Parsed command line arguments

    Returns
    -------
    LoadedData
        Container with frequencies, Z, title, metadata, and ocv_data
    """
    metadata = None
    ocv_data = None

    if args.input is None:
        # Synthetic data
        frequencies, Z = generate_synthetic_data(
            Rs=10, R0=1e5, Q0=(1e-6, 0.6),
            R1=8e5, Q1=(3e-5, 0.43), noise=0.01
        )
        title = "Synthetic data"
    else:
        # File input
        if not os.path.exists(args.input):
            raise EISAnalysisError(f"File '{args.input}' does not exist!")

        ext = os.path.splitext(args.input)[1].lower()
        try:
            if ext == '.dta':
                frequencies, Z = load_data(args.input)
                metadata = parse_dta_metadata(args.input)
                log_metadata(metadata)
                # Load OCV data if available
                ocv_data = parse_ocv_curve(args.input)
                if ocv_data is not None:
                    logger.info(f"OCV data: {len(ocv_data['time'])} points, "
                               f"duration {ocv_data['time'][-1]/60:.1f} min")
            elif ext == '.csv':
                frequencies, Z = load_csv_data(args.input)
            else:
                raise EISAnalysisError(
                    f"Unsupported format '{ext}'. Supported: .DTA (Gamry), .csv"
                )
        except EISAnalysisError:
            raise
        except Exception as e:
            raise EISAnalysisError(f"Error loading file: {e}") from e

        title = os.path.basename(args.input)

    return LoadedData(frequencies=frequencies, Z=Z, title=title, metadata=metadata, ocv_data=ocv_data)


def filter_by_frequency(
    data: LoadedData,
    args: argparse.Namespace
) -> LoadedData:
    """
    Filter data by frequency range.

    Parameters
    ----------
    data : LoadedData
        Loaded EIS data
    args : Namespace
        Parsed command line arguments

    Returns
    -------
    LoadedData
        Filtered data (or original if no filtering needed)
    """
    if args.f_min is None and args.f_max is None:
        return data

    frequencies = data.frequencies
    Z = data.Z
    original_count = len(frequencies)

    mask = np.ones(len(frequencies), dtype=bool)

    if args.f_min is not None:
        mask &= (frequencies >= args.f_min)
        logger.info(f"Applying f_min = {args.f_min} Hz")

    if args.f_max is not None:
        mask &= (frequencies <= args.f_max)
        logger.info(f"Applying f_max = {args.f_max} Hz")

    frequencies = frequencies[mask]
    Z = Z[mask]

    filtered_count = len(frequencies)
    removed_count = original_count - filtered_count
    logger.info(f"Frequency filter: {original_count} -> {filtered_count} points (removed: {removed_count})")

    if filtered_count == 0:
        raise EISAnalysisError("No data remaining after filtering! Check --f-min and --f-max.")

    return LoadedData(
        frequencies=frequencies,
        Z=Z,
        title=data.title,
        metadata=data.metadata,
        ocv_data=data.ocv_data
    )


# =============================================================================
# Analysis Functions
# =============================================================================

def run_kk_validation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Optional[plt.Figure]:
    """
    Run Kramers-Kronig validation.

    Returns
    -------
    fig : Figure or None
        KK validation figure
    """
    if args.no_kk:
        return None

    M, mu, Z_kk, residuals, fig = kramers_kronig_validation(
        frequencies, Z,
        mu_threshold=args.mu_threshold
    )
    save_figure(fig, args.save, 'kk', args.format)
    return fig


def run_rinf_estimation(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[float], Optional[plt.Figure]]:
    """
    Run R_inf estimation if --ri-fit is specified.

    Returns
    -------
    R_inf : float or None
        Estimated R_inf value
    fig : Figure or None
        R_inf fit figure
    """
    if not args.ri_fit:
        return None, None

    from eis_analysis.rinf_estimation import estimate_rinf_with_inductance

    log_separator()
    logger.info("R_inf estimation (high-frequency resistance)")
    log_separator()

    try:
        R_inf, L_fit, circuit_rl, diag_rl, fig = estimate_rinf_with_inductance(
            frequencies, Z, verbose=False, plot=True
        )

        # Print results
        method = diag_rl.get('method', 'unknown')
        logger.info(f"R_inf = {R_inf:.3f} Ohm ({diag_rl['n_points_used']} HF points)")

        if 'r_squared' in diag_rl and diag_rl['r_squared'] > 0:
            logger.info(f"  Quality: R² = {diag_rl['r_squared']:.4f}")
        if 'L_nH' in diag_rl and diag_rl['L_nH'] > 0:
            logger.info(f"  Inductance: L = {diag_rl['L_nH']:.2f} nH")

        # For Voigt fit: show R_ct and C
        if 'voigt' in method and 'R_ct' in diag_rl and 'C_nF' in diag_rl:
            logger.info(f"  Voigt params: R_ct = {diag_rl['R_ct']:.2f} Ohm, C = {diag_rl['C_nF']:.2f} nF")
            if 'f_characteristic' in diag_rl:
                logger.info(f"  Characteristic freq: f_char = {diag_rl['f_characteristic']/1e6:.3f} MHz")

        # Comparison with median
        n_avg = min(5, max(1, len(frequencies) // 10))
        high_freq_indices = np.argsort(frequencies)[-n_avg:]
        R_inf_median = np.median(Z.real[high_freq_indices])
        diff_abs = R_inf - R_inf_median
        diff_pct = (diff_abs / R_inf_median * 100) if R_inf_median != 0 else 0
        logger.info(f"  For comparison: median = {R_inf_median:.3f} Ohm (diff: {diff_abs:+.3f} Ohm, {diff_pct:+.1f}%)")

        # Warnings
        if diag_rl.get('warnings'):
            for warning in diag_rl['warnings']:
                logger.warning(f"  {warning}")

        log_separator()
        save_figure(fig, args.save, 'ri_fit', args.format)

        return R_inf, fig

    except Exception as e:
        logger.error(f"R_inf estimation failed: {e}")
        logger.debug("Traceback:", exc_info=True)
        return None, None


def run_drt_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    R_inf_computed: Optional[float],
    peak_method: str
) -> DRTResult:
    """
    Run DRT analysis.

    Returns
    -------
    DRTResult
        Container with tau, gamma, peaks, and figures
    """
    if args.no_drt:
        return DRTResult(tau=None, gamma=None, peaks_gmm=None, fig_drt=None, fig_rinf=None)

    use_auto_lambda = args.lambda_reg is None

    tau, gamma, fig_drt, peaks_gmm, fig_rinf = calculate_drt(
        frequencies, Z,
        n_tau=args.n_tau,
        lambda_reg=args.lambda_reg,
        auto_lambda=use_auto_lambda,
        normalize_rpol=args.normalize_rpol,
        peak_method=peak_method,
        use_rl_fit=False,
        use_voigt_fit=args.ri_fit,
        r_inf_preset=R_inf_computed
    )

    save_figure(fig_drt, args.save, 'drt', args.format)

    # Save R_inf figure from DRT only if not already saved via --ri-fit
    if not args.ri_fit:
        save_figure(fig_rinf, args.save, 'ri_fit', args.format)

    return DRTResult(
        tau=tau,
        gamma=gamma,
        peaks_gmm=peaks_gmm,
        fig_drt=fig_drt,
        fig_rinf=fig_rinf
    )


def run_voigt_analysis(
    drt_result: DRTResult,
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> None:
    """Run Voigt element analysis from DRT results."""
    if args.no_drt or args.no_voigt_info:
        return
    if drt_result.tau is None or drt_result.gamma is None:
        return

    from eis_analysis.fitting import analyze_voigt_elements, format_voigt_report

    try:
        voigt_info = analyze_voigt_elements(
            drt_result.tau, drt_result.gamma, frequencies, Z,
            peaks_gmm=drt_result.peaks_gmm,
            classify_terms=args.classify_terms
        )
        report = format_voigt_report(voigt_info)
        logger.info(report)

    except Exception as e:
        logger.warning(f"Voigt element analysis failed: {e}")
        logger.debug(f"Traceback: {e}", exc_info=True)


def run_circuit_fitting(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[object], Optional[plt.Figure]]:
    """
    Run equivalent circuit fitting.

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
        logger.info("For fitting, use: --circuit 'R(100) - (R(1000) | C(1e-6))' or --voigt-chain")
        return None, None

    # Default circuit for synthetic data: Rs - (R||Q) - (R||Q)
    circuit_expr = args.circuit
    if args.input is None and circuit_expr is None:
        circuit_expr = 'R() - (R() | Q()) - (R() | Q())'

    # Voigt chain fitting
    if args.voigt_chain:
        return _fit_voigt_chain(frequencies, Z, args)
    else:
        return _fit_standard_circuit(frequencies, Z, args, circuit_expr)


def _fit_voigt_chain(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace
) -> Tuple[Optional[object], Optional[plt.Figure]]:
    """Fit using Voigt chain linear method."""
    from eis_analysis.fitting.circuit import FitResult
    from eis_analysis.fitting.diagnostics import compute_fit_metrics

    log_separator()
    logger.info("Using Voigt chain linear fit")
    log_separator()

    try:
        circuit, initial_params = fit_voigt_chain_linear(
            frequencies, Z,
            n_per_decade=args.voigt_n_per_decade,
            extend_decades=args.voigt_extend_decades,
            include_L=not args.voigt_no_inductance,
            fit_type=args.voigt_fit_type,
            prune_threshold=args.voigt_prune_threshold,
            allow_negative=args.voigt_allow_negative,
            auto_optimize_M=args.voigt_auto_M,
            mu_threshold=args.voigt_mu_threshold,
            max_M=args.voigt_max_M,
            weighting=args.weighting
        )
        logger.info(f"Created circuit: {circuit}")
        logger.info(f"Number of parameters: {len(initial_params)}")
    except Exception as e:
        raise EISAnalysisError(f"Voigt chain linear fit error: {e}") from e

    log_separator()
    logger.info("Using linear result (no nonlinear optimization)")
    log_separator()

    # Compute Z_fit
    Z_fit = circuit.impedance(frequencies, initial_params)

    # Compute fit metrics using shared utility function
    fit_error_rel, fit_error_abs, quality = compute_fit_metrics(Z, Z_fit, args.weighting)

    logger.info("Fit error:")
    logger.info(f"  Relative:  {fit_error_rel:.2f}%")
    logger.info(f"  Absolute:  {fit_error_abs:.2e} Ohm")
    logger.info(f"  Quality:   {quality}")

    # Update circuit with fitted parameters
    if hasattr(circuit, 'update_params'):
        circuit.update_params(initial_params)

    # Create FitResult
    result = FitResult(
        circuit=circuit,
        params_opt=np.array(initial_params),
        params_stderr=np.zeros(len(initial_params)),
        fit_error_rel=fit_error_rel,
        fit_error_abs=fit_error_abs,
        quality=quality,
        _n_data=len(frequencies)
    )

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(Z.real, -Z.imag, 'o', label='Data', markersize=6, alpha=0.7)
    ax.plot(Z_fit.real, -Z_fit.imag, '-', label='Linear fit', linewidth=2)
    ax.set_xlabel("Z' [Ohm]")
    ax.set_ylabel("-Z'' [Ohm]")
    ax.set_title(f"Voigt chain - Linear fit (error: {fit_error_rel:.2f}%)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    plt.tight_layout()

    save_figure(fig, args.save, 'fit', args.format)

    return result, fig


def _fit_standard_circuit(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    circuit_expr: str
) -> Tuple[Optional[object], Optional[plt.Figure]]:
    """Fit using standard circuit with nonlinear optimization."""
    try:
        circuit = parse_circuit_expression(circuit_expr)
    except ValueError as e:
        raise EISAnalysisError(
            f"Circuit parsing error: {e}\n"
            "Use syntax like: --circuit 'R(100) - (R(5000) | C(1e-6))'"
        ) from e

    # Log circuit expression
    logger.info("="*60)
    logger.info("Equivalent circuit")
    logger.info("="*60)
    logger.info(f"Circuit: {circuit_expr}")

    try:
        if args.optimizer == 'de':
            # Differential Evolution global optimization
            diffevo_result, Z_fit, fig = fit_circuit_diffevo(
                circuit,
                frequencies, Z,
                strategy=args.de_strategy,
                popsize=args.de_popsize,
                maxiter=args.de_maxiter,
                tol=args.de_tol,
                workers=args.de_workers,
                weighting=args.weighting,
                verbose=True,
                use_analytic_jacobian=not args.numeric_jacobian
            )
            result = diffevo_result.best_result
        elif args.optimizer == 'multistart':
            # Multi-start optimization
            n_restarts = args.multistart if args.multistart > 0 else 16
            multistart_result, Z_fit, fig = fit_circuit_multistart(
                circuit,
                frequencies, Z,
                n_restarts=n_restarts,
                scale=args.multistart_scale,
                weighting=args.weighting,
                verbose=True,
                use_analytic_jacobian=not args.numeric_jacobian
            )
            result = multistart_result.best_result
        else:
            # Single fit (args.optimizer == 'single')
            result, Z_fit, fig = fit_equivalent_circuit(
                frequencies, Z,
                circuit,
                weighting=args.weighting,
                use_analytic_jacobian=not args.numeric_jacobian
            )

        save_figure(fig, args.save, 'fit', args.format)
        return result, fig

    except Exception as e:
        logger.error(f"Fitting error: {e}")
        logger.error("Try adjusting --circuit expression")
        logger.debug("Traceback:", exc_info=True)
        return None, None


def run_oxide_analysis(
    frequencies: NDArray,
    Z: NDArray,
    args: argparse.Namespace,
    fitted_result: Optional[object],
    metadata: Optional[dict]
) -> None:
    """Run oxide layer analysis."""
    if not args.analyze_oxide:
        return

    # Use area from metadata if available and not explicitly specified
    area_to_use = args.area
    if metadata is not None and metadata.get('area') is not None:
        if args.area == 1.0:  # Default value was not changed
            area_to_use = metadata['area']
            logger.info(f"Using area from DTA metadata: {area_to_use:.4f} cm2")
        else:
            logger.info(f"Using explicitly specified area: {area_to_use:.4f} cm2 (metadata: {metadata['area']:.4f} cm2)")

    analyze_oxide_layer(
        frequencies, Z,
        epsilon_r=args.epsilon_r,
        area_cm2=area_to_use,
        fit_result=fitted_result
    )


# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main CLI entry point."""
    # Parse arguments and setup
    args = parse_arguments()
    setup_logging(args)

    try:
        _run_analysis(args)
    except EISAnalysisError as e:
        logger.error(str(e))
        exit(1)
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        exit(130)


def _run_analysis(args: argparse.Namespace) -> None:
    """Run the full analysis pipeline."""
    # Auto-enable GMM if --classify-terms is used
    peak_method = args.peak_method
    if args.classify_terms and peak_method != 'gmm':
        logger.info("--classify-terms requires GMM peak detection")
        logger.info("Automatically enabling --peak-method gmm")
        peak_method = 'gmm'

    log_separator(60)
    logger.info(f"EIS Analysis ({get_version_string()})")
    log_separator(60)

    # Load and filter data
    data = load_eis_data(args)
    data = filter_by_frequency(data, args)

    # Visualization
    fig_vis = visualize_data(data.frequencies, data.Z, data.title)
    save_figure(fig_vis, args.save, 'nyquist_bode', args.format)

    # OCV visualization (for DTA files with OCV data)
    if args.ocv and data.ocv_data is not None:
        fig_ocv = visualize_ocv(data.ocv_data, data.title)
        save_figure(fig_ocv, args.save, 'ocv', args.format)

    # Kramers-Kronig validation
    run_kk_validation(data.frequencies, data.Z, args)

    # R_inf estimation
    R_inf_computed, _ = run_rinf_estimation(data.frequencies, data.Z, args)

    # DRT analysis
    drt_result = run_drt_analysis(data.frequencies, data.Z, args, R_inf_computed, peak_method)

    # Voigt element analysis
    run_voigt_analysis(drt_result, data.frequencies, data.Z, args)

    # Circuit fitting
    fitted_result, _ = run_circuit_fitting(data.frequencies, data.Z, args)

    # Oxide layer analysis
    run_oxide_analysis(data.frequencies, data.Z, args, fitted_result, data.metadata)

    # Display plots
    if not args.no_show:
        plt.show()

    log_separator(60)
    logger.info("Analysis complete")
    log_separator(60)


if __name__ == "__main__":
    main()
