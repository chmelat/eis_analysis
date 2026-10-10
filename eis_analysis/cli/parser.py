"""
Argument parsing for EIS CLI.

Provides structured argument parsing with logical grouping:
- Input/Output options
- DRT analysis options
- Kramers-Kronig options
- Z-HIT options
- Circuit fitting options
- Voigt chain options
- Oxide analysis options
"""

import argparse
import math
from ..version import get_version_string
from ..drt.core import DRT_WEIGHTINGS
from ..fitting.diffevo import DEFAULT_DE_STRATEGY



def _tau_extend(value: str):
    """argparse type for --tau-extend: 'auto' or a number of decades >= 0."""
    if value == 'auto':
        return value
    try:
        decades = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a number of decades or 'auto', got {value!r}")
    if not (math.isfinite(decades) and decades >= 0):
        raise argparse.ArgumentTypeError(f"must be a finite number >= 0, got {decades}")
    return decades

def _positive_float(value: str) -> float:
    """argparse type for physical quantities that must be finite and > 0."""
    number = float(value)  # argparse reports a ValueError as "invalid value"
    if not (math.isfinite(number) and number > 0):
        raise argparse.ArgumentTypeError(f"must be a finite number > 0, got {number}")
    return number

def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
    args : argparse.Namespace
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description=f'EIS analysis with DRT ({get_version_string()})',
        usage='eis [input] [options]',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  eis                                Synthetic data demo
  eis --ri-fit data.csv              Analyze CSV data file
  eis data.DTA                       Analyze Gamry DTA file
  eis data.DTA --circuit 'R()-(R()|Q())-(R()|Q())'
                                     Fit equivalent circuit
        """
    )

    # Version
    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {get_version_string()}')

    # ==========================================================================
    # Input/Output Group
    # ==========================================================================
    io_group = parser.add_argument_group('Input/Output')

    io_group.add_argument('input', nargs='?', default=None,
                          help='Input file (.DTA for Gamry, .csv for CSV). '
                               'Without argument, synthetic data is used.')

    io_group.add_argument('--f-min', type=float, default=None,
                          help='Minimum frequency [Hz] - data below will be cut off')
    io_group.add_argument('--f-max', type=float, default=None,
                          help='Maximum frequency [Hz] - data above will be cut off')

    io_group.add_argument('--save', '-s', type=str, default=None,
                          help='Save plots and fit results (JSON, CSV) with this prefix')
    io_group.add_argument('--format', '-f', type=str, default='png',
                          choices=['png', 'pdf', 'svg', 'eps'],
                          help='Output format for saved plots (default: png). '
                               'Use pdf/svg/eps for vector graphics.')
    io_group.add_argument('--no-show', action='store_true',
                          help='Do not display plots (useful with --save)')

    io_group.add_argument('--verbose', '-v', action='count', default=0,
                          help='Show debug messages on stderr')
    io_group.add_argument('--quiet', '-q', action='store_true',
                          help='Quiet mode - hide INFO messages, show only warnings and errors')

    # ==========================================================================
    # DRT Analysis Group
    # ==========================================================================
    drt_group = parser.add_argument_group('DRT Analysis')

    drt_group.add_argument('--lambda', '-l', dest='lambda_reg', type=_positive_float,
                           default=None,
                           help='Manual regularization parameter for DRT: weight of the '
                                'roughness integral of gamma against the mean squared '
                                'residual; independent of --n-tau and the frequency '
                                'range. Typical values 1e-9 to 1e-3, more for noisier '
                                'data. Without this, automatic selection (GCV + L-curve) '
                                'is used.')
    drt_group.add_argument('--drt-weighting', type=str, default='sqrt',
                           choices=DRT_WEIGHTINGS,
                           help='Weighting of the DRT least-squares term (default: sqrt, '
                                'w = 1/sqrt|Z|). modulus (1/|Z|) suits noise proportional '
                                'to |Z|; uniform is the unweighted pre-0.38 behaviour. '
                                'Independent of --weighting for circuit fitting.')
    drt_group.add_argument('--tau-extend', type=_tau_extend, default=0.0,
                           metavar='DECADES|auto',
                           help='Extend the DRT tau grid this many decades past the slow end '
                                'of the measured window (default: 0). auto extends only when '
                                'that resolves a pile-up at the slow end and the low-frequency '
                                'end is not capacitive. Peaks past the window are extrapolated.')
    drt_group.add_argument('--drt-inductance', choices=('auto', 'on', 'off'), default='auto',
                           help='Series inductance L in the DRT model (default: auto = only '
                                'when the top decade has a point with Im(Z) > 0). Without it '
                                'the DRT cannot fit an inductive high-frequency end.')
    drt_group.add_argument('--normalize-rpol', action='store_true',
                           help='Normalize gamma(tau) by R_pol so that integral gamma(tau) d(ln tau) = 1.')
    drt_group.add_argument('--n-tau', '-n', type=int, default=100,
                           help='Number of points on tau axis for DRT (default: 100)')
    drt_group.add_argument('--no-voigt-info', action='store_true',
                           help='Do not display Voigt element analysis from DRT')
    drt_group.add_argument('--no-drt', action='store_true',
                           help='Skip DRT analysis')

    # Peak detection
    drt_group.add_argument('--peak-method', type=str, default='scipy',
                           choices=['scipy', 'gmm'],
                           help='Peak detection method: scipy (fast) or gmm (robust). Default: scipy')
    drt_group.add_argument('--gmm-bic-threshold',
                           type=float,
                           default=10.0,
                           metavar='THRESHOLD',
                           help='BIC threshold for GMM peak detection. Lower values detect more peaks. '
                                'Typical range: 2-20. Default: 10.0 (conservative)')
    drt_group.add_argument('--lambda-probe', action='store_true',
                           help='Re-solve DRT at lambda*10^(+-0.5) and lambda*10^(+-1) and report '
                                'per-peak stability (stable/marginal/artifact). Helps distinguish '
                                'real relaxation processes from regularization artifacts.')

    # R_inf estimation
    drt_group.add_argument('--ri-fit', action='store_true',
                           help='Perform robust R_inf estimation before DRT analysis')

    # ==========================================================================
    # Kramers-Kronig Validation Group
    # ==========================================================================
    kk_group = parser.add_argument_group('Kramers-Kronig Validation')

    kk_group.add_argument('--no-kk', action='store_true',
                          help='Skip Kramers-Kronig validation')
    kk_group.add_argument('--mu-threshold', type=float, default=0.85,
                          help='Stopping threshold of the Lin-KK M-iteration '
                               '(default: 0.85). Lower values allow more Voigt '
                               'elements; not a data-quality criterion.')
    kk_group.add_argument('--auto-extend', action=argparse.BooleanOptionalAction,
                          default=True,
                          help='Automatically optimize extend_decades for KK validation '
                               '(minimizes pseudo chi-squared). Extends the tau grid '
                               'below the lowest frequency only, against truncation '
                               'bias on capacitive low-frequency tails; the inductive '
                               'high-frequency tail is covered by the series L. On by '
                               'default; use --no-auto-extend to disable.')
    kk_group.add_argument('--extend-decades-max', type=float, default=1.0,
                          help='Maximum extend_decades for --auto-extend search range '
                               '(searches from 0 to max, default: 1.0)')
    kk_group.add_argument('--kk-series-c', action='store_true',
                          help='Include a series capacitance in the Lin-KK model '
                               '(Schonleber add_cap). Use for blocking/capacitive '
                               'low-frequency behavior (e.g. two-electrode cells), '
                               'where imaginary residuals otherwise grow toward '
                               'low frequencies while the real fit stays good.')

    # ==========================================================================
    # Z-HIT Validation Group
    # ==========================================================================
    zhit_group = parser.add_argument_group('Z-HIT Validation')

    zhit_group.add_argument('--no-zhit', action='store_true',
                            help='Skip Z-HIT validation')

    # ==========================================================================
    # Data Quality Group
    # ==========================================================================
    # Reads the residuals of BOTH validations above, so it belongs to neither.
    quality_group = parser.add_argument_group('Data Quality')

    quality_group.add_argument('--max-residual', type=float, default=5.0,
                               help='Residual threshold for flagging an individual '
                                    'point as suspicious [%%] (default: 5.0). '
                                    'Applies to both KK and Z-HIT; a point is '
                                    'listed when either method exceeds it. '
                                    'Higher = less sensitive.')

    # ==========================================================================
    # Circuit Fitting Group
    # ==========================================================================
    fit_group = parser.add_argument_group('Circuit Fitting')

    fit_group.add_argument('--circuit', '-c', type=str, action='append', default=None,
                           help='Equivalent circuit for fitting. '
                                'Syntax: R(100) - (R(5000) | C(1e-6))  [- = series, | = parallel]. '
                                'Repeat to fit several candidates and rank them by AIC/BIC.')
    fit_group.add_argument('--weighting', type=str, default='modulus',
                           choices=['uniform', 'sqrt', 'modulus', 'proportional'],
                           help='Weighting type for fitting (default: modulus)')
    fit_group.add_argument('--fit-on', type=str, default='original',
                           choices=['original', 'zhit', 'all'],
                           help='Data the fit runs on: original (default), '
                                'zhit (the circuit fit uses the Z-HIT '
                                'reconstruction of |Z| from the phase), or all '
                                '(the reconstruction also feeds R_inf, DRT and '
                                'oxide analysis). Corrects drift of the '
                                'low-frequency modulus; needs Z-HIT validation, '
                                'so it cannot be combined with --no-zhit.')
    fit_group.add_argument('--no-fit', action='store_true',
                           help='Skip equivalent circuit fitting')
    fit_group.add_argument('--numeric-jacobian', action='store_true',
                           help='Use numeric Jacobian instead of analytic (fallback for custom elements)')

    # Optimizer selection
    fit_group.add_argument('--optimizer', type=str, default=None,
                           choices=['single', 'multistart', 'de'],
                           help='Optimizer: de (differential evolution, default), multistart, '
                                'or single (one local fit). --multistart N implies multistart.')

    # Multi-start options
    fit_group.add_argument('--multistart', type=int, default=None, metavar='N',
                           help='Number of restarts for multi-start optimization '
                                '(default: 16). Implies --optimizer multistart.')
    fit_group.add_argument('--multistart-scale', type=float, default=2.0,
                           help='Perturbation scaling in sigma units (default: 2.0)')

    # Differential Evolution options
    fit_group.add_argument('--de-strategy', type=int, default=DEFAULT_DE_STRATEGY, choices=[1, 2, 3],
                           help='DE strategy: 1=randtobest1bin, 2=best1bin, 3=rand1bin (default; '
                                'the two others converge faster but get trapped in local minima)')
    fit_group.add_argument('--de-popsize', type=int, default=15,
                           help='DE population size multiplier (default: 15)')
    fit_group.add_argument('--de-maxiter', type=int, default=1000,
                           help='DE maximum generations (default: 1000)')
    fit_group.add_argument('--de-tol', type=float, default=0.01,
                           help='DE convergence tolerance (default: 0.01)')
    fit_group.add_argument('--de-workers', type=int, default=1,
                           help='DE parallel workers (default: 1, use -1 for all CPUs)')
    fit_group.add_argument('--no-archive-check', action='store_true',
                           help='Skip the DE archive check (refining early-generation candidates '
                                'against local minima, ~1 s per fit)')

    # ==========================================================================
    # Voigt Chain Group
    # ==========================================================================
    voigt_group = parser.add_argument_group('Voigt Chain Fitting')

    voigt_group.add_argument('--voigt-chain', action='store_true',
                             help='Use automatic Voigt chain fitting via linear regression.')
    voigt_group.add_argument('--voigt-n-per-decade', type=int, default=3,
                             help='Time constants per decade for --voigt-chain (default: 3)')
    voigt_group.add_argument('--voigt-extend-decades', type=float, default=0.0,
                             help='Extend tau range by N decades (default: 0.0)')
    voigt_group.add_argument('--voigt-prune-threshold', type=float, default=0.01,
                             help='Threshold for removing small R_i (default: 0.01)')
    voigt_group.add_argument('--voigt-allow-negative', action='store_true',
                             help='Allow negative R_i values (Lin-KK style)')
    voigt_group.add_argument('--voigt-no-inductance', action='store_true',
                             help='Do not include series inductance L')
    voigt_group.add_argument('--voigt-fit-type', type=str, default='complex',
                             choices=['real', 'imag', 'complex'],
                             help='Fit type: complex (default), real, or imag')
    voigt_group.add_argument('--voigt-auto-M', action='store_true',
                             help='Auto-optimize M elements using mu metric')
    voigt_group.add_argument('--voigt-mu-threshold', type=float, default=0.85,
                             help='Stopping threshold for --voigt-auto-M (default: '
                                  '0.85; lower values allow more elements)')
    voigt_group.add_argument('--voigt-max-M', type=int, default=50,
                             help='Maximum M elements for --voigt-auto-M (default: 50)')

    # ==========================================================================
    # Oxide Analysis Group
    # ==========================================================================
    oxide_group = parser.add_argument_group('Oxide Layer Analysis')

    oxide_group.add_argument('--analyze-oxide', action='store_true',
                             help='Perform oxide layer analysis')
    oxide_group.add_argument('--local-exponent', action='store_true',
                             help='Map the local CPE exponent n(f) from the real '
                                  'part of the admittance - shows where one CPE '
                                  'fails; needs no circuit fit')
    oxide_group.add_argument('--epsilon-r', type=_positive_float, default=None,
                             help='Relative permittivity of oxide (default: 22 for ZrO2)')
    oxide_group.add_argument('--thickness', type=_positive_float, default=None,
                             help='Known oxide thickness [nm] - switches '
                                  '--analyze-oxide to permittivity estimation')
    oxide_group.add_argument('--area', type=_positive_float, default=None,
                             help='Electrode area in cm^2 (default: from DTA '
                                  'metadata, else 1.0)')
    oxide_group.add_argument('--rho-delta', type=_positive_float, default=None,
                             help='Film resistivity at the electrolyte interface '
                                  '[Ohm cm] - enables the power-law (Hirschorn-'
                                  'Orazem) thickness for a CPE')

    # ==========================================================================
    # Visualization Group
    # ==========================================================================
    vis_group = parser.add_argument_group('Visualization')

    vis_group.add_argument('--ocv', action='store_true',
                           help='Enable OCV (Open Circuit Voltage) curve visualization')

    args = parser.parse_args()

    # Resolve optimizer vs. --multistart (audit P1). Previously --multistart
    # without --optimizer multistart was silently ignored and the fit ran DE.
    if args.multistart is not None and args.multistart <= 0:
        parser.error('--multistart must be a positive integer')
    if args.optimizer is None:
        args.optimizer = 'multistart' if args.multistart is not None else 'de'
    elif args.optimizer != 'multistart' and args.multistart is not None:
        parser.error(
            f"--multistart has no effect with --optimizer {args.optimizer}; "
            "use --optimizer multistart (or drop --multistart)"
        )

    # --fit-on feeds on the Z-HIT reconstruction, so switching Z-HIT off leaves
    # it with nothing to fit. Caught here rather than mid-run, where the
    # spectrum is already loaded and the message arrives after the validation
    # sections have scrolled past.
    if args.fit_on != 'original' and args.no_zhit:
        parser.error(
            f"--fit-on {args.fit_on} needs the Z-HIT reconstruction; "
            "drop --no-zhit"
        )

    # --fit-on zhit corrects the data the circuit fit reads, and nothing else.
    # With --no-fit there is no such reader, so the flag would only print a
    # section and - if Z-HIT failed - abort a run that was never going to fit.
    # --fit-on all is different: it feeds R_inf and DRT too.
    if args.fit_on == 'zhit' and args.no_fit:
        parser.error(
            "--fit-on zhit has no effect with --no-fit; "
            "use --fit-on all to correct R_inf and DRT as well"
        )

    return args
