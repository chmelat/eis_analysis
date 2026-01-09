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
from ..version import get_version_string


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
                          help='Save plots to files with this prefix')
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

    drt_group.add_argument('--lambda', '-l', dest='lambda_reg', type=float,
                           default=None,
                           help='Manual regularization parameter for DRT. '
                                'Without this, automatic selection (GCV + L-curve) is used.')
    drt_group.add_argument('--normalize-rpol', action='store_true',
                           help='Normalize gamma(tau) by R_pol so that integral gamma(tau) d(ln tau) = 1.')
    drt_group.add_argument('--n-tau', '-n', type=int, default=100,
                           help='Number of points on tau axis for DRT (default: 100)')
    drt_group.add_argument('--no-voigt-info', action='store_true',
                           help='Do not display Voigt element analysis from DRT')
    drt_group.add_argument('--classify-terms', action='store_true',
                           help='Classify term types from DRT peaks. Enables GMM detection.')
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
                          help='mu metric threshold for Lin-KK test (default: 0.85)')
    kk_group.add_argument('--auto-extend', action='store_true',
                          help='Automatically optimize extend_decades for KK validation '
                               '(minimizes pseudo chi-squared)')
    kk_group.add_argument('--extend-decades-max', type=float, default=1.0,
                          help='Maximum extend_decades for --auto-extend search range '
                               '(searches from -max to +max, default: 1.0)')

    # ==========================================================================
    # Z-HIT Validation Group
    # ==========================================================================
    zhit_group = parser.add_argument_group('Z-HIT Validation')

    zhit_group.add_argument('--no-zhit', action='store_true',
                            help='Skip Z-HIT validation')
    zhit_group.add_argument('--zhit-optimize-offset', action='store_true',
                            help='Optimize Z-HIT offset using weighted least-squares '
                                 '(default: fixed reference point)')

    # ==========================================================================
    # Circuit Fitting Group
    # ==========================================================================
    fit_group = parser.add_argument_group('Circuit Fitting')

    fit_group.add_argument('--circuit', '-c', type=str, default=None,
                           help='Equivalent circuit for fitting. '
                                'Syntax: R(100) - (R(5000) | C(1e-6))  [- = series, | = parallel].')
    fit_group.add_argument('--weighting', type=str, default='modulus',
                           choices=['uniform', 'sqrt', 'modulus', 'proportional'],
                           help='Weighting type for fitting (default: modulus)')
    fit_group.add_argument('--no-fit', action='store_true',
                           help='Skip equivalent circuit fitting')
    fit_group.add_argument('--numeric-jacobian', action='store_true',
                           help='Use numeric Jacobian instead of analytic (fallback for custom elements)')

    # Optimizer selection
    fit_group.add_argument('--optimizer', type=str, default='de',
                           choices=['single', 'multistart', 'de'],
                           help='Optimizer: de (differential evolution, default), multistart, or single (one local fit)')

    # Multi-start options
    fit_group.add_argument('--multistart', type=int, default=0, metavar='N',
                           help='Multi-start optimization with N restarts (default: 0 = disabled)')
    fit_group.add_argument('--multistart-scale', type=float, default=2.0,
                           help='Perturbation scaling in sigma units (default: 2.0)')

    # Differential Evolution options
    fit_group.add_argument('--de-strategy', type=int, default=1, choices=[1, 2, 3],
                           help='DE strategy: 1=randtobest1bin (default), 2=best1bin, 3=rand1bin')
    fit_group.add_argument('--de-popsize', type=int, default=15,
                           help='DE population size multiplier (default: 15)')
    fit_group.add_argument('--de-maxiter', type=int, default=1000,
                           help='DE maximum generations (default: 1000)')
    fit_group.add_argument('--de-tol', type=float, default=0.01,
                           help='DE convergence tolerance (default: 0.01)')
    fit_group.add_argument('--de-workers', type=int, default=1,
                           help='DE parallel workers (default: 1, use -1 for all CPUs)')

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
                             help='mu threshold for --voigt-auto-M (default: 0.85)')
    voigt_group.add_argument('--voigt-max-M', type=int, default=50,
                             help='Maximum M elements for --voigt-auto-M (default: 50)')

    # ==========================================================================
    # Oxide Analysis Group
    # ==========================================================================
    oxide_group = parser.add_argument_group('Oxide Layer Analysis')

    oxide_group.add_argument('--analyze-oxide', action='store_true',
                             help='Perform oxide layer analysis')
    oxide_group.add_argument('--epsilon-r', type=float, default=22.0,
                             help='Relative permittivity of oxide (default: 22 for ZrO2)')
    oxide_group.add_argument('--area', type=float, default=1.0,
                             help='Electrode area in cm^2 (default: 1.0)')

    # ==========================================================================
    # Visualization Group
    # ==========================================================================
    vis_group = parser.add_argument_group('Visualization')

    vis_group.add_argument('--ocv', action='store_true',
                           help='Enable OCV (Open Circuit Voltage) curve visualization')

    return parser.parse_args()
