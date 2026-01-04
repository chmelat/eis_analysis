"""
EIS Analysis Toolkit
====================

Modular toolkit for electrochemical impedance spectroscopy (EIS) analysis
with focus on Distribution of Relaxation Times (DRT).

Modules:
- io: Data loading and synthetic data generation
- validation: Kramers-Kronig validation
- rinf_estimation: R_inf estimation with inductance compensation
- drt: Distribution of Relaxation Times analysis
- fitting: Equivalent circuit fitting with operator overloading
- analysis: Oxide layer analysis
- visualization: Plotting and diagnostics

Version is imported from eis_analysis.version (single source of truth).
"""

# Import version from single source of truth
from .version import __version__, __version_info__, get_version_string

# I/O
from .io import (
    load_data,
    load_csv_data,
    parse_dta_metadata,
    parse_ocv_curve,
    log_metadata,
    generate_synthetic_data,
)

# Validation
from .validation import (
    kramers_kronig_validation,
)

# R_inf estimation (shared utility for DRT and circuit fitting)
from .rinf_estimation import (
    estimate_rinf_with_inductance,
)

# DRT Analysis
from .drt import (
    calculate_drt,
    compute_gcv_score,
    find_optimal_lambda_gcv,
    gmm_peak_detection,
    GMM_AVAILABLE,
)

# Fitting (new operator overloading approach)
from .fitting import (
    # Circuit elements
    R, C, Q, L, W, Wo,
    # Main functions
    fit_equivalent_circuit,
    fit_circuit_multistart,  # Multi-start optimization
    fit_circuit_diffevo,  # Differential evolution optimization
    analyze_voigt_elements,
    format_voigt_report,
    FitResult,
    MultistartResult,
    DiffEvoResult,
    # Voigt chain linear fitting
    fit_voigt_chain_linear,
)

# Analysis
from .analysis import (
    analyze_oxide_layer,
)

# Visualization
from .visualization import (
    visualize_data,
    visualize_ocv,
    plot_rl_fit_diagnostics,
)

__all__ = [
    # Version info
    '__version__',
    '__version_info__',
    'get_version_string',
    # I/O
    'load_data',
    'load_csv_data',
    'parse_dta_metadata',
    'parse_ocv_curve',
    'log_metadata',
    'generate_synthetic_data',
    # Validation
    'kramers_kronig_validation',
    # DRT
    'calculate_drt',
    'compute_gcv_score',
    'find_optimal_lambda_gcv',
    'estimate_rinf_with_inductance',
    'gmm_peak_detection',
    'GMM_AVAILABLE',
    # Fitting (circuit elements)
    'R', 'C', 'Q', 'L', 'W', 'Wo',
    # Fitting (main functions)
    'fit_equivalent_circuit',
    'fit_circuit_multistart',
    'fit_circuit_diffevo',
    'analyze_voigt_elements',
    'format_voigt_report',
    'FitResult',
    'MultistartResult',
    'DiffEvoResult',
    'fit_voigt_chain_linear',
    # Analysis
    'analyze_oxide_layer',
    # Visualization
    'visualize_data',
    'visualize_ocv',
    'plot_rl_fit_diagnostics',
]
