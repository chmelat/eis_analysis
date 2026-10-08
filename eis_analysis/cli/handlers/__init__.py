"""
Analysis workflow handlers for the EIS CLI.

Each handler corresponds to a step in the analysis pipeline. The handlers are
split across submodules by pipeline stage; this package re-exports the public
``run_*`` functions so ``from eis_analysis.cli.handlers import run_*`` keeps
working.

- validation: run_kk_validation, run_zhit_validation, report_thd,
              apply_zhit_reconstruction, report_outliers, plot_validation
- rinf:       run_rinf_estimation
- drt:        run_drt_analysis, run_voigt_analysis
- fitting:    run_circuit_fitting
- oxide:      run_oxide_analysis
- local_exponent: run_local_exponent
"""

from .validation import (
    run_kk_validation,
    run_zhit_validation,
    report_thd,
    apply_zhit_reconstruction,
    report_outliers,
    plot_validation,
)
from .rinf import run_rinf_estimation
from .drt import run_drt_analysis, run_voigt_analysis
from .fitting import run_circuit_fitting
from .oxide import run_oxide_analysis
from .local_exponent import run_local_exponent

__all__ = [
    'run_kk_validation',
    'run_zhit_validation',
    'report_thd',
    'apply_zhit_reconstruction',
    'report_outliers',
    'plot_validation',
    'run_rinf_estimation',
    'run_drt_analysis',
    'run_voigt_analysis',
    'run_circuit_fitting',
    'run_oxide_analysis',
    'run_local_exponent',
]
