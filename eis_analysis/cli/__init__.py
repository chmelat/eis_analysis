"""
CLI module for EIS Analysis Toolkit.

This module provides the command-line interface components:
- logging: Custom log formatters and setup
- parser: Argument parsing
- data_handling: Data loading and filtering
- handlers: Analysis workflow handlers
- utils: Helper functions and dataclasses

The main entry point is in the root eis.py script.
"""

from .logging import setup_logging, log_separator
from .parser import parse_arguments
from .data_handling import load_eis_data, filter_by_frequency
from .handlers import (
    run_kk_validation,
    run_zhit_validation,
    run_rinf_estimation,
    run_drt_analysis,
    run_voigt_analysis,
    run_circuit_fitting,
    run_oxide_analysis,
)
from .utils import (
    EISAnalysisError,
    LoadedData,
    CLIDRTResult,
    save_figure,
    parse_circuit_expression,
)

__all__ = [
    # Logging
    'setup_logging',
    'log_separator',
    # Parser
    'parse_arguments',
    # Data handling
    'load_eis_data',
    'filter_by_frequency',
    # Handlers
    'run_kk_validation',
    'run_zhit_validation',
    'run_rinf_estimation',
    'run_drt_analysis',
    'run_voigt_analysis',
    'run_circuit_fitting',
    'run_oxide_analysis',
    # Utils
    'EISAnalysisError',
    'LoadedData',
    'CLIDRTResult',
    'save_figure',
    'parse_circuit_expression',
]
