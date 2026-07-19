"""
Custom logging configuration for EIS CLI.

Provides formatters and handlers for clean CLI output:
- INFO: no prefix (clean output)
- WARNING: "! " prefix
- ERROR: "!! " prefix
- DEBUG: "[DEBUG] " prefix
"""

import argparse
import logging
import sys

logger = logging.getLogger(__name__)


# =============================================================================
# Custom Formatters
# =============================================================================

class PrefixFormatter(logging.Formatter):
    """Level-specific prefix formatter - INFO stays clean (no prefix)."""

    PREFIXES = {
        logging.WARNING: "! ",
        logging.ERROR: "!! ",
        logging.CRITICAL: "!! ",
        logging.DEBUG: "[DEBUG] ",
    }

    def format(self, record):
        return self.PREFIXES.get(record.levelno, "") + record.getMessage()


class LevelFilter(logging.Filter):
    """Filter that accepts only specific log levels."""

    def __init__(self, levels):
        super().__init__()
        self.levels = levels if isinstance(levels, (list, tuple)) else [levels]

    def filter(self, record):
        return record.levelno in self.levels


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

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command line arguments with 'quiet' and 'verbose' attributes
    """
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
        info_handler.setFormatter(PrefixFormatter())
        root_logger.addHandler(info_handler)

    # Handler for WARNING -> stdout (prefix "!")
    warning_handler = logging.StreamHandler(sys.stdout)
    warning_handler.setLevel(logging.WARNING)
    warning_handler.addFilter(LevelFilter(logging.WARNING))
    warning_handler.setFormatter(PrefixFormatter())
    root_logger.addHandler(warning_handler)

    # Handler for ERROR/CRITICAL -> stderr (prefix "!!")
    error_handler = logging.StreamHandler(sys.stderr)
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(PrefixFormatter())
    root_logger.addHandler(error_handler)

    # Handler for DEBUG -> stderr (prefix "[DEBUG]")
    if show_debug:
        debug_handler = logging.StreamHandler(sys.stderr)
        debug_handler.setLevel(logging.DEBUG)
        debug_handler.addFilter(LevelFilter(logging.DEBUG))
        debug_handler.setFormatter(PrefixFormatter())
        root_logger.addHandler(debug_handler)


def log_separator(length: int = 50, char: str = "=") -> None:
    """
    Log a separator line for visual clarity.

    Parameters
    ----------
    length : int
        Length of separator line (default: 50)
    char : str
        Character to use for separator (default: "=")
    """
    logging.getLogger(__name__).info(char * length)
