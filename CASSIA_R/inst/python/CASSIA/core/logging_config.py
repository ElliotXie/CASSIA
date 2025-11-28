"""
Centralized logging configuration for CASSIA package.
Errors are VISIBLE by default with actionable messages.
"""
import logging
import warnings
import sys
from typing import Optional

# Package-level logger instance
_logger: Optional[logging.Logger] = None


def get_logger(name: str = "CASSIA") -> logging.Logger:
    """
    Get or create the CASSIA logger.

    Args:
        name: Logger name. Use __name__ for module-specific loggers.
              Will be prefixed with "CASSIA." automatically.

    Returns:
        logging.Logger: Configured logger instance.

    Example:
        >>> logger = get_logger(__name__)
        >>> logger.warning("Something went wrong")
        [CASSIA WARNING] Something went wrong
    """
    global _logger
    if _logger is None:
        _logger = logging.getLogger("CASSIA")
        if not _logger.handlers:
            handler = logging.StreamHandler(sys.stderr)
            formatter = logging.Formatter(
                '[CASSIA %(levelname)s] %(message)s'
            )
            handler.setFormatter(formatter)
            _logger.addHandler(handler)
            _logger.setLevel(logging.INFO)  # VISIBLE by default

    # Return child logger if name provided, otherwise root CASSIA logger
    if name and name != "CASSIA":
        return logging.getLogger(f"CASSIA.{name}")
    return _logger


def set_log_level(level: str = "INFO"):
    """
    Set CASSIA log level to control message visibility.

    Args:
        level: One of "DEBUG", "INFO", "WARNING", "ERROR", or "QUIET"
               - DEBUG: Show all messages including detailed debug info
               - INFO: Show informational messages, warnings, and errors (default)
               - WARNING: Show only warnings and errors
               - ERROR: Show only errors
               - QUIET: Suppress almost all messages

    Example:
        >>> import CASSIA
        >>> CASSIA.set_log_level("DEBUG")  # See everything
        >>> CASSIA.set_log_level("QUIET")  # Suppress messages
    """
    levels = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
        "QUIET": logging.CRITICAL  # Suppress almost everything
    }
    logger = get_logger()
    log_level = levels.get(level.upper(), logging.INFO)
    logger.setLevel(log_level)

    # Also set level on all child loggers
    for handler in logger.handlers:
        handler.setLevel(log_level)


def warn_user(message: str):
    """
    Show a warning to the user that is always visible.
    Uses Python's warnings module to ensure visibility.

    Args:
        message: Warning message to display.
    """
    warnings.warn(f"[CASSIA] {message}", UserWarning, stacklevel=2)
