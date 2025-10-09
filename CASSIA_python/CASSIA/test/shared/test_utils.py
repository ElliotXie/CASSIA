"""
Common testing utilities for CASSIA tests.

This module provides utilities for logging, timing, validation, and result saving.
"""

import logging
import time
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional, Union


def setup_logging(
    test_name: str,
    level=logging.INFO,
    log_file: Optional[Path] = None
) -> logging.Logger:
    """
    Setup logging for test with both console and file handlers.

    Args:
        test_name: Name of the test
        level: Logging level (default: INFO)
        log_file: Path to log file (optional)

    Returns:
        Configured logger

    Example:
        >>> logger = setup_logging("test_batch", level=logging.DEBUG)
        >>> logger.info("Test started")
    """
    logger = logging.getLogger(test_name)
    logger.setLevel(level)

    # Remove existing handlers
    logger.handlers = []

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler (if specified)
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(level)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_timestamp(format_string: str = "%Y%m%d_%H%M%S") -> str:
    """
    Get current timestamp as formatted string.

    Args:
        format_string: strftime format string

    Returns:
        Formatted timestamp string

    Example:
        >>> timestamp = get_timestamp()
        >>> print(timestamp)  # 20251007_143022
    """
    return datetime.now().strftime(format_string)


def log_test_start(
    logger: logging.Logger,
    config: Dict[str, Any]
) -> float:
    """
    Log test start information and return start time.

    Args:
        logger: Logger instance
        config: Test configuration dictionary

    Returns:
        Start time (seconds since epoch)

    Example:
        >>> logger = setup_logging("test")
        >>> config = {'test_name': 'batch_test', 'model': 'gemini'}
        >>> start_time = log_test_start(logger, config)
    """
    logger.info("=" * 80)
    logger.info(f"Starting test: {config.get('test_name', 'Unknown')}")
    logger.info(f"Description: {config.get('description', 'N/A')}")
    logger.info(f"Model: {config.get('model', 'N/A')}")
    logger.info(f"Provider: {config.get('provider', 'N/A')}")
    logger.info(f"Data: {config.get('data_file', 'N/A')}")
    logger.info(f"Timestamp: {get_timestamp('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 80)

    return time.time()


def log_test_end(
    logger: logging.Logger,
    start_time: float,
    success: bool = True
):
    """
    Log test completion information.

    Args:
        logger: Logger instance
        start_time: Test start time
        success: Whether test succeeded

    Example:
        >>> log_test_end(logger, start_time, success=True)
    """
    elapsed = time.time() - start_time
    status = "SUCCESS" if success else "FAILED"
    status_symbol = "✓" if success else "✗"

    logger.info("=" * 80)
    logger.info(f"Test Status: {status_symbol} {status}")
    logger.info(f"Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    logger.info(f"Timestamp: {get_timestamp('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 80)


def save_results(
    data: Union[pd.DataFrame, Dict, str, Any],
    test_name: str,
    results_dir: Path,
    prefix: str = "",
    suffix: str = "results"
) -> Path:
    """
    Save test results with timestamp.

    Args:
        data: Data to save (DataFrame, dict, or string)
        test_name: Name of test
        results_dir: Directory to save results
        prefix: Optional prefix for filename
        suffix: Suffix describing result type (default: 'results')

    Returns:
        Path to saved file

    Example:
        >>> df = pd.DataFrame({'a': [1, 2, 3]})
        >>> path = save_results(df, "test_batch", Path("results"))
    """
    # Ensure results directory exists
    results_dir.mkdir(exist_ok=True, parents=True)

    # Generate timestamp
    timestamp = get_timestamp()

    # Build filename
    prefix_str = f"{prefix}_" if prefix else ""
    filename = f"{timestamp}_{prefix_str}{suffix}"

    # Save based on data type
    if isinstance(data, pd.DataFrame):
        output_path = results_dir / f"{filename}.csv"
        data.to_csv(output_path, index=False)

    elif isinstance(data, dict):
        output_path = results_dir / f"{filename}.json"
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)

    elif isinstance(data, str):
        output_path = results_dir / f"{filename}.txt"
        with open(output_path, 'w') as f:
            f.write(data)

    else:
        # Try to convert to string
        output_path = results_dir / f"{filename}.txt"
        with open(output_path, 'w') as f:
            f.write(str(data))

    return output_path


def validate_output(
    data: pd.DataFrame,
    expected_columns: Optional[List[str]] = None,
    min_rows: int = 1,
    max_nulls: Optional[Dict[str, float]] = None
) -> bool:
    """
    Validate output DataFrame structure and content.

    Args:
        data: DataFrame to validate
        expected_columns: List of required column names
        min_rows: Minimum number of rows required
        max_nulls: Dict mapping column names to max null fraction allowed

    Returns:
        True if validation passes

    Raises:
        ValueError: If validation fails

    Example:
        >>> df = pd.DataFrame({'cluster': ['A'], 'annotation': ['T cell']})
        >>> validate_output(df, ['cluster', 'annotation'], min_rows=1)
        True
    """
    # Check it's a DataFrame
    if not isinstance(data, pd.DataFrame):
        raise ValueError(f"Expected DataFrame, got {type(data)}")

    # Check columns
    if expected_columns:
        missing_cols = set(expected_columns) - set(data.columns)
        if missing_cols:
            raise ValueError(
                f"Missing required columns: {', '.join(missing_cols)}"
            )

    # Check row count
    if len(data) < min_rows:
        raise ValueError(
            f"Expected at least {min_rows} rows, got {len(data)}"
        )

    # Check nulls
    if max_nulls:
        for col, max_null_frac in max_nulls.items():
            if col in data.columns:
                null_frac = data[col].isnull().sum() / len(data)
                if null_frac > max_null_frac:
                    raise ValueError(
                        f"Column '{col}' has {null_frac:.1%} nulls, "
                        f"max allowed is {max_null_frac:.1%}"
                    )

    return True


class Timer:
    """
    Context manager for timing code blocks.

    Example:
        >>> logger = setup_logging("test")
        >>> with Timer("Data loading", logger):
        ...     data = load_data()
    """

    def __init__(self, name: str, logger: Optional[logging.Logger] = None):
        """
        Initialize timer.

        Args:
            name: Name of the timed operation
            logger: Optional logger to log timing info
        """
        self.name = name
        self.logger = logger
        self.start_time = None
        self.elapsed = None

    def __enter__(self):
        """Start timer."""
        self.start_time = time.time()
        if self.logger:
            self.logger.info(f"Starting: {self.name}")
        return self

    def __exit__(self, *args):
        """Stop timer and log elapsed time."""
        self.elapsed = time.time() - self.start_time
        if self.logger:
            self.logger.info(
                f"Completed: {self.name} "
                f"({self.elapsed:.2f}s = {self.elapsed/60:.2f}min)"
            )


def create_test_summary(
    test_name: str,
    config: Dict[str, Any],
    results: Dict[str, Any],
    elapsed_time: float,
    success: bool
) -> Dict[str, Any]:
    """
    Create a standardized test summary dictionary.

    Args:
        test_name: Name of the test
        config: Test configuration
        results: Test results
        elapsed_time: Time taken in seconds
        success: Whether test succeeded

    Returns:
        Summary dictionary

    Example:
        >>> summary = create_test_summary(
        ...     "test_batch",
        ...     config,
        ...     {"n_clusters": 6},
        ...     120.5,
        ...     True
        ... )
    """
    return {
        'test_name': test_name,
        'timestamp': get_timestamp('%Y-%m-%d %H:%M:%S'),
        'model': config.get('model'),
        'provider': config.get('provider'),
        'data_file': config.get('data_file'),
        'elapsed_time_seconds': elapsed_time,
        'elapsed_time_minutes': elapsed_time / 60,
        'success': success,
        'results': results
    }


def print_banner(message: str, char: str = "=", width: int = 80):
    """
    Print a banner message.

    Args:
        message: Message to print
        char: Character to use for border
        width: Width of banner

    Example:
        >>> print_banner("Test Starting")
        ================================================================================
        Test Starting
        ================================================================================
    """
    print(char * width)
    print(message)
    print(char * width)
