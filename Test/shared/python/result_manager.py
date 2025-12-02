"""
CASSIA Test Suite - Results Manager
===================================
Functions for managing test results and outputs.
"""

import os
import sys
import json
import shutil
from pathlib import Path
from datetime import datetime


class TestLogger:
    """
    Tee-style logger that captures output to both console and log file.

    This class redirects stdout or stderr to write to both the terminal
    and a log file simultaneously, ensuring all console output is captured.
    """

    def __init__(self, log_path: Path, stream_type: str = 'stdout'):
        """
        Initialize the TestLogger.

        Args:
            log_path: Path to the log file
            stream_type: 'stdout' or 'stderr'
        """
        self.log_path = log_path
        self.stream_type = stream_type
        self.terminal = sys.stdout if stream_type == 'stdout' else sys.stderr
        self.log_file = open(log_path, 'a', encoding='utf-8')

        # Write header to log file
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log_file.write(f"\n{'='*60}\n")
        self.log_file.write(f"Log started: {timestamp}\n")
        self.log_file.write(f"Stream: {stream_type}\n")
        self.log_file.write(f"{'='*60}\n\n")
        self.log_file.flush()

    def write(self, message):
        """Write message to both terminal and log file."""
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()

    def flush(self):
        """Flush both terminal and log file."""
        self.terminal.flush()
        self.log_file.flush()

    def close(self):
        """Close the log file and write footer."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log_file.write(f"\n{'='*60}\n")
        self.log_file.write(f"Log ended: {timestamp}\n")
        self.log_file.write(f"{'='*60}\n")
        self.log_file.close()

    def fileno(self):
        """Return the file descriptor of the terminal (for compatibility)."""
        return self.terminal.fileno()

    def isatty(self):
        """Check if the terminal is a tty (for compatibility)."""
        return self.terminal.isatty()


def setup_logging(logs_dir: Path, log_filename: str = "test_log.txt") -> dict:
    """
    Setup logging to capture stdout/stderr to a log file.

    The log file will be created in the logs directory with a unique
    timestamped name if the default already exists.

    Args:
        logs_dir: Path to the logs directory (e.g., results['logs'])
        log_filename: Name of the log file (default: test_log.txt)

    Returns:
        dict: Contains loggers and original streams for cleanup
    """
    log_path = logs_dir / log_filename

    # Ensure unique filename if file already exists
    if log_path.exists():
        base = log_path.stem
        ext = log_path.suffix
        counter = 1
        while log_path.exists():
            log_path = results_dir / f"{base}_{counter}{ext}"
            counter += 1

    # Store original streams
    original_stdout = sys.stdout
    original_stderr = sys.stderr

    # Create loggers
    stdout_logger = TestLogger(log_path, 'stdout')
    stderr_logger = TestLogger(log_path, 'stderr')

    # Redirect streams
    sys.stdout = stdout_logger
    sys.stderr = stderr_logger

    print(f"Console output logging enabled: {log_path}")

    return {
        'stdout_logger': stdout_logger,
        'stderr_logger': stderr_logger,
        'original_stdout': original_stdout,
        'original_stderr': original_stderr,
        'log_path': log_path
    }


def cleanup_logging(logging_context: dict):
    """
    Restore original stdout/stderr and close log files.

    Args:
        logging_context: Dictionary returned by setup_logging()
    """
    if logging_context is None:
        return

    # Get components from context
    stdout_logger = logging_context.get('stdout_logger')
    stderr_logger = logging_context.get('stderr_logger')
    original_stdout = logging_context.get('original_stdout')
    original_stderr = logging_context.get('original_stderr')
    log_path = logging_context.get('log_path')

    # Print final message before closing
    if log_path:
        print(f"\nConsole output saved to: {log_path}")

    # Restore original streams
    if original_stdout:
        sys.stdout = original_stdout
    if original_stderr:
        sys.stderr = original_stderr

    # Close loggers
    if stdout_logger:
        stdout_logger.close()
    if stderr_logger:
        stderr_logger.close()


class LoggingContext:
    """
    Context manager for console output logging.

    Usage:
        with LoggingContext(logs_dir) as log_ctx:
            # All print statements will be captured
            print("This goes to both console and log file")
    """

    def __init__(self, logs_dir: Path, log_filename: str = "test_log.txt"):
        self.logs_dir = logs_dir
        self.log_filename = log_filename
        self.logging_context = None

    def __enter__(self):
        self.logging_context = setup_logging(self.logs_dir, self.log_filename)
        return self.logging_context

    def __exit__(self, exc_type, exc_val, exc_tb):
        cleanup_logging(self.logging_context)
        return False  # Don't suppress exceptions


def get_test_root() -> Path:
    """Get the root directory of the test suite."""
    return Path(__file__).parent.parent.parent


def create_results_dir(test_folder: str, mode: str = 'development') -> dict:
    """
    Create a timestamped results directory for a test with organized subfolders.

    Args:
        test_folder: Name of the test folder (e.g., '01_single_annotation')
        mode: 'installed' or 'development' (default: 'development')

    Returns:
        dict: Dictionary with paths:
            - 'base': Path to the timestamp directory
            - 'logs': Path to logs/ subdirectory
            - 'outputs': Path to outputs/ subdirectory

    Directory structure:
        results/python/{mode}/{timestamp}/outputs/
        results/python/{mode}/{timestamp}/logs/
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_base = get_test_root() / test_folder / "results"
    timestamp_dir = results_base / "python" / mode / timestamp

    # Create organized subdirectories
    logs_dir = timestamp_dir / "logs"
    outputs_dir = timestamp_dir / "outputs"

    logs_dir.mkdir(parents=True, exist_ok=True)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    return {
        'base': timestamp_dir,
        'logs': logs_dir,
        'outputs': outputs_dir
    }


def save_test_metadata(results_dir: Path, metadata: dict):
    """
    Save test metadata to JSON file.

    Args:
        results_dir: Path to results directory
        metadata: Dictionary containing test metadata
    """
    metadata_path = results_dir / "test_metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2, default=str)


def save_test_results(results_dir: Path, results: dict, filename: str = "results.json"):
    """
    Save test results to JSON file.

    Args:
        results_dir: Path to results directory
        results: Dictionary containing test results
        filename: Name of the output file
    """
    results_path = results_dir / filename
    with open(results_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)


def create_test_metadata(
    test_name: str,
    config: dict,
    duration_seconds: float,
    status: str,
    clusters_tested: list = None,
    errors: list = None
) -> dict:
    """
    Create a test metadata dictionary.

    Args:
        test_name: Name of the test
        config: Test configuration
        duration_seconds: Test duration in seconds
        status: Test status ('passed', 'failed', 'error')
        clusters_tested: List of clusters tested
        errors: List of error messages

    Returns:
        dict: Test metadata
    """
    return {
        "test_name": test_name,
        "timestamp": datetime.now().isoformat(),
        "language": "python",
        "config": {
            "model": config.get('llm', {}).get('model'),
            "provider": config.get('llm', {}).get('provider'),
            "tissue": config.get('data', {}).get('tissue'),
            "species": config.get('data', {}).get('species'),
            "n_genes": config.get('data', {}).get('n_genes'),
        },
        "duration_seconds": round(duration_seconds, 2),
        "status": status,
        "clusters_tested": clusters_tested or [],
        "errors": errors or []
    }


def cleanup_old_results(test_folder: str, keep_n: int = 10):
    """
    Remove old results directories, keeping only the most recent N.

    Args:
        test_folder: Name of the test folder
        keep_n: Number of recent results to keep
    """
    results_base = get_test_root() / test_folder / "results"
    if not results_base.exists():
        return

    # Get all timestamped directories
    result_dirs = sorted(
        [d for d in results_base.iterdir() if d.is_dir()],
        key=lambda x: x.name,
        reverse=True
    )

    # Remove old directories
    for old_dir in result_dirs[keep_n:]:
        shutil.rmtree(old_dir)
        print(f"Removed old results: {old_dir.name}")


def get_latest_results(test_folder: str) -> Path:
    """
    Get the most recent results directory for a test.

    Args:
        test_folder: Name of the test folder

    Returns:
        Path: Path to the most recent results directory, or None
    """
    results_base = get_test_root() / test_folder / "results"
    if not results_base.exists():
        return None

    result_dirs = sorted(
        [d for d in results_base.iterdir() if d.is_dir()],
        key=lambda x: x.name,
        reverse=True
    )

    return result_dirs[0] if result_dirs else None


def list_all_results(test_folder: str) -> list:
    """
    List all results directories for a test.

    Args:
        test_folder: Name of the test folder

    Returns:
        list: List of results directory paths
    """
    results_base = get_test_root() / test_folder / "results"
    if not results_base.exists():
        return []

    return sorted(
        [d for d in results_base.iterdir() if d.is_dir()],
        key=lambda x: x.name,
        reverse=True
    )
