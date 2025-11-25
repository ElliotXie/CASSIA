"""
CASSIA Test Suite - Results Manager
===================================
Functions for managing test results and outputs.
"""

import os
import json
import shutil
from pathlib import Path
from datetime import datetime


def get_test_root() -> Path:
    """Get the root directory of the test suite."""
    return Path(__file__).parent.parent.parent


def create_results_dir(test_folder: str) -> Path:
    """
    Create a timestamped results directory for a test.

    Args:
        test_folder: Name of the test folder (e.g., '01_single_annotation')

    Returns:
        Path: Path to the created results directory
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_base = get_test_root() / test_folder / "results"
    results_dir = results_base / timestamp

    results_dir.mkdir(parents=True, exist_ok=True)
    return results_dir


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
