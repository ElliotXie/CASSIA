"""
CASSIA Test 10: Symphony Compare (PIP INSTALL MODE)
=====================================================
Tests the symphonyCompare function using pip-installed CASSIA.

Usage:
    python test_symphony_compare_install.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_cluster_markers
from test_utils import (
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary,
    verify_cassia_pip_install
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA


def run_symphony_compare_test(results_dir):
    """Test Symphony Compare functionality using pip-installed CASSIA."""
    print_test_header("10 - Symphony Compare (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    data_config = config['data']

    # Results directory passed in from main()
    print(f"Results will be saved to: {results_dir['base']}")

    # Test parameters
    test_cluster = "plasma cell"
    markers = get_cluster_markers(test_cluster)
    marker_set = ", ".join(markers[:15])

    celltypes = ["Plasma cell", "B cell", "T cell"]

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    symphony_results = {}

    try:
        print(f"\n--- Test: symphonyCompare ---")
        print(f"  Tissue: {data_config.get('tissue', 'large intestine')}")
        print(f"  Species: {data_config.get('species', 'human')}")
        print(f"  Cell types to compare: {', '.join(celltypes)}")
        print(f"  Model preset: budget")

        result = CASSIA.symphonyCompare(
            tissue=data_config.get('tissue', 'large intestine'),
            celltypes=celltypes,
            marker_set=marker_set,
            species=data_config.get('species', 'human'),
            model_preset="budget",
            output_dir=str(results_dir['outputs']),
            output_basename="symphony_test",
            enable_discussion=False,
            max_discussion_rounds=0,
            consensus_threshold=0.6,
            generate_report=True,
            verbose=True
        )

        # Check result structure
        if isinstance(result, dict):
            print(f"\n--- Symphony Compare Results ---")
            print(f"  Consensus: {result.get('consensus', 'No consensus')}")
            print(f"  Confidence: {result.get('confidence', 0):.1%}")

            symphony_results = {
                'consensus': result.get('consensus'),
                'confidence': result.get('confidence'),
                'csv_file': result.get('csv_file'),
                'html_file': result.get('html_file'),
            }

            if result.get('csv_file') and Path(result['csv_file']).exists():
                status = "passed"
            else:
                status = "failed"
                errors.append("CSV file not created")
        else:
            status = "failed"
            errors.append("Unexpected result format")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="symphony_compare_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=celltypes,
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "celltypes_compared": celltypes,
        "marker_set": marker_set,
        "results": symphony_results,
        "mode": "pip_install"
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("10_symphony_compare")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_symphony_compare_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
