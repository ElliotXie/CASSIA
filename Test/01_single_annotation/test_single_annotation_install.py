"""
CASSIA Test 01: Single Cluster Annotation (PIP INSTALL MODE)
=============================================================
Tests the runCASSIA function on a single cell cluster using pip-installed CASSIA.

Usage:
    python test_single_annotation_install.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_marker_dataframe_for_cluster
from test_utils import (
    load_config,
    setup_api_keys,
    validate_annotation_result,
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


def run_single_annotation_test(results_dir):
    """Run single cluster annotation test using pip-installed CASSIA."""
    print_test_header("01 - Single Cluster Annotation (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Location: {pip_info['location']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    if not pip_info['is_pip_install']:
        print("\nWARNING: CASSIA may not be from pip install!")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Test cluster
    test_cluster = "plasma cell"
    n_genes = data_config.get('n_genes', 30)

    print(f"\nTesting single annotation for: {test_cluster}")
    print(f"Using top {n_genes} markers")

    # Get markers for the test cluster
    marker_df = get_marker_dataframe_for_cluster(test_cluster, n_genes)
    marker_list = marker_df['gene'].tolist()
    print(f"Loaded {len(marker_list)} markers")

    # Results directory passed in from main()
    print(f"Results will be saved to: {results_dir['base']}")

    # Run the test
    start_time = time.time()
    errors = []
    result = None
    validation = None

    try:
        print("\nRunning runCASSIA (pip install mode)...")
        result, conversation_history, _ = CASSIA.runCASSIA(
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            marker_list=marker_list,
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Validate result
        validation = validate_annotation_result(result)

        if validation['valid']:
            print(f"\nAnnotation Result:")
            print(f"  Main cell type: {result.get('main_cell_type')}")
            print(f"  Sub cell types: {result.get('sub_cell_types')}")
            status = "passed"
        else:
            errors = validation['errors']
            status = "failed"
            print(f"\nValidation errors: {errors}")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")

    duration = time.time() - start_time

    # Save results
    metadata = create_test_metadata(
        test_name="single_annotation_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[test_cluster],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    if result:
        save_test_results(results_dir['base'], {
            "cluster": test_cluster,
            "result": result,
            "validation": validation,
            "mode": "pip_install"
        })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("01_single_annotation")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_single_annotation_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
