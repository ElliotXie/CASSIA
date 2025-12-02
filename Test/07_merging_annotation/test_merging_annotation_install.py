"""
CASSIA Test 07: Merging Annotation (PIP INSTALL MODE)
======================================================
Tests the merge_annotations and merge_annotations_all functions
using pip-installed CASSIA.

Usage:
    python test_merging_annotation_install.py

Note: This test requires batch annotation results. If none exist, it will
run a batch annotation first.
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe, get_marker_file_path
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
    get_latest_results,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA
import pandas as pd


def run_merging_annotation_test(results_dir):
    """Test merging annotation functionality using pip-installed CASSIA."""
    print_test_header("07 - Merging Annotation (PIP INSTALL MODE)")

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
    llm_config = config['llm']
    data_config = config['data']

    # Results directory passed in from main()
    print(f"Results will be saved to: {results_dir['base']}")

    # First, we need batch results for merging
    batch_results_dir = get_latest_results("02_batch_annotation")
    batch_results_file = None

    if batch_results_dir:
        potential_file = batch_results_dir / "batch_results_full.csv"
        if potential_file.exists():
            batch_results_file = str(potential_file)
            print(f"\nUsing existing batch results: {batch_results_file}")

    # If no existing results, run batch annotation first
    if not batch_results_file:
        print("\nNo existing batch results found. Running batch annotation first...")
        marker_df = get_full_marker_dataframe()
        batch_output = str(results_dir['outputs'] / "batch_for_merge")

        CASSIA.runCASSIA_batch(
            marker=marker_df,
            output_name=batch_output,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )
        batch_results_file = f"{batch_output}_full.csv"

    # Run merging tests
    start_time = time.time()
    errors = []
    status = "error"
    merge_results = {}

    try:
        # Test 1: Single detail level merge (broad)
        print("\n--- Test 1: Single detail level merge (broad) ---")
        broad_output = str(results_dir['outputs'] / "merge_broad.csv")

        result_broad = CASSIA.merge_annotations(
            csv_path=batch_results_file,
            output_path=broad_output,
            provider=llm_config.get('provider', 'openrouter'),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            detail_level="broad",
            batch_size=10
        )

        if result_broad is not None and "Merged_Grouping_1" in result_broad.columns:
            print(f"  Broad merge: SUCCESS")
            print(f"  Output rows: {len(result_broad)}")
            merge_results['broad'] = {
                'status': 'success',
                'output_file': broad_output,
                'num_rows': len(result_broad)
            }
        else:
            errors.append("Broad merge failed - missing Merged_Grouping_1 column")
            merge_results['broad'] = {'status': 'failed'}

        # Test 2: All detail levels merge (parallel)
        print("\n--- Test 2: All detail levels merge (parallel) ---")
        all_output = str(results_dir['outputs'] / "merge_all.csv")

        result_all = CASSIA.merge_annotations_all(
            csv_path=batch_results_file,
            output_path=all_output,
            provider=llm_config.get('provider', 'openrouter'),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            batch_size=10
        )

        expected_columns = ["Merged_Grouping_1", "Merged_Grouping_2", "Merged_Grouping_3"]
        if result_all is not None and all(col in result_all.columns for col in expected_columns):
            print(f"  All levels merge: SUCCESS")
            print(f"  Output rows: {len(result_all)}")
            merge_results['all'] = {
                'status': 'success',
                'output_file': all_output,
                'num_rows': len(result_all),
                'columns': expected_columns
            }
        else:
            missing = [col for col in expected_columns if col not in (result_all.columns if result_all is not None else [])]
            errors.append(f"All levels merge failed - missing columns: {missing}")
            merge_results['all'] = {'status': 'failed'}

        # Determine overall status
        if all(r.get('status') == 'success' for r in merge_results.values()):
            status = "passed"
        elif any(r.get('status') == 'success' for r in merge_results.values()):
            status = "partial"
        else:
            status = "failed"

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="merging_annotation_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=["all"],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "batch_results_file": batch_results_file,
        "merge_results": merge_results,
        "mode": "pip_install"
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("07_merging_annotation")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_merging_annotation_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
