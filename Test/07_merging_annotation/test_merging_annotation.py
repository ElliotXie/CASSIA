"""
CASSIA Test 07: Merging Annotation
==================================
Tests the merge_annotations and merge_annotations_all functions for grouping
cell type annotations at different detail levels.

Usage:
    python test_merging_annotation.py

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
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    get_latest_results
)

# Setup CASSIA imports
setup_cassia_imports()


def run_merging_annotation_test():
    """Test merging annotation functionality."""
    print_test_header("07 - Merging Annotation")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from tools_function import runCASSIA_batch
    from merging_annotation import merge_annotations, merge_annotations_all
    import pandas as pd

    # Create results directory
    results_dir = create_results_dir("07_merging_annotation")
    print(f"Results will be saved to: {results_dir}")

    # First, we need batch results for merging
    # Check if we have existing batch results
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
        batch_output = str(results_dir / "batch_for_merge")

        runCASSIA_batch(
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
        broad_output = str(results_dir / "merge_broad.csv")

        result_broad = merge_annotations(
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
            print(f"  Sample groupings:")
            for idx, row in result_broad.head(3).iterrows():
                print(f"    {row['True Cell Type']} -> {row['Merged_Grouping_1']}")
            merge_results['broad'] = {
                'status': 'success',
                'output_file': broad_output,
                'num_rows': len(result_broad)
            }
        else:
            errors.append("Broad merge failed - missing Merged_Grouping_1 column")
            merge_results['broad'] = {'status': 'failed'}

        # Test 2: Single detail level merge (detailed)
        print("\n--- Test 2: Single detail level merge (detailed) ---")
        detailed_output = str(results_dir / "merge_detailed.csv")

        result_detailed = merge_annotations(
            csv_path=batch_results_file,
            output_path=detailed_output,
            provider=llm_config.get('provider', 'openrouter'),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            detail_level="detailed",
            batch_size=10
        )

        if result_detailed is not None and "Merged_Grouping_2" in result_detailed.columns:
            print(f"  Detailed merge: SUCCESS")
            print(f"  Output rows: {len(result_detailed)}")
            print(f"  Sample groupings:")
            for idx, row in result_detailed.head(3).iterrows():
                print(f"    {row['True Cell Type']} -> {row['Merged_Grouping_2']}")
            merge_results['detailed'] = {
                'status': 'success',
                'output_file': detailed_output,
                'num_rows': len(result_detailed)
            }
        else:
            errors.append("Detailed merge failed - missing Merged_Grouping_2 column")
            merge_results['detailed'] = {'status': 'failed'}

        # Test 3: All detail levels merge (parallel)
        print("\n--- Test 3: All detail levels merge (parallel) ---")
        all_output = str(results_dir / "merge_all.csv")

        result_all = merge_annotations_all(
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
            print(f"  Columns: {', '.join(expected_columns)}")
            print(f"\n  Sample comparison:")
            for idx, row in result_all.head(3).iterrows():
                print(f"    {row['True Cell Type']}:")
                print(f"      Broad:        {row['Merged_Grouping_1']}")
                print(f"      Detailed:     {row['Merged_Grouping_2']}")
                print(f"      Very Detailed: {row['Merged_Grouping_3']}")
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
        test_name="merging_annotation",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=["all"],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)

    save_test_results(results_dir, {
        "batch_results_file": batch_results_file,
        "merge_results": merge_results,
        "detail_levels_tested": ["broad", "detailed", "very_detailed", "all"]
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_merging_annotation_test()
    sys.exit(0 if success else 1)
