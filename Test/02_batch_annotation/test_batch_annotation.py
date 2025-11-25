"""
CASSIA Test 02: Batch Annotation
================================
Tests the runCASSIA_batch function on all 6 cell clusters.

Usage:
    python test_batch_annotation.py
"""

import sys
import time
import os
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe, get_all_clusters
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
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def run_batch_annotation_test():
    """Run batch annotation test on all clusters."""
    print_test_header("02 - Batch Annotation")

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

    # Get all clusters
    all_clusters = get_all_clusters()
    print(f"\nTesting batch annotation for {len(all_clusters)} clusters:")
    for cluster in all_clusters:
        print(f"  - {cluster}")

    # Get full marker dataframe
    marker_df = get_full_marker_dataframe()
    print(f"\nLoaded marker data: {marker_df.shape}")

    # Create results directory
    results_dir = create_results_dir("02_batch_annotation")
    output_name = str(results_dir / "batch_results")
    print(f"Results will be saved to: {results_dir}")

    # Run the test
    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")
        summary_csv = Path(f"{output_name}_summary.csv")

        if full_csv.exists():
            import pandas as pd
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(all_clusters)}")
            print(f"  Output files created:")
            print(f"    - {full_csv.name}")
            if summary_csv.exists():
                print(f"    - {summary_csv.name}")

            if clusters_annotated == len(all_clusters):
                status = "passed"
            else:
                status = "failed"
                errors.append(f"Only {clusters_annotated}/{len(all_clusters)} clusters annotated")
        else:
            status = "failed"
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")

    duration = time.time() - start_time

    # Save metadata
    metadata = create_test_metadata(
        test_name="batch_annotation",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters,
        errors=errors
    )
    save_test_metadata(results_dir, metadata)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_batch_annotation_test()
    sys.exit(0 if success else 1)
