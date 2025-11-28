"""
CASSIA Test 01: Single Cluster Annotation
==========================================
Tests the runCASSIA function on a single cell cluster.

Usage:
    python test_single_annotation.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_marker_dataframe_for_cluster, get_all_clusters
from test_utils import (
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    validate_annotation_result,
    print_test_header,
    print_test_result,
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def run_single_annotation_test():
    """Run single cluster annotation test."""
    print_test_header("01 - Single Cluster Annotation")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import runCASSIA

    # Test cluster
    test_cluster = "plasma cell"
    n_genes = data_config.get('n_genes', 30)

    print(f"\nTesting single annotation for: {test_cluster}")
    print(f"Using top {n_genes} markers")

    # Get markers for the test cluster
    marker_df = get_marker_dataframe_for_cluster(test_cluster, n_genes)
    # Extract marker list from DataFrame (runCASSIA expects list of strings, not DataFrame)
    marker_list = marker_df['gene'].tolist()
    print(f"Loaded {len(marker_list)} markers")

    # Create results directory
    results_dir = create_results_dir("01_single_annotation")
    print(f"Results will be saved to: {results_dir}")

    # Run the test
    start_time = time.time()
    errors = []
    result = None

    try:
        print("\nRunning runCASSIA...")
        result, conversation_history, _ = runCASSIA(
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
        test_name="single_annotation",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[test_cluster],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)

    if result:
        save_test_results(results_dir, {
            "cluster": test_cluster,
            "result": result,
            "validation": validation
        })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_single_annotation_test()
    sys.exit(0 if success else 1)
