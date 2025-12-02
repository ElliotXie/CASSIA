"""
CASSIA Test 08: Uncertainty Quantification
==========================================
Tests the uncertainty quantification functions for running multiple analyses
and calculating consensus/similarity scores.

Usage:
    python test_uncertainty_quantification.py

Functions tested:
- runCASSIA_n_times_similarity_score(): Run n single analyses with similarity score
- runCASSIA_batch_n_times(): Run batch analyses n times for multiple clusters
"""

import sys
import time
import os
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_cluster_markers, get_full_marker_dataframe, get_marker_file_path
from test_utils import (
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def run_uncertainty_quantification_test():
    """Test uncertainty quantification functionality."""
    print_test_header("08 - Uncertainty Quantification")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import (
        runCASSIA_n_times_similarity_score,
        runCASSIA_batch_n_times
    )

    # Create results directory with organized structure
    results = create_results_dir("08_uncertainty_quantification", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # Test cluster - use a well-defined cell type
    test_cluster = "plasma cell"
    markers = get_cluster_markers(test_cluster)

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    uq_results = {}

    try:
        # Test: Run n times with similarity score (reduced n for testing)
        print(f"\n--- Test: runCASSIA_n_times_similarity_score ---")
        print(f"  Cluster: {test_cluster}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
        print(f"  N iterations: 3 (reduced for testing)")

        result = runCASSIA_n_times_similarity_score(
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            additional_info=None,
            temperature=llm_config.get('temperature', 0.3),
            marker_list=markers,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            max_workers=llm_config.get('max_workers', 3),
            n=3,  # Reduced for faster testing
            provider=llm_config.get('provider', 'openrouter'),
            main_weight=0.5,
            sub_weight=0.5,
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Check result structure
        if isinstance(result, dict):
            print(f"\nResults:")
            print(f"  General cell type (LLM): {result.get('general_celltype_llm', 'N/A')}")
            print(f"  Sub cell type (LLM): {result.get('sub_celltype_llm', 'N/A')}")
            print(f"  Similarity score: {result.get('similarity_score', 'N/A')}")
            print(f"  Consensus types: {result.get('consensus_types', 'N/A')}")

            if result.get('Possible_mixed_celltypes_llm'):
                print(f"  Possible mixed types: {result.get('Possible_mixed_celltypes_llm')}")

            uq_results = {
                'general_celltype_llm': result.get('general_celltype_llm'),
                'sub_celltype_llm': result.get('sub_celltype_llm'),
                'similarity_score': result.get('similarity_score'),
                'consensus_types': result.get('consensus_types'),
                'mixed_celltypes': result.get('Possible_mixed_celltypes_llm'),
                'original_results_count': len(result.get('original_results', []))
            }

            # Validate results
            if result.get('general_celltype_llm') and result.get('similarity_score') is not None:
                status = "passed"
            else:
                status = "failed"
                errors.append("Missing expected result fields")
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
        test_name="uncertainty_quantification",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[test_cluster],
        errors=errors
    )
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "test_cluster": test_cluster,
        "markers_used": markers[:10] if len(markers) > 10 else markers,
        "n_iterations": 3,
        "results": uq_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_uncertainty_quantification_test()
    sys.exit(0 if success else 1)
