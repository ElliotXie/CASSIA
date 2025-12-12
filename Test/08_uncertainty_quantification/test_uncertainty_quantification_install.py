"""
CASSIA Test 08: Uncertainty Quantification (PIP INSTALL MODE)
==============================================================
Tests the uncertainty quantification functions using pip-installed CASSIA.

Usage:
    python test_uncertainty_quantification_install.py

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

from fixtures import get_cluster_markers, get_full_marker_dataframe
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


def run_uncertainty_quantification_test(results_dir):
    """Test uncertainty quantification functionality using pip-installed CASSIA."""
    print_test_header("08 - Uncertainty Quantification (PIP INSTALL MODE)")

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

    # Test cluster
    test_cluster = "plasma cell"
    markers = get_cluster_markers(test_cluster)

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    uq_results = {}

    try:
        print(f"\n--- Test: runCASSIA_n_times_similarity_score ---")
        print(f"  Cluster: {test_cluster}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
        print(f"  N iterations: 3 (reduced for testing)")

        result = CASSIA.runCASSIA_n_times_similarity_score(
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            additional_info=None,
            temperature=llm_config.get('temperature', 0.3),
            marker_list=markers,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            max_workers=llm_config.get('max_workers', 3),
            n=3,
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

            uq_results = {
                'general_celltype_llm': result.get('general_celltype_llm'),
                'sub_celltype_llm': result.get('sub_celltype_llm'),
                'similarity_score': result.get('similarity_score'),
                'consensus_types': result.get('consensus_types'),
            }

            if result.get('general_celltype_llm') and result.get('similarity_score') is not None:
                status = "passed"
            else:
                status = "failed"
                errors.append("Missing expected result fields for similarity score test")
        else:
            status = "failed"
            errors.append("Unexpected result format for similarity score test")

    except Exception as e:
        errors.append(f"Similarity score test error: {str(e)}")
        status = "error"
        print(f"\nError in similarity score test: {e}")
        import traceback
        traceback.print_exc()

    # =========================================================================
    # Test 2: runCASSIA_batch_n_times - Multiple clusters, multiple iterations
    # =========================================================================
    batch_status = "error"
    batch_results = {}
    batch_clusters_tested = []

    try:
        print(f"\n--- Test: runCASSIA_batch_n_times ---")

        # Get full marker dataframe and limit to 2 clusters for testing
        full_markers = get_full_marker_dataframe()
        test_markers = full_markers.head(2)  # First 2 clusters: monocyte, plasma cell
        batch_clusters_tested = test_markers['Broad.cell.type'].tolist()

        print(f"  Clusters: {batch_clusters_tested}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
        print(f"  N iterations: 2 (reduced for testing)")

        # Set output path for batch results
        batch_output_name = str(results_dir['outputs'] / "batch_results")

        # Run batch n times
        CASSIA.runCASSIA_batch_n_times(
            n=2,  # 2 iterations for testing
            marker=test_markers,
            output_name=batch_output_name,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            additional_info=None,
            celltype_column='Broad.cell.type',
            gene_column_name='Top.Markers',
            max_workers=3,
            batch_max_workers=2,
            provider=llm_config.get('provider', 'openrouter'),
            max_retries=1,
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Check that output files were created
        expected_files = [
            f"{batch_output_name}_1_summary.csv",
            f"{batch_output_name}_1_conversations.json",
            f"{batch_output_name}_2_summary.csv",
            f"{batch_output_name}_2_conversations.json"
        ]

        files_found = []
        files_missing = []
        for f in expected_files:
            if os.path.exists(f):
                files_found.append(os.path.basename(f))
            else:
                files_missing.append(os.path.basename(f))

        print(f"\nBatch Results:")
        print(f"  Files created: {len(files_found)}")
        for f in files_found:
            print(f"    - {f}")

        if files_missing:
            print(f"  Files missing: {len(files_missing)}")
            for f in files_missing:
                print(f"    - {f}")

        batch_results = {
            'clusters_tested': batch_clusters_tested,
            'n_iterations': 2,
            'files_created': files_found,
            'files_missing': files_missing
        }

        # Validate batch results - at least some files should be created
        if len(files_found) >= 2:  # At least 2 files (one iteration with full + summary)
            batch_status = "passed"
            print(f"\n[OK] Batch test PASSED")
        else:
            batch_status = "failed"
            errors.append(f"Batch test: Expected at least 2 output files, found {len(files_found)}")
            print(f"\n[FAIL] Batch test FAILED - insufficient output files")

    except Exception as e:
        errors.append(f"Batch test error: {str(e)}")
        batch_status = "error"
        print(f"\nError in batch test: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Combine statuses - both tests must pass for overall success
    if status == "passed" and batch_status == "passed":
        overall_status = "passed"
    elif status == "error" or batch_status == "error":
        overall_status = "error"
    else:
        overall_status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="uncertainty_quantification_install_py",
        config=config,
        duration_seconds=duration,
        status=overall_status,
        clusters_tested=[test_cluster] + batch_clusters_tested,
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    metadata['test_details'] = {
        'similarity_score_test': {
            'status': status,
            'cluster': test_cluster
        },
        'batch_n_times_test': {
            'status': batch_status,
            'clusters': batch_clusters_tested
        }
    }
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "similarity_score_test": {
            "test_cluster": test_cluster,
            "markers_used": markers[:10] if len(markers) > 10 else markers,
            "n_iterations": 3,
            "results": uq_results
        },
        "batch_n_times_test": batch_results,
        "mode": "pip_install"
    })

    # Print final result
    success = overall_status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("08_uncertainty_quantification")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_uncertainty_quantification_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
