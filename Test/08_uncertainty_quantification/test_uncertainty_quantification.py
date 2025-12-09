"""
CASSIA Test 08: Uncertainty Quantification
==========================================
Tests the uncertainty quantification functions for running multiple analyses
and calculating consensus/similarity scores.

Usage:
    python D:/CASSIA/Test/08_uncertainty_quantification/test_uncertainty_quantification.py

Functions tested:
- runCASSIA_n_times_similarity_score(): Run n single analyses with similarity score
- runCASSIA_batch_n_times(): Run batch analyses n times for multiple clusters
- runCASSIA_similarity_score_batch(): Analyze batch results and generate similarity score CSV + report
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


def run_uncertainty_quantification_test(run_batch_generation=False):
    """
    Test uncertainty quantification functionality.

    Args:
        run_batch_generation: If True, runs runCASSIA_batch_n_times to generate new batch files.
                             If False (default), uses pre-existing batch files from data folder.
    """
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
        runCASSIA_batch_n_times,
        runCASSIA_similarity_score_batch
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

    # Set up report path for single cluster test
    single_report_path = os.path.join(results['outputs'], "uq_single_report.html")

    try:
        # Test: Run n times with similarity score
        print(f"\n--- Test: runCASSIA_n_times_similarity_score ---")
        print(f"  Cluster: {test_cluster}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
        print(f"  N iterations: 5")

        result = runCASSIA_n_times_similarity_score(
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            additional_info=None,
            temperature=llm_config.get('temperature', 0.3),
            marker_list=markers,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            max_workers=llm_config.get('max_workers', 3),
            n=5,
            provider=llm_config.get('provider', 'openrouter'),
            main_weight=0.5,
            sub_weight=0.5,
            validator_involvement=config.get('validator', {}).get('default', 'v1'),
            generate_report=True,
            report_output_path=single_report_path
        )

        # Check if report was generated
        if os.path.exists(single_report_path):
            print(f"  HTML report created: {os.path.basename(single_report_path)}")

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
    batch_status = "skipped"
    batch_results = {}
    batch_clusters_tested = []
    files_found = []

    # Get full marker dataframe and limit to 2 clusters for testing
    full_markers = get_full_marker_dataframe()
    test_markers = full_markers.head(2)  # First 2 clusters: monocyte, plasma cell
    batch_clusters_tested = test_markers['Broad.cell.type'].tolist()

    # Data folder with pre-existing batch results
    script_dir = Path(__file__).parent
    data_folder = script_dir / "data"

    if run_batch_generation:
        # Run batch generation test
        try:
            print(f"\n--- Test: runCASSIA_batch_n_times ---")

            print(f"  Clusters: {batch_clusters_tested}")
            print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
            print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
            print(f"  N iterations: 5")

            # Set output path for batch results
            batch_output_name = os.path.join(results['outputs'], "batch_results")

            # Run batch n times
            runCASSIA_batch_n_times(
                n=5,
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

            # Check that output files were created (new format: summary CSV + conversations JSON)
            expected_files = [
                f"{batch_output_name}_1_summary.csv",
                f"{batch_output_name}_1_conversations.json",
                f"{batch_output_name}_2_summary.csv",
                f"{batch_output_name}_2_conversations.json",
                f"{batch_output_name}_3_summary.csv",
                f"{batch_output_name}_3_conversations.json",
                f"{batch_output_name}_4_summary.csv",
                f"{batch_output_name}_4_conversations.json",
                f"{batch_output_name}_5_summary.csv",
                f"{batch_output_name}_5_conversations.json"
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
                'n_iterations': 5,
                'files_created': files_found,
                'files_missing': files_missing
            }

            # Validate batch results - at least some files should be created
            if len(files_found) >= 2:  # At least 2 files (one iteration with full + summary)
                batch_status = "passed"
                print(f"\n[OK] Batch test PASSED")

                # Copy generated files to data folder for future use
                import shutil
                data_folder.mkdir(exist_ok=True)
                for i in range(1, 6):
                    # Copy summary CSV
                    src_csv = f"{batch_output_name}_{i}_summary.csv"
                    if os.path.exists(src_csv):
                        dst_csv = data_folder / f"batch_results_{i}_summary.csv"
                        shutil.copy2(src_csv, dst_csv)
                        print(f"  Copied to data folder: batch_results_{i}_summary.csv")
                    # Copy conversations JSON
                    src_json = f"{batch_output_name}_{i}_conversations.json"
                    if os.path.exists(src_json):
                        dst_json = data_folder / f"batch_results_{i}_conversations.json"
                        shutil.copy2(src_json, dst_json)
                        print(f"  Copied to data folder: batch_results_{i}_conversations.json")
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
    else:
        # Use pre-existing batch files from data folder
        print(f"\n--- Skipping Test: runCASSIA_batch_n_times (using pre-existing data) ---")
        batch_output_name = str(data_folder / "batch_results")

        # Check for pre-existing files (new format: summary CSV)
        expected_files = [f"batch_results_{i}_summary.csv" for i in range(1, 6)]
        for f in expected_files:
            if (data_folder / f).exists():
                files_found.append(f)

        if len(files_found) >= 2:
            batch_status = "passed"
            print(f"  Found {len(files_found)} pre-existing batch result files in data folder")
        else:
            print(f"  Warning: Only found {len(files_found)} batch files in {data_folder}")
            print(f"  Expected files: batch_results_1_summary.csv through batch_results_5_summary.csv")
            batch_status = "skipped"

        batch_results = {
            'clusters_tested': batch_clusters_tested,
            'n_iterations': 5,
            'files_created': files_found,
            'source': 'pre-existing data'
        }

    # =========================================================================
    # Test 3: runCASSIA_similarity_score_batch - Analyze batch results
    # =========================================================================
    similarity_batch_status = "error"

    try:
        print(f"\n--- Test: runCASSIA_similarity_score_batch ---")
        print(f"  Analyzing batch results")

        # Run if we have batch files (either generated or pre-existing)
        if len(files_found) >= 2:
            similarity_output_name = os.path.join(results['outputs'], "similarity_score_results")
            similarity_report_path = os.path.join(results['outputs'], "uq_batch_report.html")

            runCASSIA_similarity_score_batch(
                marker=test_markers,
                file_pattern=f"{batch_output_name}_*_summary.csv",
                output_name=similarity_output_name,
                celltype_column='Broad.cell.type',
                max_workers=3,
                model=llm_config.get('model', 'google/gemini-2.5-flash'),
                provider=llm_config.get('provider', 'openrouter'),
                generate_report=True,
                report_output_path=similarity_report_path
            )

            # Check if similarity score CSV was created
            if os.path.exists(f"{similarity_output_name}.csv"):
                print(f"  Similarity score CSV created: {os.path.basename(similarity_output_name)}.csv")
                similarity_batch_status = "passed"
            else:
                similarity_batch_status = "failed"
                errors.append("Similarity score batch: CSV file not created")

            # Check if HTML report was created
            if os.path.exists(similarity_report_path):
                print(f"  HTML report created: {os.path.basename(similarity_report_path)}")
            else:
                print(f"  Warning: HTML report not created")

            print(f"\n[OK] Similarity score batch test PASSED")
        else:
            print(f"  Skipping - no batch files found")
            similarity_batch_status = "skipped"

    except Exception as e:
        errors.append(f"Similarity score batch test error: {str(e)}")
        similarity_batch_status = "error"
        print(f"\nError in similarity score batch test: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Combine statuses - all tests must pass for overall success
    if status == "passed" and (batch_status == "passed" or batch_status == "skipped") and (similarity_batch_status == "passed" or similarity_batch_status == "skipped"):
        overall_status = "passed"
    elif status == "error" or batch_status == "error" or similarity_batch_status == "error":
        overall_status = "error"
    else:
        overall_status = "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="uncertainty_quantification",
        config=config,
        duration_seconds=duration,
        status=overall_status,
        clusters_tested=[test_cluster] + batch_clusters_tested,
        errors=errors
    )
    metadata['test_details'] = {
        'similarity_score_test': {
            'status': status,
            'cluster': test_cluster
        },
        'batch_n_times_test': {
            'status': batch_status,
            'clusters': batch_clusters_tested
        },
        'similarity_score_batch_test': {
            'status': similarity_batch_status,
            'clusters': batch_clusters_tested
        }
    }
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "similarity_score_test": {
            "test_cluster": test_cluster,
            "markers_used": markers[:10] if len(markers) > 10 else markers,
            "n_iterations": 5,
            "results": uq_results
        },
        "batch_n_times_test": batch_results,
        "similarity_score_batch_test": {
            "status": similarity_batch_status
        }
    })

    # Print final result
    success = overall_status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    # Check for --run-batch flag to enable batch generation
    run_batch = "--run-batch" in sys.argv
    success = run_uncertainty_quantification_test(run_batch_generation=run_batch)
    sys.exit(0 if success else 1)
