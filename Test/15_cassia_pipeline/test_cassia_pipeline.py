"""
CASSIA Test 15: CASSIA Pipeline
================================
Tests the runCASSIA_pipeline function, the complete end-to-end cell type
annotation orchestrator.

Usage:
    python D:\CASSIA\Test\15_cassia_pipeline\test_cassia_pipeline.py
"""

import sys
import time
import os
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import load_markers
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
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def run_cassia_pipeline_test():
    """Run CASSIA pipeline test."""
    print_test_header("15 - CASSIA Pipeline")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA pipeline function
    from CASSIA import runCASSIA_pipeline

    # Use only 2 clusters for faster testing
    test_clusters = ['monocyte', 'plasma cell']

    # Load raw Scanpy rank_genes_groups data from data folder for annotation boost
    import pandas as pd
    data_folder = Path(__file__).parent / "data"
    raw_markers_path = data_folder / "scanpy_rank_genes_df.csv"

    if not raw_markers_path.exists():
        raise FileNotFoundError(f"Raw marker file not found: {raw_markers_path}")

    # Load and filter to only the 2 test clusters
    raw_markers = pd.read_csv(raw_markers_path)
    marker_df = raw_markers[raw_markers['group'].isin(test_clusters)].copy()

    print(f"\nTesting pipeline for {len(test_clusters)} clusters:")
    for cluster in test_clusters:
        print(f"  - {cluster}")
    print(f"\nLoaded raw marker data: {marker_df.shape}")
    print(f"Data source: {raw_markers_path}")

    # Create results directory with organized structure
    results = create_results_dir("15_cassia_pipeline", get_test_mode())
    output_name = str(results['outputs'] / "pipeline_test")
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # Change to results directory so pipeline output goes there
    original_dir = os.getcwd()
    os.chdir(results['outputs'])

    # Run the test
    start_time = time.time()
    errors = []
    status = "error"
    pipeline_output_dir = None

    try:
        print("\nRunning runCASSIA_pipeline...")
        runCASSIA_pipeline(
            output_file_name=output_name,
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            marker=marker_df,
            max_workers=llm_config.get('max_workers', 3),
            annotation_model=llm_config.get('model', 'google/gemini-2.5-flash'),
            annotation_provider=llm_config.get('provider', 'openrouter'),
            score_model=llm_config.get('model', 'google/gemini-2.5-flash'),
            score_provider=llm_config.get('provider', 'openrouter'),
            annotationboost_model=llm_config.get('model', 'google/gemini-2.5-flash'),
            annotationboost_provider=llm_config.get('provider', 'openrouter'),
            score_threshold=99,
            merge_annotations=True,
            merge_model=llm_config.get('model', 'google/gemini-2.5-flash'),
            merge_provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Find the pipeline output directory (starts with CASSIA_)
        for item in os.listdir(results['outputs']):
            if item.startswith('CASSIA_') and os.path.isdir(results['outputs'] / item):
                pipeline_output_dir = results['outputs'] / item
                break

        if pipeline_output_dir and pipeline_output_dir.exists():
            print(f"\nPipeline output directory: {pipeline_output_dir.name}")

            # Check expected subdirectories (updated structure)
            annotation_report_dir = pipeline_output_dir / "01_annotation_report"
            boost_dir = pipeline_output_dir / "02_annotation_boost"
            csv_dir = pipeline_output_dir / "03_csv_files"

            checks_passed = True
            print("\nValidating output structure:")

            # Check 01_annotation_report (HTML reports)
            if annotation_report_dir.exists():
                print(f"  [OK] 01_annotation_report exists")
                html_files = list(annotation_report_dir.glob("*.html"))
                print(f"       - Contains {len(html_files)} HTML report(s)")
            else:
                print(f"  [WARN] 01_annotation_report missing")

            # Check 02_annotation_boost
            if boost_dir.exists():
                print(f"  [OK] 02_annotation_boost exists")
                boost_subdirs = [d for d in boost_dir.iterdir() if d.is_dir()]
                print(f"       - Contains {len(boost_subdirs)} cluster folder(s)")
            else:
                print(f"  [INFO] 02_annotation_boost missing (expected if all scores above threshold)")

            # Check 03_csv_files (includes FINAL_RESULTS)
            if csv_dir.exists():
                print(f"  [OK] 03_csv_files exists")
                # Check for FINAL_RESULTS.csv (pattern: *_FINAL_RESULTS.csv)
                final_results_files = list(csv_dir.glob("*_FINAL_RESULTS.csv"))
                if final_results_files:
                    final_results = final_results_files[0]
                    print(f"  [OK] FINAL_RESULTS.csv exists: {final_results.name}")
                    import pandas as pd
                    results_df = pd.read_csv(final_results)
                    print(f"       - Contains {len(results_df)} rows")
                else:
                    print(f"  [WARN] FINAL_RESULTS.csv not found in 03_csv_files")
                csv_files = list(csv_dir.glob("*.csv"))
                print(f"       - Contains {len(csv_files)} CSV file(s)")
            else:
                print(f"  [FAIL] 03_csv_files missing")
                checks_passed = False
                errors.append("03_csv_files directory missing")

            if checks_passed:
                status = "passed"
            else:
                status = "failed"
        else:
            status = "failed"
            errors.append("Pipeline output directory not created")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Change back to original directory
        os.chdir(original_dir)

    duration = time.time() - start_time

    # Save metadata
    metadata = create_test_metadata(
        test_name="cassia_pipeline",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=test_clusters,
        errors=errors
    )
    save_test_metadata(results['outputs'], metadata)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_cassia_pipeline_test()
    sys.exit(0 if success else 1)
