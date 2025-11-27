"""
CASSIA Test 15: CASSIA Pipeline
================================
Tests the runCASSIA_pipeline function, the complete end-to-end cell type
annotation orchestrator.

Usage:
    python test_cassia_pipeline.py
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
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    create_test_metadata
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
    from pipeline import runCASSIA_pipeline

    # Use only 2 clusters for faster testing
    test_clusters = ['monocyte', 'plasma cell']

    # Load full marker data and filter to test clusters
    full_df = load_markers()
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nTesting pipeline for {len(test_clusters)} clusters:")
    for cluster in test_clusters:
        print(f"  - {cluster}")
    print(f"\nLoaded marker data: {marker_df.shape}")

    # Create results directory
    results_dir = create_results_dir("15_cassia_pipeline")
    output_name = str(results_dir / "pipeline_test")
    print(f"Results will be saved to: {results_dir}")

    # Change to results directory so pipeline output goes there
    original_dir = os.getcwd()
    os.chdir(results_dir)

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
            score_threshold=75,
            merge_annotations=True,
            merge_model=llm_config.get('model', 'google/gemini-2.5-flash'),
            merge_provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Find the pipeline output directory (starts with CASSIA_)
        for item in os.listdir(results_dir):
            if item.startswith('CASSIA_') and os.path.isdir(results_dir / item):
                pipeline_output_dir = results_dir / item
                break

        if pipeline_output_dir and pipeline_output_dir.exists():
            print(f"\nPipeline output directory: {pipeline_output_dir.name}")

            # Check expected subdirectories
            annotation_dir = pipeline_output_dir / "01_annotation_results"
            reports_dir = pipeline_output_dir / "02_reports"
            boost_dir = pipeline_output_dir / "03_boost_analysis"

            checks_passed = True
            print("\nValidating output structure:")

            # Check 01_annotation_results
            if annotation_dir.exists():
                print(f"  [OK] 01_annotation_results exists")
                # Check for FINAL_RESULTS.csv
                final_results = annotation_dir / "FINAL_RESULTS.csv"
                if final_results.exists():
                    print(f"  [OK] FINAL_RESULTS.csv exists")
                    import pandas as pd
                    results_df = pd.read_csv(final_results)
                    print(f"       - Contains {len(results_df)} rows")
                else:
                    print(f"  [WARN] FINAL_RESULTS.csv not found")
            else:
                print(f"  [FAIL] 01_annotation_results missing")
                checks_passed = False
                errors.append("01_annotation_results directory missing")

            # Check 02_reports
            if reports_dir.exists():
                print(f"  [OK] 02_reports exists")
                html_files = list(reports_dir.glob("*.html"))
                print(f"       - Contains {len(html_files)} HTML report(s)")
            else:
                print(f"  [WARN] 02_reports missing (may be expected if no reports generated)")

            # Check 03_boost_analysis
            if boost_dir.exists():
                print(f"  [OK] 03_boost_analysis exists")
            else:
                print(f"  [INFO] 03_boost_analysis missing (expected if all scores above threshold)")

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
    save_test_metadata(results_dir, metadata)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_cassia_pipeline_test()
    sys.exit(0 if success else 1)
