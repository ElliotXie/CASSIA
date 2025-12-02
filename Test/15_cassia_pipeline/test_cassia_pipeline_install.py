"""
CASSIA Test 15: CASSIA Pipeline (PIP INSTALL MODE)
===================================================
Tests the runCASSIA_pipeline function using pip-installed CASSIA.

Usage:
    python test_cassia_pipeline_install.py
"""

import sys
import time
import os
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import load_markers
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
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA


def run_cassia_pipeline_test(results_dir):
    """Run CASSIA pipeline test using pip-installed CASSIA."""
    print_test_header("15 - CASSIA Pipeline (PIP INSTALL MODE)")

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

    # Get settings
    llm_config = config['llm']
    data_config = config['data']

    # Use only 2 clusters for faster testing
    test_clusters = ['monocyte', 'plasma cell']

    # Load full marker data and filter
    full_df = load_markers()
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nTesting pipeline for {len(test_clusters)} clusters:")
    for cluster in test_clusters:
        print(f"  - {cluster}")
    print(f"\nLoaded marker data: {marker_df.shape}")

    # Results directory passed in from main()
    output_name = str(results_dir['outputs'] / "pipeline_test")
    print(f"Results will be saved to: {results_dir['base']}")

    # Change to results directory
    original_dir = os.getcwd()
    os.chdir(results_dir['base'])

    # Run the test
    start_time = time.time()
    errors = []
    status = "error"
    pipeline_output_dir = None

    try:
        print("\nRunning runCASSIA_pipeline (pip install mode)...")
        CASSIA.runCASSIA_pipeline(
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

        # Find the pipeline output directory
        for item in os.listdir(results_dir['base']):
            if item.startswith('CASSIA_') and os.path.isdir(results_dir['base'] / item):
                pipeline_output_dir = results_dir['base'] / item
                break

        if pipeline_output_dir and pipeline_output_dir.exists():
            print(f"\nPipeline output directory: {pipeline_output_dir.name}")

            annotation_dir = pipeline_output_dir / "01_annotation_results"

            checks_passed = True
            print("\nValidating output structure:")

            if annotation_dir.exists():
                print(f"  [OK] 01_annotation_results exists")
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
        os.chdir(original_dir)

    duration = time.time() - start_time

    # Save metadata
    metadata = create_test_metadata(
        test_name="cassia_pipeline_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=test_clusters,
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    # Create results directory first to enable logging
    results_dir = create_results_dir("15_cassia_pipeline")

    # Setup logging to capture all console output
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_cassia_pipeline_test(results_dir)
    finally:
        # Ensure logging is cleaned up even if test fails
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
