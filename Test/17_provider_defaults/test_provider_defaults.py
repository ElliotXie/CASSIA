"""
CASSIA Test 17: Provider Defaults
==================================
Tests the overall_provider parameter and provider-specific default models
in runCASSIA_pipeline.

Usage:
    python test_provider_defaults.py

Functions tested:
- runCASSIA_pipeline() with overall_provider parameter
- get_pipeline_defaults() for each provider
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
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def run_provider_defaults_test():
    """Test provider-specific default models in CASSIA pipeline."""
    print_test_header("17 - Provider Defaults")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import runCASSIA_pipeline
    from CASSIA.core.model_settings import get_pipeline_defaults

    # Use only 2 clusters for faster testing
    test_clusters = ['monocyte', 'plasma cell']

    # Load raw Scanpy rank_genes_groups data
    import pandas as pd
    data_folder = Path(__file__).parent.parent / "16_cassia_pipeline" / "data"
    raw_markers_path = data_folder / "scanpy_rank_genes_df.csv"

    if not raw_markers_path.exists():
        raise FileNotFoundError(f"Raw marker file not found: {raw_markers_path}")

    # Load and filter to only the 2 test clusters
    raw_markers = pd.read_csv(raw_markers_path)
    marker_df = raw_markers[raw_markers['group'].isin(test_clusters)].copy()

    print(f"\nTesting pipeline with overall_provider parameter")
    print(f"Testing for {len(test_clusters)} clusters:")
    for cluster in test_clusters:
        print(f"  - {cluster}")
    print(f"\nLoaded raw marker data: {marker_df.shape}")

    # Create results directory with organized structure
    results = create_results_dir("17_provider_defaults", get_test_mode())
    output_name = str(results['outputs'] / "provider_defaults_test")
    print(f"Results will be saved to: {results['base']}")

    # Setup logging
    logging_ctx = setup_logging(results['logs'])

    # Test overall_provider parameter - using openrouter as specified
    test_provider = "openrouter"

    # Print expected defaults for the provider
    defaults = get_pipeline_defaults(test_provider)
    print(f"\n--- Testing overall_provider='{test_provider}' ---")
    print(f"Expected defaults:")
    print(f"  Annotation:       {defaults.get('annotation')}")
    print(f"  Score:            {defaults.get('score')}")
    print(f"  Merge:            {defaults.get('merge')}")
    print(f"  AnnotationBoost:  {defaults.get('annotationboost')}")

    # Change to results directory
    original_dir = os.getcwd()
    os.chdir(results['outputs'])

    # Run the test
    start_time = time.time()
    errors = []
    status = "error"
    pipeline_output_dir = None

    try:
        print("\nRunning runCASSIA_pipeline with overall_provider...")
        runCASSIA_pipeline(
            output_file_name=output_name,
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            marker=marker_df,
            max_workers=llm_config.get('max_workers', 3),
            overall_provider=test_provider,  # Key parameter being tested
            # All model params are None - should use provider defaults
            annotation_model=None,
            annotation_provider=None,
            score_model=None,
            score_provider=None,
            annotationboost_model=None,
            annotationboost_provider=None,
            merge_model=None,
            merge_provider=None,
            score_threshold=99,
            merge_annotations=True,
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )

        # Find the pipeline output directory
        for item in os.listdir(results['outputs']):
            if item.startswith('CASSIA_') and os.path.isdir(results['outputs'] / item):
                pipeline_output_dir = results['outputs'] / item
                break

        if pipeline_output_dir and pipeline_output_dir.exists():
            print(f"\nPipeline output directory: {pipeline_output_dir.name}")

            # Check expected subdirectories
            csv_dir = pipeline_output_dir / "03_csv_files"

            checks_passed = True
            print("\nValidating output structure:")

            if csv_dir.exists():
                print(f"  [OK] 03_csv_files exists")
                final_results_files = list(csv_dir.glob("*_FINAL_RESULTS.csv"))
                if final_results_files:
                    final_results = final_results_files[0]
                    print(f"  [OK] FINAL_RESULTS.csv exists: {final_results.name}")
                    results_df = pd.read_csv(final_results)
                    print(f"       - Contains {len(results_df)} rows")
                else:
                    print(f"  [WARN] FINAL_RESULTS.csv not found")
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
        os.chdir(original_dir)

    duration = time.time() - start_time

    # Save metadata
    metadata = create_test_metadata(
        test_name="provider_defaults",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=test_clusters,
        errors=errors
    )
    metadata['test_provider'] = test_provider
    metadata['expected_defaults'] = defaults
    save_test_metadata(results['outputs'], metadata)

    # Save results
    save_test_results(results['outputs'], {
        'overall_provider': test_provider,
        'expected_defaults': defaults,
        'clusters_tested': test_clusters
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_provider_defaults_test()
    sys.exit(0 if success else 1)
