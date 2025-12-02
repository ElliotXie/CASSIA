"""
CASSIA Test 06: Annotation Boost (PIP INSTALL MODE)
====================================================
Tests the runCASSIA_annotationboost function using pip-installed CASSIA.

Usage:
    python test_annotation_boost_install.py

Note: This test requires batch annotation results. If none exist, it will
run a batch annotation first.
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe


def get_local_marker_file_path():
    """Get the path to the local FindAllMarkers output file."""
    return Path(__file__).parent / "data" / "findallmarkers_output.csv"


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
    get_latest_results,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA
import pandas as pd


def run_annotation_boost_test(results_dir):
    """Test annotation boost functionality using pip-installed CASSIA."""
    print_test_header("06 - Annotation Boost (PIP INSTALL MODE)")

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

    # First, we need batch results for annotation boost
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
        batch_output = str(results_dir['outputs'] / "batch_for_boost")

        CASSIA.runCASSIA_batch(
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

    # Load the batch results to get cluster names
    batch_df = pd.read_csv(batch_results_file)
    cluster_col = 'Cluster ID' if 'Cluster ID' in batch_df.columns else ('True Cell Type' if 'True Cell Type' in batch_df.columns else batch_df.columns[0])
    available_clusters = batch_df[cluster_col].tolist()

    # Test cluster for annotation boost
    test_cluster = "plasma cell"
    if test_cluster not in available_clusters:
        test_cluster = available_clusters[0]

    print(f"\nTesting annotation boost for: {test_cluster}")

    # Get marker data path
    marker_path = str(get_local_marker_file_path())

    # Run annotation boost
    start_time = time.time()
    errors = []
    status = "error"
    boost_results = {}

    output_name = str(results_dir['outputs'] / f"boost_{test_cluster.replace(' ', '_')}")

    try:
        print(f"\nRunning annotation boost (pip install mode)...")
        print(f"  Cluster: {test_cluster}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")
        print(f"  Max iterations: 3 (reduced for testing)")

        result = CASSIA.runCASSIA_annotationboost(
            full_result_path=batch_results_file,
            marker=marker_path,
            cluster_name=test_cluster,
            major_cluster_info=f"{data_config.get('species', 'human')} {data_config.get('tissue', 'large intestine')}",
            output_name=output_name,
            num_iterations=3,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            provider=llm_config.get('provider', 'openrouter'),
            temperature=llm_config.get('temperature', 0.3),
            conversation_history_mode="final",
            search_strategy="breadth",
            report_style="per_iteration"
        )

        # Check result
        if isinstance(result, dict):
            boost_results = result

            if result.get('status') == 'success':
                print(f"\nAnnotation Boost Results:")
                print(f"  Status: {result.get('status')}")
                print(f"  Execution time: {result.get('execution_time', 0):.1f}s")

                analysis_text = result.get('analysis_text', '')
                if analysis_text and len(analysis_text) > 100:
                    status = "passed"
                else:
                    status = "failed"
                    errors.append("Analysis text is empty or too short")
            else:
                status = "failed"
                errors.append(result.get('error_message', 'Unknown error'))
        else:
            status = "failed"
            errors.append("Unexpected result format from annotation boost")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="annotation_boost_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[test_cluster],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "cluster": test_cluster,
        "batch_results_file": batch_results_file,
        "boost_results": {
            "status": boost_results.get('status'),
            "execution_time": boost_results.get('execution_time'),
        },
        "mode": "pip_install"
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("06_annotation_boost")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_annotation_boost_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
