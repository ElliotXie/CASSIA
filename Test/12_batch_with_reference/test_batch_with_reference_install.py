"""
CASSIA Test 12: Batch Annotation with Reference (PIP INSTALL MODE)
===================================================================
Tests the runCASSIA_batch_with_reference function using pip-installed CASSIA.
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe, get_all_clusters
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


def run_batch_with_reference_test(results_dir):
    print_test_header("12 - Batch Annotation with Reference (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    config = load_config()
    print_config_summary(config)
    setup_api_keys()

    llm_config = config['llm']
    data_config = config['data']

    all_clusters = get_all_clusters()
    marker_df = get_full_marker_dataframe()

    # Results directory passed in from main()
    output_name = str(results_dir['outputs'] / "batch_ref_results")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch_with_reference (pip install mode)...")
        CASSIA.runCASSIA_batch_with_reference(
            marker=marker_df,
            output_name=output_name,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1'),
            use_reference=True,
            verbose=True
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")

        if full_csv.exists():
            import pandas as pd
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            has_ref_column = 'Reference Used' in results_df.columns

            print(f"\nBatch with Reference Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(all_clusters)}")
            print(f"  Reference column present: {has_ref_column}")

            if clusters_annotated == len(all_clusters) and has_ref_column:
                status = "passed"
            else:
                status = "failed"
                if not has_ref_column:
                    errors.append("Reference Used column not found")
        else:
            status = "failed"
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")

    duration = time.time() - start_time

    metadata = create_test_metadata(
        test_name="batch_with_reference_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters,
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    metadata['use_reference'] = True
    save_test_metadata(results_dir['base'], metadata)

    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("12_batch_with_reference")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_batch_with_reference_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
