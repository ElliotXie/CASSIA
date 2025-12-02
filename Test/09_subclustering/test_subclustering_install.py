"""
CASSIA Test 09: Subclustering (PIP INSTALL MODE)
=================================================
Tests the subclustering annotation functions using pip-installed CASSIA.

Usage:
    python test_subclustering_install.py
"""

import sys
import time
from pathlib import Path
import pandas as pd

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe
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


def run_subclustering_test(results_dir):
    """Test subclustering annotation functionality using pip-installed CASSIA."""
    print_test_header("09 - Subclustering (PIP INSTALL MODE)")

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

    # Create simulated subcluster data
    marker_df = get_full_marker_dataframe()
    subcluster_df = marker_df.head(3).copy()
    subcluster_df = subcluster_df.reset_index(drop=True)

    major_cluster_info = f"{data_config.get('species', 'human')} {data_config.get('tissue', 'large intestine')} immune cells"

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    subcluster_results = {}

    try:
        print(f"\n--- Test: runCASSIA_subclusters ---")
        print(f"  Major cluster: {major_cluster_info}")
        print(f"  Subclusters to annotate: {len(subcluster_df)}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")

        output_name = str(results_dir['outputs'] / "subcluster_results")

        CASSIA.runCASSIA_subclusters(
            marker=subcluster_df,
            major_cluster_info=major_cluster_info,
            output_name=output_name,
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            provider=llm_config.get('provider', 'openrouter'),
            n_genes=data_config.get('n_genes', 30)
        )

        # Check output files
        csv_file = Path(f"{output_name}.csv")
        if csv_file.exists():
            result_df = pd.read_csv(csv_file)
            print(f"\nSubclustering Results:")
            print(f"  Output file: {csv_file.name}")
            print(f"  Subclusters annotated: {len(result_df)}")

            subcluster_results = {
                'output_file': str(csv_file),
                'num_subclusters': len(result_df),
                'columns': list(result_df.columns),
            }

            status = "passed"
        else:
            status = "failed"
            errors.append(f"Output file not created: {csv_file}")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="subclustering_install_py",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[major_cluster_info],
        errors=errors
    )
    metadata['pip_install_info'] = pip_info
    save_test_metadata(results_dir['base'], metadata)

    save_test_results(results_dir['base'], {
        "major_cluster_info": major_cluster_info,
        "num_input_subclusters": len(subcluster_df),
        "results": subcluster_results,
        "mode": "pip_install"
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


def main():
    """Main entry point with logging."""
    results_dir = create_results_dir("09_subclustering")
    logging_context = setup_logging(results_dir['logs'])

    try:
        success = run_subclustering_test(results_dir)
    finally:
        cleanup_logging(logging_context)

    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
