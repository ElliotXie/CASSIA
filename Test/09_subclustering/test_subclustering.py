"""
CASSIA Test 09: Subclustering
==============================
Tests the subclustering annotation functions for annotating subclusters
from a major cluster.

Usage:
    python test_subclustering.py

Functions tested:
- runCASSIA_subclusters(): Single-run subcluster annotation
- annotate_subclusters(): Core subcluster annotation function
"""

import sys
import time
from pathlib import Path
import pandas as pd

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe
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
    save_test_results,
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def run_subclustering_test():
    """Test subclustering annotation functionality."""
    print_test_header("09 - Subclustering")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from subclustering import runCASSIA_subclusters

    # Create results directory
    results_dir = create_results_dir("09_subclustering")
    print(f"Results will be saved to: {results_dir}")

    # For subclustering test, we'll create a simulated subcluster marker set
    # In real usage, these would come from re-clustering a major cell type cluster
    marker_df = get_full_marker_dataframe()

    # Create a subset as "subclusters" - use first 3 clusters as simulated subclusters
    subcluster_df = marker_df.head(3).copy()
    subcluster_df = subcluster_df.reset_index(drop=True)

    # Define major cluster context
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

        output_name = str(results_dir / "subcluster_results")

        runCASSIA_subclusters(
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

            if len(result_df) > 0:
                print(f"\n  Sample annotations:")
                for idx, row in result_df.head(3).iterrows():
                    print(f"    {row.get('Result ID', idx)}: {row.get('main_cell_type', 'N/A')} / {row.get('sub_cell_type', 'N/A')}")

            subcluster_results = {
                'output_file': str(csv_file),
                'num_subclusters': len(result_df),
                'columns': list(result_df.columns),
                'sample_results': result_df.head(3).to_dict('records') if len(result_df) > 0 else []
            }

            # Check for HTML report
            html_file = Path(f"{output_name}.html")
            if html_file.exists():
                print(f"  HTML report: {html_file.name}")
                subcluster_results['html_report'] = str(html_file)

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
        test_name="subclustering",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=[major_cluster_info],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)

    save_test_results(results_dir, {
        "major_cluster_info": major_cluster_info,
        "num_input_subclusters": len(subcluster_df),
        "results": subcluster_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_subclustering_test()
    sys.exit(0 if success else 1)
