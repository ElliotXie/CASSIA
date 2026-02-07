"""
CASSIA Test 09: Subclustering
==============================
Tests the subclustering annotation functions for annotating subclusters
from a major cluster, using results from RunCassia batch as context.

Usage:
    python test_subclustering.py

Prerequisites:
    - Run test 02_batch_annotation first to generate batch results

Functions tested:
- runCASSIA_subclusters(): Single-run subcluster annotation
- annotate_subclusters(): Core subcluster annotation function
"""

import sys
import time
import json
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


def find_latest_batch_results():
    """
    Find the most recent batch annotation results folder.

    Returns:
        dict with paths to summary_csv, conversations_json, and outputs_dir,
        or None if no valid results found.
    """
    batch_results_dir = Path(__file__).parent.parent / "02_batch_annotation" / "results" / "python" / "development"

    if not batch_results_dir.exists():
        print(f"  Batch results directory not found: {batch_results_dir}")
        return None

    # Find latest timestamped folder (sorted in reverse order)
    folders = sorted([f for f in batch_results_dir.iterdir() if f.is_dir()], reverse=True)

    for folder in folders:
        outputs = folder / "outputs"
        if outputs.exists():
            summary_csv = list(outputs.glob("*_summary.csv"))
            conversations_json = list(outputs.glob("*_conversations.json"))
            if summary_csv and conversations_json:
                return {
                    'summary_csv': summary_csv[0],
                    'conversations_json': conversations_json[0],
                    'outputs_dir': outputs,
                    'timestamp': folder.name
                }

    return None


def load_batch_results(batch_paths):
    """
    Load summary CSV and conversations JSON from batch results.

    Args:
        batch_paths: dict with paths to summary_csv and conversations_json

    Returns:
        tuple of (summary_df, conversations_dict)
    """
    summary_df = pd.read_csv(batch_paths['summary_csv'])

    with open(batch_paths['conversations_json'], 'r', encoding='utf-8') as f:
        conversations = json.load(f)

    return summary_df, conversations


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
    from CASSIA import runCASSIA_subclusters

    # Create results directory with organized structure
    results = create_results_dir("09_subclustering", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # Find and load batch annotation results
    print("\n--- Loading Batch Annotation Results ---")
    batch_paths = find_latest_batch_results()

    if not batch_paths:
        print("ERROR: No batch annotation results found.")
        print("Please run test 02_batch_annotation first to generate batch results.")
        cleanup_logging(logging_ctx)
        return False

    print(f"  Using batch results from: {batch_paths['timestamp']}")
    print(f"  Summary CSV: {batch_paths['summary_csv'].name}")
    print(f"  Conversations JSON: {batch_paths['conversations_json'].name}")

    # Load batch results
    summary_df, conversations = load_batch_results(batch_paths)

    # Validate batch results structure
    required_columns = ['Cluster ID', 'Predicted General Cell Type']
    missing_columns = [col for col in required_columns if col not in summary_df.columns]
    if missing_columns:
        print(f"ERROR: Batch summary CSV missing required columns: {missing_columns}")
        cleanup_logging(logging_ctx)
        return False

    # Pick first cluster from batch results for context
    cluster_row = summary_df.iloc[0]
    cluster_id = cluster_row['Cluster ID']
    predicted_cell_type = cluster_row['Predicted General Cell Type']

    print(f"\n  Selected cluster: '{cluster_id}' -> {predicted_cell_type}")
    print(f"  Conversation history available: {cluster_id in conversations}")

    # Create major_cluster_info using the annotated cell type from batch results
    major_cluster_info = f"{data_config.get('species', 'human')} {data_config.get('tissue', 'large intestine')} {predicted_cell_type}"

    # For subclustering test, we'll create a simulated subcluster marker set
    # In real usage, these would come from re-clustering a major cell type cluster
    # Load ALL clusters from the marker file (not just the first 2 for comprehensive testing)
    from fixtures import load_markers
    marker_df = load_markers()  # Load all clusters instead of just the first 2

    # Use ALL clusters as "subclusters" for comprehensive testing
    subcluster_df = marker_df.copy()

    # TEST: Recode Result ID column to numeric (0, 1, 2, ...) to verify fix works with numeric cluster IDs
    # Original cluster names are preserved in the index, but we're testing with numeric IDs
    subcluster_df.iloc[:, 0] = [str(i) for i in range(len(subcluster_df))]  # Replace cluster names with 0, 1, 2, ...

    subcluster_df = subcluster_df.reset_index(drop=True)

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    subcluster_results = {}
    batch_context = {
        'batch_timestamp': batch_paths['timestamp'],
        'source_cluster_id': cluster_id,
        'source_cell_type': predicted_cell_type,
        'has_conversation_history': cluster_id in conversations
    }

    try:
        print(f"\n--- Test: runCASSIA_subclusters ---")
        print(f"  Major cluster: {major_cluster_info}")
        print(f"  Subclusters to annotate: {len(subcluster_df)}")
        print(f"  Model: {llm_config.get('model', 'google/gemini-2.5-flash')}")
        print(f"  Provider: {llm_config.get('provider', 'openrouter')}")

        output_name = str(results['outputs'] / "subcluster_results")

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
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "major_cluster_info": major_cluster_info,
        "num_input_subclusters": len(subcluster_df),
        "batch_context": batch_context,
        "results": subcluster_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_subclustering_test()
    sys.exit(0 if success else 1)
