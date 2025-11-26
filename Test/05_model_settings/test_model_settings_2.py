"""
CASSIA Test 05b: Fuzzy Model Names in Batch Processing
=======================================================
Tests that runCASSIA_batch correctly resolves fuzzy/informal model names.

This test uses multiple fuzzy model name examples to verify the auto-resolution
feature works correctly in practical batch annotation scenarios.

Usage:
    python test_model_settings_2.py
"""

import sys
import time
import pandas as pd
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
    save_test_results,
    create_test_metadata
)

# Setup CASSIA imports
setup_cassia_imports()


def run_fuzzy_batch_test():
    """Test runCASSIA_batch with various fuzzy model names."""
    print_test_header("05b - Fuzzy Model Names in Batch")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Import CASSIA functions
    from tools_function import runCASSIA_batch

    # Create results directory
    results_dir = create_results_dir("05_model_settings")
    print(f"Results will be saved to: {results_dir}")

    start_time = time.time()
    errors = []
    test_results = {
        "fuzzy_batch_tests": []
    }

    data_config = config['data']

    # Prepare marker data - subset to 2 clusters for efficiency
    print("\n" + "="*50)
    print("Preparing marker data (2 clusters for efficiency)")
    print("="*50)

    # Load full marker data and filter to 2 clusters
    # Format: ['Broad.cell.type', 'Top.Markers'] - one row per cluster with comma-separated genes
    full_df = load_markers()
    test_clusters = ['monocyte', 'plasma cell']
    subset_marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()
    print(f"  Clusters: {', '.join(test_clusters)}")
    print(f"  Rows in marker dataframe: {len(subset_marker_df)}")

    # Define fuzzy model name test cases - one per provider
    # Each tuple: (fuzzy_name, provider, expected_resolved_contains)
    fuzzy_test_cases = [
        ("gpt", "openai", "gpt-5.1"),
        ("claude", "anthropic", "claude-sonnet-4-5"),
        ("gemini", "openrouter", "google/gemini-2.5-flash"),
    ]

    passed_count = 0
    total_tests = len(fuzzy_test_cases)

    for fuzzy_name, provider, expected_contains in fuzzy_test_cases:
        print("\n" + "="*50)
        print(f"Testing fuzzy model: '{fuzzy_name}' with provider: '{provider}'")
        print("="*50)

        test_result = {
            "fuzzy_name": fuzzy_name,
            "provider": provider,
            "expected_contains": expected_contains,
            "success": False,
            "num_clusters": 0,
            "results": []
        }

        try:
            print(f"  Running batch annotation...")
            print(f"  (Should print: Note: Resolved '{fuzzy_name}' to '...' for {provider})")

            # runCASSIA_batch doesn't return anything - it writes CSV files
            output_base = str(results_dir / f"fuzzy_batch_{fuzzy_name}_{provider}")
            runCASSIA_batch(
                marker=subset_marker_df,
                output_name=output_base,
                n_genes=20,
                model=fuzzy_name,  # <-- Fuzzy name, should auto-resolve
                temperature=0.3,
                tissue=data_config.get('tissue', 'large intestine'),
                species=data_config.get('species', 'human'),
                max_workers=2,
                provider=provider,
                validator_involvement="v1"
            )

            # Read the generated summary CSV to verify results
            summary_csv = Path(f"{output_base}_summary.csv")
            batch_success = False
            actual_model = None

            if summary_csv.exists():
                result_df = pd.read_csv(summary_csv)
                test_result["num_clusters"] = len(result_df)

                if len(result_df) == 2:
                    print(f"  [OK] Batch completed with {len(result_df)} clusters")

                    # Verify the resolved model name in output
                    if 'Model' in result_df.columns:
                        actual_model = result_df['Model'].iloc[0]
                        print(f"  Model used: {actual_model}")

                        # Check if resolved model contains expected string
                        if expected_contains in actual_model:
                            print(f"  [OK] Model correctly resolved: '{fuzzy_name}' -> '{actual_model}'")
                            batch_success = True
                        else:
                            print(f"  [X] Model NOT resolved! Expected '{expected_contains}' but got '{actual_model}'")
                            errors.append(f"Model resolution failed for '{fuzzy_name}': expected '{expected_contains}', got '{actual_model}'")

                    # Extract results from DataFrame
                    for idx, row in result_df.iterrows():
                        cluster_name = row.get('Cluster ID', idx)
                        cell_type = row.get('Predicted General Cell Type', 'N/A')
                        print(f"       {cluster_name}: {cell_type}")
                        test_result["results"].append({
                            "cluster": str(cluster_name),
                            "cell_type": cell_type
                        })
                else:
                    print(f"  [X] Batch completed but with {len(result_df)} clusters instead of 2")
                    errors.append(f"Fuzzy batch test failed for '{fuzzy_name}': got {len(result_df)} clusters")
            else:
                print(f"  [X] Summary CSV not found: {summary_csv}")
                errors.append(f"Fuzzy batch test failed for '{fuzzy_name}': CSV not generated")

            test_result["actual_model"] = actual_model
            test_result["success"] = batch_success
            if batch_success:
                passed_count += 1

        except Exception as e:
            test_result["error"] = str(e)
            errors.append(f"Fuzzy batch test error for '{fuzzy_name}': {e}")
            print(f"  [X] ERROR: {e}")

        test_results["fuzzy_batch_tests"].append(test_result)

    duration = time.time() - start_time

    # Summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)
    print(f"  Passed: {passed_count}/{total_tests}")
    print(f"  Duration: {duration:.2f}s")

    # Determine overall status
    status = "passed" if passed_count == total_tests else "failed"

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="model_settings_fuzzy_batch",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=["monocyte", "plasma cell"],
        errors=errors
    )
    save_test_metadata(results_dir, metadata)
    save_test_results(results_dir, test_results)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success


if __name__ == "__main__":
    success = run_fuzzy_batch_test()
    sys.exit(0 if success else 1)
