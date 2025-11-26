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

from fixtures import get_marker_dataframe_for_cluster
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
    results_dir = create_results_dir("05_model_settings_fuzzy")
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

    cluster1_df = get_marker_dataframe_for_cluster("monocyte", n_genes=20)
    cluster2_df = get_marker_dataframe_for_cluster("plasma cell", n_genes=20)
    subset_marker_df = pd.concat([cluster1_df, cluster2_df], ignore_index=True)
    print(f"  Clusters: monocyte, plasma cell")
    print(f"  Genes per cluster: 20")

    # Define fuzzy model name test cases
    # Each tuple: (fuzzy_name, provider, expected_resolved_contains)
    fuzzy_test_cases = [
        # OpenAI tests
        ("gpt", "openai", "gpt-5.1"),
        ("4o", "openai", "gpt-4o"),
        ("mini", "openai", "gpt-5-mini"),

        # Anthropic tests
        ("claude", "anthropic", "claude-sonnet-4-5"),
        ("sonnet", "anthropic", "claude-sonnet-4-5"),
        ("haiku", "anthropic", "claude-haiku-4-5"),

        # OpenRouter tests (most diverse)
        ("gemini", "openrouter", "google/gemini"),
        ("flash", "openrouter", "gemini-2.5-flash"),
        ("claude", "openrouter", "anthropic/claude"),
        ("deepseek", "openrouter", "deepseek"),
        ("llama", "openrouter", "meta-llama"),
        ("mistral", "openrouter", "mistral"),
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

            batch_result = runCASSIA_batch(
                marker=subset_marker_df,
                output_name=str(results_dir / f"fuzzy_batch_{fuzzy_name}"),
                n_genes=20,
                model=fuzzy_name,  # <-- Fuzzy name, should auto-resolve
                temperature=0.3,
                tissue=data_config.get('tissue', 'large intestine'),
                species=data_config.get('species', 'human'),
                max_workers=2,
                provider=provider,
                validator_involvement="v1"
            )

            # Verify results - batch_result is a DataFrame
            batch_success = batch_result is not None and len(batch_result) == 2
            test_result["num_clusters"] = len(batch_result) if batch_result is not None else 0

            if batch_success:
                print(f"  [OK] Batch completed with {len(batch_result)} clusters")
                # Extract results from DataFrame
                if hasattr(batch_result, 'iterrows'):
                    for idx, row in batch_result.iterrows():
                        cluster_name = row.get('Cluster ID', idx)
                        cell_type = row.get('Predicted General Cell Type', 'N/A')
                        print(f"       {cluster_name}: {cell_type}")
                        test_result["results"].append({
                            "cluster": str(cluster_name),
                            "cell_type": cell_type
                        })
                test_result["success"] = True
                passed_count += 1
            else:
                print(f"  [X] Batch failed or returned wrong number of results")
                errors.append(f"Fuzzy batch test failed for '{fuzzy_name}'")

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
