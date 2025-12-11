"""
CASSIA Test 19: Reasoning Effort
================================
Tests the reasoning parameter functionality.

Tests:
1. Single annotation with LOW reasoning (OpenRouter)
2. Batch annotation with GPT-5.1 (OpenAI, no reasoning)
3. Batch annotation with LOW reasoning (OpenRouter)

Usage:
    python test_reasoning_effort.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import load_markers, get_cluster_markers
from test_utils import (
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def test_single_low_reasoning():
    """Test single annotation with LOW reasoning effort via OpenRouter."""
    print("\n" + "="*60)
    print("TEST 1: Single Annotation - LOW Reasoning")
    print("         Provider: OpenRouter")
    print("="*60)

    from CASSIA import runCASSIA

    # Load monocyte markers for quick test
    marker_list = get_cluster_markers('monocyte', n_genes=15)

    print(f"\nProvider: openrouter")
    print(f"Model: openai/gpt-4o")
    print(f"Reasoning: {{'effort': 'low'}}")
    print(f"Markers: {len(marker_list)} genes")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA...")
        result, conv, extra = runCASSIA(
            marker_list=marker_list,
            model="openai/gpt-4o",
            temperature=0.3,
            tissue="blood",
            species="human",
            provider="openrouter",
            reasoning={"effort": "low"},
            validator_involvement="v1"
        )

        if result and result.get("main_cell_type"):
            print(f"\nResults:")
            print(f"  Main cell type: {result.get('main_cell_type')}")
            print(f"  Iterations: {result.get('iterations', 'N/A')}")
            status = "passed"
        else:
            status = "failed"
            errors.append("No valid result returned")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time
    print(f"\nDuration: {duration:.2f}s")
    print(f"Status: {status.upper()}")

    return status == "passed", duration, errors


def test_batch_gpt5_no_reasoning(output_dir):
    """Test batch annotation with GPT-5.1 (no reasoning) via direct OpenAI."""
    print("\n" + "="*60)
    print("TEST 2: Batch Annotation - GPT-5.1 (NO reasoning)")
    print("         Provider: OpenAI (direct)")
    print("="*60)

    from CASSIA import runCASSIA_batch
    import pandas as pd

    # Load marker data for 2 clusters
    full_df = load_markers()
    test_clusters = ['monocyte', 'plasma cell']
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nProvider: openai")
    print(f"Model: gpt-5.1")
    print(f"Reasoning: None (Chat Completions API)")
    print(f"Clusters: {test_clusters}")

    output_name = str(output_dir / "batch_gpt5_no_reasoning")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=15,
            model="gpt-5.1",
            temperature=0.3,
            tissue="blood",
            species="human",
            max_workers=2,
            provider="openai",
            # NO reasoning parameter - uses Chat Completions API
            validator_involvement="v1"
        )

        # Check output files
        summary_csv = Path(f"{output_name}_summary.csv")

        if summary_csv.exists():
            results_df = pd.read_csv(summary_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(test_clusters)}")
            print(f"  Output file: {summary_csv.name}")

            if clusters_annotated == len(test_clusters):
                status = "passed"
            else:
                status = "failed"
                errors.append(f"Only {clusters_annotated}/{len(test_clusters)} clusters annotated")
        else:
            status = "failed"
            errors.append("Output CSV file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time
    print(f"\nDuration: {duration:.2f}s")
    print(f"Status: {status.upper()}")

    return status == "passed", duration, errors


def test_batch_low_reasoning(output_dir):
    """Test batch annotation with LOW reasoning effort via OpenRouter."""
    print("\n" + "="*60)
    print("TEST 3: Batch Annotation - LOW Reasoning")
    print("         Provider: OpenRouter")
    print("="*60)

    from CASSIA import runCASSIA_batch
    import pandas as pd

    # Load marker data for 2 clusters
    full_df = load_markers()
    test_clusters = ['monocyte', 'plasma cell']
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nProvider: openrouter")
    print(f"Model: openai/gpt-4o")
    print(f"Reasoning: {{'effort': 'low'}}")
    print(f"Clusters: {test_clusters}")

    output_name = str(output_dir / "batch_low_reasoning")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=15,
            model="openai/gpt-4o",
            temperature=0.3,
            tissue="blood",
            species="human",
            max_workers=2,
            provider="openrouter",
            reasoning={"effort": "low"},
            validator_involvement="v1"
        )

        # Check output files
        summary_csv = Path(f"{output_name}_summary.csv")

        if summary_csv.exists():
            results_df = pd.read_csv(summary_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(test_clusters)}")
            print(f"  Output file: {summary_csv.name}")

            if clusters_annotated == len(test_clusters):
                status = "passed"
            else:
                status = "failed"
                errors.append(f"Only {clusters_annotated}/{len(test_clusters)} clusters annotated")
        else:
            status = "failed"
            errors.append("Output CSV file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time
    print(f"\nDuration: {duration:.2f}s")
    print(f"Status: {status.upper()}")

    return status == "passed", duration, errors


def run_reasoning_effort_test():
    """Run reasoning effort tests."""
    print_test_header("19 - Reasoning Effort")

    # Load configuration
    config = load_config()

    # Setup API keys
    setup_api_keys()

    # Create results directory
    results = create_results_dir("19_reasoning_effort", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging
    logging_ctx = setup_logging(results['logs'])

    total_start = time.time()
    all_errors = []
    test_results = {}

    # Test 1: Single annotation with LOW reasoning
    print("\n" + "#"*70)
    print("# RUNNING TEST 1: Single Annotation (LOW reasoning)")
    print("#"*70)
    success1, duration1, errors1 = test_single_low_reasoning()
    test_results['single_low'] = {'passed': success1, 'duration': duration1, 'errors': errors1}
    all_errors.extend(errors1)

    # Test 2: Batch annotation with GPT-5.1 (no reasoning)
    print("\n" + "#"*70)
    print("# RUNNING TEST 2: Batch Annotation GPT-5.1 (NO reasoning)")
    print("#"*70)
    success2, duration2, errors2 = test_batch_gpt5_no_reasoning(results['outputs'])
    test_results['batch_gpt5'] = {'passed': success2, 'duration': duration2, 'errors': errors2}
    all_errors.extend(errors2)

    # Test 3: Batch annotation with LOW reasoning
    print("\n" + "#"*70)
    print("# RUNNING TEST 3: Batch Annotation (LOW reasoning)")
    print("#"*70)
    success3, duration3, errors3 = test_batch_low_reasoning(results['outputs'])
    test_results['batch_low'] = {'passed': success3, 'duration': duration3, 'errors': errors3}
    all_errors.extend(errors3)

    total_duration = time.time() - total_start

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTest 1 (Single, LOW reasoning):      {'PASSED' if success1 else 'FAILED'} ({duration1:.2f}s)")
    print(f"Test 2 (Batch GPT-5.1, no reasoning): {'PASSED' if success2 else 'FAILED'} ({duration2:.2f}s)")
    print(f"Test 3 (Batch, LOW reasoning):        {'PASSED' if success3 else 'FAILED'} ({duration3:.2f}s)")
    print(f"\nTotal Duration: {total_duration:.2f}s")

    # Overall status
    all_passed = success1 and success2 and success3
    status = "passed" if all_passed else "failed"

    # Save metadata
    metadata = create_test_metadata(
        test_name="reasoning_effort",
        config={
            "tests": {
                "single_low": {"model": "openai/gpt-4o", "reasoning": {"effort": "low"}, "provider": "openrouter"},
                "batch_gpt5": {"model": "gpt-5.1", "reasoning": None, "provider": "openai", "type": "batch"},
                "batch_low": {"model": "openai/gpt-4o", "reasoning": {"effort": "low"}, "provider": "openrouter", "type": "batch"}
            }
        },
        duration_seconds=total_duration,
        status=status,
        errors=all_errors
    )
    save_test_metadata(results['outputs'], metadata)

    # Print final result
    print_test_result(all_passed, f"Total Duration: {total_duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return all_passed


if __name__ == "__main__":
    success = run_reasoning_effort_test()
    sys.exit(0 if success else 1)
