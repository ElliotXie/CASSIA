"""
CASSIA Test 18: Reasoning Effort
================================
Tests the reasoning parameter with runCASSIA for:
1. Direct OpenAI - GPT-4o (no reasoning) - Chat Completions API
2. Direct OpenAI - GPT-5.1 (no reasoning) - Chat Completions API
3. Direct OpenAI - GPT-5 (with reasoning) - Responses API [SKIPPED - org verification]
4. Batch Annotation - GPT-5.1 (no reasoning) - CSV/HTML outputs
5. OpenRouter GPT-5.1 - reasoning effort LOW

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


def test_openai_gpt4o_no_reasoning():
    """Test direct OpenAI with GPT-4o (no reasoning) - uses Chat Completions API."""
    print("\n" + "="*60)
    print("TEST 1: Direct OpenAI GPT-4o (NO reasoning)")
    print("         API: Chat Completions")
    print("="*60)

    from CASSIA import runCASSIA

    # Load monocyte markers for quick test
    marker_list = get_cluster_markers('monocyte', n_genes=15)

    print(f"\nProvider: openai (direct)")
    print(f"Model: gpt-4o")
    print(f"Reasoning: None (Chat Completions API)")
    print(f"Markers: {len(marker_list)} genes")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA...")
        result, conv, extra = runCASSIA(
            marker_list=marker_list,
            model="gpt-4o",
            temperature=0.3,
            tissue="blood",
            species="human",
            provider="openai",  # Direct OpenAI
            # NO reasoning parameter - uses Chat Completions API
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


def test_openai_gpt5_no_reasoning():
    """Test direct OpenAI with GPT-5.1 (no reasoning) - uses Chat Completions API."""
    print("\n" + "="*60)
    print("TEST 2: Direct OpenAI GPT-5.1 (NO reasoning)")
    print("         API: Chat Completions")
    print("="*60)

    from CASSIA import runCASSIA

    # Load monocyte markers for quick test
    marker_list = get_cluster_markers('monocyte', n_genes=15)

    print(f"\nProvider: openai (direct)")
    print(f"Model: gpt-5.1")
    print(f"Reasoning: None (Chat Completions API)")
    print(f"Markers: {len(marker_list)} genes")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA...")
        result, conv, extra = runCASSIA(
            marker_list=marker_list,
            model="gpt-5.1",
            temperature=0.3,
            tissue="blood",
            species="human",
            provider="openai",  # Direct OpenAI
            # NO reasoning parameter - should use Chat Completions API
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


def test_openai_gpt5_with_reasoning():
    """Test direct OpenAI with GPT-5 and reasoning - uses Responses API."""
    print("\n" + "="*60)
    print("TEST 3: Direct OpenAI GPT-5 (WITH reasoning)")
    print("         API: Responses API")
    print("="*60)

    from CASSIA import runCASSIA

    # Load monocyte markers for quick test
    marker_list = get_cluster_markers('monocyte', n_genes=15)

    print(f"\nProvider: openai (direct)")
    print(f"Model: gpt-5")
    print(f"Reasoning: {{'effort': 'medium'}} (Responses API)")
    print(f"Markers: {len(marker_list)} genes")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA...")
        result, conv, extra = runCASSIA(
            marker_list=marker_list,
            model="gpt-5",
            temperature=0.3,
            tissue="blood",
            species="human",
            provider="openai",  # Direct OpenAI
            reasoning={"effort": "medium"},  # Uses Responses API
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
    """Test batch annotation with GPT-5.1 (no reasoning) - generates CSV/HTML reports."""
    print("\n" + "="*60)
    print("TEST 4: Batch Annotation GPT-5.1 (NO reasoning)")
    print("         API: Chat Completions")
    print("         Output: CSV + HTML reports")
    print("="*60)

    from CASSIA import runCASSIA_batch
    import pandas as pd

    # Load marker data for 2 clusters
    full_df = load_markers()
    test_clusters = ['monocyte', 'plasma cell']
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nProvider: openai (direct)")
    print(f"Model: gpt-5.1")
    print(f"Reasoning: None (Chat Completions API)")
    print(f"Clusters: {test_clusters}")

    output_name = str(output_dir / "batch_gpt5_results")

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
            provider="openai",  # Direct OpenAI
            # NO reasoning parameter - uses Chat Completions API
            validator_involvement="v1"
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")
        summary_csv = Path(f"{output_name}_summary.csv")

        if full_csv.exists():
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(test_clusters)}")
            print(f"  Output files created:")
            print(f"    - {full_csv.name}")
            if summary_csv.exists():
                print(f"    - {summary_csv.name}")

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


def test_openrouter_reasoning_effort(output_dir, effort_level: str):
    """Test OpenRouter GPT-5.1 with specified reasoning effort level."""
    print("\n" + "="*60)
    print(f"TEST: OpenRouter GPT-5.1 - Reasoning Effort: {effort_level.upper()}")
    print("="*60)

    from CASSIA import runCASSIA_batch
    import pandas as pd

    # Load marker data for 2 clusters
    full_df = load_markers()
    test_clusters = ['monocyte', 'plasma cell']
    marker_df = full_df[full_df['Broad.cell.type'].isin(test_clusters)].copy()

    print(f"\nProvider: openrouter")
    print(f"Model: openai/gpt-5.1")
    print(f"Reasoning: {{'effort': '{effort_level}'}}")
    print(f"Clusters: {test_clusters}")

    output_name = str(output_dir / f"batch_openrouter_{effort_level}")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=15,
            model="openai/gpt-5.1",
            temperature=0.3,
            tissue="blood",
            species="human",
            max_workers=2,
            provider="openrouter",
            reasoning={"effort": effort_level},
            validator_involvement="v1"
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")
        summary_csv = Path(f"{output_name}_summary.csv")

        if full_csv.exists():
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results (effort={effort_level}):")
            print(f"  Clusters annotated: {clusters_annotated}/{len(test_clusters)}")
            print(f"  Output files:")
            print(f"    - {full_csv.name}")
            if summary_csv.exists():
                print(f"    - {summary_csv.name}")

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
    """Run all reasoning effort tests."""
    print_test_header("18 - Reasoning Effort (Direct OpenAI)")

    # Load configuration
    config = load_config()

    # Setup API keys
    setup_api_keys()

    # Create results directory
    results = create_results_dir("18_reasoning_effort", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging
    logging_ctx = setup_logging(results['logs'])

    total_start = time.time()
    all_errors = []
    test_results = {}

    # Test 1: GPT-4o without reasoning (Chat Completions API)
    print("\n" + "#"*70)
    print("# RUNNING TEST 1: Direct OpenAI GPT-4o (Chat Completions API)")
    print("#"*70)
    success1, duration1, errors1 = test_openai_gpt4o_no_reasoning()
    test_results['gpt4o'] = {'passed': success1, 'duration': duration1, 'errors': errors1}
    all_errors.extend(errors1)

    # Test 2: GPT-5.1 without reasoning (Chat Completions API)
    print("\n" + "#"*70)
    print("# RUNNING TEST 2: Direct OpenAI GPT-5.1 (Chat Completions API)")
    print("#"*70)
    success2, duration2, errors2 = test_openai_gpt5_no_reasoning()
    test_results['gpt5_no_reasoning'] = {'passed': success2, 'duration': duration2, 'errors': errors2}
    all_errors.extend(errors2)

    # Test 3: GPT-5 with reasoning (Responses API) - skip for now due to org verification
    print("\n" + "#"*70)
    print("# SKIPPING TEST 3: Direct OpenAI GPT-5 (Responses API)")
    print("# Reason: Requires organization verification")
    print("#"*70)
    success3 = True  # Skip this test
    duration3 = 0

    # Test 4: Batch annotation with GPT-5.1 (generates CSV/HTML)
    print("\n" + "#"*70)
    print("# RUNNING TEST 4: Batch Annotation GPT-5.1 (Chat Completions API)")
    print("#"*70)
    success4, duration4, errors4 = test_batch_gpt5_no_reasoning(results['outputs'])
    test_results['batch_gpt5'] = {'passed': success4, 'duration': duration4, 'errors': errors4}
    all_errors.extend(errors4)

    # Test 5: OpenRouter GPT-5.1 with reasoning effort LOW
    print("\n" + "#"*70)
    print("# RUNNING TEST 5: OpenRouter GPT-5.1 (Reasoning: LOW)")
    print("#"*70)
    success5, duration5, errors5 = test_openrouter_reasoning_effort(results['outputs'], "low")
    test_results['openrouter_low'] = {'passed': success5, 'duration': duration5, 'errors': errors5}
    all_errors.extend(errors5)

    total_duration = time.time() - total_start

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nTest 1 (GPT-4o, Chat Completions):      {'PASSED' if success1 else 'FAILED'} ({duration1:.2f}s)")
    print(f"Test 2 (GPT-5.1, Chat Completions):     {'PASSED' if success2 else 'FAILED'} ({duration2:.2f}s)")
    print(f"Test 3 (GPT-5, Responses API):          SKIPPED (org verification)")
    print(f"Test 4 (Batch GPT-5.1, no reasoning):   {'PASSED' if success4 else 'FAILED'} ({duration4:.2f}s)")
    print(f"Test 5 (OpenRouter, effort=LOW):        {'PASSED' if success5 else 'FAILED'} ({duration5:.2f}s)")
    print(f"\nTotal Duration: {total_duration:.2f}s")

    # Overall status - based on all active tests
    all_passed = success1 and success2 and success4 and success5
    status = "passed" if all_passed else "failed"

    # Save metadata
    metadata = create_test_metadata(
        test_name="reasoning_effort_comparison",
        config={
            "tests": {
                "gpt4o": {"model": "gpt-4o", "reasoning": None, "api": "Chat Completions"},
                "gpt5_no_reasoning": {"model": "gpt-5.1", "reasoning": None, "api": "Chat Completions"},
                "gpt5_with_reasoning": {"model": "gpt-5", "reasoning": {"effort": "medium"}, "api": "Responses", "status": "skipped"},
                "batch_gpt5": {"model": "gpt-5.1", "reasoning": None, "api": "Chat Completions", "type": "batch"},
                "openrouter_low": {"model": "openai/gpt-5.1", "reasoning": {"effort": "low"}, "provider": "openrouter"}
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
