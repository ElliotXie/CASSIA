"""
CASSIA Test 24: Bugfix Validation
==================================
Tests for fixes reported by user tjli:

1. Fix #1: OpenRouter scoring should not send `reasoning` param for Anthropic models (400 error)
2. Fix #3: Subcluster report key_markers should not be null
3. Fix #4: Subcluster functions should support additional_context parameter

Usage:
    python test_bugfix_validation.py
"""

import sys
import time
from pathlib import Path
import pandas as pd

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import load_markers
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


def test_fix1_reasoning_filter():
    """
    Fix #1: Verify reasoning parameter is only sent to models that support it.
    Tests that OpenRouter calls with Anthropic Claude models do not include reasoning.
    """
    print("\n" + "=" * 60)
    print("TEST 1: Reasoning param filtered for non-supporting models")
    print("=" * 60)

    from CASSIA.core.llm_utils import call_llm

    errors = []
    status = "error"

    try:
        # Test: call_llm with an Anthropic model and reasoning param via OpenRouter
        # This should NOT raise a 400 error because the fix filters out reasoning
        # for non-supporting models. We do a simple annotation call.
        print("\n  Testing: Anthropic model via OpenRouter with reasoning param...")
        print("  Model: anthropic/claude-sonnet-4-5")
        print("  Reasoning: {'effort': 'low'}")

        result = call_llm(
            prompt="What cell type expresses CD3D, CD3E, CD4, IL7R? Reply in one sentence.",
            provider="openrouter",
            model="anthropic/claude-sonnet-4-5",
            max_tokens=200,
            temperature=0,
            reasoning={"effort": "low"}  # This should be silently ignored for Claude
        )

        if result and len(result) > 0:
            print(f"  Result: {result[:100]}...")
            print("  PASS: No 400 error - reasoning param correctly filtered")
            status = "passed"
        else:
            errors.append("Empty result returned")
            status = "failed"

    except Exception as e:
        error_msg = str(e)
        if "400" in error_msg:
            errors.append(f"400 error still occurring - fix not working: {error_msg}")
            status = "failed"
        else:
            errors.append(error_msg)
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def test_fix3_subcluster_markers():
    """
    Fix #3: Verify subcluster results CSV has non-empty key_markers column.
    """
    print("\n" + "=" * 60)
    print("TEST 2: Subcluster key_markers populated (not null)")
    print("=" * 60)

    from CASSIA import runCASSIA_subclusters

    config = load_config()
    llm_config = config['llm']
    data_config = config['data']

    errors = []
    status = "error"

    try:
        # Load markers and use first 2 clusters as subclusters
        marker_df = load_markers()
        subcluster_df = marker_df.head(2).copy()
        subcluster_df.iloc[:, 0] = [str(i) for i in range(len(subcluster_df))]
        subcluster_df = subcluster_df.reset_index(drop=True)

        # Store expected markers for validation
        expected_markers = {str(i): str(subcluster_df.iloc[i, 1]) for i in range(len(subcluster_df))}

        major_cluster_info = "human large intestine immune cells"

        results_dir = create_results_dir("24_bugfix_validation", get_test_mode())
        output_name = str(results_dir['outputs'] / "subcluster_markers_test")

        print(f"\n  Subclusters: {len(subcluster_df)}")
        print(f"  Model: {llm_config.get('model')}")
        print(f"  Provider: {llm_config.get('provider')}")

        runCASSIA_subclusters(
            marker=subcluster_df,
            major_cluster_info=major_cluster_info,
            output_name=output_name,
            model=llm_config.get('model'),
            temperature=llm_config.get('temperature', 0.3),
            provider=llm_config.get('provider', 'openrouter'),
            n_genes=data_config.get('n_genes', 30)
        )

        csv_file = Path(f"{output_name}.csv")
        if csv_file.exists():
            result_df = pd.read_csv(csv_file)
            print(f"\n  Output: {csv_file.name}")
            print(f"  Rows: {len(result_df)}")

            # Check key_markers column
            if 'key_markers' not in result_df.columns:
                errors.append("key_markers column missing from output")
                status = "failed"
            else:
                null_markers = result_df['key_markers'].isna() | (result_df['key_markers'] == '') | (result_df['key_markers'] == 'null')
                null_count = null_markers.sum()

                for idx, row in result_df.iterrows():
                    markers_val = row.get('key_markers', '')
                    rid = str(row.get('Result ID', idx))
                    preview = str(markers_val)[:80] if pd.notna(markers_val) else 'NULL'
                    print(f"    Cluster {rid} markers: {preview}...")

                if null_count == 0:
                    print(f"\n  PASS: All {len(result_df)} clusters have non-null key_markers")
                    status = "passed"
                else:
                    errors.append(f"{null_count}/{len(result_df)} clusters have null/empty key_markers")
                    status = "failed"
        else:
            errors.append(f"Output file not created: {csv_file}")
            status = "failed"

    except Exception as e:
        errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def test_fix4_additional_context():
    """
    Fix #4: Verify subcluster functions accept and use additional_context parameter.
    """
    print("\n" + "=" * 60)
    print("TEST 3: Subcluster additional_context parameter support")
    print("=" * 60)

    from CASSIA import runCASSIA_subclusters

    config = load_config()
    llm_config = config['llm']
    data_config = config['data']

    errors = []
    status = "error"

    try:
        # Load markers and use first 2 clusters as subclusters
        marker_df = load_markers()
        subcluster_df = marker_df.head(2).copy()
        subcluster_df.iloc[:, 0] = [str(i) for i in range(len(subcluster_df))]
        subcluster_df = subcluster_df.reset_index(drop=True)

        major_cluster_info = "human large intestine T cells"
        additional_context = "These subclusters are from colonic lamina propria T cells in an IBD patient sample. The tissue shows signs of active inflammation."

        results_dir = create_results_dir("24_bugfix_validation", get_test_mode())
        output_name = str(results_dir['outputs'] / "subcluster_context_test")

        print(f"\n  Subclusters: {len(subcluster_df)}")
        print(f"  additional_context: {additional_context[:60]}...")
        print(f"  Model: {llm_config.get('model')}")

        runCASSIA_subclusters(
            marker=subcluster_df,
            major_cluster_info=major_cluster_info,
            output_name=output_name,
            model=llm_config.get('model'),
            temperature=llm_config.get('temperature', 0.3),
            provider=llm_config.get('provider', 'openrouter'),
            n_genes=data_config.get('n_genes', 30),
            additional_context=additional_context  # NEW PARAMETER
        )

        csv_file = Path(f"{output_name}.csv")
        if csv_file.exists():
            result_df = pd.read_csv(csv_file)
            print(f"\n  Output: {csv_file.name}")
            print(f"  Rows: {len(result_df)}")

            for idx, row in result_df.iterrows():
                rid = row.get('Result ID', idx)
                ct1 = row.get('main_cell_type', 'N/A')
                ct2 = row.get('sub_cell_type', 'N/A')
                print(f"    Cluster {rid}: {ct1} / {ct2}")

            if len(result_df) > 0:
                print(f"\n  PASS: additional_context accepted, {len(result_df)} clusters annotated")
                status = "passed"
            else:
                errors.append("No results returned")
                status = "failed"
        else:
            errors.append(f"Output file not created: {csv_file}")
            status = "failed"

    except TypeError as e:
        if "additional_context" in str(e):
            errors.append(f"additional_context parameter not accepted: {e}")
            status = "failed"
        else:
            errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()
    except Exception as e:
        errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def run_bugfix_validation_test():
    """Run all bugfix validation tests."""
    print_test_header("24 - Bugfix Validation (tjli report)")

    # Load configuration and setup
    config = load_config()
    setup_api_keys()

    total_start = time.time()
    all_errors = []
    test_results = {}

    # Test 1: Reasoning filter
    print("\n" + "#" * 70)
    print("# FIX #1: OpenRouter reasoning param filter")
    print("#" * 70)
    success1, errors1 = test_fix1_reasoning_filter()
    test_results['fix1_reasoning'] = {'passed': success1, 'errors': errors1}
    all_errors.extend(errors1)

    # Test 2: Subcluster markers
    print("\n" + "#" * 70)
    print("# FIX #3: Subcluster key_markers populated")
    print("#" * 70)
    success2, errors2 = test_fix3_subcluster_markers()
    test_results['fix3_markers'] = {'passed': success2, 'errors': errors2}
    all_errors.extend(errors2)

    # Test 3: Additional context
    print("\n" + "#" * 70)
    print("# FIX #4: Subcluster additional_context support")
    print("#" * 70)
    success3, errors3 = test_fix4_additional_context()
    test_results['fix4_context'] = {'passed': success3, 'errors': errors3}
    all_errors.extend(errors3)

    total_duration = time.time() - total_start

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  Fix #1 (Reasoning filter):      {'PASSED' if success1 else 'FAILED'}")
    print(f"  Fix #3 (Subcluster markers):     {'PASSED' if success2 else 'FAILED'}")
    print(f"  Fix #4 (additional_context):     {'PASSED' if success3 else 'FAILED'}")
    print(f"\n  Total Duration: {total_duration:.2f}s")

    if all_errors:
        print(f"\n  Errors:")
        for err in all_errors:
            print(f"    - {err[:100]}")

    all_passed = success1 and success2 and success3
    print_test_result(all_passed, f"Total Duration: {total_duration:.2f}s")

    return all_passed


if __name__ == "__main__":
    success = run_bugfix_validation_test()
    sys.exit(0 if success else 1)
