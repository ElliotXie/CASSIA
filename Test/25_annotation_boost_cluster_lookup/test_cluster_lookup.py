"""
CASSIA Test 25: Annotation Boost Cluster Lookup (Issue #18)
============================================================
Tests for the fix of GitHub Issue #18:
  "Annotation Boost fails with 'Cluster not found in conversations JSON'"

Two bugs fixed:
1. JSON key type mismatch: cluster_name is int (from pandas) but JSON keys are strings
2. NoneType error: call_llm returns None → 'completion_marker in conversation' crashes

These tests are unit tests that do NOT require API calls.

Usage:
    python test_cluster_lookup.py
"""

import sys
import os
import json
import time
import tempfile
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    setup_cassia_imports,
    print_test_header,
    print_test_result,
)

# Setup CASSIA imports
setup_cassia_imports()


def test_json_key_lookup_string_cluster():
    """
    Test that prepare_analysis_data finds cluster when cluster_name is a string
    and JSON key is also a string. (Baseline - should always work)
    """
    print("\n" + "=" * 60)
    print("TEST 1: String cluster_name matches string JSON key")
    print("=" * 60)

    from CASSIA.agents.annotation_boost.annotation_boost import prepare_analysis_data

    errors = []
    status = "error"

    try:
        # Create a temporary conversations JSON with string keys
        conversations = {
            "0": {
                "annotations": ["Cluster 0 is T cells based on CD3D, CD3E expression."],
                "validations": [],
                "formatting": "",
                "scoring": ""
            },
            "1": {
                "annotations": ["Cluster 1 is B cells based on CD19, MS4A1 expression."],
                "validations": [],
                "formatting": "",
                "scoring": ""
            }
        }

        # Create temp summary CSV
        import pandas as pd
        summary_data = pd.DataFrame({
            "Cluster ID": ["0", "1"],
            "Predicted General Cell Type": ["T cells", "B cells"],
            "Predicted Detailed Cell Type": ["CD4+ T cells", "Naive B cells"],
            "Possible Mixed Cell Types": ["", ""],
            "Marker Number": [10, 10],
            "Marker List": ["CD3D,CD3E,CD4", "CD19,MS4A1,CD79A"],
            "Iterations": [1, 1],
            "Model": ["test", "test"],
            "Provider": ["test", "test"],
            "Tissue": ["blood", "blood"],
            "Species": ["human", "human"]
        })

        # Create temp marker CSV
        marker_data = pd.DataFrame({
            "cluster": ["0", "0", "0", "1", "1", "1"],
            "gene": ["CD3D", "CD3E", "CD4", "CD19", "MS4A1", "CD79A"],
            "p_val_adj": [0.001] * 6,
            "avg_log2FC": [2.0] * 6
        })

        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = os.path.join(tmpdir, "test_summary.csv")
            marker_path = os.path.join(tmpdir, "test_markers.csv")
            json_path = os.path.join(tmpdir, "test_conversations.json")

            summary_data.to_csv(summary_path, index=False)
            marker_data.to_csv(marker_path, index=False)
            with open(json_path, 'w') as f:
                json.dump(conversations, f)

            # When CSV has string "0", pd.read_csv may read it as int 0.
            # Use the actual type that prepare_analysis_data will see after CSV read.
            # Read back to check the actual type:
            check_df = pd.read_csv(summary_path)
            actual_cluster = check_df["Cluster ID"].iloc[0]
            print(f"  Cluster ID type after CSV round-trip: {type(actual_cluster).__name__} = {actual_cluster!r}")

            _, _, _, annotation_history = prepare_analysis_data(
                full_result_path=summary_path,
                marker_path=marker_path,
                cluster_name=actual_cluster,
                conversation_history_mode="full",
                conversations_json_path=json_path,
                species="human",
                auto_convert_ids=False
            )

            if annotation_history and "T cells" in annotation_history:
                print("  PASS: String key found annotation history correctly")
                status = "passed"
            else:
                errors.append(f"Expected annotation history with 'T cells', got: {annotation_history!r}")
                status = "failed"

    except Exception as e:
        errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def test_json_key_lookup_int_cluster():
    """
    Test that prepare_analysis_data finds cluster when cluster_name is an int
    but JSON key is a string. This is the core bug from Issue #18.
    """
    print("\n" + "=" * 60)
    print("TEST 2: Integer cluster_name matches string JSON key (Issue #18)")
    print("=" * 60)

    from CASSIA.agents.annotation_boost.annotation_boost import prepare_analysis_data

    errors = []
    status = "error"

    try:
        # JSON always has string keys
        conversations = {
            "0": {
                "annotations": ["Cluster 0 is T cells based on CD3D, CD3E expression."],
                "validations": [],
                "formatting": "",
                "scoring": ""
            }
        }

        import pandas as pd
        # Use integer cluster IDs (as pandas would produce from numeric CSV data)
        summary_data = pd.DataFrame({
            "Cluster ID": [0, 1],  # integers, matching how pandas reads numeric CSV
            "Predicted General Cell Type": ["T cells", "B cells"],
            "Predicted Detailed Cell Type": ["CD4+ T cells", "Naive B cells"],
            "Possible Mixed Cell Types": ["", ""],
            "Marker Number": [10, 10],
            "Marker List": ["CD3D,CD3E,CD4", "CD19,MS4A1,CD79A"],
            "Iterations": [1, 1],
            "Model": ["test", "test"],
            "Provider": ["test", "test"],
            "Tissue": ["blood", "blood"],
            "Species": ["human", "human"]
        })

        marker_data = pd.DataFrame({
            "cluster": [0, 0, 0],
            "gene": ["CD3D", "CD3E", "CD4"],
            "p_val_adj": [0.001] * 3,
            "avg_log2FC": [2.0] * 3
        })

        with tempfile.TemporaryDirectory() as tmpdir:
            summary_path = os.path.join(tmpdir, "test_summary.csv")
            marker_path = os.path.join(tmpdir, "test_markers.csv")
            json_path = os.path.join(tmpdir, "test_conversations.json")

            summary_data.to_csv(summary_path, index=False)
            marker_data.to_csv(marker_path, index=False)
            with open(json_path, 'w') as f:
                json.dump(conversations, f)

            # Call with INTEGER cluster_name 0 — before the fix, this would fail
            # with "Cluster 0 not found in conversations JSON"
            _, _, _, annotation_history = prepare_analysis_data(
                full_result_path=summary_path,
                marker_path=marker_path,
                cluster_name=0,  # INTEGER — this is the bug trigger
                conversation_history_mode="full",
                conversations_json_path=json_path,
                species="human",
                auto_convert_ids=False
            )

            if annotation_history and "T cells" in annotation_history:
                print("  PASS: Integer cluster_name 0 matched string JSON key '0'")
                status = "passed"
            else:
                errors.append(f"Integer key 0 did not match string JSON key '0'. Got: {annotation_history!r}")
                status = "failed"

    except Exception as e:
        errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def test_llm_returns_none():
    """
    Test that iterative_marker_analysis handles None return from call_llm
    without crashing with 'argument of type NoneType is not iterable'.
    """
    print("\n" + "=" * 60)
    print("TEST 3: call_llm returns None — no TypeError crash (Issue #18)")
    print("=" * 60)

    errors = []
    status = "error"

    try:
        from unittest.mock import patch, MagicMock
        from CASSIA.agents.annotation_boost.annotation_boost import iterative_marker_analysis
        import pandas as pd

        marker_data = pd.DataFrame({
            "cluster": [0, 0],
            "gene": ["CD3D", "CD3E"],
            "p_val_adj": [0.001, 0.001],
            "avg_log2FC": [2.0, 2.0]
        })

        call_count = {"n": 0}

        def mock_call_llm(*args, **kwargs):
            call_count["n"] += 1
            if call_count["n"] == 1:
                # First call returns None (simulates API failure)
                return None
            else:
                # Second call returns a valid completion with the completion marker
                return "FINAL ANNOTATION\nT cells based on CD3D, CD3E markers."

        with patch('CASSIA.agents.annotation_boost.annotation_boost.call_llm', side_effect=mock_call_llm):
            try:
                result, messages = iterative_marker_analysis(
                    major_cluster_info="blood human immune cells",
                    marker=marker_data,
                    comma_separated_genes="CD3D, CD3E",
                    annotation_history="",
                    provider="openrouter",
                    model="test-model",
                    temperature=0.3,
                    num_iterations=3
                )

                if result is not None:
                    print(f"  Result obtained after {call_count['n']} calls (1st was None)")
                    print(f"  PASS: No TypeError crash when call_llm returns None")
                    status = "passed"
                else:
                    # Even if result is None after all iterations, the point is it didn't crash
                    print(f"  PASS: Function completed without TypeError (result={result!r})")
                    status = "passed"

            except TypeError as e:
                if "NoneType" in str(e):
                    errors.append(f"TypeError still occurs: {e}")
                    status = "failed"
                else:
                    raise

    except Exception as e:
        errors.append(str(e))
        print(f"\n  Error: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n  Status: {status.upper()}")
    return status == "passed", errors


def run_cluster_lookup_test():
    """Run all cluster lookup bugfix tests."""
    print_test_header("25 - Annotation Boost Cluster Lookup (Issue #18)")

    total_start = time.time()
    all_errors = []
    test_results = {}

    # Test 1: String key baseline
    print("\n" + "#" * 70)
    print("# BUG FIX 1: String cluster_name lookup (baseline)")
    print("#" * 70)
    success1, errors1 = test_json_key_lookup_string_cluster()
    test_results['string_key'] = {'passed': success1, 'errors': errors1}
    all_errors.extend(errors1)

    # Test 2: Integer key — the core Issue #18 bug
    print("\n" + "#" * 70)
    print("# BUG FIX 2: Integer cluster_name → string JSON key (Issue #18 core)")
    print("#" * 70)
    success2, errors2 = test_json_key_lookup_int_cluster()
    test_results['int_key'] = {'passed': success2, 'errors': errors2}
    all_errors.extend(errors2)

    # Test 3: NoneType crash
    print("\n" + "#" * 70)
    print("# BUG FIX 3: call_llm returns None — no crash (Issue #18 secondary)")
    print("#" * 70)
    success3, errors3 = test_llm_returns_none()
    test_results['none_crash'] = {'passed': success3, 'errors': errors3}
    all_errors.extend(errors3)

    total_duration = time.time() - total_start

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\n  Test 1 (String key lookup):       {'PASSED' if success1 else 'FAILED'}")
    print(f"  Test 2 (Int→String key fix):      {'PASSED' if success2 else 'FAILED'}")
    print(f"  Test 3 (None return handling):     {'PASSED' if success3 else 'FAILED'}")
    print(f"\n  Total Duration: {total_duration:.2f}s")

    if all_errors:
        print(f"\n  Errors:")
        for err in all_errors:
            print(f"    - {err[:100]}")

    all_passed = success1 and success2 and success3
    print_test_result(all_passed, f"Total Duration: {total_duration:.2f}s")

    return all_passed


if __name__ == "__main__":
    success = run_cluster_lookup_test()
    sys.exit(0 if success else 1)
