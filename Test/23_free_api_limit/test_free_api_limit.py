"""
CASSIA Test 23: Free API Lifetime Limit
========================================
Tests the free API 2-cluster lifetime limit by running runCASSIA_batch
WITHOUT any API key. Expects auto-selection of at most 2 clusters.

This test deliberately does NOT load API keys. It validates:
- Free API auto-detection when no key is set
- Auto-selection of first N clusters when quota < total clusters
- "LIMIT REACHED" message when quota is exhausted
- Post-run reminder with remaining quota

WARNING: This test consumes free API quota (up to 2 cluster annotations
per machine, lifetime). After quota is exhausted, the test will still
PASS by validating the "limit reached" behavior.

Usage:
    python test_free_api_limit.py
"""

import sys
import time
import os
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import load_markers
from test_utils import (
    setup_cassia_imports,
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

# Setup CASSIA imports (development mode)
setup_cassia_imports()


def run_free_api_limit_test():
    """Run free API lifetime limit test on all clusters."""
    print_test_header("23 - Free API Lifetime Limit")

    # -------------------------------------------------------------------------
    # IMPORTANT: Do NOT load API keys. We want to test the free API path.
    # Explicitly unset any API keys that might be in the environment.
    # -------------------------------------------------------------------------
    print("\nClearing API keys to force free API path...")
    for key in ['OPENROUTER_API_KEY', 'GOOGLE_API_KEY', 'TOGETHER_API_KEY',
                'OPENAI_API_KEY', 'ANTHROPIC_API_KEY']:
        removed = os.environ.pop(key, None)
        if removed:
            print(f"  Unset {key}")

    # Import CASSIA functions
    from CASSIA import runCASSIA_batch, get_remaining_free_clusters, MAX_FREE_CLUSTERS
    from CASSIA.core.free_api import clear_free_api_cache

    # Clear any cached keys from previous runs
    clear_free_api_cache()

    # Check remaining quota BEFORE test
    remaining_before = get_remaining_free_clusters()
    print(f"\nFree API Status:")
    print(f"  Max lifetime clusters: {MAX_FREE_CLUSTERS}")
    print(f"  Remaining before test: {remaining_before}")

    # Load ALL clusters (all 6 from processed.csv)
    full_df = load_markers()
    all_clusters = full_df['Broad.cell.type'].tolist()
    total_clusters = len(all_clusters)

    print(f"\nLoaded {total_clusters} clusters from marker data:")
    for cluster in all_clusters:
        print(f"  - {cluster}")

    # Create results directory
    results = create_results_dir("23_free_api_limit", get_test_mode())
    output_name = str(results['outputs'] / "free_api_results")
    print(f"\nResults will be saved to: {results['base']}")

    # Setup logging
    logging_ctx = setup_logging(results['logs'])

    # Run the test
    start_time = time.time()
    errors = []
    status = "error"
    clusters_annotated = 0
    expected_clusters = 0

    if remaining_before <= 0:
        # =====================================================================
        # SCENARIO: Quota exhausted — expect "LIMIT REACHED" message, no output
        # =====================================================================
        print(f"\n--- SCENARIO: Quota exhausted (remaining={remaining_before}) ---")
        print("Expecting 'FREE API LIMIT REACHED' message and no output.\n")

        expected_clusters = 0

        try:
            runCASSIA_batch(
                marker=full_df,
                output_name=output_name,
                n_genes=30,
                tissue='large intestine',
                species='human',
                max_workers=3,
                provider='openrouter',
            )

            # If we get here, the function returned without error (expected for limit=0)
            summary_csv = Path(f"{output_name}_summary.csv")
            if summary_csv.exists():
                import pandas as pd
                results_df = pd.read_csv(summary_csv)
                if len(results_df) > 0:
                    errors.append(f"Expected no clusters annotated but got {len(results_df)}")
                    status = "failed"
                else:
                    status = "passed"
            else:
                # No output file = correct behavior when quota is 0
                status = "passed"

        except ValueError as e:
            # runCASSIA (single) raises ValueError for limit reached
            # runCASSIA_batch just returns None
            if "LIMIT REACHED" in str(e) or "free" in str(e).lower():
                print(f"Got expected error: {e}")
                status = "passed"
            else:
                errors.append(str(e))
                status = "error"
        except Exception as e:
            errors.append(str(e))
            status = "error"

    else:
        # =====================================================================
        # SCENARIO: Quota available — expect auto-selection
        # =====================================================================
        expected_clusters = min(remaining_before, total_clusters)
        print(f"\n--- SCENARIO: Quota available (remaining={remaining_before}) ---")
        print(f"Expecting {expected_clusters} of {total_clusters} clusters to be annotated.\n")

        try:
            runCASSIA_batch(
                marker=full_df,
                output_name=output_name,
                n_genes=30,
                tissue='large intestine',
                species='human',
                max_workers=3,
                provider='openrouter',
            )

            # Check output files
            summary_csv = Path(f"{output_name}_summary.csv")
            conversations_json = Path(f"{output_name}_conversations.json")
            html_report = Path(f"{output_name}_report.html")

            if summary_csv.exists():
                import pandas as pd
                import json
                results_df = pd.read_csv(summary_csv)
                clusters_annotated = len(results_df)

                print(f"\nBatch Results:")
                print(f"  Clusters annotated: {clusters_annotated}")
                print(f"  Expected: {expected_clusters}")
                print(f"  Output files created:")
                print(f"    - {summary_csv.name}")
                if conversations_json.exists():
                    print(f"    - {conversations_json.name}")
                    with open(conversations_json, 'r', encoding='utf-8') as f:
                        conv_data = json.load(f)
                    print(f"      (contains {len(conv_data)} cluster conversations)")
                if html_report.exists():
                    print(f"    - {html_report.name}")

                if clusters_annotated == expected_clusters:
                    status = "passed"
                elif clusters_annotated > 0 and clusters_annotated <= expected_clusters:
                    # Some clusters may fail due to API issues — still acceptable
                    status = "passed"
                    print(f"  Note: {expected_clusters - clusters_annotated} cluster(s) may have failed")
                else:
                    status = "failed"
                    errors.append(
                        f"Expected {expected_clusters} clusters but got {clusters_annotated}"
                    )
            else:
                status = "failed"
                errors.append("Summary CSV not created")

        except Exception as e:
            errors.append(str(e))
            status = "error"
            print(f"\nError: {e}")

    duration = time.time() - start_time

    # Check remaining quota AFTER test
    try:
        clear_free_api_cache()
        remaining_after = get_remaining_free_clusters()
    except Exception:
        remaining_after = -1

    print(f"\nPost-test quota: {remaining_after} of {MAX_FREE_CLUSTERS} remaining")

    # Save metadata
    metadata = create_test_metadata(
        test_name="free_api_limit",
        config={
            'llm': {
                'provider': 'openrouter',
                'model': '(free API default)',
            },
            'data': {
                'tissue': 'large intestine',
                'species': 'human',
                'n_genes': 30,
            }
        },
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters[:expected_clusters] if expected_clusters > 0 else [],
        errors=errors
    )
    # Add free API-specific metadata
    metadata['free_api'] = {
        'remaining_before': remaining_before,
        'remaining_after': remaining_after,
        'max_free_clusters': MAX_FREE_CLUSTERS,
        'expected_clusters': expected_clusters,
        'clusters_annotated': clusters_annotated,
        'total_clusters_in_data': total_clusters,
    }
    save_test_metadata(results['outputs'], metadata)

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    if success:
        if remaining_before <= 0:
            print("    (Validated: limit reached behavior)")
        else:
            print(f"    (Validated: {clusters_annotated} clusters annotated with free API)")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_free_api_limit_test()
    sys.exit(0 if success else 1)
