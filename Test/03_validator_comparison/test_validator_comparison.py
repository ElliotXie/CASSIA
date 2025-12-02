"""
CASSIA Test 03: Validator Comparison
====================================
Compares v0 (strict) and v1 (moderate) validators.

Usage:
    python test_validator_comparison.py
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_full_marker_dataframe, get_all_clusters
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


def run_validator_comparison_test():
    """Compare v0 and v1 validators."""
    print_test_header("03 - Validator Comparison")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import runCASSIA_batch
    import pandas as pd

    # Get marker data
    marker_df = get_full_marker_dataframe()
    all_clusters = get_all_clusters()

    print(f"\nComparing validators on {len(all_clusters)} clusters")

    # Create results directory with organized structure
    results = create_results_dir("03_validator_comparison", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # Run with both validators
    validators = ["v0", "v1"]
    validator_results = {}
    errors = []

    start_time = time.time()

    for validator in validators:
        print(f"\n{'='*40}")
        print(f"Testing {validator} validator ({'strict' if validator == 'v0' else 'moderate'})")
        print(f"{'='*40}")

        output_name = str(results['outputs'] / f"results_{validator}")

        try:
            runCASSIA_batch(
                marker=marker_df,
                output_name=output_name,
                n_genes=data_config.get('n_genes', 30),
                model=llm_config.get('model', 'google/gemini-2.5-flash'),
                temperature=llm_config.get('temperature', 0.3),
                tissue=data_config.get('tissue', 'large intestine'),
                species=data_config.get('species', 'human'),
                max_workers=llm_config.get('max_workers', 3),
                provider=llm_config.get('provider', 'openrouter'),
                validator_involvement=validator
            )

            # Load results
            full_csv = Path(f"{output_name}_full.csv")
            if full_csv.exists():
                df = pd.read_csv(full_csv)
                validator_results[validator] = {
                    "clusters_annotated": len(df),
                    "results_file": str(full_csv)
                }
                print(f"  Annotated {len(df)} clusters")
            else:
                errors.append(f"{validator}: Output file not created")

        except Exception as e:
            errors.append(f"{validator}: {str(e)}")
            print(f"  Error: {e}")

    duration = time.time() - start_time

    # Compare results
    print(f"\n{'='*40}")
    print("Comparison Results")
    print(f"{'='*40}")

    comparison = {}
    if "v0" in validator_results and "v1" in validator_results:
        v0_count = validator_results["v0"]["clusters_annotated"]
        v1_count = validator_results["v1"]["clusters_annotated"]
        comparison = {
            "v0_clusters": v0_count,
            "v1_clusters": v1_count,
            "both_complete": v0_count == len(all_clusters) and v1_count == len(all_clusters)
        }
        print(f"  v0 (strict): {v0_count} clusters")
        print(f"  v1 (moderate): {v1_count} clusters")
        status = "passed" if comparison["both_complete"] else "failed"
    else:
        status = "failed"
        print("  Could not compare - one or both validators failed")

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="validator_comparison",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters,
        errors=errors
    )
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "validator_results": validator_results,
        "comparison": comparison
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_validator_comparison_test()
    sys.exit(0 if success else 1)
