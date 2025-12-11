"""
CASSIA Test 04: Quality Scoring
===============================
Tests the runCASSIA_score_batch function to score annotation results.

Usage:
    python test_quality_scoring.py

Note: This test requires batch annotation results from Test 02.
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
    get_latest_results,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def run_quality_scoring_test():
    """Test quality scoring functionality."""
    print_test_header("04 - Quality Scoring")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import runCASSIA_batch, runCASSIA_score_batch
    import pandas as pd

    # Create results directory with organized structure
    results = create_results_dir("04_quality_scoring", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # First, we need batch results to score
    # Check if we have existing batch results from Test 02
    batch_results_dir = get_latest_results("02_batch_annotation")
    input_file = None

    if batch_results_dir:
        potential_file = batch_results_dir / "batch_results_summary.csv"
        if potential_file.exists():
            input_file = str(potential_file)
            print(f"\nUsing existing batch results: {input_file}")

    # If no existing results, run a quick batch annotation
    if not input_file:
        print("\nNo existing batch results found. Running batch annotation first...")
        marker_df = get_full_marker_dataframe()
        batch_output = str(results['outputs'] / "batch_for_scoring")

        runCASSIA_batch(
            marker=marker_df,
            output_name=batch_output,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1')
        )
        input_file = f"{batch_output}_summary.csv"

    # Run quality scoring
    start_time = time.time()
    errors = []
    status = "error"
    scoring_results = {}

    output_file = str(results['outputs'] / "scored_results.csv")

    try:
        print(f"\nRunning quality scoring...")
        print(f"Input: {input_file}")
        print(f"Output: {output_file}")

        runCASSIA_score_batch(
            input_file=input_file,
            output_file=output_file,
            max_workers=llm_config.get('max_workers', 3),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            provider=llm_config.get('provider', 'openrouter')
        )

        # Check results
        if Path(output_file).exists():
            scored_df = pd.read_csv(output_file)

            # Check for score column
            score_col = None
            for col in ['score', 'Score', 'SCORE']:
                if col in scored_df.columns:
                    score_col = col
                    break

            if score_col:
                scores = scored_df[score_col].dropna()
                scoring_results = {
                    "clusters_scored": len(scores),
                    "avg_score": float(scores.mean()),
                    "min_score": float(scores.min()),
                    "max_score": float(scores.max()),
                    "scores": scores.tolist()
                }

                print(f"\nScoring Results:")
                print(f"  Clusters scored: {len(scores)}")
                print(f"  Average score: {scores.mean():.1f}")
                print(f"  Score range: {scores.min():.0f} - {scores.max():.0f}")

                # Validate scores are in expected range
                if scores.min() >= 0 and scores.max() <= 100:
                    status = "passed"
                else:
                    status = "failed"
                    errors.append("Scores outside expected range (0-100)")
            else:
                status = "failed"
                errors.append("No score column found in output")
        else:
            status = "failed"
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")

    duration = time.time() - start_time

    # Save metadata and results
    all_clusters = get_all_clusters()
    metadata = create_test_metadata(
        test_name="quality_scoring",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters,
        errors=errors
    )
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "input_file": input_file,
        "output_file": output_file,
        "scoring_results": scoring_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_quality_scoring_test()
    sys.exit(0 if success else 1)
