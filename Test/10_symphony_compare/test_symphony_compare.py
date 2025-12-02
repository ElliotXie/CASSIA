"""
CASSIA Test 10: Symphony Compare
=================================
Tests the symphonyCompare function for multi-model cell type comparison
with AI consensus building.

Usage:
    python test_symphony_compare.py

Functions tested:
- symphonyCompare(): Orchestrate multiple AI models to compare cell types
"""

import sys
import time
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from fixtures import get_cluster_markers
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


def run_symphony_compare_test():
    """Test Symphony Compare functionality."""
    print_test_header("10 - Symphony Compare")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    # Get LLM settings
    llm_config = config['llm']
    data_config = config['data']

    # Import CASSIA functions
    from CASSIA import symphonyCompare

    # Create results directory with organized structure
    results = create_results_dir("10_symphony_compare", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging to capture console output
    logging_ctx = setup_logging(results['logs'])

    # Test parameters - compare cell types with a marker set
    test_cluster = "plasma cell"
    markers = get_cluster_markers(test_cluster)
    marker_set = ", ".join(markers[:15])  # Use top 15 markers

    # Cell types to compare
    celltypes = ["Plasma cell", "B cell", "T cell"]

    # Run tests
    start_time = time.time()
    errors = []
    status = "error"
    symphony_results = {}

    try:
        print(f"\n--- Test: symphonyCompare ---")
        print(f"  Tissue: {data_config.get('tissue', 'large intestine')}")
        print(f"  Species: {data_config.get('species', 'human')}")
        print(f"  Cell types to compare: {', '.join(celltypes)}")
        print(f"  Marker set: {marker_set[:50]}...")
        print(f"  Model preset: budget (cost-effective for testing)")
        print(f"  Discussion: Disabled (for faster testing)")

        result = symphonyCompare(
            tissue=data_config.get('tissue', 'large intestine'),
            celltypes=celltypes,
            marker_set=marker_set,
            species=data_config.get('species', 'human'),
            model_preset="budget",  # Use budget preset for cost-effective testing
            output_dir=str(results['outputs']),
            output_basename="symphony_test",
            enable_discussion=False,  # Disable discussion for faster testing
            max_discussion_rounds=0,
            consensus_threshold=0.6,
            generate_report=True,
            verbose=True
        )

        # Check result structure
        if isinstance(result, dict):
            print(f"\n--- Symphony Compare Results ---")
            print(f"  Consensus: {result.get('consensus', 'No consensus')}")
            print(f"  Confidence: {result.get('confidence', 0):.1%}")
            print(f"  CSV file: {Path(result.get('csv_file', '')).name if result.get('csv_file') else 'N/A'}")
            print(f"  HTML report: {Path(result.get('html_file', '')).name if result.get('html_file') else 'N/A'}")

            # Summary statistics
            summary = result.get('summary', {})
            if summary:
                print(f"\n  Summary:")
                print(f"    Models used: {summary.get('models_used', 'N/A')}")
                print(f"    Total rounds: {summary.get('total_rounds', 'N/A')}")
                print(f"    Consensus reached: {summary.get('consensus_reached', False)}")

            # Cell type scores
            celltype_scores = summary.get('celltype_scores', {})
            if celltype_scores:
                print(f"\n  Cell Type Scores:")
                for ct, scores in celltype_scores.items():
                    print(f"    {ct}: mean={scores.get('mean', 0):.1f}, range=[{scores.get('min', 0):.1f}-{scores.get('max', 0):.1f}]")

            symphony_results = {
                'consensus': result.get('consensus'),
                'confidence': result.get('confidence'),
                'csv_file': result.get('csv_file'),
                'html_file': result.get('html_file'),
                'summary': summary,
                'num_results': len(result.get('results', []))
            }

            # Validate results
            if result.get('csv_file') and Path(result['csv_file']).exists():
                status = "passed"
            else:
                status = "failed"
                errors.append("CSV file not created or consensus not reached")
        else:
            status = "failed"
            errors.append("Unexpected result format")

    except Exception as e:
        errors.append(str(e))
        status = "error"
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time

    # Save metadata and results
    metadata = create_test_metadata(
        test_name="symphony_compare",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=celltypes,
        errors=errors
    )
    save_test_metadata(results['outputs'], metadata)

    save_test_results(results['outputs'], {
        "celltypes_compared": celltypes,
        "marker_set": marker_set,
        "model_preset": "budget",
        "results": symphony_results
    })

    # Print final result
    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    # Cleanup logging
    cleanup_logging(logging_ctx)

    return success


if __name__ == "__main__":
    success = run_symphony_compare_test()
    sys.exit(0 if success else 1)
