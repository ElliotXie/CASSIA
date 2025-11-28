"""
CASSIA Test 12: Batch Annotation with Reference Retrieval
=========================================================
Tests the runCASSIA_batch_with_reference function with use_reference=True.
Verifies the two-step ReAct reference workflow.
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
    print_config_summary
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    create_test_metadata
)

setup_cassia_imports()

def run_batch_with_reference_test():
    print_test_header("12 - Batch Annotation with Reference")

    config = load_config()
    print_config_summary(config)
    setup_api_keys()

    llm_config = config['llm']
    data_config = config['data']

    from CASSIA import runCASSIA_batch_with_reference

    all_clusters = get_all_clusters()
    marker_df = get_full_marker_dataframe()

    results_dir = create_results_dir("12_batch_with_reference")
    output_name = str(results_dir / "batch_ref_results")

    start_time = time.time()
    errors = []
    status = "error"

    try:
        print("\nRunning runCASSIA_batch_with_reference with use_reference=True...")
        runCASSIA_batch_with_reference(
            marker=marker_df,
            output_name=output_name,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'large intestine'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1'),
            use_reference=True,  # KEY: Enable reference retrieval
            verbose=True
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")

        if full_csv.exists():
            import pandas as pd
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            # Verify reference columns exist
            has_ref_column = 'Reference Used' in results_df.columns
            refs_used = results_df['Reference Used'].value_counts() if has_ref_column else {}

            print(f"\nBatch with Reference Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(all_clusters)}")
            print(f"  Reference column present: {has_ref_column}")
            if has_ref_column:
                print(f"  Reference usage: {dict(refs_used)}")

            if clusters_annotated == len(all_clusters) and has_ref_column:
                status = "passed"
            else:
                status = "failed"
                if not has_ref_column:
                    errors.append("Reference Used column not found")
        else:
            status = "failed"
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")

    duration = time.time() - start_time

    metadata = create_test_metadata(
        test_name="batch_with_reference",
        config=config,
        duration_seconds=duration,
        status=status,
        clusters_tested=all_clusters,
        errors=errors
    )
    metadata['use_reference'] = True  # Additional metadata
    save_test_metadata(results_dir, metadata)

    success = status == "passed"
    print_test_result(success, f"Duration: {duration:.2f}s")

    return success

if __name__ == "__main__":
    success = run_batch_with_reference_test()
    sys.exit(0 if success else 1)
