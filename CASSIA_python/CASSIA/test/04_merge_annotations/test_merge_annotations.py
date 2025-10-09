#!/usr/bin/env python3
"""
Test: merge_annotations_all
Description: Test annotation merging at multiple granularity levels
Expected Runtime: 3-5 minutes

This test validates the annotation merging functionality which consolidates
cell type annotations at three different granularity levels:
1. Broad (Merged_Grouping_1): General lineage categories
2. Detailed (Merged_Grouping_2): Intermediate specificity
3. Very Detailed (Merged_Grouping_3): Normalized specific annotations

Usage:
    python test_merge_annotations.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import shutil

# Add parent directories to path for imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Import CASSIA
import CASSIA

# Import shared utilities
shared_dir = Path(__file__).parent.parent / "shared"
sys.path.insert(0, str(shared_dir))

from test_config import load_config, validate_api_key
from test_utils import (
    setup_logging,
    log_test_start,
    log_test_end,
    save_results,
    Timer,
    get_timestamp
)
from sample_data import SampleDataLoader


def create_prerequisite_batch_results(logger, config, marker_data):
    """
    Create the prerequisite batch annotation results needed for merging.

    Returns:
        str: Path to the batch results CSV
    """
    logger.info("\n" + "=" * 80)
    logger.info("CREATING PREREQUISITE BATCH RESULTS")
    logger.info("=" * 80)

    batch_output = config['prerequisites']['batch_results_file'].replace('_full.csv', '')

    logger.info(f"Running batch annotation to create prerequisite data...")
    logger.info(f"Output name: {batch_output}")

    with Timer("Prerequisite batch annotation", logger):
        CASSIA.runCASSIA_batch(
            marker=marker_data,
            output_name=batch_output,
            model=config['model'],
            provider=config['provider'],
            tissue="large intestine",
            species="human",
            max_workers=4,
            n_genes=50
        )

    batch_csv = f"{batch_output}_full.csv"

    if not os.path.exists(batch_csv):
        raise FileNotFoundError(f"Batch results not created: {batch_csv}")

    logger.info(f"✓ Created prerequisite batch results: {batch_csv}")
    return batch_csv


def validate_merged_columns(merged_df, logger, config):
    """Validate that all expected merged annotation columns exist."""

    logger.info("Validating merged annotation columns:")

    expected_cols = config['validation']['expected_merged_columns']
    missing_cols = []

    for col in expected_cols:
        if col in merged_df.columns:
            non_null = merged_df[col].notna().sum()
            logger.info(f"  ✓ Found {col}: {non_null}/{len(merged_df)} non-null values")
        else:
            missing_cols.append(col)
            logger.error(f"  ✗ Missing column: {col}")

    if missing_cols:
        raise ValueError(f"Missing merged annotation columns: {missing_cols}")

    return True


def analyze_grouping_hierarchy(merged_df, logger, config):
    """Analyze and validate the hierarchy of grouping granularities."""

    logger.info("\nAnalyzing grouping hierarchy:")

    grouping_cols = [
        'Merged_Grouping_1',  # Broad
        'Merged_Grouping_2',  # Detailed
        'Merged_Grouping_3'   # Very detailed
    ]

    # Check that each level exists
    for col in grouping_cols:
        if col not in merged_df.columns:
            logger.warning(f"  ⚠ Missing column: {col}")
            return False

    # Count unique values at each level
    unique_counts = {}
    for col in grouping_cols:
        unique_counts[col] = merged_df[col].nunique()
        logger.info(f"  - {col}: {unique_counts[col]} unique groups")

    # Validate hierarchy: Grouping_1 should have fewest unique values (most general)
    # and Grouping_3 should have most (most specific)
    if config['validation']['check_grouping_hierarchy']:
        if unique_counts['Merged_Grouping_1'] <= unique_counts['Merged_Grouping_2'] <= unique_counts['Merged_Grouping_3']:
            logger.info("  ✓ Hierarchy validation passed: Broad → Detailed → Very Detailed")
        else:
            logger.warning("  ⚠ Hierarchy may not follow expected pattern")
            logger.warning(f"    Expected: Grouping_1 ≤ Grouping_2 ≤ Grouping_3")
            logger.warning(f"    Got: {unique_counts['Merged_Grouping_1']} ≤ {unique_counts['Merged_Grouping_2']} ≤ {unique_counts['Merged_Grouping_3']}")

    return True


def show_grouping_examples(merged_df, logger, config):
    """Display example groupings at each level."""

    if not config['reporting']['show_grouping_examples']:
        return

    logger.info("\nExample groupings by level:")
    logger.info("-" * 80)

    # Show first few rows with all grouping levels
    display_cols = ['True Cell Type', 'Predicted Main Cell Type']
    grouping_cols = ['Merged_Grouping_1', 'Merged_Grouping_2', 'Merged_Grouping_3']

    for col in grouping_cols:
        if col in merged_df.columns:
            display_cols.append(col)

    for idx, row in merged_df.head(5).iterrows():
        logger.info(f"\nCluster: {row['True Cell Type']}")
        logger.info(f"  Predicted: {row['Predicted Main Cell Type']}")

        if 'Merged_Grouping_1' in row.index:
            logger.info(f"  Broad (G1): {row['Merged_Grouping_1']}")
        if 'Merged_Grouping_2' in row.index:
            logger.info(f"  Detailed (G2): {row['Merged_Grouping_2']}")
        if 'Merged_Grouping_3' in row.index:
            logger.info(f"  Very Detailed (G3): {row['Merged_Grouping_3']}")


def main():
    """Main test function."""
    # Setup
    test_dir = Path(__file__).parent
    config = load_config(test_dir / "config.json")

    # Setup logging with file output
    timestamp = get_timestamp()
    log_file = test_dir / "results" / f"{timestamp}_test_log.txt"
    logger = setup_logging(config['test_name'], log_file=log_file)

    logger.info("=" * 80)
    logger.info("CASSIA Test: merge_annotations_all (Multi-Level Merging)")
    logger.info("=" * 80)

    # Validate API key
    try:
        api_key = validate_api_key(config['provider'])
        CASSIA.set_api_key(api_key, provider=config['provider'])
        logger.info(f"✓ API key validated for provider: {config['provider']}")
    except EnvironmentError as e:
        logger.error(str(e))
        return 1

    # Load data
    logger.info(f"\nLoading data: {config['data_file']}")
    try:
        loader = SampleDataLoader()
        dataset_name = config['data_file'].replace('.csv', '')
        marker_data = loader.load_processed() if dataset_name == 'processed' else loader.load_unprocessed()

        logger.info(f"✓ Loaded {len(marker_data)} rows from {config['data_file']}")

        # Get cluster info
        cluster_col = marker_data.columns[1]
        clusters = marker_data[cluster_col].unique()
        logger.info(f"  - Number of clusters: {len(clusters)}")
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}", exc_info=True)
        return 1

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    try:
        # Step 1: Create prerequisite batch results
        if config['prerequisites']['requires_batch_results']:
            batch_csv = create_prerequisite_batch_results(logger, config, marker_data)
        else:
            batch_csv = config['prerequisites']['batch_results_file']
            if not os.path.exists(batch_csv):
                raise FileNotFoundError(f"Batch results file not found: {batch_csv}")

        # Step 2: Run annotation merging
        logger.info("\n" + "=" * 80)
        logger.info("RUNNING ANNOTATION MERGING")
        logger.info("=" * 80)
        logger.info(f"\nMerge configuration:")
        logger.info(f"  - Input file: {batch_csv}")
        logger.info(f"  - Batch size: {config['method_params']['batch_size']}")
        logger.info(f"  - Context: {config['method_params']['additional_context']}")
        logger.info(f"  - Levels: Broad, Detailed, Very Detailed (3 parallel calls)")

        # Import merge_annotations_all
        try:
            from CASSIA.merging_annotation import merge_annotations_all
        except ImportError:
            from merging_annotation import merge_annotations_all

        with Timer("Annotation merging (all levels)", logger):
            merged_df = merge_annotations_all(
                csv_path=batch_csv,
                output_path=None,  # We'll save it ourselves with timestamp
                provider=config['provider'],
                model=config['model'],
                api_key=api_key,
                additional_context=config['method_params']['additional_context'],
                batch_size=config['method_params']['batch_size']
            )

        logger.info(f"\n✓ Annotation merging completed")
        logger.info(f"  - Result shape: {merged_df.shape}")
        logger.info(f"  - Columns added: {[col for col in merged_df.columns if 'Merged_Grouping' in col]}")

        # Step 3: Validate merged columns
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATING MERGED ANNOTATIONS")
        logger.info("=" * 80)

        validate_merged_columns(merged_df, logger, config)

        # Step 4: Analyze grouping hierarchy
        analyze_grouping_hierarchy(merged_df, logger, config)

        # Step 5: Show examples
        show_grouping_examples(merged_df, logger, config)

        # Step 6: Archive results
        logger.info("\n" + "=" * 80)
        logger.info("ARCHIVING RESULTS")
        logger.info("=" * 80)

        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Save merged results with timestamp
        merged_path = save_results(
            merged_df,
            config['test_name'],
            results_dir,
            prefix="merged",
            suffix="annotations"
        )
        logger.info(f"  ✓ Archived merged annotations: {merged_path.name}")

        # Clean up batch results
        try:
            if os.path.exists(batch_csv):
                os.remove(batch_csv)
            batch_summary = batch_csv.replace('_full.csv', '_summary.csv')
            if os.path.exists(batch_summary):
                os.remove(batch_summary)
            logger.info(f"  ✓ Cleaned up prerequisite files")
        except:
            pass

        # Create test summary
        summary = {
            'total_clusters': len(merged_df),
            'merged_columns': [col for col in merged_df.columns if 'Merged_Grouping' in col],
            'unique_groups': {
                'Merged_Grouping_1': merged_df['Merged_Grouping_1'].nunique() if 'Merged_Grouping_1' in merged_df.columns else 0,
                'Merged_Grouping_2': merged_df['Merged_Grouping_2'].nunique() if 'Merged_Grouping_2' in merged_df.columns else 0,
                'Merged_Grouping_3': merged_df['Merged_Grouping_3'].nunique() if 'Merged_Grouping_3' in merged_df.columns else 0
            }
        }

        summary_path = save_results(
            summary,
            config['test_name'],
            results_dir,
            prefix="merge",
            suffix="summary"
        )
        logger.info(f"  ✓ Saved test summary: {summary_path.name}")

        # Log final summary
        logger.info("\n" + "=" * 80)
        logger.info("MERGING RESULTS SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Clusters processed: {len(merged_df)}")
        logger.info(f"Merging levels completed: 3")
        logger.info(f"Unique groups per level:")
        logger.info(f"  - Broad (Grouping_1): {summary['unique_groups']['Merged_Grouping_1']}")
        logger.info(f"  - Detailed (Grouping_2): {summary['unique_groups']['Merged_Grouping_2']}")
        logger.info(f"  - Very Detailed (Grouping_3): {summary['unique_groups']['Merged_Grouping_3']}")

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ ANNOTATION MERGING TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Results directory: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ Merging test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
