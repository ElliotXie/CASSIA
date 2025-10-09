#!/usr/bin/env python3
"""
Test: runCASSIA_subclusters
Description: Test hierarchical subcluster annotation functionality
Expected Runtime: 3-5 minutes

This test validates the subclustering functionality which annotates
sub-populations within a major cluster type.

Usage:
    python test_subclustering.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample subcluster marker data
"""

import os
import sys
from pathlib import Path
import pandas as pd
import glob

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


def create_synthetic_subcluster_data(logger):
    """
    Create synthetic subcluster marker data for testing.

    Returns:
        pd.DataFrame: Subcluster marker data
    """
    logger.info("Creating synthetic subcluster marker data...")

    # Create synthetic CD8 T cell subcluster data
    subclusters = pd.DataFrame({
        'Subcluster': [
            'Subcluster_0',
            'Subcluster_1',
            'Subcluster_2',
            'Subcluster_3'
        ],
        'Top_Markers': [
            'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK, LEF1, SELL',
            'LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3, PDCD1',
            'GZMK, GZMH, PRF1, NKG7, CCR7, CD27, GZMA, GNLY',
            'WFDC2, CEACAM7, CLDN8, PPARG, MKI67, PCNA, TOP2A'
        ]
    })

    logger.info(f"  ✓ Created {len(subclusters)} synthetic subclusters")
    for idx, row in subclusters.iterrows():
        logger.info(f"    - {row['Subcluster']}: {row['Top_Markers'][:50]}...")

    return subclusters


def validate_subcluster_results(result_df, logger, config):
    """Validate subcluster annotation results."""

    logger.info("Validating subcluster results:")

    # Check expected columns
    expected_cols = config['validation']['expected_columns']
    missing_cols = set(expected_cols) - set(result_df.columns)

    if missing_cols:
        logger.error(f"  ✗ Missing columns: {missing_cols}")
        raise ValueError(f"Missing required columns: {missing_cols}")

    logger.info(f"  ✓ All expected columns present: {expected_cols}")

    # Check for non-null cell type values
    if config['validation']['check_cell_type_values']:
        null_main = result_df['main_cell_type'].isna().sum()
        null_sub = result_df['sub_cell_type'].isna().sum()

        if null_main > 0:
            logger.warning(f"  ⚠ {null_main} null main_cell_type values")
        else:
            logger.info(f"  ✓ All main_cell_type values non-null")

        if null_sub > 0:
            logger.warning(f"  ⚠ {null_sub} null sub_cell_type values")
        else:
            logger.info(f"  ✓ All sub_cell_type values non-null")

    # Check for reasons/explanations
    if config['validation']['check_reason_provided']:
        if 'reason' in result_df.columns:
            has_reason = result_df['reason'].notna().sum()
            logger.info(f"  ✓ Reasons provided: {has_reason}/{len(result_df)} subclusters")
        else:
            logger.warning(f"  ⚠ 'reason' column not found")

    return True


def display_subcluster_results(result_df, logger):
    """Display annotated subcluster results."""

    logger.info("\nSubcluster annotation results:")
    logger.info("-" * 80)

    for idx, row in result_df.iterrows():
        logger.info(f"\n{row.get('Result ID', idx+1)}. Subcluster annotation:")
        logger.info(f"   Main type: {row['main_cell_type']}")
        logger.info(f"   Sub type: {row['sub_cell_type']}")

        if 'key_markers' in row.index and pd.notna(row['key_markers']):
            logger.info(f"   Key markers: {row['key_markers'][:60]}...")

        if 'reason' in row.index and pd.notna(row['reason']) and row['reason']:
            reason_preview = row['reason'][:100] + "..." if len(str(row['reason'])) > 100 else row['reason']
            logger.info(f"   Reason: {reason_preview}")


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
    logger.info("CASSIA Test: runCASSIA_subclusters (Hierarchical Annotation)")
    logger.info("=" * 80)

    # Validate API key
    try:
        api_key = validate_api_key(config['provider'])
        CASSIA.set_api_key(api_key, provider=config['provider'])
        logger.info(f"✓ API key validated for provider: {config['provider']}")
    except EnvironmentError as e:
        logger.error(str(e))
        return 1

    # Create/load subcluster data
    logger.info(f"\nPreparing subcluster data:")
    try:
        # Try to load from data directory first
        data_dir = Path(__file__).parent.parent.parent / "data"
        subcluster_file = data_dir / "subcluster_results.csv"

        if subcluster_file.exists():
            logger.info(f"Loading subcluster data from: {subcluster_file}")
            subcluster_data = pd.read_csv(subcluster_file)
            logger.info(f"✓ Loaded {len(subcluster_data)} subclusters from file")
        else:
            logger.info("Subcluster data file not found, creating synthetic data")
            subcluster_data = create_synthetic_subcluster_data(logger)

    except Exception as e:
        logger.warning(f"Could not load data file: {e}")
        logger.info("Creating synthetic subcluster data...")
        subcluster_data = create_synthetic_subcluster_data(logger)

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    try:
        # Run subcluster annotation
        logger.info("\n" + "=" * 80)
        logger.info("RUNNING SUBCLUSTER ANNOTATION")
        logger.info("=" * 80)
        logger.info(f"\nSubcluster configuration:")
        logger.info(f"  - Major cluster: {config['method_params']['major_cluster_info']}")
        logger.info(f"  - Number of subclusters: {len(subcluster_data)}")
        logger.info(f"  - Model: {config['model']}")
        logger.info(f"  - Number of genes: {config['method_params']['n_genes']}")

        output_name = config['method_params']['output_name']

        with Timer("Subcluster annotation", logger):
            # Import subclustering module
            from CASSIA.subclustering import runCASSIA_subclusters

            runCASSIA_subclusters(
                marker=subcluster_data,
                major_cluster_info=config['method_params']['major_cluster_info'],
                output_name=output_name,
                model=config['model'],
                temperature=config['temperature'],
                provider=config['provider'],
                n_genes=config['method_params']['n_genes']
            )

        # Read results
        result_csv = f"{output_name}.csv" if not output_name.endswith('.csv') else output_name

        if not os.path.exists(result_csv):
            raise FileNotFoundError(f"Subcluster results not created: {result_csv}")

        result_df = pd.read_csv(result_csv)
        logger.info(f"\n✓ Subcluster annotation completed")
        logger.info(f"  - Results shape: {result_df.shape}")
        logger.info(f"  - Output file: {result_csv}")

        # Validate results
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATING SUBCLUSTER RESULTS")
        logger.info("=" * 80)

        validate_subcluster_results(result_df, logger, config)

        # Display results
        display_subcluster_results(result_df, logger)

        # Check for HTML report
        if config['reporting']['check_html_report']:
            html_file = f"{output_name}.html" if not output_name.endswith('.csv') else output_name.replace('.csv', '.html')

            if os.path.exists(html_file):
                logger.info(f"\n✓ HTML report generated: {html_file}")
            else:
                logger.warning(f"\n⚠ HTML report not found: {html_file}")

        # Archive results
        logger.info("\n" + "=" * 80)
        logger.info("ARCHIVING RESULTS")
        logger.info("=" * 80)

        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Save results with timestamp
        archived_path = save_results(
            result_df,
            config['test_name'],
            results_dir,
            prefix="subcluster",
            suffix="annotations"
        )
        logger.info(f"  ✓ Archived results: {archived_path.name}")

        # Copy HTML report if exists
        if os.path.exists(html_file):
            import shutil
            html_dest = results_dir / f"{timestamp}_subcluster_report.html"
            shutil.copy2(html_file, html_dest)
            logger.info(f"  ✓ Archived HTML report: {html_dest.name}")

        # Clean up original files
        try:
            if os.path.exists(result_csv):
                os.remove(result_csv)
            if os.path.exists(html_file):
                os.remove(html_file)
            logger.info(f"  ✓ Cleaned up working files")
        except:
            pass

        # Create test summary
        summary = {
            'major_cluster': config['method_params']['major_cluster_info'],
            'n_subclusters': len(result_df),
            'main_cell_types': result_df['main_cell_type'].tolist(),
            'sub_cell_types': result_df['sub_cell_type'].tolist()
        }

        summary_path = save_results(
            summary,
            config['test_name'],
            results_dir,
            prefix="subcluster",
            suffix="summary"
        )
        logger.info(f"  ✓ Saved test summary: {summary_path.name}")

        # Log final summary
        logger.info("\n" + "=" * 80)
        logger.info("SUBCLUSTER ANNOTATION SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Major cluster: {config['method_params']['major_cluster_info']}")
        logger.info(f"Subclusters annotated: {len(result_df)}")
        logger.info(f"Unique main types: {result_df['main_cell_type'].nunique()}")
        logger.info(f"Unique sub types: {result_df['sub_cell_type'].nunique()}")

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ SUBCLUSTER ANNOTATION TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Results directory: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ Subcluster test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
