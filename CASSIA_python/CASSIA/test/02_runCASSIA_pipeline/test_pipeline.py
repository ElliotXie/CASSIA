#!/usr/bin/env python3
"""
Test: runCASSIA_pipeline
Description: Test full end-to-end CASSIA pipeline
Expected Runtime: 10-15 minutes

This test validates the complete CASSIA pipeline including:
- Batch annotation
- Quality scoring
- Annotation boost for low-scoring clusters
- Annotation merging
- HTML report generation

Usage:
    python test_pipeline.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import glob
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
from sample_data import load_sample_data


def find_pipeline_output_dir(base_name, tissue, species):
    """Find the pipeline output directory with timestamp."""
    pattern = f"CASSIA_{tissue}_{species}_*".replace(" ", "_")
    matching_dirs = glob.glob(pattern)

    if matching_dirs:
        # Return the most recent one (sorted by name, which includes timestamp)
        return sorted(matching_dirs)[-1]
    return None


def validate_pipeline_structure(output_dir, logger):
    """Validate that all expected directories and files were created."""
    expected_dirs = [
        "01_annotation_results",
        "02_reports",
        "03_boost_analysis"
    ]

    for expected_dir in expected_dirs:
        dir_path = os.path.join(output_dir, expected_dir)
        if os.path.exists(dir_path):
            logger.info(f"  ✓ Found directory: {expected_dir}")
        else:
            raise FileNotFoundError(f"Expected directory not found: {expected_dir}")

    # Check for final results file
    final_results_pattern = os.path.join(output_dir, "01_annotation_results", "*_FINAL_RESULTS.csv")
    final_results = glob.glob(final_results_pattern)

    if final_results:
        logger.info(f"  ✓ Found final results: {os.path.basename(final_results[0])}")
        return final_results[0]
    else:
        raise FileNotFoundError("Final results CSV not found")


def validate_reports(output_dir, logger):
    """Check that HTML reports were generated."""
    reports_dir = os.path.join(output_dir, "02_reports")
    html_files = glob.glob(os.path.join(reports_dir, "*.html"))

    if html_files:
        logger.info(f"  ✓ Found {len(html_files)} HTML report(s)")
        for html_file in html_files:
            logger.info(f"    - {os.path.basename(html_file)}")
        return True
    else:
        logger.warning("  ⚠ No HTML reports found")
        return False


def validate_boost_results(output_dir, logger):
    """Check if boost analysis was performed for low-scoring clusters."""
    boost_dir = os.path.join(output_dir, "03_boost_analysis")

    if not os.path.exists(boost_dir):
        logger.info("  ℹ No boost analysis directory (all scores may be high)")
        return True

    # Check for subdirectories (one per boosted cluster)
    boost_subdirs = [d for d in os.listdir(boost_dir) if os.path.isdir(os.path.join(boost_dir, d))]

    if boost_subdirs:
        logger.info(f"  ✓ Found boost analysis for {len(boost_subdirs)} cluster(s)")
        for subdir in boost_subdirs:
            logger.info(f"    - {subdir}")
        return True
    else:
        logger.info("  ℹ No boost analysis performed (all scores above threshold)")
        return True


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
    logger.info("CASSIA Test: runCASSIA_pipeline (Full End-to-End)")
    logger.info("=" * 80)

    # Validate API key
    try:
        api_key = validate_api_key(config['annotation_provider'])
        CASSIA.set_api_key(api_key, provider=config['annotation_provider'])
        logger.info(f"✓ API key validated for provider: {config['annotation_provider']}")
    except EnvironmentError as e:
        logger.error(str(e))
        return 1

    # Load data
    logger.info(f"\nLoading data: {config['data_file']}")
    try:
        dataset_name = config['data_file'].replace('.csv', '')
        marker_data = load_sample_data(dataset_name)
        logger.info(f"✓ Loaded {len(marker_data)} rows from {config['data_file']}")

        if len(marker_data.columns) >= 2:
            cluster_col = marker_data.columns[1]
            n_clusters = marker_data[cluster_col].nunique()
            logger.info(f"  - Number of clusters: {n_clusters}")
            logger.info(f"  - Cluster names: {marker_data[cluster_col].unique().tolist()}")
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}", exc_info=True)
        return 1

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    # Store original directory to cleanup later
    original_dir = os.getcwd()

    try:
        # Run CASSIA pipeline
        logger.info("\n" + "=" * 80)
        logger.info("STARTING FULL CASSIA PIPELINE")
        logger.info("=" * 80)
        logger.info("\nPipeline stages:")
        logger.info("  1. Batch annotation")
        logger.info("  2. Quality scoring")
        logger.info("  3. Annotation boost (for low-scoring clusters)")
        logger.info("  4. Annotation merging")
        logger.info("  5. Report generation")
        logger.info("")

        with Timer("Complete CASSIA pipeline", logger):
            CASSIA.runCASSIA_pipeline(
                marker=marker_data,
                annotation_model=config['annotation_model'],
                annotation_provider=config['annotation_provider'],
                score_model=config['score_model'],
                score_provider=config['score_provider'],
                annotationboost_model=config['annotationboost_model'],
                annotationboost_provider=config['annotationboost_provider'],
                merge_model=config['merge_model'],
                merge_provider=config['merge_provider'],
                **config['method_params']
            )

        # Find the output directory
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATING PIPELINE OUTPUT")
        logger.info("=" * 80)

        output_dir = find_pipeline_output_dir(
            config['method_params']['output_file_name'],
            config['method_params']['tissue'],
            config['method_params']['species']
        )

        if not output_dir:
            raise FileNotFoundError("Pipeline output directory not found")

        logger.info(f"\n✓ Found pipeline output directory: {output_dir}")

        # Validate directory structure
        logger.info("\nValidating directory structure:")
        final_results_path = validate_pipeline_structure(output_dir, logger)

        # Read and validate final results
        logger.info("\nValidating final results:")
        final_df = pd.read_csv(final_results_path)
        logger.info(f"  ✓ Final results shape: {final_df.shape}")

        if config['validation']['check_output_format']:
            expected_cols = config['validation']['expected_final_columns']
            missing_cols = set(expected_cols) - set(final_df.columns)
            if missing_cols:
                raise ValueError(f"Missing columns: {missing_cols}")
            logger.info(f"  ✓ All expected columns present")

        if config['validation']['check_score_range']:
            if 'Score' in final_df.columns:
                scores = final_df['Score'].dropna()
                min_score = scores.min()
                max_score = scores.max()
                logger.info(f"  ✓ Score range: {min_score} - {max_score}")

                if min_score < 0 or max_score > 100:
                    logger.warning(f"  ⚠ Scores outside expected range [0-100]")

        # Validate reports
        if config['reporting']['check_html_reports']:
            logger.info("\nValidating HTML reports:")
            validate_reports(output_dir, logger)

        # Validate boost results
        if config['reporting']['check_boost_results']:
            logger.info("\nValidating boost analysis:")
            validate_boost_results(output_dir, logger)

        # Copy final results to test results directory for archiving
        logger.info("\nArchiving results to test directory:")
        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Save final results with timestamp
        archived_path = save_results(
            final_df,
            config['test_name'],
            results_dir,
            prefix="pipeline",
            suffix="final"
        )
        logger.info(f"  ✓ Archived final results: {archived_path}")

        # Also save a summary of pipeline structure
        summary = {
            'output_directory': output_dir,
            'final_results_file': final_results_path,
            'total_clusters': len(final_df),
            'score_range': f"{final_df['Score'].min():.1f} - {final_df['Score'].max():.1f}" if 'Score' in final_df.columns else 'N/A',
            'html_reports': len(glob.glob(os.path.join(output_dir, "02_reports", "*.html"))),
            'boost_clusters': len([d for d in os.listdir(os.path.join(output_dir, "03_boost_analysis")) if os.path.isdir(os.path.join(output_dir, "03_boost_analysis", d))]) if os.path.exists(os.path.join(output_dir, "03_boost_analysis")) else 0
        }

        summary_path = save_results(
            summary,
            config['test_name'],
            results_dir,
            prefix="pipeline",
            suffix="summary"
        )
        logger.info(f"  ✓ Saved pipeline summary: {summary_path}")

        # Log summary statistics
        logger.info("\n" + "=" * 80)
        logger.info("PIPELINE RESULTS SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Total clusters processed: {len(final_df)}")

        if 'Score' in final_df.columns:
            avg_score = final_df['Score'].mean()
            logger.info(f"Average quality score: {avg_score:.2f}")

            threshold = config['method_params'].get('score_threshold', 75)
            low_scores = final_df[final_df['Score'] < threshold]
            logger.info(f"Clusters below threshold ({threshold}): {len(low_scores)}")

        if 'Merged_Annotation_Broad' in final_df.columns:
            logger.info("✓ Merged annotations generated")

        # Show sample results
        logger.info("\nSample pipeline results:")
        display_cols = ['True Cell Type', 'Predicted Main Cell Type']
        if 'Score' in final_df.columns:
            display_cols.append('Score')
        if 'Merged_Annotation_Broad' in final_df.columns:
            display_cols.append('Merged_Annotation_Broad')

        for idx, row in final_df.head(3).iterrows():
            result_str = " | ".join([f"{col}: {row[col]}" for col in display_cols if col in row.index])
            logger.info(f"  {result_str}")

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ PIPELINE TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Pipeline output: {output_dir}")
        logger.info(f"📁 Test results: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ Pipeline test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1

    finally:
        # Ensure we're back in the original directory
        os.chdir(original_dir)


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
