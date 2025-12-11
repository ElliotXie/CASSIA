"""
Complete CASSIA analysis pipeline.

This module provides the main pipeline function that orchestrates the entire
cell type annotation workflow including annotation, scoring, merging, and
report generation.
"""

import os
import json
import datetime
import pandas as pd


def runCASSIA_pipeline(
    output_file_name: str,
    tissue: str,
    species: str,
    marker,  # Can be DataFrame or file path string
    max_workers: int = 4,
    annotation_model: str = "meta-llama/llama-4-maverick",
    annotation_provider: str = "openrouter",
    score_model: str = "google/gemini-2.5-pro-preview-03-25",
    score_provider: str = "openrouter",
    annotationboost_model: str = "google/gemini-2.5-flash-preview",
    annotationboost_provider: str = "openrouter",
    score_threshold: float = 75,
    additional_info: str = "None",
    max_retries: int = 1,
    merge_annotations: bool = True,
    merge_model: str = "deepseek/deepseek-chat-v3-0324",
    merge_provider: str = "openrouter",
    conversation_history_mode: str = "final",
    ranking_method: str = "avg_log2FC",
    ascending: bool = None,
    report_style: str = "per_iteration",
    validator_involvement: str = "v1",
    output_dir: str = None,
    validate_api_keys_before_start: bool = True,
    auto_convert_ids: bool = True
):
    """
    Run the complete cell analysis pipeline including annotation, scoring, and report generation.

    Args:
        output_file_name (str): Base name for output files
        tissue (str): Tissue type being analyzed
        species (str): Species being analyzed
        marker: Marker data (pandas DataFrame or path to CSV file)
        max_workers (int): Maximum number of concurrent workers
        annotation_model (str): Model to use for initial annotation
        annotation_provider (str): Provider for initial annotation
        score_model (str): Model to use for scoring
        score_provider (str): Provider for scoring
        annotationboost_model (str): Model to use for boosting low-scoring annotations
        annotationboost_provider (str): Provider for boosting low-scoring annotations
        score_threshold (float): Threshold for identifying low-scoring clusters
        additional_info (str): Additional information for analysis
        max_retries (int): Maximum number of retries for failed analyses
        merge_annotations (bool): Whether to merge annotations from LLM
        merge_model (str): Model to use for merging annotations
        merge_provider (str): Provider to use for merging annotations
        conversation_history_mode (str): Mode for extracting conversation history ("full", "final", or "none")
        ranking_method (str): Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending (bool): Sort direction (None uses default for each method)
        report_style (str): Style of report generation ("per_iteration" or "total_summary")
        validator_involvement (str): Validator involvement level
        output_dir (str): Directory where the output folder will be created. If None, uses current working directory.
        validate_api_keys_before_start (bool): If True, validates the API key before starting.
            Fails fast with clear error message if the key is invalid. Default: True.
        auto_convert_ids (bool): Automatically convert Ensembl/Entrez gene IDs to gene symbols.
            If True (default), detects and converts IDs in the marker data before processing.
            Requires the 'mygene' package to be installed for conversion.
    """
    # Import dependencies here to avoid circular imports
    try:
        from CASSIA.engine.tools_function import runCASSIA_batch
        from CASSIA.evaluation.scoring import runCASSIA_score_batch
        from CASSIA.reports.generate_batch_report import generate_batch_html_report_from_data
        from CASSIA.agents.annotation_boost.annotation_boost import runCASSIA_annotationboost
    except ImportError:
        try:
            from ..engine.tools_function import runCASSIA_batch
            from ..evaluation.scoring import runCASSIA_score_batch
            from ..reports.generate_batch_report import generate_batch_html_report_from_data
            from ..agents.annotation_boost.annotation_boost import runCASSIA_annotationboost
        except ImportError:
            from tools_function import runCASSIA_batch
            from scoring import runCASSIA_score_batch
            from generate_batch_report import generate_batch_html_report_from_data
            from annotation_boost import runCASSIA_annotationboost

    # Import validation function
    try:
        from CASSIA.core.validation import validate_runCASSIA_pipeline_inputs
    except ImportError:
        try:
            from ..core.validation import validate_runCASSIA_pipeline_inputs
        except ImportError:
            from validation import validate_runCASSIA_pipeline_inputs

    # Validate all inputs early (fail-fast)
    validate_runCASSIA_pipeline_inputs(
        output_file_name=output_file_name,
        tissue=tissue,
        species=species,
        marker=marker,
        max_workers=max_workers,
        max_retries=max_retries,
        score_threshold=score_threshold,
        conversation_history_mode=conversation_history_mode,
        report_style=report_style
    )

    # Note: API key validation is handled by runCASSIA_batch (called first)
    # via the validate_api_key_before_start parameter

    # Determine base directory for output
    if output_dir is not None:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
        base_dir = output_dir
    else:
        base_dir = "."  # Current working directory (backward compatible)

    # Create a main folder based on tissue and species for organizing reports
    main_folder_name = f"CASSIA_Pipeline_{tissue}_{species}"
    main_folder_name = "".join(c for c in main_folder_name if c.isalnum() or c in (' ', '-', '_')).strip()
    main_folder_name = main_folder_name.replace(' ', '_')

    # Remove .csv extension if present
    if output_file_name.lower().endswith('.csv'):
        output_file_name = output_file_name[:-4]  # Remove last 4 characters (.csv)

    # Extract just the filename (in case an absolute path was provided)
    # This ensures internal folder paths work correctly
    output_base_name = os.path.basename(output_file_name)

    # Add timestamp to prevent overwriting existing folders with the same name
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    main_folder_name = f"{main_folder_name}_{timestamp}"

    # Create the main folder inside base_dir
    main_folder_path = os.path.join(base_dir, main_folder_name)
    if not os.path.exists(main_folder_path):
        os.makedirs(main_folder_path)
        print(f"Created main folder: {main_folder_path}")

    # Create organized subfolders according to user's specifications
    reports_folder = os.path.join(main_folder_path, "01_annotation_report")  # HTML report
    boost_folder = os.path.join(main_folder_path, "02_annotation_boost")  # All annotation boost related results
    csv_folder = os.path.join(main_folder_path, "03_csv_files")  # All CSV files

    # Create all subfolders
    for folder in [reports_folder, boost_folder, csv_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"Created subfolder: {folder}")

    # Define derived file names with folder paths (use output_base_name for internal paths)
    # All CSV files go to the csv_folder
    raw_summary_csv = os.path.join(csv_folder, f"{output_base_name}_summary.csv")
    raw_conversations_json = os.path.join(csv_folder, f"{output_base_name}_conversations.json")
    raw_sorted_csv = os.path.join(csv_folder, f"{output_base_name}_sorted_summary.csv")
    score_file_name = os.path.join(csv_folder, f"{output_base_name}_scored.csv")
    merged_annotation_file = os.path.join(csv_folder, f"{output_base_name}_merged.csv")

    # Reports go to the reports_folder - ALL HTML reports should be in this folder
    report_base_name = os.path.join(reports_folder, f"{output_base_name}")

    # First annotation output uses original output_file_name (may be absolute path from user)
    annotation_output = output_base_name

    print("\n=== Starting cell type analysis ===")
    # Run initial cell type analysis
    # API key validation happens here via validate_api_key_before_start
    runCASSIA_batch(
        marker=marker,
        output_name=annotation_output,
        model=annotation_model,
        tissue=tissue,
        species=species,
        additional_info=additional_info,
        provider=annotation_provider,
        max_workers=max_workers,
        max_retries=max_retries,
        ranking_method=ranking_method,
        ascending=ascending,
        validator_involvement=validator_involvement,
        validate_api_key_before_start=validate_api_keys_before_start,
        auto_convert_ids=auto_convert_ids
    )
    print("✓ Cell type analysis completed")

    # Copy the generated files to the organized folders
    original_summary_csv = annotation_output + "_summary.csv"
    original_conversations_json = annotation_output + "_conversations.json"

    # Copy the files if they exist
    if os.path.exists(original_summary_csv):
        df_summary = pd.read_csv(original_summary_csv)
        df_summary.to_csv(raw_summary_csv, index=False)
        print(f"Copied summary results to {raw_summary_csv}")
    if os.path.exists(original_conversations_json):
        import shutil
        shutil.copy2(original_conversations_json, raw_conversations_json)
        print(f"Copied conversations JSON to {raw_conversations_json}")

    # Verify batch output was created before proceeding
    if not os.path.exists(raw_summary_csv):
        raise RuntimeError(
            f"Batch annotation failed to create output file: {raw_summary_csv}\n"
            "Cannot proceed with scoring."
        )

    # Check if file has actual data (not just headers)
    df_check = pd.read_csv(raw_summary_csv)
    if len(df_check) == 0:
        raise RuntimeError(
            f"Batch annotation produced empty results in: {raw_summary_csv}\n"
            "Cannot proceed with scoring."
        )
    print(f"Verified batch output: {len(df_check)} clusters annotated")

    # Merge annotations if requested
    if merge_annotations:
        print("\n=== Starting annotation merging ===")

        # Import the merge_annotations function dynamically
        try:
            try:
                from CASSIA.agents.merging.merging_annotation import merge_annotations_all
            except ImportError:
                try:
                    from ..agents.merging.merging_annotation import merge_annotations_all
                except ImportError:
                    from merging_annotation import merge_annotations_all

            # Sort the CSV file by Cluster ID before merging to ensure consistent order
            print("Sorting CSV by Cluster ID before merging...")
            df = pd.read_csv(raw_summary_csv)
            df = df.sort_values(by=['Cluster ID'])
            df.to_csv(raw_sorted_csv, index=False)

            # Run the merging process on the sorted CSV
            merge_annotations_all(
                csv_path=raw_sorted_csv,
                output_path=merged_annotation_file,
                provider=merge_provider,
                model=merge_model,
                additional_context=f"These are cell clusters from {species} {tissue}. {additional_info}"
            )
            print(f"✓ Annotations merged and saved to {merged_annotation_file}")
        except Exception as e:
            print(f"Warning: Annotation merging failed: {str(e)}")
            print("Continuing without merged annotations (this is an optional feature)...")

    print("\n=== Starting scoring process ===")
    # Run scoring (generate_report=False because pipeline handles its own report)
    runCASSIA_score_batch(
        input_file=raw_summary_csv,
        output_file=score_file_name,
        max_workers=max_workers,
        model=score_model,
        provider=score_provider,
        max_retries=max_retries,
        generate_report=False,
        conversations_json_path=raw_conversations_json
    )
    print("✓ Scoring process completed")

    print("\n=== Creating final combined results ===")
    # Create final combined CSV with all results
    try:
        # Read the scored file (which has all the original data plus scores)
        final_df = pd.read_csv(score_file_name)

        # If merged annotations exist, add merged columns
        if os.path.exists(merged_annotation_file):
            merged_df = pd.read_csv(merged_annotation_file)
            # Merge on 'Cluster ID' to add merged annotation columns
            if 'Cluster ID' in merged_df.columns:
                # Keep only the merged columns (not duplicating existing ones)
                merge_columns = [col for col in merged_df.columns if col not in final_df.columns or col == 'Cluster ID']
                final_df = final_df.merge(merged_df[merge_columns], on='Cluster ID', how='left')

        # Sort the final results by Cluster ID
        final_df = final_df.sort_values(by=['Cluster ID'])

        # Save the final combined results
        final_combined_file = os.path.join(csv_folder, f"{output_base_name}_FINAL_RESULTS.csv")
        final_df.to_csv(final_combined_file, index=False)
        print(f"✓ Final combined results saved to {final_combined_file}")

    except Exception as e:
        raise RuntimeError(f"Failed to create final combined results: {str(e)}")

    print("\n=== Generating main reports ===")
    # Read final combined CSV (includes merged groupings) and convert to list of dicts
    final_df = pd.read_csv(final_combined_file)
    rows_data = final_df.to_dict('records')

    # Load conversation history from JSON and add to rows_data
    if os.path.exists(raw_conversations_json):
        try:
            with open(raw_conversations_json, 'r', encoding='utf-8') as f:
                conversations_data = json.load(f)
            print(f"Loaded conversation history for {len(conversations_data)} clusters from JSON")

            # Add conversation history to each row (as structured dict for HTML generation)
            for row in rows_data:
                cluster_id = str(row.get('Cluster ID', ''))
                if cluster_id in conversations_data:
                    # Pass structured dict directly - generate_batch_report handles it
                    row['Conversation History'] = conversations_data[cluster_id]
        except Exception as e:
            print(f"Warning: Could not load conversations JSON: {e}")

    # Generate the HTML report (report_base_name already includes reports_folder path)
    report_output_path = f"{report_base_name}_report.html"
    generate_batch_html_report_from_data(
        rows=rows_data,
        output_path=report_output_path,
        report_title=f"CASSIA Pipeline Analysis - {tissue} ({species})"
    )
    print(f"✓ Generated report: {report_output_path}")

    # Clean up the batch HTML report (generated by runCASSIA_batch in current directory)
    batch_report = f"{annotation_output}_report.html"
    if os.path.exists(batch_report):
        try:
            os.remove(batch_report)
            print(f"Cleaned up redundant batch report: {batch_report}")
        except Exception as e:
            print(f"Warning: Could not remove batch report: {e}")

    print("✓ Main reports generated")

    print("\n=== Analyzing low-scoring clusters ===")
    # Handle low-scoring clusters
    df = pd.read_csv(score_file_name)
    low_score_clusters = df[df['Score'] < score_threshold]['Cluster ID'].tolist()

    print(f"Found {len(low_score_clusters)} clusters with scores below {score_threshold}:")
    print(low_score_clusters)

    if low_score_clusters:
        print("\n=== Starting boost annotation for low-scoring clusters ===")

        # Create boosted reports list - we will NOT generate a combined report
        for cluster in low_score_clusters:
            print(f"Processing low score cluster: {cluster}")

            # Keep the original cluster name for data lookup
            original_cluster_name = cluster

            # Sanitize the cluster name only for file naming purposes
            sanitized_cluster_name = "".join(c for c in str(cluster) if c.isalnum() or c in (' ', '-', '_')).strip()

            # Create individual folder for this cluster's boost analysis
            cluster_boost_folder = os.path.join(boost_folder, sanitized_cluster_name)
            if not os.path.exists(cluster_boost_folder):
                os.makedirs(cluster_boost_folder)

            # Define output name for the cluster boost report
            cluster_output_name = os.path.join(cluster_boost_folder, f"{output_base_name}_{sanitized_cluster_name}_boosted")

            # Use the original name for data lookup
            try:
                # major_cluster_info should be simple user-provided information like "human large intestine"
                # NOT complex data extracted from CSV
                major_cluster_info = f"{species} {tissue}"

                # Run annotation boost - use original cluster name for data lookup, but sanitized name for output file
                runCASSIA_annotationboost(
                    full_result_path=raw_summary_csv,
                    marker=marker,
                    cluster_name=original_cluster_name,
                    major_cluster_info=major_cluster_info,
                    output_name=cluster_output_name,
                    num_iterations=5,
                    model=annotationboost_model,
                    provider=annotationboost_provider,
                    temperature=0,
                    conversation_history_mode=conversation_history_mode,
                    report_style=report_style,
                    conversations_json_path=raw_conversations_json,
                    species=species,
                    auto_convert_ids=auto_convert_ids
                )
            except IndexError:
                print(f"Error in pipeline: No data found for cluster: {original_cluster_name}")
            except Exception as e:
                print(f"Error in pipeline processing cluster {original_cluster_name}: {str(e)}")

        print("✓ Boost annotation completed")

    # Try to clean up the original files in the root directory
    try:
        for file_to_remove in [original_summary_csv, original_conversations_json, annotation_output + "_report.html"]:
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)
                print(f"Removed original file: {file_to_remove}")
    except Exception as e:
        print(f"Warning: Could not remove some temporary files: {str(e)}")

    print("\n=== Cell type analysis pipeline completed ===")
    print(f"All results have been organized in the '{main_folder_path}' folder:")
    print(f"  📊 MAIN RESULTS: {final_combined_file}")
    print(f"  📁 HTML Report: {reports_folder}")
    print(f"  🔍 Annotation Boost Results: {boost_folder}")
    print(f"  📂 CSV Files: {csv_folder}")
    print(f"\n✅ Your final results are in: {os.path.basename(final_combined_file)}")
