"""
Complete CASSIA analysis pipeline.

This module provides the main pipeline function that orchestrates the entire
cell type annotation workflow including annotation, scoring, merging, and
report generation.
"""

import os
import datetime
import pandas as pd
import re


def extract_conversation_from_html(html_path):
    """
    Extract formatted conversation history from batch HTML report.

    The batch HTML report contains properly formatted annotation analysis
    with line breaks preserved. This function extracts that content for
    use in the final pipeline report.

    Args:
        html_path: Path to the batch HTML report

    Returns:
        Dictionary mapping cluster_id -> formatted conversation history HTML
    """
    if not os.path.exists(html_path):
        return {}

    with open(html_path, 'r', encoding='utf-8') as f:
        html_content = f.read()

    conversations = {}

    # Find all modal content sections using regex
    # Pattern: <div class="modal-content" id="modal-content-X">...</div>
    modal_pattern = r'<div class="modal-content" id="modal-content-(\d+)">(.*?)</div>\s*</div>\s*</div>\s*</div>'

    # Find cluster IDs from modal titles
    # Pattern: <h2 class="modal-title">CLUSTER_ID</h2>
    for match in re.finditer(r'<div class="modal-content" id="modal-content-\d+">', html_content):
        start_pos = match.start()
        # Find the end of this modal (look for the next modal or end of modals container)

        # Extract cluster ID from modal title
        title_match = re.search(r'<h2 class="modal-title">([^<]+)</h2>', html_content[start_pos:start_pos+2000])
        if not title_match:
            continue
        cluster_id = title_match.group(1).strip()

        # Extract annotation section content
        # Look for: <div class="modal-section annotation-section">...<div class="section-content">CONTENT</div>
        section_start = html_content.find('<div class="modal-section annotation-section">', start_pos)
        if section_start == -1 or section_start > start_pos + 50000:
            continue

        content_start = html_content.find('<div class="section-content">', section_start)
        if content_start == -1 or content_start > section_start + 5000:
            continue

        content_start += len('<div class="section-content">')

        # Find the closing </div> for section-content
        # Need to handle nested divs properly
        depth = 1
        pos = content_start
        while depth > 0 and pos < len(html_content):
            next_open = html_content.find('<div', pos)
            next_close = html_content.find('</div>', pos)

            if next_close == -1:
                break

            if next_open != -1 and next_open < next_close:
                depth += 1
                pos = next_open + 4
            else:
                depth -= 1
                if depth == 0:
                    content_end = next_close
                else:
                    pos = next_close + 6

        if depth == 0:
            annotation_content = html_content[content_start:content_end].strip()
            conversations[cluster_id] = annotation_content

    return conversations


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
    validate_api_keys_before_start: bool = False
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
        validate_api_keys_before_start (bool): If True, validates all required API keys before starting the pipeline.
            Fails fast with clear error messages if any keys are invalid. Default: False.
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

    # Validate API keys before starting pipeline (if requested)
    if validate_api_keys_before_start:
        try:
            from CASSIA.core.api_validation import validate_api_keys
        except ImportError:
            try:
                from ..core.api_validation import validate_api_keys
            except ImportError:
                from api_validation import validate_api_keys

        print("\n=== Validating API Keys ===")

        # Collect all providers that will be used
        providers_to_check = set()
        if annotation_provider:
            providers_to_check.add(annotation_provider)
        if score_provider:
            providers_to_check.add(score_provider)
        if annotationboost_provider:
            providers_to_check.add(annotationboost_provider)
        if merge_annotations and merge_provider:
            providers_to_check.add(merge_provider)

        # Remove None values
        providers_to_check.discard(None)

        # Validate each provider
        validation_failed = False
        for provider in providers_to_check:
            # Skip custom HTTP URLs (can't validate easily)
            if provider.startswith("http"):
                print(f"Skipping validation for custom provider: {provider}")
                continue

            is_valid = validate_api_keys(provider, verbose=True)
            if not is_valid:
                validation_failed = True

        if validation_failed:
            raise ValueError(
                "API key validation failed for one or more providers. "
                "Please check the error messages above and fix your API keys. "
                "You can set API keys with: CASSIA.set_api_key(provider, 'your-key')"
            )

        print("✓ All API keys validated successfully\n")

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
    raw_full_csv = os.path.join(csv_folder, f"{output_base_name}_full.csv")
    raw_summary_csv = os.path.join(csv_folder, f"{output_base_name}_summary.csv")
    raw_sorted_csv = os.path.join(csv_folder, f"{output_base_name}_sorted_full.csv")
    score_file_name = os.path.join(csv_folder, f"{output_base_name}_scored.csv")
    merged_annotation_file = os.path.join(csv_folder, f"{output_base_name}_merged.csv")

    # Reports go to the reports_folder - ALL HTML reports should be in this folder
    report_base_name = os.path.join(reports_folder, f"{output_base_name}")

    # First annotation output uses original output_file_name (may be absolute path from user)
    annotation_output = output_base_name

    print("\n=== Starting cell type analysis ===")
    # Run initial cell type analysis
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
        validator_involvement=validator_involvement
    )
    print("✓ Cell type analysis completed")

    # Copy the generated files to the organized folders
    original_full_csv = annotation_output + "_full.csv"
    original_summary_csv = annotation_output + "_summary.csv"

    # Copy the files if they exist
    if os.path.exists(original_full_csv):
        # Read and write instead of just copying to ensure compatibility
        df_full = pd.read_csv(original_full_csv)
        df_full.to_csv(raw_full_csv, index=False)
        print(f"Copied full results to {raw_full_csv}")
    if os.path.exists(original_summary_csv):
        df_summary = pd.read_csv(original_summary_csv)
        df_summary.to_csv(raw_summary_csv, index=False)
        print(f"Copied summary results to {raw_summary_csv}")

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
            df = pd.read_csv(raw_full_csv)
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
            print(f"! Error during annotation merging: {str(e)}")

    print("\n=== Starting scoring process ===")
    # Run scoring (generate_report=False because pipeline handles its own report)
    runCASSIA_score_batch(
        input_file=raw_full_csv,
        output_file=score_file_name,
        max_workers=max_workers,
        model=score_model,
        provider=score_provider,
        max_retries=max_retries,
        generate_report=False
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
        print(f"Warning: Could not create final combined results: {str(e)}")
        final_combined_file = score_file_name  # Fallback to scored file

    print("\n=== Generating main reports ===")
    # Read final combined CSV (includes merged groupings) and convert to list of dicts
    final_df = pd.read_csv(final_combined_file)
    rows_data = final_df.to_dict('records')

    # Extract formatted conversation history from batch HTML report (preserves line breaks)
    # The batch HTML report has properly formatted annotation analysis that we want to reuse
    batch_report = f"{annotation_output}_report.html"
    formatted_conversations = extract_conversation_from_html(batch_report)

    if formatted_conversations:
        print(f"Extracted formatted conversation history for {len(formatted_conversations)} clusters from batch report")
        # Replace conversation history with formatted version from HTML
        for row in rows_data:
            cluster_id = str(row.get('Cluster ID', ''))
            if cluster_id in formatted_conversations:
                # Store the pre-formatted HTML content
                row['_formatted_annotation_html'] = formatted_conversations[cluster_id]

    # Generate the HTML report (report_base_name already includes reports_folder path)
    report_output_path = f"{report_base_name}_report.html"
    generate_batch_html_report_from_data(
        rows=rows_data,
        output_path=report_output_path,
        report_title=f"CASSIA Pipeline Analysis - {tissue} ({species})"
    )
    print(f"✓ Generated report: {report_output_path}")

    # Clean up the batch HTML report (generated by runCASSIA_batch in current directory)
    # Now safe to delete since we've extracted the formatted content
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
                # NOTE: Using the raw_full_csv path to ensure the CSV can be found
                runCASSIA_annotationboost(
                    full_result_path=raw_full_csv,  # This is in the annotation_results_folder
                    marker=marker,
                    cluster_name=original_cluster_name,
                    major_cluster_info=major_cluster_info,
                    output_name=cluster_output_name,
                    num_iterations=5,
                    model=annotationboost_model,
                    provider=annotationboost_provider,
                    temperature=0,
                    conversation_history_mode=conversation_history_mode,
                    report_style=report_style
                )
            except IndexError:
                print(f"Error in pipeline: No data found for cluster: {original_cluster_name}")
            except Exception as e:
                print(f"Error in pipeline processing cluster {original_cluster_name}: {str(e)}")

        print("✓ Boost annotation completed")

    # Try to clean up the original files in the root directory
    try:
        for file_to_remove in [original_full_csv, original_summary_csv, annotation_output + "_sorted_full.csv"]:
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
