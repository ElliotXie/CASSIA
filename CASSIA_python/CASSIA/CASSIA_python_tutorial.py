#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CASSIA Analysis Tutorial

This Python script demonstrates a complete workflow using CASSIA for cell type annotation 
of single-cell RNA sequencing data. We'll analyze an intestinal cell dataset containing 
six distinct populations:

1. monocyte
2. plasma cells
3. cd8-positive, alpha-beta t cell
4. transit amplifying cell of large intestine
5. intestinal enteroendocrine cell
6. intestinal crypt stem cell
"""

# --------------------- Setup and Environment Preparation ---------------------

# Add current directory to path for proper imports
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Direct imports from local files, not from the installed package
from tools_function import *
from main_function_code import *
import pandas as pd
import numpy as np
import argparse
import re

# Import the new unified modules for annotation boost
try:
    from annotation_boost import (
        iterative_marker_analysis,
        runCASSIA_annotationboost,
        runCASSIA_annotationboost_additional_task
    )
    from llm_utils import call_llm
    print("Successfully imported unified annotation boost modules")
except ImportError as e:
    print(f"Note: Could not import unified modules: {str(e)}")
    print("Using original implementations from tools_function.py")
    # These will be provided by tools_function import

# Setup configuration variables
script_dir = os.path.dirname(os.path.abspath(__file__))
output_name = "intestine_detailed"
model_name = "google/gemini-2.5-flash-preview"
provider = "openrouter"
tissue = "large intestine"
species = "human"

# Load marker data (using relative file paths instead of the builtin loadmarker function)
def load_marker_data():
    """Load marker data from the CASSIA data directory, with column name compatibility handling."""
    processed_markers = pd.read_csv(os.path.join(script_dir, "data", "processed.csv"))
    unprocessed_markers = pd.read_csv(os.path.join(script_dir, "data", "unprocessed.csv"))
    subcluster_results = pd.read_csv(os.path.join(script_dir, "data", "subcluster_results.csv"))
    
    # Remove 'Unnamed: 0' column if it exists, as it's redundant with the 'gene' column
    for df in [processed_markers, unprocessed_markers, subcluster_results]:
        if 'Unnamed: 0' in df.columns:
            df.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Handle potential column name differences - ensure we have the required columns
    for df in [processed_markers, unprocessed_markers, subcluster_results]:
        # Check if avg_log2FC is missing but logFC is present (or similar variations)
        if 'avg_log2FC' not in df.columns:
            # Try alternative column names
            for alt_col in ['logFC', 'log2FC', 'Log2_fold_change', 'log2FoldChange']:
                if alt_col in df.columns:
                    df['avg_log2FC'] = df[alt_col]
                    break
        
        # Ensure all required columns exist, if not, add them with default values
        required_cols = ['avg_log2FC', 'pct.1', 'pct.2', 'p_val', 'p_val_adj'] 
        for col in required_cols:
            if col not in df.columns:
                df[col] = 0.0  # Add default value
    
    return processed_markers, unprocessed_markers, subcluster_results

# List available marker files
def list_available_markers():
    data_dir = os.path.join(script_dir, "data")
    return [f.replace('.csv', '') for f in os.listdir(data_dir) if f.endswith('.csv')]

# --------------------- Step 1: Full Pipeline ---------------------
def run_full_pipeline(marker_data):
    """
    This is the main function that runs the entire CASSIA pipeline in one go.
    If you just want to run the complete analysis without the modular approach,
    you can simply call this function with your marker data.
    
    Example:
        processed, unprocessed, subcluster = load_marker_data()
        run_full_pipeline(unprocessed)
    """
    print("\n=== Running Full CASSIA Pipeline ===")
    # Just call runCASSIA_pipeline directly if you want to run the whole pipeline at once
    
    # First, run the main pipeline without the annotation boost for low-scoring clusters
    runCASSIA_pipeline(
        output_file_name = "FastAnalysisResults",
        tissue = tissue,
        species = species,
        marker_path = marker_data,
        max_workers = 6,  # Matches the number of clusters in dataset
        annotation_model = model_name,
        annotation_provider = provider,
        score_model = model_name,
        score_provider = provider,
        score_threshold = 97,
        annotationboost_model = model_name,
        annotationboost_provider = provider,
        merge_annotations = True,  # Enable the built-in merging functionality
        merge_model = model_name   # Use the same model for merging
    )
    

# --------------------- Step 2: Batch Analysis Only ---------------------
def run_batch_analysis(marker_data):
    print("\n=== Running Batch Analysis Only ===")
    runCASSIA_batch(
        marker = marker_data,
        output_name = output_name,
        model = model_name,
        tissue = tissue,
        species = species,
        max_workers = 6,  # Matching cluster count
        n_genes = 50,
        additional_info = None,
        provider = provider
    )

# --------------------- Custom Merging Function ---------------------
def run_custom_merge(input_csv):
    """
    Run the merging process in a linear (non-parallel) way to avoid potential issues.
    
    This is an alternative implementation to the built-in merging in runCASSIA_pipeline.
    The key differences are:
    
    1. This implementation processes all clusters in a single LLM call instead of batching
    2. No parallelization is used, avoiding potential race conditions
    3. The results are sorted by True Cell Type before processing for consistency
    4. This uses direct access to the internal functions of the merging module
    
    Use this approach if you're experiencing issues with the default merging process
    in the runCASSIA_pipeline function.
    """
    from merging_annotation import call_llm, _create_annotation_prompt, _parse_llm_response
    
    print("\n=== Running Custom Linear Merging ===")
    
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Sort by True Cell Type to ensure consistent order
    df = df.sort_values(by=['True Cell Type'])
    
    # Setup columns
    cluster_col = 'True Cell Type'
    general_col = 'Predicted Main Cell Type'
    subtype_col = 'Predicted Sub Cell Types'
    
    # Process each row individually
    results = {}
    
    additional_context = f"These are cell clusters from {species} {tissue}."
    
    # Create one big prompt with all cells - avoids batch processing entirely
    prompt = _create_annotation_prompt(df, cluster_col, general_col, subtype_col, additional_context, detail_level="broad")
    
    # Call LLM once with all data
    print(f"Calling LLM to process {len(df)} clusters...")
    llm_response = call_llm(
        prompt=prompt,
        provider=provider,
        model=model_name,
        temperature=0
    )
    
    # Parse the response
    merged_annotations = _parse_llm_response(llm_response, df.index)
    
    # Add a new column for merged annotations
    df['Merged Annotation'] = df.index.map(lambda idx: merged_annotations.get(idx, "No annotation"))
    
    # Save the merged results
    output_path = input_csv.replace(".csv", "_merged.csv")
    df.to_csv(output_path, index=False)
    print(f"Merged annotations saved to {output_path}")
    
    return df

# --------------------- Step 3: Quality Scoring ---------------------
def run_quality_scoring(input_csv, output_csv=None):
    print("\n=== Running Quality Scoring ===")
    if output_csv is None:
        output_csv = input_csv.replace(".csv", "_scored.csv")
        
    runCASSIA_score_batch(
        input_file = input_csv,
        output_file = output_csv,
        max_workers = 6,
        model = model_name,
        provider = provider
    )
    return output_csv

# --------------------- Step 4: Generate Report ---------------------
def generate_report(scored_csv, report_name=None):
    print("\n=== Generating Report ===")
    if report_name is None:
        report_name = scored_csv.replace("_scored.csv", "_report")
        
    runCASSIA_generate_score_report(
        csv_path = scored_csv,
        index_name = report_name
    )

# --------------------- Step 5: Uncertainty Quantification ---------------------
def run_uncertainty_quantification(marker_data):
    print("\n=== Running Uncertainty Quantification ===")
    # Run multiple iterations
    iteration_results = runCASSIA_batch_n_times(
        n=2,
        marker=marker_data,
        output_name=output_name + "_Uncertainty",
        model=model_name,
        provider=provider,
        tissue=tissue,
        species=species,
        max_workers=6,
        batch_max_workers=3  # Conservative setting for API rate limits
    )

    # Calculate similarity scores
    similarity_scores = runCASSIA_similarity_score_batch(
        marker=marker_data,
        file_pattern=output_name + "_Uncertainty_*_full.csv",
        output_name="intestine_uncertainty",
        max_workers=6,
        model=model_name,
        provider=provider,
        main_weight=0.5,
        sub_weight=0.5
    )

# --------------------- Step 6: Annotation Boost ---------------------
def run_annotation_boost(marker_data, full_csv, cluster_name="monocyte", provider_test=None, debug_mode=False, test_genes=None, conversation_history_mode="final"):
    """
    Run annotation boost for a specific cluster.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Name of the cluster to analyze (default: "monocyte")
        provider_test: Optional provider to test (default: uses global provider)
        debug_mode: Enable debug mode for diagnostics
        test_genes: List of test genes to check in the marker data
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    """
    print(f"\n=== Running Annotation Boost for {cluster_name} ===")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # Import debug utilities if in debug mode
    if debug_mode:
        try:
            from CASSIA.debug_genes import run_gene_diagnostics
            print("Successfully imported debug utilities")
            
            # Define the test genes if not specified
            if test_genes is None:
                test_genes = ["CD133", "CD9", "ChAT", "DCLK1", "EDNRB", "ERBB3", "FABP7", "GFAP", "KIT", "LGR5", "NGFR", "NKX2-2", "NOS1", "OLIG2", "PGP9.5", "PROM1", "RET", "S100B", "SOX9", "UCHL1", "VIP"]
            
            # Generate a test conversation
            test_conversation = f"""
            Based on the marker genes, I would like to check some additional genes to confirm this cell type:
            <check_genes>{', '.join(test_genes[:10])}</check_genes>
            
            Let's also check these additional markers:
            <check_genes>{', '.join(test_genes[10:])}</check_genes>
            """
            
            # Run diagnostics
            print("\n=== Running Gene Extraction Diagnostics ===")
            print(f"Testing with marker data: {marker_data.shape}")
            
            try:
                # Try normal import first
                from debug_genes import run_gene_diagnostics
            except ImportError:
                try:
                    # Try relative import as fallback
                    from CASSIA.debug_genes import run_gene_diagnostics
                except ImportError:
                    raise ImportError("Could not import debug_genes module. Make sure it's in the correct directory.")
            
            # Run full diagnostics
            run_gene_diagnostics(marker_data, test_conversation, test_genes)
            
        except ImportError as e:
            print(f"Could not import debug utilities: {e}")
    
    # Make sure the CSV file exists
    if not os.path.exists(full_csv):
        print(f"Error: File not found: {full_csv}")
        print(f"Current working directory: {os.getcwd()}")
        print("Please provide the correct path to the CSV file")
        return
    
    # Read the CSV to ensure the cluster name matches exactly what's in the file
    try:
        df = pd.read_csv(full_csv)
        print(f"Successfully loaded {full_csv} with {len(df)} rows")
    except Exception as e:
        print(f"Error reading CSV file: {str(e)}")
        return
    
    # Check if the True Cell Type column exists
    if 'True Cell Type' not in df.columns:
        print(f"Error: 'True Cell Type' column not found in {full_csv}")
        print(f"Available columns: {df.columns.tolist()}")
        return
    
    # Check if the cluster exists in the dataframe
    if cluster_name not in df['True Cell Type'].values:
        # Try to find the closest match (case insensitive)
        matches = df[df['True Cell Type'].str.lower() == cluster_name.lower()]
        if not matches.empty:
            # Use the exact case/format from the file
            cluster_name = matches.iloc[0]['True Cell Type']
            print(f"Using exact cluster name from file: '{cluster_name}'")
        else:
            print(f"Warning: Cluster '{cluster_name}' not found in {full_csv}")
            print(f"Available clusters: {df['True Cell Type'].tolist()}")
            return
    
    # Create a sanitized version for the output filename
    output_filename = f"{cluster_name.replace(',', '')}_annotationboost"
        
    try:
        # Now run with the exact cluster name
        print(f"Running annotation boost with {test_provider} provider")
        print(f"Using conversation history mode: {conversation_history_mode}")
        result = runCASSIA_annotationboost(
            full_result_path = full_csv,
            marker = marker_data,
            output_name = output_filename,  # Remove comma for filename only
            cluster_name = cluster_name,  # Use original format for data lookup
            major_cluster_info = f"{species.title()} {tissue.title()}",
            num_iterations = 5,
            model = model_name,
            provider = test_provider,
            conversation_history_mode = conversation_history_mode
        )
        
        # Check if the result is successful
        if isinstance(result, dict):
            status = result.get('status', 'unknown')
            if status == 'success':
                print(f"\n✅ Successfully completed annotation boost for {cluster_name}")
                print(f"Results saved to:")
                for key in ['formatted_report_path', 'raw_report_path']:
                    if key in result:
                        print(f"  - {key}: {result[key]}")
                print(f"Execution time: {result.get('execution_time', 0):.2f} seconds")
            elif status in ['error', 'partial_error', 'critical_error']:
                print(f"\n❌ Error in annotation boost: {result.get('error_message', 'Unknown error')}")
                # Try to report any partial results if available
                for key in ['formatted_report_path', 'raw_report_path']:
                    if key in result and result[key]:
                        print(f"  - Partial {key}: {result[key]}")
            else:
                print(f"\n⚠️ Unknown result status: {status}")
        else:
            # Handle old style return value (for backward compatibility)
            print(f"Successfully completed annotation boost for {cluster_name}")
            print(f"Results saved with prefix: {output_filename}")
        
        # Offer to open the reports if they were generated
        if isinstance(result, dict) and result.get('status') == 'success':
            if 'formatted_report_path' in result and os.path.exists(result['formatted_report_path']):
                print(f"\nTo view the formatted report, open: {result['formatted_report_path']}")
            if 'raw_report_path' in result and result['raw_report_path'] and os.path.exists(result['raw_report_path']):
                print(f"To view the raw conversation report, open: {result['raw_report_path']}")
        
        print(f"Successfully completed annotation boost for {cluster_name}")
    except Exception as e:
        print(f"Error in run_annotation_boost: {str(e)}")
        import traceback
        traceback.print_exc()
        print("Available clusters in the CSV file:")
        print(df['True Cell Type'].tolist())

# --------------------- Step 7: Compare Celltypes ---------------------
def run_celltype_comparison():
    print("\n=== Running Cell Type Comparison ===")
    # The marker here are copied from CASSIA's previous results
    marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

    compareCelltypes(
        tissue = tissue,
        celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
        marker_set = marker,
        species = species,
        output_file = "plasama_cell_subtype"
    )

# --------------------- Step 8: Subclustering ---------------------
def run_subclustering(subcluster_data):
    print("\n=== Running Subclustering Analysis ===")
    runCASSIA_subclusters(
        marker = subcluster_data,
        major_cluster_info = "cd8 t cell",
        output_name = "subclustering_results",
        model = model_name,
        provider = provider
    )

    # Run multiple iterations for subclustering
    runCASSIA_n_subcluster(
        n=5, 
        marker=subcluster_data,
        major_cluster_info="cd8 t cell", 
        base_output_name="subclustering_results_n",
        model=model_name,
        temperature=0,
        provider=provider,
        max_workers=5,
        n_genes=50
    )

    # Calculate similarity scores for subclusters
    runCASSIA_similarity_score_batch(
        marker = subcluster_data,
        file_pattern = "subclustering_results_n_*.csv",
        output_name = "subclustering_uncertainty",
        max_workers = 6,
        model = model_name,
        provider = provider,
        main_weight = 0.5,
        sub_weight = 0.5
    )

# --------------------- Step 9: Annotation Boost with Additional Task ---------------------
def run_annotation_boost_with_task(marker_data, full_csv, cluster_name=None, additional_task=None, provider_test=None, conversation_history_mode="final"):
    """
    Run annotation boost with an additional task for a cluster.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Optional name of the cluster to analyze (default: "cd8-positive, alpha-beta t cell")
        additional_task: Optional task description (default: infer cell state and function)
        provider_test: Optional provider to test (default: uses global provider)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    """
    # Default cluster to analyze if not specified
    if cluster_name is None:
        cluster_name = "cd8-positive, alpha-beta t cell"
        
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
        
    print(f"\n=== Running Annotation Boost with Additional Task for {cluster_name} ===")
    
    # Make sure the CSV file exists
    if not os.path.exists(full_csv):
        print(f"Error: File not found: {full_csv}")
        print(f"Current working directory: {os.getcwd()}")
        print("Please provide the correct path to the CSV file")
        return
    
    # Read the CSV to ensure the cluster name matches exactly what's in the file
    try:
        df = pd.read_csv(full_csv)
        print(f"Successfully loaded {full_csv} with {len(df)} rows")
    except Exception as e:
        print(f"Error reading CSV file: {str(e)}")
        return
    
    # Check if the True Cell Type column exists
    if 'True Cell Type' not in df.columns:
        print(f"Error: 'True Cell Type' column not found in {full_csv}")
        print(f"Available columns: {df.columns.tolist()}")
        return
    
    # Check if the cluster exists in the dataframe
    if cluster_name not in df['True Cell Type'].values:
        # Try to find the closest match (case insensitive)
        matches = df[df['True Cell Type'].str.lower() == cluster_name.lower()]
        if not matches.empty:
            # Use the exact case/format from the file
            cluster_name = matches.iloc[0]['True Cell Type']
            print(f"Using exact cluster name from file: '{cluster_name}'")
        else:
            print(f"Warning: Cluster '{cluster_name}' not found in {full_csv}")
            print(f"Available clusters: {df['True Cell Type'].tolist()}")
            return
    
    # Create a sanitized version for the output filename
    # Replace commas, spaces and other problematic characters
    safe_cluster = re.sub(r'[,\s+]', '_', cluster_name)
    output_filename = f"{output_name}_{safe_cluster}_boosted"
    
    try:
        # Define the additional task if not specified
        if additional_task is None:
            additional_task = "infer the state and function of this cell cluster, and determine if it shows signs of exhaustion or activation"
        
        print(f"Additional task: {additional_task}")
        print(f"Running annotation boost with {test_provider} provider")
        print(f"Using conversation history mode: {conversation_history_mode}")
        
        # Call the annotation boost function with the exact cluster name from the file
        result = runCASSIA_annotationboost_additional_task(
            full_result_path = full_csv,
            marker = marker_data,
            output_name = output_filename,
            cluster_name = cluster_name,  # Use original cluster name with comma
            major_cluster_info = f"{species.title()} {tissue.title()}",
            num_iterations = 5,
            model = model_name,
            provider = test_provider,
            additional_task = additional_task,
            conversation_history_mode = conversation_history_mode
        )
        
        # Check if the result is successful
        if isinstance(result, dict):
            status = result.get('status', 'unknown')
            if status == 'success':
                print(f"\n✅ Successfully completed annotation boost for {cluster_name}")
                print(f"Results saved to:")
                for key in ['summary_report_path', 'raw_report_path', 'formatted_report_path']:
                    if key in result:
                        print(f"  - {key}: {result[key]}")
                print(f"Execution time: {result.get('execution_time', 0):.2f} seconds")
            elif status in ['error', 'partial_error', 'critical_error']:
                print(f"\n❌ Error in annotation boost: {result.get('error_message', 'Unknown error')}")
                # Try to report any partial results if available
                for key in ['summary_report_path', 'raw_report_path', 'formatted_report_path']:
                    if key in result and result[key]:
                        print(f"  - Partial {key}: {result[key]}")
            else:
                print(f"\n⚠️ Unknown result status: {status}")
        else:
            # Handle old style return value (for backward compatibility)
            print(f"Successfully completed annotation boost for {cluster_name}")
            print(f"Results saved with prefix: {output_filename}")
        
        # Offer to open the reports if they were generated
        if isinstance(result, dict) and result.get('status') == 'success':
            if 'formatted_report_path' in result and os.path.exists(result['formatted_report_path']):
                print(f"\nTo view the formatted report, open: {result['formatted_report_path']}")
            if 'raw_report_path' in result and result['raw_report_path'] and os.path.exists(result['raw_report_path']):
                print(f"To view the raw conversation report, open: {result['raw_report_path']}")
    
    except Exception as e:
        print(f"Error in run_annotation_boost_with_task: {str(e)}")
        import traceback
        traceback.print_exc()
        print("\nAvailable clusters in the CSV file:")
        print(df['True Cell Type'].tolist())

# --------------------- New: Test All Annotation Boost Providers ---------------------
def test_annotation_boost_providers(marker_data, full_csv, cluster_name="monocyte", conversation_history_mode="final"):
    """
    Test the annotation boost functionality with different providers.
    This function helps validate that the unified annotation boost implementation 
    works correctly with all supported providers.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Name of the cluster to analyze (default: "monocyte")
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    """
    print("\n=== Testing Annotation Boost with Multiple Providers ===")
    print(f"Using conversation history mode: {conversation_history_mode}")
    
    # Test with OpenAI provider
    print("\n----- Testing with OpenAI provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="openai", conversation_history_mode=conversation_history_mode)
    except Exception as e:
        print(f"Error with OpenAI provider: {str(e)}")
    
    # Test with Anthropic provider
    print("\n----- Testing with Anthropic provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="anthropic", conversation_history_mode=conversation_history_mode)
    except Exception as e:
        print(f"Error with Anthropic provider: {str(e)}")
    
    # Test with OpenRouter provider
    print("\n----- Testing with OpenRouter provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="openrouter", conversation_history_mode=conversation_history_mode)
    except Exception as e:
        print(f"Error with OpenRouter provider: {str(e)}")
    
    # Test additional task functionality with OpenRouter
    print("\n----- Testing Annotation Boost with Additional Task -----")
    try:
        run_annotation_boost_with_task(
            marker_data, 
            full_csv, 
            cluster_name, 
            additional_task="check if this cell type expresses cancer markers", 
            provider_test="openrouter",
            conversation_history_mode=conversation_history_mode
        )
    except Exception as e:
        print(f"Error with additional task: {str(e)}")

def setup_api_keys():
    """Setup API keys for various providers from environment variables."""
    import os
    
    # Check if API keys are already set in environment
    api_key_openai = os.environ.get('OPENAI_API_KEY', '')
    api_key_anthropic = os.environ.get('ANTHROPIC_API_KEY', '')
    api_key_openrouter = os.environ.get('OPENROUTER_API_KEY', '')
    
    # Only prompt if keys are not already set
    if not api_key_openai and not api_key_anthropic and not api_key_openrouter:
        print("No API keys found in environment variables.")
        print("CASSIA requires at least one API key to function.")
        print("You can set these in your environment or enter them now.")
        
        # Prompt for OpenAI key if not set
        if not api_key_openai:
            key = input("Enter your OpenAI API key (or press Enter to skip): ")
            if key.strip():
                os.environ['OPENAI_API_KEY'] = key
                set_openai_api_key(key)
        
        # Prompt for Anthropic key if not set
        if not api_key_anthropic:
            key = input("Enter your Anthropic API key (or press Enter to skip): ")
            if key.strip():
                os.environ['ANTHROPIC_API_KEY'] = key
                set_anthropic_api_key(key)
        
        # Prompt for OpenRouter key if not set
        if not api_key_openrouter:
            key = input("Enter your OpenRouter API key (or press Enter to skip): ")
            if key.strip():
                os.environ['OPENROUTER_API_KEY'] = key
                set_openrouter_api_key(key)
    else:
        # Set the keys in CASSIA even if they're already in the environment
        if api_key_openai:
            set_openai_api_key(api_key_openai)
        if api_key_anthropic:
            set_anthropic_api_key(api_key_anthropic)
        if api_key_openrouter:
            set_openrouter_api_key(api_key_openrouter)
            
    print("API key setup complete")

def main():
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description='Run CASSIA analysis pipelines')
    parser.add_argument('--step', type=str, default='all',
                      help='Which step to run: all, batch, merge, score, report, uncertainty, boost, compare, subcluster, boost_task, test_boost, debug_genes')
    parser.add_argument('--input_csv', type=str, default=None,
                      help='Input CSV file for steps that require it (merge, score, report, boost)')
    parser.add_argument('--cluster', type=str, default='monocyte',
                      help='Cluster name for annotation boost')
    parser.add_argument('--task', type=str, default=None,
                      help='Additional task for annotation boost with task, e.g., "check if this is a cancer cell"')
    parser.add_argument('--provider', type=str, default=None,
                      help='Provider to use for API calls (openai, anthropic, openrouter)')
    parser.add_argument('--test_genes', type=str, default=None,
                      help='Comma-separated list of genes to test for the debug_genes step')
    parser.add_argument('--history_mode', type=str, default="final",
                      help='Conversation history mode for annotation boost: "full", "final", or "none"')
    args = parser.parse_args()
    
    # Override default provider if specified in command line
    global provider
    if args.provider:
        provider = args.provider
        print(f"Using provider specified in command line: {provider}")
    
    # Setup API keys first
    setup_api_keys()
    
    # Load marker data
    processed, unprocessed, subcluster = load_marker_data()
    print(f"Successfully loaded marker data with {unprocessed.shape[0]} genes and {unprocessed.shape[1]} columns")
    
    # Print available markers
    if args.step == 'list_markers':
        print("Available markers:", list_available_markers())
        return
    
    # Get the default input CSV if not provided
    input_csv = args.input_csv or r"C:\Users\ellio\OneDrive - UW-Madison\CASSIA+\CASSIA_large_intestine_human_20250513_225204\TEST2_full.csv"
    
    # Run the selected step
    if args.step == 'all':
        # This will use the built-in merging process in runCASSIA_pipeline
        try:
            run_full_pipeline(unprocessed)
        except Exception as e:
            print(f"Error in pipeline: {str(e)}")
    elif args.step == 'batch':
        run_batch_analysis(unprocessed)
    elif args.step == 'merge':
        # This uses our custom non-parallel merging implementation
        # If you're having issues with the built-in merging, use this approach instead
        run_custom_merge(input_csv)
    elif args.step == 'score':
        run_quality_scoring(input_csv)
    elif args.step == 'report':
        input_csv = args.input_csv or f"{output_name}_scored.csv"
        generate_report(input_csv)
    elif args.step == 'uncertainty':
        run_uncertainty_quantification(unprocessed)
    elif args.step == 'boost':
        run_annotation_boost(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode)
    elif args.step == 'compare':
        run_celltype_comparison()
    elif args.step == 'subcluster':
        run_subclustering(subcluster)
    elif args.step == 'boost_task':
        try:
            run_annotation_boost_with_task(unprocessed, input_csv, args.cluster, args.task, conversation_history_mode=args.history_mode)
        except Exception as e:
            print(f"Error in annotation boost with task: {str(e)}")
    elif args.step == 'test_boost':
        # New option to test the unified annotation boost functionality
        try:
            test_annotation_boost_providers(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode)
        except Exception as e:
            print(f"Error testing annotation boost: {str(e)}")
    elif args.step == 'debug_genes':
        # New option to run gene extraction diagnostics
        try:
            # Parse test genes if provided
            test_genes = None
            if args.test_genes:
                test_genes = [g.strip() for g in args.test_genes.split(',')]
            
            print(f"=== Running Gene Extraction Diagnostics ===")
            print(f"Testing with marker data: {unprocessed.shape}")
            
            try:
                # Try normal import first
                from debug_genes import run_gene_diagnostics
            except ImportError:
                try:
                    # Try relative import as fallback
                    from CASSIA.debug_genes import run_gene_diagnostics
                except ImportError:
                    raise ImportError("Could not import debug_genes module. Make sure it's in the correct directory.")
            
            # Create a sample conversation with the problematic genes
            if test_genes is None:
                test_genes = ["CD133", "CD9", "ChAT", "DCLK1", "EDNRB", "ERBB3", "FABP7", "GFAP", "KIT", "LGR5", "NGFR", "NKX2-2", "NOS1", "OLIG2", "PGP9.5", "PROM1", "RET", "S100B", "SOX9", "UCHL1", "VIP"]
            
            test_conversation = f"""
            Based on the marker genes, I would like to check some additional genes to confirm this cell type:
            <check_genes>{', '.join(test_genes[:10])}</check_genes>
            
            Let's also check these additional markers:
            <check_genes>{', '.join(test_genes[10:])}</check_genes>
            """
            
            # Run full diagnostics
            run_gene_diagnostics(unprocessed, test_conversation, test_genes)
            
            # Test with a specific cluster
            if args.cluster:
                run_annotation_boost(unprocessed, input_csv, args.cluster, debug_mode=True, test_genes=test_genes)
            
        except Exception as e:
            print(f"Error in gene extraction diagnostics: {str(e)}")
            import traceback
            traceback.print_exc()
    else:
        print(f"Unknown step: {args.step}")
        print("Available steps: all, batch, merge, score, report, uncertainty, boost, compare, subcluster, boost_task, test_boost, debug_genes")


if __name__ == "__main__":
    # Import necessary modules
    import os
    import sys
    import pandas as pd
    import time
    import csv
    import argparse
    from concurrent.futures import ThreadPoolExecutor, as_completed
    
    # Global configuration
    species = "human"
    tissue = "large intestine"
    model_name = "google/gemini-2.5-flash-preview"  # Using specified model
    output_name = f"CASSIA_{tissue.replace(' ', '_')}_{species}"

    
    # Run the main function
    main() 


    # Example commands to test each case in a Windows terminal (Command Prompt or PowerShell):
    # python CASSIA_python_tutorial.py --step all
    # python CASSIA_python_tutorial.py --step batch
    # python CASSIA_python_tutorial.py --step merge
    # python CASSIA_python_tutorial.py --step score
    # python CASSIA_python_tutorial.py --step report
    # python CASSIA_python_tutorial.py --step uncertainty
    # python CASSIA_python_tutorial.py --step boost
    # python CASSIA_python_tutorial.py --step compare
    # python CASSIA_python_tutorial.py --step subcluster
    # python CASSIA_python_tutorial.py --step boost_task
    # python CASSIA_python_tutorial.py --step test_boost
    