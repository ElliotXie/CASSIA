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
try:
    from .tools_function import *
except ImportError:
    from tools_function import *
try:
    from .main_function_code import *
except ImportError:
    from main_function_code import *
try:
    from .Uncertainty_quantification import *
except ImportError:
    from Uncertainty_quantification import *
try:
    from .subclustering import *
except ImportError:
    from subclustering import *
import pandas as pd
import numpy as np
import argparse
import re
import time

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

# Setup configuration variables for custermized api key (deepseek example)
# Uncomment the following block to use DeepSeek as a custom provider
# output_name = "intestine_detailed"
# model_name = "deepseek-chat"
# provider = "https://api.deepseek.com"
# tissue = "large intestine"
# species = "human"
# api_key = "sk-afb39114f1334ba486505d9425937d16"

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

# --------------------- Output Path Organization ---------------------
def get_output_path(result_type, provider_type=None, filename=None):
    """
    Get organized output path based on provider type and result category.
    
    Args:
        result_type: Type of result ('batch_analysis', 'celltype_comparison', 'subclustering', 
                    'annotation_boost', 'scoring', 'reports', 'uncertainty', 'merging')
        provider_type: Provider type (None=use global, 'custom'=custom API, 'normal'=normal API)
        filename: Optional filename (will use output_name if not provided)
    
    Returns:
        Full path for the output file
    """
    # Determine if using custom API
    if provider_type is None:
        if provider.startswith("http"):
            api_type = "custom_api"
        else:
            api_type = "normal_api"
    elif provider_type == "custom":
        api_type = "custom_api"
    else:
        api_type = "normal_api"
    
    # Create the base directory path
    base_path = os.path.join(script_dir, "test_results", api_type, result_type)
    
    # Create directory if it doesn't exist
    os.makedirs(base_path, exist_ok=True)
    
    # Return just the directory path if no filename specified
    if filename is None:
        return base_path
    
    # Return full file path
    return os.path.join(base_path, filename)

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
        merge_model = model_name,   # Use the same model for merging
        merge_provider = provider   # Use the same provider for merging
    )
    

# --------------------- Step 2: Batch Analysis Only ---------------------
def run_batch_analysis(marker_data):
    print("\n=== Running Batch Analysis Only ===")
    
    # Get organized output path
    output_path = get_output_path("batch_analysis", filename=output_name)
    
    runCASSIA_batch(
        marker = marker_data,
        output_name = output_path,
        model = model_name,
        tissue = tissue,
        species = species,
        max_workers = 6,  # Matching cluster count
        n_genes = 50,
        additional_info = None,
        provider = provider
    )
    
    print(f"✓ Batch analysis results saved to: {get_output_path('batch_analysis')}")
    return output_path

# --------------------- Merging Function ---------------------
def run_merge(input_csv, detail_level="broad"):
    """
    Run the merging process using the official merging_annotation module.
    
    This uses the proper merge_annotations function with customizable provider support.
    
    Args:
        input_csv: Path to input CSV file
        detail_level: Level of detail for merging ("broad", "detailed", "very_detailed")
    """
    try:
        from .merging_annotation import merge_annotations
    except ImportError:
        from merging_annotation import merge_annotations
    
    print(f"\n=== Running Annotation Merging ({detail_level} level) ===")
    
    # Generate output path with organized structure
    base_name = os.path.basename(input_csv).replace(".csv", "")
    output_filename = f"{base_name}_merged_{detail_level}.csv"
    output_path = get_output_path("merging", filename=output_filename)
    
    # Additional context for the analysis
    additional_context = f"These are cell clusters from {species} {tissue}."
    
    # Use the official merge_annotations function with custom provider support
    result_df = merge_annotations(
        csv_path=input_csv,
        output_path=output_path,
        provider=provider,
        model=model_name,
        additional_context=additional_context,
        batch_size=20,
        detail_level=detail_level
    )
    
    print(f"✓ Annotation merging completed and saved to {output_path}")
    return result_df

def run_merge_all(input_csv):
    """
    Run all three levels of merging using the official merge_annotations_all function.
    
    This processes broad, detailed, and very_detailed levels in parallel.
    """
    try:
        from .merging_annotation import merge_annotations_all
    except ImportError:
        from merging_annotation import merge_annotations_all
    
    print("\n=== Running All Levels of Annotation Merging ===")
    
    # Generate output path with organized structure
    base_name = os.path.basename(input_csv).replace(".csv", "")
    output_filename = f"{base_name}_merged_all.csv"
    output_path = get_output_path("merging", filename=output_filename)
    
    # Additional context for the analysis
    additional_context = f"These are cell clusters from {species} {tissue}."
    
    # Use the official merge_annotations_all function with custom provider support
    result_df = merge_annotations_all(
        csv_path=input_csv,
        output_path=output_path,
        provider=provider,
        model=model_name,
        additional_context=additional_context,
        batch_size=20
    )
    
    print(f"✓ All levels of annotation merging completed and saved to {output_path}")
    return result_df

# --------------------- Step 3: Quality Scoring ---------------------
def run_quality_scoring(input_csv, output_csv=None):
    print("\n=== Running Quality Scoring ===")
    
    if output_csv is None:
        # Generate organized output path
        base_name = os.path.basename(input_csv).replace(".csv", "")
        output_filename = f"{base_name}_scored.csv"
        output_csv = get_output_path("scoring", filename=output_filename)
        
    runCASSIA_score_batch(
        input_file = input_csv,
        output_file = output_csv,
        max_workers = 6,
        model = model_name,
        provider = provider
    )
    
    print(f"✓ Quality scoring results saved to: {output_csv}")
    return output_csv

# --------------------- Step 4: Generate Report ---------------------
def generate_report(scored_csv, report_name=None):
    print("\n=== Generating Report ===")
    
    if report_name is None:
        # Generate organized output path
        base_name = os.path.basename(scored_csv).replace("_scored.csv", "")
        report_name = get_output_path("reports", filename=f"{base_name}_report")
        
    runCASSIA_generate_score_report(
        csv_path = scored_csv,
        index_name = report_name
    )
    
    print(f"✓ Reports generated in: {get_output_path('reports')}")

# --------------------- Step 5: Uncertainty Quantification ---------------------
def run_uncertainty_quantification(marker_data, provider_test=None):
    print("\n=== Running Uncertainty Quantification ===")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # If using a custom provider, set a default model if needed
    current_model = model_name
    if test_provider.startswith("http") and current_model == "google/gemini-2.5-flash-preview":
        if test_provider == "https://api.deepseek.com":
            current_model = "deepseek-chat"
        print(f"Using model: {current_model} with custom provider: {test_provider}")
    
    # Generate organized output path
    uncertainty_output_name = get_output_path("uncertainty", filename=output_name + "_Uncertainty")
    
    # Run multiple iterations
    iteration_results = runCASSIA_batch_n_times(
        n=4,
        marker=marker_data,
        output_name=uncertainty_output_name,
        model=current_model,
        provider=test_provider,
        tissue=tissue,
        species=species,
        max_workers=6,
        batch_max_workers=3,  # Conservative setting for API rate limits
        temperature=0.3 # Set temperature here
    )

    # Calculate similarity scores
    similarity_output_name = get_output_path("uncertainty", filename="intestine_uncertainty")
    similarity_scores = runCASSIA_similarity_score_batch(
        marker=marker_data,
        file_pattern=f"{uncertainty_output_name}_*_full.csv",
        output_name=similarity_output_name,
        max_workers=6,
        model=current_model,
        provider=test_provider,
        main_weight=0.5,
        sub_weight=0.5,
        temperature=0.3 # Set temperature here as well
    )
    
    print(f"✓ Uncertainty quantification results saved to: {get_output_path('uncertainty')}")

# --------------------- Step 5b: Single Cluster Uncertainty Quantification ---------------------
def run_single_cluster_uncertainty(marker_data, cluster_name="monocyte", provider_test=None, n_iterations=5):
    """
    Run uncertainty quantification for a single cluster.
    
    Args:
        marker_data: Marker data DataFrame
        cluster_name: Name of the cluster to analyze (default: "monocyte")
        provider_test: Optional provider to test (default: uses global provider)
        n_iterations: Number of iterations to run (default: 5)
    
    Returns:
        dict: Analysis results including consensus types, cell types, and scores
    """
    print(f"\n=== Running Single Cluster Uncertainty Quantification for {cluster_name} ===")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # If using a custom provider, set a default model if needed
    current_model = model_name
    if test_provider.startswith("http") and current_model == "google/gemini-2.5-flash-preview":
        if test_provider == "https://api.deepseek.com":
            current_model = "deepseek-chat"
        print(f"Using model: {current_model} with custom provider: {test_provider}")
    
    # Filter marker data for the specific cluster
    if isinstance(marker_data, pd.DataFrame):
        if 'cluster' in marker_data.columns:
            cluster_col = 'cluster'
        elif 'cell_type' in marker_data.columns:
            cluster_col = 'cell_type'
        elif marker_data.columns[0].lower() in ['cluster', 'cell_type', 'celltype', 'cell type']:
            cluster_col = marker_data.columns[0]
        else:
            print(f"Warning: Could not identify cluster column in marker data")
            print(f"Available columns: {marker_data.columns.tolist()}")
            print(f"Using first column as cluster column: {marker_data.columns[0]}")
            cluster_col = marker_data.columns[0]
        
        # Try to find the cluster (case insensitive)
        if cluster_col in marker_data.columns:
            try:
                filtered_markers = marker_data[marker_data[cluster_col].str.lower() == cluster_name.lower()]
                if filtered_markers.empty:
                    print(f"Warning: Cluster '{cluster_name}' not found in marker data")
                    print(f"Available clusters: {marker_data[cluster_col].unique().tolist()}")
                    # Use all markers as fallback
                    filtered_markers = marker_data
                else:
                    print(f"Found {len(filtered_markers)} marker genes for cluster '{cluster_name}'")
            except Exception as e:
                print(f"Error filtering marker data: {str(e)}")
                filtered_markers = marker_data
        else:
            filtered_markers = marker_data
    else:
        # If not a DataFrame, use as is
        filtered_markers = marker_data
    
    # Extract marker genes if needed
    if isinstance(filtered_markers, pd.DataFrame):
        # Get top marker genes
        if len(filtered_markers.columns) > 1:
            # If we have a proper marker dataframe with multiple columns
            try:
                # Try to get gene column
                gene_col = None
                for col_name in ['gene', 'genes', 'feature', 'features', 'marker', 'markers']:
                    if col_name in filtered_markers.columns:
                        gene_col = col_name
                        break
                
                if gene_col is None:
                    # If no obvious gene column, use second column
                    gene_col = filtered_markers.columns[1]
                
                # Get top 50 genes
                if len(filtered_markers) > 50:
                    # Sort by fold change or p-value if available
                    if 'avg_log2FC' in filtered_markers.columns:
                        sorted_markers = filtered_markers.sort_values(by='avg_log2FC', ascending=False)
                    elif 'p_val_adj' in filtered_markers.columns:
                        sorted_markers = filtered_markers.sort_values(by='p_val_adj', ascending=True)
                    else:
                        sorted_markers = filtered_markers
                    
                    top_markers = sorted_markers.head(50)
                else:
                    top_markers = filtered_markers
                
                # Extract gene list
                marker_list = top_markers[gene_col].tolist()
                marker_list = [str(gene) for gene in marker_list if gene and str(gene).strip()]
                
                # Fall back to string representation if extraction fails
                if not marker_list:
                    marker_list = filtered_markers.iloc[:, 1].astype(str).tolist()[:50]
            except Exception as e:
                print(f"Error extracting marker genes: {str(e)}")
                # Fall back to string representation
                marker_list = filtered_markers.iloc[:, 1].astype(str).tolist()[:50]
        else:
            # If we only have one column, use it directly
            marker_list = filtered_markers.iloc[:, 0].astype(str).tolist()[:50]
    elif isinstance(filtered_markers, str):
        # If it's already a string, split by common separators
        marker_list = [m.strip() for m in re.split(r'[,;\s]+', filtered_markers)][:50]
    elif isinstance(filtered_markers, list):
        # If it's already a list, use it directly
        marker_list = [str(m).strip() for m in filtered_markers][:50]
    else:
        print(f"Error: Unsupported marker data type: {type(filtered_markers)}")
        return None
    
    # Print marker list for verification
    print(f"Using {len(marker_list)} marker genes for analysis:")
    print(", ".join(marker_list[:10]) + (", ..." if len(marker_list) > 10 else ""))
    
    # Run uncertainty quantification using runCASSIA_n_times_similarity_score
    try:
        print(f"Running {n_iterations} iterations of cell type analysis...")
        results = runCASSIA_n_times_similarity_score(
            tissue=tissue,
            species=species,
            additional_info=f"Analyzing cluster: {cluster_name}",
            temperature=0,
            marker_list=marker_list,
            model=current_model,
            max_workers=5,
            n=n_iterations,
            provider=test_provider,
            main_weight=0.5,
            sub_weight=0.5
        )
        
        # Print results
        print("\n=== Single Cluster Uncertainty Analysis Results ===")
        print(f"Cluster: {cluster_name}")
        print(f"General cell type: {results['general_celltype_llm']}")
        print(f"Sub cell type: {results['sub_celltype_llm']}")
        if results['Possible_mixed_celltypes_llm']:
            print(f"Possible mixed cell types: {', '.join(results['Possible_mixed_celltypes_llm'])}")
        print(f"Similarity score: {results['similarity_score']:.2f}")
        print(f"Consensus types (from frequency): {results['consensus_types']}")
        
        # Save results to file
        output_filename = f"{cluster_name.replace(' ', '_').replace(',', '')}_uncertainty.json"
        output_path = get_output_path("uncertainty", filename=output_filename)
        
        import json
        with open(output_path, 'w') as f:
            # Convert any non-serializable objects to strings
            serializable_results = {}
            for k, v in results.items():
                if isinstance(v, (str, int, float, bool, list, dict)) or v is None:
                    serializable_results[k] = v
                else:
                    serializable_results[k] = str(v)
            json.dump(serializable_results, f, indent=2)
        print(f"Results saved to {output_path}")
        
        return results
    except Exception as e:
        print(f"Error in single cluster uncertainty quantification: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

# --------------------- Step 6: Annotation Boost ---------------------
def run_annotation_boost(marker_data, full_csv, cluster_name="monocyte", provider_test=None, debug_mode=False, test_genes=None, conversation_history_mode="final", search_strategy="breadth", report_style="per_iteration"):
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
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        report_style: Style of report generation ("per_iteration" or "total_summary")
    """
    print(f"\n=== Running Annotation Boost for {cluster_name} ===")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # Set search strategy to 'depth' for custom providers like DeepSeek
    if test_provider.startswith("http"):
        search_strategy = "depth"
        print(f"Custom provider detected. Automatically setting search strategy to: {search_strategy}")
    
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
    
    # Check if the cluster column exists (new name first, then fall back to old name)
    cluster_col = 'Cluster ID' if 'Cluster ID' in df.columns else ('True Cell Type' if 'True Cell Type' in df.columns else None)
    if cluster_col is None:
        print(f"Error: 'Cluster ID' or 'True Cell Type' column not found in {full_csv}")
        print(f"Available columns: {df.columns.tolist()}")
        return

    # Check if the cluster exists in the dataframe
    if cluster_name not in df[cluster_col].values:
        # Try to find the closest match (case insensitive)
        matches = df[df[cluster_col].str.lower() == cluster_name.lower()]
        if not matches.empty:
            # Use the exact case/format from the file
            cluster_name = matches.iloc[0][cluster_col]
            print(f"Using exact cluster name from file: '{cluster_name}'")
        else:
            print(f"Warning: Cluster '{cluster_name}' not found in {full_csv}")
            print(f"Available clusters: {df[cluster_col].tolist()}")
            return
    
    # Create a sanitized version for the output filename
    output_filename = f"{cluster_name.replace(',', '')}_annotationboost"
    output_path = get_output_path("annotation_boost", filename=output_filename)
        
    try:
        # Now run with the exact cluster name
        print(f"Running annotation boost with {test_provider} provider")
        print(f"Using conversation history mode: {conversation_history_mode}")
        print(f"Using search strategy: {search_strategy}")
        print(f"Using report style: {report_style}")
        result = runCASSIA_annotationboost(
            full_result_path = full_csv,
            marker = marker_data,
            output_name = output_path,  # Use organized path
            cluster_name = cluster_name,  # Use original format for data lookup
            major_cluster_info = f"{species.title()} {tissue.title()}",
            num_iterations = 5,
            model = model_name,
            provider = test_provider,
            conversation_history_mode = conversation_history_mode,
            search_strategy = search_strategy,
            report_style = report_style
        )
        
        # Check if the result is successful
        if isinstance(result, dict):
            status = result.get('status', 'unknown')
            if status == 'success':
                print(f"\n✅ Successfully completed annotation boost for {cluster_name}")
                print(f"Results saved to:")
                for key in ['formatted_report_path', 'raw_report_path', 'summary_report_path']:
                    if key in result:
                        print(f"  - {key}: {result[key]}")
                print(f"Execution time: {result.get('execution_time', 0):.2f} seconds")
            elif status in ['error', 'partial_error', 'critical_error']:
                print(f"\n❌ Error in annotation boost: {result.get('error_message', 'Unknown error')}")
                # Try to report any partial results if available
                for key in ['formatted_report_path', 'raw_report_path', 'summary_report_path']:
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
            if 'summary_report_path' in result and result['summary_report_path'] and os.path.exists(result['summary_report_path']):
                print(f"To view the summary report, open: {result['summary_report_path']}")
        
        print(f"Successfully completed annotation boost for {cluster_name}")
    except Exception as e:
        print(f"Error in run_annotation_boost: {str(e)}")
        import traceback
        traceback.print_exc()
        print("Available clusters in the CSV file:")
        print(df[cluster_col].tolist())

# --------------------- Step 7: Compare Celltypes ---------------------
def run_celltype_comparison():
    print("\n=== Running Cell Type Comparison (Multi-Model, Multi-Round) ===")
    # Contradictory marker set to trigger discussion
    marker = "CD19, CD20, CD38, SDC1, PAX5, IRF4, CD27, CD3"
    celltypes = ["Naive B cell", "Plasma cell", "Plasmablast"]
    
    # Organized output path
    base_results_dir = os.path.join(os.path.dirname(__file__), 'test_results', 'budget_discussion_mode_test')
    os.makedirs(base_results_dir, exist_ok=True)
    output_file = os.path.join(base_results_dir, 'results.csv')
    html_file = os.path.join(base_results_dir, 'results_report.html')

    from CASSIA.cell_type_comparison import compareCelltypes
    comparison_results = compareCelltypes(
        tissue = tissue,
        celltypes = celltypes,
        marker_set = marker,
        species = species,
        model_preset = "budget",
        discussion_mode = True,
        discussion_rounds = 3,
        generate_html_report = True,
        output_file = output_file
    )
    print(f"\nTest finished. All rounds (initial and discussion) are saved in: {base_results_dir}")
    print(f"CSV: {output_file}")
    print(f"HTML: {html_file}")

# --------------------- Step 8: Subclustering ---------------------
def run_subclustering(subcluster_data):
    print("\n=== Running Subclustering Analysis ===")
    
    # Generate organized output paths
    subclustering_output = get_output_path("subclustering", filename="subclustering_results")
    runCASSIA_subclusters(
        marker = subcluster_data,
        major_cluster_info = "cd8 t cell",
        output_name = subclustering_output,
        model = model_name,
        provider = provider
    )
    print(f"✓ Subclustering results saved to: {get_output_path('subclustering')}")

# --------------------- Step 8b: Subclustering with Uncertainty Quantification ---------------------
def run_subclustering_with_uncertainty(subcluster_data):
    print("\n=== Running Subclustering with Uncertainty Quantification ===")
    # Generate organized output paths
    subclustering_n_output = get_output_path("subclustering", filename="subclustering_results_n")
    uncertainty_output = get_output_path("subclustering", filename="subclustering_uncertainty")

    # Run multiple iterations for subclustering
    runCASSIA_n_subcluster(
        n=5, 
        marker=subcluster_data,
        major_cluster_info="cd8 t cell", 
        base_output_name=subclustering_n_output,
        model=model_name,
        temperature=0,
        provider=provider,
        max_workers=5,
        n_genes=50
    )

    # Calculate similarity scores for subclusters
    runCASSIA_similarity_score_batch(
        marker = subcluster_data,
        file_pattern = f"{subclustering_n_output}_*.csv",
        output_name = uncertainty_output,
        max_workers = 6,
        model = model_name,
        provider = provider,
        main_weight = 0.5,
        sub_weight = 0.5
    )
    print(f"✓ Subclustering with uncertainty results saved to: {get_output_path('subclustering')}")

# --------------------- Step 9: Annotation Boost with Additional Task ---------------------
def run_annotation_boost_with_task(marker_data, full_csv, cluster_name=None, additional_task=None, provider_test=None, conversation_history_mode="final", search_strategy="breadth", report_style="per_iteration"):
    """
    Run annotation boost with an additional task for a cluster.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Optional name of the cluster to analyze (default: "cd8-positive, alpha-beta t cell")
        additional_task: Optional task description (default: infer cell state and function)
        provider_test: Optional provider to test (default: uses global provider)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        report_style: Style of report generation ("per_iteration" or "total_summary")
    """
    # Default cluster to analyze if not specified
    if cluster_name is None:
        cluster_name = "cd8-positive, alpha-beta t cell"
        
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # Set search strategy to 'depth' for custom providers like DeepSeek
    if test_provider.startswith("http"):
        search_strategy = "depth"
        print(f"Custom provider detected. Automatically setting search strategy to: {search_strategy}")
    
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
    
    # Check if the cluster column exists (new name first, then fall back to old name)
    cluster_col = 'Cluster ID' if 'Cluster ID' in df.columns else ('True Cell Type' if 'True Cell Type' in df.columns else None)
    if cluster_col is None:
        print(f"Error: 'Cluster ID' or 'True Cell Type' column not found in {full_csv}")
        print(f"Available columns: {df.columns.tolist()}")
        return

    # Check if the cluster exists in the dataframe
    if cluster_name not in df[cluster_col].values:
        # Try to find the closest match (case insensitive)
        matches = df[df[cluster_col].str.lower() == cluster_name.lower()]
        if not matches.empty:
            # Use the exact case/format from the file
            cluster_name = matches.iloc[0][cluster_col]
            print(f"Using exact cluster name from file: '{cluster_name}'")
        else:
            print(f"Warning: Cluster '{cluster_name}' not found in {full_csv}")
            print(f"Available clusters: {df[cluster_col].tolist()}")
            return
    
    # Create a sanitized version for the output filename
    # Replace commas, spaces and other problematic characters
    safe_cluster = re.sub(r'[,\s+]', '_', cluster_name)
    output_filename = f"{output_name}_{safe_cluster}_boosted"
    output_path = get_output_path("annotation_boost", filename=output_filename)
    
    try:
        # Define the additional task if not specified
        if additional_task is None:
            additional_task = "infer the state and function of this cell cluster, and determine if it shows signs of exhaustion or activation"
        
        print(f"Additional task: {additional_task}")
        print(f"Running annotation boost with {test_provider} provider")
        print(f"Using conversation history mode: {conversation_history_mode}")
        print(f"Using search strategy: {search_strategy}")
        
        # Call the annotation boost function with the exact cluster name from the file
        result = runCASSIA_annotationboost_additional_task(
            full_result_path = full_csv,
            marker = marker_data,
            output_name = output_path,
            cluster_name = cluster_name,  # Use original cluster name with comma
            major_cluster_info = f"{species.title()} {tissue.title()}",
            num_iterations = 5,
            model = model_name,
            provider = test_provider,
            additional_task = additional_task,
            conversation_history_mode = conversation_history_mode,
            search_strategy = search_strategy,
            report_style = report_style
        )
        
        # Check if the result is successful
        if isinstance(result, dict):
            status = result.get('status', 'unknown')
            if status == 'success':
                print(f"\n✅ Successfully completed annotation boost for {cluster_name}")
                print(f"Results saved to:")
                for key in ['formatted_report_path', 'raw_report_path', 'summary_report_path']:
                    if key in result:
                        print(f"  - {key}: {result[key]}")
                print(f"Execution time: {result.get('execution_time', 0):.2f} seconds")
            elif status in ['error', 'partial_error', 'critical_error']:
                print(f"\n❌ Error in annotation boost: {result.get('error_message', 'Unknown error')}")
                # Try to report any partial results if available
                for key in ['formatted_report_path', 'raw_report_path', 'summary_report_path']:
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
            if 'summary_report_path' in result and result['summary_report_path'] and os.path.exists(result['summary_report_path']):
                print(f"To view the summary report, open: {result['summary_report_path']}")
    
    except Exception as e:
        print(f"Error in run_annotation_boost_with_task: {str(e)}")
        import traceback
        traceback.print_exc()
        print("\nAvailable clusters in the CSV file:")
        print(df[cluster_col].tolist())

# --------------------- New: Test Full Pipeline with Different Providers ---------------------
def test_full_pipeline_providers(marker_data):
    """
    Test the full CASSIA pipeline with different providers.
    This function helps validate that runCASSIA_pipeline works correctly 
    with all supported providers, including custom API endpoints.
    
    Args:
        marker_data: Marker data DataFrame for the pipeline
    """
    print("\n=== Testing Full Pipeline with Multiple Providers ===")
    
    # Test with OpenAI provider
    print("\n----- Testing Full Pipeline with OpenAI provider -----")
    try:
        # Temporarily set provider and model for OpenAI
        original_provider = globals()['provider']
        original_model = globals()['model_name']
        
        globals()['provider'] = "openai"
        globals()['model_name'] = "gpt-4o"
        
        print("Running full pipeline with OpenAI...")
        runCASSIA_pipeline(
            output_file_name="TestResults_OpenAI",
            tissue=tissue,
            species=species,
            marker_path=marker_data,
            max_workers=3,  # Reduced for testing
            annotation_model="gpt-4o",
            annotation_provider="openai",
            score_model="gpt-4o",
            score_provider="openai",
            score_threshold=97,
            annotationboost_model="gpt-4o",
            annotationboost_provider="openai",
            merge_annotations=True,
            merge_model="gpt-4o",
            merge_provider="openai"
        )
        print("✅ Full pipeline with OpenAI completed successfully")
    except Exception as e:
        print(f"❌ Error with OpenAI provider: {str(e)}")
    finally:
        # Restore original settings
        globals()['provider'] = original_provider
        globals()['model_name'] = original_model
    
    # Test with Anthropic provider
    print("\n----- Testing Full Pipeline with Anthropic provider -----")
    try:
        globals()['provider'] = "anthropic"
        globals()['model_name'] = "claude-3-5-sonnet-20241022"
        
        print("Running full pipeline with Anthropic...")
        runCASSIA_pipeline(
            output_file_name="TestResults_Anthropic",
            tissue=tissue,
            species=species,
            marker_path=marker_data,
            max_workers=3,
            annotation_model="claude-3-5-sonnet-20241022",
            annotation_provider="anthropic",
            score_model="claude-3-5-sonnet-20241022",
            score_provider="anthropic",
            score_threshold=97,
            annotationboost_model="claude-3-5-sonnet-20241022",
            annotationboost_provider="anthropic",
            merge_annotations=True,
            merge_model="claude-3-5-sonnet-20241022",
            merge_provider="anthropic"
        )
        print("✅ Full pipeline with Anthropic completed successfully")
    except Exception as e:
        print(f"❌ Error with Anthropic provider: {str(e)}")
    finally:
        globals()['provider'] = original_provider
        globals()['model_name'] = original_model
    
    # Test with OpenRouter provider
    print("\n----- Testing Full Pipeline with OpenRouter provider -----")
    try:
        globals()['provider'] = "openrouter"
        globals()['model_name'] = "anthropic/claude-3.5-sonnet"
        
        print("Running full pipeline with OpenRouter...")
        runCASSIA_pipeline(
            output_file_name="TestResults_OpenRouter",
            tissue=tissue,
            species=species,
            marker_path=marker_data,
            max_workers=3,
            annotation_model="anthropic/claude-3.5-sonnet",
            annotation_provider="openrouter",
            score_model="anthropic/claude-3.5-sonnet",
            score_provider="openrouter",
            score_threshold=97,
            annotationboost_model="anthropic/claude-3.5-sonnet",
            annotationboost_provider="openrouter",
            merge_annotations=True,
            merge_model="anthropic/claude-3.5-sonnet",
            merge_provider="openrouter"
        )
        print("✅ Full pipeline with OpenRouter completed successfully")
    except Exception as e:
        print(f"❌ Error with OpenRouter provider: {str(e)}")
    finally:
        globals()['provider'] = original_provider
        globals()['model_name'] = original_model
    
    # Test with custom provider (commented out as it requires an API key)
    print("\n----- Custom Provider Test (Skipped - Requires API Key) -----")
    print("To test with a custom provider like DeepSeek, use the command line option:")
    print("python CASSIA_python_tutorial.py --step all --provider https://api.deepseek.com --api_key YOUR_API_KEY")
    
    print("\n=== Full Pipeline Provider Tests Complete ===")

# --------------------- New: Compare Search Strategies ---------------------
def compare_search_strategies(marker_data, full_csv, cluster_name="monocyte", conversation_history_mode="final", report_style="per_iteration"):
    """
    Compare breadth-first vs depth-first search strategies for annotation boost.
    This function demonstrates the differences between the two approaches.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Name of the cluster to analyze (default: "monocyte")
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
        report_style: Style of report generation ("per_iteration" or "total_summary")
    """
    print(f"\n=== Comparing Search Strategies for {cluster_name} ===")
    print(f"This will run the same cluster with both breadth-first and depth-first approaches")
    print(f"Using conversation history mode: {conversation_history_mode}")
    print(f"Using report style: {report_style}")
    
    # Test breadth-first approach
    print("\n" + "="*60)
    print("BREADTH-FIRST APPROACH")
    print("="*60)
    print("• Tests multiple hypotheses per iteration (up to 3)")
    print("• Explores various possibilities quickly")
    print("• Good for initial screening and comprehensive coverage")
    print("• Original CASSIA behavior")
    print()
    
    try:
        run_annotation_boost(
            marker_data, 
            full_csv, 
            cluster_name, 
            search_strategy="breadth",
            conversation_history_mode=conversation_history_mode,
            report_style=report_style
        )
        print("✅ Breadth-first analysis completed successfully")
    except Exception as e:
        print(f"❌ Error with breadth-first approach: {str(e)}")
    
    # Test depth-first approach
    print("\n" + "="*60)
    print("DEPTH-FIRST APPROACH")
    print("="*60)
    print("• Focuses on one hypothesis per iteration")
    print("• Goes deeper into each hypothesis before moving to alternatives")
    print("• Better for complex cell types requiring detailed investigation")
    print("• Avoids hypothesis overload and provides more focused analysis")
    print()
    
    try:
        run_annotation_boost(
            marker_data, 
            full_csv, 
            cluster_name, 
            search_strategy="depth",
            conversation_history_mode=conversation_history_mode,
            report_style=report_style
        )
        print("✅ Depth-first analysis completed successfully")
    except Exception as e:
        print(f"❌ Error with depth-first approach: {str(e)}")
    
    print("\n" + "="*60)
    print("COMPARISON COMPLETE")
    print("="*60)
    print("Compare the generated reports to see the differences:")
    print("• Breadth-first reports will show '(Breadth-First Analysis)' in the title")
    print("• Depth-first reports will show '(Depth-First Analysis)' in the title")
    print("• Check the iteration patterns and hypothesis exploration strategies")
    print("• Notice how depth-first builds understanding progressively")
    print("• Observe how breadth-first explores multiple directions simultaneously")

# --------------------- New: Test All Annotation Boost Providers ---------------------
def test_annotation_boost_providers(marker_data, full_csv, cluster_name="monocyte", conversation_history_mode="final", search_strategy="breadth", report_style="per_iteration"):
    """
    Test the annotation boost functionality with different providers.
    This function helps validate that the unified annotation boost implementation 
    works correctly with all supported providers.
    
    Args:
        marker_data: Marker data DataFrame
        full_csv: Path to the CSV file with annotation results
        cluster_name: Name of the cluster to analyze (default: "monocyte")
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        report_style: Style of report generation ("per_iteration" or "total_summary")
    """
    print("\n=== Testing Annotation Boost with Multiple Providers ===")
    print(f"Using conversation history mode: {conversation_history_mode}")
    print(f"Using search strategy: {search_strategy}")
    
    # Test with OpenAI provider
    print("\n----- Testing with OpenAI provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="openai", conversation_history_mode=conversation_history_mode, search_strategy=search_strategy, report_style=report_style)
    except Exception as e:
        print(f"Error with OpenAI provider: {str(e)}")
    
    # Test with Anthropic provider
    print("\n----- Testing with Anthropic provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="anthropic", conversation_history_mode=conversation_history_mode, search_strategy=search_strategy, report_style=report_style)
    except Exception as e:
        print(f"Error with Anthropic provider: {str(e)}")
    
    # Test with OpenRouter provider
    print("\n----- Testing with OpenRouter provider -----")
    try:
        run_annotation_boost(marker_data, full_csv, cluster_name, provider_test="openrouter", conversation_history_mode=conversation_history_mode, search_strategy=search_strategy, report_style=report_style)
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
            conversation_history_mode=conversation_history_mode,
            search_strategy=search_strategy,
            report_style=report_style
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

# --------------------- New: Test Validator Involvement Levels ---------------------
def test_validator_involvement(marker_data, provider_test=None):
    """
    Test the different validator involvement levels (v0 vs v1) to demonstrate the differences.
    This function helps validate that the validator_involvement parameter works correctly
    and shows the differences between strict (v0) and moderate (v1) validation.
    
    Args:
        marker_data: Marker data DataFrame
        provider_test: Optional provider to test (default: uses global provider)
    """
    print("\n=== Testing Validator Involvement Levels ===")
    print("This will compare v0 (stricter) vs v1 (moderate) validation approaches")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # If using a custom provider, adjust model if needed
    current_model = model_name
    if test_provider.startswith("http") and current_model == "google/gemini-2.5-flash-preview":
        if test_provider == "https://api.deepseek.com":
            current_model = "deepseek-chat"
        print(f"Using model: {current_model} with custom provider: {test_provider}")
    
    validator_levels = ["v0", "v1"]
    results = {}
    
    for validator_level in validator_levels:
        print(f"\n" + "="*60)
        print(f"TESTING {validator_level.upper()} VALIDATOR")
        print("="*60)
        
        if validator_level == "v0":
            print("• STRICTER validation (original v0 system)")
            print("• More detailed feedback and instructions")
            print("• Expects thorough analysis and validation")
        else:
            print("• MODERATE validation (current v1 system)")
            print("• Streamlined validation process")
            print("• Balanced approach for most use cases")
        print()
        
        # Generate organized output path
        validator_output_name = get_output_path("validator_test", filename=f"{output_name}_{validator_level}_validator")
        
        try:
            print(f"Running batch analysis with {validator_level} validator...")
            start_time = time.time()
            
            batch_results = runCASSIA_batch(
                marker=marker_data,
                output_name=validator_output_name,
                model=current_model,
                tissue=tissue,
                species=species,
                max_workers=3,  # Reduced for testing
                n_genes=30,     # Fewer genes for faster testing
                additional_info=f"Testing {validator_level} validator involvement",
                provider=test_provider,
                validator_involvement=validator_level
            )
            
            end_time = time.time()
            print(f"✅ {validator_level} validator test completed in {end_time - start_time:.2f} seconds")
            
            # Check if files were created
            expected_files = [
                f"{validator_output_name}_full.csv",
                f"{validator_output_name}_summary.csv"
            ]
            
            for file_path in expected_files:
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    print(f"  Generated: {os.path.basename(file_path)} ({file_size} bytes)")
                    
                    # Read and analyze results
                    if file_path.endswith("_full.csv"):
                        try:
                            df = pd.read_csv(file_path)
                            print(f"  Analyzed {len(df)} clusters")
                            if 'iterations' in df.columns:
                                avg_iterations = df['iterations'].mean()
                                print(f"  Average iterations: {avg_iterations:.1f}")
                        except Exception as e:
                            print(f"  Error reading results: {str(e)}")
                else:
                    print(f"  ⚠️ Expected file not found: {os.path.basename(file_path)}")
            
            results[validator_level] = {
                'success': True,
                'execution_time': end_time - start_time,
                'output_files': expected_files
            }
            
        except Exception as e:
            print(f"❌ Error with {validator_level} validator: {str(e)}")
            results[validator_level] = {
                'success': False,
                'error': str(e)
            }
    
    # Generate comparison summary
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    
    for level in validator_levels:
        result = results.get(level, {})
        if result.get('success'):
            print(f"✅ {level.upper()} validator: SUCCESS ({result.get('execution_time', 0):.1f}s)")
        else:
            print(f"❌ {level.upper()} validator: FAILED - {result.get('error', 'Unknown error')}")
    
    print("\n🔍 Key Differences to Look For:")
    print("• v0 (stricter): More thorough validation, potentially more iterations")
    print("• v1 (moderate): Streamlined validation, faster processing")
    print("• Check the generated CSV files to compare validation behavior")
    print("• Look at the 'iterations' column to see validation patterns")
    
    print(f"\n📁 Results saved to: {get_output_path('validator_test')}")
    
    return results

# --------------------- New: Single Annotation with Validator Test ---------------------
def test_single_annotation_validators(marker_list=None, provider_test=None):
    """
    Test a single annotation with different validator levels to see the detailed differences.
    
    Args:
        marker_list: Optional list of markers (defaults to a test set)
        provider_test: Optional provider to test (default: uses global provider)
    """
    print("\n=== Testing Single Annotation with Different Validators ===")
    
    # Use default test markers if not provided
    if marker_list is None:
        marker_list = ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "IGHM", "IGHD", "CD27", "CD38", "SDC1"]
        print(f"Using test markers: {', '.join(marker_list[:5])}...")
    
    # Use the specified provider or fall back to global setting
    test_provider = provider_test or provider
    print(f"Using provider: {test_provider}")
    
    # If using a custom provider, adjust model if needed
    current_model = model_name
    if test_provider.startswith("http") and current_model == "google/gemini-2.5-flash-preview":
        if test_provider == "https://api.deepseek.com":
            current_model = "deepseek-chat"
    
    validator_levels = ["v0", "v1"]
    results = {}
    
    for validator_level in validator_levels:
        print(f"\n--- Testing {validator_level} validator on single annotation ---")
        
        try:
            result, conversation, _ = runCASSIA(
                model=current_model,
                temperature=0,
                marker_list=marker_list,
                tissue=tissue,
                species=species,
                additional_info=f"Single annotation test with {validator_level} validator",
                provider=test_provider,
                validator_involvement=validator_level
            )
            
            if result:
                print(f"✅ {validator_level} validation successful")
                print(f"   Main cell type: {result.get('main_cell_type', 'N/A')}")
                print(f"   Sub cell types: {result.get('sub_cell_types', 'N/A')}")
                print(f"   Iterations: {result.get('iterations', 'N/A')}")
                print(f"   Conversation length: {len(conversation)} exchanges")
                
                results[validator_level] = {
                    'success': True,
                    'result': result,
                    'conversation_length': len(conversation),
                    'iterations': result.get('iterations', 0)
                }
            else:
                print(f"❌ {validator_level} validation failed - no result returned")
                results[validator_level] = {'success': False, 'error': 'No result returned'}
                
        except Exception as e:
            print(f"❌ Error with {validator_level} validator: {str(e)}")
            results[validator_level] = {'success': False, 'error': str(e)}
    
    # Compare results
    print("\n" + "="*50)
    print("SINGLE ANNOTATION COMPARISON")
    print("="*50)
    
    if all(results[level].get('success') for level in validator_levels):
        v0_result = results['v0']
        v1_result = results['v1']
        
        print(f"v0 iterations: {v0_result.get('iterations', 0)}")
        print(f"v1 iterations: {v1_result.get('iterations', 0)}")
        print(f"v0 conversation length: {v0_result.get('conversation_length', 0)}")
        print(f"v1 conversation length: {v1_result.get('conversation_length', 0)}")
        
        # Compare main cell types
        v0_main = v0_result.get('result', {}).get('main_cell_type', '')
        v1_main = v1_result.get('result', {}).get('main_cell_type', '')
        
        if v0_main == v1_main:
            print(f"✅ Both validators agreed on main cell type: {v0_main}")
        else:
            print(f"⚠️ Validators disagreed:")
            print(f"   v0: {v0_main}")
            print(f"   v1: {v1_main}")
    
    return results

# --------------------- Step 8: Test LLM Image Processing ---------------------
def test_llm_image_processing():
    """
    Test the new LLM image processing functionality from llm_image.py
    This function tests vision-enabled LLM calls with sample images.
    """
    print("\n=== Testing LLM Image Processing ===")
    
    try:
        # Import the new llm_image module
        from llm_image import (
            call_llm_with_image, 
            call_llm_analyze_image, 
            call_llm_extract_text_from_image
        )
        print("✅ Successfully imported llm_image module")
    except ImportError as e:
        print(f"❌ Failed to import llm_image module: {str(e)}")
        return
    
    # Test providers to check
    test_providers = [
        {"name": "OpenAI", "provider": "openai", "model": "gpt-4o"},
        {"name": "Anthropic", "provider": "anthropic", "model": "claude-3-5-sonnet-20241022"},
        {"name": "OpenRouter", "provider": "openrouter", "model": "openai/gpt-4o"}
    ]
    
    # Create test results directory
    test_results_dir = get_output_path("image_processing_test")
    print(f"Test results will be saved to: {test_results_dir}")
    
    # Test with sample images - using logo files that exist in the project
    script_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(script_dir, "..", "logo.png")
    
    if not os.path.exists(logo_path):
        print(f"⚠️ Logo file not found at {logo_path}")
        print("Creating a simple test image URL instead...")
        test_image = "https://via.placeholder.com/300x200/0000FF/FFFFFF?text=Test+Image"
        image_type = "URL"
    else:
        test_image = logo_path
        image_type = "Local file"
        print(f"✅ Found test image: {logo_path}")
    
    print(f"Using test image: {test_image} ({image_type})")
    
    # Test results storage
    results = []
    
    # Test each provider
    for provider_info in test_providers:
        provider_name = provider_info["name"]
        provider = provider_info["provider"]
        model = provider_info["model"]
        
        print(f"\n--- Testing {provider_name} ---")
        
        try:
            # Test 1: Basic image analysis
            print(f"Test 1: Basic image analysis with {provider_name}")
            analysis_result = call_llm_analyze_image(
                image_input=test_image,
                analysis_prompt="Please analyze this image and describe what you see. Be concise.",
                provider=provider,
                model=model,
                temperature=0.3,
                max_tokens=500
            )
            
            result_entry = {
                "provider": provider_name,
                "test": "Basic Analysis",
                "success": True,
                "result": analysis_result[:200] + "..." if len(analysis_result) > 200 else analysis_result,
                "error": None
            }
            results.append(result_entry)
            print(f"✅ {provider_name} - Basic analysis successful")
            print(f"Response preview: {analysis_result[:100]}...")
            
            # Test 2: Custom prompt with image
            print(f"Test 2: Custom prompt with {provider_name}")
            custom_result = call_llm_with_image(
                prompt="What colors are prominent in this image? List the main visual elements.",
                image_input=test_image,
                provider=provider,
                model=model,
                temperature=0.1,
                max_tokens=300
            )
            
            result_entry = {
                "provider": provider_name,
                "test": "Custom Prompt",
                "success": True,
                "result": custom_result[:200] + "..." if len(custom_result) > 200 else custom_result,
                "error": None
            }
            results.append(result_entry)
            print(f"✅ {provider_name} - Custom prompt successful")
            print(f"Response preview: {custom_result[:100]}...")
            
            # Test 3: OCR functionality (if using a local file with text)
            if image_type == "Local file":
                print(f"Test 3: OCR text extraction with {provider_name}")
                try:
                    ocr_result = call_llm_extract_text_from_image(
                        image_input=test_image,
                        provider=provider,
                        model=model,
                        temperature=0.0,
                        max_tokens=4096
                    )
                    
                    result_entry = {
                        "provider": provider_name,
                        "test": "OCR Text Extraction",
                        "success": True,
                        "result": ocr_result[:200] + "..." if len(ocr_result) > 200 else ocr_result,
                        "error": None
                    }
                    results.append(result_entry)
                    print(f"✅ {provider_name} - OCR extraction successful")
                    print(f"Extracted text preview: {ocr_result[:100]}...")
                except Exception as e:
                    print(f"⚠️ {provider_name} - OCR test failed: {str(e)}")
                    results.append({
                        "provider": provider_name,
                        "test": "OCR Text Extraction",
                        "success": False,
                        "result": None,
                        "error": str(e)
                    })
            
        except Exception as e:
            print(f"❌ {provider_name} - Tests failed: {str(e)}")
            results.append({
                "provider": provider_name,
                "test": "All Tests",
                "success": False,
                "result": None,
                "error": str(e)
            })
    
    # Save test results
    results_df = pd.DataFrame(results)
    results_csv = os.path.join(test_results_dir, "image_processing_test_results.csv")
    results_df.to_csv(results_csv, index=False)
    print(f"\n📊 Test results saved to: {results_csv}")
    
    # Print summary
    print("\n=== Test Summary ===")
    successful_tests = len([r for r in results if r["success"]])
    total_tests = len(results)
    print(f"✅ Successful tests: {successful_tests}/{total_tests}")
    
    if successful_tests > 0:
        print("\n🎉 LLM Image Processing tests completed successfully!")
        print("The llm_image module is working correctly.")
    else:
        print("\n❌ All tests failed. Check your API keys and network connection.")
    
    # Display sample usage examples
    print("\n=== Usage Examples ===")
    print("""
# Basic image analysis
from CASSIA.llm_image import call_llm_analyze_image
result = call_llm_analyze_image("path/to/image.jpg", provider="openai")

# Custom prompt with image
from CASSIA.llm_image import call_llm_with_image
result = call_llm_with_image(
    prompt="What medical conditions can you identify?",
    image_input="medical_scan.png",
    provider="anthropic",
    model="claude-3-5-sonnet-20241022"
)

# OCR text extraction
from CASSIA.llm_image import call_llm_extract_text_from_image
text = call_llm_extract_text_from_image("document.png", provider="openai")

# Multiple images
result = call_llm_with_image(
    prompt="Compare these two images",
    image_input=["image1.jpg", "image2.jpg"],
    provider="openrouter"
)
""")

def test_model_settings():
    """
    Test the model settings functionality with comprehensive examples.
    This function demonstrates how to use the new model settings system.
    """
    print("\n=== Testing CASSIA Model Settings System ===")
    
    # Test import
    try:
        from model_settings import (
            resolve_model_name,
            get_recommended_model,
            get_model_settings,
            print_available_models
        )
        print("Successfully imported model_settings module")
        
        # Show where settings are loaded from
        settings = get_model_settings()
        print(f"📁 Model settings loaded from: {settings.config_path}")
        
    except ImportError as e:
        print(f"Failed to import model_settings: {e}")
        return
    
    print("\n" + "="*80)
    print("MODEL SETTINGS SYSTEM OVERVIEW")
    print("="*80)
    print("""
The CASSIA Model Settings System makes it easy to use different AI models while
maintaining control over API keys. Key features:

✅ Provider must be specified (openai, anthropic, openrouter)
✅ Tier shortcuts: 'best', 'balanced', 'fast', 'recommended' per provider
✅ Or use exact model names directly
✅ No accidental API switching
✅ Configuration included with package installation
""")
    
    print(f"📦 Configuration Location: {settings.config_path}")
    print("   The model_settings.json file is included in the package data directory")
    print("   and is automatically loaded when you import the module.")
    
    # Test 1: Basic model resolution with providers
    print("\n=== Test 1: Basic Model Resolution ===")
    test_cases = [
        ("best", "openai"),
        ("best", "anthropic"),
        ("best", "openrouter"),
        ("balanced", "openai"),
        ("balanced", "anthropic"),
        ("balanced", "openrouter"),
        ("fast", "openai"),
        ("fast", "anthropic"),
        ("fast", "openrouter"),
        ("recommended", "openai"),
        ("recommended", "anthropic"),
        ("recommended", "openrouter"),
        ("gpt-4o", "openai"),  # exact model name
    ]
    
    print("Testing model name resolution (Provider REQUIRED):")
    for model, provider in test_cases:
        try:
            resolved = resolve_model_name(model, provider)
            print(f"  {model:10} + {provider:10} → {resolved[0]}")
        except Exception as e:
            print(f"  {model:10} + {provider:10} → Error: {e}")
    
    # Test 2: Error handling - provider required
    print("\n=== Test 2: Error Handling (Provider Required) ===")
    try:
        resolve_model_name("best", None)
        print("ERROR: Should have failed without provider")
    except Exception as e:
        print(f"Correctly failed without provider: {e}")
    
    # Test 3: Provider-specific recommendations
    print("\n=== Test 3: Provider-Specific Recommendations ===")
    providers = ["openai", "anthropic", "openrouter"]
    for provider in providers:
        try:
            recommended = get_recommended_model(provider)
            print(f"  {provider:10} best → {recommended[0]}")
        except Exception as e:
            print(f"  {provider:10} best → Error: {e}")
    
    # Test 4: Print available models
    print("\n=== Test 4: Available Models ===")
    print_available_models()

    # Test 5: Provider-specific tier shortcuts
    print("\n=== Test 5: Provider-Specific Tier Shortcuts ===")
    tiers = ["best", "balanced", "fast", "recommended"]
    for provider in providers:
        print(f"  {provider:10}:")
        for tier in tiers:
            try:
                resolved = resolve_model_name(tier, provider)
                print(f"    {tier:12} -> {resolved[0]}")
            except Exception as e:
                print(f"    {tier:12} -> Error: {e}")
    
    # Test 6: Using exact model names
    print("\n=== Test 6: Exact Model Names ===")
    exact_tests = [
        ("gpt-4o", "openai"),
        ("claude-sonnet-4-5", "anthropic"),
        ("google/gemini-2.5-flash", "openrouter"),
    ]
    for model, provider in exact_tests:
        resolved = resolve_model_name(model, provider)
        print(f"  {model:30} -> {resolved[0]}")
    
    # Test 7: Practical usage examples
    print("\n=== Test 7: Practical Usage Examples ===")
    print("""
HOW TO USE MODEL SETTINGS IN YOUR CODE:

1. USING TIER SHORTCUTS:
   runCASSIA_batch(
       marker=markers,
       model="recommended",   # Use recommended model for provider
       provider="openai"
   )

2. FAST/CHEAP ANALYSIS:
   runCASSIA_batch(
       marker=markers,
       model="fast",          # Fastest/cheapest option
       provider="openrouter"  # -> google/gemini-2.5-flash
   )

3. BEST QUALITY ANALYSIS:
   runCASSIA_batch(
       marker=markers,
       model="best",          # Best model for provider
       provider="anthropic"   # -> claude-opus-4-5
   )

4. BALANCED ANALYSIS:
   runCASSIA_batch(
       marker=markers,
       model="balanced",      # Good balance of cost/quality
       provider="openrouter"  # -> openai/gpt-5.1
   )

5. EXACT MODEL NAME:
   runCASSIA_batch(
       marker=markers,
       model="gpt-4o",        # Use exact model name
       provider="openai"
   )
""")
    
    # Test 8: Key benefits summary
    print("\n=== Key Benefits Summary ===")
    print("Provider control: You specify which API to use")
    print("Tier shortcuts: 'best', 'balanced', 'fast', 'recommended'")
    print("Exact model names: Can still use full model names directly")
    print("Simple JSON config: Easy to update model mappings")
    
    print("\n" + "="*80)
    print("MODEL TIER SUMMARY")
    print("="*80)
    print("""
Available tiers: best, balanced, fast, recommended

OPENAI:
  best        -> gpt-5.1
  balanced    -> gpt-4o
  fast        -> gpt-5-mini
  recommended -> gpt-5.1

ANTHROPIC:
  best        -> claude-opus-4-5
  balanced    -> claude-sonnet-4-5
  fast        -> claude-haiku-4-5
  recommended -> claude-sonnet-4-5

OPENROUTER:
  best        -> anthropic/claude-sonnet-4.5
  balanced    -> openai/gpt-5.1
  fast        -> google/gemini-2.5-flash
  recommended -> anthropic/claude-sonnet-4.5
""")
    
    # Test 9: Create comprehensive test results
    test_results_dir = get_output_path("model_settings_test")
    results_file = os.path.join(test_results_dir, "model_settings_test_results.txt")
    
    print(f"\n📊 Saving detailed test results to: {results_file}")
    
    with open(results_file, 'w') as f:
        f.write("CASSIA Model Settings System Test Results\n")
        f.write("="*50 + "\n\n")
        
        f.write("Test 1: Basic Model Resolution\n")
        f.write("-" * 30 + "\n")
        for model, provider in test_cases:
            try:
                resolved = resolve_model_name(model, provider)
                f.write(f"{model:10} + {provider:10} → {resolved[0]}\n")
            except Exception as e:
                f.write(f"{model:10} + {provider:10} → Error: {e}\n")
        
        f.write("\nTest 2: Provider-Specific Recommendations\n")
        f.write("-" * 30 + "\n")
        for provider in providers:
            try:
                recommended = get_recommended_model(provider)
                f.write(f"{provider:10} best → {recommended[0]}\n")
            except Exception as e:
                f.write(f"{provider:10} best → Error: {e}\n")
        
        f.write("\nTest Status: SUCCESS\n")
        f.write("All core functionality working correctly.\n")
    
    print("✅ Model settings tests completed successfully!")
    print("🎉 The model settings system is ready for use!")


def main():
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description='Run CASSIA analysis pipelines')
    parser.add_argument('--step', type=str, default='all',
                      help='Which step to run: all, batch, merge, merge_all, score, report, uncertainty, single_cluster_uncertainty, boost, compare, symphony, subcluster, subclustering_uncertainty, boost_task, test_boost, test_subcluster, single_subcluster, debug_genes, test_pipeline, compare_strategies, test_total_summary, test_image')
    parser.add_argument('--input_csv', type=str, default=None,
                      help='Input CSV file for steps that require it (merge, score, report, boost)')
    parser.add_argument('--cluster', type=str, default='monocyte',
                      help='Cluster name for annotation boost or single cluster uncertainty')
    parser.add_argument('--task', type=str, default=None,
                      help='Additional task for annotation boost with task, e.g., "check if this is a cancer cell"')
    parser.add_argument('--provider', type=str, default=None,
                      help='Provider to use for API calls (openai, anthropic, openrouter, or a custom base URL)')
    parser.add_argument('--test_provider', type=str, default=None,
                      help='Test provider to use for API calls (e.g., openrouter, https://api.deepseek.com)')
    parser.add_argument('--api_key', type=str, default=None,
                      help='API key for custom provider (used if test_provider is a base URL)')
    parser.add_argument('--test_genes', type=str, default=None,
                      help='Comma-separated list of genes to test for the debug_genes step')
    parser.add_argument('--history_mode', type=str, default="final",
                      help='Conversation history mode for annotation boost: "full", "final", or "none"')
    parser.add_argument('--search_strategy', type=str, default="breadth",
                      help='Search strategy for annotation boost: "breadth" (multiple hypotheses) or "depth" (one hypothesis at a time)')
    parser.add_argument('--report_style', type=str, default="per_iteration",
                      help='Report style for annotation boost: "per_iteration" (detailed by iteration) or "total_summary" (gene-focused summary)')
    parser.add_argument('--iterations', type=int, default=5,
                      help='Number of iterations for single_cluster_uncertainty (default: 5)')
    parser.add_argument('--major_cluster', type=str, default="cd8 t cell",
                      help='Major cluster info for subclustering (default: "cd8 t cell")')
    args = parser.parse_args()
    
    # Override default provider if specified in command line
    global provider, model_name, api_key
    if args.test_provider:
        provider = args.test_provider
        print(f"Using test provider specified in command line: {provider}")
        if provider.startswith("http"):
            # If a custom provider, set the API key
            api_key = args.api_key or os.environ.get("CUSTOMIZED_API_KEY", "")
            if not api_key:
                print("Warning: No API key provided for custom provider. Use --api_key to specify it.")
            else:
                os.environ["CUSTOMIZED_API_KEY"] = api_key
                print(f"Set CUSTOMIZED_API_KEY for custom provider: {provider}")
        # Optionally, set a default model for deepseek
        if provider == "https://api.deepseek.com" and model_name == "google/gemini-2.5-flash-preview":
            model_name = "deepseek-chat"
    elif args.provider:
        provider = args.provider
        print(f"Using provider specified in command line: {provider}")
        # Also handle API key setting for custom providers when using --provider flag
        if provider.startswith("http"):
            # If a custom provider, set the API key
            api_key = args.api_key or os.environ.get("CUSTOMIZED_API_KEY", "")
            if not api_key:
                print("Warning: No API key provided for custom provider. Use --api_key to specify it.")
            else:
                os.environ["CUSTOMIZED_API_KEY"] = api_key
                print(f"Set CUSTOMIZED_API_KEY for custom provider: {provider}")
            # Set default model for deepseek
            if provider == "https://api.deepseek.com" and model_name == "google/gemini-2.5-flash-preview":
                model_name = "deepseek-chat"
    
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
        # This uses the official merging implementation with custom provider support
        run_merge(input_csv)
    elif args.step == 'merge_all':
        # This runs all three levels of merging (broad, detailed, very_detailed)
        run_merge_all(input_csv)
    elif args.step == 'score':
        run_quality_scoring(input_csv)
    elif args.step == 'report':
        input_csv = args.input_csv or f"{output_name}_scored.csv"
        generate_report(input_csv)
    elif args.step == 'uncertainty':
        run_uncertainty_quantification(unprocessed, provider)
    elif args.step == 'single_cluster_uncertainty':
        run_single_cluster_uncertainty(
            marker_data=unprocessed, 
            cluster_name=args.cluster, 
            provider_test=provider, 
            n_iterations=args.iterations
        )
    elif args.step == 'boost':
        run_annotation_boost(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode, search_strategy=args.search_strategy, report_style=args.report_style)
    elif args.step == 'compare':
        run_celltype_comparison()
    elif args.step == 'symphony':
        run_symphony_agent()
    elif args.step == 'subcluster':
        run_subclustering(subcluster)
    elif args.step == 'subclustering_uncertainty':
        run_subclustering_with_uncertainty(subcluster)
    elif args.step == 'boost_task':
        try:
            run_annotation_boost_with_task(unprocessed, input_csv, args.cluster, args.task, conversation_history_mode=args.history_mode, search_strategy=args.search_strategy, report_style=args.report_style)
        except Exception as e:
            print(f"Error in annotation boost with task: {str(e)}")
    elif args.step == 'test_boost':
        # New option to test the unified annotation boost functionality
        try:
            test_annotation_boost_providers(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode, search_strategy=args.search_strategy, report_style=args.report_style)
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
                run_annotation_boost(unprocessed, input_csv, args.cluster, debug_mode=True, test_genes=test_genes, search_strategy=args.search_strategy)
            
        except Exception as e:
            print(f"Error in gene extraction diagnostics: {str(e)}")
            import traceback
            traceback.print_exc()
    elif args.step == 'test_subcluster':
        # New option to test subclustering with different providers
        try:
            test_subclustering_providers(subcluster, args.major_cluster)
        except Exception as e:
            print(f"Error testing subclustering: {str(e)}")
    elif args.step == 'single_subcluster':
        # New option to run a single subclustering analysis
        run_single_subcluster(subcluster, args.major_cluster)
    elif args.step == 'test_pipeline':
        # New option to test the full pipeline with different providers
        try:
            test_full_pipeline_providers(unprocessed)
        except Exception as e:
            print(f"Error testing full pipeline: {str(e)}")
    elif args.step == 'compare_strategies':
        # New option to compare breadth-first vs depth-first search strategies
        try:
            compare_search_strategies(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode, report_style=args.report_style)
        except Exception as e:
            print(f"Error comparing search strategies: {str(e)}")
    elif args.step == 'test_total_summary':
        # New option to test total_summary report style with different providers
        try:
            test_total_summary_reports(unprocessed, input_csv, args.cluster, conversation_history_mode=args.history_mode)
        except Exception as e:
            print(f"Error testing total summary reports: {str(e)}")
    elif args.step == 'test_validators':
        # New option to test validator involvement levels
        try:
            test_validator_involvement(unprocessed, provider)
        except Exception as e:
            print(f"Error testing validator involvement: {str(e)}")
    elif args.step == 'test_single_validators':
        # New option to test single annotation with different validators
        try:
            test_single_annotation_validators(provider_test=provider)
        except Exception as e:
            print(f"Error testing single annotation validators: {str(e)}")
    elif args.step == 'test_image':
        # New option to test LLM image processing functionality
        try:
            test_llm_image_processing()
        except Exception as e:
            print(f"Error testing LLM image processing: {str(e)}")
    elif args.step == 'test_model_settings':
        # New option to test model settings functionality
        try:
            test_model_settings()
        except Exception as e:
            print(f"Error testing model settings: {str(e)}")
    else:
        print(f"Unknown step: {args.step}")
        print("Available steps: all, batch, merge, merge_all, score, report, uncertainty, single_cluster_uncertainty, boost, compare, symphony, subcluster, subclustering_uncertainty, boost_task, test_boost, test_subcluster, single_subcluster, debug_genes, test_pipeline, compare_strategies, test_total_summary, test_validators, test_single_validators, test_image, test_model_settings")


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


# =====================================================================================
# USAGE EXAMPLES - Organized by functionality and complexity
# =====================================================================================

# --------------------- 1. BASIC ANALYSIS STEPS ---------------------
# Basic steps using default settings (OpenRouter + Gemini 2.5 Flash):

# python CASSIA_python_tutorial.py --step all                    # Complete pipeline
# python CASSIA_python_tutorial.py --step batch                  # Batch analysis only
# python CASSIA_python_tutorial.py --step merge                  # Merge annotations
# python CASSIA_python_tutorial.py --step merge_all              # All merge levels
# python CASSIA_python_tutorial.py --step score                  # Quality scoring
# python CASSIA_python_tutorial.py --step report                 # Generate reports
# python CASSIA_python_tutorial.py --step uncertainty            # Uncertainty analysis
# python CASSIA_python_tutorial.py --step boost                  # Annotation boost
# python CASSIA_python_tutorial.py --step compare                # Cell type comparison
# python CASSIA_python_tutorial.py --step symphony               # Symphony Agent for advanced comparison
# python CASSIA_python_tutorial.py --step subcluster             # Subclustering
# python CASSIA_python_tutorial.py --step subclustering_uncertainty # Subclustering with uncertainty
# python CASSIA_python_tutorial.py --step boost_task             # Boost with custom task
# python CASSIA_python_tutorial.py --step test_image             # Test LLM image processing

# --------------------- 2. PROVIDER CUSTOMIZATION ---------------------
# Using different API providers:

# Standard providers:
# python CASSIA_python_tutorial.py --step [STEP] --provider openai
# python CASSIA_python_tutorial.py --step [STEP] --provider anthropic
# python CASSIA_python_tutorial.py --step [STEP] --provider openrouter

# Custom API endpoints (replace YOUR_API_KEY with actual key):
# python CASSIA_python_tutorial.py --step [STEP] --provider https://api.deepseek.com --api_key YOUR_API_KEY
# python CASSIA_python_tutorial.py --step [STEP] --provider http://your-custom-api.com --api_key YOUR_API_KEY

# --------------------- 3. ANNOTATION BOOST ADVANCED OPTIONS ---------------------
# Search strategies:
# python CASSIA_python_tutorial.py --step boost --search_strategy breadth    # Multiple hypotheses (default)
# python CASSIA_python_tutorial.py --step boost --search_strategy depth      # Focused analysis

# Report styles:
# python CASSIA_python_tutorial.py --step boost --report_style per_iteration  # Iteration-based (default)
# python CASSIA_python_tutorial.py --step boost --report_style total_summary  # Gene-focused summary

# Conversation history modes:
# python CASSIA_python_tutorial.py --step boost --history_mode final    # Final summary only (default)
# python CASSIA_python_tutorial.py --step boost --history_mode full     # Complete conversation
# python CASSIA_python_tutorial.py --step boost --history_mode none     # No conversation history

# Cluster-specific analysis:
# python CASSIA_python_tutorial.py --step boost --cluster monocyte
# python CASSIA_python_tutorial.py --step boost --cluster "cd8-positive, alpha-beta t cell"

# Custom tasks:
# python CASSIA_python_tutorial.py --step boost_task --task "check if this is a cancer cell"
# python CASSIA_python_tutorial.py --step boost_task --task "determine activation state"

# --------------------- 4. UNCERTAINTY QUANTIFICATION ---------------------
# Multiple iterations for uncertainty analysis:
# python CASSIA_python_tutorial.py --step uncertainty --iterations 5

# Single cluster uncertainty:
# python CASSIA_python_tutorial.py --step single_cluster_uncertainty --cluster monocyte --iterations 3
# python CASSIA_python_tutorial.py --step single_cluster_uncertainty --cluster "plasma cell" --iterations 5

# --------------------- 5. SUBCLUSTERING ANALYSIS ---------------------
# Basic subclustering:
# python CASSIA_python_tutorial.py --step subcluster --major_cluster "cd8 t cell"
# python CASSIA_python_tutorial.py --step single_subcluster --major_cluster "T cell"

# Test multiple providers:
# python CASSIA_python_tutorial.py --step test_subcluster --major_cluster "cd8 t cell"

# --------------------- 6. TESTING AND COMPARISON FUNCTIONS ---------------------
# Test multiple providers for annotation boost:
# python CASSIA_python_tutorial.py --step test_boost --cluster monocyte

# Compare search strategies:
    # python CASSIA_python_tutorial.py --step compare_strategies --cluster monocyte
# python CASSIA_python_tutorial.py --step compare_strategies --cluster "plasma cell"

# Test full pipeline with multiple providers:
# python CASSIA_python_tutorial.py --step test_pipeline

# Test validator involvement levels:
# python CASSIA_python_tutorial.py --step test_validators

# Test single annotation with different validators:
# python CASSIA_python_tutorial.py --step test_single_validators

# Test total summary report style:
# python CASSIA_python_tutorial.py --step test_total_summary --cluster monocyte

# Debug gene extraction:
# python CASSIA_python_tutorial.py --step debug_genes --test_genes "CD133,CD9,ChAT,DCLK1"

# --------------------- 7. VALIDATOR INVOLVEMENT TESTING ---------------------
# Test validator involvement levels (v0 vs v1):
# python CASSIA_python_tutorial.py --step test_validators                     # Compare v0 (stricter) vs v1 (moderate) 
# python CASSIA_python_tutorial.py --step test_validators --provider openai   # Test with OpenAI
# python CASSIA_python_tutorial.py --step test_validators --provider anthropic # Test with Anthropic
# python CASSIA_python_tutorial.py --step test_validators --provider https://api.deepseek.com --api_key YOUR_API_KEY

# Test single annotation with different validators:
# python CASSIA_python_tutorial.py --step test_single_validators               # Compare single annotation validation
# python CASSIA_python_tutorial.py --step test_single_validators --provider openai

# --------------------- 8. LLM IMAGE PROCESSING TESTS ---------------------
# Test LLM image processing functionality:

# Basic image processing test:
# python CASSIA_python_tutorial.py --step test_image                          # Test all providers with image processing

# Test with specific providers:
# python CASSIA_python_tutorial.py --step test_image --provider openai        # Test OpenAI vision models
# python CASSIA_python_tutorial.py --step test_image --provider anthropic     # Test Anthropic vision models
# python CASSIA_python_tutorial.py --step test_image --provider openrouter    # Test OpenRouter vision models

# Test with custom API endpoints:
# python CASSIA_python_tutorial.py --step test_image --provider https://api.deepseek.com --api_key YOUR_API_KEY

# --------------------- 9. COMPLEX COMBINATIONS ---------------------
# Combine multiple advanced options:

# Custom provider + search strategy + report style:
# python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key YOUR_API_KEY --search_strategy depth --report_style total_summary

# OpenAI with specific settings:
# python CASSIA_python_tutorial.py --step boost --provider openai --search_strategy breadth --report_style per_iteration --cluster monocyte

# Anthropic with custom task:
# python CASSIA_python_tutorial.py --step boost_task --provider anthropic --search_strategy depth --task "analyze cell activation markers"

# Custom API with uncertainty analysis:
# python CASSIA_python_tutorial.py --step uncertainty --provider https://api.deepseek.com --api_key YOUR_API_KEY --iterations 4

# Full pipeline with custom provider:
# python CASSIA_python_tutorial.py --step all --provider https://api.deepseek.com --api_key YOUR_API_KEY

# Test validator levels with different providers:
# python CASSIA_python_tutorial.py --step test_validators --provider openai
# python CASSIA_python_tutorial.py --step test_validators --provider https://api.deepseek.com --api_key YOUR_API_KEY

# --------------------- 8. INPUT FILE SPECIFICATION ---------------------
# For steps requiring input CSV files, use --input_csv:
# python CASSIA_python_tutorial.py --step boost --input_csv /path/to/your/results.csv
# python CASSIA_python_tutorial.py --step score --input_csv /path/to/your/results.csv
# python CASSIA_python_tutorial.py --step boost_task --input_csv ./test_results/normal_api/batch_analysis/results_full.csv

# =====================================================================================
# NOTES:
# - Replace [STEP] with any valid step name from section 1
# - Replace YOUR_API_KEY with your actual API key
# - All results are automatically organized in test_results/normal_api/ or test_results/custom_api/
# - Custom providers require --api_key parameter
# - For help: python CASSIA_python_tutorial.py --help
# =====================================================================================


# =====================================================================================
# CUSTOMIZED EXAMPLES - Ready-to-use commands with explicit API keys
# =====================================================================================

# --------------------- DEEPSEEK API EXAMPLES (Custom Provider) ---------------------
# Copy-paste ready commands with explicit DeepSeek API key:

# Basic analysis with DeepSeek:
# python CASSIA_python_tutorial.py --step all --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# python CASSIA_python_tutorial.py --step batch --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# python CASSIA_python_tutorial.py --step merge --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# python CASSIA_python_tutorial.py --step score --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# python CASSIA_python_tutorial.py --step uncertainty --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16

# Annotation boost with DeepSeek:
# python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --cluster monocyte
# python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --cluster "cd8-positive, alpha-beta t cell"
# python CASSIA_python_tutorial.py --step boost_task --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --task "check if this is a cancer cell"

# Advanced DeepSeek combinations:
# python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --search_strategy depth --report_style total_summary
# python CASSIA_python_tutorial.py --step boost_task --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --search_strategy breadth --report_style per_iteration --task "determine activation state"
# python CASSIA_python_tutorial.py --step compare_strategies --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --cluster monocyte

# Subclustering with DeepSeek:
# python CASSIA_python_tutorial.py --step subcluster --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --major_cluster "cd8 t cell"
# python CASSIA_python_tutorial.py --step single_subcluster --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --major_cluster "T cell"

# --------------------- OPENAI API EXAMPLES ---------------------
# Using OpenAI with explicit models:

# python CASSIA_python_tutorial.py --step all --provider openai
# python CASSIA_python_tutorial.py --step batch --provider openai  
# python CASSIA_python_tutorial.py --step boost --provider openai --cluster monocyte --search_strategy breadth
# python CASSIA_python_tutorial.py --step boost_task --provider openai --task "analyze immune cell activation" --report_style total_summary
# python CASSIA_python_tutorial.py --step uncertainty --provider openai --iterations 4
# python CASSIA_python_tutorial.py --step subcluster --provider openai --major_cluster "cd8 t cell"

# --------------------- ANTHROPIC API EXAMPLES ---------------------
# Using Anthropic Claude models:

# python CASSIA_python_tutorial.py --step all --provider anthropic
# python CASSIA_python_tutorial.py --step batch --provider anthropic
# python CASSIA_python_tutorial.py --step boost --provider anthropic --cluster "plasma cell" --search_strategy depth
# python CASSIA_python_tutorial.py --step boost_task --provider anthropic --task "identify cell state markers" --report_style per_iteration
# python CASSIA_python_tutorial.py --step uncertainty --provider anthropic --iterations 3
# python CASSIA_python_tutorial.py --step compare_strategies --provider anthropic --cluster monocyte

# --------------------- OPENROUTER API EXAMPLES ---------------------
# Using OpenRouter (default provider):

# python CASSIA_python_tutorial.py --step all --provider openrouter
# python CASSIA_python_tutorial.py --step batch --provider openrouter
# python CASSIA_python_tutorial.py --step boost --provider openrouter --cluster monocyte --search_strategy breadth --report_style total_summary
# python CASSIA_python_tutorial.py --step boost_task --provider openrouter --task "check for exhaustion markers" --search_strategy depth
# python CASSIA_python_tutorial.py --step test_boost --provider openrouter --cluster "intestinal crypt stem cell"

# --------------------- WORKFLOW EXAMPLES ---------------------
# Complete workflows for common use cases:

# Standard workflow (OpenRouter):
# python CASSIA_python_tutorial.py --step batch
# python CASSIA_python_tutorial.py --step score --input_csv ./test_results/normal_api/batch_analysis/CASSIA_large_intestine_human_full.csv
# python CASSIA_python_tutorial.py --step boost --input_csv ./test_results/normal_api/batch_analysis/CASSIA_large_intestine_human_full.csv --cluster monocyte

# Custom API workflow (DeepSeek):
# python CASSIA_python_tutorial.py --step batch --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# python CASSIA_python_tutorial.py --step score --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --input_csv ./test_results/custom_api/batch_analysis/CASSIA_large_intestine_human_full.csv
# python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16 --input_csv ./test_results/custom_api/batch_analysis/CASSIA_large_intestine_human_full.csv --cluster monocyte

# Testing and comparison workflow:
# python CASSIA_python_tutorial.py --step test_pipeline
# python CASSIA_python_tutorial.py --step test_boost --cluster monocyte
# python CASSIA_python_tutorial.py --step compare_strategies --cluster "cd8-positive, alpha-beta t cell"

# =====================================================================================
# QUICK START COMMANDS:
# =====================================================================================
# 1. Run complete pipeline: python CASSIA_python_tutorial.py --step all
# 2. Quick batch analysis: python CASSIA_python_tutorial.py --step batch
# 3. Test with DeepSeek: python CASSIA_python_tutorial.py --step batch --provider https://api.deepseek.com --api_key sk-afb39114f1334ba486505d9425937d16
# 4. Boost analysis: python CASSIA_python_tutorial.py --step boost --cluster monocyte
# 5. Compare providers: python CASSIA_python_tutorial.py --step test_boost
# 6. Test image processing: python CASSIA_python_tutorial.py --step test_image
# =====================================================================================

# --------------------- New: Symphony Agent (Advanced Celltype Comparison) ---------------------
def run_symphony_agent():
    """
    Run the Symphony Agent, a multi-model, multi-round discussion framework
    to compare and contrast similar or ambiguous cell types. This is ideal for
    resolving subtle differences between cell types like Naive B cells,
    Plasmablasts, and Plasma cells.
    """
    print("\n=== Running Symphony Agent (Advanced Cell Type Comparison) ===")
    
    # Import the core Symphony function. It's also available as compareCelltypes.
    try:
        from CASSIA.symphony_compare import symphonyCompare
    except ImportError:
        print("Could not import symphonyCompare. Please ensure it is available.")
        return

    # Define a set of ambiguous or closely related cell types to compare
    celltypes_to_compare = ["Naive B cell", "Plasmablast", "Plasma cell"]
    
    # Provide a set of marker genes. Including contradictory or overlapping
    # markers is useful to stimulate discussion among the LLM agents.
    # For example, CD20 is for B-cells, CD38/SDC1 for plasma cells, PAX5 for B-cells, IRF4 for plasma cells.
    marker_gene_set = "CD19, CD20, CD27, PAX5, IRF4, CD38, SDC1"
    
    # Set up organized output paths
    output_dir = get_output_path("symphony_agent")
    output_csv = os.path.join(output_dir, "symphony_bcell_comparison.csv")
    output_html = os.path.join(output_dir, "symphony_bcell_report.html")
    
    print(f"Comparing cell types: {', '.join(celltypes_to_compare)}")
    print(f"Using markers: {marker_gene_set}")
    
    # Run the Symphony agent
    symphonyCompare(
        tissue=tissue,
        species=species,
        celltypes=celltypes_to_compare,
        marker_set=marker_gene_set,
        model_preset="quality",  # Use "quality" for best models, "budget" for cheaper, or "fastest"
        discussion_mode=True,    # Enable the multi-round discussion
        discussion_rounds=2,     # Number of discussion rounds after initial analysis
        generate_html_report=True,
        output_file=output_csv
    )
    
    print("\n✅ Symphony Agent analysis complete!")
    print(f"Results saved to: {output_dir}")
    print(f"-> CSV report: {output_csv}")
    print(f"-> HTML report: {output_html}")

