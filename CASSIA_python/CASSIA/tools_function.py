import pandas as pd
import json
import re
import csv
import os
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from openai import OpenAI
try:
    from .main_function_code import *
except ImportError:
    from main_function_code import *

import requests
import threading
import numpy as np
from importlib import resources
import datetime
import shutil
from collections import Counter

# Suppress httpx and API client logs to reduce noise during batch operations
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("openai").setLevel(logging.WARNING)
logging.getLogger("anthropic").setLevel(logging.WARNING)

try:
    from .llm_utils import *
except ImportError:
    from llm_utils import *

try:
    from .model_settings import resolve_model_name, get_recommended_model
except ImportError:
    from model_settings import resolve_model_name, get_recommended_model

try:
    from .merging_annotation import *
except ImportError:
    from merging_annotation import *

try:
    from .cell_type_comparison import compareCelltypes
except ImportError:
    from cell_type_comparison import compareCelltypes

def set_openai_api_key(api_key):
    os.environ["OPENAI_API_KEY"] = api_key

def set_anthropic_api_key(api_key):
    """Set the Anthropic API key in environment variables."""
    os.environ["ANTHROPIC_API_KEY"] = api_key

def set_openrouter_api_key(api_key):
    os.environ["OPENROUTER_API_KEY"] = api_key

def set_api_key(api_key, provider="openai"):
    """
    Set the API key for the specified provider in environment variables.
    
    Args:
        api_key (str): The API key to set
        provider (str): The provider to set the key for ('openai', 'anthropic', 'openrouter', or a custom base URL)
    
    Raises:
        ValueError: If provider is not recognized
    """
    if provider.lower() == "openai":
        os.environ["OPENAI_API_KEY"] = api_key
    elif provider.lower() == "anthropic":
        os.environ["ANTHROPIC_API_KEY"] = api_key
    elif provider.lower() == "openrouter":
        os.environ["OPENROUTER_API_KEY"] = api_key
    elif provider.lower().startswith("http"):
        os.environ["CUSTERMIZED_API_KEY"] = api_key
    else:
        raise ValueError("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)")
    
def split_markers(marker_string):
    # First, try splitting by comma and space
    markers = re.split(r',\s*', marker_string)
    
    # If that results in only one marker, try splitting by comma only
    if len(markers) == 1:
        markers = marker_string.split(',')
    
    # If still only one marker, try splitting by space
    if len(markers) == 1:
        markers = marker_string.split()
    
    # Remove any empty strings
    markers = [m.strip() for m in markers if m.strip()]
    
    return markers


def _validate_ranking_parameters(df, ranking_method, ascending):
    """Validate ranking method and column existence"""
    valid_methods = ["avg_log2FC", "p_val_adj", "pct_diff", "Score"]
    if ranking_method not in valid_methods:
        raise ValueError(f"ranking_method must be one of {valid_methods}")
    
    if ranking_method == "Score" and "Score" not in df.columns:
        available_cols = [col for col in df.columns if col.lower() in ['score', 'scores']]
        if available_cols:
            suggestion = f". Did you mean '{available_cols[0]}'?"
        else:
            suggestion = ". Available numeric columns: " + ", ".join(df.select_dtypes(include=[np.number]).columns.tolist())
        raise ValueError(f"Column 'Score' not found in DataFrame{suggestion}")
    
    if ranking_method == "p_val_adj" and "p_val_adj" not in df.columns:
        raise ValueError("Column 'p_val_adj' not found in DataFrame")


def _prepare_ranking_column(df, ranking_method):
    """Prepare the ranking column, calculating if necessary"""
    df_copy = df.copy()
    
    if ranking_method == "pct_diff":
        if "pct.1" not in df_copy.columns or "pct.2" not in df_copy.columns:
            raise ValueError("Columns 'pct.1' and 'pct.2' required for pct_diff ranking")
        df_copy["pct_diff"] = df_copy["pct.1"] - df_copy["pct.2"]
        return df_copy, "pct_diff"
    
    return df_copy, ranking_method


def _get_sort_direction(ranking_method, ascending):
    """Get the sort direction for the ranking method"""
    DEFAULT_SORT_DIRECTIONS = {
        "avg_log2FC": False,    # Higher is better
        "p_val_adj": True,      # Lower p-value is better
        "pct_diff": False,      # Higher difference is better
        "Score": False          # Higher score is better (default)
    }
    
    if ascending is not None:
        return ascending
    return DEFAULT_SORT_DIRECTIONS.get(ranking_method, False)


def get_top_markers(df, n_genes=10, format_type=None, ranking_method="avg_log2FC", ascending=None):
    """
    Get top markers from either Seurat or Scanpy differential expression results.
    
    Args:
        df: Either a pandas DataFrame (Seurat format) or dictionary (Scanpy format)
        n_genes: Number of top genes to return per cluster
        format_type: Either 'seurat', 'scanpy', or None (auto-detect)
        ranking_method: Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending: Sort direction (None uses default for each method)
    
    Returns:
        pandas DataFrame with cluster and marker columns
    """
    # Auto-detect format if not specified
    if format_type is None:
        if 'names' in df and 'scores' in df:
            format_type = 'scanpy'
        else:
            format_type = 'seurat'
    
    if format_type == 'scanpy':
        # Process Scanpy format
        clusters = df['names'].dtype.names
        top_markers = []
        
        for cluster in clusters:
            # Get data for this cluster
            genes = df['names'][cluster]
            logfc = df['logfoldchanges'][cluster]
            pvals_adj = df['pvals_adj'][cluster]
            pcts = df['pcts'][cluster]  # Get percentage information
            
            # Create temporary DataFrame for sorting
            cluster_df = pd.DataFrame({
                'gene': genes,
                'avg_log2FC': logfc,
                'p_val_adj': pvals_adj,
                'pct.1': pcts,  # Assuming this represents pct.1
                'pct.2': np.zeros_like(pcts)  # May need to adjust based on data structure
            })
            
            # Filter for significant upregulated genes with PCT threshold
            mask = (pvals_adj < 0.05) & (logfc > 0.25) & (pcts >= 0.1)  # Add PCT filter
            filtered_df = cluster_df[mask]
            
            if not filtered_df.empty:
                # Validate parameters and prepare ranking
                _validate_ranking_parameters(filtered_df, ranking_method, ascending)
                df_prepared, sort_column = _prepare_ranking_column(filtered_df, ranking_method)
                sort_ascending = _get_sort_direction(ranking_method, ascending)
                
                # Sort and get top n genes
                top_genes_df = (df_prepared
                               .sort_values(sort_column, ascending=sort_ascending, na_position='last')
                               .head(n_genes))
                valid_genes = top_genes_df['gene'].values
                
                # Join genes with commas
                markers = ','.join(valid_genes)
                top_markers.append({
                    'cluster': cluster,
                    'markers': markers
                })
        
        return pd.DataFrame(top_markers)
    
    else:  # Seurat format
        # Convert string \'inf\' and \'-inf\' to numeric values first
        df['avg_log2FC'] = pd.to_numeric(df['avg_log2FC'].replace({'inf': np.inf, '-inf': -np.inf}), errors='coerce')
        
        # Replace inf and -inf values with max and min non-inf values respectively
        max_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).max()
        min_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).min()
        df['avg_log2FC'] = df['avg_log2FC'].replace([np.inf, -np.inf], [max_non_inf, min_non_inf])
        
        # Filter by adjusted p-value, positive log2FC, and PCT
        df_filtered = df[
            (df['p_val_adj'] < 0.05) & 
            (df['avg_log2FC'] > 0.25) &
            ((df['pct.1'] >=0.1) | (df['pct.2'] >=0.1))  # Add PCT filter
        ].copy()
        
        # Validate parameters and prepare ranking
        _validate_ranking_parameters(df_filtered, ranking_method, ascending)
        df_prepared, sort_column = _prepare_ranking_column(df_filtered, ranking_method)
        sort_ascending = _get_sort_direction(ranking_method, ascending)
        
        # Sort within each cluster by specified method and get top n genes
        top_markers = []
        
        for cluster in df_prepared['cluster'].unique():
            cluster_data = df_prepared[df_prepared['cluster'] == cluster]
            # Sort by specified column and direction, then take top n
            top_n = (cluster_data
                    .sort_values(sort_column, ascending=sort_ascending, na_position='last')
                    .head(n_genes))
            
            # Handle NaN values warning
            nan_count = cluster_data[sort_column].isna().sum()
            if nan_count > 0:
                print(f"Warning: {nan_count} NaN values found in {sort_column} column for cluster {cluster}, they will be placed at the end")
            
            top_markers.append(top_n)
        
        # Combine all results
        if top_markers:
            top_markers = pd.concat(top_markers, ignore_index=True)
            
            # Create markers column by concatenating genes in order
            result = (top_markers
                     .groupby('cluster',observed=True)
                     .agg({'gene': lambda x: ','.join(x)})
                     .rename(columns={'gene': 'markers'})
                     .reset_index())
            
            return result
        else:
            return pd.DataFrame(columns=['cluster', 'markers'])
        

def check_formatted_output(structured_output):
    return 'main_cell_type' in structured_output and 'sub_cell_types' in structured_output


def rerun_formatting_agent(agent, full_conversation_history):
    full_text = "\n\n".join([f"{role}: {message}" for role, message in full_conversation_history])
    formatted_result = agent(full_text, "user")
    return extract_json_from_reply(formatted_result)


def safe_get(dict_obj, *keys):
    """Safely get nested dictionary values"""
    for key in keys:
        if isinstance(dict_obj, dict) and key in dict_obj:
            dict_obj = dict_obj[key]
        else:
            return None
    return dict_obj


def natural_sort_key(cell_type):
    """
    Create a sort key that handles numeric cluster names properly.

    Handles various formats:
    - Pure numbers: "0", "1", "10" → sorted as integers 0, 1, 10 (priority 0)
    - "cluster X": "cluster 0", "cluster 1", "cluster 10" → sorted by X numerically (priority 1)
    - "Cluster X": case-insensitive
    - Other text: sorted alphabetically (priority 2)

    Args:
        cell_type (str): The cell type or cluster name

    Returns:
        tuple: (sort_priority, numeric_value, string_value) for proper sorting
    """
    import re

    if not cell_type or not isinstance(cell_type, str):
        return (3, 0, str(cell_type))  # Non-string values go last

    cell_type_str = str(cell_type).strip()

    # Try to parse as pure integer
    try:
        return (0, int(cell_type_str), "")  # Pure numbers have priority 0
    except ValueError:
        pass

    # Try to extract number from "cluster X" or "Cluster X" format (case-insensitive)
    cluster_match = re.match(r'^cluster\s+(\d+)$', cell_type_str, re.IGNORECASE)
    if cluster_match:
        cluster_num = int(cluster_match.group(1))
        return (1, cluster_num, "")  # Cluster numbers have priority 1 (after pure numbers)

    # For any other text, sort alphabetically (priority 2)
    return (2, 0, cell_type_str.lower())


def clean_conversation_history(history_text):
    """
    Clean conversation history for safe CSV storage while preserving full content.

    Args:
        history_text (str): Raw conversation history text

    Returns:
        str: Cleaned text safe for CSV storage (no truncation)
    """
    if not history_text:
        return ""

    # Replace newlines with spaces (prevents row breaks in CSV/Excel)
    cleaned = history_text.replace('\n', ' ').replace('\r', ' ')

    # Collapse multiple spaces into single spaces
    cleaned = ' '.join(cleaned.split())

    # Double quotes will be handled by csv.writer's automatic escaping
    # No need to replace quotes - csv module handles this correctly

    # Return full content without truncation
    return cleaned


def write_csv(filename, headers, row_data):
    import os
    
    # Make sure the directory exists
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(row_data)


def runCASSIA(model="google/gemini-2.5-flash-preview", temperature=0, marker_list=None, tissue="lung", species="human", additional_info=None, provider="openrouter", validator_involvement="v1"):
    """
    Wrapper function to run cell type analysis using OpenAI, Anthropic, OpenRouter, or a custom OpenAI-compatible provider.
    
    Args:
        model (str): Model name to use
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for the analysis
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a base URL)
    
    Returns:
        tuple: (analysis_result, conversation_history)
    """
    if provider.lower() == "openai":
        return run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
    elif provider.lower() == "anthropic":
        return run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
    elif provider.lower() == "openrouter":
        return run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
    elif provider.lower().startswith("http"):
        api_key = os.environ.get("CUSTERMIZED_API_KEY")
        if not api_key:
            raise ValueError("CUSTERMIZED_API_KEY environment variable is not set. Please call set_api_key with your API key and provider (base URL).")
        return run_cell_type_analysis_custom(
            base_url=provider,
            api_key=api_key,
            model=model,
            temperature=temperature,
            marker_list=marker_list,
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            validator_involvement=validator_involvement
        )
    else:
        raise ValueError("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)")





def runCASSIA_batch(marker, output_name="cell_type_analysis_results.json", n_genes=50, model="google/gemini-2.5-flash-preview", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column=None, gene_column_name=None, max_workers=10, provider="openrouter", max_retries=1, ranking_method="avg_log2FC", ascending=None, validator_involvement="v1"):
    """
    Run cell type analysis on multiple clusters in parallel.
    
    Args:
        marker: Input marker data (pandas DataFrame or path to CSV file)
        output_name (str): Base name for output files  
        n_genes (int): Number of top genes to extract per cluster
        model (str): Model name to use for analysis
        temperature (float): Temperature parameter for the model
        tissue (str): Tissue type being analyzed
        species (str): Species being analyzed
        additional_info (str): Additional information for analysis
        celltype_column (str): Column name containing cell type names (default: first column)
        gene_column_name (str): Column name containing gene markers (default: second column)
        max_workers (int): Maximum number of parallel workers
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or base URL)
        max_retries (int): Maximum number of retries for failed analyses
        ranking_method (str): Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending (bool): Sort direction (None uses default for each method)
    
    Returns:
        dict: Results dictionary containing analysis results for each cell type
    """
    # Load the dataframe

    if isinstance(marker, pd.DataFrame):
        df = marker.copy()
    elif isinstance(marker, str):
        df = pd.read_csv(marker)
    else:
        raise ValueError("marker must be either a pandas DataFrame or a string path to a CSV file")
    
    # If dataframe has only two columns, assume it's already processed
    if len(df.columns) == 2:
        print("Using input dataframe directly as it appears to be pre-processed (2 columns)")
    else:
        print("Processing input dataframe to get top markers")
        df = get_top_markers(df, n_genes=n_genes, ranking_method=ranking_method, ascending=ascending)
    
    # If celltype_column is not specified, use the first column
    if celltype_column is None:
        celltype_column = df.columns[0]
    
    # If gene_column_name is not specified, use the second column
    if gene_column_name is None:
        gene_column_name = df.columns[1]
    
    
    # Choose the appropriate analysis function based on provider
    analysis_function = runCASSIA
    
    def analyze_cell_type(cell_type, marker_list):
        print(f"\nAnalyzing {cell_type}...")
        for attempt in range(max_retries + 1):
            try:
                result, conversation_history = runCASSIA(
                    model=model,
                    temperature=temperature,
                    marker_list=marker_list,
                    tissue=tissue,
                    species=species,
                    additional_info=additional_info,
                    provider=provider,
                    validator_involvement=validator_involvement
                )
                # Add the number of markers and marker list to the result
                result['num_markers'] = len(marker_list)
                result['marker_list'] = marker_list
                print(f"Analysis for {cell_type} completed.")
                return cell_type, result, conversation_history
            except Exception as exc:
                # Don't retry authentication errors
                if "401" in str(exc) or "API key" in str(exc) or "authentication" in str(exc).lower():
                    print(f'{cell_type} generated an exception: {exc}')
                    print(f'This appears to be an API authentication error. Please check your API key.')
                    raise exc
                
                # For other errors, retry if attempts remain
                if attempt < max_retries:
                    print(f'{cell_type} generated an exception: {exc}')
                    print(f'Retrying analysis for {cell_type} (attempt {attempt + 2}/{max_retries + 1})...')
                else:
                    print(f'{cell_type} failed after {max_retries + 1} attempts with error: {exc}')
                    raise exc

    results = {}
    
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_celltype = {executor.submit(analyze_cell_type, row[celltype_column], split_markers(row[gene_column_name])): row[celltype_column] for _, row in df.iterrows()}
        
        # Process completed tasks
        for future in as_completed(future_to_celltype):
            cell_type = future_to_celltype[future]
            try:
                cell_type, result, conversation_history = future.result()
                if result:
                    results[cell_type] = {
                        "analysis_result": result,
                        "conversation_history": conversation_history,
                        "iterations": result.get("iterations", 1)
                    }
            except Exception as exc:
                print(f'{cell_type} failed: {exc}')
    

    
    print(f"All analyses completed. Results saved to '{output_name}'.")
    
    # Prepare data for both CSV files

    full_data = []
    summary_data = []

    for true_cell_type, details in results.items():
        main_cell_type = safe_get(details, 'analysis_result', 'main_cell_type')
        sub_cell_types = ', '.join(safe_get(details, 'analysis_result', 'sub_cell_types') or [])
        possible_mixed_cell_types = ', '.join(safe_get(details, 'analysis_result', 'possible_mixed_cell_types') or [])
        marker_number = safe_get(details, 'analysis_result', 'num_markers')
        marker_list = ', '.join(safe_get(details, 'analysis_result', 'marker_list') or [])
        iterations = safe_get(details, 'analysis_result', 'iterations')
        
        # Process and clean conversation history for safe CSV storage
        raw_conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in safe_get(details, 'conversation_history') or []])
        conversation_history = clean_conversation_history(raw_conversation_history)  # Clean while preserving full content
        
        full_data.append([
            true_cell_type, 
            main_cell_type, 
            sub_cell_types, 
            possible_mixed_cell_types,
            marker_number, 
            marker_list,
            iterations,
            model,
            provider,
            tissue,
            species,
            additional_info or "N/A",
            conversation_history
        ])
        summary_data.append([
            true_cell_type, 
            main_cell_type, 
            sub_cell_types, 
            possible_mixed_cell_types, 
            marker_list,
            iterations,
            model,
            provider,
            tissue,
            species
        ])

    # Generate output filenames based on input JSON filename
    base_name = os.path.splitext(output_name)[0]
    full_csv_name = f"{base_name}_full.csv"
    summary_csv_name = f"{base_name}_summary.csv"

    # Make sure the output directory exists
    output_dir = os.path.dirname(full_csv_name)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Sort the data by True Cell Type with natural/numeric ordering
    # This ensures "cluster 1", "cluster 2", "cluster 10" (not "cluster 1", "cluster 10", "cluster 2")
    full_data.sort(key=lambda x: natural_sort_key(x[0]))  # Sort by first column (True Cell Type) numerically
    summary_data.sort(key=lambda x: natural_sort_key(x[0]))  # Sort by first column (True Cell Type) numerically

    # Write the full data CSV with updated headers
    write_csv(full_csv_name, 
              ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types', 
               'Possible Mixed Cell Types', 'Marker Number', 'Marker List', 'Iterations',
               'Model', 'Provider', 'Tissue', 'Species', 'Additional Info', 'Conversation History'],
              full_data)

    # Write the summary data CSV with updated headers
    write_csv(summary_csv_name,
              ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types',
               'Possible Mixed Cell Types', 'Marker List', 'Iterations', 'Model', 'Provider',
               'Tissue', 'Species'],
              summary_data)

    print(f"Two CSV files have been created:")
    print(f"1. {full_csv_name} (full data)")
    print(f"2. {summary_csv_name} (summary data)")



############### Scoring annotation #################

def prompt_creator_score(major_cluster_info, marker, annotation_history):
    prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and rate single-cell annotation results, focusing on their correctness and ability to capture the overall picture of the data. You will provide a score from 0 to 100 and justify your rating.

Here are the single-cell annotation results to evaluate:



<marker>
{marker}
</marker>

<Cluster Origin>
{major_cluster_info}
</Cluster Origin>

<annotation_history>
{annotation_history}
</annotation_history>

Carefully analyze these results, paying particular attention to the following aspects:
1. Correctness of the annotations
2. Balanced consideration of multiple markers rather than over-focusing on a specific one
3. Ability to capture the general picture of the cell populations

When evaluating, consider:
- Are the annotations scientifically accurate?
- Is there a good balance in the use of different markers?
- Does the annotation provide a comprehensive view of the cell types present?
- Are there any obvious misclassifications or oversights?
- Did it consider the rank of the marker? marker appear first is more important.

Provide your analysis in the following format:
1. Start with a <reasoning> tag, where you explain your evaluation of the annotation results. Discuss the strengths and weaknesses you've identified, referring to specific examples from the results where possible.
2. After your reasoning, use a <score> tag to provide a numerical score from 0 to 100, where 0 represents completely incorrect or unusable results, and 100 represents perfect annotation that captures all aspects of the data correctly.

Your response should look like this:

<reasoning>
[Your detailed analysis and justification here]
</reasoning>

<score>[Your numerical score between 0 and 100]</score>

Remember, the focus is on correctness and the ability to see the general picture, rather than the structure of the results. Be critical but fair in your assessment.
    """
    return prompt


def extract_score_and_reasoning(text):
    """
    Extract both score and reasoning from annotation text.
    
    Args:
        text (str): Text containing score and reasoning between XML-like tags
        
    Returns:
        tuple: (score, reasoning_text) where score is int or None and reasoning_text is str or None
        
    Example:
        >>> score, reasoning = extract_score_and_reasoning("<reasoning>Good analysis</reasoning><score>85</score>")
        >>> print(f"Score: {score}, Reasoning: {reasoning[:20]}...")
        Score: 85, Reasoning: Good analysis...
    """
    try:
        # Initialize results
        score = None
        reasoning = None
        
        # Extract score - try multiple patterns
        score_patterns = [
            r'<score>(\d+)</score>',  # Original format
            r'Score:\s*(\d+)',        # "Score: 85"
            r'score:\s*(\d+)',        # "score: 85"  
            r'(\d+)/100',             # "85/100"
            r'(\d+)\s*out\s*of\s*100', # "85 out of 100"
            r'rating.*?(\d+)',        # "rating of 85"
            r'(\d+)%'                 # "85%"
        ]
        
        for pattern in score_patterns:
            score_match = re.search(pattern, text, re.IGNORECASE)
            if score_match:
                score = int(score_match.group(1))
                break
        
        # Extract reasoning - try multiple patterns
        reasoning_patterns = [
            r'<reasoning>(.*?)</reasoning>',  # Original format
            r'Reasoning:\s*(.*?)(?=Score:|$)',  # "Reasoning: ..." until "Score:" or end
            r'reasoning:\s*(.*?)(?=score:|$)',  # lowercase version
            r'Analysis:\s*(.*?)(?=Score:|$)',   # "Analysis: ..."
            r'Evaluation:\s*(.*?)(?=Score:|$)' # "Evaluation: ..."
        ]
        
        for pattern in reasoning_patterns:
            reasoning_match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
            if reasoning_match:
                reasoning = reasoning_match.group(1).strip()
                break
        
        # If no specific reasoning found, use the entire text as reasoning
        if reasoning is None and text.strip():
            reasoning = text.strip()
        
        return score, reasoning
        
    except Exception as e:
        print(f"Error extracting data: {str(e)}")
        return None, None

def score_single_analysis(major_cluster_info, marker, annotation_history, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Score a single cell type annotation analysis.
    
    Args:
        major_cluster_info (str): Information about species and tissue
        marker (str): Comma-separated list of marker genes
        annotation_history (str): History of annotation conversation
        model (str): Model to use (e.g., "gpt-4" for OpenAI or "claude-3-5-sonnet-20241022" for Anthropic)
        provider (str): AI provider to use ('openai', 'anthropic', or 'openrouter')
        
    Returns:
        tuple: (score, reasoning) where score is int and reasoning is str
    """
    prompt = prompt_creator_score(major_cluster_info, marker, annotation_history)
    
    # Add explicit max_tokens to ensure responses aren't truncated
    response = call_llm(
        prompt=prompt, 
        provider=provider, 
        model=model, 
        max_tokens=4096  # Maximum tokens allowed for most models
    )
    
    score, reasoning = extract_score_and_reasoning(response)
    return score, reasoning



def process_single_row(row_data, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Process a single row of data.
    
    Args:
        row_data (tuple): (idx, row) containing index and row data
        model (str): Model to use
        provider (str): AI provider to use ('openai' or 'anthropic')
        
    Returns:
        tuple: (idx, score, reasoning)
    """
    idx, row = row_data
    
    try:
        major_cluster_info = f"{row['Species']} {row['Tissue']}"
        
        # Handle both 'Marker List' and 'Marker.List' column names
        marker_column_options = ['Marker List', 'Marker.List', 'marker_list', 'Marker_List']
        marker = None
        for col in marker_column_options:
            if col in row:
                marker = row[col]
                break
        if marker is None:
            raise KeyError(f"Could not find marker column. Available columns: {list(row.index)}")
        
        # Handle both 'Conversation History' and 'Conversation.History' column names
        history_column_options = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History']
        annotation_history = None
        for col in history_column_options:
            if col in row:
                annotation_history = row[col]
                break
        if annotation_history is None:
            raise KeyError(f"Could not find conversation history column. Available columns: {list(row.index)}")
        
        # Try up to 3 times for a valid score if we get None
        score, reasoning = None, None
        max_retries_for_none = 3
        retry_count = 0
        
        while score is None and retry_count < max_retries_for_none:
            if retry_count > 0:
                print(f"Retry {retry_count}/{max_retries_for_none} for row {idx + 1} due to None score")
            
            score, reasoning = score_single_analysis(
                major_cluster_info, 
                marker, 
                annotation_history,
                model=model,
                provider=provider
            )
            
            if score is not None:
                break
                
            retry_count += 1

        print(f"Processed row {idx + 1}: Score = {score}")
        return (idx, score, reasoning)
        
    except Exception as e:
        print(f"Error processing row {idx + 1}: {str(e)}")
        return (idx, None, f"Error: {str(e)}")


def score_annotation_batch(results_file_path, output_file_path=None, max_workers=4, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Process and score all rows in a results CSV file in parallel.
    
    Args:
        results_file_path (str): Path to the results CSV file
        output_file_path (str, optional): Path to save the updated results
        max_workers (int): Maximum number of parallel threads
        model (str): Model to use
        provider (str): AI provider to use ('openai' or 'anthropic')
        
    Returns:
        pd.DataFrame: Original results with added score and reasoning columns
    """
    # Read results file
    results = pd.read_csv(results_file_path)
    
    # Initialize new columns if they don't exist
    if 'Score' not in results.columns:
        results['Score'] = None
    if 'Scoring_Reasoning' not in results.columns:
        results['Scoring_Reasoning'] = None
    
    # Create a list of unscored rows to process
    rows_to_process = [
        (idx, row) for idx, row in results.iterrows() 
        if pd.isna(row['Score'])
    ]
    
    if not rows_to_process:
        print("All rows already scored!")
        return results
    
    # Set up a lock for DataFrame updates
    df_lock = threading.Lock()
    
    # Process rows in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_row = {
            executor.submit(
                process_single_row, 
                row_data,
                model=model,
                provider=provider
            ): row_data[0] 
            for row_data in rows_to_process
        }
        
        # Process completed jobs
        for future in as_completed(future_to_row):
            idx, score, reasoning = future.result()
            
            # Safely update DataFrame
            with df_lock:
                results.loc[idx, 'Score'] = score
                results.loc[idx, 'Scoring_Reasoning'] = reasoning
                
                # Save intermediate results
                if output_file_path is None:
                    output_file_path = results_file_path.replace('.csv', '_scored.csv')
                results.to_csv(output_file_path, index=False)
    
    return results

def runCASSIA_score_batch(input_file, output_file=None, max_workers=4, model="deepseek/deepseek-chat-v3-0324", provider="openrouter", max_retries=1):
    """
    Run scoring with progress updates.
    
    Args:
        input_file (str): Path to input CSV file (with or without .csv extension)
        output_file (str, optional): Path to output CSV file (with or without .csv extension)
        max_workers (int): Maximum number of parallel workers
        model (str): Model to use
        provider (str): AI provider to use ('openai' or 'anthropic')
        max_retries (int): Maximum number of retries for failed analyses
        
    Returns:
        pd.DataFrame: Results DataFrame with scores
    """
    # Add .csv extension if not present
    if not input_file.lower().endswith('.csv'):
        input_file = input_file + '.csv'
    
    if output_file and not output_file.lower().endswith('.csv'):
        output_file = output_file + '.csv'
    
    print(f"Starting scoring process with {max_workers} workers using {provider} ({model})...")
    
    try:
        # Read the input file
        results = pd.read_csv(input_file)
        
        # Initialize new columns if they don't exist
        if 'Score' not in results.columns:
            results['Score'] = None
        if 'Scoring_Reasoning' not in results.columns:
            results['Scoring_Reasoning'] = None
        
        # Create a list of unscored rows to process
        rows_to_process = [
            (idx, row) for idx, row in results.iterrows() 
            if pd.isna(row['Score'])
        ]
        
        if not rows_to_process:
            print("All rows already scored!")
            return results
        
        # Set up a lock for DataFrame updates
        df_lock = threading.Lock()
        
        # Define a function that includes retry logic
        def process_with_retry(row_data):
            idx, row = row_data
            for attempt in range(max_retries + 1):
                try:
                    return process_single_row(row_data, model=model, provider=provider)
                except Exception as exc:
                    # Don't retry authentication errors
                    if "401" in str(exc) or "API key" in str(exc) or "authentication" in str(exc).lower():
                        print(f'⚠️  Row {idx + 1} API ERROR: {exc}')
                        print(f'⚠️  This appears to be an API authentication error. Please check your API key.')
                        # Return error info instead of raising
                        return idx, None, f"API error: {str(exc)}"
                    
                    # For other errors, retry if attempts remain
                    if attempt < max_retries:
                        print(f'⚠️  Row {idx + 1} ERROR: {exc}')
                        print(f'🔄 RETRYING row {idx + 1} (exception retry {attempt + 1}/{max_retries})...')
                    else:
                        print(f'❌ Row {idx + 1} FAILED after {max_retries + 1} attempts with error: {exc}')
                        # Return error info instead of raising
                        return idx, None, f"Failed after {max_retries + 1} attempts: {str(exc)}"
        
        # Process rows in parallel
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_row = {
                executor.submit(process_with_retry, row_data): row_data[0] 
                for row_data in rows_to_process
            }
            
            # Process completed jobs
            for future in as_completed(future_to_row):
                try:
                    idx, score, reasoning = future.result()
                    
                    # Safely update DataFrame
                    with df_lock:
                        results.loc[idx, 'Score'] = score
                        results.loc[idx, 'Scoring_Reasoning'] = reasoning
                        
                        # Save intermediate results if output file is specified
                        if output_file:
                            results.to_csv(output_file, index=False)
                        else:
                            output_file = input_file.replace('.csv', '_scored.csv')
                            results.to_csv(output_file, index=False)
                except Exception as exc:
                    print(f"Failed to process row: {exc}")
        
        # Print summary statistics
        total_rows = len(results)
        scored_rows = results['Score'].notna().sum()
        print(f"\nScoring completed!")
        print(f"\nSummary:")
        print(f"Total rows: {total_rows}")
        print(f"Successfully scored: {scored_rows}")
        print(f"Failed/Skipped: {total_rows - scored_rows}")
        
        
    except Exception as e:
        print(f"Error in runCASSIA_score_batch: {str(e)}")
        raise




################### Validator plus ########################

def runCASSIA_annotationboost_additional_task(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="claude-3-5-sonnet-20241022",
    additional_task="",
    provider=None,
    temperature=0,
    conversation_history_mode="final",
    report_style="per_iteration"
):
    """
    Generate a detailed HTML report for cell type analysis of a specific cluster.
    
    Args:
        full_result_path (str): Path to the full results CSV file
        marker_path (str): Path to the marker genes CSV file
        cluster_name (str): Name of the cluster to analyze
        major_cluster_info (str): General information about the dataset (e.g., "Human PBMC")
        output_name (str): Name of the output HTML file
        num_iterations (int): Number of iterations for marker analysis (default=5)
        model (str): Model to use for analysis (default="claude-3-5-sonnet-20241022")
        additional_task (str): Additional analysis task to perform
        provider (str): AI provider to use (default=None, will be inferred from model name)
        temperature (float): Sampling temperature (0-1)
        conversation_history_mode (str): Mode for extracting conversation history ("full", "final", or "none")
        report_style (str): Style of report generation ("per_iteration" or "total_summary")
        
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    # Import here to avoid circular imports
    try:
        # Try normal import first
        from annotation_boost import runCASSIA_annotationboost_additional_task as run_annotationboost_additional_task
    except ImportError:
        try:
            # Try relative import as fallback
            try:
                from .annotation_boost import runCASSIA_annotationboost_additional_task as run_annotationboost_additional_task
            except ImportError:
                from annotation_boost import runCASSIA_annotationboost_additional_task as run_annotationboost_additional_task
        except ImportError:
            raise ImportError("Could not import annotation_boost module.")
    
    # Determine provider based on model name if not provided
    if provider is None:
        provider = "anthropic"
        if "gpt" in model.lower():
            provider = "openai"
        elif not model.startswith("claude"):
            provider = "openrouter"
    
    return run_annotationboost_additional_task(
        full_result_path=full_result_path,
        marker=marker,
        cluster_name=cluster_name,
        major_cluster_info=major_cluster_info,
        output_name=output_name,
        num_iterations=num_iterations,
        model=model,
        provider=provider,
        additional_task=additional_task,
        temperature=temperature,
        conversation_history_mode=conversation_history_mode,
        report_style=report_style
    )






def runCASSIA_annotationboost(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="google/gemini-2.5-flash-preview",
    provider="openrouter",
    temperature=0,
    conversation_history_mode="final",
    report_style="per_iteration"
):
    """
    Wrapper function to generate cell type analysis report using either OpenAI or Anthropic models.
    
    Args:
        full_result_path (str): Path to the full results CSV file
        marker (str): Path to the marker genes CSV file
        cluster_name (str): Name of the cluster to analyze
        major_cluster_info (str): General information about the dataset (e.g., "Human PBMC")
        output_name (str): Name of the output HTML file
        num_iterations (int): Number of iterations for marker analysis (default=5)
        model (str): Model to use for analysis 
            - OpenAI options: "gpt-4", "gpt-3.5-turbo", etc.
            - Anthropic options: "claude-3-opus-20240229", "claude-3-sonnet-20240229", etc.
        provider (str): AI provider to use ('openai' or 'anthropic' or 'openrouter')
        temperature (float): Sampling temperature (0-1)
        conversation_history_mode (str): Mode for extracting conversation history ("full", "final", or "none")
        report_style (str): Style of report generation ("per_iteration" or "total_summary")
    
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    # Import here to avoid circular imports
    try:
        # Try normal import first
        from annotation_boost import runCASSIA_annotationboost as run_annotationboost
    except ImportError:
        try:
            # Try relative import as fallback
            try:
                from .annotation_boost import runCASSIA_annotationboost as run_annotationboost
            except ImportError:
                from annotation_boost import runCASSIA_annotationboost as run_annotationboost
        except ImportError:
            raise ImportError("Could not import annotation_boost module.")
    
    return run_annotationboost(
        full_result_path=full_result_path,
        marker=marker,
        cluster_name=cluster_name,
        major_cluster_info=major_cluster_info,
        output_name=output_name,
        num_iterations=num_iterations,
        model=model,
        provider=provider,
        temperature=temperature,
        conversation_history_mode=conversation_history_mode,
        report_style=report_style
    )



def generate_html_report(analysis_text):
    # Split the text into sections based on agents
    sections = analysis_text.split(" | ")
    
    # HTML template with CSS styling - note the double curly braces for CSS
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ 
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif; 
                max-width: 1200px; 
                margin: 0 auto; 
                padding: 20px; 
                background-color: #f0f2f5;
                line-height: 1.6;
            }}
            .container {{ 
                background-color: white; 
                padding: 40px; 
                border-radius: 16px; 
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }}
            .agent-section {{ 
                margin-bottom: 35px; 
                padding: 25px; 
                border-radius: 12px; 
                transition: all 0.3s ease;
            }}
            .agent-section:hover {{
                transform: translateY(-2px);
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }}
            .final-annotation {{ 
                background-color: #f0f7ff; 
                border-left: 5px solid #2196f3; 
            }}
            .validator {{ 
                background-color: #f0fdf4; 
                border-left: 5px solid #22c55e; 
            }}
            .formatting {{ 
                background: linear-gradient(145deg, #fff7ed, #ffe4c4);
                border-left: 5px solid #f97316; 
                box-shadow: 0 4px 15px rgba(249, 115, 22, 0.1);
            }}
            h2 {{ 
                color: #1a2b3c; 
                margin-top: 0; 
                font-size: 1.5rem;
                font-weight: 600;
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            ul {{ 
                margin: 15px 0; 
                padding-left: 20px; 
            }}
            pre {{ 
                background-color: #f8fafc; 
                padding: 20px; 
                border-radius: 8px; 
                overflow-x: auto;
                font-family: 'Consolas', 'Monaco', monospace;
                font-size: 0.9rem;
                line-height: 1.5;
            }}
            .validation-result {{ 
                font-weight: 600; 
                color: #16a34a; 
                padding: 12px 20px;
                background-color: #dcfce7; 
                border-radius: 8px; 
                display: inline-block;
                margin: 10px 0;
            }}
            br {{ 
                margin-bottom: 8px; 
            }}
            p {{
                margin: 12px 0;
                color: #374151;
            }}
            .summary-content {{
                display: flex;
                flex-direction: column;
                gap: 24px;
            }}
            .summary-item {{
                display: flex;
                flex-direction: column;
                gap: 8px;
                background: rgba(255, 255, 255, 0.7);
                padding: 16px;
                border-radius: 12px;
                backdrop-filter: blur(8px);
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
            }}
            .summary-label {{
                font-weight: 600;
                color: #c2410c;
                font-size: 0.95rem;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }}
            .summary-value {{
                color: #1f2937;
                font-size: 1.1rem;
                padding: 8px 16px;
                background-color: rgba(255, 255, 255, 0.9);
                border-radius: 8px;
                display: inline-block;
                box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            }}
            .summary-list {{
                margin: 0;
                padding-left: 24px;
                list-style-type: none;
            }}
            .summary-list li {{
                color: #1f2937;
                padding: 8px 0;
                position: relative;
            }}
            .summary-list li:before {{
                content: "•";
                color: #f97316;
                font-weight: bold;
                position: absolute;
                left: -20px;
            }}
            .report-header {{
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }}
            
            .report-title {{
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }}
            
            .report-subtitle {{
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }}
            .scoring {{ 
                background: linear-gradient(145deg, #f0fdf4, #dcfce7);
                border-left: 5px solid #22c55e;
                box-shadow: 0 4px 15px rgba(34, 197, 94, 0.1);
            }}
            .scoring-content {{
                display: flex;
                flex-direction: column;
                gap: 16px;
                color: #1f2937;
                line-height: 1.8;
            }}
            .scoring-content br + br {{
                content: "";
                display: block;
                margin: 12px 0;
            }}
            .empty-list {{
                color: #6b7280;
                font-style: italic;
            }}
            .error-message {{
                color: #dc2626;
                padding: 12px;
                background-color: #fef2f2;
                border-radius: 6px;
                border-left: 4px solid #dc2626;
            }}
            .score-badge {{
                background: linear-gradient(135deg, #22c55e, #16a34a);
                color: white;
                padding: 8px 16px;
                border-radius: 12px;
                font-size: 1.5rem;
                font-weight: 700;
                display: inline-block;
                margin: 12px 0;
                box-shadow: 0 4px 12px rgba(34, 197, 94, 0.2);
                position: relative;
                top: -10px;
            }}
            .score-badge::before {{
                content: "Score:";
                font-size: 0.9rem;
                font-weight: 500;
                margin-right: 8px;
                opacity: 0.9;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <h1 class="report-title">CASSIA Analysis Report</h1>
                <p class="report-subtitle">Comprehensive Cell Type Analysis and Annotation</p>
            </div>
            {0}
        </div>
    </body>
    </html>
    """
    
    content = []
    
    # Process each section
    for section in sections:
        if section.startswith("Final Annotation Agent:"):
            annotation_content = section.replace("Final Annotation Agent:", "").strip()
            content.append("""
                <div class="agent-section final-annotation">
                    <h2>🔍 Final Annotation Analysis</h2>
                    {0}
                </div>
            """.format(annotation_content.replace('\n', '<br>')))
            
        elif section.startswith("Coupling Validator:"):
            validator_content = section.replace("Coupling Validator:", "").strip()
            validation_result = '<div class="validation-result">✅ VALIDATION PASSED</div>' if "VALIDATION PASSED" in validator_content else ""
            
            content.append("""
                <div class="agent-section validator">
                    <h2>✓ Validation Check</h2>
                    {0}
                    {1}
                </div>
            """.format(validation_result, validator_content.replace('\n', '<br>')))
            
        elif section.startswith("Formatting Agent:"):
            try:
                import json
                # Get the content after "Formatting Agent:"
                json_text = section.replace("Formatting Agent:", "").strip()
                
                # Since the JSON is consistently formatted with newlines,
                # we can find where it ends (the last '}' followed by a newline or end of string)
                json_end = json_text.rfind('}')
                if json_end != -1:
                    json_content = json_text[:json_end + 1]
                    data = json.loads(json_content)
                    
                    # Process the data...
                    main_cell_type = data.get('main_cell_type', 'Not specified')
                    sub_cell_types = data.get('sub_cell_types', [])
                    mixed_types = data.get('possible_mixed_cell_types', [])
                    num_markers = data.get('num_markers', 'Not specified')
                    
                    # Format the content...
                    formatted_content = f"""
                        <div class="summary-content">
                            <div class="summary-item">
                                <span class="summary-label">Main Cell Type:</span>
                                <span class="summary-value">{main_cell_type}</span>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Sub Cell Types:</span>
                                <ul class="summary-list">
                                    {"".join(f'<li>{item}</li>' for item in sub_cell_types) if sub_cell_types 
                                     else '<li class="empty-list">No sub cell types identified</li>'}
                                </ul>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Possible Mixed Cell Types:</span>
                                <ul class="summary-list">
                                    {"".join(f'<li>{item}</li>' for item in mixed_types) if mixed_types 
                                     else '<li class="empty-list">No mixed cell types identified</li>'}
                                </ul>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Number of Markers:</span>
                                <span class="summary-value">{num_markers}</span>
                            </div>
                        </div>
                    """
                    
                    content.append(f"""
                        <div class="agent-section formatting">
                            <h2>📋 Summary</h2>
                            {formatted_content}
                        </div>
                    """)
                else:
                    raise ValueError("Could not find JSON content")
                    
            except Exception as e:
                content.append(f"""
                    <div class="agent-section formatting">
                        <h2>📋 Summary</h2>
                        <p class="error-message">Error formatting data: {str(e)}</p>
                    </div>
                """)
        elif section.startswith("Scoring Agent:"):
            try:
                # Get the content after "Scoring Agent:"
                scoring_text = section.split("Scoring Agent:", 1)[1].strip()
                
                # Split the score from the main text
                main_text, score = scoring_text.rsplit("Score:", 1)
                score = score.strip()
                
                content.append(r"""
                    <div class="agent-section scoring">
                        <h2>🎯 Quality Assessment</h2>
                        <div class="score-badge">{0}</div>
                        <div class="scoring-content">
                            {1}
                        </div>
                    </div>
                """.format(score, main_text.replace('\n', '<br>')))
            except Exception as e:
                content.append(r"""
                    <div class="agent-section scoring">
                        <h2>🎯 Quality Assessment</h2>
                        <p class="error-message">Error formatting scoring data: {0}</p>
                    </div>
                """.format(str(e)))
    
    # Combine all sections
    final_html = html_template.format(''.join(content))
    return final_html



def process_single_report(text, score_reasoning, score):
    combined = (
        f"{text}\n"
        f" | Scoring Agent: {score_reasoning}\n"
        f"Score: {score}"
    )
    return generate_html_report(combined)


def generate_index_page(report_files):
    index_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ 
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif; 
                max-width: 1200px; 
                margin: 0 auto; 
                padding: 20px; 
                background-color: #f0f2f5;
                line-height: 1.6;
            }}
            .container {{ 
                background-color: white; 
                padding: 40px; 
                border-radius: 16px; 
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }}
            .report-list {{
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
                gap: 20px;
                padding: 20px 0;
            }}
            .report-link {{
                background: white;
                padding: 20px;
                border-radius: 12px;
                text-decoration: none;
                color: #1a2b3c;
                border: 1px solid #e5e7eb;
                transition: all 0.3s ease;
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            .report-link:hover {{
                transform: translateY(-2px);
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                border-color: #f97316;
            }}
            .report-icon {{
                font-size: 24px;
            }}
            .report-header {{
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }}
            .index-title {{
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }}
            .index-subtitle {{
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <h1 class="index-title">CASSIA Reports Summary</h1>
                <p class="index-subtitle">Select a report to view detailed analysis</p>
            </div>
            <div class="report-list">
                {0}
            </div>
        </div>
    </body>
    </html>
    """
    
    # Generate links for each report
    links = []
    for filename in sorted(report_files):
        display_name = filename.replace('report_', '').replace('.html', '')
        links.append(f'<a href="{filename}" class="report-link"><span class="report-icon">📊</span>{display_name}</a>')
    
    return index_template.format('\n'.join(links))

def runCASSIA_generate_score_report(csv_path, index_name="CASSIA_reports_summary"):
    """
    Generate HTML reports from a scored CSV file and create an index page.
    
    Args:
        csv_path (str): Path to the CSV file containing the score results
        index_name (str): Base name for the index file (without .html extension)
    """
    # Read the CSV file
    report = pd.read_csv(csv_path)
    report_files = []
    
    # Determine output folder (same folder as the CSV file)
    output_folder = os.path.dirname(csv_path)
    if not output_folder:
        output_folder = "."
    
    # Process each row
    for index, row in report.iterrows():
        # Get the first column value for the filename
        filename = str(row.iloc[0]).strip()
        filename = "".join(c for c in filename if c.isalnum() or c in (' ', '-', '_')).strip()
        
        # Handle both 'Conversation History' and 'Conversation.History' column names
        history_column_options = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History']
        text = None
        for col in history_column_options:
            if col in row:
                text = row[col]
                break
        if text is None:
            raise KeyError(f"Could not find conversation history column. Available columns: {list(row.index)}")
        
        # Handle both 'Scoring_Reasoning' and 'Scoring.Reasoning' column names
        reasoning_column_options = ['Scoring_Reasoning', 'Scoring.Reasoning', 'scoring_reasoning', 'Scoring_reasoning']
        score_reasoning = None
        for col in reasoning_column_options:
            if col in row:
                score_reasoning = row[col]
                break
        if score_reasoning is None:
            raise KeyError(f"Could not find scoring reasoning column. Available columns: {list(row.index)}")
        
        score = row["Score"]
        
        # Generate HTML for this row
        html_content = process_single_report(text, score_reasoning, score)
        
        # Save using the first column value as filename in the output folder
        output_path = os.path.join(output_folder, f"report_{filename}.html")
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        
        # Store just the filename for the index (not the full path)
        report_files.append(os.path.basename(output_path))
        print(f"Report saved to {output_path}")
    
    # Generate and save index page in the same folder
    index_html = generate_index_page(report_files)
    index_filename = os.path.join(output_folder, f"{os.path.basename(index_name)}.html")
    with open(index_filename, "w", encoding="utf-8") as f:
        f.write(index_html)
    print(f"Index page saved to {index_filename}")


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
    validator_involvement: str = "v1"
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
    """
    # Create a main folder based on tissue and species for organizing reports
    main_folder_name = f"CASSIA_{tissue}_{species}"
    main_folder_name = "".join(c for c in main_folder_name if c.isalnum() or c in (' ', '-', '_')).strip()
    main_folder_name = main_folder_name.replace(' ', '_')
    
    # Remove .csv extension if present
    if output_file_name.lower().endswith('.csv'):
        output_file_name = output_file_name[:-4] # Remove last 4 characters (.csv)

    # Add timestamp to prevent overwriting existing folders with the same name
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    main_folder_name = f"{main_folder_name}_{timestamp}"
    
    # Create the main folder if it doesn't exist
    if not os.path.exists(main_folder_name):
        os.makedirs(main_folder_name)
        print(f"Created main folder: {main_folder_name}")
        
    # Create organized subfolders according to user's specifications
    annotation_results_folder = os.path.join(main_folder_name, "01_annotation_results")  # All CSV files
    reports_folder = os.path.join(main_folder_name, "02_reports")  # All HTML reports except annotation boost
    boost_folder = os.path.join(main_folder_name, "03_boost_analysis")   # All annotation boost related results
    
    # Create all subfolders
    for folder in [annotation_results_folder, reports_folder, boost_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print(f"Created subfolder: {folder}")
    
    # Define derived file names with folder paths
    # All CSV files go to the annotation_results_folder
    raw_full_csv = os.path.join(annotation_results_folder, f"{output_file_name}_full.csv")
    raw_summary_csv = os.path.join(annotation_results_folder, f"{output_file_name}_summary.csv")
    raw_sorted_csv = os.path.join(annotation_results_folder, f"{output_file_name}_sorted_full.csv")
    score_file_name = os.path.join(annotation_results_folder, f"{output_file_name}_scored.csv")
    merged_annotation_file = os.path.join(annotation_results_folder, f"{output_file_name}_merged.csv")
    
    # Reports go to the reports_folder - ALL HTML reports should be in this folder
    report_base_name = os.path.join(reports_folder, f"{output_file_name}")
    
    # First annotation output is in the current directory but will be moved later
    annotation_output = output_file_name

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
        summary_csv = raw_summary_csv
        
        # Import the merge_annotations function dynamically
        try:
            try:
                from .merging_annotation import merge_annotations_all
            except ImportError:
                from merging_annotation import merge_annotations_all
            
            # Sort the CSV file by True Cell Type before merging to ensure consistent order
            print("Sorting CSV by True Cell Type before merging...")
            df = pd.read_csv(raw_full_csv)
            df = df.sort_values(by=['True Cell Type'])
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
    # Run scoring
    runCASSIA_score_batch(
        input_file=raw_full_csv,
        output_file=score_file_name,
        max_workers=max_workers,
        model=score_model,
        provider=score_provider,
        max_retries=max_retries
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
            # Merge on 'True Cell Type' to add merged annotation columns
            if 'True Cell Type' in merged_df.columns:
                # Keep only the merged columns (not duplicating existing ones)
                merge_columns = [col for col in merged_df.columns if col not in final_df.columns or col == 'True Cell Type']
                final_df = final_df.merge(merged_df[merge_columns], on='True Cell Type', how='left')
        
        # Sort the final results by True Cell Type
        final_df = final_df.sort_values(by=['True Cell Type'])
        
        # Save the final combined results
        final_combined_file = os.path.join(annotation_results_folder, f"{output_file_name}_FINAL_RESULTS.csv")
        final_df.to_csv(final_combined_file, index=False)
        print(f"✓ Final combined results saved to {final_combined_file}")
        
    except Exception as e:
        print(f"Warning: Could not create final combined results: {str(e)}")
        final_combined_file = score_file_name  # Fallback to scored file

    print("\n=== Generating main reports ===")
    # Process reports - ensure they go to reports_folder
    runCASSIA_generate_score_report(
        csv_path=score_file_name,
        index_name=report_base_name  # This will create reports in the reports_folder
    )
    
    # Move any HTML files from annotation_results_folder to reports_folder
    for file in os.listdir(annotation_results_folder):
        if file.endswith('.html'):
            src_path = os.path.join(annotation_results_folder, file)
            dst_path = os.path.join(reports_folder, file)
            try:
                shutil.copy2(src_path, dst_path)
                os.remove(src_path)  # Remove from original location after copying
                print(f"Moved HTML report {file} to reports folder")
            except Exception as e:
                print(f"Error moving HTML file {file}: {str(e)}")
    
    print("✓ Main reports generated")

    print("\n=== Analyzing low-scoring clusters ===")
    # Handle low-scoring clusters
    df = pd.read_csv(score_file_name)
    low_score_clusters = df[df['Score'] < score_threshold]['True Cell Type'].tolist()

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
            cluster_output_name = os.path.join(cluster_boost_folder, f"{output_file_name}_{sanitized_cluster_name}_boosted")
            
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
    
    print("\n=== Organizing intermediate files ===")
    # Create intermediate files folder
    intermediate_folder = os.path.join(annotation_results_folder, "intermediate_files")
    if not os.path.exists(intermediate_folder):
        os.makedirs(intermediate_folder)
    
    # List of intermediate files to move
    intermediate_files = [
        raw_full_csv,
        raw_summary_csv, 
        raw_sorted_csv,
        score_file_name,
        merged_annotation_file
    ]
    
    # Move intermediate files to intermediate folder
    for file_path in intermediate_files:
        if os.path.exists(file_path):
            try:
                filename = os.path.basename(file_path)
                destination = os.path.join(intermediate_folder, filename)
                shutil.move(file_path, destination)
                print(f"Moved {filename} to intermediate_files folder")
            except Exception as e:
                print(f"Warning: Could not move {filename}: {str(e)}")
    
    print("✓ Intermediate files organized")
    
    # Try to clean up the original files in the root directory
    try:
        for file_to_remove in [original_full_csv, original_summary_csv, annotation_output + "_sorted_full.csv"]:
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)
                print(f"Removed original file: {file_to_remove}")
    except Exception as e:
        print(f"Warning: Could not remove some temporary files: {str(e)}")
    
    print("\n=== Cell type analysis pipeline completed ===")
    print(f"All results have been organized in the '{main_folder_name}' folder:")
    print(f"  📊 MAIN RESULTS: {final_combined_file}")
    print(f"  📁 HTML Reports: {reports_folder}")
    print(f"  🔍 Annotation Boost Results: {boost_folder}")
    print(f"  📂 Intermediate Files: {intermediate_folder}")
    print(f"\n✅ Your final results are in: {os.path.basename(final_combined_file)}")


def loadmarker(marker_type="processed"):
    """
    Load built-in marker files.
    
    Args:
        marker_type (str): Type of markers to load. Options:
            - "processed": For processed marker data
            - "unprocessed": For raw unprocessed marker data
            - "subcluster_results": For subcluster analysis results
    
    Returns:
        pandas.DataFrame: Marker data
    
    Raises:
        ValueError: If marker_type is not recognized
    """
    marker_files = {
        "processed": "processed.csv",
        "unprocessed": "unprocessed.csv",
        "subcluster_results": "subcluster_results.csv"
    }
    
    if marker_type not in marker_files:
        raise ValueError(f"Unknown marker type: {marker_type}. Available types: {list(marker_files.keys())}")
    
    filename = marker_files[marker_type]
    
    try:
        # Using importlib.resources for Python 3.7+
        with resources.path('CASSIA.data', filename) as file_path:
            return pd.read_csv(file_path)
    except Exception as e:
        raise Exception(f"Error loading marker file: {str(e)}")

def list_available_markers():
    """List all available built-in marker sets."""
    try:
        with resources.path('CASSIA.data', '') as data_path:
            marker_files = [f for f in os.listdir(data_path) if f.endswith('.csv')]
        return [f.replace('.csv', '') for f in marker_files]
    except Exception as e:
        raise Exception(f"Error listing marker files: {str(e)}")