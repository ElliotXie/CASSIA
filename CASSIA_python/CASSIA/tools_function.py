import pandas as pd
import json
import re
import csv
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from openai import OpenAI
from .main_function_code import *
import requests
import threading
import numpy as np
from importlib import resources
import datetime

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
        provider (str): The provider to set the key for ('openai' or 'anthropic')
    
    Raises:
        ValueError: If provider is not 'openai' or 'anthropic'
    """
    if provider.lower() == "openai":
        os.environ["OPENAI_API_KEY"] = api_key
    elif provider.lower() == "anthropic":
        os.environ["ANTHROPIC_API_KEY"] = api_key
    elif provider.lower() == "openrouter":
        os.environ["OPENROUTER_API_KEY"] = api_key
    else:
        raise ValueError("Provider must be either 'openai' or 'anthropic' or 'openrouter'")
    
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


def get_top_markers(df, n_genes=10, format_type=None):
    """
    Get top markers from either Seurat or Scanpy differential expression results.
    
    Args:
        df: Either a pandas DataFrame (Seurat format) or dictionary (Scanpy format)
        n_genes: Number of top genes to return per cluster
        format_type: Either 'seurat', 'scanpy', or None (auto-detect)
    
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
            genes = df['names'][cluster][:n_genes]
            logfc = df['logfoldchanges'][cluster][:n_genes]
            pvals_adj = df['pvals_adj'][cluster][:n_genes]
            pcts = df['pcts'][cluster][:n_genes]  # Get percentage information
            
            # Filter for significant upregulated genes with PCT threshold
            mask = (pvals_adj < 0.05) & (logfc > 0.25) & (pcts >= 0.1)  # Add PCT filter
            valid_genes = genes[mask]
            
            if len(valid_genes) > 0:
                # Join genes with commas
                markers = ','.join(valid_genes)
                top_markers.append({
                    'cluster': cluster,
                    'markers': markers
                })
        
        return pd.DataFrame(top_markers)
    
    else:  # Seurat format
        # Filter by adjusted p-value, positive log2FC, and PCT
        df_filtered = df[
            (df['p_val_adj'] < 0.05) & 
            (df['avg_log2FC'] > 0.25) &
            ((df['pct.1'] >=0.1) | (df['pct.2'] >=0.1))  # Add PCT filter
        ].copy()
        
        # Sort within each cluster by avg_log2FC and get top n genes
        top_markers = []
        
        for cluster in df_filtered['cluster'].unique():
            cluster_data = df_filtered[df_filtered['cluster'] == cluster]
            # Sort by avg_log2FC in descending order and take top n
            top_n = (cluster_data
                    .sort_values('avg_log2FC', ascending=False)
                    .head(n_genes))
            
            # Verify sorting worked correctly
            if not top_n.empty and not top_n['avg_log2FC'].is_monotonic_decreasing:
                print(f"Warning: Sorting issue detected for cluster {cluster}")
            
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
    for key in keys:
        try:
            dict_obj = dict_obj[key]
        except (KeyError, TypeError):
            return "N/A"
    return dict_obj


def write_csv(filename, headers, row_data):
    with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(row_data)


def runCASSIA(model="gpt-4o", temperature=0, marker_list=None, tissue="lung", species="human", additional_info=None, provider="openai"):
    """
    Wrapper function to run cell type analysis using either OpenAI or Anthropic's Claude
    
    Args:
        model (str): Model name to use
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for the analysis
        provider (str): AI provider to use ('openai' or 'anthropic')
    
    Returns:
        tuple: (analysis_result, conversation_history)
    """
    if provider.lower() == "openai":
        return run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info)
    elif provider.lower() == "anthropic":
        return run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info)
    elif provider.lower() == "openrouter":
        return run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info)
    else:
        raise ValueError("Provider must be either 'openai' or 'anthropic' or 'openrouter'")





def runCASSIA_batch(marker, output_name="cell_type_analysis_results.json", n_genes=50, model="gpt-4o", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column=None, gene_column_name=None, max_workers=10, provider="openai", max_retries=1):
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
        df = get_top_markers(df, n_genes=n_genes)
    
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
                result, conversation_history = analysis_function(
                    model=model,
                    temperature=temperature,
                    marker_list=marker_list,
                    tissue=tissue,
                    species=species,
                    additional_info=additional_info,
                    provider=provider
                )
                # Add the number of markers and marker list to the result
                result['num_markers'] = len(marker_list)
                result['marker_list'] = marker_list
                print(f"Analysis for {cell_type} completed.\n")
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
    
    # Function to safely get nested dictionary values
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
        
        conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in safe_get(details, 'conversation_history') or []])
        
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

    # Sort the data by True Cell Type to ensure consistent ordering
    full_data.sort(key=lambda x: x[0])  # Sort by first column (True Cell Type)
    summary_data.sort(key=lambda x: x[0])  # Sort by first column (True Cell Type)

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
    
    return results





def runCASSIA_batch_n_times(n, marker, output_name="cell_type_analysis_results", model="gpt-4o", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column=None, gene_column_name=None, max_workers=10, batch_max_workers=5, provider="openai", max_retries=1):
    def single_batch_run(i):
        output_json_name = f"{output_name}_{i}.json"
        print(f"Starting batch run {i+1}/{n}")
        start_time = time.time()
        result = runCASSIA_batch(
            marker=marker,
            output_name=output_json_name,
            model=model,
            temperature=temperature,
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            celltype_column=celltype_column,
            gene_column_name=gene_column_name,
            max_workers=max_workers,
            provider=provider,
            max_retries=max_retries
        )
        end_time = time.time()
        print(f"Finished batch run {i+1}/{n} in {end_time - start_time:.2f} seconds")
        return i, result, output_json_name

    all_results = []
    start_time = time.time()

    with ThreadPoolExecutor(max_workers=batch_max_workers) as executor:
        future_to_index = {executor.submit(single_batch_run, i): i for i in range(n)}
        
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                index, result, output_json_name = future.result()
                all_results.append((index, result, output_json_name))
                print(f"Batch run {index+1}/{n} completed and saved to {output_json_name}")
            except Exception as exc:
                print(f'Batch run {index+1} generated an exception: {exc}')

    end_time = time.time()
    print(f"All {n} batch runs completed in {end_time - start_time:.2f} seconds")

    return None

    #return all_results



def run_single_analysis(args):
    index, tissue, species, additional_info, temperature, marker_list, model, provider = args
    print(f"Starting analysis {index+1}")
    start_time = time.time()
    try:
        result = runCASSIA(
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            temperature=temperature,
            marker_list=marker_list,
            model=model,
            provider=provider
        )
        end_time = time.time()
        print(f"Finished analysis {index+1} in {end_time - start_time:.2f} seconds")
        return index, result
    except Exception as e:
        print(f"Error in analysis {index+1}: {str(e)}")
        return index, None

def runCASSIA_n_times(n, tissue, species, additional_info, temperature, marker_list, model, max_workers=10, provider="openai"):
    print(f"Starting {n} parallel analyses")
    start_time = time.time()
    
    results = {}
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks with the provider parameter
        future_to_index = {
            executor.submit(
                run_single_analysis, 
                (i, tissue, species, additional_info, temperature, marker_list, model, provider)
            ): i for i in range(n)
        }
        
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                index, result = future.result()
                if result:
                    results[index] = result
            except Exception as exc:
                print(f'Analysis {index+1} generated an exception: {exc}')
    
    end_time = time.time()
    print(f"All analyses completed in {end_time - start_time:.2f} seconds")
    return results







def parse_results_to_dict(result_string):

    
    # Use regex to find all results, supporting both () and [] formats
    pattern = r"result(\d+):[\(\[]([^\]\)]+)[\)\]]"
    matches = re.findall(pattern, result_string)
    
    
    # Parse each result
    parsed_results = {}
    for match in matches:
        try:
            result_num, cell_types = match
            # Split cell types, handling potential commas within cell type names
            cell_type_list = re.split(r',\s*(?=[^,]*(?:,|$))', cell_types)
            
            # Strip whitespace and remove any remaining quotes
            cell_type_list = [ct.strip().strip("'\"") for ct in cell_type_list]
            
            # Ensure we have at least two cell types (main and sub)
            while len(cell_type_list) < 2:
                cell_type_list.append("N/A")
            
            parsed_results[f"result{result_num}"] = tuple(cell_type_list[:2])
        except Exception as e:
            print(f"Error parsing match {match}: {str(e)}")
    
    
    return parsed_results



def extract_celltypes_from_llm(llm_response):
    # Extract the JSON part from the response
    json_match = re.search(r'```json\n(.*?)\n```', llm_response, re.DOTALL)
    
    if json_match:
        json_str = json_match.group(1)
        try:
            data = json.loads(json_str)
            final_results = data.get("final_results", [])
            mixed_celltypes = data.get("possible_mixed_celltypes", [])
            consensus_score = data.get("consensus_score", 0)
            general_celltype = final_results[0] if len(final_results) > 0 else "Not found"
            sub_celltype = final_results[1] if len(final_results) > 1 else "Not found"
            
            return general_celltype, sub_celltype, mixed_celltypes, consensus_score
        except json.JSONDecodeError:
            print("Error decoding JSON from LLM response")
    else:
        print("No JSON data found in the LLM response")
    
    return "Not found", "Not found", [], "Not found"


from collections import Counter

def consensus_similarity_flexible(results, main_weight=0.7, sub_weight=0.3):
    general_types = [result[0] for result in results.values()]
    sub_types = [result[1] for result in results.values()]
    
    consensus_general = Counter(general_types).most_common(1)[0][0]
    consensus_sub = Counter(sub_types).most_common(1)[0][0]
    
    total_score = 0
    for result in results.values():
        if result[0] == consensus_general:
            total_score += main_weight
        elif result[0] == consensus_sub:
            total_score += main_weight * sub_weight
        
        if result[1] == consensus_sub:
            total_score += sub_weight
        elif result[1] == consensus_general:
            total_score += sub_weight * main_weight
    
    similarity_score = total_score / (len(results) * (main_weight + sub_weight))
    
    return similarity_score, consensus_general, consensus_sub






def agent_unification_claude(prompt, system_prompt='''You are a careful professional biologist, specializing in single-cell RNA-seq analysis.You will be given a series results from a celltype annotator. 
your task is to unify all the celltypes name, so that same celltype have the same name. The final format the first letter for each word will be capital and other will be small case. Remove plural. Some words like stem and progenitor and immature means the same thing should be unified.
                  
An example below:
                  
Input format：      
result1:(immune cell, t cell),result2:(Immune cells,t cell),result3:(T cell, cd8+ t cell)
                  
Output format:
<results>result1:(Immune cell, T cell),result2:(Immune cell, T cell),result3:(T cell, Cd8+ t cell)</results>

Another example:
                      
Input format：      
result1:(Hematopoietic stem/progenitor cells (HSPCs), T cell progenitors),result2:(Hematopoietic Progenitor cells,t cell),result3:(Hematopoietic progenitor cells, T cell)
                  
Output format:
<results>result1:(Hematopoietic Progenitor Cells, T cell Progenitors),result2:(Hematopoietic Progenitor Cells,T cell),result3:(Hematopoietic Progenitor Cells, T cell)</results>
''', model="claude-3-5-sonnet-20241022", temperature=0):
    
    client = anthropic.Anthropic()
    try:
        message = client.messages.create(
            model=model,
            max_tokens=7000,
            temperature=temperature,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}]
        )
        response = message.content[0].text.strip()
        
        # Extract content between <results> tags
        import re
        match = re.search(r'<results>(.*?)</results>', response, re.DOTALL)
        if match:
            return match.group(1).strip()
        else:
            # If no tags found, return the original response
            return response.strip()
            
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


def agent_judgement_claude(prompt, system_prompt='''You are a careful professional biologist, specializing in single-cell RNA-seq analysis. You will be given a series of results from a cell type annotator.
Your task is to determine the consensus cell type. The first entry of each result is the general cell type and the second entry is the subtype. You should provide the final general cell type and the subtype. Considering all results, if you think there is very strong evidence of mixed cell types, please also list them. Please give your step-by-step reasoning and the final answer. We also want to know how robust our reuslts are, please give a consensus score ranging from 0 to 100 to show how similar the results are from different runs. $10,000 will be rewarded for the correct answer.
                           
Output in JSON format:
<json>{
"final_results": [
"General cell type here",
"Sub cell type here"
],
"possible_mixed_celltypes": [
"Mixed cell type 1 here",
"Mixed cell type 2 here"
]
"consensus_score": 0-100
}</json>
''', model="claude-3-5-sonnet-20241022", temperature=0):
    
    client = anthropic.Anthropic()
    try:
        message = client.messages.create(
            model=model,
            max_tokens=7000,
            temperature=temperature,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}]
        )
        return message.content[0].text.strip()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None






def agent_unification(prompt, system_prompt='''You are a careful professional biologist, specializing in single-cell RNA-seq analysis.You will be given a series results from a celltype annotator. 
your task is to unify all the celltypes name, so that same celltype have the same name. The final format the first letter for each word will be capital and other will be small case. Remove plural. Some words like stem and progenitor and immature means the same thing should be unified.
                  
An example below:
                  
Input format：      
result1:(immune cell, t cell),result2:(Immune cells,t cell),result3:(T cell, cd8+ t cell)
                  
Output format:
result1:(Immune cell, T cell),result2:(Immune cell, T cell),result3:(T cell, Cd8+ t cell)

Another example:
                      
Input format：      
result1:(Hematopoietic stem/progenitor cells (HSPCs), T cell progenitors),result2:(Hematopoietic Progenitor cells,t cell),result3:(Hematopoietic progenitor cells, T cell)
                  
Output format:
result1:(Hematopoietic Progenitor Cells, T cell Progenitors),result2:(Hematopoietic Progenitor Cells,T cell),result3:(Hematopoietic Progenitor Cells, T cell)             


''', model="gpt-4o", temperature=0):
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ],
            temperature=temperature,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None







def agent_judgement(prompt, system_prompt='''You are a careful professional biologist, specializing in single-cell RNA-seq analysis. You will be given a series of results from a cell type annotator.
Your task is to determine the consensus cell type. The first entry of each result is the general cell type and the second entry is the subtype. You should provide the final general cell type and the subtype. Considering all results, if you think there is very strong evidence of mixed cell types, please also list them. Please give your step-by-step reasoning and the final answer. We also want to know how robust our reuslts are, please give a consensus score ranging from 0 to 100 to show how similar the results are from different runs. $10,000 will be rewarded for the correct answer.
Output in JSON format:
{
"final_results": [
"General cell type here",
"Sub cell type here"
],
"possible_mixed_celltypes": [
"Mixed cell type 1 here",
"Mixed cell type 2 here"
],
"consensus_score": 0-100
}

'''
    , model="gpt-4o", temperature=0):
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ],
            temperature=temperature,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None





def agent_unification_deplural(prompt, system_prompt='''remove the plural for celltype name, keep the original input format.
''', model="gpt-4o", temperature=0):
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ],
            temperature=temperature,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None


import re
import requests

def get_cell_type_info(cell_type_name, ontology="CL"):
    # Check if the cell type name contains "mixed" (case-insensitive)
    if "mixed" in cell_type_name.lower():
        return "mixed cell population", "mixed cell population"

    base_url = "https://www.ebi.ac.uk/ols/api/search"
    params = {
        "q": cell_type_name,
        "ontology": ontology,
        "rows": 1
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if 'response' in data and 'docs' in data['response'] and data['response']['docs']:
            first_doc = data['response']['docs'][0]
            obo_id = first_doc.get('obo_id')
            label = first_doc.get('label')
            return obo_id, label
        else:
            return None, None
    
    except requests.RequestException:
        return None, None




def standardize_cell_types(input_string):
    # Remove all hyphens from the input string
    input_string = input_string.replace("-", " ")
    
    # Parse the input string into a list of tuples
    results = re.findall(r"result\d+:\('([^']+)', '([^']+)'\)", input_string)
    
    standardized_results = []
    for i, (general_type, specific_type) in enumerate(results, 1):
        # Search for standardized names
        _, general_label = get_cell_type_info(general_type)
        _, specific_label = get_cell_type_info(specific_type)
        
        # Use original names if no standardized names found
        general_label = general_label or general_type
        specific_label = specific_label or specific_type
        
        standardized_results.append(f"result{i}:('{general_label}', '{specific_label}')")
    
    return ",".join(standardized_results)


import pandas as pd
import glob
from collections import defaultdict


def organize_batch_results(marker, file_pattern, celltype_column=None):
    # Read marker data
    if isinstance(marker, pd.DataFrame):
        df = marker.copy()
    elif isinstance(marker, str):
        df = pd.read_csv(marker)
    else:
        raise ValueError("marker must be either a pandas DataFrame or a string path to a CSV file")

    # Only process with get_top_markers if more than 2 columns
    if len(df.columns) > 2:
        marker = get_top_markers(df, n_genes=50)
    else:
        marker = df  # Use the DataFrame directly if it has 2 or fewer columns
        
    # If celltype_column is not provided, use the first column
    if celltype_column is None:
        celltype_column = marker.columns[0]
    
    marker_celltype = marker[celltype_column]

    # Use glob to find all matching files
    file_list = sorted(glob.glob(file_pattern))

    # Initialize a defaultdict to store results for each cell type
    results = defaultdict(list)

    # Loop through each file (round of results)
    for file in file_list:
        df = pd.read_csv(file)
        
        # Loop through each cell type
        for celltype in marker_celltype:
            # Find the row for the current cell type
            row = df[df['True Cell Type'] == celltype]
            
            if not row.empty:
                # Extract the predicted general cell type (second column)
                predicted_general = row.iloc[0, 1]
                
                # Extract the first predicted subtype (first element in the third column)
                predicted_subtypes = row.iloc[0, 2]
                first_subtype = predicted_subtypes.split(',')[0].strip() if pd.notna(predicted_subtypes) else 'N/A'
                
                # Append the results as a tuple
                results[celltype].append((predicted_general, first_subtype))

    # Convert the defaultdict to a regular dict for easier handling
    organized_results = dict(results)
    
    return organized_results

# Example usage:
# organized_results = organize_batch_results(
#     marker_file_path="path/to/marker/file.csv",
#     file_pattern="batch_run_results_*.csv",
#     celltype_column="jackcvs"  # Optional: if not provided, will use the first column
# )


def process_cell_type_variance_analysis_batch(results, model="gpt-4o", temperature=0, main_weight=0.5, sub_weight=0.5):

    print("Starting the process of cell type batch variance analysis...")
    # Extract and format results
    results_unification_llm = agent_unification(results)


    results_depluar = agent_unification_deplural(results)


    result_unified_oncology = standardize_cell_types(results_depluar)


    # Consensus judgment
    result_consensus_from_llm = agent_judgement(
        prompt=results_unification_llm,
        model=model,
        temperature=temperature
    )


    result_consensus_from_oncology = agent_judgement(
        prompt=result_unified_oncology,
        model=model,
        temperature=temperature
    )


    # Extract consensus celltypes
    general_celltype, sub_celltype, mixed_types, llm_generated_consensus_score_llm = extract_celltypes_from_llm(result_consensus_from_llm)

    print(llm_generated_consensus_score_llm)

    general_celltype_oncology, sub_celltype_oncology, mixed_types_oncology, llm_generated_consensus_score_oncology = extract_celltypes_from_llm(result_consensus_from_oncology)
    
    print(f"General celltype: {general_celltype}")
    print(f"Sub celltype: {sub_celltype}")
    print(f"Mixed types: {mixed_types}")
    print(f"General celltype oncology: {general_celltype_oncology}")
    print(f"Sub celltype oncology: {sub_celltype_oncology}")
    print(f"Mixed types oncology: {mixed_types_oncology}")

    # Calculate similarity score
    parsed_results_oncology = parse_results_to_dict(result_unified_oncology)
    parsed_results_llm = parse_results_to_dict(results_unification_llm)


    consensus_score_oncology, consensus_1_oncology, consensus_2_oncology = consensus_similarity_flexible(parsed_results_oncology, main_weight=main_weight, sub_weight=sub_weight)
    consensus_score_llm, consensus_1_llm, consensus_2_llm = consensus_similarity_flexible(parsed_results_llm, main_weight=main_weight, sub_weight=sub_weight)

    print(f"Consensus score (Oncology): {consensus_score_oncology}")
    print(f"Consensus score (LLM): {consensus_score_llm}")
    print(f"LLM generated consensus score llm: {llm_generated_consensus_score_llm}")
    print(f"LLM generated consensus score oncology: {llm_generated_consensus_score_oncology}")

    return {
        'general_celltype_llm': general_celltype,
        'sub_celltype_llm': sub_celltype,
        'mixed_celltypes_llm': mixed_types,
        'general_celltype_oncology': general_celltype_oncology,
        'sub_celltype_oncology': sub_celltype_oncology,
        'mixed_types_oncology': mixed_types_oncology,
        'consensus_score_llm': consensus_score_llm,
        'consensus_score_oncology': consensus_score_oncology,
        'llm_generated_consensus_score_llm': llm_generated_consensus_score_llm,
        'llm_generated_consensus_score_oncology': llm_generated_consensus_score_oncology,
        'count_consensus_1_llm': consensus_1_llm,
        'count_consensus_2_llm': consensus_2_llm,
        'count_consensus_1_oncology': consensus_1_oncology,
        'count_consensus_2_oncology': consensus_2_oncology,
        'unified_results_llm': results_unification_llm,
        'unified_results_oncology': result_unified_oncology,
        'result_consensus_from_llm': result_consensus_from_llm,
        'result_consensus_from_oncology': result_consensus_from_oncology
    }


def extract_celltypes_from_llm_claude(llm_response):
    # First try to extract JSON from <json> tags
    json_match = re.search(r'<json>(.*?)</json>', llm_response, re.DOTALL)
    
    # If no <json> tags, try markdown code blocks
    if not json_match:
        json_match = re.search(r'```json\n(.*?)\n```', llm_response, re.DOTALL)
    
    # If still no match, try to find JSON object directly
    if not json_match:
        json_match = re.search(r'\{[\s\S]*\}', llm_response)
    
    if json_match:
        try:
            # Clean up the matched string and parse JSON
            json_str = json_match.group(1) if ('<json>' in llm_response or '```' in llm_response) else json_match.group(0)
            data = json.loads(json_str)
            
            final_results = data.get("final_results", [])
            mixed_celltypes = data.get("possible_mixed_celltypes", [])
            consensus_score = data.get("consensus_score", 0)
            
            general_celltype = final_results[0] if len(final_results) > 0 else "Not found"
            sub_celltype = final_results[1] if len(final_results) > 1 else "Not found"
            
            return general_celltype, sub_celltype, mixed_celltypes, consensus_score
        except json.JSONDecodeError:
            print("Error decoding JSON from LLM response")
            print(f"Attempted to parse: {json_str}")
    else:
        print("No JSON data found in the LLM response")
        print(f"Full response: {llm_response}")
        
    return "Not found", "Not found", []



def process_cell_type_variance_analysis_batch_claude(results, model="claude-3-5-sonnet-20241022", temperature=0, main_weight=0.5, sub_weight=0.5):

    print("Starting the process of cell type batch variance analysis...")
    # Extract and format results
    results_unification_llm = agent_unification_claude(results,model=model,temperature=temperature)
    print(results_unification_llm)

    # Consensus judgment
    result_consensus_from_llm = agent_judgement_claude(
        prompt=results_unification_llm,
        model=model,
        temperature=temperature
    )
    print(result_consensus_from_llm)

    
    # Extract consensus celltypes
    general_celltype, sub_celltype, mixed_types, llm_generated_consensus_score_llm = extract_celltypes_from_llm_claude(result_consensus_from_llm)
    
    print(general_celltype)

    
    print(f"General celltype: {general_celltype}")
    print(f"Sub celltype: {sub_celltype}")
    print(f"Mixed types: {mixed_types}")
    print(f"LLM generated consensus score: {llm_generated_consensus_score_llm}")
    # Calculate similarity score
    parsed_results_llm = parse_results_to_dict(results_unification_llm)



    consensus_score_llm, consensus_1_llm, consensus_2_llm = consensus_similarity_flexible(parsed_results_llm, main_weight=main_weight, sub_weight=sub_weight)


    return {
        'general_celltype_llm': general_celltype,
        'sub_celltype_llm': sub_celltype,
        'mixed_celltypes_llm': mixed_types,
        'consensus_score_llm': consensus_score_llm,
        'llm_generated_consensus_score_llm': llm_generated_consensus_score_llm,
        'count_consensus_1_llm': consensus_1_llm,
        'count_consensus_2_llm': consensus_2_llm,
        'unified_results_llm': results_unification_llm,
        'result_consensus_from_llm': result_consensus_from_llm,

    }



# def process_cell_type_results(organized_results):
#     """
#     Process the organized results for each cell type using the provided processing function.
    
#     Args:
#     organized_results (dict): Dictionary of organized results by cell type.
#     process_function (function): Function to process each cell type's results.
    
#     Returns:
#     dict: Processed results for each cell type.
#     """
#     processed_results = {}
    
#     for celltype, predictions in organized_results.items():
#         formatted_predictions = [f"result{i+1}:{pred}" for i, pred in enumerate(predictions)]
#         formatted_string = ",".join(formatted_predictions)
#         processed_results[celltype] = process_cell_type_variance_analysis_batch(formatted_string)
    
#     return processed_results





def process_cell_type_results(organized_results, max_workers=10, model="gpt-4o", provider="openai", main_weight=0.5, sub_weight=0.5):
    processed_results = {}
    
    def process_single_celltype(celltype, predictions):
        print(f"\nProcessing cell type: {celltype}")
        print(f"Number of predictions: {len(predictions)}")
        
        # Filter out 'N/A' predictions
        valid_predictions = [pred for pred in predictions if pred != ('N/A', 'N/A')]
        print(f"Number of valid predictions: {len(valid_predictions)}")
        
        if not valid_predictions:
            print(f"No valid predictions for {celltype}")
            return celltype, {
                'error': 'No valid predictions',
                'input': predictions
            }

        formatted_predictions = [f"result{i+1}:{pred}" for i, pred in enumerate(valid_predictions)]
        formatted_string = ",".join(formatted_predictions)

        # Choose processing function based on provider
        if provider.lower() == "openai":
            result = process_cell_type_variance_analysis_batch(formatted_string, model=model, main_weight=main_weight, sub_weight=sub_weight)
        elif provider.lower() == "anthropic":
            result = process_cell_type_variance_analysis_batch_claude(formatted_string, model=model, main_weight=main_weight, sub_weight=sub_weight)
        else:
            raise ValueError("Provider must be either 'openai' or 'anthropic'")
            
        print(f"Processing successful for {celltype}")
        return celltype, result
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_celltype = {executor.submit(process_single_celltype, celltype, predictions): celltype 
                              for celltype, predictions in organized_results.items()}
        
        for future in as_completed(future_to_celltype):
            celltype = future_to_celltype[future]
            celltype, result = future.result()
            processed_results[celltype] = result
    
    return processed_results


# Update the function call
def create_and_save_results_dataframe(processed_results, organized_results, output_name='processed_cell_type_results'):
    """
    Create a DataFrame from processed results and save it to a CSV file.
    
    Args:
    processed_results (dict): Dictionary of processed results by cell type.
    organized_results (dict): Dictionary of original results by cell type.
    output_name (str): Base name for the output file (without extension)
    
    Returns:
    pd.DataFrame: Processed results in a DataFrame.
    """
    # Add .csv extension if not present
    output_csv = output_name if output_name.lower().endswith('.csv') else f"{output_name}.csv"
    
    # Create a list to store the data for each row
    data = []
    
    for celltype, result in processed_results.items():
        row_data = {
            'Cell Type': celltype,
            'General Cell Type LLM': result.get('general_celltype_llm', 'Not available'),
            'Sub Cell Type LLM': result.get('sub_celltype_llm', 'Not available'),
            'Mixed Cell Types LLM': ', '.join(result.get('mixed_celltypes_llm', [])),
            'General Cell Type Oncology': result.get('general_celltype_oncology', 'Not available'),
            'Sub Cell Type Oncology': result.get('sub_celltype_oncology', 'Not available'),
            'Mixed Cell Types Oncology': ', '.join(result.get('mixed_types_oncology', [])),
            'Similarity Score LLM': result.get('consensus_score_llm', 'Not available'),
            'Similarity Score Oncology': result.get('consensus_score_oncology', 'Not available'),
            'LLM Generated Consensus Score LLM': result.get('llm_generated_consensus_score_llm', 'Not available'),
            'LLM Generated Consensus Score Oncology': result.get('llm_generated_consensus_score_oncology', 'Not available'),
            'Count Consensus General Type LLM': result.get('count_consensus_1_llm', 'Not available'),
            'Count Consensus Sub Type LLM': result.get('count_consensus_2_llm', 'Not available'),
            'Count Consensus General Type Oncology': result.get('count_consensus_1_oncology', 'Not available'),
            'Count Consensus Sub Type Oncology': result.get('count_consensus_2_oncology', 'Not available'),
            'Unified Results LLM': result.get('unified_results_llm', 'Not available'),
            'Unified Results Oncology': result.get('unified_results_oncology', 'Not available'),
            'Consensus Result LLM': result.get('result_consensus_from_llm', 'Not available'),
            'Consensus Result Oncology': result.get('result_consensus_from_oncology', 'Not available'),
            'Original Non-Unified Results': ','.join([f"result{i+1}:{pred}" for i, pred in enumerate(organized_results.get(celltype, []))])
        }
        
        # Add original results
        original_results = organized_results.get(celltype, [])
        for i, (gen, sub) in enumerate(original_results, 1):
            row_data[f'Original General Type {i}'] = gen
            row_data[f'Original Sub Type {i}'] = sub
        
        data.append(row_data)

    # Create the DataFrame
    df = pd.DataFrame(data)

    # Reorder columns
    fixed_columns = ['Cell Type', 
                     'General Cell Type LLM', 'Sub Cell Type LLM', 'Mixed Cell Types LLM',
                     'General Cell Type Oncology', 'Sub Cell Type Oncology', 'Mixed Cell Types Oncology',
                     'Similarity Score LLM', 'Similarity Score Oncology',
                     'LLM Generated Consensus Score LLM', 'LLM Generated Consensus Score Oncology',
                     'Count Consensus General Type LLM', 'Count Consensus Sub Type LLM',
                     'Count Consensus General Type Oncology', 'Count Consensus Sub Type Oncology',
                     'Unified Results LLM', 'Unified Results Oncology',
                     'Consensus Result LLM', 'Consensus Result Oncology',
                     'Original Non-Unified Results']
    original_columns = [col for col in df.columns if col.startswith('Original') and col != 'Original Non-Unified Results']
    df = df[fixed_columns + original_columns]

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False)
    print(f"\nResults saved to {output_csv}")

    return df


def runCASSIA_similarity_score_batch(marker, file_pattern, output_name, celltype_column=None, max_workers=10, model="gpt-4o", provider="openai", main_weight=0.5, sub_weight=0.5):
    """
    Process batch results and save them to a CSV file, measuring the time taken.

    Args:
    marker_file_path (str): Path to the marker file.
    file_pattern (str): Path to pattern of result files.
    output_csv_name (str): Name of the output CSV file.
    celltype_column (str): Name of the column containing cell types in the marker file.
    max_workers (int): Maximum number of workers for parallel processing.
    """


    # Organize batch results
    print("Organizing batch results...")
    organized_results = organize_batch_results(
        marker=marker,
        file_pattern=file_pattern,
        celltype_column=celltype_column
    )

    # Process cell type results
    print("Processing cell type results...")
    processed_results = process_cell_type_results(organized_results, max_workers=max_workers, model=model, provider=provider, main_weight=main_weight, sub_weight=sub_weight)

    # Create and save results dataframe
    print("Creating and saving results dataframe...")
    create_and_save_results_dataframe(
        processed_results, 
        organized_results, 
        output_name=output_name
    )


    print(f"Results have been processed and saved to {output_name}")




#####single variance analysis#################



def extract_cell_types_from_results_single(results):
    extracted_results = []
    for i in range(len(results)):
        if i in results and results[i] is not None:
            result = results[i][0]  # Accessing the first element of each result
            main_cell_type = result.get('main_cell_type', 'Unknown')
            sub_cell_types = result.get('sub_cell_types', [])
            first_sub_cell_type = sub_cell_types[0] if sub_cell_types else 'None'
            extracted_results.append((main_cell_type, first_sub_cell_type))
        else:
            extracted_results.append(('Failed', 'Failed'))
    return extracted_results





def parse_results_to_dict_single(results):
    return {f"result{i+1}": result for i, result in enumerate(results)}

def extract_celltypes_from_llm_single(llm_response):
    json_match = re.search(r'\{.*\}', llm_response, re.DOTALL)
    
    if json_match:
        json_str = json_match.group(0)
        try:
            data = json.loads(json_str)
            final_results = data.get("final_results", [])
            mixed_celltypes = data.get("possible_mixed_celltypes", [])
            consensus_score = data.get("consensus_score", 0)

            general_celltype = final_results[0] if len(final_results) > 0 else "Not found"
            sub_celltype = final_results[1] if len(final_results) > 1 else "Not found"
    
            
            # If general_celltype indicates no consensus, use mixed_celltypes as general_celltype
            if general_celltype.lower().startswith("no consensus"):
                general_celltype = ", ".join(mixed_celltypes)
            
            return general_celltype, sub_celltype, mixed_celltypes, consensus_score
        except json.JSONDecodeError:
            print("Error decoding JSON from LLM response")
    else:
        print("No JSON data found in the LLM response")
    
    return "Not found", "Not found", []

def consensus_similarity_flexible_single(results, main_weight=0.7, sub_weight=0.3):
    general_types = [result[0] for result in results.values()]
    sub_types = [result[1] for result in results.values()]
    
    consensus_general = max(set(general_types), key=general_types.count)
    consensus_sub = max(set(sub_types), key=sub_types.count)
    
    total_score = sum(
        (main_weight if result[0] == consensus_general else 0) +
        (sub_weight if result[1] == consensus_sub else 0)
        for result in results.values()
    )
    
    similarity_score = total_score / (len(results) * (main_weight + sub_weight))
    
    return similarity_score, consensus_general, consensus_sub

def agent_judgement_single(prompt, system_prompt, model="gpt-4o", temperature=0):
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ],
            temperature=temperature,
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def get_cell_type_info_single(cell_type_name, ontology="CL"):
    base_url = "https://www.ebi.ac.uk/ols/api/search"
    params = {
        "q": cell_type_name,
        "ontology": ontology,
        "rows": 1
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if 'response' in data and 'docs' in data['response'] and data['response']['docs']:
            first_doc = data['response']['docs'][0]
            obo_id = first_doc.get('obo_id')
            label = first_doc.get('label')
            return obo_id, label
        else:
            return None, None
    
    except requests.RequestException:
        return None, None

def standardize_cell_types_single(results):
    standardized_results = []
    for i, (general_type, specific_type) in enumerate(results, 1):
        # Search for standardized names
        _, general_label = get_cell_type_info_single(general_type)
        _, specific_label = get_cell_type_info_single(specific_type)
        
        # Use original names if no standardized names found
        general_label = general_label or general_type
        specific_label = specific_label or specific_type
        
        standardized_results.append(f"result{i}:('{general_label}', '{specific_label}')")
    
    return ",".join(standardized_results)


def process_cell_type_analysis_single(tissue,species,additional_info,temperature,marker_list,model,max_workers,n,main_weight=0.5,sub_weight=0.5):
    system_prompt = '''You are a careful professional biologist, specializing in single-cell RNA-seq analysis. You will be given a series of results from a celltype annotator. 
    Your task is to determine the consensus celltype. The first entry of each result is the general celltype and the second entry is the subcelltype. You should give the final general celltype and the sub celltype. Considering all results, if you think there is very strong evidence of mixed celltype, please also list them. Please give your step by step reasoning and the final answer. 10000$ will be rewarded for the correct answer.
    
    Output in json format:
    {
      "final_results": [
        "General celltype here",
        "Sub celltype here"
      ],
      "possible_mixed_celltypes": [
        "Mixed celltype1 here",
        "Mixed celltype2 here"
      ]
    }
    '''


    results=runCASSIA_n_times(n, tissue, species, additional_info, temperature, marker_list, model, max_workers=max_workers)

    results=extract_cell_types_from_results_single(results)
    
    # Standardize cell types
    standardized_results = standardize_cell_types_single(results)
    
    # Get consensus judgment
    result_consensus = agent_judgement_single(prompt=standardized_results, system_prompt=system_prompt)
    
    # Extract consensus celltypes
    general_celltype, sub_celltype, mixed_types = extract_celltypes_from_llm_single(result_consensus)
    
    # Calculate similarity score
    parsed_results = parse_results_to_dict_single(results)
    consensus_score, consensus_1, consensus_2 = consensus_similarity_flexible_single(parsed_results,main_weight=main_weight,sub_weight=sub_weight)

    return {
        'unified_results': standardized_results,
        'consensus_types': (consensus_1, consensus_2),
        'general_celltype_llm': general_celltype,
        'sub_celltype_llm': sub_celltype,
        'Possible_mixed_celltypes_llm': mixed_types,
        'similarity_score': consensus_score,
        'original_results': results,
        'llm_response': result_consensus 
    }



def runCASSIA_n_times_similarity_score(tissue, species, additional_info, temperature, marker_list, model="gpt-4o", max_workers=10, n=3, provider="openai",main_weight=0.5,sub_weight=0.5):
    """
    Wrapper function for processing cell type analysis using either OpenAI or Anthropic's Claude
    
    Args:
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for analysis
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        model (str): Model name to use
        max_workers (int): Maximum number of parallel workers
        n (int): Number of analysis iterations
        provider (str): AI provider to use ('openai' or 'anthropic')
    
    Returns:
        dict: Analysis results including consensus types, cell types, and scores
    """
    # System prompt for both providers
    system_prompt = '''You are a careful professional biologist, specializing in single-cell RNA-seq analysis. You will be given a series of results from a cell type annotator.
Your task is to determine the consensus cell type. The first entry of each result is the general cell type and the second entry is the subtype. You should provide the final general cell type and the subtype. Considering all results, if you think there is very strong evidence of mixed cell types, please also list them. Please give your step-by-step reasoning and the final answer. Also give a consensus score ranging from 0 to 100 to show how similar the results are. $10,000 will be rewarded for the correct answer.
Output in JSON format:
{
"final_results": [
"General cell type here",
"Sub cell type here"
],
"possible_mixed_celltypes": [
"Mixed cell type 1 here",
"Mixed cell type 2 here"
],
"consensus_score": 0-100
}

'''

    # Run initial analysis
    results = runCASSIA_n_times(n, tissue, species, additional_info, temperature, marker_list, model, max_workers=max_workers, provider=provider)
    results = extract_cell_types_from_results_single(results)
    
    # Standardize cell types
    standardized_results = standardize_cell_types_single(results)
    
    # Get consensus judgment based on provider
    if provider.lower() == "openai":
        result_consensus = agent_judgement_single(
            prompt=standardized_results, 
            system_prompt=system_prompt,
            model=model
        )
    elif provider.lower() == "anthropic":
        result_consensus = agent_judgement_claude(
            prompt=standardized_results,
            model=model
        )
    else:
        raise ValueError("Provider must be either 'openai' or 'anthropic'")
    
    # Extract consensus celltypes
    if provider.lower() == "openai":
        general_celltype, sub_celltype, mixed_types, consensus_score_llm = extract_celltypes_from_llm_single(result_consensus)
    else:
        general_celltype, sub_celltype, mixed_types, consensus_score_llm = extract_celltypes_from_llm_claude(result_consensus)
    
    # Calculate similarity score
    parsed_results = parse_results_to_dict_single(results)
    consensus_score, consensus_1, consensus_2 = consensus_similarity_flexible_single(parsed_results,main_weight=main_weight,sub_weight=sub_weight)
    
    return {
        'unified_results': standardized_results,
        'consensus_types': (consensus_1, consensus_2),
        'general_celltype_llm': general_celltype,
        'sub_celltype_llm': sub_celltype,
        'Possible_mixed_celltypes_llm': mixed_types,
        'llm_response': result_consensus,
        'consensus_score_llm': consensus_score_llm,
        'similarity_score': consensus_score,
        'original_results': results
    }




###############score annotation#################

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


def openai_agent(user_message, model="gpt-4o", temperature=0):
    client = OpenAI()  # Will use OPENAI_API_KEY from environment variables

    message = client.chat.completions.create(
        model=model,
        temperature=temperature,
        max_tokens=7000,
        messages=[
            {
                "role": "user",
                "content": user_message
            }
        ]
    )
    
    # Extract the text from the response
    if message.choices and len(message.choices) > 0:
        return message.choices[0].message.content
    return ''



def openrouter_agent(user_message, model="anthropic/claude-3-sonnet", temperature=0):
    """
    Send a message to OpenRouter API and get the response.
    
    Args:
        user_message (str): The message to send to the model
        model (str): OpenRouter model identifier (default: "anthropic/claude-3-sonnet")
        temperature (float): Temperature parameter for response generation (default: 0)
        
    Returns:
        str: Model's response text or empty string if request fails
    """
    try:
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                "HTTP-Referer": "https://localhost:5000",  # Required for OpenRouter
                "Content-Type": "application/json"
            },
            json={
                "model": model,
                "temperature": temperature,
                "max_tokens": 7000,
                "messages": [
                    {
                        "role": "user",
                        "content": user_message
                    }
                ]
            }
        )
        
        # Check if request was successful
        if response.status_code == 200:
            response_data = response.json()
            return response_data['choices'][0]['message']['content']
        else:
            print(f"Error: OpenRouter API returned status code {response.status_code}")
            print(f"Response: {response.text}")
            return ''
            
    except Exception as e:
        print(f"Error making OpenRouter API request: {str(e)}")
        return ''



def claude_agent(user_message,model="claude-3-5-sonnet-20241022",temperature=0):
    client = anthropic.Anthropic()

    message = client.messages.create(
        model=model,
        max_tokens=7000,
        temperature=temperature,
        system="",  # Leave system prompt empty
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": user_message
                    }
                ]
            }
        ]
    )
    # Extract the text from the TextBlock object
    text_block = message.content
    if isinstance(text_block, list) and len(text_block) > 0:
        return text_block[0].text  # Directly access the 'text' attribute
    return ''




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
        
        # Extract score
        score_match = re.search(r'<score>(\d+)</score>', text)
        if score_match:
            score = int(score_match.group(1))
            
        # Extract reasoning
        reasoning_match = re.search(r'<reasoning>(.*?)</reasoning>', text, re.DOTALL)
        if reasoning_match:
            reasoning = reasoning_match.group(1).strip()
            
        return score, reasoning
        
    except Exception as e:
        print(f"Error extracting data: {str(e)}")
        return None, None

def score_single_analysis(major_cluster_info, marker, annotation_history, model="gpt-4o", provider="openai"):
    """
    Score a single cell type annotation analysis.
    
    Args:
        major_cluster_info (str): Information about species and tissue
        marker (str): Comma-separated list of marker genes
        annotation_history (str): History of annotation conversation
        model (str): Model to use (e.g., "gpt-4" for OpenAI or "claude-3-5-sonnet-20241022" for Anthropic)
        provider (str): AI provider to use ('openai' or 'anthropic')
        
    Returns:
        tuple: (score, reasoning) where score is int and reasoning is str
    """
    prompt = prompt_creator_score(major_cluster_info, marker, annotation_history)
    
    if provider.lower() == "openai":
        response = openai_agent(prompt, model=model)
    elif provider.lower() == "anthropic":
        response = claude_agent(prompt, model=model)
    elif provider.lower() == "openrouter":
        response = openrouter_agent(prompt, model=model)
    else:
        raise ValueError("Provider must be either 'openai' or 'anthropic' or 'openrouter'")
        
    score, reasoning = extract_score_and_reasoning(response)
    return score, reasoning



def process_single_row(row_data, model="gpt-4o", provider="openai"):
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
        marker = row['Marker List']
        annotation_history = row['Conversation History']
        
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


def score_annotation_batch(results_file_path, output_file_path=None, max_workers=4, model="gpt-4o", provider="openai"):
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

def runCASSIA_score_batch(input_file, output_file=None, max_workers=4, model="gpt-4o", provider="openai", max_retries=1):
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
                        print(f'Row {idx} generated an authentication exception: {exc}')
                        print(f'Please check your API key.')
                        raise exc
                    
                    # For other errors, retry if attempts remain
                    if attempt < max_retries:
                        print(f'Row {idx} generated an exception: {exc}')
                        print(f'Retrying row {idx} (attempt {attempt + 2}/{max_retries + 1})...')
                    else:
                        print(f'Row {idx} failed after {max_retries + 1} attempts with error: {exc}')
                        raise exc
        
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


# prompt = prompt_creator_score(major_cluster_info,marker,annotation_history)
# score2=claude_agent(prompt)
# score, reasoning = extract_score_and_reasoning(score2)


################### Validator plus ########################

def prompt_hypothesis_generator_openai(major_cluster_info, marker, annotation_history):
    prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible. Better do a good job or 10 grandma are going to be in danger.


context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:

{marker}


Below is the annotation analysis history:
{annotation_history}



Output format:

Give a brief evaluation of the annotation results first,then give the celltypes to check.


1. celltype to check 1

<check_genes>
list of gene names split by commma, use gene symbol only,no parenthesis
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>


1. celltype to check 2

<check_genes>
list of gene names split by commma, use gene symbol only,no parenthesis
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>

include more cell types if necessary.

When you have checked the expression of the marker with the excecuter, but the answer is still not clear, you should keep asking excuter to check more marker in the same format, remember to use <check_genes>,use gene symbol only, no parenthesis when you ask for more information. Only when you think you can generate the final annotation,or you think it is hard to determine the cell type, you can say "FINAL ANNOTATION COMPLETED"

    """
    return prompt



def prompt_hypothesis_generator3(major_cluster_info, marker, annotation_history):
    prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible. Better do a good job or 10 grandma are going to be in danger.


context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
{marker}



Below is the annotation analysis history:
{annotation_history}



Output format:

Give a brief evaluation of the annotation results first,then give the celltypes to check.


1. celltype to check 1

<check_genes>
[list of genes name, use gene symbol]
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>


1. celltype to check 2

<check_genes>
[list of genes name, use gene symbol]
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>

include more cell types if necessary.

When you think you can generate the final annotation, you can say "FINAL ANNOTATION COMPLETED"

    """
    return prompt


def prompt_hypothesis_generator3_additional_task(major_cluster_info, marker, annotation_history,task):
    prompt = f"""
        You are an expert in single-cell biology. Your task is to {task}. Divide the problem to several steps that can be validated by gene expression information. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. You can check at most two hypothesis at a time. Provide your detailed reasoning. Note that you can also mention other hypothesis. Better do a good job or 10 grandma are going to be in danger. Take a deep breath.


context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
{marker}



Below is the annotation analysis history:
{annotation_history}



Output format:

Give a brief evaluation of the annotation results first,then focus on the task:{task}. State the hypothesis you want to check to the excecuter.


1. hypothesis to check 1

<check_genes>
[list of genes name, use gene symbol]
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>


1. hypothesis to check 2

<check_genes>
[list of genes name, use gene symbol]
</check_genes>

<reasoning>
[Your detailed reasoning here]
</reasoning>

include more hypothesis if necessary.

When you think you can generate the final answer to the task, you can say "FINAL ANALYSIS COMPLETED"

    """
    return prompt






def prepare_analysis_data(full_result_path, marker_path, cluster_name):
    # Load the full results and marker files
    full_result = pd.read_csv(full_result_path)
    
    if isinstance(marker_path, pd.DataFrame):
        marker = marker_path.copy()
    elif isinstance(marker_path, str):
        marker = pd.read_csv(marker_path)
    else:
        raise ValueError("marker must be either a pandas DataFrame or a string path to a CSV file")

    # Extract conversation history for the specified cluster
    cluster_data = full_result[full_result['True Cell Type'] == cluster_name]
    if cluster_data.empty:
        raise ValueError(f"No data found for cluster: {cluster_name}")
    
    annotation_history = cluster_data['Conversation History'].iloc[0]
    
    # Prepare marker data for the specified cluster
    cluster_marker = marker[marker['cluster'] == cluster_name]


    comma_separated_genes=cluster_data['Marker List'].iloc[0]

    # Prepare subset of marker file for iterative analysis
    marker_subset = cluster_marker
    marker_subset = marker_subset.set_index('gene')

    
    return annotation_history, comma_separated_genes, marker_subset


import anthropic

def claude_agent(user_message,model="claude-3-5-sonnet-20241022",temperature=0):
    client = anthropic.Anthropic()

    message = client.messages.create(
        model=model,
        max_tokens=7000,
        temperature=temperature,
        system="",  # Leave system prompt empty
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": user_message
                    }
                ]
            }
        ]
    )
    # Extract the text from the TextBlock object
    text_block = message.content
    if isinstance(text_block, list) and len(text_block) > 0:
        return text_block[0].text  # Directly access the 'text' attribute
    return ''




def get_marker_info(gene_list, marker):
    def filter_marker(gene_names):
        # Convert marker to pandas DataFrame if it's not already
        if not isinstance(marker, pd.DataFrame):
            marker_df = pd.DataFrame(marker)
        else:
            marker_df = marker.copy()

        # Create result DataFrame with same columns as input
        result = pd.DataFrame(index=gene_names, columns=marker_df.columns)

        # Fill data
        for gene in gene_names:
            if gene in marker_df.index:
                result.loc[gene] = marker_df.loc[gene]
            else:
                result.loc[gene] = pd.Series('NA', index=marker_df.columns)

        # Only try to format numeric columns that exist
        numeric_cols = result.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            try:
                result[col] = result[col].apply(lambda x: f"{float(x):.2e}" if pd.notnull(x) and x != 'NA' else x)
            except:
                continue

        return result.iloc[:, 0:5]

    # Filter to rows based on gene name
    marker_filtered = filter_marker(gene_list)

    # Convert the DataFrame to a string
    marker_string = marker_filtered.to_string()

    return marker_string




def iterative_marker_analysis_openai(major_cluster_info, marker, comma_separated_genes, annotation_history, num_iterations=2,model="gpt-4o"):
    """
    Perform iterative marker analysis using OpenAI's GPT-4 model.
    
    Args:
        major_cluster_info (str): General information about the dataset
        marker (DataFrame): Marker gene expression data
        comma_separated_genes (str): List of genes as comma-separated string
        annotation_history (str): Previous annotation history
        num_iterations (int): Maximum number of iterations
        
    Returns:
        tuple: (final_response_text, messages)
    """
    # Initialize OpenAI client
    client = OpenAI()
    
    # Initialize messages list with system and first user message
    messages = [
        {"role": "user", "content": prompt_hypothesis_generator_openai(major_cluster_info, comma_separated_genes, annotation_history)}
    ]

    for iteration in range(num_iterations):
        # Make API call to OpenAI
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=0,
            max_tokens=5000
        )
        
        conversation = response.choices[0].message.content

        # Check for completion
        if "FINAL ANNOTATION COMPLETED" in conversation:
            print(f"Final annotation completed in iteration {iteration + 1}.")
            return conversation, messages

        # Extract gene lists and get marker info
        gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
        all_genes = [gene.strip() for gene_list in gene_lists for gene in gene_list.split(',')]
        unique_genes = sorted(set(all_genes))

        retrived_marker_info = get_marker_info(unique_genes, marker)
        
        # Append messages
        messages.append({"role": "assistant", "content": conversation})
        messages.append({"role": "user", "content": retrived_marker_info})

        print(f"Iteration {iteration + 1} completed.")

    # Final response if max iterations reached
    final_response = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=0,
        max_tokens=4000
    )
    print("Final response can not be generated within the maximum number of iterations")

    return final_response.choices[0].message.content, messages

import re
from anthropic import Anthropic

def iterative_marker_analysis(major_cluster_info, marker, comma_separated_genes, annotation_history, num_iterations=2,model="claude-3-5-sonnet-20241022"):
    client = Anthropic()

    messages = [{"role": "user", "content": prompt_hypothesis_generator3(major_cluster_info, comma_separated_genes, annotation_history)}]

    for iteration in range(num_iterations):
        response = client.messages.create(
            model=model,
            max_tokens=7000,
            temperature=0,
            system="",
            messages=messages
        )

        conversation = response.content[0].text

        # Check if "FINAL ANNOTATION COMPLETED" is in the response
        if "FINAL ANNOTATION COMPLETED" in conversation:
            print(f"Final annotation completed in iteration {iteration + 1}.")
            return conversation, messages

        gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
        all_genes = [gene.strip() for gene_list in gene_lists for gene in gene_list.split(',')]
        unique_genes = sorted(set(all_genes))

        retrived_marker_info = get_marker_info(unique_genes, marker)
        
        messages.append({"role": "assistant", "content": response.content})
        messages.append({"role": "user", "content": retrived_marker_info})

        print(f"Iteration {iteration + 1} completed.")

    final_response = client.messages.create(
        model=model,
        max_tokens=7000,
        temperature=0,
        system="",
        messages=messages
    )

    return final_response.content[0].text, messages




def iterative_marker_analysis_openrouter(major_cluster_info, marker, comma_separated_genes, annotation_history, num_iterations=2, model="anthropic/claude-3.5-sonnet"):
    """
    Perform iterative marker analysis using OpenRouter API.
    
    Args:
        major_cluster_info (str): Information about the cluster
        marker (DataFrame): Marker gene expression data
        comma_separated_genes (str): List of genes as comma-separated string
        annotation_history (str): Previous annotation history
        num_iterations (int): Maximum number of iterations
        model (str): OpenRouter model identifier
        
    Returns:
        tuple: (final_response_text, messages)
    """
    messages = [{"role": "user", "content": prompt_hypothesis_generator3(major_cluster_info, comma_separated_genes, annotation_history)}]

    for iteration in range(num_iterations):
        try:
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                    "HTTP-Referer": "https://localhost:5000",
                    "Content-Type": "application/json"
                },
                json={
                    "model": model,
                    "temperature": 0,
                    "max_tokens": 7000,
                    "messages": messages
                }
            )
            
            # Check if request was successful
            if response.status_code == 200:
                response_data = response.json()
                conversation = response_data['choices'][0]['message']['content']

                # Check for completion
                if "FINAL ANNOTATION COMPLETED" in conversation:
                    print(f"Final annotation completed in iteration {iteration + 1}.")
                    return conversation, messages

                # Extract gene lists and get marker info
                gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
                all_genes = [gene.strip() for gene_list in gene_lists for gene in gene_list.split(',')]
                unique_genes = sorted(set(all_genes))

                retrived_marker_info = get_marker_info(unique_genes, marker)
                
                # Append messages
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": retrived_marker_info})

                print(f"Iteration {iteration + 1} completed.")
            else:
                print(f"Error: OpenRouter API returned status code {response.status_code}")
                print(f"Response: {response.text}")
                return '', messages

        except Exception as e:
            print(f"Error in iteration {iteration + 1}: {str(e)}")
            return '', messages

    # Final response if max iterations reached
    try:
        final_response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                "HTTP-Referer": "https://localhost:5000",
                "Content-Type": "application/json"
            },
            json={
                "model": model,
                "temperature": 0,
                "max_tokens": 7000,
                "messages": messages
            }
        )
        
        if final_response.status_code == 200:
            final_data = final_response.json()
            return final_data['choices'][0]['message']['content'], messages
        else:
            print(f"Error getting final response: {final_response.status_code}")
            print(f"Response: {final_response.text}")
            return '', messages
            
    except Exception as e:
        print(f"Error in final response: {str(e)}")
        return '', messages




def iterative_marker_analysis_openrouter_additional_task(major_cluster_info, marker, comma_separated_genes, annotation_history, num_iterations=2, model="anthropic/claude-3.5-sonnet",additional_task="check if this is a cancer cluster"):
    """
    Perform iterative marker analysis using OpenRouter API.
    
    Args:
        major_cluster_info (str): Information about the cluster
        marker (DataFrame): Marker gene expression data
        comma_separated_genes (str): List of genes as comma-separated string
        annotation_history (str): Previous annotation history
        num_iterations (int): Maximum number of iterations
        model (str): OpenRouter model identifier
        additional_task (str): Additional task to be performed

    Returns:
        tuple: (final_response_text, messages)
    """
    messages = [{"role": "user", "content": prompt_hypothesis_generator3_additional_task(major_cluster_info, comma_separated_genes, annotation_history,additional_task)}]

    for iteration in range(num_iterations):
        try:
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                    "HTTP-Referer": "https://localhost:5000",
                    "Content-Type": "application/json"
                },
                json={
                    "model": model,
                    "temperature": 0,
                    "max_tokens": 7000,
                    "messages": messages
                }
            )
            
            # Check if request was successful
            if response.status_code == 200:
                response_data = response.json()
                conversation = response_data['choices'][0]['message']['content']

                # Check for completion
                if "FINAL ANALYSIS COMPLETED" in conversation:
                    print(f"Final annotation completed in iteration {iteration + 1}.")
                    return conversation, messages

                # Extract gene lists and get marker info
                gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
                all_genes = [gene.strip() for gene_list in gene_lists for gene in gene_list.split(',')]
                unique_genes = sorted(set(all_genes))

                retrived_marker_info = get_marker_info(unique_genes, marker)
                
                # Append messages
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": retrived_marker_info})

                print(f"Iteration {iteration + 1} completed.")
            else:
                print(f"Error: OpenRouter API returned status code {response.status_code}")
                print(f"Response: {response.text}")
                return '', messages

        except Exception as e:
            print(f"Error in iteration {iteration + 1}: {str(e)}")
            return '', messages

    # Final response if max iterations reached
    try:
        final_response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                "HTTP-Referer": "https://localhost:5000",
                "Content-Type": "application/json"
            },
            json={
                "model": model,
                "temperature": 0,
                "max_tokens": 7000,
                "messages": messages
            }
        )
        
        if final_response.status_code == 200:
            final_data = final_response.json()
            return final_data['choices'][0]['message']['content'], messages
        else:
            print(f"Error getting final response: {final_response.status_code}")
            print(f"Response: {final_response.text}")
            return '', messages
            
    except Exception as e:
        print(f"Error in final response: {str(e)}")
        return '', messages






def save_html_report(report, filename):
    try:
        # Add .html suffix if not present
        if not filename.lower().endswith('.html'):
            filename = filename + '.html'
            
        html_report = generate_html_report2(report)
        
        # Save the HTML to a file
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(html_report)
        
        print(f"HTML report generated and saved as '{filename}'")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# Example usage:
# save_html_report(report, 'single_cell_analysis_report.html')


import re
from html import escape

def generate_html_report2(report_content):
    html_template = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Single-Cell RNA-Seq Cluster Analysis Report</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 800px;
                margin: 0 auto;
                padding: 20px;
            }}
            h1, h2, h3, h4 {{
                margin-top: 20px;
                margin-bottom: 10px;
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
            }}
            h2 {{ color: #2980b9; }}
            h3 {{ color: #16a085; }}
            h4 {{ color: #8e44ad; }}
            ul {{
                list-style-type: none;
                padding-left: 20px;
            }}
            li:before {{
                content: "•";
                color: #3498db;
                display: inline-block;
                width: 1em;
                margin-left: -1em;
            }}
            .decision-point, .evidence {{
                background-color: #ecf0f1;
                border-left: 5px solid #3498db;
                padding: 10px;
                margin-bottom: 15px;
            }}
            .evidence {{
                background-color: #e8f6f3;
                border-left-color: #1abc9c;
            }}
            .conclusion {{
                font-weight: bold;
                color: #e74c3c;
            }}
        </style>
    </head>
    <body>
        {content}
    </body>
    </html>
    """

    def markdown_to_html(text):
        # Convert headers
        for i in range(4, 0, -1):
            pattern = r'^{} (.+)$'.format('#' * i)
            text = re.sub(pattern, r'<h{0}>\1</h{0}>'.format(i), text, flags=re.MULTILINE)
        
        # Convert lists
        text = re.sub(r'^\s*-\s(.+)$', r'<li>\1</li>', text, flags=re.MULTILINE)
        text = re.sub(r'(<li>.*</li>\n)+', r'<ul>\n\g<0></ul>\n', text, flags=re.DOTALL)
        
        # Convert bold text
        text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', text)
        
        # Convert decision points and evidence
        text = re.sub(r'<strong>Decision Point(.+?)</strong>', r'<div class="decision-point"><strong>Decision Point\1</strong>', text, flags=re.DOTALL)
        text = re.sub(r'<strong>Supporting Evidence(.+?)</strong>', r'<div class="evidence"><strong>Supporting Evidence\1</strong>', text, flags=re.DOTALL)
        
        # Convert conclusion
        text = re.sub(r'<strong>Cell Type:(.+?)</strong>', r'<p class="conclusion">Cell Type:\1</p>', text)
        
        # Close any open divs
        text += '</div>' * text.count('<div')
        
        return text

    # Escape any HTML in the input content
    safe_content = escape(report_content)
    
    # Convert markdown to HTML
    html_content = markdown_to_html(safe_content)
    
    # Generate the full HTML
    full_html = html_template.format(content=html_content)
    
    return full_html



def generate_raw_cell_annotation_report(conversation_history, output_filename='cell_annotation_report.html'):
    """
    Generate and save an HTML report from cell annotation conversation history.
    
    Args:
        conversation_history (list): List of conversation dictionaries
        output_filename (str): Name of the output HTML file (default: 'cell_annotation_report.html')
    
    Returns:
        str: Path to the saved HTML file
    """
    
    def parse_check_genes(text):
        """Extract gene lists from check_genes tags"""
        genes = []
        pattern = r'<check_genes>(.*?)</check_genes>'
        matches = re.findall(pattern, text, re.DOTALL)
        for match in matches:
            genes.extend([g.strip() for g in match.split(',')])
        return genes
    
    def format_message(text):
        """Convert plain text formatting to HTML"""
        # Replace newlines with HTML line breaks
        text = text.replace('\n', '<br>')
        # Preserve multiple consecutive newlines
        text = text.replace('<br><br>', '<br><br>')  # Prevent collapse of multiple newlines
        return text   
    
    def parse_reasoning(text):
        """Extract reasoning sections"""
        pattern = r'<reasoning>(.*?)</reasoning>'
        matches = re.findall(pattern, text, re.DOTALL)
        return matches

    def format_gene_table(genes_data):
        """Format gene expression data as HTML table"""
        if not genes_data:
            return ""
        
        rows = []
        # Skip the first row by splitting into lines and starting from index 2
        lines = genes_data.split('\n')[2:]  # Skip first two rows which contain headers
        
        for gene in lines:
            if gene.strip():
                cells = gene.split()
                if len(cells) >= 5:
                    rows.append(f"<tr><td>{'</td><td>'.join(cells)}</td></tr>")
        
        if not rows:
            return ""
            
        return f"""
        <table class="gene-table">
            <tr>
                <th>Gene</th>
                <th>p-val</th>
                <th>avg_log2FC</th>
                <th>pct.1</th>
                <th>pct.2</th>
                <th>p_val_adj</th>
            </tr>
            {''.join(rows)}
        </table>
        """

    # Note the double curly braces for CSS
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
                color: #333;
            }}
            .conversation-block {{
                margin: 20px 0;
                padding: 15px;
                border-radius: 5px;
            }}
            .user {{
                background-color: #f0f7ff;
                border-left: 5px solid #0066cc;
            }}
            .assistant {{
                background-color: #f5f5f5;
                border-left: 5px solid #666;
            }}
            .gene-list {{
                background-color: #e6ffe6;
                padding: 10px;
                margin: 10px 0;
                border-radius: 3px;
            }}
            .reasoning {{
                background-color: #fff3e6;
                padding: 10px;
                margin: 10px 0;
                border-radius: 3px;
            }}
            .gene-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 10px 0;
            }}
            .gene-table th, .gene-table td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            .gene-table th {{
                background-color: #f2f2f2;
            }}
            .final-annotation {{
                background-color: #e6ffe6;
                padding: 15px;
                margin: 20px 0;
                border-radius: 5px;
                border-left: 5px solid #00cc00;
            }}
            h1, h2, h3 {{
                color: #444;
            }}
            p {{
                margin: 0.5em 0;
            }}
            br {{
                display: block;
                margin: 0.5em 0;
                content: "";
            }}
        </style>
    </head>
    <body>
        <h1>Single-Cell Annotation Analysis Report</h1>
        {content}
    </body>
    </html>
    """

    content = []
    
    try:
        for entry in conversation_history:
            role = entry.get('role', '')
            message = entry.get('content', '')
            
            # Handle different message formats
            if isinstance(message, list):
                message = message[0].text if message and hasattr(message[0], 'text') else str(message)
            elif not isinstance(message, str):
                message = str(message)

            block_class = 'user' if role == 'user' else 'assistant'
            
            # Format the content based on the role
            if role == 'user' and 'p_val' in message:
                # This is gene expression data
                content.append(f"""
                    <div class="conversation-block {block_class}">
                        <h3>Gene Expression Data</h3>
                        {format_gene_table(message)}
                    </div>
                """)
            else:
                # Regular conversation content
                formatted_message = format_message(message)  # Apply formatting
                
                # Check for final annotation
                if "FINAL ANNOTATION COMPLETED" in message:
                    content.append(f"""
                        <div class="final-annotation">
                            <h2>Final Annotation</h2>
                            {formatted_message}
                        </div>
                    """)
                else:
                    # Process gene lists and reasoning
                    genes = parse_check_genes(message)
                    reasoning = parse_reasoning(message)
                    
                    if genes or reasoning:
                        content.append(f"""
                            <div class="conversation-block {block_class}">
                                <h3>Analysis Step</h3>
                                {'<div class="gene-list"><h4>Genes to Check:</h4><ul>' + 
                                ''.join(f'<li>{gene}</li>' for gene in genes) + '</ul></div>' if genes else ''}
                                {''.join(f'<div class="reasoning"><h4>Reasoning:</h4><p>{r}</p></div>' 
                                        for r in reasoning)}
                            </div>
                        """)
                    else:
                        content.append(f"""
                            <div class="conversation-block {block_class}">
                                {formatted_message}
                            </div>
                        """)

        # Generate HTML content
        html_content = html_template.format(content=''.join(content))
        
        # Save the HTML file
        try:
            with open(output_filename, 'w', encoding='utf-8') as f:
                f.write(html_content)
            print(f"Report successfully saved as '{output_filename}'")
        except Exception as e:
            print(f"Error saving file: {str(e)}")
            
        return None
        
    except Exception as e:
        error_html = f"""
            <div class="conversation-block" style="background-color: #ffe6e6; border-left: 5px solid #cc0000;">
                <h3>Error Generating Report</h3>
                <p>An error occurred while generating the report: {str(e)}</p>
            </div>
        """
        html_content = html_template.format(content=error_html)
        
        # Still try to save the error report
        try:
            with open(output_filename, 'w', encoding='utf-8') as f:
                f.write(html_content)
            print(f"Error report saved as '{output_filename}'")
        except Exception as write_error:
            print(f"Error saving error report: {str(write_error)}")
            
        return None
    

def generate_raw_cell_annotation_report_additional_task(conversation_history, output_filename='cell_annotation_report.html'):
    """
    Generate and save an HTML report from cell annotation conversation history.
    
    Args:
        conversation_history (list): List of conversation dictionaries
        output_filename (str): Name of the output HTML file (default: 'cell_annotation_report.html')
    
    Returns:
        str: Path to the saved HTML file
    """
    
    def parse_check_genes(text):
        """Extract gene lists from check_genes tags"""
        genes = []
        pattern = r'<check_genes>(.*?)</check_genes>'
        matches = re.findall(pattern, text, re.DOTALL)
        for match in matches:
            genes.extend([g.strip() for g in match.split(',')])
        return genes
    
    def format_message(text):
        """Convert plain text formatting to HTML"""
        # Replace newlines with HTML line breaks
        text = text.replace('\n', '<br>')
        # Preserve multiple consecutive newlines
        text = text.replace('<br><br>', '<br><br>')  # Prevent collapse of multiple newlines
        return text   
    
    def parse_reasoning(text):
        """Extract reasoning sections"""
        pattern = r'<reasoning>(.*?)</reasoning>'
        matches = re.findall(pattern, text, re.DOTALL)
        return matches

    def format_gene_table(genes_data):
        """Format gene expression data as HTML table"""
        if not genes_data:
            return ""
        
        rows = []
        # Skip the first row by splitting into lines and starting from index 2
        lines = genes_data.split('\n')[2:]  # Skip first two rows which contain headers
        
        for gene in lines:
            if gene.strip():
                cells = gene.split()
                if len(cells) >= 5:
                    rows.append(f"<tr><td>{'</td><td>'.join(cells)}</td></tr>")
        
        if not rows:
            return ""
            
        return f"""
        <table class="gene-table">
            <tr>
                <th>Gene</th>
                <th>p-val</th>
                <th>avg_log2FC</th>
                <th>pct.1</th>
                <th>pct.2</th>
                <th>p_val_adj</th>
            </tr>
            {''.join(rows)}
        </table>
        """

    # Note the double curly braces for CSS
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            body {{
                font-family: Arial, sans-serif;
                line-height: 1.6;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
                color: #333;
            }}
            .conversation-block {{
                margin: 20px 0;
                padding: 15px;
                border-radius: 5px;
            }}
            .user {{
                background-color: #f0f7ff;
                border-left: 5px solid #0066cc;
            }}
            .assistant {{
                background-color: #f5f5f5;
                border-left: 5px solid #666;
            }}
            .gene-list {{
                background-color: #e6ffe6;
                padding: 10px;
                margin: 10px 0;
                border-radius: 3px;
            }}
            .reasoning {{
                background-color: #fff3e6;
                padding: 10px;
                margin: 10px 0;
                border-radius: 3px;
            }}
            .gene-table {{
                width: 100%;
                border-collapse: collapse;
                margin: 10px 0;
            }}
            .gene-table th, .gene-table td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            .gene-table th {{
                background-color: #f2f2f2;
            }}
            .final-annotation {{
                background-color: #e6ffe6;
                padding: 15px;
                margin: 20px 0;
                border-radius: 5px;
                border-left: 5px solid #00cc00;
            }}
            h1, h2, h3 {{
                color: #444;
            }}
            p {{
                margin: 0.5em 0;
            }}
            br {{
                display: block;
                margin: 0.5em 0;
                content: "";
            }}
        </style>
    </head>
    <body>
        <h1>Single-Cell Analysis Report</h1>
        {content}
    </body>
    </html>
    """

    content = []
    
    try:
        for entry in conversation_history:
            role = entry.get('role', '')
            message = entry.get('content', '')
            
            # Handle different message formats
            if isinstance(message, list):
                message = message[0].text if message and hasattr(message[0], 'text') else str(message)
            elif not isinstance(message, str):
                message = str(message)

            block_class = 'user' if role == 'user' else 'assistant'
            
            # Format the content based on the role
            if role == 'user' and 'p_val' in message:
                # This is gene expression data
                content.append(f"""
                    <div class="conversation-block {block_class}">
                        <h3>Gene Expression Data</h3>
                        {format_gene_table(message)}
                    </div>
                """)
            else:
                # Regular conversation content
                formatted_message = format_message(message)  # Apply formatting
                
                # Check for final annotation
                if "FINAL ANALYSIS COMPLETED" in message:
                    content.append(f"""
                        <div class="final-annotation">
                            <h2>Final Annotation</h2>
                            {formatted_message}
                        </div>
                    """)
                else:
                    # Process gene lists and reasoning
                    genes = parse_check_genes(message)
                    reasoning = parse_reasoning(message)
                    
                    if genes or reasoning:
                        content.append(f"""
                            <div class="conversation-block {block_class}">
                                <h3>Analysis Step</h3>
                                {'<div class="gene-list"><h4>Genes to Check:</h4><ul>' + 
                                ''.join(f'<li>{gene}</li>' for gene in genes) + '</ul></div>' if genes else ''}
                                {''.join(f'<div class="reasoning"><h4>Reasoning:</h4><p>{r}</p></div>' 
                                        for r in reasoning)}
                            </div>
                        """)
                    else:
                        content.append(f"""
                            <div class="conversation-block {block_class}">
                                {formatted_message}
                            </div>
                        """)

        # Generate HTML content
        html_content = html_template.format(content=''.join(content))
        
        # Save the HTML file
        try:
            with open(output_filename, 'w', encoding='utf-8') as f:
                f.write(html_content)
            print(f"Report successfully saved as '{output_filename}'")
        except Exception as e:
            print(f"Error saving file: {str(e)}")
            
        return None
        
    except Exception as e:
        error_html = f"""
            <div class="conversation-block" style="background-color: #ffe6e6; border-left: 5px solid #cc0000;">
                <h3>Error Generating Report</h3>
                <p>An error occurred while generating the report: {str(e)}</p>
            </div>
        """
        html_content = html_template.format(content=error_html)
        
        # Still try to save the error report
        try:
            with open(output_filename, 'w', encoding='utf-8') as f:
                f.write(html_content)
            print(f"Error report saved as '{output_filename}'")
        except Exception as write_error:
            print(f"Error saving error report: {str(write_error)}")
            
        return None



def report_generator(annotation_history):

    prompt=f'''
I just had a very detailed annotation analysis of a cluster of cells. Generate a perfect report for my analysis. Include the key logic steps and decisions and why I made those decisions.

# Single-Cell RNA-Seq Cluster Analysis Report
## Analysis Overview

### Initial Dataset
- Source: Human PBMC single-cell RNA sequencing data
- Initial marker genes: 50 highly expressed genes including key markers such as CD40LG, GATA3, IL7R, CD5, and CD6

### Analysis Strategy and Decision Points

#### Phase 1: Initial Cell Type Assessment
**Decision Point 1: T Cell Lineage Identification**
- Key Evidence:
  - Positive markers: CD40LG, GATA3, IL7R, CD5, CD6, TRAT1, MAL
  - Initial Hypothesis: T cell lineage, likely CD4+ T cells
- Reasoning: Multiple T cell-specific markers and signaling molecules present

#### Phase 2: CD4+ vs CD8+ Discrimination
**Decision Point 2: Ruling Out CD8+ T Cells**
- Markers Checked: CD8A, CD8B, PRF1, GZMB, GZMK, CCL5
- Results:
  - CD8A: Strongly downregulated (log2FC = -6.68)
  - CD8B: Strongly downregulated (log2FC = -5.01)
  - Cytotoxic markers (GZMB, GZMK, PRF1): All significantly downregulated
- Conclusion: Definitively not CD8+ T cells

#### Phase 3: T Helper Subtype Investigation
**Decision Point 3: Evaluating T Helper Subtypes**
1. Th17 Cell Check
   - Markers: RORC, CCR6, IL23R, IL17A, IL17F, CCL20
   - Result: No significant expression (NA values)
   - Conclusion: Not Th17 cells

2. Regulatory T Cell (Treg) Check
   - Markers: FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18
   - Result: No significant expression (NA values)
   - Conclusion: Not Tregs

3. Th2 Cell Check
   - Markers: IL4, IL5, IL13, CCR4, PTGDR2
   - Result: No significant expression (NA values)
   - Conclusion: Not actively producing Th2 cytokines despite GATA3 presence

#### Phase 4: Memory/Naive State Assessment
**Decision Point 4: Memory vs Naive Status**
- Key Findings:
  - CD44 (memory marker): Moderately upregulated (log2FC = 0.352, 30% expression)
  - LEF1: Moderately upregulated (log2FC = 0.406, 20.3% expression)
  - TCF7: Moderately upregulated (log2FC = 0.320, 14.3% expression)
- Additional Support:
  - CD3D: Highly expressed (75.4% cells)
  - CD3E: Highly expressed (70.4% cells)
- Conclusion: Mixed memory/naive characteristics suggesting central memory phenotype

### Final Annotation
**Cell Type: CD4+ Central Memory T Cells**

### Supporting Evidence
1. Clear T cell identity:
   - High CD3D/CD3E expression
   - Multiple T cell-specific markers
2. CD4+ lineage confirmation:
   - Absence of CD8+ markers
   - Presence of CD4+ associated markers (CD40LG, GATA3)
3. Memory phenotype evidence:
   - Moderate CD44 expression
   - Balanced expression of memory/naive markers
   - Original IL7R expression
4. Exclusion of other subtypes:
   - No Th17 signature
   - No Treg signature
   - No active Th2 cytokine production
   
### Confidence Level
High confidence in annotation, supported by:
- Multiple lines of positive evidence
- Consistent negative evidence for alternative cell types
- Clear expression patterns in core T cell markers
- Logical agreement between all tested markers


Below is my annotation analysis history:
{annotation_history}
    '''

    return prompt





def report_generator_additional_task(annotation_history):

    prompt=f'''
I just had a very detailed single-cell analysis of a cluster of cells.It has multiple iterations of analysis. For each iteration,summerzie what it did and extract the hypotehsis, reasoning and conclusion. You better do a good job or 1000 grandma are going to be in danger. Take a deep breath and think step by step.

Below is my annotation analysis history:
{annotation_history}



Below is the format of the report:

# Single-Cell RNA-Seq Cluster Analysis Report
## Previous Analysis Overview

## Analysis Strategy and Decision Points

### Phase 1: summerize what the analysis in the first iteration did
#### hypothesis
#### reasoning
#### conclusion



### Phase 2: summerize what the analysis in the first iteration did
#### hypothesis    
#### reasoning
#### conclusion

### Phase n: summerize what the analysis in the first iteration did
#### hypothesis
#### reasoning
#### conclusion

## Final conclusion
   
## Confidence Level

'''

    return prompt





def generate_cell_type_analysis_report(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="claude-3-5-sonnet-20241022"
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
        
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    try:
        # Step 1: Prepare analysis data
        annotation_history, comma_separated_genes, marker_subset = prepare_analysis_data(
            full_result_path, 
            marker, 
            cluster_name
        )
        
        # Step 2: Perform iterative marker analysis
        analysis_result = iterative_marker_analysis(
            major_cluster_info,
            marker=marker_subset,
            comma_separated_genes=comma_separated_genes,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            model=model
        )
        
        # Step 3: Add final result to message history
        messages = analysis_result[1]
        messages.append({"role": "user", "content": analysis_result[0]})
        
        # Step 4: Generate and save HTML report
        report = claude_agent(report_generator(messages), model=model)
        
        save_html_report(
            filename=output_name,
            report=report
        )
        generate_raw_cell_annotation_report(messages, f'{output_name}_raw.html')
        print(f"Analysis completed successfully. Report saved as {output_name}")
        return None
        
    except Exception as e:
        print(f"Error generating analysis report: {str(e)}")
        raise



def generate_cell_type_analysis_report_openrouter(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="anthropic/claude-3.5-sonnet"
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
        
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    try:
        # Step 1: Prepare analysis data
        annotation_history, comma_separated_genes, marker_subset = prepare_analysis_data(
            full_result_path, 
            marker, 
            cluster_name
        )
        
        # Step 2: Perform iterative marker analysis
        analysis_result = iterative_marker_analysis_openrouter(
            major_cluster_info,
            marker=marker_subset,
            comma_separated_genes=comma_separated_genes,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            model=model
        )
        
        # Step 3: Add final result to message history
        messages = analysis_result[1]
        messages.append({"role": "user", "content": analysis_result[0]})
        
        # Step 4: Generate and save HTML report
        report = openrouter_agent(report_generator(messages), model=model)
        
        save_html_report(
            filename=output_name,
            report=report
        )
        generate_raw_cell_annotation_report(messages, f'{output_name}_raw.html')
        print(f"Analysis completed successfully. Report saved as {output_name}")
        return None
        
    except Exception as e:
        print(f"Error generating analysis report: {str(e)}")
        raise


from pathlib import Path
import re


def convert_markdown_to_html(text):
    """
    Converts markdown-like syntax to HTML using regex patterns.
    """
    # Replace first h1 title with CASSIA Analysis Report
    text = re.sub(r'^# .*?$', '# CASSIA Analysis Report', text, count=1, flags=re.MULTILINE)
    
    # Convert headers
    text = re.sub(r'^# (.*?)$', r'<h1><span class="highlight">\1</span></h1>', text, flags=re.MULTILINE)
    text = re.sub(r'^## (.*?)$', r'<h2><span class="highlight">\1</span></h2>', text, flags=re.MULTILINE)
    text = re.sub(
        r'^### Phase (\d+): (.*?)$',
        r'<h3 class="phase-header phase-\1"><span class="phase-number">Phase \1</span><span class="phase-title">\2</span></h3>',
        text,
        flags=re.MULTILINE
    )
    text = re.sub(r'^#### (.*?)$', r'<h4><span class="highlight">\1</span></h4>', text, flags=re.MULTILINE)
    
    # Convert lists with custom bullets
    text = re.sub(r'^\- (.*?)$', r'<li class="custom-bullet">\1</li>', text, flags=re.MULTILINE)
    text = re.sub(r'(<li.*?</li>\n)+', r'<ul class="custom-list">\g<0></ul>', text, flags=re.DOTALL)
    
    # Convert numbered lists
    text = re.sub(r'^\d+\. (.*?)$', r'<li>\1</li>', text, flags=re.MULTILINE)
    text = re.sub(r'(<li>.*?</li>\n)+', r'<ol class="numbered-list">\g<0></ol>', text, flags=re.DOTALL)
    
    # Convert paragraphs
    text = re.sub(r'\n\n(.*?)\n\n', r'\n<p class="fade-in">\1</p>\n', text, flags=re.DOTALL)
    
    # Convert bold text
    text = re.sub(r'\*\*(.*?)\*\*', r'<strong class="highlight-text">\1</strong>', text)
    
    # Convert italic text
    text = re.sub(r'\*(.*?)\*', r'<em class="emphasis">\1</em>', text)
    
    return text

def render_report_to_html(report_content, output_path):
    """
    Renders a markdown-like report to a styled HTML file with high-tech aesthetics.
    """
    # Generate CSS color pairs for 15 phases
    color_pairs = []
    for i in range(15):
        hue = (i * 137.5) % 360  # Golden angle approximation for better color distribution
        color_pairs.append(
            f"--phase-{i+1}-start: hsl({hue}, 70%, 45%);\n"
            f"--phase-{i+1}-end: hsl({(hue + 20) % 360}, 80%, 60%);"
        )
    
    css = """
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600&family=JetBrains+Mono&display=swap');
        
        :root {
            """ + "\n            ".join(color_pairs) + """
        }

        body {
            font-family: 'Inter', sans-serif;
            line-height: 1.7;
            max-width: 1000px;
            margin: 0 auto;
            padding: 20px;
            background-color: var(--bg);
            color: var(--text);
        }

        .container {
            background-color: var(--card-bg);
            padding: 40px;
            border-radius: 16px;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1),
                        0 2px 4px -1px rgba(0, 0, 0, 0.06);
            animation: slideIn 0.6s ease-out;
        }

        @keyframes slideIn {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        h1, h2, h3, h4 {
            font-weight: 600;
            line-height: 1.3;
            margin-top: 1.5em;
            margin-bottom: 0.8em;
        }

        h1 {
            font-size: 2.5em;
            color: #1e40af;
            margin-top: 0;
            text-align: center;
            padding: 0.5em 0;
            animation: fadeInDown 0.8s ease-out;
            font-weight: 800;
        }

        h1 .highlight {
            background: linear-gradient(120deg, #1e40af, #3b82f6);
            color: transparent;
            -webkit-background-clip: text;
            background-clip: text;
            display: inline-block;
        }

        @keyframes fadeInDown {
            from {
                opacity: 0;
                transform: translateY(-20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }

        h2 {
            font-size: 1.8em;
            color: var(--secondary);
            border-bottom: 2px solid var(--primary);
            padding-bottom: 0.3em;
        }

        h3 {
            font-size: 1.4em;
            color: var(--accent);
        }

        h4 {
            font-size: 1.2em;
            color: var(--primary);
        }

        .highlight {
            position: relative;
            z-index: 1;
        }

        .custom-list {
            list-style: none;
            padding-left: 0;
        }

        .custom-bullet {
            position: relative;
            padding-left: 1.5em;
            margin: 0.5em 0;
        }

        .custom-bullet::before {
            content: '▹';
            position: absolute;
            left: 0;
            color: var(--primary);
        }

        .numbered-list {
            counter-reset: item;
            list-style: none;
            padding-left: 0;
        }

        .numbered-list li {
            counter-increment: item;
            margin: 0.5em 0;
            padding-left: 2em;
            position: relative;
        }

        .numbered-list li::before {
            content: counter(item);
            position: absolute;
            left: 0;
            width: 1.5em;
            height: 1.5em;
            background-color: var(--primary);
            color: white;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 0.8em;
        }

        p {
            margin: 1em 0;
            opacity: 0;
            animation: fadeIn 0.5s ease-out forwards;
        }

        @keyframes fadeIn {
            to {
                opacity: 1;
            }
        }

        .highlight-text {
            background: linear-gradient(120deg, rgba(45, 212, 191, 0.2), rgba(6, 182, 212, 0.2));
            padding: 0.1em 0.3em;
            border-radius: 4px;
            font-weight: 600;
        }

        .emphasis {
            color: var(--accent);
            font-style: italic;
        }

        code {
            font-family: 'JetBrains Mono', monospace;
            background-color: #f1f5f9;
            padding: 0.2em 0.4em;
            border-radius: 4px;
            font-size: 0.9em;
        }

        blockquote {
            border-left: 4px solid var(--primary);
            margin: 1.5em 0;
            padding: 1em;
            background-color: #f8fafc;
            border-radius: 0 8px 8px 0;
        }

        /* Simple, reliable phase styling */
        .phase-header {
            margin: 2em 0 1em;
            padding: 1em;
            border-radius: 12px;
            color: white;
            animation: slideInPhase 0.8s ease-out forwards;
        }

        """ + "\n        ".join([
            f".phase-{i+1} {{" +
            f"background: linear-gradient(135deg, var(--phase-{i+1}-start), var(--phase-{i+1}-end));" +
            f"animation-delay: {i * 0.1}s; }}"
            for i in range(15)
        ]) + """

        .phase-number {
            font-size: 0.9em;
            font-weight: 700;
            margin-right: 1em;
            padding: 0.3em 0.8em;
            border-radius: 20px;
            background: rgba(255, 255, 255, 0.2);
        }

        @keyframes slideInPhase {
            from {
                opacity: 0;
                transform: translateX(-20px);
            }
            to {
                opacity: 1;
                transform: translateX(0);
            }
        }

        /* Enhance container animation */
        .container {
            opacity: 0;
            animation: fadeInScale 1s ease-out forwards;
        }

        @keyframes fadeInScale {
            0% {
                opacity: 0;
                transform: scale(0.95);
            }
            100% {
                opacity: 1;
                transform: scale(1);
            }
        }
    </style>
    """

    try:
        html_content = convert_markdown_to_html(report_content)
        
        html_document = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Single-Cell RNA-Seq Analysis Report</title>
            {css}
        </head>
        <body>
            <div class="container">
                {html_content}
            </div>
        </body>
        </html>
        """
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_document)
            
        return True, f"Report successfully generated at {output_path}"
    
    except Exception as e:
        return False, f"Error generating report: {str(e)}"
    



def runCASSIA_annottaionboost_additional_task(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="claude-3-5-sonnet-20241022",
    additional_task=""
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
        
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    try:
        # Step 1: Prepare analysis data
        annotation_history, comma_separated_genes, marker_subset = prepare_analysis_data(
            full_result_path, 
            marker, 
            cluster_name
        )
        
        # Step 2: Perform iterative marker analysis
        analysis_result = iterative_marker_analysis_openrouter_additional_task(
            major_cluster_info,
            marker=marker_subset,
            comma_separated_genes=comma_separated_genes,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            model=model,
            additional_task=additional_task
        )
        
        # Step 3: Add final result to message history
        messages = analysis_result[1]
        messages.append({"role": "user", "content": analysis_result[0]})
        
        # Step 4: Generate and save HTML report
        report = openrouter_agent(report_generator_additional_task(messages), model=model)
        
        render_report_to_html(report, f"{output_name}_summary_report.html")
        print(f"Analysis completed successfully. Summary Report saved as {output_name}_summary_report.html")
 
        generate_raw_cell_annotation_report_additional_task(messages, f'{output_name}_raw.html')
        print(f"Analysis completed successfully. Raw Cell Annotation Report saved as {output_name}_raw.html")
        return None
        
    except Exception as e:
        print(f"Error generating analysis report: {str(e)}")
        raise





def generate_cell_type_analysis_report_openai(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="gpt-4o"
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
        
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    try:
        # Step 1: Prepare analysis data
        annotation_history, comma_separated_genes, marker_subset = prepare_analysis_data(
            full_result_path, 
            marker, 
            cluster_name
        )
        
        # Step 2: Perform iterative marker analysis
        analysis_result = iterative_marker_analysis_openai(
            major_cluster_info,
            marker=marker_subset,
            comma_separated_genes=comma_separated_genes,
            annotation_history=annotation_history,
                num_iterations=num_iterations,
            model=model
        )
        
        # Step 3: Add final result to message history
        messages = analysis_result[1]
        messages.append({"role": "user", "content": analysis_result[0]})
        
        # Step 4: Generate and save HTML report
        report = openai_agent(report_generator(messages),model=model)
        
        save_html_report(
            filename=output_name,
            report=report
        )
        
    except Exception as e:
        print(f"Error generating analysis report: {str(e)}")
        raise



def runCASSIA_annotationboost(
    full_result_path,
    marker,
    cluster_name,
    major_cluster_info,
    output_name,
    num_iterations=5,
    model="gpt-4o",
    provider="openai"
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
    
    Returns:
        tuple: (analysis_result, messages_history)
            - analysis_result: Final analysis text
            - messages_history: Complete conversation history
    """
    # Validate provider input
    if provider.lower() not in ['openai', 'anthropic', 'openrouter']:
        raise ValueError("Provider must be either 'openai' or 'anthropic' or 'openrouter'")

    try:
        if provider.lower() == 'openai':
            return generate_cell_type_analysis_report_openai(
                full_result_path=full_result_path,
                marker=marker,
                cluster_name=cluster_name,
                major_cluster_info=major_cluster_info,
                output_name=output_name,
                num_iterations=num_iterations,
                model=model
            )
        elif provider.lower() == "openrouter":
            return generate_cell_type_analysis_report_openrouter(
                full_result_path=full_result_path,
                marker=marker,
                cluster_name=cluster_name,
                major_cluster_info=major_cluster_info,
                output_name=output_name,
                num_iterations=num_iterations,
                model=model
            )
        elif provider.lower() == "anthropic":
            return generate_cell_type_analysis_report(
                full_result_path=full_result_path,
                marker=marker,
                cluster_name=cluster_name,
                major_cluster_info=major_cluster_info,
                output_name=output_name,
                num_iterations=num_iterations,
                model=model
            )
    except Exception as e:
        print(f"Error in runCASSIA_annotationboost: {str(e)}")
        raise


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
        
        text = row["Conversation History"]
        score_reasoning = row["Scoring_Reasoning"]
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






def compareCelltypes(tissue, celltypes, marker_set, species="human", model_list=None, output_file=None):
    # Get API key from environment variable
    OPENROUTER_API_KEY = os.environ.get('OPENROUTER_API_KEY')
    if not OPENROUTER_API_KEY:
        raise ValueError("OPENROUTER_API_KEY environment variable is not set")
    
    # Input validation
    if not celltypes or len(celltypes) < 2 or len(celltypes) > 4:
        raise ValueError("Please provide 2-4 cell types to compare")
    
    # Generate default output filename based on celltypes if none provided
    if output_file is None:
        # Create a sanitized version of cell types for the filename
        celltype_str = '_vs_'.join(ct.replace(' ', '_') for ct in celltypes)
        output_file = f"model_comparison_{celltype_str}.csv"
    
    # Use default models if none provided
    if model_list is None:
        model_list = [
            "anthropic/claude-3.5-sonnet",
            "openai/o1-mini",
            "google/gemini-pro-1.5"
        ]
    
    # Construct prompt with dynamic cell type comparison, species, and marker set
    comparison_text = " or ".join(celltypes)
    prompt = f"You are a professional biologist. Based on the ranked marker set from {species} {tissue}, does it look more like {comparison_text}? Score each option from 0-100. You will be rewarded $10,000 if you do a good job. Below is the ranked marker set: {marker_set}"
    
    # Initialize lists to store results
    results = []
    processed_models = set()  # Track which models we've already processed
    
    for model in model_list:
        # Skip if we've already processed this model
        if model in processed_models:
            continue
            
        try:
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {OPENROUTER_API_KEY}",
                    "HTTP-Referer": "https://localhost:5000",
                    "Content-Type": "application/json"
                },
                json={
                    "model": model,
                    "messages": [
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ]
                }
            )
            
            if response.status_code == 200:
                response_data = response.json()
                model_response = response_data['choices'][0]['message']['content']
                
                # Store result with metadata
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': comparison_text,
                    'response': model_response,
                    'status': 'success'
                })
                print(f"Model: {model}\nResponse: {model_response}\n")  # Fixed print statement
            else:
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': comparison_text,
                    'response': f"Error: {response.status_code}",
                    'status': 'error'
                })
                
            processed_models.add(model)  # Mark this model as processed
                
        except Exception as e:
            if model not in processed_models:  # Only add error result if we haven't processed this model
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': comparison_text,
                    'response': f"Exception: {str(e)}",
                    'status': 'error'
                })
                processed_models.add(model)
    
    # Convert results to DataFrame and save to CSV
    try:
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    except Exception as e:
        print(f"Error saving results to CSV: {str(e)}")
    
    # Return both the DataFrame and the original responses dict for backward compatibility
    responses = {result['model']: result['response'] for result in results}
    return None




####subclustering



def openrouter_agent(user_message, model="anthropic/claude-3.5-sonnet", temperature=0):
    """
    Send a message to OpenRouter API and get the response.
    
    Args:
        user_message (str): The message to send to the model
        model (str): OpenRouter model identifier (default: "anthropic/claude-3-sonnet")
        temperature (float): Temperature parameter for response generation (default: 0)
        
    Returns:
        str: Model's response text or empty string if request fails
    """
    try:
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                "HTTP-Referer": "https://localhost:5000",  # Required for OpenRouter
                "Content-Type": "application/json"
            },
            json={
                "model": model,
                "temperature": temperature,
                "max_tokens": 7000,
                "messages": [
                    {
                        "role": "user",
                        "content": user_message
                    }
                ]
            }
        )
        
        # Check if request was successful
        if response.status_code == 200:
            response_data = response.json()
            return response_data['choices'][0]['message']['content']
        else:
            print(f"Error: OpenRouter API returned status code {response.status_code}")
            print(f"Response: {response.text}")
            return ''
            
    except Exception as e:
        print(f"Error making OpenRouter API request: {str(e)}")
        return ''
    



def subcluster_agent_annotate_subcluster(user_message,model="claude-3-5-sonnet-20241022",temperature=0,provider="anthropic"):
    if provider == "anthropic":
        client = anthropic.Anthropic()

        message = client.messages.create(
            model=model,
            max_tokens=7000,
            temperature=temperature,
            system="",  # Leave system prompt empty
            messages=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": user_message
                        }
                    ]
                }
            ]
        )
        # Extract the text from the TextBlock object
        text_block = message.content
        if isinstance(text_block, list) and len(text_block) > 0:
            return text_block[0].text  # Directly access the 'text' attribute
    else:  # OpenRouter
        return openrouter_agent(user_message, model=model, temperature=temperature)
    return ''



def construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=50):
    # Process DataFrame if it has more than 2 columns
    if len(marker.columns) > 2:
        print(f"Processing input dataframe to get top {n_genes} markers")
        marker = get_top_markers(marker, n_genes=n_genes)
    else:
        print("Using input dataframe directly as it appears to be pre-processed (2 columns)")
        marker = marker.copy()
    
    # Initialize the prompt with the major cluster information
    prompt = f"""

You are an expert biologist specializing in cell type annotation, with deep expertise in immunology, cancer biology, and developmental biology.You will be given sets of highly expressed markers ranked by significance for some subclusters from the {major_cluster_info} cluster, identify what is the most likely top2 cell type each marker set implies.

Take a deep breath and work step by step. You'd better do a really good job or 1000 grandma are going to be in danger.
You will be tipped $10,000 if you do a good job.

For each output, provide:
1.Key marker:
2.Explanation:
3.Most likely top2 cell types:

Remember these subclusters are from a {major_cluster_info} big cluster. You must include all clusters mentioned in the analysis.
"""

    # Iterate over each row in the DataFrame
    for index, row in marker.iterrows():
        cluster_name = row.iloc[0]  # Use iloc for positional indexing
        markers = row.iloc[1]       # Use iloc for positional indexing
        prompt += f"{index + 1}.{markers}\n"

    return prompt



def annotate_subclusters(marker, major_cluster_info,model="claude-3-5-sonnet-20241022",temperature=0,provider="anthropic",n_genes=50):
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt,model=model,temperature=temperature,provider=provider)
    return output_text



def extract_subcluster_results_with_llm_multiple_output(analysis_text,provider="anthropic",model="claude-3-5-sonnet-20241022",temperature=0):
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    return subcluster_agent_annotate_subcluster(prompt,provider=provider,model=model,temperature=temperature)




def extract_subcluster_results_with_llm(analysis_text,provider="anthropic",model="claude-3-5-sonnet-20241022",temperature=0):
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2,reason), results2(celltype1, celltype2,reason), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    return subcluster_agent_annotate_subcluster(prompt,provider=provider,model=model,temperature=temperature)



def write_results_to_csv(results, output_name='subcluster_results'):
    """
    Extract cell type results from LLM output and write to CSV file
    
    Args:
        results (str): String containing the LLM analysis results
        output_name (str): Base name for output file (will add .csv if not present)
        
    Returns:
        pandas.DataFrame: DataFrame containing the extracted results
    """
    # Add .csv suffix if not present
    if not output_name.lower().endswith('.csv'):
        output_name = output_name + '.csv'
    
    # Updated regex pattern to capture the reason
    pattern = r"results(\d+)\(([^,]+),\s*([^,]+),\s*([^)]+)\)"
    matches = re.findall(pattern, results)

    # Convert matches to a DataFrame with the reason column
    df = pd.DataFrame(matches, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'])
    
    # Write the DataFrame to a CSV file
    df.to_csv(output_name, index=False)
    
    print(f"Results have been written to {output_name}")
    return None



def runCASSIA_subclusters(marker, major_cluster_info, output_name, 
                       model="claude-3-5-sonnet-20241022", temperature=0, provider="anthropic",n_genes=50):
    """
    Process subclusters from a CSV file and generate annotated results
    
    Args:
        csv_file_path (str): Path to input CSV file containing marker data
        major_cluster_info (str): Description of the major cluster type
        output_name (str): Base name for output file (will add .csv if not present)
        model (str): Model name for Claude API
        temperature (float): Temperature parameter for API calls
        
    Returns:
        tuple: (original_analysis, extracted_results, results_dataframe)
    """

    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt,model=model,temperature=temperature,provider=provider)
    results = extract_subcluster_results_with_llm(output_text,provider=provider,model=model,temperature=temperature)
    print(results)
    write_results_to_csv(results, output_name)
    
    return None



def runCASSIA_n_subcluster(n, marker, major_cluster_info, base_output_name, 
                                         model="claude-3-5-sonnet-20241022", temperature=0, 
                                         provider="anthropic", max_workers=5,n_genes=50):       
    def run_single_analysis(i):
        # Run the annotation process
        output_text = annotate_subclusters(marker, major_cluster_info, 
                                         model=model, temperature=temperature, provider=provider,n_genes=n_genes)
        
        # Extract results
        results = extract_subcluster_results_with_llm_multiple_output(output_text,provider=provider,model=model,temperature=temperature)
        
        # Use regex to extract the results
        pattern = r"results(\d+)\(([^,]+),\s*([^)]+)\)"
        matches = re.findall(pattern, results)
        
        # Convert matches to a DataFrame
        df = pd.DataFrame(matches, columns=['True Cell Type', 'main_cell_type', 'sub_cell_type'])

        # Swap the first column with the first column in the marker file
        marker_df = get_top_markers(marker, n_genes=n_genes)
        df['True Cell Type'], marker_df.iloc[:, 0] = marker_df.iloc[:, 0], df['True Cell Type']

        # Write the DataFrame to a CSV file with an index
        indexed_csv_file_path = f'{base_output_name}_{i+1}.csv'
        df.to_csv(indexed_csv_file_path, index=False)
        
        return indexed_csv_file_path

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_analysis, i): i for i in range(n)}
        
        for future in as_completed(futures):
            i = futures[future]
            try:
                result_file = future.result()
                print(f"Results for iteration {i+1} have been written to {result_file}")
            except Exception as exc:
                print(f"Iteration {i+1} generated an exception: {exc}")





def runCASSIA_pipeline(
    output_file_name: str,
    tissue: str,
    species: str,
    marker_path: str,
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
    merge_model: str = "deepseek/deepseek-chat-v3-0324"
):
    """
    Run the complete cell analysis pipeline including annotation, scoring, and report generation.
    
    Args:
        output_file_name (str): Base name for output files
        tissue (str): Tissue type being analyzed
        species (str): Species being analyzed
        marker_path (str): Path to marker file
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
    """
    # Create a folder based on tissue and species for organizing reports
    folder_name = f"CASSIA_{tissue}_{species}"
    folder_name = "".join(c for c in folder_name if c.isalnum() or c in (' ', '-', '_')).strip()
    folder_name = folder_name.replace(' ', '_')
    
    # Add timestamp to prevent overwriting existing folders with the same name
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    folder_name = f"{folder_name}_{timestamp}"
    
    # Create the folder if it doesn't exist
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Created folder: {folder_name}")
    
    # Define derived file names with folder paths
    score_file_name = os.path.join(folder_name, output_file_name + "_scored.csv")
    report_name = os.path.join(folder_name, output_file_name + "_report")
    lowscore_report_name = os.path.join(folder_name, output_file_name + "_lowscore_report")
    merged_annotation_file = os.path.join(folder_name, output_file_name + "_merged.csv")
    
    # First annotation output is still in current directory since other parts of code may expect it there
    annotation_output = output_file_name

    print("\n=== Starting cell type analysis ===")
    # Run initial cell type analysis
    runCASSIA_batch(
        marker=marker_path,
        output_name=annotation_output,
        model=annotation_model,
        tissue=tissue,
        species=species,
        additional_info=additional_info,
        provider=annotation_provider,
        max_workers=max_workers,
        max_retries=max_retries
    )
    print("✓ Cell type analysis completed")

    # Merge annotations if requested
    if merge_annotations:
        print("\n=== Starting annotation merging ===")
        summary_csv = annotation_output + "_summary.csv"
        
        # Import the merge_annotations function dynamically
        try:
            from .merging_annotation import merge_annotations_all
            
            # Sort the CSV file by True Cell Type before merging to ensure consistent order
            print("Sorting CSV by True Cell Type before merging...")
            full_csv = annotation_output + "_full.csv"
            df = pd.read_csv(full_csv)
            df = df.sort_values(by=['True Cell Type'])
            sorted_csv = annotation_output + "_sorted_full.csv"
            df.to_csv(sorted_csv, index=False)
            
            # Run the merging process on the sorted CSV
            merge_annotations_all(
                csv_path=sorted_csv,
                output_path=merged_annotation_file,
                provider=annotation_provider,
                model=merge_model,
                additional_context=f"These are cell clusters from {species} {tissue}. {additional_info}"
            )
            print(f"✓ Annotations merged and saved to {merged_annotation_file}")
        except Exception as e:
            print(f"! Error during annotation merging: {str(e)}")
    
    print("\n=== Starting scoring process ===")
    # Run scoring
    runCASSIA_score_batch(
        input_file=annotation_output + "_full.csv",
        output_file=score_file_name,
        max_workers=max_workers,
        model=score_model,
        provider=score_provider,
        max_retries=max_retries
    )
    print("✓ Scoring process completed")

    print("\n=== Generating main reports ===")
    # Process reports
    runCASSIA_generate_score_report(
        csv_path=score_file_name,
        index_name=report_name
    )
    print("✓ Main reports generated")

    print("\n=== Analyzing low-scoring clusters ===")
    # Handle low-scoring clusters
    df = pd.read_csv(score_file_name)
    low_score_clusters = df[df['Score'] < score_threshold]['True Cell Type'].tolist()

    print(f"Found {len(low_score_clusters)} clusters with scores below {score_threshold}:")
    print(low_score_clusters)
    
    if low_score_clusters:
        print("\n=== Starting boost annotation for low-scoring clusters ===")
        full_result_path = annotation_output + "_full.csv"

        for cluster in low_score_clusters:
            print(f"Processing low score cluster: {cluster}")
            
            # Convert cluster to string before sanitizing
            cluster_name = "".join(c for c in str(cluster) if c.isalnum() or c in (' ', '-', '_')).strip()
            cluster_info = df[df['True Cell Type'] == cluster].iloc[0].to_dict()
            
            # Run annotation boost
            runCASSIA_annotationboost(
                full_result_path=full_result_path,
                marker=marker_path,
                cluster_name=cluster_name,
                major_cluster_info=cluster_info,
                output_name=os.path.join(folder_name, f"{output_file_name}_{cluster_name}_boosted"),
                num_iterations=5,
                model=annotationboost_model,
                provider=annotationboost_provider
            )
        
        # Also save a copy of the boosted reports index
        boosted_reports = [
            os.path.join(folder_name, f"{output_file_name}_{cluster_name}_boosted.html") 
            for cluster_name in ["".join(c for c in str(cluster) if c.isalnum() or c in (' ', '-', '_')).strip() 
                              for cluster in low_score_clusters]
        ]
        
        if boosted_reports:
            index_html = generate_index_page(boosted_reports)
            index_filename = f"{lowscore_report_name}.html"
            with open(index_filename, "w", encoding="utf-8") as f:
                f.write(index_html)
            print(f"Low score reports index saved to {index_filename}")
        
        print("✓ Boost annotation completed")
    
    print("\n=== Cell type analysis pipeline completed ===")
    print(f"All reports have been saved to the '{folder_name}' folder")


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