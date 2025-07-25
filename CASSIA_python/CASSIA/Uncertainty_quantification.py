from concurrent.futures import ThreadPoolExecutor, as_completed
try:
    from .tools_function import *
except ImportError:
    from tools_function import *

import time

try:
    from .llm_utils import *
except ImportError:
    from llm_utils import *
import re



def runCASSIA_batch_n_times(n, marker, output_name="cell_type_analysis_results", model="google/gemini-2.5-flash-preview", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column=None, gene_column_name=None, max_workers=10, batch_max_workers=5, provider="openrouter", max_retries=1, validator_involvement="v1"):
    """
    Run multiple batch cell type analyses in parallel.
    
    Args:
        n (int): Number of batch analyses to run
        marker: DataFrame or path to CSV file containing marker data
        output_name (str): Base name for output files
        model (str): Model name to use
        temperature (float): Temperature parameter for the model
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for analysis
        celltype_column (str): Name of column containing cell types
        gene_column_name (str): Name of column containing gene lists
        max_workers (int): Maximum number of workers for each batch
        batch_max_workers (int): Maximum number of concurrent batch runs
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a custom URL)
        max_retries (int): Maximum number of retries for failed analyses
        validator_involvement (str): Validator involvement level ('v0' for high involvement, 'v1' for moderate involvement)
    
    Returns:
        None: Results are saved to files
    """
    def single_batch_run(i):
        output_json_name = f"{output_name}_{i+1}.json"
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
            max_retries=max_retries,
            validator_involvement=validator_involvement
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
    index, tissue, species, additional_info, temperature, marker_list, model, provider, validator_involvement = args
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
            provider=provider,
            validator_involvement=validator_involvement
        )
        end_time = time.time()
        print(f"Finished analysis {index+1} in {end_time - start_time:.2f} seconds")
        return index, result
    except Exception as e:
        print(f"Error in analysis {index+1}: {str(e)}")
        return index, None





def runCASSIA_n_times(n, tissue, species, additional_info, temperature, marker_list, model, max_workers=10, provider="openrouter", validator_involvement="v1"):
    """
    Run multiple cell type analyses in parallel.
    
    Args:
        n (int): Number of analyses to run
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for analysis
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        model (str): Model name to use
        max_workers (int): Maximum number of parallel workers
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a custom URL)
    
    Returns:
        dict: Dictionary of analysis results indexed by iteration number
    """
    print(f"Starting {n} parallel analyses")
    start_time = time.time()
    
    results = {}
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks with the provider parameter
        future_to_index = {
            executor.submit(
                run_single_analysis, 
                (i, tissue, species, additional_info, temperature, marker_list, model, provider, validator_involvement)
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
    """
    Robust parsing function that tries multiple regex patterns to handle different model outputs.
    """
    
    # Define multiple patterns to try, in order of preference
    patterns = [
        # Original strict pattern: result1:(cell_type, subtype) or result1:[cell_type, subtype]
        r"result(\d+):[\(\[]([^\]\)]+)[\)\]]",
        
        # Case insensitive with optional spaces: Result 1: (cell_type, subtype)
        r"(?i)result\s*(\d+)\s*:\s*[\(\[]([^\]\)]+)[\)\]]",
        
        # Simple numbered format: 1. (cell_type, subtype) or 1: (cell_type, subtype)
        r"(\d+)[\.\:]\s*[\(\[]([^\]\)]+)[\)\]]",
        
        # Even more flexible: any number followed by parentheses/brackets
        r"(\d+)[^\(\[]*[\(\[]([^\]\)]+)[\)\]]",
        
        # Pattern for results without explicit numbering but with parentheses
        r"[\(\[]([^,\]\)]+),\s*([^\]\)]+)[\)\]]",
        
        # Very loose pattern: find any content in parentheses with comma
        r"[\(\[]([^,\]\)]{3,}),\s*([^\]\)]{2,})[\)\]]"
    ]
    
    parsed_results = {}
    
    # Try each pattern until we find matches
    for pattern_idx, pattern in enumerate(patterns):
        matches = re.findall(pattern, result_string, re.IGNORECASE)
        
        if matches:
            print(f"Successfully parsed using pattern {pattern_idx + 1}: {pattern}")
            break
    
    # If no numbered patterns worked, try the non-numbered ones differently
    if not matches and pattern_idx >= 4:  # If we reached the non-numbered patterns
        # For non-numbered patterns, create sequential numbering
        sequential_matches = re.findall(patterns[4], result_string, re.IGNORECASE)
        if sequential_matches:
            matches = [(str(i+1), cell_types) for i, cell_types in enumerate(sequential_matches)]
            print(f"Created sequential numbering from pattern matches")
        else:
            sequential_matches = re.findall(patterns[5], result_string, re.IGNORECASE)
            if sequential_matches:
                matches = [(str(i+1), cell_types) for i, cell_types in enumerate(sequential_matches)]
                print(f"Created sequential numbering from loose pattern matches")
    
    # Parse each result
    for match in matches:
        try:
            if len(match) == 2:
                result_num, cell_types = match
            else:
                # Handle case where match might have different structure
                result_num = str(len(parsed_results) + 1)
                cell_types = str(match)
            
            # Handle different cell type formats
            if isinstance(cell_types, tuple):
                # If cell_types is already a tuple (from non-numbered patterns)
                cell_type_list = [str(ct).strip().strip("'\"") for ct in cell_types]
            else:
                # Split cell types, handling potential commas within cell type names
                # Try different splitting approaches
                if ',' in cell_types:
                    # Split by comma and clean up
                    cell_type_list = [ct.strip().strip("'\"") for ct in cell_types.split(',')]
                else:
                    # If no comma, try to split by other delimiters
                    for delimiter in [';', '|', ' and ', ' & ']:
                        if delimiter in cell_types:
                            cell_type_list = [ct.strip().strip("'\"") for ct in cell_types.split(delimiter)]
                            break
                    else:
                        # If no delimiter found, treat as single cell type
                        cell_type_list = [cell_types.strip().strip("'\"")]
            
            # Clean up cell types - remove common prefixes/suffixes that models might add
            cleaned_cell_types = []
            for ct in cell_type_list:
                ct = ct.strip()
                # Remove common model artifacts
                ct = re.sub(r'^(cell type:?\s*|type:?\s*)', '', ct, flags=re.IGNORECASE)
                ct = re.sub(r'\s*(cell type|type)\s*$', '', ct, flags=re.IGNORECASE)
                if ct:  # Only add non-empty cell types
                    cleaned_cell_types.append(ct)
            
            # Ensure we have at least two cell types (main and sub)
            while len(cleaned_cell_types) < 2:
                cleaned_cell_types.append("N/A")
            
            parsed_results[f"result{result_num}"] = tuple(cleaned_cell_types[:2])
            
        except Exception as e:
            print(f"Error parsing match {match}: {str(e)}")
            continue
    
    # If still no results, try one last desperate attempt
    if not parsed_results:
        print("No patterns matched. Attempting fallback parsing...")
        # Look for any two consecutive words that might be cell types
        fallback_pattern = r'\b([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s*,?\s*([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\b'
        fallback_matches = re.findall(fallback_pattern, result_string)
        
        for i, (general, specific) in enumerate(fallback_matches[:5]):  # Limit to first 5 matches
            if general.lower() != specific.lower():  # Avoid duplicates
                parsed_results[f"result{i+1}"] = (general, specific)
        
        if parsed_results:
            print(f"Fallback parsing found {len(parsed_results)} potential cell type pairs")
    
    if not parsed_results:
        print("Warning: No results found to parse from input:")
        print(f"Input string: {result_string[:200]}..." if len(result_string) > 200 else result_string)
    
    return parsed_results



def extract_celltypes_from_llm(llm_response, provider="openai", single_analysis=False):
    """
    Unified function to extract cell type information from LLM responses.
    Handles differences in response format between different providers.
    
    Args:
        llm_response: The text response from the LLM
        provider: The provider that generated the response ("openai", "anthropic", "openrouter", or a custom URL)
        single_analysis: Whether this is from a single analysis (different return format)
        
    Returns:
        If single_analysis=False:
            general_celltype, sub_celltype, mixed_celltypes, consensus_score
        If single_analysis=True:
            general_celltype, sub_celltype, mixed_celltypes
    """
    import json
    import re
    
    # Default return values
    default_not_found = "Not found"
    default_score = 0
    
    # Try multiple extraction methods for all providers
    json_match = None
    
    # First try to extract JSON from <json> tags (common for Claude and some custom APIs)
    if not json_match:
        json_match = re.search(r'<json>(.*?)</json>', llm_response, re.DOTALL)
        
    # Try json> format (sometimes Claude drops the opening <)
        if not json_match:
            json_match = re.search(r'json>(.*?)</json>', llm_response, re.DOTALL)
    
    # Try markdown code blocks (common for OpenAI and some custom APIs)
    if not json_match:
        json_match = re.search(r'```json\n(.*?)\n```', llm_response, re.DOTALL)
        
    # Try code blocks with different formatting
    if not json_match:
        json_match = re.search(r'```\s*json\s*(.*?)\s*```', llm_response, re.DOTALL)
    
    # If still no match, try to find JSON object directly
    if not json_match:
        json_match = re.search(r'\{[\s\S]*\}', llm_response)
    
    # Process the match if found
    if json_match:
        try:
            # Extract the JSON content based on the matched pattern
            if '<json>' in llm_response or 'json>' in llm_response:
                json_str = json_match.group(1)
            elif '```json' in llm_response or '```' in llm_response:
                json_str = json_match.group(1)
            else:
                json_str = json_match.group(0)
                
            # Parse the JSON
            data = json.loads(json_str)
            
            final_results = data.get("final_results", [])
            mixed_celltypes = data.get("possible_mixed_celltypes", [])
            consensus_score = data.get("consensus_score", default_score)
            
            general_celltype = final_results[0] if len(final_results) > 0 else default_not_found
            sub_celltype = final_results[1] if len(final_results) > 1 else default_not_found
            
            # For single analysis, if general_celltype indicates no consensus, use mixed_celltypes
            if single_analysis and general_celltype.lower().startswith("no consensus"):
                general_celltype = ", ".join(mixed_celltypes)
            
            if single_analysis:
                return general_celltype, sub_celltype, mixed_celltypes
            else:
                return general_celltype, sub_celltype, mixed_celltypes, consensus_score
        except json.JSONDecodeError:
            print("Error decoding JSON from LLM response")
            print(f"Attempted to parse: {json_str if 'json_str' in locals() else 'No JSON found'}")
    else:
        print("No JSON data found in the LLM response")
        print(f"Full response: {llm_response[:500]}...")  # Print first 500 chars to avoid overwhelming output
    
    if single_analysis:
        return default_not_found, default_not_found, []
    else:
        return default_not_found, default_not_found, [], default_score


def consensus_similarity_flexible(results, main_weight=0.7, sub_weight=0.3):
    # Check if results is empty
    if not results:
        print("Warning: No results to calculate consensus similarity")
        return 0.0, "Unknown", "Unknown"
        
    general_types = [result[0] for result in results.values()]
    sub_types = [result[1] for result in results.values()]
    
    # Check if general_types or sub_types is empty
    if not general_types or not sub_types:
        print("Warning: Empty general_types or sub_types in consensus calculation")
        return 0.0, "Unknown", "Unknown"
    
    # Use Counter to get most common types
    from collections import Counter
    general_counter = Counter(general_types)
    sub_counter = Counter(sub_types)
    
    # Check if Counter is empty
    if not general_counter or not sub_counter:
        print("Warning: Empty counter in consensus calculation")
        return 0.0, "Unknown", "Unknown"
    
    consensus_general = general_counter.most_common(1)[0][0]
    consensus_sub = sub_counter.most_common(1)[0][0]
    
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






def agent_unification(prompt, system_prompt='''You are a careful professional biologist, specializing in single-cell RNA-seq analysis. You will be given a series of results from a celltype annotator. 

Your task is to unify all the celltype names, so that same celltype have the same name. Follow these rules EXACTLY:

1. The first letter of each word should be CAPITAL and others should be lowercase
2. Remove plurals (e.g., "cells" becomes "cell")  
3. Unify similar terms: "stem", "progenitor", and "immature" mean the same thing
4. Keep the exact same format as the input - only change the cell type names
5. ALWAYS wrap your output in <results></results> tags

CRITICAL: You MUST preserve the exact input format. If input is "result1:(cell1, cell2)" then output MUST be "result1:(Cell1, Cell2)"

Examples:

Input format:      
result1:(immune cell, t cell),result2:(Immune cells,t cell),result3:(T cell, cd8+ t cell)
                  
Output format:
<results>result1:(Immune Cell, T Cell),result2:(Immune Cell, T Cell),result3:(T Cell, Cd8+ T Cell)</results>

Another example:
                      
Input format:      
result1:(Hematopoietic stem/progenitor cells (HSPCs), T cell progenitors),result2:(Hematopoietic Progenitor cells,t cell),result3:(Hematopoietic progenitor cells, T cell)
                  
Output format:
<results>result1:(Hematopoietic Progenitor Cell, T Cell Progenitor),result2:(Hematopoietic Progenitor Cell, T Cell),result3:(Hematopoietic Progenitor Cell, T Cell)</results>

REMEMBER: Always use <results></results> tags around your output and maintain the exact same structure as the input.
''', model=None, provider="openai", temperature=0, deplural_only=False):
    """
    Unified function to call any LLM provider for cell type name unification.
    
    Args:
        prompt: The unification prompt containing cell type names to unify
        system_prompt: Instructions for the LLM
        model: Model to use (defaults to provider's default if None)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom URL)
        temperature: Temperature for generation (0-1)
        deplural_only: If True, only removes plurals without other unification
    
    Returns:
        The unified cell type names as a string
    """
    # If deplural_only is True, override the system prompt
    if deplural_only:
        system_prompt = "Remove the plural for celltype name, keep the original input format."
    
    # Select default model based on provider if not specified
    if model is None:
        if provider == "openai":
            model = "gpt-4o"
        elif provider == "anthropic":
            model = "claude-3-5-sonnet-20241022"
        elif provider == "openrouter":
            model = "anthropic/claude-3.5-sonnet"
        elif provider.startswith("http"):
            # For custom API endpoints, use a default model if none specified
            model = model or "deepseek-chat"
    
    # Call LLM using the unified function
    result = call_llm(
        prompt=prompt,
        provider=provider,
        model=model,
        temperature=temperature,
        system_prompt=system_prompt
    )
    
    # Process the result for all providers to extract content from <results> tags
    # This ensures consistent behavior regardless of provider
    if not deplural_only:
        import re
        results_match = re.search(r'<results>(.*?)</results>', result, re.DOTALL)
        if results_match:
            return results_match.group(1)
    
    return result







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
    , model=None, provider="openai", temperature=0, single_analysis=False, use_json_tags=False):
    """
    Unified function to call any LLM provider for cell type judgment.
    
    Args:
        prompt: The prompt containing cell type info for judgment
        system_prompt: Instructions for the LLM
        model: Model to use (defaults to provider's default if None)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom URL)
        temperature: Temperature for generation (0-1)
        single_analysis: Whether this is a single analysis judgment 
        use_json_tags: Whether to add <json></json> tags around the output
        
    Returns:
        The judgment result as a string
    """
    # Select default model based on provider if not specified
    if model is None:
        if provider == "openai":
            model = "gpt-4o"
        elif provider == "anthropic":
            model = "claude-3-5-sonnet-20241022"
        elif provider == "openrouter":
            model = "anthropic/claude-3.5-sonnet"
        elif provider.startswith("http"):
            # For custom API endpoints, use a default model if none specified
            model = model or "deepseek-chat"
            
    # Modify system prompt for all providers to include JSON tags
    # This helps standardize the output format across all providers
    if use_json_tags or True:  # Always use JSON tags for consistency
        # Replace the JSON format instruction to include tags
        if "Output in JSON format:" in system_prompt:
            system_prompt = system_prompt.replace(
                "Output in JSON format:", 
                "Output in JSON format with tags:"
            ).replace(
                "consensus_score\": 0-100\n}", 
                "consensus_score\": 0-100\n}</json>"
            ).replace(
                "{\n", 
                "<json>{\n"
            )
    
    # Call LLM using the unified function
    return call_llm(
        prompt=prompt,
        provider=provider,
        model=model,
        temperature=temperature,
        system_prompt=system_prompt
    )





def agent_unification_deplural(prompt, model="gpt-4o", provider="openai", temperature=0):
    """
    Wrapper around agent_unification that only removes plurals from cell type names.
    
    Args:
        prompt: The prompt containing cell type names to depluralize
        model: Model to use
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom URL)
        temperature: Temperature for generation (0-1)
        
    Returns:
        The depluralized cell type names as a string
    """
    return agent_unification(
        prompt=prompt,
        model=model,
        provider=provider,
        temperature=temperature,
        deplural_only=True
    )



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
    """
    Robust standardization function that tries multiple parsing approaches.
    If parsing fails, it returns the original input_string.
    """
    # Remove all hyphens from the input string
    current_input_string = input_string.replace("-", " ")
    
    results = []
    parsing_successful = False

    # Try multiple regex patterns to extract results
    patterns = [
        # Original strict pattern with single quotes
        r"result\d+:\('([^']+)', '([^']+)'\)",
        # Pattern with double quotes
        r"result\d+:\(\"([^\"]+)\", \"([^\"]+)\"\)",
        # Pattern without quotes
        r"result\d+:\(([^,]+), ([^)]+)\)",
        # Case insensitive result pattern
        r"(?i)result\s*\d+\s*:\s*[\(\[]([^,\]\)]+),\s*([^\]\)]+)[\)\]]",
        # Very flexible pattern
        r"(?i)\d+[^\(\[]*[\(\[]([^,\]\)]+),\s*([^\]\)]+)[\)\]]"
    ]
    
    # Try each pattern until we find matches
    for pattern_idx, pattern in enumerate(patterns):
        results = re.findall(pattern, current_input_string)
        if results:
            print(f"Standardization: Successfully parsed using pattern {pattern_idx + 1}")
            parsing_successful = True
            break
    
    # If no pattern worked, try to use the robust parser
    if not parsing_successful:
        print(f"Standardization: No direct patterns matched, trying robust parser...")
        parsed_dict = parse_results_to_dict(current_input_string)
        
        if parsed_dict:
            # Convert parsed dictionary back to tuples
            results = [(cell_types[0], cell_types[1]) for cell_types in parsed_dict.values()]
            print(f"Standardization: Extracted {len(results)} results using robust parser")
            parsing_successful = True

    # If parsing failed after all attempts
    if not parsing_successful: # or not results
        print("Warning: No results found to standardize.")
        print(f"Input string for standardization was: {input_string[:200]}..." if len(input_string) > 200 else input_string)
        return input_string # Return the original input_string if parsing failed
    
    standardized_results_list = []
    for i, (general_type, specific_type) in enumerate(results, 1):
        # Clean up the cell types before ontology lookup
        general_type = general_type.strip().strip("'\"")
        specific_type = specific_type.strip().strip("'\"")
        
        # Search for standardized names
        _, general_label = get_cell_type_info(general_type)
        _, specific_label = get_cell_type_info(specific_type)
        
        # Use original names if no standardized names found
        general_label = general_label or general_type
        specific_label = specific_label or specific_type
        
        standardized_results_list.append(f"result{i}:('{general_label}', '{specific_label}')")
    
    return ",".join(standardized_results_list)


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


def process_cell_type_variance_analysis_batch(results, model=None, provider="openai", temperature=0.0, main_weight=0.5, sub_weight=0.5, include_ontology=True):
    """
    Unified function to process cell type variance analysis for batch results.
    
    Args:
        results: The formatted batch results string
        model: Model to use (defaults to provider's default if None)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom URL)
        temperature: Temperature for generation (0-1)
        main_weight: Weight for the main cell type in similarity calculation
        sub_weight: Weight for the sub cell type in similarity calculation
        include_ontology: Whether to include ontology-based standardization (additional analysis)
        
    Returns:
        Dictionary containing processed results
    """
    
    # Default model based on provider if not specified
    if model is None:
        if provider == "openai":
            model = "gpt-4o"
        elif provider == "anthropic":
            model = "claude-3-5-sonnet-20241022"
        elif provider == "openrouter":
            model = "anthropic/claude-3.5-sonnet"
        elif provider.startswith("http"):
            # For custom API endpoints, use a default model if none specified
            model = model or "deepseek-chat"
    
    # Extract and format results using the unified agent_unification function
    try:
        results_unification_llm = agent_unification(
            prompt=results,
            model=model,
            provider=provider,
            temperature=temperature
        )
    except Exception as e:
        print(f"Error in agent_unification: {str(e)}")
        results_unification_llm = results  # Fall back to original results
    
    # Process with ontology if requested
    ontology_results_available = False
    if include_ontology:
        try:
            # Depluralize using unified function
            results_depluar = agent_unification_deplural(
                prompt=results, # This 'results' is the original full batch string
                model=model,
                provider=provider,
                temperature=temperature
            )
            
            # Attempt to standardize cell types
            result_unified_oncology = standardize_cell_types(results_depluar)
            results_depluar_retry = None # Initialize results_depluar_retry to None
            
            # Check if standardization failed (returned input string unchanged)
            # and if the input string was non-trivial
            if result_unified_oncology == results_depluar and results_depluar and results_depluar.strip():
                print("Standardization failed on first attempt. Attempting re-generation of depluralized string and retrying standardization...")
                # Re-run depluralization
                results_depluar_retry = agent_unification_deplural(
                    prompt=results, # Use the same original prompt
                    model=model,
                    provider=provider,
                    temperature=temperature
                )
                # Retry standardization with the new depluralized string
                result_unified_oncology = standardize_cell_types(results_depluar_retry)
                # If it still fails, result_unified_oncology will be results_depluar_retry. 
                # The subsequent logic will handle it as a failed standardization.

            # Check if any actual standardization occurred or if it's still the raw depluralized string
            if result_unified_oncology != results_depluar and \
               (results_depluar_retry is None or result_unified_oncology != results_depluar_retry):
                 ontology_results_available = True
            else:
                 # This means even after a potential retry, standardization didn't change the string,
                 # or the initial standardization worked but produced the same string (e.g., already standardized)
                 # or it failed and returned the input string.
                 # We should consider it not truly "available" if it's identical to the unstandardized depluralized form
                 # unless results_depluar itself was already in the final standardized format.
                 # A simple check is if it contains "result1:(" which is our target format.
                 if "result1:(" in result_unified_oncology:
                     ontology_results_available = True
                 else:
                     ontology_results_available = False
                     print("Ontology standardization did not produce a parseable result or made no changes.")

        except Exception as e:
            print(f"Error in ontology processing: {str(e)}")
            result_unified_oncology = results  # Fall back to original results
    
    # Consensus judgment using unified function - always use JSON tags for consistency
    try:
        result_consensus_from_llm = agent_judgement(
            prompt=results_unification_llm,
            model=model,
            provider=provider,
            temperature=temperature,
            use_json_tags=True  # Always use JSON tags for all providers
        )
    except Exception as e:
        print(f"Error in agent_judgement: {str(e)}")
        result_consensus_from_llm = "{}"  # Empty JSON as fallback
    
    # Extract consensus celltypes using unified function
    try:
        general_celltype, sub_celltype, mixed_types, llm_generated_consensus_score_llm = extract_celltypes_from_llm(
            llm_response=result_consensus_from_llm,
            provider=provider
        )
    except Exception as e:
        print(f"Error extracting celltypes: {str(e)}")
        general_celltype, sub_celltype, mixed_types, llm_generated_consensus_score_llm = "Unknown", "Unknown", [], 0
    
    # Calculate similarity score
    try:
        parsed_results_llm = parse_results_to_dict(results_unification_llm)
        consensus_score_llm, consensus_1_llm, consensus_2_llm = consensus_similarity_flexible(
            parsed_results_llm, 
            main_weight=main_weight, 
            sub_weight=sub_weight
        )
    except Exception as e:
        print(f"Error calculating similarity score: {str(e)}")
        consensus_score_llm, consensus_1_llm, consensus_2_llm = 0, "Unknown", "Unknown"
    
    # Build result dictionary
    result = {
        'general_celltype_llm': general_celltype,
        'sub_celltype_llm': sub_celltype,
        'mixed_celltypes_llm': mixed_types,
        'consensus_score_llm': consensus_score_llm,
        'llm_generated_consensus_score_llm': llm_generated_consensus_score_llm,
        'count_consensus_1_llm': consensus_1_llm,
        'count_consensus_2_llm': consensus_2_llm,
        'unified_results_llm': results_unification_llm,
        'result_consensus_from_llm': result_consensus_from_llm
    }
    
    # Add ontology-based results if included and available
    if include_ontology and ontology_results_available:
        try:
            # Get consensus from ontology-standardized results - always use JSON tags
            result_consensus_from_oncology = agent_judgement(
                prompt=result_unified_oncology,
                model=model,
                provider=provider,
                temperature=temperature,
                use_json_tags=True  # Always use JSON tags for all providers
            )
            
            # Extract ontology-based consensus
            general_celltype_oncology, sub_celltype_oncology, mixed_types_oncology, llm_generated_consensus_score_oncology = extract_celltypes_from_llm(
                llm_response=result_consensus_from_oncology,
                provider=provider
            )
            
            # Calculate ontology-based similarity score
            parsed_results_oncology = parse_results_to_dict(result_unified_oncology)
            consensus_score_oncology, consensus_1_oncology, consensus_2_oncology = consensus_similarity_flexible(
                parsed_results_oncology, 
                main_weight=main_weight, 
                sub_weight=sub_weight
            )
            
            # Add ontology-based results to the result dictionary
            result.update({
                'general_celltype_oncology': general_celltype_oncology,
                'sub_celltype_oncology': sub_celltype_oncology,
                'mixed_types_oncology': mixed_types_oncology,
                'consensus_score_oncology': consensus_score_oncology,
                'llm_generated_consensus_score_oncology': llm_generated_consensus_score_oncology,
                'count_consensus_1_oncology': consensus_1_oncology,
                'count_consensus_2_oncology': consensus_2_oncology,
                'unified_results_oncology': result_unified_oncology,
                'result_consensus_from_oncology': result_consensus_from_oncology
            })
        except Exception as e:
            print(f"Error in ontology-based consensus processing: {str(e)}")
            # Add default values for ontology results
            result.update({
                'general_celltype_oncology': "Error in processing",
                'sub_celltype_oncology': "Error in processing",
                'mixed_types_oncology': [],
                'consensus_score_oncology': 0,
                'llm_generated_consensus_score_oncology': 0,
                'count_consensus_1_oncology': "Unknown",
                'count_consensus_2_oncology': "Unknown",
                'unified_results_oncology': result_unified_oncology,
                'result_consensus_from_oncology': "{}"
            })
    
    return result


def process_cell_type_results(organized_results, max_workers=10, model="google/gemini-2.5-flash-preview", provider="openrouter", main_weight=0.5, sub_weight=0.5, temperature=0.0):
    processed_results = {}
    
    def process_single_celltype(celltype, predictions):
        # Filter out 'N/A' predictions
        valid_predictions = [pred for pred in predictions if pred != ('N/A', 'N/A')]
        
        if not valid_predictions:
            return celltype, {
                'error': 'No valid predictions',
                'input': predictions
            }

        formatted_predictions = [f"result{i+1}:{pred}" for i, pred in enumerate(valid_predictions)]
        formatted_string = ",".join(formatted_predictions)

        # Use only the unified function
        result = process_cell_type_variance_analysis_batch(
            formatted_string,
            model=model,
            provider=provider,
            main_weight=main_weight,
            sub_weight=sub_weight,
            temperature=temperature
        )
        
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
    
    # Check if processed_results is empty
    if not processed_results:
        # Create an empty DataFrame with minimal columns
        df = pd.DataFrame(columns=['Cell Type', 'Status'])
        df.to_csv(output_csv, index=False)
        print(f"Empty results saved to {output_csv}")
        return df
    
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

    # Check if DataFrame is empty
    if df.empty:
        df.to_csv(output_csv, index=False)
        print(f"Empty results saved to {output_csv}")
        return df
        
    # Define expected columns
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
    
    # Filter fixed columns to only include those that exist in the DataFrame
    available_fixed_columns = [col for col in fixed_columns if col in df.columns]
    
    # Get original columns
    original_columns = [col for col in df.columns if col.startswith('Original') and col != 'Original Non-Unified Results']
    
    # Reorder columns if possible
    if available_fixed_columns:
        df = df[available_fixed_columns + original_columns]

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False)
    print(f"Results saved to {output_csv}")

    return df


def runCASSIA_similarity_score_batch(marker, file_pattern, output_name, celltype_column=None, max_workers=10, model="google/gemini-2.5-flash-preview", provider="openrouter", main_weight=0.5, sub_weight=0.5, temperature=0.0):
    """
    Process batch results and save them to a CSV file, measuring the time taken.

    Args:
    marker_file_path (str): Path to the marker file.
    file_pattern (str): Path to pattern of result files.
    output_csv_name (str): Name of the output CSV file.
    celltype_column (str): Name of the column containing cell types in the marker file.
    max_workers (int): Maximum number of workers for parallel processing.
    model (str): Model name to use for LLM calls.
    provider (str): LLM provider ("openai", "anthropic", "openrouter", or a custom URL).
    main_weight (float): Weight for the main cell type in similarity calculation.
    sub_weight (float): Weight for the sub cell type in similarity calculation.
    temperature (float): Temperature for the LLM calls.
    """

    # Organize batch results
    organized_results = organize_batch_results(
        marker=marker,
        file_pattern=file_pattern,
        celltype_column=celltype_column
    )

    # Process cell type results
    processed_results = process_cell_type_results(organized_results, max_workers=max_workers, model=model, provider=provider, main_weight=main_weight, sub_weight=sub_weight, temperature=temperature)

    # Create and save results dataframe
    create_and_save_results_dataframe(
        processed_results, 
        organized_results, 
        output_name=output_name
    )

    print(f"Similarity analysis completed: {output_name}")



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

def extract_celltypes_from_llm_single(llm_response, provider="openai"):
    """
    Wrapper around the unified extract_celltypes_from_llm function for single analysis.
    
    Args:
        llm_response: The text response from the LLM
        provider: The provider that generated the response ("openai", "anthropic", "openrouter", or a custom URL)
        
    Returns:
        general_celltype, sub_celltype, mixed_celltypes
    """
    return extract_celltypes_from_llm(
        llm_response=llm_response,
        provider=provider,
        single_analysis=True
    )

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

def agent_judgement_single(prompt, system_prompt, model="gpt-4o", provider="openai", temperature=0):
    """
    Wrapper around agent_judgement for single analysis operations.
    
    Args:
        prompt: The prompt containing cell type info for judgment
        system_prompt: Instructions for the LLM
        model: Model to use (defaults to gpt-4o for OpenAI)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom URL)
        temperature: Temperature for generation (0-1)
        
    Returns:
        The judgment result as a string
    """
    return agent_judgement(
        prompt=prompt, 
        system_prompt=system_prompt,
        model=model,
        provider=provider,
        temperature=temperature,
        single_analysis=True,
        use_json_tags=True  # Always use JSON tags for all providers
    )

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


def runCASSIA_n_times_similarity_score(tissue, species, additional_info, temperature, marker_list, model="google/gemini-2.5-flash-preview", max_workers=10, n=3, provider="openrouter",main_weight=0.5,sub_weight=0.5, validator_involvement="v1"):
    """
    Wrapper function for processing cell type analysis using any supported provider.
    
    Args:
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for analysis
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        model (str): Model name to use
        max_workers (int): Maximum number of parallel workers
        n (int): Number of analysis iterations
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a custom URL)
        main_weight (float): Weight for main cell type in similarity calculation
        sub_weight (float): Weight for sub cell type in similarity calculation
    
    Returns:
        dict: Analysis results including consensus types, cell types, and scores
    """
    # System prompt for all providers
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
    results = runCASSIA_n_times(n, tissue, species, additional_info, temperature, marker_list, model, max_workers=max_workers, provider=provider, validator_involvement=validator_involvement)
    results = extract_cell_types_from_results_single(results)
    
    # Standardize cell types
    standardized_results = standardize_cell_types_single(results)
    
    # Get consensus judgment using the unified function
    result_consensus = agent_judgement_single(
        prompt=standardized_results,
        system_prompt=system_prompt,
        model=model,
        provider=provider,
        temperature=0.3
    )
    
    # Extract consensus celltypes using the unified function
    general_celltype, sub_celltype, mixed_types = extract_celltypes_from_llm_single(
        result_consensus, provider=provider
    )
    
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
        'consensus_score_llm': None,  # Using None instead of consensus_score_llm which isn't returned by extract_celltypes_from_llm_single
        'similarity_score': consensus_score,
        'original_results': results
    }

