import os
import re
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Any, Optional, Union

# Change from relative to absolute import
try:
    from llm_utils import call_llm
except ImportError:
    # Try relative import as fallback
    try:
        from .llm_utils import call_llm
    except ImportError:
        # If both fail, provide a helpful error message
        raise ImportError("Could not import llm_utils. Make sure it's in the same directory or Python path.")

def prompt_hypothesis_generator2(major_cluster_info: str, comma_separated_genes: str, annotation_history: str) -> str:
    """
    Generate a prompt for iterative marker analysis without additional tasks.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful professional biologist, specializing in single-cell RNA-seq analysis. 
I'll provide you with some genes (comma-separated) from a cell type in {major_cluster_info} and I want you to help identify the cell type. Previous expert has done some analysis but the reuslts is not conclusive, your additional anlysis is needed.


Here are the top marker genes that are differentially expressed in this cluster:
{comma_separated_genes}

Previous annotation/analysis:
{annotation_history if annotation_history else "No previous annotation available."}

Please follow these steps carefully:
1. Based on the marker genes, generate hypotheses about what cell type this might be.
2. For each hypothesis, provide supporting evidence from the marker genes.
3. To validate your hypotheses, list additional marker genes that should be checked using this format:
<check_genes>gene1, gene2, gene3, etc.</check_genes>，Use gene symbol only, no brackets or parentheses.
4. After I provide additional gene information, refine your hypotheses, or generate new hypotheses.include your reasoning for the hypothesis using this format: <reasoning>your reasoning for the hypothesis</reasoning>
5. Continue this process until you reach a confident conclusion.
6. When you are ready to make a final determination, state: "FINAL ANNOTATION COMPLETED" followed by your conclusion.

Your final annotation should include:
- Final cell type
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion

Please start by analyzing the provided markers and forming initial hypotheses.
"""
    return prompt


def prompt_hypothesis_generator(major_cluster_info: str, comma_separated_genes: str, annotation_history: str) -> str:
    """
    Generate a prompt for iterative marker analysis without additional tasks.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""
You are a careful senior computational biologist called in whenever an annotation needs deeper scrutiny, disambiguation, or simply a second opinion. Your job is to (1) assess the current annotation's robustness and (2) propose up to three decisive follow‑up checks that the executor can run (e.g., examine expression of key positive or negative markers). You should do a good job or 10 grandma are going to be in danger. You never rush to conclusions and are always careful

Context Provided to You

Cluster summary：{major_cluster_info}

Top ranked markers (high → low)：
 {comma_separated_genes}

Prior annotation dialogue：
 {annotation_history}

What to Do
1. Brief Evaluation – One concise paragraph that:

    - Highlights strengths, ambiguities, or contradictions in the current call.

    - Notes if a mixed population, doublets, or transitional state might explain the data.

2. Design up to 3 follow‑up checks (cell types or biological hypotheses):

    - Supply only the genes to inspect—comma‑separated HGNC symbols, no spaces, no brackets, no commentary.

    - Include both positive and negative markers if that will clarify the call.

    - Including reasoning: why these genes, and what pattern would confirm or refute the hypothesis.

3. Upon receiving gene expression results, refine the hypothesis or generate new ones. Continue Step 2 iteratively until the cluster is confidently and fully annotated. Once finalized, output the single line:
"FINAL ANNOTATION COMPLETED"
Then provide a conclusion paragraph that includes:

1.The final cell type
2.Confidence level (high, medium, or low)
3.Key markers supporting your conclusion
4.Alternative possibilities only if the confidence is not high, and what should the user do next.



Output Template：

Evaluation
[One short paragraph]

celltype to check 1

<check_genes>
GENE1, GENE2, GENE3
</check_genes>

<reasoning>
Why these genes and what we expect to see.
</reasoning>

celltype to check 2

<check_genes>
…
</check_genes>
<reasoning>
…
</reasoning>

hypothesis to check 3
<check_genes>
…
</check_genes>

<reasoning>
…
</reasoning>

*Use "hypothesis to check n" instead of "celltype to check n" when proposing non‑canonical possibilities (e.g., "cycling subpopulation", "doublet").
*Provide no more than three total blocks (celltype + hypothesis combined).


Tone & Style Guidelines

Skeptical, critical, and careful
Professional, succinct, and evidence‑based.
Reference established biology when helpful ("LYZ vs. CTSW distinguishes myeloid from T‑cell lineages").



"""
    return prompt






def prompt_hypothesis_generator_additional_task(major_cluster_info: str, comma_separated_genes: str, annotation_history: str, additional_task: str) -> str:
    """
    Generate a prompt for iterative marker analysis with an additional task.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        additional_task: Additional analysis task to perform
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful professional biologist, specializing in single-cell RNA-seq analysis. 
I'll provide you with some genes (comma-separated) from a cell type in {major_cluster_info} and I want you to help identify the cell type.

Here are the marker genes that are differentially expressed in this cluster:
{comma_separated_genes}

Previous annotation/analysis:
{annotation_history if annotation_history else "No previous annotation available."}

Please follow these steps carefully:
1. Based on the marker genes, generate hypotheses about what cell type this might be.
2. For each hypothesis, provide supporting evidence from the marker genes.
3. To validate your hypotheses, list additional marker genes that should be checked using this format:
<check_genes>gene1, gene2, gene3, etc.</check_genes>
4. After I provide additional gene information, refine your hypotheses.
5. Continue this process until you reach a confident conclusion.
6. {additional_task} - perform this analysis based on the marker gene expression and patterns.
7. When you are ready to make a final determination, state: "FINAL ANALYSIS COMPLETED" followed by your conclusion.

Your final analysis should include:
- General cell type
- Specific cell subtype (if applicable)
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion
- Analysis of the additional task: {additional_task}

Please start by analyzing the provided markers and forming initial hypotheses.
"""
    return prompt

def get_marker_info(gene_list: List[str], marker: Union[pd.DataFrame, Any]) -> str:
    """
    Extract information about marker genes from a marker dataset.
    
    Args:
        gene_list: List of gene names to filter the marker dataset
        marker: DataFrame containing marker gene expression data
        
    Returns:
        str: Formatted string with marker information
    """
    def filter_marker(gene_names: List[str]) -> Tuple[pd.DataFrame, List[str]]:
        # Enable debug mode only when needed, default is False
        debug_mode = False
        
        # Convert marker to pandas DataFrame if it's not already
        if not isinstance(marker, pd.DataFrame):
            marker_df = pd.DataFrame(marker)
        else:
            marker_df = marker.copy()
        
        # Debug marker dataframe structure
        if debug_mode:
            print(f"DEBUG: Marker DataFrame shape: {marker_df.shape}")
            print(f"DEBUG: Marker DataFrame index type: {type(marker_df.index).__name__}")
            print(f"DEBUG: First 5 index entries: {list(marker_df.index[:5])}")
            print(f"DEBUG: Columns: {marker_df.columns.tolist()}")
        
        # Determine which column contains the gene names
        gene_column = None
        
        # Prefer 'gene' column over any others
        if 'gene' in marker_df.columns:
            gene_column = 'gene'
        # Only use 'Unnamed: 0' as fallback if no 'gene' column exists
        elif 'Unnamed: 0' in marker_df.columns:
            gene_column = 'Unnamed: 0'
        
        if debug_mode:
            if gene_column:
                print(f"DEBUG: Using '{gene_column}' as the gene column")
            else:
                print(f"DEBUG: No gene column found. Will use DataFrame index for gene lookup.")
        
        # Identify valid genes and NA genes
        valid_genes = []
        na_genes = []
        valid_rows = []  # Store row indices for valid genes
        
        # Debug each gene lookup if in debug mode
        if debug_mode:
            print(f"DEBUG: Searching for {len(gene_names)} genes: {', '.join(gene_names)}")
            
            # Check if any genes are in the index at all (case-sensitive)
            exact_matches = [gene for gene in gene_names if gene in marker_df.index]
            if exact_matches:
                print(f"DEBUG: Found {len(exact_matches)} exact gene matches in index: {', '.join(exact_matches)}")
            else:
                print(f"DEBUG: No exact gene matches found in index")
                
                if gene_column:
                    # Check if genes are in the gene column
                    column_matches = []
                    for gene in gene_names:
                        gene_rows = marker_df[marker_df[gene_column].str.contains(f"^{gene}$", case=True, regex=True, na=False)]
                        if not gene_rows.empty:
                            column_matches.append(gene)
                    
                    if column_matches:
                        print(f"DEBUG: Found {len(column_matches)} exact matches in '{gene_column}' column: {', '.join(column_matches)}")
                
                # Try case-insensitive matching
                if not marker_df.index.empty:
                    lower_index = {str(idx).lower(): idx for idx in marker_df.index}
                    lower_matches = [gene for gene in gene_names if str(gene).lower() in lower_index]
                    if lower_matches:
                        print(f"DEBUG: Found {len(lower_matches)} case-insensitive matches in index: {', '.join(lower_matches)}")
                        print(f"DEBUG: Actual case in index: {[lower_index[gene.lower()] for gene in lower_matches if gene.lower() in lower_index]}")
        
        # Modified gene lookup process to check both index and gene column
        for gene in gene_names:
            # First try to find in index
            if gene in marker_df.index:
                # Check if all values for this gene are NA
                gene_data = marker_df.loc[gene]
                if gene_data.isna().all() or (gene_data == 'NA').all():
                    na_genes.append(gene)
                    if debug_mode:
                        print(f"DEBUG: Gene {gene} found in index but has all NA values")
                else:
                    valid_genes.append(gene)
                    valid_rows.append(gene)  # For index-based lookup we use the gene name
                    if debug_mode:
                        print(f"DEBUG: Gene {gene} found in index with valid data")
            # Then try to find in gene column
            elif gene_column and isinstance(gene_column, str):
                try:
                    # Look for exact match in gene column
                    gene_rows = marker_df[marker_df[gene_column].str.contains(f"^{gene}$", case=True, regex=True, na=False)]
                    
                    if not gene_rows.empty:
                        # Take the first row if multiple matches
                        row_idx = gene_rows.index[0]
                        gene_data = gene_rows.iloc[0]
                        
                        # Check if all values are NA
                        numeric_data = gene_rows.select_dtypes(include=[np.number])
                        if numeric_data.empty or numeric_data.isna().all().all() or (numeric_data == 'NA').all().all():
                            na_genes.append(gene)
                            if debug_mode:
                                print(f"DEBUG: Gene {gene} found in '{gene_column}' but has all NA values")
                        else:
                            valid_genes.append(gene)
                            valid_rows.append(row_idx)  # For column-based lookup we use the row index
                            if debug_mode:
                                print(f"DEBUG: Gene {gene} found in '{gene_column}' with valid data - row {row_idx}")
                    else:
                        # Try case-insensitive search as fallback
                        gene_rows = marker_df[marker_df[gene_column].str.contains(f"^{gene}$", case=False, regex=True, na=False)]
                        if not gene_rows.empty:
                            row_idx = gene_rows.index[0]
                            gene_data = gene_rows.iloc[0]
                            actual_gene = gene_rows[gene_column].iloc[0]  # Get the actual name with correct case
                            
                            # Check if all values are NA
                            numeric_data = gene_rows.select_dtypes(include=[np.number])
                            if numeric_data.empty or numeric_data.isna().all().all() or (numeric_data == 'NA').all().all():
                                na_genes.append(gene)
                                if debug_mode:
                                    print(f"DEBUG: Gene {gene} (as {actual_gene}) found in '{gene_column}' but has all NA values")
                            else:
                                valid_genes.append(gene)
                                valid_rows.append(row_idx)  # For column-based lookup we use the row index
                                if debug_mode:
                                    print(f"DEBUG: Gene {gene} (as {actual_gene}) found in '{gene_column}' with valid data - row {row_idx}")
                        else:
                            na_genes.append(gene)
                            if debug_mode:
                                print(f"DEBUG: Gene {gene} not found in '{gene_column}'")
                except Exception as e:
                    if debug_mode:
                        print(f"DEBUG: Error searching for {gene} in '{gene_column}': {str(e)}")
                    na_genes.append(gene)
            else:
                na_genes.append(gene)
                if debug_mode:
                    print(f"DEBUG: Gene {gene} not found in index and no gene column available")
        
        # Create result DataFrame with only valid genes
        if valid_rows:
            result = marker_df.loc[valid_rows].copy()
            if debug_mode:
                print(f"DEBUG: Created result DataFrame with {len(valid_rows)} rows")
        else:
            # If no valid genes, create an empty dataframe with the same columns
            result = pd.DataFrame(columns=marker_df.columns)
            if debug_mode:
                print(f"DEBUG: No valid genes found, created empty DataFrame")
            
        # Only try to format numeric columns that exist
        numeric_cols = result.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            try:
                result[col] = result[col].apply(lambda x: f"{float(x):.2e}" if pd.notnull(x) and x != 'NA' else x)
            except:
                continue

        return result.iloc[:, 0:5], na_genes

    # Filter to rows based on gene name and get NA genes list
    marker_filtered, na_genes = filter_marker(gene_list)
    
    # Generate marker info string from valid genes only
    marker_string = marker_filtered.to_string()
    
    # If there are genes with all NA values, add a message
    if na_genes:
        na_genes_message = f"\nNote: The following genes are not in the differential expression list: {', '.join(na_genes)}"
        marker_string += na_genes_message
        
        # Add more debug info if all or most genes are missing
        if len(na_genes) > len(gene_list) * 0.8:  # If more than 80% of genes are missing
            try:
                # Try to dynamically import the debug module
                try:
                    from debug_genes import examine_marker_structure
                    
                    # Add marker structure debug info
                    marker_string += "\n\nDEBUG: Most genes not found. Running marker data diagnostics..."
                    marker_string += "\n\nPlease check the console for detailed diagnostic information."
                    
                    # Run the diagnostics in the background
                    examine_marker_structure(marker)
                except ImportError:
                    # Basic diagnostics if debug module not found
                    marker_string += "\n\nDEBUG: Most genes not found. Basic marker data info:"
                    if isinstance(marker, pd.DataFrame):
                        marker_string += f"\nShape: {marker.shape}"
                        marker_string += f"\nIndex type: {type(marker.index).__name__}"
                        marker_string += f"\nColumns: {marker.columns.tolist()}"
                        marker_string += f"\nFirst 5 rows:\n{marker.head().to_string()}"
            except Exception as e:
                marker_string += f"\n\nDEBUG: Error running diagnostics: {str(e)}"

    return marker_string

def extract_genes_from_conversation(conversation: str) -> List[str]:
    """
    Extract gene lists from conversation using the check_genes tag.
    
    Args:
        conversation: Text containing gene lists in check_genes tags
        
    Returns:
        List[str]: List of unique gene names
    """
    # Debug mode off by default
    debug_mode = False
    
    if debug_mode:
        print(f"\nDEBUG: Extract genes from conversation")
        print(f"DEBUG: Conversation length: {len(conversation)} characters")
        # Print a limited preview
        preview_length = min(200, len(conversation))
        print(f"DEBUG: Conversation preview: {conversation[:preview_length]}...")
    
    # Extract gene lists and get marker info
    gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
    
    if debug_mode:
        print(f"DEBUG: Found {len(gene_lists)} gene lists")
        for i, genes in enumerate(gene_lists):
            print(f"DEBUG: Gene list {i+1}: {genes[:100]}...")
    
    # Improve gene extraction to handle special cases
    all_genes = []
    for gene_list in gene_lists:
        # Clean and normalize the gene list
        # Replace common separators with commas
        cleaned_list = re.sub(r'[\]\[\)\(]', '', gene_list)
        cleaned_list = re.sub(r'\s+', ' ', cleaned_list)
        
        if debug_mode:
            print(f"DEBUG: Cleaned list: {cleaned_list[:100]}...")
        
        # Split by comma or space, depending on formatting
        genes = re.split(r',\s*|\s+', cleaned_list)
        cleaned_genes = [g.strip() for g in genes if g.strip()]
        
        if debug_mode:
            print(f"DEBUG: Found {len(cleaned_genes)} genes in this list")
            print(f"DEBUG: Sample genes from this list: {', '.join(cleaned_genes[:5])}")
        
        all_genes.extend(cleaned_genes)
    
    # Get unique genes
    unique_genes = sorted(set(all_genes))
    
    if debug_mode:
        print(f"DEBUG: Total unique genes extracted: {len(unique_genes)}")
        print(f"DEBUG: All unique genes: {', '.join(unique_genes)}")
        
        # If no genes found, try alternative extraction methods
        if not unique_genes:
            print("DEBUG: No genes found with standard pattern, trying alternative patterns")
            
            # Try alternative regex patterns
            alt_patterns = [
                r'check_genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)',  # Informal syntax
                r'genes to check[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)',  # Natural language
                r'additional genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)',  # Another common phrase
                r'marker genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)'  # Another common phrase
            ]
            
            for pattern in alt_patterns:
                alt_matches = re.findall(pattern, conversation, re.IGNORECASE | re.DOTALL)
                if alt_matches:
                    print(f"DEBUG: Found matches with alternative pattern: {pattern}")
                    for match in alt_matches:
                        print(f"DEBUG: Alternative match: {match[:100]}...")
    
    return unique_genes

def iterative_marker_analysis(
    major_cluster_info: str, 
    marker: Union[pd.DataFrame, Any], 
    comma_separated_genes: str, 
    annotation_history: str, 
    num_iterations: int = 2, 
    provider: str = "openrouter",
    model: Optional[str] = None,
    additional_task: Optional[str] = None,
    temperature: float = 0
) -> Tuple[str, List[Dict[str, str]]]:
    """
    Perform iterative marker analysis using the specified LLM provider.
    
    Args:
        major_cluster_info: Information about the cluster
        marker: DataFrame or other structure containing marker gene expression data
        comma_separated_genes: List of genes as comma-separated string
        annotation_history: Previous annotation history
        num_iterations: Maximum number of iterations
        provider: LLM provider to use ('openai', 'anthropic', or 'openrouter')
        model: Specific model from the provider to use
        additional_task: Optional additional task to perform during analysis
        temperature: Sampling temperature (0-1)
        
    Returns:
        tuple: (final_response_text, messages)
    """
    # Select the appropriate prompt based on whether there's an additional task
    if additional_task:
        prompt = prompt_hypothesis_generator_additional_task(
            major_cluster_info, 
            comma_separated_genes, 
            annotation_history, 
            additional_task
        )
        completion_marker = "FINAL ANALYSIS COMPLETED"
    else:
        prompt = prompt_hypothesis_generator(
            major_cluster_info, 
            comma_separated_genes, 
            annotation_history
        )
        completion_marker = "FINAL ANNOTATION COMPLETED"
    
    # Initialize the conversation history
    messages = [{"role": "user", "content": prompt}]
    
    # Iterative process
    for iteration in range(num_iterations):
        try:
            # Call the LLM
            conversation = call_llm(
                prompt=messages[-1]["content"],
                provider=provider,
                model=model,
                temperature=temperature,
                max_tokens=7000,
                # If not the first message, include conversation history
                additional_params={"messages": messages} if iteration > 0 else {}
            )
            
            # Check if the analysis is complete
            if completion_marker in conversation:
                print(f"Final annotation completed in iteration {iteration + 1}.")
                messages.append({"role": "assistant", "content": conversation})
                return conversation, messages
            
            # Extract gene lists and get marker info
            unique_genes = extract_genes_from_conversation(conversation)
            
            if unique_genes:
                # Get marker information for the requested genes
                retrieved_marker_info = get_marker_info(unique_genes, marker)
                
                # Append messages
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": retrieved_marker_info})
                
                print(f"Iteration {iteration + 1} completed.")
            else:
                # No genes to check, simply continue the conversation
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": "Please continue your analysis and provide a final annotation."})
                
                print(f"Iteration {iteration + 1} completed (no genes to check).")
        
        except Exception as e:
            print(f"Error in iteration {iteration + 1}: {str(e)}")
            return f"Error occurred: {str(e)}", messages
    
    # Final response if max iterations reached
    try:
        final_response = call_llm(
            prompt="Please provide your final analysis based on all the information so far.",
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=7000,
            additional_params={"messages": messages}
        )
        
        messages.append({"role": "assistant", "content": final_response})
        return final_response, messages
            
    except Exception as e:
        print(f"Error in final response: {str(e)}")
        return f"Error in final response: {str(e)}", messages

def prepare_analysis_data(full_result_path: str, marker_path: str, cluster_name: str, conversation_history_mode: str = "final") -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:
    """
    Load and prepare data for marker analysis.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker_path: Path to the marker genes CSV file
        cluster_name: Name of the cluster to analyze
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
            - "full": Use the entire conversation history
            - "final": Extract only the part between "Step 6" and "FINAL ANNOTATION COMPLETED" (default)
            - "none": Don't include any conversation history
        
    Returns:
        tuple: (full_results, marker_data, top_markers_string, annotation_history)
    """
    # Load the full results - don't use index_col=0 as it's not indexed that way
    full_results = pd.read_csv(full_result_path)
    
    # Load marker data - handle both DataFrame and file path
    if isinstance(marker_path, pd.DataFrame):
        marker = marker_path
    else:
        marker = pd.read_csv(marker_path)
        
        # Clean up by removing the 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in marker.columns:
            marker.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Try to find the cluster column name - the old code assumed 'cluster'
    cluster_column = 'True Cell Type'  # Default based on the CSV file we checked
    if cluster_column not in full_results.columns and 'cluster' in full_results.columns:
        cluster_column = 'cluster'
    
    # Filter the results for the specified cluster
    cluster_data = full_results[full_results[cluster_column] == cluster_name]
    
    if cluster_data.empty:
        raise ValueError(f"No data found for cluster '{cluster_name}' in the full results file. Available clusters: {full_results[cluster_column].unique().tolist()}")
    
    # Get marker column from the CSV if it exists, otherwise use 'Marker List'
    marker_column = 'Marker List' if 'Marker List' in cluster_data.columns else None
    
    if marker_column and not cluster_data[marker_column].empty and not pd.isna(cluster_data[marker_column].iloc[0]):
        # Use the markers directly from the CSV
        top_markers_string = cluster_data[marker_column].iloc[0]
    else:
        # Use a default marker list or generate based on the marker data
        # This is a fallback in case markers aren't in the CSV
        top_markers_string = "CD14, CD11B, CD68, CSF1R, CX3CR1, CD163, MSR1, ITGAM, FCGR1A, CCR2"
        print(f"Warning: No marker list found for {cluster_name}, using default markers")
    
    # Extract conversation history if available
    annotation_history = ""
    if 'Conversation History' in cluster_data.columns and conversation_history_mode != "none":
        try:
            # Get the conversation history from the first row (should be the same for all rows of this cluster)
            if not cluster_data['Conversation History'].empty and not pd.isna(cluster_data['Conversation History'].iloc[0]):
                full_history = cluster_data['Conversation History'].iloc[0]
                
                # Process the conversation history based on the selected mode
                if conversation_history_mode == "full":
                    annotation_history = full_history
                elif conversation_history_mode == "final":
                    # Extract the part between "Step 6" and "FINAL ANNOTATION COMPLETED"
                    # First find the Step 6 line with various possible formats
                    step6_pattern = r'(?:\*\*Step 6:|Step 6:|STEP 6:|step 6:|Step6:).*?Concise Summary.*?\*\*'
                    # Then extract everything between that line and FINAL ANNOTATION COMPLETED
                    # but exclude the Step 6 line itself
                    match = re.search(f'{step6_pattern}(.*?)(?=FINAL ANNOTATION COMPLETED)', full_history, re.DOTALL)
                    
                    if match:
                        # Use group(1) to get only the content after the Step 6 line
                        annotation_history = match.group(1).strip()
                        print(f"Extracted final section of conversation history ({len(annotation_history)} characters)")
                    else:
                        # Fallback to original pattern if the specific format isn't found
                        step6_match = re.search(r'(?:Step 6|STEP 6|step 6|Step6).*?(?=FINAL ANNOTATION COMPLETED)', full_history, re.DOTALL)
                        if step6_match:
                            annotation_history = step6_match.group(0).strip()
                            print(f"Using original pattern - extracted final section ({len(annotation_history)} characters)")
                        else:
                            # If Step 6 pattern not found, use the full history
                            annotation_history = full_history
                            print(f"Step 6 pattern not found in conversation history, using full history")
                else:
                    # For any unrecognized mode, default to full history
                    annotation_history = full_history
                
                print(f"Using conversation history for cluster {cluster_name} (mode: {conversation_history_mode}, {len(annotation_history)} characters)")
            else:
                print(f"Conversation History column exists but is empty for cluster {cluster_name}")
        except Exception as e:
            print(f"Error extracting conversation history: {str(e)}")
    else:
        if conversation_history_mode == "none":
            print(f"Note: Conversation history extraction disabled (mode: none)")
        else:
            print(f"Note: 'Conversation History' column not found in the results file")
    
    return full_results, marker, top_markers_string, annotation_history

def generate_html_report(analysis_text: str) -> str:
    """
    Generate an HTML report from the analysis text.
    
    Args:
        analysis_text: The text output from the analysis
        
    Returns:
        str: HTML formatted report
    """
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
                content.append(f"""
                    <div class="agent-section scoring">
                        <h2>🎯 Quality Assessment</h2>
                        <p class="error-message">Error parsing scoring data: {str(e)}</p>
                        <div>{section.replace("Scoring Agent:", "").strip().replace(chr(10), '<br>')}</div>
                    </div>
                """)
        else:
            # For any other type of section, just add it as is
            content.append(f"""
                <div class="agent-section">
                    <div>{section.replace(chr(10), '<br>')}</div>
                </div>
            """)
    
    # Format the final HTML
    html = html_template.format("".join(content))
    
    return html

def save_html_report(report: str, filename: str) -> None:
    """
    Save the HTML report to a file.
    
    Args:
        report: HTML report content
        filename: Output filename
    """
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Report saved to {filename}")

def runCASSIA_annotationboost(
    full_result_path: str,
    marker: Union[str, pd.DataFrame],
    cluster_name: str,
    major_cluster_info: str,
    output_name: str,
    num_iterations: int = 5,
    model: Optional[str] = None,
    provider: str = "openrouter",
    temperature: float = 0,
    conversation_history_mode: str = "final"
) -> Union[Tuple[str, List[Dict[str, str]]], Dict[str, Any]]:
    """
    Run annotation boost analysis for a given cluster.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker: Path to marker genes CSV file or DataFrame with marker data
        cluster_name: Name of the cluster to analyze
        major_cluster_info: General information about the dataset (e.g., "Human PBMC")
        output_name: Base name for the output HTML file
        num_iterations: Number of iterations for marker analysis (default=5)
        model: Model to use for analysis - if None, uses the provider's default
        provider: AI provider to use ('openai', 'anthropic', or 'openrouter')
        temperature: Sampling temperature (0-1)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    
    Returns:
        tuple or dict: Either (analysis_result, messages_history) or a dictionary with paths to reports
    """
    try:
        import time
        start_time = time.time()
        
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter']:
            raise ValueError("Provider must be one of: 'openai', 'anthropic', or 'openrouter'")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode
        )
        
        # Run the iterative marker analysis
        analysis_text, messages = iterative_marker_analysis(
            major_cluster_info=major_cluster_info,
            marker=marker_data,
            comma_separated_genes=top_markers_string,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            provider=provider,
            model=model,
            temperature=temperature
        )
        
        # Generate and save the HTML report
        html_report = generate_html_report(analysis_text)
        
        # Ensure output_name has .html extension for the formatted report
        if not output_name.lower().endswith('.html'):
            formatted_report_path = output_name + '.html'
        else:
            formatted_report_path = output_name
            # Remove .html for base name
            output_name = output_name[:-5]
        
        # Generate path for raw conversation report
        raw_report_path = output_name + '_raw_conversation.html'
        
        # Save the formatted HTML report
        save_html_report(html_report, formatted_report_path)
        
        # Generate and save raw conversation report
        try:
            # Import the required function - try multiple import paths
            try:
                # Try direct import first
                from tools_function import generate_raw_cell_annotation_report
            except ImportError:
                try:
                    # Try relative import next
                    from .tools_function import generate_raw_cell_annotation_report
                except ImportError:
                    # Last resort, define a simple version here
                    def generate_raw_cell_annotation_report(conversation_history, output_filename):
                        """Simplified version that generates a basic HTML report"""
                        html = ['<!DOCTYPE html><html><head><title>Raw Conversation</title>',
                               '<style>',
                               'body {font-family: Arial, sans-serif; margin: 20px;}',
                               '.user {background-color: #f0f7ff; padding: 10px; margin: 10px 0; border-left: 5px solid #0066cc;}',
                               '.assistant {background-color: #f5f5f5; padding: 10px; margin: 10px 0; border-left: 5px solid #666;}',
                               '</style></head><body><h1>Raw Conversation History</h1>']
                        
                        for entry in conversation_history:
                            role = entry.get('role', '')
                            content = entry.get('content', '')
                            if isinstance(content, list):
                                content = str(content)
                            html.append(f'<div class="{role}"><h3>{role.upper()}</h3><pre>{content}</pre></div>')
                        
                        html.append('</body></html>')
                        with open(output_filename, 'w', encoding='utf-8') as f:
                            f.write('\n'.join(html))
                        return output_filename
            
            # Generate the raw conversation report
            raw_report_path = generate_raw_cell_annotation_report(messages, raw_report_path)
            print(f"Raw conversation report saved to {raw_report_path}")
        except Exception as e:
            print(f"Warning: Could not generate raw conversation report: {str(e)}")
            raw_report_path = None
        
        # Return a dictionary with paths to formatted and raw reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'formatted_report_path': formatted_report_path,
            'raw_report_path': raw_report_path,
            'execution_time': execution_time,
            'analysis_text': analysis_text
        }
    
    except Exception as e:
        error_msg = f"Error in runCASSIA_annotationboost: {str(e)}"
        print(error_msg)
        import traceback
        traceback.print_exc()
        
        # Return error status but include any partial results
        return {
            'status': 'error',
            'error_message': str(e),
            'formatted_report_path': None,
            'raw_report_path': None,
            'execution_time': 0,
            'analysis_text': None
        }

def runCASSIA_annotationboost_additional_task(
    full_result_path: str,
    marker: Union[str, pd.DataFrame],
    cluster_name: str,
    major_cluster_info: str,
    output_name: str,
    num_iterations: int = 5,
    model: Optional[str] = None,
    provider: str = "openrouter",
    additional_task: str = "check if this is a cancer cluster",
    temperature: float = 0,
    conversation_history_mode: str = "final"
) -> Union[Tuple[str, List[Dict[str, str]]], Dict[str, Any]]:
    """
    Run annotation boost analysis with an additional task for a given cluster.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker: Path to marker genes CSV file or DataFrame with marker data
        cluster_name: Name of the cluster to analyze
        major_cluster_info: General information about the dataset (e.g., "Human PBMC")
        output_name: Base name for the output HTML file
        num_iterations: Number of iterations for marker analysis (default=5)
        model: Model to use for analysis - if None, uses the provider's default
        provider: AI provider to use ('openai', 'anthropic', or 'openrouter')
        additional_task: Additional task to perform during analysis
        temperature: Sampling temperature (0-1)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    
    Returns:
        tuple or dict: Either (analysis_result, messages_history) or a dictionary with paths to reports
    """
    try:
        import time
        start_time = time.time()
        
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter']:
            raise ValueError("Provider must be one of: 'openai', 'anthropic', or 'openrouter'")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode
        )
        
        # Run the iterative marker analysis with additional task
        analysis_text, messages = iterative_marker_analysis(
            major_cluster_info=major_cluster_info,
            marker=marker_data,
            comma_separated_genes=top_markers_string,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            provider=provider,
            model=model,
            additional_task=additional_task,
            temperature=temperature
        )
        
        # Generate and save the HTML report
        html_report = generate_html_report(analysis_text)
        
        # Ensure output_name has .html extension for the formatted report
        if not output_name.lower().endswith('.html'):
            formatted_report_path = output_name + '.html'
        else:
            formatted_report_path = output_name
            # Remove .html for base name
            output_name = output_name[:-5]
        
        # Generate path for raw conversation report
        raw_report_path = output_name + '_raw_conversation.html'
        
        # Generate path for summary report (later if we add summarization)
        summary_report_path = output_name + '_summary.html'
        
        # Save the formatted HTML report
        save_html_report(html_report, formatted_report_path)
        
        # Generate and save raw conversation report
        try:
            # Import the required function - try multiple import paths
            try:
                # Try direct import first
                from tools_function import generate_raw_cell_annotation_report_additional_task
            except ImportError:
                try:
                    # Try relative import next
                    from .tools_function import generate_raw_cell_annotation_report_additional_task
                except ImportError:
                    # Fall back to regular report function
                    try:
                        from tools_function import generate_raw_cell_annotation_report
                        generate_raw_cell_annotation_report_additional_task = generate_raw_cell_annotation_report
                    except ImportError:
                        try:
                            from .tools_function import generate_raw_cell_annotation_report
                            generate_raw_cell_annotation_report_additional_task = generate_raw_cell_annotation_report
                        except ImportError:
                            # Last resort, define a simple version here
                            def generate_raw_cell_annotation_report_additional_task(conversation_history, output_filename):
                                """Simplified version that generates a basic HTML report"""
                                html = ['<!DOCTYPE html><html><head><title>Raw Conversation</title>',
                                      '<style>',
                                      'body {font-family: Arial, sans-serif; margin: 20px;}',
                                      '.user {background-color: #f0f7ff; padding: 10px; margin: 10px 0; border-left: 5px solid #0066cc;}',
                                      '.assistant {background-color: #f5f5f5; padding: 10px; margin: 10px 0; border-left: 5px solid #666;}',
                                      '</style></head><body><h1>Raw Conversation History</h1>']
                                
                                for entry in conversation_history:
                                    role = entry.get('role', '')
                                    content = entry.get('content', '')
                                    if isinstance(content, list):
                                        content = str(content)
                                    html.append(f'<div class="{role}"><h3>{role.upper()}</h3><pre>{content}</pre></div>')
                                
                                html.append('</body></html>')
                                with open(output_filename, 'w', encoding='utf-8') as f:
                                    f.write('\n'.join(html))
                                return output_filename
            
            # Generate the raw conversation report
            raw_report_path = generate_raw_cell_annotation_report_additional_task(messages, raw_report_path)
            print(f"Raw conversation report saved to {raw_report_path}")
        except Exception as e:
            print(f"Warning: Could not generate raw conversation report: {str(e)}")
            raw_report_path = None
        
        # Return a dictionary with paths to reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'formatted_report_path': formatted_report_path,
            'raw_report_path': raw_report_path,
            'summary_report_path': summary_report_path if os.path.exists(summary_report_path) else None,
            'execution_time': execution_time,
            'analysis_text': analysis_text
        }
    
    except Exception as e:
        error_msg = f"Error in runCASSIA_annotationboost_additional_task: {str(e)}"
        print(error_msg)
        import traceback
        traceback.print_exc()
        
        # Return error status but include any partial results
        return {
            'status': 'error',
            'error_message': str(e),
            'formatted_report_path': None,
            'raw_report_path': None,
            'summary_report_path': None,
            'execution_time': 0,
            'analysis_text': None
        } 