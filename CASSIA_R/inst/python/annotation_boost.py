import os
import re
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Any, Optional, Union

# Change from relative to absolute import
from llm_utils import call_llm


def prompt_hypothesis_generator2(
    major_cluster_info: str,
    comma_separated_genes: str,
    annotation_history: str) -> str:
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


def prompt_hypothesis_generator(
    major_cluster_info: str,
    comma_separated_genes: str,
    annotation_history: str) -> str:
    """
    Generate a prompt for iterative marker analysis without additional tasks.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful senior computational biologist called in whenever an annotation needs deeper scrutiny, disambiguation, or simply a second opinion. Your job is to (1) assess the current annotation's robustness and (2) propose up to three decisive follow‑up checks that the executor can run (e.g., examine expression of key positive or negative markers). You should do a good job or 10 grandma are going to be in danger. You never rush to conclusions and are always careful

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

2. Design up to 3 follow‑up checks (cell types or biological hypotheses):

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
*For each hypothesis check no more than 7 genes.


Tone & Style Guidelines

Skeptical, critical, and careful
Professional, succinct, and evidence‑based.
Reference established biology when helpful ("LYZ vs. CTSW distinguishes myeloid from T‑cell lineages").


"""
    return prompt


def prompt_hypothesis_generator_additional_task(
    major_cluster_info: str,
    comma_separated_genes: str,
    annotation_history: str,
    additional_task: str) -> str:
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


def get_marker_info(gene_list: List[str],
                    marker: Union[pd.DataFrame, Any]) -> str:
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
            print(
                f"DEBUG: Marker DataFrame index type: {type(marker_df.index).__name__}")
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
                print(
                    f"DEBUG: No gene column found. Will use DataFrame index for gene lookup.")

        # Identify valid genes and NA genes
        valid_genes = []
        na_genes = []
        valid_rows = []  # Store row indices for valid genes
        gene_name_map = {}  # Map of row index to gene name for output

        if debug_mode:
            print(
                f"DEBUG: Searching for {len(gene_names)} genes: {', '.join(gene_names)}")

            # Check if any genes are in the index at all (case-sensitive)
            exact_matches = [
                gene for gene in gene_names if gene in marker_df.index]
            if exact_matches:
                print(
                    f"DEBUG: Found {len(exact_matches)} exact gene matches in index: {', '.join(exact_matches)}")
            else:
                print(f"DEBUG: No exact gene matches found in index")

            if gene_column:
                # Check if genes are in the gene column
                column_matches = []
                for gene in gene_names:
                    gene_rows = marker_df[marker_df[gene_column].str.contains(
                        f"^{gene}$", case=True, regex=True, na=False)]
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
                gene_name_map[gene] = gene  # Map the index to gene name
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
                        # gene_data = gene_rows.iloc[0] # This line seems unused, consider removing if truly not needed

                        # Check if all values are NA
                        numeric_data = gene_rows.select_dtypes(include=[np.number])
                        if numeric_data.empty or numeric_data.isna().all().all() or (numeric_data == 'NA').all().all():
                            na_genes.append(gene)
                            if debug_mode:
                                print(f"DEBUG: Gene {gene} found in '{gene_column}' but has all NA values")
                        else:
                            valid_genes.append(gene)
                            valid_rows.append(row_idx)  # For column-based lookup we use the row index
                            gene_name_map[row_idx] = gene  # Map the index to gene name
                            if debug_mode:
                                print(f"DEBUG: Gene {gene} found in '{gene_column}' with valid data - row {row_idx}")
                    else:
                        # Try case-insensitive search as fallback
                        gene_rows = marker_df[marker_df[gene_column].str.contains(f"^{gene}$", case=False, regex=True, na=False)]
                        if not gene_rows.empty:
                            row_idx = gene_rows.index[0]
                            # gene_data = gene_rows.iloc[0] # This line seems unused
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
                                gene_name_map[row_idx] = gene  # Map the index to gene name
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
            
            # Add gene name as a column if it's not the index
            if result.index.name != 'gene' and 'gene' not in result.columns:
                result['gene'] = [gene_name_map.get(idx, str(idx)) for idx in result.index]
                # Move gene column to first position
                cols = result.columns.tolist()
                cols.insert(0, cols.pop(cols.index('gene')))
                result = result[cols]

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
            except BaseException:
                continue

        return result, na_genes

    # Filter to rows based on gene name and get NA genes list
    marker_filtered, na_genes = filter_marker(gene_list)
    
    # Ensure gene names are visible in output by making a copy with gene as first column if needed
    if 'gene' not in marker_filtered.columns and marker_filtered.index.name != 'gene':
        # Add the index as a column named 'gene'
        output_df = marker_filtered.reset_index()
        # Rename index column to 'gene' if it has no name
        if output_df.columns[0] == 'index':
            output_df = output_df.rename(columns={'index': 'gene'})
    else:
        output_df = marker_filtered
    
    # Remove the 'cluster' column if it exists - it's not needed in the gene expression output
    if 'cluster' in output_df.columns:
        output_df = output_df.drop(columns=['cluster'])
    
    # Ensure 'gene' column is the first column
    if 'gene' in output_df.columns and list(output_df.columns).index('gene') > 0:
        cols = list(output_df.columns)
        cols.remove('gene')
        cols.insert(0, 'gene')
        output_df = output_df[cols]
        
    # Generate marker info string from valid genes only - don't show the row indices
    marker_string = output_df.to_string(index=False)
    
    # If there are genes with all NA values, add a message
    if na_genes:
        na_genes_message = f"\nNote: The following genes are not in the differential expression list: {', '.join(na_genes)}"
        marker_string += na_genes_message

        # Add more debug info if all or most genes are missing
        if len(na_genes) > len(gene_list) * 0.8:  # If more than 80% of genes are missing
            try:
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


def prepare_analysis_data(full_result_path: str, marker_path: str, cluster_name: str,
                          conversation_history_mode: str = "final") -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:
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
    try:
        full_results = pd.read_csv(full_result_path)
    except BaseException:
        full_results = pd.read_csv(full_result_path, index_col=0)
    
    # Load marker data - handle both DataFrame and file path
    if isinstance(marker_path, pd.DataFrame):
        marker = marker_path
    else:
        try:
            marker = pd.read_csv(marker_path)
        except BaseException:
            marker = pd.read_csv(marker_path, index_col=0)
        
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
        raise ValueError(
            f"No data found for cluster '{cluster_name}' in the full results file. Available clusters: {full_results[cluster_column].unique().tolist()}")
    
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
                
                print(
                    f"Using conversation history for cluster {cluster_name} (mode: {conversation_history_mode}, {len(annotation_history)} characters)")
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


def generate_summary_report(conversation_history: List[Dict[str, str]], output_filename: str) -> str:
    """
    Generate a summarized report from the raw conversation history.
    
    Args:
        conversation_history: List of conversation messages
        output_filename: Path to save the summary report
        
    Returns:
        str: Path to the saved HTML report
    """
    try:
        # Extract content from conversation history, alternating between assistant and user
        full_conversation = ""
        for msg in conversation_history:
            role = msg.get('role', '')
            content = msg.get('content', '')
            if isinstance(content, list):
                content = str(content)
            full_conversation += f"\n## {role.upper()}\n{content}\n"

        # Craft the prompt for the LLM to generate a summary using simple tags
        prompt = f"""You are a specialized scientific report generator focusing on gene expression analysis.
        I will provide you a raw conversation history from a cell type annotation tool called CASSIA, which conducts iterative gene expression analysis.

        Your task is to generate a structured, concise summary report that highlights the key findings, hypotheses tested, and conclusions drawn.

        Use the following simple tag format exactly in your response:

        <OVERVIEW>
        Brief overview of what was analyzed and how many iterations were performed.
        Include the total number of genes examined across all iterations.
        </OVERVIEW>

        <INITIAL_ASSESSMENT>
        Summarize the initial cell type hypothesis and evidence from the first evaluation.
        </INITIAL_ASSESSMENT>

        <ITERATION_1>
        <HYPOTHESES>
        Clear explanation of what hypotheses were being tested in the first iteration and WHY.
        Format multiple hypotheses as numbered points (1., 2., 3.) each on a new line.
        </HYPOTHESES>

        <GENES_CHECKED>
        List the specific genes checked in this iteration (comma-separated).
        </GENES_CHECKED>

        <KEY_FINDINGS>
        Concise summary of the key results from gene expression analysis and what was learned.
        </KEY_FINDINGS>
        </ITERATION_1>

        # Repeat for each additional iteration (ITERATION_2, ITERATION_3, etc.)

        <FINAL_ANNOTATION>
        The conclusive cell type annotation, exactly as determined in the conversation.
        Include the confidence level and main supporting evidence.
        </FINAL_ANNOTATION>

        <MARKER_SUMMARY>
        List the most important marker genes that defined this cell type, organized by their functional groups.
        </MARKER_SUMMARY>

        <RECOMMENDATIONS>
        Any suggested next steps or validation approaches mentioned in the conversation.
        </RECOMMENDATIONS>

        Follow these guidelines:
        1. Maintain scientific precision while making the report accessible
        2. Include exact gene names with proper capitalization
        3. Keep your summary factual and based strictly on the conversation
        4. Use the exact tags as shown above
        5. Make sure to separate each iteration clearly
        6. Make the final cell type annotation stand out prominently

        Here's the conversation history to summarize:

        {full_conversation}
        """

        # Use the call_llm function to generate the summary
        summary = call_llm(
            prompt=prompt,
            provider="openrouter",  # Using OpenRouter as the provider
            model="google/gemini-2.5-flash-preview",  # Using Gemini 2.5 Flash model
            temperature=0.3,      # Low temperature for more consistent output
            max_tokens=4000
        )

        # Get the base name without extension
        base_filename = os.path.splitext(output_filename)[0]

        # Create filenames for raw and HTML versions
        raw_output_filename = base_filename + "_raw.txt"
        html_output_filename = base_filename + ".html"

        # Save the raw summary to a file
        with open(raw_output_filename, 'w', encoding='utf-8') as f:
            f.write(summary)
        print(f"Raw summary saved to {raw_output_filename}")

        # Convert to HTML and save
        html_path = format_summary_to_html(summary, html_output_filename)
        print(f"HTML summary saved to {html_path}")

        # Return the HTML file path
        return html_output_filename

    except Exception as e:
        error_msg = f"Error generating summary report: {str(e)}"
        print(error_msg)
        # Save error message to file so there's still an output
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(f"<error>{error_msg}</error>")
        return output_filename


def format_summary_to_html(summary_text: str, output_filename: str) -> str:
    """
    Convert the tagged summary into a properly formatted HTML report.
    
    Args:
        summary_text: Text with tags like <OVERVIEW>, <ITERATION_1>, etc.
        output_filename: Path to save the HTML report
        
    Returns:
        str: Path to the saved HTML report
    """
    try:
        import re
        
        # Helper function to format hypotheses with better separation
        def format_hypotheses(text: str) -> str:
            """
            Format hypotheses text to make numbered points stand out better.
            Detects numbered lists (1., 2., 3.) and formats them with proper styling.
            
            Args:
                text: The raw hypotheses text
                
            Returns:
                str: Formatted HTML for the hypotheses
            """
            # Handle case where there's no content
            if not text or text == "No information available":
                return text
                
            # Check if text contains numbered points (1. 2. 3. etc.)
            numbered_points = re.findall(r'(?:^|\s)(\d+\.)\s+([^0-9\.].*?)(?=\s+\d+\.\s+|\s*$)', text, re.DOTALL)
            
            # If we found numbered points, format them as a list
            if numbered_points:
                result = '<div class="hypothesis-list">'
                for num, content in numbered_points:
                    result += f'<div class="hypothesis-item"><span class="hypothesis-number">{num}</span> {content.strip()}</div>'
                result += '</div>'
                return result
                
            # Alternative approach: look for numbers at beginning of paragraphs
            paragraphs = text.split('\n')
            if any(re.match(r'^\s*\d+[\.\)\-]', p) for p in paragraphs if p.strip()):
                result = '<div class="hypothesis-list">'
                for para in paragraphs:
                    if not para.strip():
                        continue
                    if re.match(r'^\s*\d+[\.\)\-]', para):
                        num, content = re.match(r'^\s*(\d+[\.\)\-])\s*(.*)', para).groups()
                        result += f'<div class="hypothesis-item"><span class="hypothesis-number">{num}</span> {content.strip()}</div>'
                    else:
                        result += f'<div class="hypothesis-item-continued">{para.strip()}</div>'
                result += '</div>'
                return result
            
            # If the above patterns don't match, do a more aggressive search for numbered items
            numbered_items = re.split(r'(?:^|\s)(\d+[\.\)\-])(?=\s)', text)
            if len(numbered_items) > 2:  # We have at least one number + content
                result = '<div class="hypothesis-list">'
                # Skip the first empty item
                for i in range(1, len(numbered_items), 2):
                    if i+1 < len(numbered_items):
                        num = numbered_items[i]
                        content = numbered_items[i+1].strip()
                        result += f'<div class="hypothesis-item"><span class="hypothesis-number">{num}</span> {content}</div>'
                result += '</div>'
                return result
            
            # If no patterns match, just return the original text
            return text
            
        # Extract sections using regex
        sections = {}
        
        # Define the sections to extract
        section_patterns = {
            'overview': r'<OVERVIEW>\s*([\s\S]*?)\s*</OVERVIEW>',
            'initial_assessment': r'<INITIAL_ASSESSMENT>\s*([\s\S]*?)\s*</INITIAL_ASSESSMENT>',
            'final_annotation': r'<FINAL_ANNOTATION>\s*([\s\S]*?)\s*</FINAL_ANNOTATION>',
            'marker_summary': r'<MARKER_SUMMARY>\s*([\s\S]*?)\s*</MARKER_SUMMARY>',
            'recommendations': r'<RECOMMENDATIONS>\s*([\s\S]*?)\s*</RECOMMENDATIONS>',
        }
        
        # Extract each section
        for key, pattern in section_patterns.items():
            match = re.search(pattern, summary_text)
            if match:
                sections[key] = match.group(1).strip()
            else:
                sections[key] = "No information available"
        
        # Extract iterations
        iterations = []
        iter_pattern = r'<ITERATION_(\d+)>\s*([\s\S]*?)\s*</ITERATION_\1>'
        for match in re.finditer(iter_pattern, summary_text):
            iter_num = match.group(1)
            iter_content = match.group(2)
            
            # Extract subsections within each iteration
            hypotheses = re.search(r'<HYPOTHESES>\s*([\s\S]*?)\s*</HYPOTHESES>', iter_content)
            genes_checked = re.search(r'<GENES_CHECKED>\s*([\s\S]*?)\s*</GENES_CHECKED>', iter_content)
            key_findings = re.search(r'<KEY_FINDINGS>\s*([\s\S]*?)\s*</KEY_FINDINGS>', iter_content)
            
            iterations.append({
                'number': iter_num,
                'hypotheses': hypotheses.group(1).strip() if hypotheses else "No information available",
                'genes_checked': genes_checked.group(1).strip() if genes_checked else "No information available",
                'key_findings': key_findings.group(1).strip() if key_findings else "No information available"
            })
        
        # Sort iterations by number
        iterations.sort(key=lambda x: int(x['number']))
        
        # HTML template with CSS styling
        html = f"""
    <!DOCTYPE html>
    <html>
    <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>CASSIA Cell Type Annotation Summary</title>
        <style>
                    :root {{
                        --primary-color: #2563eb;
                        --secondary-color: #0891b2;
                        --accent-color: #4f46e5;
                        --light-bg: #f3f4f6;
                        --border-color: #e5e7eb;
                        --text-color: #1f2937;
                        --text-light: #6b7280;
                    }}
                    
            body {{ 
                        font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                line-height: 1.6;
                        color: var(--text-color);
                        background-color: #ffffff;
                        margin: 0;
                        padding: 0;
            }}
                    
            .container {{ 
                        max-width: 900px;
                        margin: 0 auto;
                        padding: 2rem;
                }}
                    
                    header {{
                        text-align: center;
                        margin-bottom: 2.5rem;
                        border-bottom: 2px solid var(--border-color);
                        padding-bottom: 1rem;
                    }}
                    
                    h1 {{
                        color: var(--primary-color);
                        font-size: 2.5rem;
                        margin-bottom: 0.5rem;
                    }}
                    
                    .subtitle {{
                        color: var(--text-light);
                        font-size: 1.2rem;
                        margin-top: 0;
                    }}
                    
                    section {{
                        margin-bottom: 2.5rem;
                        background-color: white;
                        border-radius: 0.5rem;
                        box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
                        padding: 1.5rem;
                        border-left: 4px solid var(--primary-color);
                    }}
                    
                    h2 {{
                        color: var(--primary-color);
                        font-size: 1.8rem;
                        margin-top: 0;
                        border-bottom: 1px solid var(--border-color);
                        padding-bottom: 0.5rem;
                    }}
                    
                    h3 {{
                        color: var(--secondary-color);
                        font-size: 1.3rem;
                        margin: 1.5rem 0 0.5rem;
                    }}
                    
                    ul, ol {{
                        padding-left: 1.5rem;
                    }}
                    
                    li {{
                        margin-bottom: 0.5rem;
                    }}
                    
                    .final-annotation {{
                        background-color: #ecfdf5;
                        border-left: 4px solid #10b981;
                    }}
                    
                    .final-annotation h2 {{
                        color: #10b981;
                    }}
                    
                    .iteration {{
                        margin-bottom: 1.5rem;
                        background-color: var(--light-bg);
                        border-radius: 0.5rem;
                        padding: 1.5rem;
                        border-left: 4px solid var(--secondary-color);
                    }}
                    
                    .iteration-title {{
                        font-size: 1.5rem;
                        color: var(--secondary-color);
                        margin-top: 0;
                        border-bottom: 1px solid var(--border-color);
                        padding-bottom: 0.5rem;
                    }}
                    
                    .gene-list {{
                        display: flex;
                        flex-wrap: wrap;
                        gap: 0.5rem;
                        margin: 1rem 0;
                    }}
                    
                    .gene-badge {{
                        background-color: var(--accent-color);
                        color: white;
                        padding: 0.3rem 0.7rem;
                        border-radius: 1rem;
                        font-size: 0.9rem;
                        font-weight: 500;
                    }}
                    
                    .sub-section {{
                        background-color: white;
                        border-radius: 0.3rem;
                        padding: 1rem;
                        margin-top: 1rem;
                    }}
                    
                    code {{
                        font-family: ui-monospace, monospace;
                        font-size: 0.9em;
                    }}
                    
                    .marker-category {{
                        margin-bottom: 1rem;
                    }}
                    
                    .marker-category-title {{
                        font-weight: 600;
                        margin-bottom: 0.5rem;
                        color: var(--secondary-color);
                    }}
                    
                    .hypothesis-list {{
                        display: flex;
                        flex-direction: column;
                        gap: 0.8rem;
                    }}
                    
                    .hypothesis-item {{
                        padding-left: 1.5rem;
                        position: relative;
                        margin-bottom: 0.5rem;
                    }}
                    
                    .hypothesis-number {{
                        position: absolute;
                        left: 0;
                        font-weight: 600;
                        color: var(--accent-color);
                    }}
                    
                    .hypothesis-item-continued {{
                        padding-left: 1.5rem;
                        margin-top: -0.3rem;
                        color: var(--text-light);
            }}
        </style>
    </head>
    <body>
        <div class="container">
                    <header>
                        <h1>CASSIA Cell Type Annotation Summary</h1>
                        <p class="subtitle">Single-cell RNA-seq Analysis Report</p>
                    </header>
                    
                    <section>
                        <h2>Overview</h2>
            <div class="content">
                            {sections['overview']}
            </div>
                    </section>
                    
                    <section>
                        <h2>Initial Assessment</h2>
                        <div class="content">
                            {sections['initial_assessment']}
                        </div>
                    </section>
            """
            
        # Add iterations
        for iteration in iterations:
            # Format genes checked as badges
            genes = iteration['genes_checked']
            gene_badges = ""
            if genes and genes != "No information available":
                gene_list = [g.strip() for g in re.split(r'[,\s]+', genes) if g.strip()]
                gene_badges = '<div class="gene-list">' + ''.join([f'<span class="gene-badge">{gene}</span>' for gene in gene_list]) + '</div>'
            
            html += f"""
                <section>
                    <h2>Iteration {iteration['number']}</h2>
                    
                    <div class="sub-section">
                        <h3>Hypotheses</h3>
                        <div class="content">
                            {format_hypotheses(iteration['hypotheses'])}
                        </div>
                    </div>
                    
                    <div class="sub-section">
                        <h3>Genes Checked</h3>
                        {gene_badges}
                    </div>
                    
                    <div class="sub-section">
                        <h3>Key Findings</h3>
                        <div class="content">
                            {iteration['key_findings']}
                        </div>
                    </div>
                </section>
            """
        
        # Add final annotation (highlighted)
        html += f"""
                <section class="final-annotation">
                    <h2>Final Annotation</h2>
                    <div class="content">
                        {sections['final_annotation']}
                    </div>
                </section>
        """
        
        # Add marker summary
        html += f"""
                <section>
                    <h2>Key Marker Genes</h2>
                    <div class="content">
                        {sections['marker_summary']}
                    </div>
                </section>
        """
        
        # Add recommendations if available
        if sections['recommendations'] and sections['recommendations'] != "No information available":
            html += f"""
                <section>
                    <h2>Recommendations</h2>
                    <div class="content">
                        {sections['recommendations']}
                    </div>
                </section>
            """
        
        # Close the HTML
        html += """
        </div>
    </body>
    </html>
    """
    
        # Save the HTML report
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"HTML report saved to {output_filename}")
        return output_filename
    
    except Exception as e:
        error_msg = f"Error formatting summary to HTML: {str(e)}"
        print(error_msg)
        # Save error message to file so there's still an output
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(f"<html><body><h1>Error</h1><p>{error_msg}</p></body></html>")
        return output_filename


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
        
        # Generate paths for reports
        if not output_name.lower().endswith('.html'):
            raw_report_path = output_name + '_raw_conversation.html'
        else:
            # Remove .html for base name
            output_name = output_name[:-5]
            raw_report_path = output_name + '_raw_conversation.html'

        # Generate path for summary report
        summary_report_path = output_name + '_summary.html'

        # Skip the first message which contains the prompt
        conversation_without_prompt = messages[1:] if len(messages) > 1 else messages
        
        # Generate the HTML reports
        try:
            # Generate the raw conversation report
            html_report = generate_html_report(analysis_text) if 'generate_html_report' in globals(
            ) else "<html><body><h1>Raw Output</h1><pre>" + analysis_text + "</pre></body></html>"
            save_html_report(html_report, raw_report_path)
            print(f"Raw conversation report saved to {raw_report_path}")

            # Generate the summary report
            try:
                summary_report_path = generate_summary_report(conversation_without_prompt, summary_report_path)
                print(f"Summary report saved to {summary_report_path}")
            except Exception as e:
                print(f"Warning: Could not generate summary report: {str(e)}")
                summary_report_path = None
        except Exception as e:
            print(f"Warning: Could not generate raw conversation report: {str(e)}")
            raw_report_path = None
            summary_report_path = None

        # Return a dictionary with paths to reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'raw_report_path': raw_report_path,
            'summary_report_path': summary_report_path,
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
            'raw_report_path': None,
            'summary_report_path': None,
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
        
        # Generate paths for reports
        if not output_name.lower().endswith('.html'):
            raw_report_path = output_name + '_raw_conversation.html'
        else:
            # Remove .html for base name
            output_name = output_name[:-5]
            raw_report_path = output_name + '_raw_conversation.html'

        # Generate path for summary report
        summary_report_path = output_name + '_summary.html'

        # Skip the first message which contains the prompt
        conversation_without_prompt = messages[1:] if len(messages) > 1 else messages
        
        # Generate the HTML reports
        try:
            # Generate the raw conversation report
            html_report = generate_html_report(analysis_text) if 'generate_html_report' in globals(
            ) else "<html><body><h1>Raw Output</h1><pre>" + analysis_text + "</pre></body></html>"
            save_html_report(html_report, raw_report_path)
            print(f"Raw conversation report saved to {raw_report_path}")

            # Generate the summary report
            try:
                summary_report_path = generate_summary_report(conversation_without_prompt, summary_report_path)
                print(f"Summary report saved to {summary_report_path}")
            except Exception as e:
                print(f"Warning: Could not generate summary report: {str(e)}")
                summary_report_path = None
        except Exception as e:
            print(f"Warning: Could not generate raw conversation report: {str(e)}")
            raw_report_path = None
            summary_report_path = None

        # Return a dictionary with paths to reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'raw_report_path': raw_report_path,
            'summary_report_path': summary_report_path,
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
            'raw_report_path': None,
            'summary_report_path': None,
            'execution_time': 0,
            'analysis_text': None
        }
