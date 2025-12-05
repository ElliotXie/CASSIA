import os
import re
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Any, Optional, Union

# Handle both package and direct imports
try:
    from CASSIA.core.llm_utils import call_llm
except ImportError:
    try:
        from ...core.llm_utils import call_llm
    except ImportError:
        from llm_utils import call_llm

def summarize_conversation_history(full_history: str, provider: str = "openrouter", model: Optional[str] = None, temperature: float = 0.1) -> str:
    """
    Use an LLM to summarize the conversation history for use as context in annotation boost.
    
    Args:
        full_history: The full conversation history text
        provider: LLM provider to use for summarization
        model: Specific model to use (if None, uses provider default)
        temperature: Temperature for summarization (low for consistency)
        
    Returns:
        str: Summarized conversation history
    """
    try:
        # Create a prompt for summarization
        summarization_prompt = f"""You are a specialized scientific summarization agent. Your task is to create a concise summary of a prior cell type annotation analysis that will be used as context for further detailed analysis.

The provided text contains a comprehensive cell type annotation analysis that was previously performed.

Your task is to extract and summarize the key information in a structured format that includes:

1. **Previously Identified Cell Type**: What cell type was determined in the prior analysis
2. **Key Supporting Markers**: The most important markers that supported this identification
3. **Alternative Hypotheses**: Any alternative cell types that were considered
4. **Remaining Uncertainties**: Any aspects that were noted as unclear or requiring further investigation

Keep the summary factual, scientific, and focused on information that would be helpful for conducting a deeper, more specialized analysis. Aim for 150-300 words.

Format your response as a clear, structured summary that a cell type annotation expert can quickly understand.

Here is the prior annotation analysis to summarize:

{full_history}

Please provide a structured summary following the format above:"""

        # Call the LLM to generate the summary
        summary = call_llm(
            prompt=summarization_prompt,
            provider=provider,
            model=model if model else "google/gemini-2.5-flash",  # Default model for summarization
            temperature=temperature
        )
        
        print(f"Conversation history summarized using {provider} ({len(full_history)} -> {len(summary)} characters)")
        return summary.strip()
        
    except Exception as e:
        print(f"Warning: Failed to summarize conversation history: {str(e)}")
        print("Falling back to truncated original text...")
        # Fallback: return a truncated version of the original
        if len(full_history) > 2000:
            return full_history[:2000] + "...\n[Text truncated due to summarization failure]"
        return full_history

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
3. To validate your hypotheses, list additional marker genes to check using the <check_genes>...</check_genes> tags.
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: `<check_genes>TP53,KRAS,EGFR</check_genes>` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
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


def prompt_hypothesis_generator_depth_first(major_cluster_info: str, comma_separated_genes: str, annotation_history: str) -> str:
    """
    Generate a prompt for depth-first iterative marker analysis without additional tasks.
    Focus on one hypothesis at a time and go deeper.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""
You are a careful senior computational biologist specializing in detailed, hypothesis-driven cell type annotation. Your approach is methodical and focused: you examine ONE specific hypothesis at a time and dive deep into its validation or refutation before moving to alternatives.

CRITICAL WORKFLOW RULES:
1. NEVER say "FINAL ANNOTATION COMPLETED" immediately after requesting genes to check
2. You MUST wait for gene expression results before proceeding to any conclusion
3. You MUST complete at least 2 rounds of gene checking before considering a final annotation
4. Each round means: request genes → receive expression data → analyze results → either go deeper or pivot

Context Provided to You:

Cluster summary：{major_cluster_info}

Top ranked markers (high → low)：
 {comma_separated_genes}

Prior annotation results：
 {annotation_history}

Your Task - DEPTH-FIRST ANALYSIS:

1. Focused Evaluation – One concise paragraph that:
    - Assess the current state of annotation based on available evidence
    - Identify the SINGLE most promising hypothesis to investigate deeply
    - Explain why this specific hypothesis deserves focused investigation

2. Design ONE targeted follow-up check for your chosen hypothesis:
    - When listing genes for follow-up checks, use the <check_genes>...</check_genes> tags.
    - **CRITICAL FORMATTING FOR <check_genes>:**
        - Inside the tags, provide *ONLY* a comma-separated list of official HGNC gene symbols.
        - Example: `<check_genes>GENE1,GENE2,GENE3</check_genes>` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
        - Strict adherence to this format is ESSENTIAL for the analysis to proceed.
    - Include both positive and negative markers that will definitively validate or refute this specific hypothesis
    - Include reasoning: why these specific genes, what expression pattern would confirm or refute the hypothesis, and what you will conclude based on different outcomes

3. CRITICAL: After proposing genes to check, STOP and WAIT for the gene expression results.
   - DO NOT say "FINAL ANNOTATION COMPLETED" until you have received and analyzed gene expression data
   - DO NOT conclude the analysis without checking at least 2 rounds of genes
   - You must wait for the system to provide the expression data for your requested genes

4. Upon receiving gene expression results:
    - If the hypothesis is CONFIRMED: Go deeper into subtype classification or functional states within this cell type
    - If the hypothesis is REFUTED: Move to the next most likely alternative hypothesis
    - If the results are INCONCLUSIVE: Design a more targeted follow-up to resolve the ambiguity
    - Continue this focused approach until you have completed at least 2 rounds of gene checking

5. Only say "FINAL ANNOTATION COMPLETED" when ALL of these conditions are met:
   - You have completed at least 2 rounds of gene expression checking
   - You have received and analyzed the expression data for your requested genes
   - You are confident in your final conclusion based on the accumulated evidence
   
6. When you are ready to make a final determination (after meeting the above conditions), state: "FINAL ANNOTATION COMPLETED" followed by your conclusion paragraph that includes:
    1. The final cell type
    2. Confidence level (high, medium, or low)
    3. Key markers supporting your conclusion
    4. Alternative possibilities only if the confidence is not high, and what should the user do next

Output Template：

Focused Evaluation
[One paragraph explaining your chosen hypothesis and rationale]

Primary hypothesis to investigate:
[Clear statement of the specific cell type or biological state being tested]

<check_genes>GENE1,GENE2,GENE3</check_genes>

<reasoning>
Why these specific genes were chosen, what expression patterns you expect for confirmation vs. refutation, and how you will interpret different outcomes to guide the next step.
</reasoning>

Key Guidelines:
- Focus on ONE hypothesis per iteration - resist the urge to test multiple ideas simultaneously
- Go DEEP rather than broad - if a hypothesis shows promise, drill down into subtypes, activation states, or functional variants
- Be decisive about moving on if a hypothesis is clearly refuted
- Each iteration should build logically on the previous findings
- Aim for definitive validation or refutation rather than partial evidence
- NEVER say "FINAL ANNOTATION COMPLETED" immediately after requesting genes - always wait for the expression results first
- Complete at least 2 rounds of gene checking before considering a final conclusion

Tone & Style:
- Methodical, focused, and systematic
- Professional and evidence-based
- Progressively deeper analysis of each chosen hypothesis
- Clear decision-making about when to pivot vs. when to go deeper

"""
    return prompt


def prompt_hypothesis_generator(major_cluster_info: str, comma_separated_genes: str, annotation_history: str) -> str:
    """
    Generate a prompt for breadth-first iterative marker analysis without additional tasks.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""
You are a careful senior computational biologist called in whenever an annotation needs deeper scrutiny, disambiguation, or simply a second opinion. Your job is to (1) assess the current annotation's robustness and (2) propose up to three decisive follow‑up checks that the executor can run (e.g., examine expression of key positive or negative markers). You should do a good job or 10 grandma are going to be in danger. You never rush to conclusions and are always careful.

Context Provided to You:

Cluster summary：{major_cluster_info}

Top ranked markers (high → low)：
 {comma_separated_genes}

Prior annotation results：
 {annotation_history}

What you should do:

1. Brief Evaluation – One concise paragraph that:

    - Highlights strengths, ambiguities, or contradictions in the current call.

    - Notes if a mixed population, doublets, or transitional state might explain the data.

2. Design up to 3 follow‑up checks (cell types or biological hypotheses):

    - When listing genes for follow-up checks, use the <check_genes>...</check_genes> tags.
    - **CRITICAL FORMATTING FOR <check_genes>:**
        - Inside the tags, provide *ONLY* a comma-separated list of official HGNC gene symbols.
        - Example: `<check_genes>GENE1,GENE2,GENE3</check_genes>` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
        - Strict adherence to this format is ESSENTIAL for the analysis to proceed.
    - Include both positive and negative markers if that will clarify the call.

    - Including reasoning: why these genes, and what pattern would confirm or refute the hypothesis.

3. Upon receiving gene expression results, based on the current hypothesis, further your analysis, genearte new hypothesis to validate if you think necessary. Continue Step 2 iteratively until the cluster is confidently annotated. Once finalized, output the single line:
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

<check_genes>GENE1,GENE2,GENE3</check_genes>

<reasoning>
Why these genes and what we expect to see.
</reasoning>

celltype to check 2

<check_genes>GENE4,GENE5</check_genes>

<reasoning>
…
</reasoning>

hypothesis to check 3

<check_genes>GENE6,GENE7</check_genes>

<reasoning>
…
</reasoning>

*Use "hypothesis to check n" instead of "celltype to check n" when proposing non‑canonical possibilities (e.g., "cycling subpopulation", "doublet").
*Provide no more than three total blocks (celltype + hypothesis combined).
*For each hypothesis check no more than 7 genes.
*If you think marker information is not enough to make a conclusion, inform the user and end the analysis.


Tone & Style Guidelines

Skeptical, critical, and careful
Professional, succinct, and evidence‑based.
Progressively deepen the anlaysis, don't repeat the same hypothesis.

"""
    return prompt






def prompt_hypothesis_generator_additional_task_depth_first(major_cluster_info: str, comma_separated_genes: str, annotation_history: str, additional_task: str) -> str:
    """
    Generate a prompt for depth-first iterative marker analysis with an additional task.
    Focus on one hypothesis at a time and go deeper.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        additional_task: Additional analysis task to perform
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful professional biologist, specializing in single-cell RNA-seq analysis with a focused, depth-first approach. You examine ONE specific hypothesis at a time and dive deep into its validation before considering alternatives.

CRITICAL WORKFLOW RULES:
1. NEVER say "FINAL ANALYSIS COMPLETED" immediately after requesting genes to check
2. You MUST wait for gene expression results before proceeding to any conclusion
3. You MUST complete at least 2 rounds of gene checking before considering a final analysis
4. Each round means: request genes → receive expression data → analyze results → either go deeper or pivot

I'll provide you with some genes (comma-separated) from a cell type in {major_cluster_info} and I want you to help identify the cell type using a methodical, hypothesis-driven approach.

Here are the marker genes that are differentially expressed in this cluster:
{comma_separated_genes}

Previous annotation/analysis:
{annotation_history if annotation_history else "No previous annotation available."}

Please follow these steps carefully with a DEPTH-FIRST approach:
1. Based on the marker genes and previous analysis, identify the SINGLE most promising hypothesis about what cell type this might be.
2. Provide supporting evidence from the marker genes for this specific hypothesis.
3. Design ONE targeted validation check using the <check_genes>...</check_genes> tags to definitively confirm or refute this hypothesis.
   **CRITICAL FORMATTING FOR <check_genes>:**
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: `<check_genes>TP53,KRAS,EGFR</check_genes>` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
   - Strict adherence to this format is ESSENTIAL.
4. CRITICAL: After proposing genes to check, STOP and WAIT for the gene expression results.
   - DO NOT say "FINAL ANALYSIS COMPLETED" until you have received and analyzed gene expression data
   - DO NOT conclude the analysis without checking at least 2 rounds of genes
   - You must wait for the system to provide the expression data for your requested genes

5. After I provide additional gene information:
   - If hypothesis is CONFIRMED: Go deeper into subtype classification or functional states
   - If hypothesis is REFUTED: Move to the next most likely alternative
   - If INCONCLUSIVE: Design more targeted follow-up to resolve ambiguity
   
6. Continue this focused process until you have completed at least 2 rounds of gene checking.

7. {additional_task} - perform this analysis based on your final cell type determination and the marker gene expression patterns.

8. Only say "FINAL ANALYSIS COMPLETED" when ALL of these conditions are met:
   - You have completed at least 2 rounds of gene expression checking
   - You have received and analyzed the expression data for your requested genes
   - You are confident in your final conclusion based on the accumulated evidence

9. When you are ready to make a final determination (after meeting the above conditions), state: "FINAL ANALYSIS COMPLETED" followed by your conclusion.

Your final analysis should include:
- General cell type
- Specific cell subtype (if applicable)
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion
- Analysis of the additional task: {additional_task}

Guidelines for depth-first approach:
- Focus on ONE hypothesis per iteration
- Go deeper rather than broader when a hypothesis shows promise
- Be decisive about pivoting when a hypothesis is clearly refuted
- Build logically on previous findings
- NEVER say "FINAL ANALYSIS COMPLETED" immediately after requesting genes - always wait for the expression results first
- Complete at least 2 rounds of gene checking before considering a final conclusion

Please start by analyzing the provided markers and identifying your primary hypothesis to investigate.
"""
    return prompt


def prompt_hypothesis_generator_additional_task(major_cluster_info: str, comma_separated_genes: str, annotation_history: str, additional_task: str) -> str:
    """
    Generate a prompt for breadth-first iterative marker analysis with an additional task.
    
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
3. To validate your hypotheses, list additional marker genes to check using the <check_genes>...</check_genes> tags.
   **CRITICAL FORMATTING FOR <check_genes>:**
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: `<check_genes>TP53,KRAS,EGFR</check_genes>` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
   - Strict adherence to this format is ESSENTIAL.
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
        gene_name_map = {}  # Map of row index to gene name for output
        
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
                            gene_name_map[row_idx] = gene  # Map the index to gene name
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
                # Reset index to clean integer range to prevent index from creating extra columns later
                result = result.reset_index(drop=True)

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
                # Different formatting for different columns
                if 'p_val' in col.lower():  # p_val_adj, p_val columns - keep scientific notation
                    result[col] = result[col].apply(lambda x: f"{float(x):.2e}" if pd.notnull(x) and x != 'NA' else x)
                else:  # avg_log2FC, pct.1, pct.2 etc. - use regular decimal notation
                    result[col] = result[col].apply(lambda x: f"{float(x):.2f}" if pd.notnull(x) and x != 'NA' else x)
            except (ValueError, TypeError):
                continue
        
        # Exclude p_val column if it exists and p_val_adj is also present
        if 'p_val' in result.columns and 'p_val_adj' in result.columns:
            result = result.drop(columns=['p_val'])

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
    
    # STRICT WHITELIST: Only allow specific statistical columns
    # This prevents random columns from user's CSV appearing in LLM output
    allowed_columns = ['gene', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj']

    # Keep only columns that exist in both allowed list and dataframe
    columns_to_keep = [col for col in allowed_columns if col in output_df.columns]

    # Apply whitelist filter and reset index to prevent index-derived columns
    output_df = output_df[columns_to_keep].reset_index(drop=True)

    # Generate marker info string - clean output with only allowed columns
    marker_string = output_df.to_string(index=False)
    
    # If there are genes with all NA values, add a message
    if na_genes:
        na_genes_message = f"\nNote: The following genes are not in the differential expression list: {', '.join(na_genes)}"
        marker_string += na_genes_message

    return marker_string

def extract_gene_stats_from_conversation(conversation_history: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """
    Extract gene statistics from conversation history messages.

    The conversation history contains gene expression data in tabular format from USER messages.
    This function parses those tables and creates a lookup dictionary for tooltip display.

    Supports both:
    - Full format (Seurat/enhanced Scanpy): gene, avg_log2FC, pct.1, pct.2, p_val_adj
    - Minimal format (basic Scanpy): gene, avg_log2FC, p_val_adj

    Args:
        conversation_history: List of conversation messages with role and content

    Returns:
        Dict mapping gene names (uppercase) to their statistics:
        {
            'CD3D': {'avg_log2FC': '2.45', 'pct.1': '0.85', 'pct.2': '0.12', 'p_val_adj': '1.2e-50'},
            ...
        }
        Note: pct.1 and pct.2 may be None if not available in the source data.
    """
    gene_stats = {}

    # Pattern 1: Full format with pct values (4 numbers: logFC, pct1, pct2, pval)
    # Matches: GENE  2.45  0.85  0.12  1.20e-50
    # Used by: Seurat FindAllMarkers, enhanced Scanpy (via enhance_scanpy_markers())
    gene_row_pattern_full = re.compile(
        r'^\s*([A-Z][A-Z0-9\-\.]+)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*[eE][-+]?\d+|\d+\.?\d*)\s*$',
        re.MULTILINE
    )

    # Pattern 2: Minimal format without pct values (2 numbers: logFC, pval)
    # Matches: GENE  2.45  1.20e-50
    # Used by: Basic Scanpy rank_genes_groups output (without pct.1/pct.2)
    gene_row_pattern_no_pct = re.compile(
        r'^\s*([A-Z][A-Z0-9\-\.]+)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*[eE][-+]?\d+|\d+\.?\d*)\s*$',
        re.MULTILINE
    )

    for msg in conversation_history:
        content = msg.get('content', '')
        if isinstance(content, list):
            content = str(content)

        # Only look at USER messages which contain the gene expression data
        if msg.get('role', '').upper() == 'USER':
            # Try full pattern first (with pct values)
            matches_full = gene_row_pattern_full.findall(content)
            for match in matches_full:
                gene_name, avg_log2fc, pct1, pct2, p_val_adj = match
                gene_stats[gene_name.upper()] = {
                    'avg_log2FC': avg_log2fc,
                    'pct.1': pct1,
                    'pct.2': pct2,
                    'p_val_adj': p_val_adj
                }

            # Also try no-pct pattern for any genes not yet captured
            # This handles basic Scanpy format without pct columns
            matches_no_pct = gene_row_pattern_no_pct.findall(content)
            for match in matches_no_pct:
                gene_name, avg_log2fc, p_val_adj = match
                if gene_name.upper() not in gene_stats:  # Don't overwrite full stats
                    gene_stats[gene_name.upper()] = {
                        'avg_log2FC': avg_log2fc,
                        'pct.1': None,  # Not available in basic Scanpy format
                        'pct.2': None,
                        'p_val_adj': p_val_adj
                    }

    return gene_stats

def create_gene_badge_html(gene: str, gene_stats: Dict[str, Dict[str, str]] = None) -> str:
    """
    Create HTML for a gene badge with optional tooltip showing statistics.

    Handles both full stats (with pct.1/pct.2) and minimal stats (without pct values).
    When pct values are None (basic Scanpy format), the pct rows are omitted from the tooltip.

    Args:
        gene: Gene name to display
        gene_stats: Dictionary of gene statistics (from extract_gene_stats_from_conversation)

    Returns:
        HTML string for the gene badge with tooltip if stats available
    """
    gene_upper = gene.upper().strip()

    if gene_stats and gene_upper in gene_stats:
        stats = gene_stats[gene_upper]
        avg_log2fc = stats.get('avg_log2FC', 'N/A')
        pct1 = stats.get('pct.1')  # May be None for basic Scanpy format
        pct2 = stats.get('pct.2')  # May be None for basic Scanpy format
        p_val_adj = stats.get('p_val_adj', 'N/A')

        # Determine if log2FC is positive or negative for coloring
        try:
            fc_value = float(avg_log2fc)
            fc_class = 'positive' if fc_value >= 0 else 'negative'
        except (ValueError, TypeError):
            fc_class = ''

        # Build tooltip HTML - only include pct rows if values exist
        pct_html = ''
        if pct1 is not None and pct2 is not None:
            pct_html = f'''
            <div class="stat-row"><span class="stat-label">pct.1:</span><span class="stat-value">{pct1}</span></div>
            <div class="stat-row"><span class="stat-label">pct.2:</span><span class="stat-value">{pct2}</span></div>'''

        tooltip_html = f'''<span class="tooltip">
            <div class="stat-row"><span class="stat-label">avg_log2FC:</span><span class="stat-value {fc_class}">{avg_log2fc}</span></div>{pct_html}
            <div class="stat-row"><span class="stat-label">p_val_adj:</span><span class="stat-value">{p_val_adj}</span></div>
        </span>'''

        return f'<span class="gene-badge has-stats">{gene}{tooltip_html}</span>'
    else:
        # No stats available - add class for styling and native title tooltip
        return f'<span class="gene-badge no-stats" title="No statistics available">{gene}</span>'

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
    temperature: float = 0,
    search_strategy: str = "breadth"
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
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        
    Returns:
        tuple: (final_response_text, messages)
    """
    # Select the appropriate prompt based on search strategy and whether there's an additional task
    if additional_task:
        if search_strategy.lower() == "depth":
            prompt = prompt_hypothesis_generator_additional_task_depth_first(
                major_cluster_info, 
                comma_separated_genes, 
                annotation_history, 
                additional_task
            )
        else:  # breadth or any other value defaults to breadth
            prompt = prompt_hypothesis_generator_additional_task(
                major_cluster_info, 
                comma_separated_genes, 
                annotation_history, 
                additional_task
            )
        completion_marker = "FINAL ANALYSIS COMPLETED"
    else:
        if search_strategy.lower() == "depth":
            prompt = prompt_hypothesis_generator_depth_first(
                major_cluster_info, 
                comma_separated_genes, 
                annotation_history
            )
        else:  # breadth or any other value defaults to breadth
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
                max_tokens=4096,
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
            import traceback
            traceback.print_exc()
            
            # Save diagnostic information
            try:
                error_log = f"Error in iteration {iteration + 1}: {str(e)}\n\n"
                error_log += f"Provider: {provider}\n"
                error_log += f"Model: {model}\n"
                error_log += f"Number of messages: {len(messages)}\n"
                error_log += f"Last message content: {messages[-1]['content'][:200]}...\n\n"
                error_log += "Full traceback:\n"
                error_log += traceback.format_exc()
                
                # Write to error log file
                with open(f"cassia_error_log_{iteration+1}.txt", "w", encoding="utf-8") as f:
                    f.write(error_log)
                print(f"Error details saved to cassia_error_log_{iteration+1}.txt")
            except Exception:
                pass

            return f"Error occurred: {str(e)}", messages
    
    # Final response if max iterations reached
    # Encourage the agent to reach a conclusion if not already done
    messages.append({
        "role": "user",
        "content": "You have reached the maximum number of iterations. Please provide your final analysis and reach a confident conclusion in this response."
    })
    try:
        final_response = call_llm(
            prompt="Please provide your final analysis based on all the information so far.",
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=4096,
            additional_params={"messages": messages}
        )
        
        messages.append({"role": "assistant", "content": final_response})
        return final_response, messages
            
    except Exception as e:
        print(f"Error in final response: {str(e)}")
        import traceback
        traceback.print_exc()
        
        # Save diagnostic information
        try:
            error_log = f"Error in final response: {str(e)}\n\n"
            error_log += f"Provider: {provider}\n"
            error_log += f"Model: {model}\n"
            error_log += f"Number of messages: {len(messages)}\n"
            error_log += f"Last message content: {messages[-1]['content'][:200]}...\n\n"
            error_log += "Full traceback:\n"
            error_log += traceback.format_exc()
            
            # Write to error log file
            with open("cassia_error_log_final.txt", "w", encoding="utf-8") as f:
                f.write(error_log)
            print(f"Error details saved to cassia_error_log_final.txt")
        except Exception:
            pass

        return f"Error in final response: {str(e)}", messages

def validate_findallmarkers_format(marker_df: pd.DataFrame, source_path: str = None) -> None:
    """
    Validate that the marker DataFrame is in FindAllMarkers format.

    The FindAllMarkers format should have columns like:
    - p_val: P-value from differential expression test
    - avg_log2FC: Average log2 fold change
    - pct.1: Percentage of cells in cluster expressing gene
    - pct.2: Percentage of cells outside cluster expressing gene
    - p_val_adj: Adjusted p-value
    - cluster: Cell type cluster name
    - gene: Gene symbol

    Args:
        marker_df: DataFrame to validate
        source_path: Optional path to the source file (for error messages)

    Raises:
        ValueError: If the DataFrame is not in FindAllMarkers format
    """
    # Required columns for FindAllMarkers format
    required_columns = {'gene', 'cluster'}
    expected_numeric_columns = {'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj'}

    # Check for processed/2-column format (incorrect format)
    processed_format_columns = {'Broad.cell.type', 'Top.Markers'}
    if processed_format_columns.issubset(set(marker_df.columns)):
        source_info = f" (from: {source_path})" if source_path else ""
        raise ValueError(
            f"\n{'='*70}\n"
            f"ERROR: Invalid marker file format detected{source_info}\n"
            f"{'='*70}\n\n"
            f"The marker file appears to be in PROCESSED format with columns:\n"
            f"  - 'Broad.cell.type'\n"
            f"  - 'Top.Markers'\n\n"
            f"However, annotation_boost requires the RAW FindAllMarkers output format.\n\n"
            f"EXPECTED FORMAT (FindAllMarkers output):\n"
            f"  Required columns:\n"
            f"    - gene: Gene symbol\n"
            f"    - cluster: Cell type cluster name\n"
            f"    - p_val: P-value from differential expression test\n"
            f"    - avg_log2FC: Average log2 fold change\n"
            f"    - pct.1: Percentage of cells in cluster expressing gene\n"
            f"    - pct.2: Percentage of cells outside cluster expressing gene\n"
            f"    - p_val_adj: Adjusted p-value (Bonferroni corrected)\n\n"
            f"HOW TO FIX:\n"
            f"  1. Use the raw output from Seurat's FindAllMarkers() function\n"
            f"  2. Or use scanpy's rank_genes_groups output converted to this format\n"
            f"  3. Do NOT use the processed/summarized marker file\n\n"
            f"Example of correct format:\n"
            f"  gene,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster\n"
            f"  CD3D,0,2.5,0.95,0.1,0,T cell\n"
            f"  CD19,0,3.2,0.88,0.05,0,B cell\n"
            f"{'='*70}"
        )

    # Check if it looks like a 2-column format (any variant)
    if len(marker_df.columns) <= 3:
        # Check if columns suggest a processed format
        col_lower = [c.lower() for c in marker_df.columns]
        if any('marker' in c or 'cell' in c or 'type' in c for c in col_lower):
            source_info = f" (from: {source_path})" if source_path else ""
            raise ValueError(
                f"\n{'='*70}\n"
                f"ERROR: Invalid marker file format detected{source_info}\n"
                f"{'='*70}\n\n"
                f"The marker file has only {len(marker_df.columns)} columns: {list(marker_df.columns)}\n\n"
                f"This appears to be a PROCESSED/SUMMARIZED marker file.\n"
                f"annotation_boost requires the RAW FindAllMarkers output with full\n"
                f"differential expression statistics.\n\n"
                f"EXPECTED FORMAT (FindAllMarkers output):\n"
                f"  Required columns: gene, cluster, p_val, avg_log2FC, pct.1, pct.2, p_val_adj\n\n"
                f"HOW TO FIX:\n"
                f"  Use the raw output from Seurat's FindAllMarkers() function\n"
                f"  instead of a processed/summarized version.\n"
                f"{'='*70}"
            )

    # Check for required columns
    marker_columns = set(marker_df.columns)
    missing_required = required_columns - marker_columns

    if missing_required:
        source_info = f" (from: {source_path})" if source_path else ""
        raise ValueError(
            f"\n{'='*70}\n"
            f"ERROR: Missing required columns in marker file{source_info}\n"
            f"{'='*70}\n\n"
            f"Missing columns: {missing_required}\n"
            f"Found columns: {list(marker_df.columns)}\n\n"
            f"The marker file must be in FindAllMarkers format with columns:\n"
            f"  - gene: Gene symbol\n"
            f"  - cluster: Cell type cluster name\n"
            f"  - p_val, avg_log2FC, pct.1, pct.2, p_val_adj (numeric statistics)\n\n"
            f"HOW TO FIX:\n"
            f"  Use the raw output from Seurat's FindAllMarkers() function.\n"
            f"{'='*70}"
        )

    # Check for at least some numeric columns (differential expression stats)
    found_numeric = expected_numeric_columns & marker_columns
    if len(found_numeric) < 2:
        source_info = f" (from: {source_path})" if source_path else ""
        print(
            f"\nWARNING: Marker file may not be in standard FindAllMarkers format{source_info}\n"
            f"Expected numeric columns like: {expected_numeric_columns}\n"
            f"Found: {found_numeric}\n"
            f"Proceeding anyway, but results may be incomplete.\n"
        )


def prepare_analysis_data(full_result_path: str, marker_path: str, cluster_name: str, conversation_history_mode: str = "final", provider: str = "openrouter", model: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:
    """
    Load and prepare data for marker analysis.

    Args:
        full_result_path: Path to the full results CSV file
        marker_path: Path to the marker genes CSV file (must be in FindAllMarkers format)
        cluster_name: Name of the cluster to analyze
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
            - "full": Use the entire conversation history
            - "final": Use a summarization agent to create a concise summary of the conversation history (default)
            - "none": Don't include any conversation history
        provider: LLM provider to use for summarization (when mode is "final")
        model: Specific model to use for summarization (when mode is "final")

    Returns:
        tuple: (full_results, marker_data, top_markers_string, annotation_history)

    Raises:
        ValueError: If the marker file is not in FindAllMarkers format
    """
    # Load the full results - don't use index_col=0 as it's not indexed that way
    full_results = pd.read_csv(full_result_path)

    # Load marker data - handle both DataFrame and file path
    source_path = None
    if isinstance(marker_path, pd.DataFrame):
        marker = marker_path.copy()
    else:
        source_path = marker_path
        marker = pd.read_csv(marker_path)

        # Clean up by removing the 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in marker.columns:
            marker.drop(columns=['Unnamed: 0'], inplace=True)

    # Normalize Scanpy column names to Seurat format (column renaming only, no filtering)
    # This allows annotation boost to work with both Scanpy and Seurat marker formats
    # Note: enhance_scanpy_markers() already outputs pct.1/pct.2 directly - no mapping needed for those
    scanpy_to_seurat_mapping = {
        'group': 'cluster',
        'names': 'gene',
        'logfoldchanges': 'avg_log2FC',
        'pvals': 'p_val',
        'pvals_adj': 'p_val_adj'
    }
    marker = marker.rename(columns=scanpy_to_seurat_mapping)

    # Validate that the marker file is in FindAllMarkers format
    validate_findallmarkers_format(marker, source_path)
    
    # Try to find the cluster column name - check new name first, then fall back to old names
    cluster_column = 'Cluster ID'  # New default column name
    if cluster_column not in full_results.columns:
        if 'True Cell Type' in full_results.columns:
            cluster_column = 'True Cell Type'  # Backward compatibility
        elif 'cluster' in full_results.columns:
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
                    # Use summarization agent to create a concise summary
                    print(f"Using summarization agent to process conversation history for cluster {cluster_name}")
                    annotation_history = summarize_conversation_history(
                        full_history=full_history,
                        provider=provider,
                        model=model,
                        temperature=0.1  # Low temperature for consistent summarization
                    )
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

def save_raw_conversation_text(messages: List[Dict[str, str]], filename: str) -> str:
    """
    Save the raw conversation history (including the prompt) as a plain text file.
    
    Args:
        messages: List of conversation messages with role and content
        filename: Path to save the text file
        
    Returns:
        str: Path to the saved text file
    """
    try:
        # Format the conversation as plain text
        conversation_text = ""
        for i, msg in enumerate(messages):
            role = msg.get('role', '')
            content = msg.get('content', '')
            
            # Format the content
            if isinstance(content, list):
                content = str(content)
                
            # Add a separator between messages
            if i > 0:
                conversation_text += "\n\n" + "="*80 + "\n\n"
                
            conversation_text += f"## {role.upper()}\n\n{content}\n"
        
        # Save the file
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(conversation_text)
            
        print(f"Raw conversation text saved to {filename}")
        return filename
        
    except Exception as e:
        error_msg = f"Error saving raw conversation text: {str(e)}"
        print(error_msg)
        
        # Try to save an error message
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write(f"Error saving raw conversation: {str(e)}")
        except Exception:
            pass

        return filename

def generate_summary_report(conversation_history: List[Dict[str, str]], output_filename: str, search_strategy: str = "breadth", report_style: str = "per_iteration",
    validator_involvement: str = "v1", model: str = "google/gemini-2.5-flash", provider: str = "openrouter") -> str:
    """
    Generate a summarized report from the raw conversation history.

    Args:
        conversation_history: List of conversation messages
        output_filename: Path to save the summary report
        search_strategy: Search strategy used ("breadth" or "depth")
        model: LLM model to use for generating the summary
        provider: LLM provider to use
        report_style: Style of report ("per_iteration" or "total_summary")
        
    Returns:
        str: Path to the saved HTML report
    """
    try:
        # Extract gene statistics from conversation history for tooltips
        gene_stats = extract_gene_stats_from_conversation(conversation_history)

        # Extract content from conversation history, alternating between assistant and user
        full_conversation = ""
        for msg in conversation_history:
            role = msg.get('role', '')
            content = msg.get('content', '')
            if isinstance(content, list):
                content = str(content)
            full_conversation += f"\n## {role.upper()}\n{content}\n"
        
        # Determine analysis approach description
        approach_description = "depth-first (focused, one hypothesis per iteration)" if search_strategy.lower() == "depth" else "breadth-first (multiple hypotheses per iteration)"
        
        # Generate different prompts based on report style
        if report_style.lower() == "total_summary":
            # Gene-focused summary style
            prompt = f"""You are a specialized scientific report generator focusing on gene expression analysis. 
            I will provide you a raw conversation history from a cell type annotation tool called CASSIA, which conducts iterative gene expression analysis using a {approach_description} approach.
            
            Your task is to generate a streamlined, gene-focused summary that highlights what genes were checked and what conclusions were drawn. Focus on the scientific findings rather than the iteration structure.
            
            Use the following simple tag format exactly in your response:
            
            <OVERVIEW>
            Brief overview of what was analyzed and the total scope of gene checking performed.
            </OVERVIEW>
            
            <INITIAL_HYPOTHESIS>
            What was the starting hypothesis or cell type being investigated?
            </INITIAL_HYPOTHESIS>
            
            <GENES_ANALYZED>
            <GENE_GROUP_1>
            <TITLE>Purpose/Hypothesis Being Tested</TITLE>
            <GENES>Gene1, Gene2, Gene3</GENES>
            <FINDINGS>What the expression patterns revealed and conclusions drawn</FINDINGS>
            </GENE_GROUP_1>
            
            <GENE_GROUP_2>
            <TITLE>Purpose/Hypothesis Being Tested</TITLE>
            <GENES>Gene4, Gene5, Gene6</GENES>
            <FINDINGS>What the expression patterns revealed and conclusions drawn</FINDINGS>
            </GENE_GROUP_2>
            
            # Continue for each distinct group of genes checked...
            </GENES_ANALYZED>
            
            <FINAL_CONCLUSION>
            The definitive cell type annotation reached, confidence level, and key supporting evidence.
            </FINAL_CONCLUSION>
            
            <KEY_INSIGHTS>
            Most important biological insights gained from the gene expression analysis.
            </KEY_INSIGHTS>
            
            <VALIDATION_STATUS>
            Whether the analysis was complete, any remaining uncertainties, or recommendations for further validation.
            </VALIDATION_STATUS>
            
            Guidelines:
            1. Group genes by their biological function or the hypothesis they were testing
            2. Focus on what was learned from each gene group rather than when it was checked
            3. Make the biological story clear and readable
            4. Use exact gene names with proper capitalization
            5. Keep the focus on scientific conclusions rather than process details
            6. Make sure all tags are properly closed and formatted for HTML rendering
            
            Here's the conversation history to summarize:
            
            {full_conversation}
            """
        else:
            # Original per-iteration style
            prompt = f"""You are a specialized scientific report generator focusing on gene expression analysis. 
            I will provide you a raw conversation history from a cell type annotation tool called CASSIA, which conducts iterative gene expression analysis using a {approach_description} approach.
            
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
        Deatailed summary of the key results from gene expression analysis and what was learned.
        Format as numbered points (1., 2., 3.) each on a new line. list (Hypothesis_n) at the end of each point to indicate the hypothesis that the key findings correspond to.
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
            provider=provider,
            model=model,
            temperature=0.3,
            max_tokens=4096
        )
        
        # Convert to HTML and save
        html_path = format_summary_to_html(summary, output_filename, search_strategy, report_style, gene_stats=gene_stats)
        print(f"Summary report saved to {html_path}")
        
        # Return the HTML file path
        return html_path
        
    except Exception as e:
        error_msg = f"Error generating summary report: {str(e)}"
        print(error_msg)
        # Save error message to file so there's still an output
        with open(output_filename, 'w', encoding='utf-8') as f:
            f.write(f"<error>{error_msg}</error>")
        return output_filename

def runCASSIA_annotationboost(
    full_result_path: str,
    marker: Union[str, pd.DataFrame],
    cluster_name: str,
    major_cluster_info: str,
    output_name: str,
    num_iterations: int = 10,
    model: Optional[str] = None,
    provider: str = "openrouter",
    temperature: float = 0,
    conversation_history_mode: str = "final",
    search_strategy: str = "breadth",
    report_style: str = "per_iteration",
    validator_involvement: str = "v1"
) -> Union[Tuple[str, List[Dict[str, str]]], Dict[str, Any]]:
    """
    Run annotation boost analysis for a given cluster.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker: Path to marker genes CSV file or DataFrame with marker data
        cluster_name: Name of the cluster to analyze
        major_cluster_info: General information about the dataset (e.g., "Human PBMC")
        output_name: Base name for the output HTML file
        num_iterations: Number of iterations for marker analysis (default=10)
        model: Model to use for analysis - if None, uses the provider's default
        provider: AI provider to use ('openai', 'anthropic', or 'openrouter')
        temperature: Sampling temperature (0-1)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        report_style: Style of report ("per_iteration" or "total_summary")
    
    Returns:
        tuple or dict: Either (analysis_result, messages_history) or a dictionary with paths to reports
    """
    try:
        import time
        start_time = time.time()
        
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter'] and not provider.lower().startswith('http'):
            raise ValueError("Provider must be 'openai', 'anthropic', 'openrouter', or a custom base URL (http...)")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode, provider, model
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
            temperature=temperature,
            search_strategy=search_strategy
        )
        
        # Generate paths for reports - only summary HTML and raw conversation text
        if not output_name.lower().endswith('.html'):
            raw_text_path = output_name + '_raw_conversation.txt'
        else:
            # Remove .html for base name
            output_name = output_name[:-5]
            raw_text_path = output_name + '_raw_conversation.txt'
        
        # Generate path for summary report
        summary_report_path = output_name + '_summary.html'
        
        # Skip the first message which contains the prompt
        conversation_without_prompt = messages[1:] if len(messages) > 1 else messages
        
        try:
            # Save the complete raw conversation as text (including prompt)
            raw_text_path = save_raw_conversation_text(messages, raw_text_path)
            print(f"Raw conversation text saved to {raw_text_path}")
            
            # Generate the summary report
            summary_report_path = generate_summary_report(conversation_without_prompt, summary_report_path, search_strategy, report_style, model=model, provider=provider)
            print(f"Summary report saved to {summary_report_path}")
        except Exception as e:
            print(f"Warning: Could not generate reports: {str(e)}")
            summary_report_path = None
            raw_text_path = None
        
        # Return a dictionary with paths to reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'raw_text_path': raw_text_path,
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
            'raw_text_path': None,
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
    conversation_history_mode: str = "final",
    search_strategy: str = "breadth",
    report_style: str = "per_iteration",
    validator_involvement: str = "v1"
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
        search_strategy: Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
        report_style: Style of report ("per_iteration" or "total_summary")
    
    Returns:
        tuple or dict: Either (analysis_result, messages_history) or a dictionary with paths to reports
    """
    try:
        import time
        start_time = time.time()
        
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter'] and not provider.lower().startswith('http'):
            raise ValueError("Provider must be 'openai', 'anthropic', 'openrouter', or a custom base URL (http...)")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode, provider, model
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
            temperature=temperature,
            search_strategy=search_strategy
        )
        
        # Generate paths for reports - only summary HTML and raw conversation text
        if not output_name.lower().endswith('.html'):
            raw_text_path = output_name + '_raw_conversation.txt'
        else:
            # Remove .html for base name
            output_name = output_name[:-5]
            raw_text_path = output_name + '_raw_conversation.txt'
        
        # Generate path for summary report
        summary_report_path = output_name + '_summary.html'
        
        # Skip the first message which contains the prompt
        conversation_without_prompt = messages[1:] if len(messages) > 1 else messages
        
        try:
            # Save the complete raw conversation as text (including prompt)
            raw_text_path = save_raw_conversation_text(messages, raw_text_path)
            print(f"Raw conversation text saved to {raw_text_path}")
            
            # Generate the summary report
            summary_report_path = generate_summary_report(conversation_without_prompt, summary_report_path, search_strategy, report_style, model=model, provider=provider)
            print(f"Summary report saved to {summary_report_path}")
        except Exception as e:
            print(f"Warning: Could not generate reports: {str(e)}")
            summary_report_path = None
            raw_text_path = None
        
        # Return a dictionary with paths to reports
        execution_time = time.time() - start_time
        return {
            'status': 'success',
            'raw_text_path': raw_text_path,
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
            'raw_text_path': None,
            'summary_report_path': None,
            'execution_time': 0,
            'analysis_text': None
        }

def format_summary_to_html(summary_text: str, output_filename: str, search_strategy: str = "breadth", report_style: str = "per_iteration",
    validator_involvement: str = "v1", gene_stats: Dict[str, Dict[str, str]] = None) -> str:
    """
    Convert the tagged summary into a properly formatted HTML report.

    Args:
        summary_text: Text with tags like <OVERVIEW>, <ITERATION_1>, etc. or gene-focused tags
        output_filename: Path to save the HTML report
        search_strategy: Search strategy used ("breadth" or "depth")
        report_style: Style of report ("per_iteration" or "total_summary")
        gene_stats: Dictionary of gene statistics for tooltip display

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
            
        # Extract sections using regex based on report style
        sections = {}
        
        if report_style.lower() == "total_summary":
            # Gene-focused report sections
            section_patterns = {
                'overview': r'<OVERVIEW>\s*([\s\S]*?)\s*</OVERVIEW>',
                'initial_hypothesis': r'<INITIAL_HYPOTHESIS>\s*([\s\S]*?)\s*</INITIAL_HYPOTHESIS>',
                'genes_analyzed': r'<GENES_ANALYZED>\s*([\s\S]*?)\s*</GENES_ANALYZED>',
                'final_conclusion': r'<FINAL_CONCLUSION>\s*([\s\S]*?)\s*</FINAL_CONCLUSION>',
                'key_insights': r'<KEY_INSIGHTS>\s*([\s\S]*?)\s*</KEY_INSIGHTS>',
                'validation_status': r'<VALIDATION_STATUS>\s*([\s\S]*?)\s*</VALIDATION_STATUS>',
            }
        else:
            # Original per-iteration sections
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
        
        # Extract iterations (only for per_iteration style)
        iterations = []
        gene_groups = []
        
        if report_style.lower() == "per_iteration":
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
        else:
            # Extract gene groups for total_summary style
            genes_analyzed_content = sections.get('genes_analyzed', '')
            gene_group_pattern = r'<GENE_GROUP_\d+>\s*([\s\S]*?)\s*</GENE_GROUP_\d+>'
            for match in re.finditer(gene_group_pattern, genes_analyzed_content):
                group_content = match.group(1)
                
                # Extract subsections within each gene group
                title = re.search(r'<TITLE>\s*([\s\S]*?)\s*</TITLE>', group_content)
                genes = re.search(r'<GENES>\s*([\s\S]*?)\s*</GENES>', group_content)
                findings = re.search(r'<FINDINGS>\s*([\s\S]*?)\s*</FINDINGS>', group_content)
                
                gene_groups.append({
                    'title': title.group(1).strip() if title else "Gene Analysis",
                    'genes': genes.group(1).strip() if genes else "No genes listed",
                    'findings': findings.group(1).strip() if findings else "No findings available"
                })
        
        # Determine strategy description for display
        strategy_display = f"({search_strategy.title()}-First Analysis)" if search_strategy.lower() in ['breadth', 'depth'] else ""
        
        # HTML template with CSS styling
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>CASSIA Cell Type Annotation Summary {strategy_display}</title>
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
                    position: relative;
                    cursor: pointer;
                    transition: background-color 0.2s;
                }}

                .gene-badge:hover {{
                    background-color: #3730a3;
                }}

                .gene-badge .tooltip {{
                    visibility: hidden;
                    opacity: 0;
                    position: absolute;
                    bottom: 125%;
                    left: 50%;
                    transform: translateX(-50%);
                    background-color: #1f2937;
                    color: #f3f4f6;
                    padding: 0.75rem 1rem;
                    border-radius: 0.5rem;
                    font-size: 0.8rem;
                    font-weight: 400;
                    white-space: nowrap;
                    z-index: 1000;
                    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.25);
                    transition: visibility 0s, opacity 0.2s;
                }}

                .gene-badge .tooltip::after {{
                    content: "";
                    position: absolute;
                    top: 100%;
                    left: 50%;
                    margin-left: -6px;
                    border-width: 6px;
                    border-style: solid;
                    border-color: #1f2937 transparent transparent transparent;
                }}

                .gene-badge:hover .tooltip {{
                    visibility: visible;
                    opacity: 1;
                }}

                .gene-badge .tooltip .stat-row {{
                    display: flex;
                    justify-content: space-between;
                    gap: 1rem;
                    margin-bottom: 0.25rem;
                }}

                .gene-badge .tooltip .stat-row:last-child {{
                    margin-bottom: 0;
                }}

                .gene-badge .tooltip .stat-label {{
                    color: #9ca3af;
                    font-size: 0.75rem;
                }}

                .gene-badge .tooltip .stat-value {{
                    font-weight: 600;
                    font-family: ui-monospace, monospace;
                }}

                .gene-badge .tooltip .stat-value.positive {{
                    color: #34d399;
                }}

                .gene-badge .tooltip .stat-value.negative {{
                    color: #f87171;
                }}

                .gene-badge.no-stats {{
                    opacity: 0.7;
                }}

                .gene-badge.has-stats {{
                    cursor: help;
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

                /* Sidebar Navigation */
                .sidebar {{
                    position: fixed;
                    top: 100px;
                    right: 20px;
                    width: 200px;
                    max-height: calc(100vh - 140px);
                    overflow-y: auto;
                    background: white;
                    border-radius: 8px;
                    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                    padding: 1rem;
                    z-index: 100;
                }}

                .sidebar h4 {{
                    font-size: 0.85rem;
                    color: var(--text-light);
                    margin: 0 0 0.75rem 0;
                    text-transform: uppercase;
                    letter-spacing: 0.05em;
                }}

                .sidebar ul {{
                    list-style: none;
                    padding: 0;
                    margin: 0;
                }}

                .sidebar li {{
                    margin-bottom: 0.4rem;
                }}

                .sidebar a {{
                    color: var(--text-color);
                    text-decoration: none;
                    font-size: 0.85rem;
                    display: block;
                    padding: 0.3rem 0.5rem;
                    border-radius: 4px;
                    transition: background-color 0.2s;
                }}

                .sidebar a:hover {{
                    background-color: var(--light-bg);
                    color: var(--primary-color);
                }}

                .sidebar a.final {{
                    color: #10b981;
                    font-weight: 600;
                }}

                /* Adjust body for sidebar */
                @media (min-width: 1200px) {{
                    body {{
                        margin-right: 240px;
                    }}
                }}

                @media (max-width: 1199px) {{
                    .sidebar {{
                        display: none;
                    }}
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <header>
                    <h1>CASSIA Cell Type Annotation Summary {strategy_display}</h1>
                    <p class="subtitle">Single-cell RNA-seq Analysis Report</p>
                </header>
                
                <section id="overview">
                    <h2>Overview</h2>
                    <div class="content">
                        {sections['overview']}
                    </div>
                </section>

                <section id="initial-assessment">
                    <h2>{('Initial Hypothesis' if report_style.lower() == 'total_summary' else 'Initial Assessment')}</h2>
                    <div class="content">
                        {sections.get('initial_hypothesis' if report_style.lower() == 'total_summary' else 'initial_assessment', 'No information available')}
                    </div>
                </section>
        """
        
        # Add content based on report style
        if report_style.lower() == "total_summary":
            # Add gene groups for total summary style
            for i, group in enumerate(gene_groups, 1):
                # Format genes as badges with tooltip stats
                genes = group['genes']
                gene_badges = ""
                if genes and genes != "No genes listed":
                    gene_list = [g.strip() for g in re.split(r'[,\s]+', genes) if g.strip()]
                    gene_badges = '<div class="gene-list">' + ''.join([create_gene_badge_html(gene, gene_stats) for gene in gene_list]) + '</div>'
                
                html += f"""
                    <section id="gene-analysis-{i}">
                        <h2>Gene Analysis {i}: {group['title']}</h2>

                        <div class="sub-section">
                            <h3>Genes Analyzed</h3>
                            {gene_badges}
                        </div>

                        <div class="sub-section">
                            <h3>Findings & Conclusions</h3>
                            <div class="content">
                                {group['findings']}
                            </div>
                        </div>
                    </section>
                """
        else:
            # Add iterations for per-iteration style
            for iteration in iterations:
                # Format genes checked as badges with tooltip stats
                genes = iteration['genes_checked']
                gene_badges = ""
                if genes and genes != "No information available":
                    gene_list = [g.strip() for g in re.split(r'[,\s]+', genes) if g.strip()]
                    gene_badges = '<div class="gene-list">' + ''.join([create_gene_badge_html(gene, gene_stats) for gene in gene_list]) + '</div>'
                
                html += f"""
                    <section id="iteration-{iteration['number']}">
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
                                {format_hypotheses(iteration['key_findings'])}
                            </div>
                        </div>
                    </section>
                """
        
        # Add final annotation (highlighted)
        final_section_key = 'final_conclusion' if report_style.lower() == 'total_summary' else 'final_annotation'
        final_section_title = 'Final Conclusion' if report_style.lower() == 'total_summary' else 'Final Annotation'
        html += f"""
                <section id="final-annotation" class="final-annotation">
                    <h2>{final_section_title}</h2>
                    <div class="content">
                        {sections.get(final_section_key, 'No information available')}
                    </div>
                </section>
        """

        if report_style.lower() == "total_summary":
            # Add sections specific to total summary style
            if sections.get('key_insights') and sections['key_insights'] != "No information available":
                html += f"""
                    <section id="key-insights">
                        <h2>Key Insights</h2>
                        <div class="content">
                            {sections['key_insights']}
                        </div>
                    </section>
                """

            if sections.get('validation_status') and sections['validation_status'] != "No information available":
                html += f"""
                    <section id="validation-status">
                        <h2>Validation Status</h2>
                        <div class="content">
                            {sections['validation_status']}
                        </div>
                    </section>
                """
        else:
            # Add sections specific to per-iteration style
            # Add marker summary
            if sections.get('marker_summary') and sections['marker_summary'] != "No information available":
                html += f"""
                    <section id="marker-genes">
                        <h2>Key Marker Genes</h2>
                        <div class="content">
                            {sections['marker_summary']}
                        </div>
                    </section>
                """

            # Add recommendations if available
            if sections.get('recommendations') and sections['recommendations'] != "No information available":
                html += f"""
                    <section id="recommendations">
                        <h2>Recommendations</h2>
                        <div class="content">
                            {sections['recommendations']}
                        </div>
                    </section>
                """
        
        # Build sidebar navigation
        sidebar_links = [
            '<li><a href="#overview">Overview</a></li>',
            f'<li><a href="#initial-assessment">{"Initial Hypothesis" if report_style.lower() == "total_summary" else "Initial Assessment"}</a></li>'
        ]

        if report_style.lower() == "total_summary":
            # Add gene analysis links
            for i in range(1, len(gene_groups) + 1):
                sidebar_links.append(f'<li><a href="#gene-analysis-{i}">Gene Analysis {i}</a></li>')
        else:
            # Add iteration links
            for iteration in iterations:
                sidebar_links.append(f'<li><a href="#iteration-{iteration["number"]}">Iteration {iteration["number"]}</a></li>')

        # Add final annotation link
        final_title = "Final Conclusion" if report_style.lower() == "total_summary" else "Final Annotation"
        sidebar_links.append(f'<li><a href="#final-annotation" class="final">{final_title}</a></li>')

        # Add optional section links
        if report_style.lower() == "total_summary":
            if sections.get('key_insights') and sections['key_insights'] != "No information available":
                sidebar_links.append('<li><a href="#key-insights">Key Insights</a></li>')
            if sections.get('validation_status') and sections['validation_status'] != "No information available":
                sidebar_links.append('<li><a href="#validation-status">Validation Status</a></li>')
        else:
            if sections.get('marker_summary') and sections['marker_summary'] != "No information available":
                sidebar_links.append('<li><a href="#marker-genes">Key Marker Genes</a></li>')
            if sections.get('recommendations') and sections['recommendations'] != "No information available":
                sidebar_links.append('<li><a href="#recommendations">Recommendations</a></li>')

        sidebar_html = f"""
            <nav class="sidebar">
                <h4>Contents</h4>
                <ul>
                    {''.join(sidebar_links)}
                </ul>
            </nav>
        """

        # Add smooth scrolling JavaScript
        smooth_scroll_js = """
            <script>
                document.querySelectorAll('.sidebar a').forEach(link => {
                    link.addEventListener('click', function(e) {
                        e.preventDefault();
                        const target = document.querySelector(this.getAttribute('href'));
                        if (target) {
                            target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                        }
                    });
                });
            </script>
        """

        # Close the HTML
        html += f"""
            </div>
            {sidebar_html}
            {smooth_scroll_js}
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