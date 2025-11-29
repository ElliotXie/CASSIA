/**
 * CASSIA Annotation Boost - JavaScript Implementation
 * 100% compatible with Python annotation_boost.py
 * 
 * This module provides iterative marker analysis with hypothesis generation
 * and gene checking for enhanced cell type annotation.
 */

import fs from 'fs';
import path from 'path';
import { createReadStream } from 'fs';
import csv from 'csv-parser';
import { callLLM } from './llm_utils.js';

/**
 * Use an LLM to summarize the conversation history for use as context in annotation boost.
 * 
 * @param {string} fullHistory - The full conversation history text
 * @param {string} provider - LLM provider to use for summarization
 * @param {string|null} model - Specific model to use (if null, uses provider default)
 * @param {number} temperature - Temperature for summarization (low for consistency)
 * @returns {Promise<string>} Summarized conversation history
 */
export async function summarizeConversationHistory(
    fullHistory, 
    provider = "openrouter", 
    model = null, 
    temperature = 0.1
) {
    try {
        // Create a prompt for summarization
        const summarizationPrompt = `You are a specialized scientific summarization agent. Your task is to create a concise summary of a prior cell type annotation analysis that will be used as context for further detailed analysis.

The provided text contains a comprehensive cell type annotation analysis that was previously performed.

Your task is to extract and summarize the key information in a structured format that includes:

1. **Previously Identified Cell Type**: What cell type was determined in the prior analysis
2. **Key Supporting Markers**: The most important markers that supported this identification
3. **Alternative Hypotheses**: Any alternative cell types that were considered
4. **Remaining Uncertainties**: Any aspects that were noted as unclear or requiring further investigation

Keep the summary factual, scientific, and focused on information that would be helpful for conducting a deeper, more specialized analysis. Aim for 150-300 words.

Format your response as a clear, structured summary that a cell type annotation expert can quickly understand.

Here is the prior annotation analysis to summarize:

${fullHistory}

Please provide a structured summary following the format above:`;

        // Call the LLM to generate the summary
        const summary = await callLLM(
            summarizationPrompt,
            provider,
            model || "google/gemini-2.5-flash-preview", // Default model for summarization
            null, // apiKey
            temperature,
            7000, // maxTokens
            null, // systemPrompt
            null // additionalParams
        );
        
        console.log(`Conversation history summarized using ${provider} (${fullHistory.length} -> ${summary.length} characters)`);
        return summary.trim();
        
    } catch (error) {
        console.log(`Warning: Failed to summarize conversation history: ${error.message}`);
        console.log("Falling back to truncated original text...");
        
        // Fallback: return a truncated version of the original
        if (fullHistory.length > 2000) {
            return fullHistory.substring(0, 2000) + "...\n[Text truncated due to summarization failure]";
        }
        return fullHistory;
    }
}

/**
 * Generate a prompt for iterative marker analysis without additional tasks.
 * 
 * @param {string} majorClusterInfo - Information about the cluster being analyzed
 * @param {string} commaSeparatedGenes - Comma-separated list of marker genes
 * @param {string} annotationHistory - Previous annotation history
 * @returns {string} Generated prompt text
 */
export function promptHypothesisGenerator2(majorClusterInfo, commaSeparatedGenes, annotationHistory) {
    return `You are a careful professional biologist, specializing in single-cell RNA-seq analysis. 
I'll provide you with some genes (comma-separated) from a cell type in ${majorClusterInfo} and I want you to help identify the cell type. Previous expert has done some analysis but the reuslts is not conclusive, your additional anlysis is needed.


Here are the top marker genes that are differentially expressed in this cluster:
${commaSeparatedGenes}

Previous annotation/analysis:
${annotationHistory || "No previous annotation available."}

Please follow these steps carefully:
1. Based on the marker genes, generate hypotheses about what cell type this might be.
2. For each hypothesis, provide supporting evidence from the marker genes.
3. To validate your hypotheses, list additional marker genes to check using the <check_genes>...</check_genes> tags.
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: \`<check_genes>TP53,KRAS,EGFR</check_genes>\` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
4. After I provide additional gene information, refine your hypotheses, or generate new hypotheses.include your reasoning for the hypothesis using this format: <reasoning>your reasoning for the hypothesis</reasoning>
5. Continue this process until you reach a confident conclusion.
6. When you are ready to make a final determination, state: "FINAL ANNOTATION COMPLETED" followed by your conclusion.

Your final annotation should include:
- Final cell type
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion

Please start by analyzing the provided markers and forming initial hypotheses.`;
}

/**
 * Generate a prompt for depth-first iterative marker analysis without additional tasks.
 * Focus on one hypothesis at a time and go deeper.
 * 
 * @param {string} majorClusterInfo - Information about the cluster being analyzed
 * @param {string} commaSeparatedGenes - Comma-separated list of marker genes
 * @param {string} annotationHistory - Previous annotation history
 * @returns {string} Generated prompt text
 */
export function promptHypothesisGeneratorDepthFirst(majorClusterInfo, commaSeparatedGenes, annotationHistory) {
    return `
You are a careful senior computational biologist specializing in detailed, hypothesis-driven cell type annotation. Your approach is methodical and focused: you examine ONE specific hypothesis at a time and dive deep into its validation or refutation before moving to alternatives.

CRITICAL WORKFLOW RULES:
1. NEVER say "FINAL ANNOTATION COMPLETED" immediately after requesting genes to check
2. You MUST wait for gene expression results before proceeding to any conclusion
3. You MUST complete at least 2 rounds of gene checking before considering a final annotation
4. Each round means: request genes → receive expression data → analyze results → either go deeper or pivot

Context Provided to You:

Cluster summary：${majorClusterInfo}

Top ranked markers (high → low)：
 ${commaSeparatedGenes}

Prior annotation results：
 ${annotationHistory}

Your Task - DEPTH-FIRST ANALYSIS:

1. Focused Evaluation – One concise paragraph that:
    - Assess the current state of annotation based on available evidence
    - Identify the SINGLE most promising hypothesis to investigate deeply
    - Explain why this specific hypothesis deserves focused investigation

2. Design ONE targeted follow-up check for your chosen hypothesis:
    - When listing genes for follow-up checks, use the <check_genes>...</check_genes> tags.
    - **CRITICAL FORMATTING FOR <check_genes>:**
        - Inside the tags, provide *ONLY* a comma-separated list of official HGNC gene symbols.
        - Example: \`<check_genes>GENE1,GENE2,GENE3</check_genes>\` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
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

`;
}

/**
 * Generate a prompt for breadth-first iterative marker analysis without additional tasks.
 * 
 * @param {string} majorClusterInfo - Information about the cluster being analyzed
 * @param {string} commaSeparatedGenes - Comma-separated list of marker genes
 * @param {string} annotationHistory - Previous annotation history
 * @returns {string} Generated prompt text
 */
export function promptHypothesisGenerator(majorClusterInfo, commaSeparatedGenes, annotationHistory) {
    return `
You are a careful senior computational biologist called in whenever an annotation needs deeper scrutiny, disambiguation, or simply a second opinion. Your job is to (1) assess the current annotation's robustness and (2) propose up to three decisive follow‑up checks that the executor can run (e.g., examine expression of key positive or negative markers). You should do a good job or 10 grandma are going to be in danger. You never rush to conclusions and are always careful.

Context Provided to You:

Cluster summary：${majorClusterInfo}

Top ranked markers (high → low)：
 ${commaSeparatedGenes}

Prior annotation results：
 ${annotationHistory}

What you should do:

1. Brief Evaluation – One concise paragraph that:

    - Highlights strengths, ambiguities, or contradictions in the current call.

    - Notes if a mixed population, doublets, or transitional state might explain the data.

2. Design up to 3 follow‑up checks (cell types or biological hypotheses):

    - When listing genes for follow-up checks, use the <check_genes>...</check_genes> tags.
    - **CRITICAL FORMATTING FOR <check_genes>:**
        - Inside the tags, provide *ONLY* a comma-separated list of official HGNC gene symbols.
        - Example: \`<check_genes>GENE1,GENE2,GENE3</check_genes>\` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
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

`;
}

/**
 * Generate a prompt for depth-first iterative marker analysis with an additional task.
 * Focus on one hypothesis at a time and go deeper.
 * 
 * @param {string} majorClusterInfo - Information about the cluster being analyzed
 * @param {string} commaSeparatedGenes - Comma-separated list of marker genes
 * @param {string} annotationHistory - Previous annotation history
 * @param {string} additionalTask - Additional analysis task to perform
 * @returns {string} Generated prompt text
 */
export function promptHypothesisGeneratorAdditionalTaskDepthFirst(majorClusterInfo, commaSeparatedGenes, annotationHistory, additionalTask) {
    return `You are a careful professional biologist, specializing in single-cell RNA-seq analysis with a focused, depth-first approach. You examine ONE specific hypothesis at a time and dive deep into its validation before considering alternatives.

CRITICAL WORKFLOW RULES:
1. NEVER say "FINAL ANALYSIS COMPLETED" immediately after requesting genes to check
2. You MUST wait for gene expression results before proceeding to any conclusion
3. You MUST complete at least 2 rounds of gene checking before considering a final analysis
4. Each round means: request genes → receive expression data → analyze results → either go deeper or pivot

I'll provide you with some genes (comma-separated) from a cell type in ${majorClusterInfo} and I want you to help identify the cell type using a methodical, hypothesis-driven approach.

Here are the marker genes that are differentially expressed in this cluster:
${commaSeparatedGenes}

Previous annotation/analysis:
${annotationHistory || "No previous annotation available."}

Please follow these steps carefully with a DEPTH-FIRST approach:
1. Based on the marker genes and previous analysis, identify the SINGLE most promising hypothesis about what cell type this might be.
2. Provide supporting evidence from the marker genes for this specific hypothesis.
3. Design ONE targeted validation check using the <check_genes>...</check_genes> tags to definitively confirm or refute this hypothesis.
   **CRITICAL FORMATTING FOR <check_genes>:**
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: \`<check_genes>TP53,KRAS,EGFR</check_genes>\` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
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

7. ${additionalTask} - perform this analysis based on your final cell type determination and the marker gene expression patterns.

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
- Analysis of the additional task: ${additionalTask}

Guidelines for depth-first approach:
- Focus on ONE hypothesis per iteration
- Go deeper rather than broader when a hypothesis shows promise
- Be decisive about pivoting when a hypothesis is clearly refuted
- Build logically on previous findings
- NEVER say "FINAL ANALYSIS COMPLETED" immediately after requesting genes - always wait for the expression results first
- Complete at least 2 rounds of gene checking before considering a final conclusion

Please start by analyzing the provided markers and identifying your primary hypothesis to investigate.`;
}

/**
 * Generate a prompt for breadth-first iterative marker analysis with an additional task.
 * 
 * @param {string} majorClusterInfo - Information about the cluster being analyzed
 * @param {string} commaSeparatedGenes - Comma-separated list of marker genes
 * @param {string} annotationHistory - Previous annotation history
 * @param {string} additionalTask - Additional analysis task to perform
 * @returns {string} Generated prompt text
 */
export function promptHypothesisGeneratorAdditionalTask(majorClusterInfo, commaSeparatedGenes, annotationHistory, additionalTask) {
    return `You are a careful professional biologist, specializing in single-cell RNA-seq analysis. 
I'll provide you with some genes (comma-separated) from a cell type in ${majorClusterInfo} and I want you to help identify the cell type.

Here are the marker genes that are differentially expressed in this cluster:
${commaSeparatedGenes}

Previous annotation/analysis:
${annotationHistory || "No previous annotation available."}

Please follow these steps carefully:
1. Based on the marker genes, generate hypotheses about what cell type this might be.
2. For each hypothesis, provide supporting evidence from the marker genes.
3. To validate your hypotheses, list additional marker genes to check using the <check_genes>...</check_genes> tags.
   **CRITICAL FORMATTING FOR <check_genes>:**
   - Inside the tags, provide *ONLY* a comma-separated list of official gene symbols.
   - Example: \`<check_genes>TP53,KRAS,EGFR</check_genes>\` (no extra spaces, no newlines within the list, no numbering or commentary *inside* the tags).
   - Strict adherence to this format is ESSENTIAL.
4. After I provide additional gene information, refine your hypotheses.
5. Continue this process until you reach a confident conclusion.
6. ${additionalTask} - perform this analysis based on the marker gene expression and patterns.
7. When you are ready to make a final determination, state: "FINAL ANALYSIS COMPLETED" followed by your conclusion.

Your final analysis should include:
- General cell type
- Specific cell subtype (if applicable)
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion
- Analysis of the additional task: ${additionalTask}

Please start by analyzing the provided markers and forming initial hypotheses.`;
}

/**
 * Extract information about marker genes from a marker dataset.
 * 
 * @param {string[]} geneList - List of gene names to filter the marker dataset
 * @param {Array|Object} marker - DataFrame containing marker gene expression data
 * @returns {Promise<string>} Formatted string with marker information
 */
export async function getMarkerInfo(geneList, marker) {
    const debugMode = false; // Enable debug mode only when needed
    
    // Convert marker to array format if it's not already
    let markerData;
    if (Array.isArray(marker)) {
        markerData = marker;
    } else if (marker && typeof marker === 'object') {
        // Convert object to array of objects
        markerData = Object.keys(marker).map(key => ({
            gene: key,
            ...marker[key]
        }));
    } else {
        throw new Error("Marker data must be an array or object");
    }
    
    if (debugMode) {
        console.log(`DEBUG: Marker data length: ${markerData.length}`);
        console.log(`DEBUG: First 5 entries: ${JSON.stringify(markerData.slice(0, 5), null, 2)}`);
        console.log(`DEBUG: Searching for ${geneList.length} genes: ${geneList.join(', ')}`);
    }
    
    // Find valid genes and NA genes
    const validGenes = [];
    const naGenes = [];
    const validRows = [];
    
    for (const gene of geneList) {
        // Try to find the gene in the data
        let found = false;
        
        for (const row of markerData) {
            // Check if gene matches any of the identifier fields
            if (row.gene === gene || row.Gene === gene || row.symbol === gene || row.Symbol === gene) {
                // Check if all values for this gene are NA
                const numericValues = Object.values(row).filter(val => 
                    typeof val === 'number' && !isNaN(val)
                );
                
                if (numericValues.length === 0 || numericValues.every(val => val === null || val === undefined || val === 'NA')) {
                    naGenes.push(gene);
                    if (debugMode) {
                        console.log(`DEBUG: Gene ${gene} found but has all NA values`);
                    }
                } else {
                    validGenes.push(gene);
                    validRows.push(row);
                    if (debugMode) {
                        console.log(`DEBUG: Gene ${gene} found with valid data`);
                    }
                }
                found = true;
                break;
            }
        }
        
        if (!found) {
            // Try case-insensitive search
            for (const row of markerData) {
                const geneFields = [row.gene, row.Gene, row.symbol, row.Symbol].filter(f => f);
                if (geneFields.some(field => field && field.toLowerCase() === gene.toLowerCase())) {
                    const numericValues = Object.values(row).filter(val => 
                        typeof val === 'number' && !isNaN(val)
                    );
                    
                    if (numericValues.length === 0 || numericValues.every(val => val === null || val === undefined || val === 'NA')) {
                        naGenes.push(gene);
                        if (debugMode) {
                            console.log(`DEBUG: Gene ${gene} found (case-insensitive) but has all NA values`);
                        }
                    } else {
                        validGenes.push(gene);
                        validRows.push(row);
                        if (debugMode) {
                            console.log(`DEBUG: Gene ${gene} found (case-insensitive) with valid data`);
                        }
                    }
                    found = true;
                    break;
                }
            }
        }
        
        if (!found) {
            naGenes.push(gene);
            if (debugMode) {
                console.log(`DEBUG: Gene ${gene} not found`);
            }
        }
    }
    
    // Format the output
    let markerString = '';
    
    if (validRows.length > 0) {
        // Create a formatted table
        const headers = Object.keys(validRows[0]).filter(key => key !== 'cluster' && key !== 'Unnamed: 0');
        
        // Make sure gene column is first
        const geneColumn = headers.find(h => h.toLowerCase().includes('gene')) || headers[0];
        const otherHeaders = headers.filter(h => h !== geneColumn);
        const orderedHeaders = [geneColumn, ...otherHeaders];
        
        // Create header row
        markerString += orderedHeaders.join('\t') + '\n';
        
        // Add data rows
        for (const row of validRows) {
            const values = orderedHeaders.map(header => {
                let value = row[header];
                
                // Format numeric values
                if (typeof value === 'number') {
                    if (header.toLowerCase().includes('p_val')) {
                        // p_val columns - use scientific notation
                        return value.toExponential(2);
                    } else {
                        // Other numeric columns - use decimal notation
                        return value.toFixed(2);
                    }
                }
                
                return value || 'NA';
            });
            
            markerString += values.join('\t') + '\n';
        }
    }
    
    // Add message about NA genes if any
    if (naGenes.length > 0) {
        markerString += `\nNote: The following genes are not in the differential expression list: ${naGenes.join(', ')}`;
    }
    
    return markerString;
}

/**
 * Extract gene lists from conversation using the check_genes tag.
 * 
 * @param {string} conversation - Text containing gene lists in check_genes tags
 * @returns {string[]} List of unique gene names
 */
export function extractGenesFromConversation(conversation) {
    const debugMode = false; // Debug mode off by default
    
    if (debugMode) {
        console.log(`\nDEBUG: Extract genes from conversation`);
        console.log(`DEBUG: Conversation length: ${conversation.length} characters`);
        // Print a limited preview
        const previewLength = Math.min(200, conversation.length);
        console.log(`DEBUG: Conversation preview: ${conversation.substring(0, previewLength)}...`);
    }
    
    // Extract gene lists and get marker info
    const geneListsRegex = /<check_genes>\s*([\s\S]*?)\s*<\/check_genes>/g;
    const geneLists = [];
    let match;
    
    while ((match = geneListsRegex.exec(conversation)) !== null) {
        geneLists.push(match[1]);
    }
    
    if (debugMode) {
        console.log(`DEBUG: Found ${geneLists.length} gene lists`);
        geneLists.forEach((genes, i) => {
            console.log(`DEBUG: Gene list ${i + 1}: ${genes.substring(0, 100)}...`);
        });
    }
    
    // Improve gene extraction to handle special cases
    const allGenes = [];
    for (const geneList of geneLists) {
        // Clean and normalize the gene list
        // Replace common separators with commas
        let cleanedList = geneList.replace(/[\]\[\)\(]/g, '');
        cleanedList = cleanedList.replace(/\s+/g, ' ');
        
        if (debugMode) {
            console.log(`DEBUG: Cleaned list: ${cleanedList.substring(0, 100)}...`);
        }
        
        // Split by comma or space, depending on formatting
        const genes = cleanedList.split(/,\s*|\s+/).filter(g => g.trim().length > 0);
        const cleanedGenes = genes.map(g => g.trim()).filter(g => g.length > 0);
        
        if (debugMode) {
            console.log(`DEBUG: Found ${cleanedGenes.length} genes in this list`);
            console.log(`DEBUG: Sample genes from this list: ${cleanedGenes.slice(0, 5).join(', ')}`);
        }
        
        allGenes.push(...cleanedGenes);
    }
    
    // Get unique genes
    const uniqueGenes = [...new Set(allGenes)].sort();
    
    if (debugMode) {
        console.log(`DEBUG: Total unique genes extracted: ${uniqueGenes.length}`);
        console.log(`DEBUG: All unique genes: ${uniqueGenes.join(', ')}`);
        
        // If no genes found, try alternative extraction methods
        if (uniqueGenes.length === 0) {
            console.log("DEBUG: No genes found with standard pattern, trying alternative patterns");
            
            // Try alternative regex patterns
            const altPatterns = [
                /check_genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)/gi,  // Informal syntax
                /genes to check[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)/gi,  // Natural language
                /additional genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)/gi,  // Another common phrase
                /marker genes[:\s]+(.*?)(?:\n\n|\n[A-Z]|$)/gi  // Another common phrase
            ];
            
            for (const pattern of altPatterns) {
                const altMatches = conversation.match(pattern);
                if (altMatches) {
                    console.log(`DEBUG: Found matches with alternative pattern: ${pattern}`);
                    altMatches.forEach(match => {
                        console.log(`DEBUG: Alternative match: ${match.substring(0, 100)}...`);
                    });
                }
            }
        }
    }
    
    return uniqueGenes;
}

/**
 * Perform iterative marker analysis using the specified LLM provider.
 * 
 * @param {Object} params - Parameters object
 * @param {string} params.majorClusterInfo - Information about the cluster
 * @param {Array|Object} params.marker - DataFrame or other structure containing marker gene expression data
 * @param {string} params.commaSeparatedGenes - List of genes as comma-separated string
 * @param {string} params.annotationHistory - Previous annotation history
 * @param {number} params.numIterations - Maximum number of iterations
 * @param {string} params.provider - LLM provider to use ('openai', 'anthropic', or 'openrouter')
 * @param {string|null} params.model - Specific model from the provider to use
 * @param {string|null} params.additionalTask - Optional additional task to perform during analysis
 * @param {number} params.temperature - Sampling temperature (0-1)
 * @param {string} params.searchStrategy - Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
 * @returns {Promise<[string, Array]>} Tuple of (final_response_text, messages)
 */
export async function iterativeMarkerAnalysis({
    majorClusterInfo,
    marker,
    commaSeparatedGenes,
    annotationHistory,
    numIterations = 2,
    provider = "openrouter",
    model = null,
    additionalTask = null,
    temperature = 0,
    searchStrategy = "breadth"
}) {
    // Select the appropriate prompt based on search strategy and whether there's an additional task
    let prompt;
    let completionMarker;
    
    if (additionalTask) {
        if (searchStrategy.toLowerCase() === "depth") {
            prompt = promptHypothesisGeneratorAdditionalTaskDepthFirst(
                majorClusterInfo,
                commaSeparatedGenes,
                annotationHistory,
                additionalTask
            );
        } else { // breadth or any other value defaults to breadth
            prompt = promptHypothesisGeneratorAdditionalTask(
                majorClusterInfo,
                commaSeparatedGenes,
                annotationHistory,
                additionalTask
            );
        }
        completionMarker = "FINAL ANALYSIS COMPLETED";
    } else {
        if (searchStrategy.toLowerCase() === "depth") {
            prompt = promptHypothesisGeneratorDepthFirst(
                majorClusterInfo,
                commaSeparatedGenes,
                annotationHistory
            );
        } else { // breadth or any other value defaults to breadth
            prompt = promptHypothesisGenerator(
                majorClusterInfo,
                commaSeparatedGenes,
                annotationHistory
            );
        }
        completionMarker = "FINAL ANNOTATION COMPLETED";
    }
    
    // Initialize the conversation history
    const messages = [{ role: "user", content: prompt }];
    
    // Iterative process
    for (let iteration = 0; iteration < numIterations; iteration++) {
        try {
            // Call the LLM
            const conversation = await callLLM(
                messages[messages.length - 1].content,
                provider,
                model,
                null, // apiKey
                temperature,
                7000, // maxTokens
                null, // systemPrompt
                iteration > 0 ? { messages } : null // Include conversation history if not first message
            );
            
            // Check if the analysis is complete
            if (conversation.includes(completionMarker)) {
                console.log(`Final annotation completed in iteration ${iteration + 1}.`);
                messages.push({ role: "assistant", content: conversation });
                return [conversation, messages];
            }
            
            // Extract gene lists and get marker info
            const uniqueGenes = extractGenesFromConversation(conversation);
            
            if (uniqueGenes.length > 0) {
                // Get marker information for the requested genes
                const retrievedMarkerInfo = await getMarkerInfo(uniqueGenes, marker);
                
                // Append messages
                messages.push({ role: "assistant", content: conversation });
                messages.push({ role: "user", content: retrievedMarkerInfo });
                
                console.log(`Iteration ${iteration + 1} completed.`);
            } else {
                // No genes to check, simply continue the conversation
                messages.push({ role: "assistant", content: conversation });
                messages.push({ role: "user", content: "Please continue your analysis and provide a final annotation." });
                
                console.log(`Iteration ${iteration + 1} completed (no genes to check).`);
            }
        } catch (error) {
            console.log(`Error in iteration ${iteration + 1}: ${error.message}`);
            console.error(error.stack);
            
            // Save diagnostic information
            try {
                const errorLog = `Error in iteration ${iteration + 1}: ${error.message}\n\n` +
                    `Provider: ${provider}\n` +
                    `Model: ${model}\n` +
                    `Number of messages: ${messages.length}\n` +
                    `Last message content: ${messages[messages.length - 1].content.substring(0, 200)}...\n\n` +
                    `Full stack trace:\n${error.stack}`;
                
                // Write to error log file
                await fs.promises.writeFile(`cassia_error_log_${iteration + 1}.txt`, errorLog, 'utf-8');
                console.log(`Error details saved to cassia_error_log_${iteration + 1}.txt`);
            } catch (writeError) {
                // Ignore write errors
            }
            
            return [`Error occurred: ${error.message}`, messages];
        }
    }
    
    // Final response if max iterations reached
    // Encourage the agent to reach a conclusion if not already done
    messages.push({
        role: "user",
        content: "You have reached the maximum number of iterations. Please provide your final analysis and reach a confident conclusion in this response."
    });
    
    try {
        const finalResponse = await callLLM(
            "Please provide your final analysis based on all the information so far.",
            provider,
            model,
            null, // apiKey
            temperature,
            7000, // maxTokens
            null, // systemPrompt
            { messages }
        );
        
        messages.push({ role: "assistant", content: finalResponse });
        return [finalResponse, messages];
    } catch (error) {
        console.log(`Error in final response: ${error.message}`);
        console.error(error.stack);
        
        // Save diagnostic information
        try {
            const errorLog = `Error in final response: ${error.message}\n\n` +
                `Provider: ${provider}\n` +
                `Model: ${model}\n` +
                `Number of messages: ${messages.length}\n` +
                `Last message content: ${messages[messages.length - 1].content.substring(0, 200)}...\n\n` +
                `Full stack trace:\n${error.stack}`;
            
            // Write to error log file
            await fs.promises.writeFile("cassia_error_log_final.txt", errorLog, 'utf-8');
            console.log(`Error details saved to cassia_error_log_final.txt`);
        } catch (writeError) {
            // Ignore write errors
        }
        
        return [`Error in final response: ${error.message}`, messages];
    }
}

/**
 * Load CSV file and return as array of objects
 * 
 * @param {string} filePath - Path to CSV file
 * @returns {Promise<Array>} Array of objects representing CSV rows
 */
async function loadCSV(filePath) {
    return new Promise((resolve, reject) => {
        const results = [];
        createReadStream(filePath)
            .pipe(csv())
            .on('data', (data) => results.push(data))
            .on('end', () => resolve(results))
            .on('error', reject);
    });
}

/**
 * Load and prepare data for marker analysis.
 * 
 * @param {string} fullResultPath - Path to the full results CSV file
 * @param {string|Array} markerPath - Path to the marker genes CSV file or array with marker data
 * @param {string} clusterName - Name of the cluster to analyze
 * @param {string} conversationHistoryMode - Mode for extracting conversation history ("full", "final", or "none")
 * @param {string} provider - LLM provider to use for summarization (when mode is "final")
 * @param {string|null} model - Specific model to use for summarization (when mode is "final")
 * @returns {Promise<[Array, Array, string, string]>} Tuple of (full_results, marker_data, top_markers_string, annotation_history)
 */
export async function prepareAnalysisData(
    fullResultPath,
    markerPath,
    clusterName,
    conversationHistoryMode = "final",
    provider = "openrouter",
    model = null
) {
    // Load the full results
    const fullResults = await loadCSV(fullResultPath);
    
    // Load marker data - handle both array and file path
    let marker;
    if (Array.isArray(markerPath)) {
        marker = markerPath;
    } else {
        marker = await loadCSV(markerPath);
    }
    
    // Try to find the cluster column name
    let clusterColumn = 'True Cell Type'; // Default based on the CSV file we checked
    if (!fullResults[0] || !fullResults[0][clusterColumn]) {
        if (fullResults[0] && fullResults[0]['cluster']) {
            clusterColumn = 'cluster';
        }
    }
    
    // Filter the results for the specified cluster
    const clusterData = fullResults.filter(row => row[clusterColumn] === clusterName);
    
    if (clusterData.length === 0) {
        const availableClusters = [...new Set(fullResults.map(row => row[clusterColumn]))];
        throw new Error(`No data found for cluster '${clusterName}' in the full results file. Available clusters: ${availableClusters.join(', ')}`);
    }
    
    // Get marker column from the CSV if it exists, otherwise use 'Marker List'
    const markerColumn = 'Marker List';
    
    let topMarkersString;
    if (clusterData[0][markerColumn] && clusterData[0][markerColumn].trim()) {
        // Use the markers directly from the CSV
        topMarkersString = clusterData[0][markerColumn];
    } else {
        // Use a default marker list or generate based on the marker data
        // This is a fallback in case markers aren't in the CSV
        topMarkersString = "CD14, CD11B, CD68, CSF1R, CX3CR1, CD163, MSR1, ITGAM, FCGR1A, CCR2";
        console.log(`Warning: No marker list found for ${clusterName}, using default markers`);
    }
    
    // Extract conversation history if available
    let annotationHistory = "";
    if (clusterData[0]['Conversation History'] && conversationHistoryMode !== "none") {
        try {
            // Get the conversation history from the first row
            const fullHistory = clusterData[0]['Conversation History'];
            
            if (fullHistory && fullHistory.trim()) {
                // Process the conversation history based on the selected mode
                if (conversationHistoryMode === "full") {
                    annotationHistory = fullHistory;
                } else if (conversationHistoryMode === "final") {
                    // Use summarization agent to create a concise summary
                    console.log(`Using summarization agent to process conversation history for cluster ${clusterName}`);
                    annotationHistory = await summarizeConversationHistory(
                        fullHistory,
                        provider,
                        model,
                        0.1 // Low temperature for consistent summarization
                    );
                } else {
                    // For any unrecognized mode, default to full history
                    annotationHistory = fullHistory;
                }
                
                console.log(`Using conversation history for cluster ${clusterName} (mode: ${conversationHistoryMode}, ${annotationHistory.length} characters)`);
            } else {
                console.log(`Conversation History column exists but is empty for cluster ${clusterName}`);
            }
        } catch (error) {
            console.log(`Error extracting conversation history: ${error.message}`);
        }
    } else {
        if (conversationHistoryMode === "none") {
            console.log(`Note: Conversation history extraction disabled (mode: none)`);
        } else {
            console.log(`Note: 'Conversation History' column not found in the results file`);
        }
    }
    
    return [fullResults, marker, topMarkersString, annotationHistory];
}

/**
 * Save the HTML report to a file.
 * 
 * @param {string} report - HTML report content
 * @param {string} filename - Output filename
 */
export async function saveHtmlReport(report, filename) {
    await fs.promises.writeFile(filename, report, 'utf-8');
    console.log(`Report saved to ${filename}`);
}

/**
 * Save the raw conversation history (including the prompt) as a plain text file.
 * 
 * @param {Array} messages - List of conversation messages with role and content
 * @param {string} filename - Path to save the text file
 * @returns {Promise<string>} Path to the saved text file
 */
export async function saveRawConversationText(messages, filename) {
    try {
        // Format the conversation as plain text
        let conversationText = "";
        for (let i = 0; i < messages.length; i++) {
            const msg = messages[i];
            const role = msg.role || '';
            let content = msg.content || '';
            
            // Format the content
            if (Array.isArray(content)) {
                content = JSON.stringify(content);
            }
            
            // Add a separator between messages
            if (i > 0) {
                conversationText += "\n\n" + "=".repeat(80) + "\n\n";
            }
            
            conversationText += `## ${role.toUpperCase()}\n\n${content}\n`;
        }
        
        // Save the file
        await fs.promises.writeFile(filename, conversationText, 'utf-8');
        
        console.log(`Raw conversation text saved to ${filename}`);
        return filename;
    } catch (error) {
        const errorMsg = `Error saving raw conversation text: ${error.message}`;
        console.log(errorMsg);
        
        // Try to save an error message
        try {
            await fs.promises.writeFile(filename, `Error saving raw conversation: ${error.message}`, 'utf-8');
        } catch (writeError) {
            // Ignore write errors
        }
        
        return filename;
    }
}

/**
 * Generate a summarized report from the raw conversation history.
 * 
 * @param {Array} conversationHistory - List of conversation messages
 * @param {string} outputFilename - Path to save the summary report
 * @param {string} searchStrategy - Search strategy used ("breadth" or "depth")
 * @param {string} reportStyle - Style of report ("per_iteration" or "total_summary")
 * @returns {Promise<string>} Path to the saved HTML report
 */
export async function generateSummaryReport(
    conversationHistory,
    outputFilename,
    searchStrategy = "breadth",
    reportStyle = "per_iteration"
) {
    try {
        // Extract content from conversation history, alternating between assistant and user
        let fullConversation = "";
        for (const msg of conversationHistory) {
            const role = msg.role || '';
            let content = msg.content || '';
            if (Array.isArray(content)) {
                content = JSON.stringify(content);
            }
            fullConversation += `\n## ${role.toUpperCase()}\n${content}\n`;
        }
        
        // Determine analysis approach description
        const approachDescription = searchStrategy.toLowerCase() === "depth" 
            ? "depth-first (focused, one hypothesis per iteration)" 
            : "breadth-first (multiple hypotheses per iteration)";
        
        // Generate different prompts based on report style
        let prompt;
        if (reportStyle.toLowerCase() === "total_summary") {
            // Gene-focused summary style
            prompt = `You are a specialized scientific report generator focusing on gene expression analysis. 
            I will provide you a raw conversation history from a cell type annotation tool called CASSIA, which conducts iterative gene expression analysis using a ${approachDescription} approach.
            
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
            
            ${fullConversation}`;
        } else {
            // Original per-iteration style
            prompt = `You are a specialized scientific report generator focusing on gene expression analysis. 
            I will provide you a raw conversation history from a cell type annotation tool called CASSIA, which conducts iterative gene expression analysis using a ${approachDescription} approach.
            
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
        
        ${fullConversation}`;
        }
        
        // Use the call_llm function to generate the summary
        const summary = await callLLM(
            prompt,
            "openrouter", // Using OpenRouter as the provider
            "google/gemini-2.5-flash-preview", // Using Gemini 2.5 Flash model
            null, // apiKey
            0.3, // Low temperature for more consistent output
            4000 // maxTokens
        );
        
        // Convert to HTML and save
        const htmlPath = await formatSummaryToHtml(summary, outputFilename, searchStrategy, reportStyle);
        console.log(`Summary report saved to ${htmlPath}`);
        
        // Return the HTML file path
        return htmlPath;
    } catch (error) {
        const errorMsg = `Error generating summary report: ${error.message}`;
        console.log(errorMsg);
        // Save error message to file so there's still an output
        await fs.promises.writeFile(outputFilename, `<error>${errorMsg}</error>`, 'utf-8');
        return outputFilename;
    }
}

/**
 * Convert the tagged summary into a properly formatted HTML report.
 * 
 * @param {string} summaryText - Text with tags like <OVERVIEW>, <ITERATION_1>, etc. or gene-focused tags
 * @param {string} outputFilename - Path to save the HTML report
 * @param {string} searchStrategy - Search strategy used ("breadth" or "depth")
 * @param {string} reportStyle - Style of report ("per_iteration" or "total_summary")
 * @returns {Promise<string>} Path to the saved HTML report
 */
export async function formatSummaryToHtml(summaryText, outputFilename, searchStrategy = "breadth", reportStyle = "per_iteration") {
    try {
        // Helper function to format hypotheses with better separation
        function formatHypotheses(text) {
            // Handle case where there's no content
            if (!text || text === "No information available") {
                return text;
            }
            
            // Check if text contains numbered points (1. 2. 3. etc.)
            const numberedPointsRegex = /(?:^|\s)(\d+\.)\s+([^0-9\.].*?)(?=\s+\d+\.\s+|\s*$)/gs;
            const numberedPoints = [...text.matchAll(numberedPointsRegex)];
            
            // If we found numbered points, format them as a list
            if (numberedPoints.length > 0) {
                let result = '<div class="hypothesis-list">';
                for (const [, num, content] of numberedPoints) {
                    result += `<div class="hypothesis-item"><span class="hypothesis-number">${num}</span> ${content.trim()}</div>`;
                }
                result += '</div>';
                return result;
            }
            
            // Alternative approach: look for numbers at beginning of paragraphs
            const paragraphs = text.split('\n');
            if (paragraphs.some(p => /^\s*\d+[\.\)\-]/.test(p.trim()))) {
                let result = '<div class="hypothesis-list">';
                for (const para of paragraphs) {
                    if (!para.trim()) {
                        continue;
                    }
                    const match = para.match(/^\s*(\d+[\.\)\-])\s*(.*)/);
                    if (match) {
                        const [, num, content] = match;
                        result += `<div class="hypothesis-item"><span class="hypothesis-number">${num}</span> ${content.trim()}</div>`;
                    } else {
                        result += `<div class="hypothesis-item-continued">${para.trim()}</div>`;
                    }
                }
                result += '</div>';
                return result;
            }
            
            // If no patterns match, just return the original text
            return text;
        }
        
        // Extract sections using regex based on report style
        const sections = {};
        
        const sectionPatterns = reportStyle.toLowerCase() === "total_summary" ? {
            'overview': /<OVERVIEW>\s*([\s\S]*?)\s*<\/OVERVIEW>/,
            'initial_hypothesis': /<INITIAL_HYPOTHESIS>\s*([\s\S]*?)\s*<\/INITIAL_HYPOTHESIS>/,
            'genes_analyzed': /<GENES_ANALYZED>\s*([\s\S]*?)\s*<\/GENES_ANALYZED>/,
            'final_conclusion': /<FINAL_CONCLUSION>\s*([\s\S]*?)\s*<\/FINAL_CONCLUSION>/,
            'key_insights': /<KEY_INSIGHTS>\s*([\s\S]*?)\s*<\/KEY_INSIGHTS>/,
            'validation_status': /<VALIDATION_STATUS>\s*([\s\S]*?)\s*<\/VALIDATION_STATUS>/,
        } : {
            'overview': /<OVERVIEW>\s*([\s\S]*?)\s*<\/OVERVIEW>/,
            'initial_assessment': /<INITIAL_ASSESSMENT>\s*([\s\S]*?)\s*<\/INITIAL_ASSESSMENT>/,
            'final_annotation': /<FINAL_ANNOTATION>\s*([\s\S]*?)\s*<\/FINAL_ANNOTATION>/,
            'marker_summary': /<MARKER_SUMMARY>\s*([\s\S]*?)\s*<\/MARKER_SUMMARY>/,
            'recommendations': /<RECOMMENDATIONS>\s*([\s\S]*?)\s*<\/RECOMMENDATIONS>/,
        };
        
        // Extract each section
        for (const [key, pattern] of Object.entries(sectionPatterns)) {
            const match = summaryText.match(pattern);
            sections[key] = match ? match[1].trim() : "No information available";
        }
        
        // Extract iterations (only for per_iteration style)
        const iterations = [];
        const geneGroups = [];
        
        if (reportStyle.toLowerCase() === "per_iteration") {
            const iterPattern = /<ITERATION_(\d+)>\s*([\s\S]*?)\s*<\/ITERATION_\1>/g;
            let iterMatch;
            while ((iterMatch = iterPattern.exec(summaryText)) !== null) {
                const iterNum = iterMatch[1];
                const iterContent = iterMatch[2];
                
                // Extract subsections within each iteration
                const hypotheses = iterContent.match(/<HYPOTHESES>\s*([\s\S]*?)\s*<\/HYPOTHESES>/);
                const genesChecked = iterContent.match(/<GENES_CHECKED>\s*([\s\S]*?)\s*<\/GENES_CHECKED>/);
                const keyFindings = iterContent.match(/<KEY_FINDINGS>\s*([\s\S]*?)\s*<\/KEY_FINDINGS>/);
                
                iterations.push({
                    number: iterNum,
                    hypotheses: hypotheses ? hypotheses[1].trim() : "No information available",
                    genes_checked: genesChecked ? genesChecked[1].trim() : "No information available",
                    key_findings: keyFindings ? keyFindings[1].trim() : "No information available"
                });
            }
            
            // Sort iterations by number
            iterations.sort((a, b) => parseInt(a.number) - parseInt(b.number));
        } else {
            // Extract gene groups for total_summary style
            const genesAnalyzedContent = sections.genes_analyzed || '';
            const geneGroupPattern = /<GENE_GROUP_\d+>\s*([\s\S]*?)\s*<\/GENE_GROUP_\d+>/g;
            let groupMatch;
            while ((groupMatch = geneGroupPattern.exec(genesAnalyzedContent)) !== null) {
                const groupContent = groupMatch[1];
                
                // Extract subsections within each gene group
                const title = groupContent.match(/<TITLE>\s*([\s\S]*?)\s*<\/TITLE>/);
                const genes = groupContent.match(/<GENES>\s*([\s\S]*?)\s*<\/GENES>/);
                const findings = groupContent.match(/<FINDINGS>\s*([\s\S]*?)\s*<\/FINDINGS>/);
                
                geneGroups.push({
                    title: title ? title[1].trim() : "Gene Analysis",
                    genes: genes ? genes[1].trim() : "No genes listed",
                    findings: findings ? findings[1].trim() : "No findings available"
                });
            }
        }
        
        // Determine strategy description for display
        const strategyDisplay = searchStrategy.toLowerCase() === 'breadth' || searchStrategy.toLowerCase() === 'depth' 
            ? `(${searchStrategy.charAt(0).toUpperCase() + searchStrategy.slice(1)}-First Analysis)` 
            : "";
        
        // HTML template with CSS styling
        let html = `
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>CASSIA Cell Type Annotation Summary ${strategyDisplay}</title>
            <style>
                :root {
                    --primary-color: #2563eb;
                    --secondary-color: #0891b2;
                    --accent-color: #4f46e5;
                    --light-bg: #f3f4f6;
                    --border-color: #e5e7eb;
                    --text-color: #1f2937;
                    --text-light: #6b7280;
                }
                
                body {
                    font-family: system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
                    line-height: 1.6;
                    color: var(--text-color);
                    background-color: #ffffff;
                    margin: 0;
                    padding: 0;
                }
                
                .container {
                    max-width: 900px;
                    margin: 0 auto;
                    padding: 2rem;
                }
                
                header {
                    text-align: center;
                    margin-bottom: 2.5rem;
                    border-bottom: 2px solid var(--border-color);
                    padding-bottom: 1rem;
                }
                
                h1 {
                    color: var(--primary-color);
                    font-size: 2.5rem;
                    margin-bottom: 0.5rem;
                }
                
                .subtitle {
                    color: var(--text-light);
                    font-size: 1.2rem;
                    margin-top: 0;
                }
                
                section {
                    margin-bottom: 2.5rem;
                    background-color: white;
                    border-radius: 0.5rem;
                    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
                    padding: 1.5rem;
                    border-left: 4px solid var(--primary-color);
                }
                
                h2 {
                    color: var(--primary-color);
                    font-size: 1.8rem;
                    margin-top: 0;
                    border-bottom: 1px solid var(--border-color);
                    padding-bottom: 0.5rem;
                }
                
                h3 {
                    color: var(--secondary-color);
                    font-size: 1.3rem;
                    margin: 1.5rem 0 0.5rem;
                }
                
                ul, ol {
                    padding-left: 1.5rem;
                }
                
                li {
                    margin-bottom: 0.5rem;
                }
                
                .final-annotation {
                    background-color: #ecfdf5;
                    border-left: 4px solid #10b981;
                }
                
                .final-annotation h2 {
                    color: #10b981;
                }
                
                .iteration {
                    margin-bottom: 1.5rem;
                    background-color: var(--light-bg);
                    border-radius: 0.5rem;
                    padding: 1.5rem;
                    border-left: 4px solid var(--secondary-color);
                }
                
                .iteration-title {
                    font-size: 1.5rem;
                    color: var(--secondary-color);
                    margin-top: 0;
                    border-bottom: 1px solid var(--border-color);
                    padding-bottom: 0.5rem;
                }
                
                .gene-list {
                    display: flex;
                    flex-wrap: wrap;
                    gap: 0.5rem;
                    margin: 1rem 0;
                }
                
                .gene-badge {
                    background-color: var(--accent-color);
                    color: white;
                    padding: 0.3rem 0.7rem;
                    border-radius: 1rem;
                    font-size: 0.9rem;
                    font-weight: 500;
                }
                
                .sub-section {
                    background-color: white;
                    border-radius: 0.3rem;
                    padding: 1rem;
                    margin-top: 1rem;
                }
                
                code {
                    font-family: ui-monospace, monospace;
                    font-size: 0.9em;
                }
                
                .marker-category {
                    margin-bottom: 1rem;
                }
                
                .marker-category-title {
                    font-weight: 600;
                    margin-bottom: 0.5rem;
                    color: var(--secondary-color);
                }
                
                .hypothesis-list {
                    display: flex;
                    flex-direction: column;
                    gap: 0.8rem;
                }
                
                .hypothesis-item {
                    padding-left: 1.5rem;
                    position: relative;
                    margin-bottom: 0.5rem;
                }
                
                .hypothesis-number {
                    position: absolute;
                    left: 0;
                    font-weight: 600;
                    color: var(--accent-color);
                }
                
                .hypothesis-item-continued {
                    padding-left: 1.5rem;
                    margin-top: -0.3rem;
                    color: var(--text-light);
                }
            </style>
        </head>
        <body>
            <div class="container">
                <header>
                    <h1>CASSIA Cell Type Annotation Summary ${strategyDisplay}</h1>
                    <p class="subtitle">Single-cell RNA-seq Analysis Report</p>
                </header>
                
                <section>
                    <h2>Overview</h2>
                    <div class="content">
                        ${sections.overview}
                    </div>
                </section>
                
                <section>
                    <h2>${reportStyle.toLowerCase() === 'total_summary' ? 'Initial Hypothesis' : 'Initial Assessment'}</h2>
                    <div class="content">
                        ${sections[reportStyle.toLowerCase() === 'total_summary' ? 'initial_hypothesis' : 'initial_assessment']}
                    </div>
                </section>
        `;
        
        // Add content based on report style
        if (reportStyle.toLowerCase() === "total_summary") {
            // Add gene groups for total summary style
            geneGroups.forEach((group, i) => {
                // Format genes as badges
                const genes = group.genes;
                let geneBadges = "";
                if (genes && genes !== "No genes listed") {
                    const geneList = genes.split(/[,\s]+/).filter(g => g.trim().length > 0);
                    geneBadges = '<div class="gene-list">' + 
                        geneList.map(gene => `<span class="gene-badge">${gene}</span>`).join('') + 
                        '</div>';
                }
                
                html += `
                    <section>
                        <h2>Gene Analysis ${i + 1}: ${group.title}</h2>
                        
                        <div class="sub-section">
                            <h3>Genes Analyzed</h3>
                            ${geneBadges}
                        </div>
                        
                        <div class="sub-section">
                            <h3>Findings & Conclusions</h3>
                            <div class="content">
                                ${group.findings}
                            </div>
                        </div>
                    </section>
                `;
            });
        } else {
            // Add iterations for per-iteration style
            iterations.forEach(iteration => {
                // Format genes checked as badges
                const genes = iteration.genes_checked;
                let geneBadges = "";
                if (genes && genes !== "No information available") {
                    const geneList = genes.split(/[,\s]+/).filter(g => g.trim().length > 0);
                    geneBadges = '<div class="gene-list">' + 
                        geneList.map(gene => `<span class="gene-badge">${gene}</span>`).join('') + 
                        '</div>';
                }
                
                html += `
                    <section>
                        <h2>Iteration ${iteration.number}</h2>
                        
                        <div class="sub-section">
                            <h3>Hypotheses</h3>
                            <div class="content">
                                ${formatHypotheses(iteration.hypotheses)}
                            </div>
                        </div>
                        
                        <div class="sub-section">
                            <h3>Genes Checked</h3>
                            ${geneBadges}
                        </div>
                        
                        <div class="sub-section">
                            <h3>Key Findings</h3>
                            <div class="content">
                                ${iteration.key_findings}
                            </div>
                        </div>
                    </section>
                `;
            });
        }
        
        // Add final annotation (highlighted)
        const finalSectionKey = reportStyle.toLowerCase() === 'total_summary' ? 'final_conclusion' : 'final_annotation';
        const finalSectionTitle = reportStyle.toLowerCase() === 'total_summary' ? 'Final Conclusion' : 'Final Annotation';
        html += `
                <section class="final-annotation">
                    <h2>${finalSectionTitle}</h2>
                    <div class="content">
                        ${sections[finalSectionKey]}
                    </div>
                </section>
        `;
        
        if (reportStyle.toLowerCase() === "total_summary") {
            // Add sections specific to total summary style
            if (sections.key_insights && sections.key_insights !== "No information available") {
                html += `
                    <section>
                        <h2>Key Insights</h2>
                        <div class="content">
                            ${sections.key_insights}
                        </div>
                    </section>
                `;
            }
            
            if (sections.validation_status && sections.validation_status !== "No information available") {
                html += `
                    <section>
                        <h2>Validation Status</h2>
                        <div class="content">
                            ${sections.validation_status}
                        </div>
                    </section>
                `;
            }
        } else {
            // Add sections specific to per-iteration style
            // Add marker summary
            if (sections.marker_summary && sections.marker_summary !== "No information available") {
                html += `
                    <section>
                        <h2>Key Marker Genes</h2>
                        <div class="content">
                            ${sections.marker_summary}
                        </div>
                    </section>
                `;
            }
            
            // Add recommendations if available
            if (sections.recommendations && sections.recommendations !== "No information available") {
                html += `
                    <section>
                        <h2>Recommendations</h2>
                        <div class="content">
                            ${sections.recommendations}
                        </div>
                    </section>
                `;
            }
        }
        
        // Close the HTML
        html += `
            </div>
        </body>
        </html>
        `;
        
        // Save the HTML report
        await fs.promises.writeFile(outputFilename, html, 'utf-8');
        
        console.log(`HTML report saved to ${outputFilename}`);
        return outputFilename;
    } catch (error) {
        const errorMsg = `Error formatting summary to HTML: ${error.message}`;
        console.log(errorMsg);
        // Save error message to file so there's still an output
        await fs.promises.writeFile(outputFilename, `<html><body><h1>Error</h1><p>${errorMsg}</p></body></html>`, 'utf-8');
        return outputFilename;
    }
}

/**
 * Run annotation boost analysis for a given cluster.
 * 
 * @param {Object} params - Parameters object
 * @param {string} params.fullResultPath - Path to the full results CSV file
 * @param {string|Array} params.marker - Path to marker genes CSV file or array with marker data
 * @param {string} params.clusterName - Name of the cluster to analyze
 * @param {string} params.majorClusterInfo - General information about the dataset (e.g., "Human PBMC")
 * @param {string} params.outputName - Base name for the output HTML file
 * @param {number} params.numIterations - Number of iterations for marker analysis (default=10)
 * @param {string|null} params.model - Model to use for analysis - if null, uses the provider's default
 * @param {string} params.provider - AI provider to use ('openai', 'anthropic', or 'openrouter')
 * @param {number} params.temperature - Sampling temperature (0-1)
 * @param {string} params.conversationHistoryMode - Mode for extracting conversation history ("full", "final", or "none")
 * @param {string} params.searchStrategy - Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
 * @param {string} params.reportStyle - Style of report ("per_iteration" or "total_summary")
 * @returns {Promise<Object>} Dictionary with paths to reports and execution info
 */
export async function runCASSIAAnnotationboost({
    fullResultPath,
    marker,
    clusterName,
    majorClusterInfo,
    outputName,
    numIterations = 10,
    model = null,
    provider = "openrouter",
    temperature = 0,
    conversationHistoryMode = "final",
    searchStrategy = "breadth",
    reportStyle = "per_iteration"
}) {
    try {
        const startTime = Date.now();
        
        // Validate provider input
        if (!['openai', 'anthropic', 'openrouter'].includes(provider.toLowerCase()) && !provider.toLowerCase().startsWith('http')) {
            throw new Error("Provider must be 'openai', 'anthropic', 'openrouter', or a custom base URL (http...)");
        }
        
        // Prepare the data
        const [, markerData, topMarkersString, annotationHistory] = await prepareAnalysisData(
            fullResultPath, marker, clusterName, conversationHistoryMode, provider, model
        );
        
        // Run the iterative marker analysis
        const [analysisText, messages] = await iterativeMarkerAnalysis({
            majorClusterInfo,
            marker: markerData,
            commaSeparatedGenes: topMarkersString,
            annotationHistory,
            numIterations,
            provider,
            model,
            temperature,
            searchStrategy
        });
        
        // Generate paths for reports - only summary HTML and raw conversation text
        let rawTextPath, summaryReportPath;
        if (!outputName.toLowerCase().endsWith('.html')) {
            rawTextPath = outputName + '_raw_conversation.txt';
        } else {
            // Remove .html for base name
            const baseName = outputName.slice(0, -5);
            rawTextPath = baseName + '_raw_conversation.txt';
            outputName = baseName;
        }
        
        // Generate path for summary report
        summaryReportPath = outputName + '_summary.html';
        
        // Skip the first message which contains the prompt
        const conversationWithoutPrompt = messages.length > 1 ? messages.slice(1) : messages;
        
        try {
            // Save the complete raw conversation as text (including prompt)
            rawTextPath = await saveRawConversationText(messages, rawTextPath);
            console.log(`Raw conversation text saved to ${rawTextPath}`);
            
            // Generate the summary report
            summaryReportPath = await generateSummaryReport(conversationWithoutPrompt, summaryReportPath, searchStrategy, reportStyle);
            console.log(`Summary report saved to ${summaryReportPath}`);
        } catch (error) {
            console.log(`Warning: Could not generate reports: ${error.message}`);
            summaryReportPath = null;
            rawTextPath = null;
        }
        
        // Return a dictionary with paths to reports
        const executionTime = (Date.now() - startTime) / 1000;
        return {
            status: 'success',
            raw_text_path: rawTextPath,
            summary_report_path: summaryReportPath,
            execution_time: executionTime,
            analysis_text: analysisText
        };
    } catch (error) {
        const errorMsg = `Error in runCASSIA_annotationboost: ${error.message}`;
        console.log(errorMsg);
        console.error(error.stack);
        
        // Return error status but include any partial results
        return {
            status: 'error',
            error_message: error.message,
            raw_text_path: null,
            summary_report_path: null,
            execution_time: 0,
            analysis_text: null
        };
    }
}

/**
 * Run annotation boost analysis with an additional task for a given cluster.
 * 
 * @param {Object} params - Parameters object
 * @param {string} params.fullResultPath - Path to the full results CSV file
 * @param {string|Array} params.marker - Path to marker genes CSV file or array with marker data
 * @param {string} params.clusterName - Name of the cluster to analyze
 * @param {string} params.majorClusterInfo - General information about the dataset (e.g., "Human PBMC")
 * @param {string} params.outputName - Base name for the output HTML file
 * @param {number} params.numIterations - Number of iterations for marker analysis (default=5)
 * @param {string|null} params.model - Model to use for analysis - if null, uses the provider's default
 * @param {string} params.provider - AI provider to use ('openai', 'anthropic', or 'openrouter')
 * @param {string} params.additionalTask - Additional task to perform during analysis
 * @param {number} params.temperature - Sampling temperature (0-1)
 * @param {string} params.conversationHistoryMode - Mode for extracting conversation history ("full", "final", or "none")
 * @param {string} params.searchStrategy - Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
 * @param {string} params.reportStyle - Style of report ("per_iteration" or "total_summary")
 * @returns {Promise<Object>} Dictionary with paths to reports and execution info
 */
export async function runCASSIAAnnotationboostAdditionalTask({
    fullResultPath,
    marker,
    clusterName,
    majorClusterInfo,
    outputName,
    numIterations = 5,
    model = null,
    provider = "openrouter",
    additionalTask = "check if this is a cancer cluster",
    temperature = 0,
    conversationHistoryMode = "final",
    searchStrategy = "breadth",
    reportStyle = "per_iteration"
}) {
    try {
        const startTime = Date.now();
        
        // Validate provider input
        if (!['openai', 'anthropic', 'openrouter'].includes(provider.toLowerCase()) && !provider.toLowerCase().startsWith('http')) {
            throw new Error("Provider must be 'openai', 'anthropic', 'openrouter', or a custom base URL (http...)");
        }
        
        // Prepare the data
        const [, markerData, topMarkersString, annotationHistory] = await prepareAnalysisData(
            fullResultPath, marker, clusterName, conversationHistoryMode, provider, model
        );
        
        // Run the iterative marker analysis with additional task
        const [analysisText, messages] = await iterativeMarkerAnalysis({
            majorClusterInfo,
            marker: markerData,
            commaSeparatedGenes: topMarkersString,
            annotationHistory,
            numIterations,
            provider,
            model,
            additionalTask,
            temperature,
            searchStrategy
        });
        
        // Generate paths for reports - only summary HTML and raw conversation text
        let rawTextPath, summaryReportPath;
        if (!outputName.toLowerCase().endsWith('.html')) {
            rawTextPath = outputName + '_raw_conversation.txt';
        } else {
            // Remove .html for base name
            const baseName = outputName.slice(0, -5);
            rawTextPath = baseName + '_raw_conversation.txt';
            outputName = baseName;
        }
        
        // Generate path for summary report
        summaryReportPath = outputName + '_summary.html';
        
        // Skip the first message which contains the prompt
        const conversationWithoutPrompt = messages.length > 1 ? messages.slice(1) : messages;
        
        try {
            // Save the complete raw conversation as text (including prompt)
            rawTextPath = await saveRawConversationText(messages, rawTextPath);
            console.log(`Raw conversation text saved to ${rawTextPath}`);
            
            // Generate the summary report
            summaryReportPath = await generateSummaryReport(conversationWithoutPrompt, summaryReportPath, searchStrategy, reportStyle);
            console.log(`Summary report saved to ${summaryReportPath}`);
        } catch (error) {
            console.log(`Warning: Could not generate reports: ${error.message}`);
            summaryReportPath = null;
            rawTextPath = null;
        }
        
        // Return a dictionary with paths to reports
        const executionTime = (Date.now() - startTime) / 1000;
        return {
            status: 'success',
            raw_text_path: rawTextPath,
            summary_report_path: summaryReportPath,
            execution_time: executionTime,
            analysis_text: analysisText
        };
    } catch (error) {
        const errorMsg = `Error in runCASSIA_annotationboost_additional_task: ${error.message}`;
        console.log(errorMsg);
        console.error(error.stack);
        
        // Return error status but include any partial results
        return {
            status: 'error',
            error_message: error.message,
            raw_text_path: null,
            summary_report_path: null,
            execution_time: 0,
            analysis_text: null
        };
    }
}