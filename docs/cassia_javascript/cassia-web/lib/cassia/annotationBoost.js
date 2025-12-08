/**
 * CASSIA Annotation Boost - JavaScript Implementation for Browser
 * This module provides iterative marker analysis with hypothesis generation
 * and gene checking for enhanced cell type annotation.
 */
import { callLLM } from './llm_utils.js';
import { parseCSV } from '../utils/csv-parser.js';

// ----------------- CSV Processing Functions -----------------

/**
 * Extract top marker genes from marker data (FindAllMarkers/scanpy output)
 * @param {Array} markerData - Array of marker gene objects from CSV
 * @param {number} topN - Number of top genes to extract (default: 20)
 * @returns {string} Comma-separated list of top marker genes
 */
export function extractTopMarkerGenes(markerData, topN = 20) {
    if (!markerData || markerData.length === 0) {
        throw new Error('No marker data provided');
    }
    
    console.log(`🧬 Extracting top ${topN} marker genes from ${markerData.length} total markers`);
    
    // Try to find gene column with common names
    const geneColumnOptions = ['gene', 'Gene', 'GENE', 'gene_name', 'Gene_name', 'symbol', 'Symbol', 'SYMBOL'];
    let geneColumn = null;
    
    // Check which gene column exists
    for (const col of geneColumnOptions) {
        if (markerData[0] && col in markerData[0]) {
            geneColumn = col;
            break;
        }
    }
    
    if (!geneColumn) {
        // Fallback: use first column if no standard gene column found
        const firstRow = markerData[0];
        const columns = Object.keys(firstRow);
        geneColumn = columns[0];
        console.log(`⚠️ No standard gene column found, using first column: "${geneColumn}"`);
    }
    
    console.log(`🎯 Using gene column: "${geneColumn}"`);
    
    // Extract top genes (assumes data is already sorted by significance)
    const topGenes = markerData
        .slice(0, topN)
        .map(row => row[geneColumn])
        .filter(gene => gene && gene.trim() !== '') // Remove empty genes
        .join(', ');
    
    console.log(`✅ Extracted genes: ${topGenes}`);
    return topGenes;
}

/**
 * Extract conversation history from runCASSIA batch full.csv for a specific cluster
 * @param {string|Array} csvData - CSV content as string or array of objects
 * @param {string} clusterName - Name of the cluster to extract conversation for
 * @param {string} clusterColumnName - Name of the column containing cluster identifiers (default: 'True Cell Type')
 * @returns {string} Conversation history for the specified cluster
 */
export function extractConversationForCluster(csvData, clusterName, clusterColumnName = 'True Cell Type') {
    try {
        // Parse CSV if it's a string using the robust parser
        let rows;
        if (typeof csvData === 'string') {
            rows = parseCSV(csvData, { debug: false });
        } else if (Array.isArray(csvData)) {
            rows = csvData;
        } else {
            throw new Error('csvData must be either a string or array');
        }
        
        if (rows.length === 0) {
            throw new Error('No data found in CSV');
        }
        
        // Find conversation history column with comprehensive options
        const historyColumnOptions = [
            'Conversation History', 'Conversation.History', 'conversation_history', 
            'Conversation_History', 'Conversation-History', 'ConversationHistory',
            'conversationhistory', 'conversation history', 'CONVERSATION HISTORY',
            'annotation_history', 'Annotation History', 'Annotation.History',
            'annotation_conversation', 'Annotation Conversation',
            'conversation', 'Conversation', 'history', 'History',
            'chat_history', 'Chat History', 'chat-history', 'ChatHistory'
        ];
        
        let historyColumn = null;
        // First try exact match
        for (const col of historyColumnOptions) {
            if (col in rows[0]) {
                historyColumn = col;
                break;
            }
        }
        
        // If no exact match, try case-insensitive and normalized match
        if (!historyColumn) {
            const availableColumns = Object.keys(rows[0]);
            for (const availCol of availableColumns) {
                const normalizedAvail = availCol.toLowerCase().replace(/[\s._-]+/g, '');
                for (const option of historyColumnOptions) {
                    const normalizedOption = option.toLowerCase().replace(/[\s._-]+/g, '');
                    if (normalizedAvail === normalizedOption || 
                        (normalizedAvail.includes('conversation') && normalizedAvail.includes('history')) ||
                        normalizedAvail.includes('conversationhistory') ||
                        normalizedAvail.includes('annotationhistory')) {
                        historyColumn = availCol;
                        break;
                    }
                }
                if (historyColumn) break;
            }
        }
        
        if (!historyColumn) {
            const availableColumns = Object.keys(rows[0]);
            throw new Error(`No conversation history column found. Available columns: ${availableColumns.join(', ')}`);
        }
        
        // Use the specified cluster column or try to find it if not specified
        let clusterColumn = clusterColumnName;
        
        // Verify the cluster column exists
        if (!(clusterColumn in rows[0])) {
            // Try to find alternative cluster column names
            const clusterColumnOptions = [
                'True Cell Type', 'True.Cell.Type', 'Cluster', 'cluster', 
                'Cell Type', 'Cell.Type', 'celltype', 'CellType',
                'Predicted Main Cell Type', 'Main Cell Type', 'main_cell_type',
                'Annotation', 'annotation', 'Label', 'label'
            ];
            
            let foundClusterColumn = null;
            for (const col of clusterColumnOptions) {
                if (col in rows[0]) {
                    foundClusterColumn = col;
                    break;
                }
            }
            
            if (foundClusterColumn) {
                clusterColumn = foundClusterColumn;
                console.log(`⚠️ Specified cluster column "${clusterColumnName}" not found. Using "${foundClusterColumn}" instead.`);
            } else {
                const availableColumns = Object.keys(rows[0]);
                throw new Error(`Cluster column "${clusterColumnName}" not found. Available columns: ${availableColumns.join(', ')}`);
            }
        }
        
        console.log(`🔍 Using cluster column: "${clusterColumn}", history column: "${historyColumn}"`);
        
        // Find the row for the specified cluster
        const targetRow = rows.find(row => {
            const rowCluster = row[clusterColumn];
            return rowCluster && rowCluster.toLowerCase().trim() === clusterName.toLowerCase().trim();
        });
        
        if (!targetRow) {
            const availableClusters = rows.map(row => row[clusterColumn]).filter(c => c).slice(0, 10);
            throw new Error(`Cluster "${clusterName}" not found in column "${clusterColumn}". Available clusters (first 10): ${availableClusters.join(', ')}`);
        }
        
        const conversationHistory = targetRow[historyColumn];
        if (!conversationHistory || conversationHistory.trim() === '') {
            throw new Error(`No conversation history found for cluster "${clusterName}" in column "${historyColumn}"`);
        }
        
        console.log(`✅ Extracted conversation history for "${clusterName}": ${conversationHistory.length} characters`);
        return conversationHistory.trim();
        
    } catch (error) {
        console.error('Error extracting conversation for cluster:', error);
        throw error;
    }
}

/**
 * Get available clusters from runCASSIA batch CSV
 * @param {string|Array} csvData - CSV content as string or array of objects
 * @param {string} clusterColumnName - Name of the column containing cluster identifiers (default: 'True Cell Type')
 * @returns {Array} Array of cluster names
 */
export function getAvailableClusters(csvData, clusterColumnName = 'True Cell Type') {
    try {
        let rows;
        if (typeof csvData === 'string') {
            rows = parseCSV(csvData, { debug: false });
        } else if (Array.isArray(csvData)) {
            rows = csvData;
        } else {
            return [];
        }
        
        if (rows.length === 0) return [];
        
        // Use the specified cluster column or try to find it if not specified
        let clusterColumn = clusterColumnName;
        
        // Verify the cluster column exists
        if (!(clusterColumn in rows[0])) {
            const clusterColumnOptions = [
                'True Cell Type', 'True.Cell.Type', 'Cluster', 'cluster', 
                'Cell Type', 'Cell.Type', 'celltype', 'CellType',
                'Predicted Main Cell Type', 'Main Cell Type', 'main_cell_type',
                'Annotation', 'annotation', 'Label', 'label'
            ];
            
            let foundClusterColumn = null;
            for (const col of clusterColumnOptions) {
                if (col in rows[0]) {
                    foundClusterColumn = col;
                    break;
                }
            }
            
            if (foundClusterColumn) {
                clusterColumn = foundClusterColumn;
            } else {
                return [];
            }
        }
        
        return rows.map(row => row[clusterColumn]).filter(c => c && c.trim() !== '');
        
    } catch (error) {
        console.error('Error getting available clusters:', error);
        return [];
    }
}

/**
 * Get available column names from CSV data
 * @param {string|Array} csvData - CSV content as string or array of objects
 * @returns {Array} Array of column names
 */
export function getAvailableColumns(csvData) {
    try {
        let rows;
        if (typeof csvData === 'string') {
            rows = parseCSV(csvData, { debug: false });
        } else if (Array.isArray(csvData)) {
            rows = csvData;
        } else {
            return [];
        }
        
        if (rows.length === 0) return [];
        
        return Object.keys(rows[0]);
        
    } catch (error) {
        console.error('Error getting available columns:', error);
        return [];
    }
}

/**
 * Extract marker genes for a specific cluster from runCASSIA CSV data
 * @param {string|Array} csvData - CSV content as string or array of objects
 * @param {string} clusterName - Name of the cluster to extract markers for
 * @param {string} clusterColumnName - Name of the column containing cluster identifiers (default: 'True Cell Type')
 * @returns {string} Comma-separated list of marker genes for the cluster
 */
export function extractMarkerGenesForCluster(csvData, clusterName, clusterColumnName = 'True Cell Type') {
    try {
        // Parse CSV if it's a string using the robust parser
        let rows;
        if (typeof csvData === 'string') {
            rows = parseCSV(csvData, { debug: false });
        } else if (Array.isArray(csvData)) {
            rows = csvData;
        } else {
            throw new Error('csvData must be either a string or array');
        }
        
        if (rows.length === 0) {
            throw new Error('No data found in CSV');
        }
        
        // Find the row for the specified cluster
        const targetRow = rows.find(row => {
            const rowCluster = row[clusterColumnName];
            return rowCluster && rowCluster.toLowerCase().trim() === clusterName.toLowerCase().trim();
        });
        
        if (!targetRow) {
            const availableClusters = rows.map(row => row[clusterColumnName]).filter(c => c).slice(0, 10);
            throw new Error(`Cluster "${clusterName}" not found in column "${clusterColumnName}". Available clusters (first 10): ${availableClusters.join(', ')}`);
        }
        
        // Extract marker genes from "Marker List" column
        const markerColumnOptions = [
            'Marker List', 'Marker.List', 'marker_list', 'Marker_List', 'Marker-List',
            'MarkerList', 'markerlist', 'marker list', 'MARKER LIST',
            'markers', 'Markers', 'marker', 'Marker', 'genes', 'Genes'
        ];
        
        let markerColumn = null;
        for (const col of markerColumnOptions) {
            if (col in targetRow && targetRow[col] && targetRow[col].trim() !== '') {
                markerColumn = col;
                break;
            }
        }
        
        if (!markerColumn) {
            const availableColumns = Object.keys(targetRow);
            throw new Error(`No marker column found for cluster "${clusterName}". Available columns: ${availableColumns.join(', ')}`);
        }
        
        const markerGenes = targetRow[markerColumn];
        if (!markerGenes || markerGenes.trim() === '') {
            throw new Error(`No marker genes found for cluster "${clusterName}" in column "${markerColumn}"`);
        }
        
        console.log(`✅ Extracted marker genes for "${clusterName}" from column "${markerColumn}": ${markerGenes.length} characters`);
        return markerGenes.trim();
        
    } catch (error) {
        console.error('Error extracting marker genes for cluster:', error);
        throw error;
    }
}

/**
 * Use an LLM to summarize the conversation history for use as context in annotation boost.
 */
export async function summarizeConversationHistory(
    fullHistory,
    apiKey,
    provider = "openrouter",
    model = null,
    temperature = 0.1,
    reasoningEffort = null
) {
    try {
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
        const summary = await callLLM(
            summarizationPrompt,
            provider,
            model || "google/gemini-2.5-flash", // Match Python's default model
            apiKey,
            temperature
            // Python doesn't specify max_tokens for this call - uses default
        );
        return summary.trim();
    } catch (error) {
        console.warn(`Warning: Failed to summarize conversation history: ${error.message}`);
        if (fullHistory.length > 2000) {
            return fullHistory.substring(0, 2000) + "...\n[Text truncated due to summarization failure]";
        }
        return fullHistory;
    }
}


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
 * Get marker information for a list of genes.
 * This is a browser-version, operating on a data object instead of a file.
 */
export async function getMarkerInfo(geneList, markerData) {
    if (!markerData || markerData.length === 0) {
        return "Marker data not provided or is empty.";
    }

    // Find gene column with comprehensive options
    const headers = Object.keys(markerData[0]);
    let geneColumn = null;
    
    // Prefer 'gene' column over any others
    if (headers.includes('gene')) {
        geneColumn = 'gene';
    } else if (headers.includes('Unnamed: 0')) {
        geneColumn = 'Unnamed: 0';
    } else {
        const geneColumnOptions = ['gene', 'Gene', 'GENE', 'gene_name', 'Gene_name', 'symbol', 'Symbol', 'SYMBOL', 'genes', 'Genes'];
        geneColumn = headers.find(h => geneColumnOptions.includes(h));
    }

    if (!geneColumn) {
        return `Could not find a 'gene' column in the marker data. Available columns: ${headers.join(', ')}`;
    }

    // Track valid genes and NA genes
    const validGenes = [];
    const naGenes = [];
    const resultRows = [];

    for (const gene of geneList) {
        // Try exact match first
        let foundRow = markerData.find(row => row[geneColumn] === gene);
        
        // Try case-insensitive match if exact match fails
        if (!foundRow) {
            foundRow = markerData.find(row => row[geneColumn] && row[geneColumn].toLowerCase() === gene.toLowerCase());
        }

        if (foundRow) {
            // Check if all values are NA or empty
            const numericColumns = Object.keys(foundRow).filter(key => 
                key !== geneColumn && 
                key !== 'cluster' && 
                key !== 'Unnamed: 0' &&
                !isNaN(parseFloat(foundRow[key])) && 
                foundRow[key] !== '' && 
                foundRow[key] !== 'NA'
            );
            
            if (numericColumns.length === 0) {
                naGenes.push(gene);
            } else {
                validGenes.push(gene);
                // Create a clean row for output - exclude unnecessary columns
                const cleanRow = { gene: gene };
                for (const col of Object.keys(foundRow)) {
                    // Skip unwanted columns: gene columns, cluster, rownames, and p_val (keep p_val_adj)
                    if (col !== 'cluster' && 
                        col !== 'Unnamed: 0' && 
                        col !== geneColumn && 
                        col !== 'gene' &&
                        col !== 'Gene' &&
                        col !== 'GENE' &&
                        col !== 'p_val' && // Remove regular p-value, keep adjusted
                        !col.toLowerCase().includes('rownames')) {
                        
                        let value = foundRow[col];
                        // Format numeric values
                        if (!isNaN(parseFloat(value)) && value !== '' && value !== 'NA') {
                            if (col.toLowerCase().includes('p_val_adj') || col.toLowerCase().includes('padj')) {
                                value = parseFloat(value).toExponential(2);
                            } else {
                                value = parseFloat(value).toFixed(2);
                            }
                        }
                        cleanRow[col] = value;
                    }
                }
                resultRows.push(cleanRow);
            }
        } else {
            naGenes.push(gene);
        }
    }

    // Generate formatted output
    let resultString = "";
    
    if (resultRows.length > 0) {
        // Create table-like output
        const columns = Object.keys(resultRows[0]);
        const columnWidths = {};
        
        // Calculate column widths
        for (const col of columns) {
            columnWidths[col] = Math.max(col.length, ...resultRows.map(row => String(row[col]).length));
        }
        
        // Format as table
        resultString = columns.map(col => col.padEnd(columnWidths[col])).join('  ') + '\n';
        for (const row of resultRows) {
            resultString += columns.map(col => String(row[col]).padEnd(columnWidths[col])).join('  ') + '\n';
        }
    }

    // Add note about missing genes
    if (naGenes.length > 0) {
        resultString += `\nNote: The following genes are not in the differential expression list: ${naGenes.join(', ')}`;
    }

    return resultString || "No gene expression data found for the requested genes.";
}

/**
 * Extracts comma-separated gene lists from <check_genes> tags in a conversation.
 */
export function extractGenesFromConversation(conversation) {
    const regex = /<check_genes>([\s\S]*?)<\/check_genes>/g;
    const matches = [...conversation.matchAll(regex)];
    const allGenes = new Set();

    for (const match of matches) {
        const geneString = match[1];
        geneString.split(',')
            .map(g => g.trim())
            .filter(g => g)
            .forEach(g => allGenes.add(g));
    }

    return Array.from(allGenes);
}

/**
 * Perform iterative marker analysis using the specified LLM provider.
 * 
 * @param {string} majorClusterInfo - Information about the cluster
 * @param {any[]} marker - DataFrame or other structure containing marker gene expression data
 * @param {string} commaSeparatedGenes - List of genes as comma-separated string
 * @param {string} annotationHistory - Previous annotation history
 * @param {number} numIterations - Maximum number of iterations
 * @param {string} provider - LLM provider to use ('openai', 'anthropic', or 'openrouter')
 * @param {string} model - Specific model from the provider to use
 * @param {string} additionalTask - Optional additional task to perform during analysis
 * @param {number} temperature - Sampling temperature (0-1)
 * @param {string} searchStrategy - Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time)
 * @param {string} apiKey - API key for browser version (additional parameter)
 * @returns {string} The final conversation text
 */
export async function iterativeMarkerAnalysis(
    majorClusterInfo,
    marker,
    commaSeparatedGenes,
    annotationHistory,
    numIterations = 2,
    provider = "openrouter",
    model = null,
    additionalTask = null,
    temperature = 0,
    searchStrategy = "breadth",
    apiKey, // Additional parameter for browser version
    reasoningEffort = null // Reasoning effort for supported models
) {

    // Select the appropriate prompt based on search strategy and whether there's an additional task
    let prompt;
    if (additionalTask) {
        // TODO: Add additional task prompts if needed (not implemented in current Python version)
        if (searchStrategy.toLowerCase() === "depth") {
            prompt = promptHypothesisGeneratorDepthFirst(majorClusterInfo, commaSeparatedGenes, annotationHistory);
        } else {
            prompt = promptHypothesisGenerator(majorClusterInfo, commaSeparatedGenes, annotationHistory);
        }
        // completion_marker = "FINAL ANALYSIS COMPLETED"  // Would be used for additional tasks
    } else {
        if (searchStrategy.toLowerCase() === "depth") {
            prompt = promptHypothesisGeneratorDepthFirst(majorClusterInfo, commaSeparatedGenes, annotationHistory);
        } else {
            prompt = promptHypothesisGenerator(majorClusterInfo, commaSeparatedGenes, annotationHistory);
        }
    }
    
    const completionMarker = "FINAL ANNOTATION COMPLETED";
    
    // Initialize the conversation history
    let conversation = `System: Starting iterative analysis for ${majorClusterInfo}.\n`;
    let messages = [{"role": "user", "content": prompt}];
    
    // Iterative process
    for (let iteration = 0; iteration < numIterations; iteration++) {
        try {
            // Call the LLM
            // Only include conversation history after first iteration (matching Python implementation)
            console.log(`🔄 Iteration ${iteration + 1}: Sending ${messages.length} messages to LLM for context`);
            const llmResponse = await callLLM(
                messages[messages.length - 1].content,
                provider,
                model,
                apiKey,
                temperature,
                4096, // Match Python's max_tokens
                null, // systemPrompt
                // Only include conversation history if not the first iteration (matching Python)
                iteration > 0 ? { messages: messages } : null,
                reasoningEffort && reasoningEffort !== 'none' ? { effort: reasoningEffort } : null
            );
            
            conversation += `\n--- Iteration ${iteration + 1} ---\n`;
            conversation += `\nExpert Analyst:\n${llmResponse}\n`;
            
            // Check if the analysis is complete
            if (llmResponse.includes(completionMarker)) {
                console.log(`🎯 Final annotation completed in iteration ${iteration + 1}. Total messages: ${messages.length + 1}`);
                messages.push({"role": "assistant", "content": llmResponse});
                return { conversation, messages };
            }
            
            // Extract gene lists and get marker info
            const uniqueGenes = extractGenesFromConversation(llmResponse);
            
            if (uniqueGenes.length > 0) {
                // Get marker information for the requested genes
                const retrievedMarkerInfo = await getMarkerInfo(uniqueGenes, marker);
                
                conversation += `\nSystem:\n${retrievedMarkerInfo}\n`;
                
                // Append messages for conversation history
                messages.push({"role": "assistant", "content": llmResponse});
                messages.push({"role": "user", "content": retrievedMarkerInfo});
                
                console.log(`✅ Iteration ${iteration + 1} completed. Total messages in history: ${messages.length}`);
            } else {
                // No genes to check, simply continue the conversation
                conversation += "\nSystem: No genes requested. Please continue your analysis and provide a final annotation.\n";
                
                messages.push({"role": "assistant", "content": llmResponse});
                messages.push({"role": "user", "content": "Please continue your analysis and provide a final annotation."});
                
                console.log(`✅ Iteration ${iteration + 1} completed (no genes to check). Total messages in history: ${messages.length}`);
            }
            
        } catch (error) {
            console.error(`Error in iteration ${iteration + 1}:`, error);
            conversation += `\nSystem: Error occurred in iteration ${iteration + 1}: ${error.message}\n`;
            return { conversation, messages };
        }
    }
    
    // Final response if max iterations reached
    // Encourage the agent to reach a conclusion if not already done
    try {
        const finalUserMessage = "You have reached the maximum number of iterations. Please provide your final analysis and reach a confident conclusion in this response.";
        conversation += `\nSystem: ${finalUserMessage}\n`;

        messages.push({"role": "user", "content": finalUserMessage});

        // Match Python: use fixed prompt text for the actual LLM call
        const finalResponse = await callLLM(
            "Please provide your final analysis based on all the information so far.", // Match Python's exact prompt
            provider,
            model,
            apiKey,
            temperature,
            4096, // Match Python's max_tokens
            null, // systemPrompt
            { messages: messages }, // Include full conversation history
            reasoningEffort && reasoningEffort !== 'none' ? { effort: reasoningEffort } : null
        );
        
        conversation += `\nExpert Analyst:\n${finalResponse}\n`;
        messages.push({"role": "assistant", "content": finalResponse});
        
        return { conversation, messages };
        
    } catch (error) {
        console.error(`Error in final response:`, error);
        conversation += `\nSystem: Error in final response: ${error.message}\n`;
        return { conversation, messages };
    }
}

/**
 * Generate a structured summary report using LLM, matching Python implementation
 * @param {Array} messages - Array of conversation messages with role and content
 * @param {string} searchStrategy - Search strategy used ("breadth" or "depth")
 * @param {string} reportStyle - Style of report ("per_iteration" or "total_summary")
 * @param {string} provider - LLM provider
 * @param {string} model - LLM model
 * @param {string} apiKey - API key
 * @returns {Promise<string>} HTML report content
 */
export async function generateSummaryReport(messages, searchStrategy = "breadth", reportStyle = "per_iteration", provider = "openrouter", model = "anthropic/claude-3.5-sonnet", apiKey, reasoningEffort = null) {
    try {
        // Ensure messages is an array
        if (!Array.isArray(messages)) {
            console.error('Messages is not an array:', typeof messages, messages);
            throw new Error('Messages parameter must be an array');
        }
        
        // Extract content from conversation history, alternating between assistant and user
        let fullConversation = "";
        for (const msg of messages) {
            const role = msg.role || '';
            const content = msg.content || '';
            fullConversation += `\n## ${role.toUpperCase()}\n${content}\n`;
        }
        
        // Determine analysis approach description
        const approachDescription = searchStrategy.toLowerCase() === "depth" 
            ? "depth-first (focused, one hypothesis per iteration)" 
            : "breadth-first (multiple hypotheses per iteration)";
        
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
        
        // Use the LLM to generate the summary
        const summary = await callLLM(
            prompt,
            provider,
            model,
            apiKey,
            0.3, // Low temperature for consistent output
            4096, // Match Python's max_tokens
            null, // systemPrompt
            null, // additionalParams
            reasoningEffort && reasoningEffort !== 'none' ? { effort: reasoningEffort } : null
        );
        
        // Convert to HTML and return
        return formatSummaryToHTML(summary, searchStrategy, reportStyle);
        
    } catch (error) {
        console.error('Error generating summary report:', error);
        // Return a basic conversation display as fallback
        return generateBasicConversationReport(messages);
    }
}

/**
 * Generate a basic conversation report as fallback
 */
function generateBasicConversationReport(messages) {
    // Ensure messages is an array
    if (!Array.isArray(messages)) {
        console.error('Messages is not an array:', typeof messages, messages);
        return `<html><body><h1>Error</h1><p>Invalid conversation data format</p></body></html>`;
    }
    
    let html = `
    <html>
    <head>
        <title>Annotation Boost Report</title>
        <style>
            body { font-family: sans-serif; }
            .message { margin-bottom: 1em; padding: 1em; border-radius: 5px; }
            .assistant { background-color: #e1f5fe; border-left: 5px solid #0288d1; }
            .user { background-color: #f1f8e9; border-left: 5px solid #689f38; }
            pre { white-space: pre-wrap; word-wrap: break-word; }
        </style>
    </head>
    <body>
        <h1>Annotation Boost Conversation</h1>
    `;

    for (const msg of messages) {
        if (!msg.content || !msg.content.trim()) continue;
        const className = msg.role === 'assistant' ? 'assistant' : 'user';
        const speaker = msg.role === 'assistant' ? 'Expert Analyst' : 'System';
        
        html += `<div class="message ${className}"><strong>${speaker}</strong><pre>${msg.content.trim()}</pre></div>`;
    }

    html += '</body></html>';
    return html;
}

/**
 * Convert tagged summary to professional HTML report
 */
function formatSummaryToHTML(summaryText, searchStrategy = "breadth", reportStyle = "per_iteration") {
    try {
        // Helper function to format hypotheses with better separation (matching Python implementation)
        function formatHypotheses(text) {
            if (!text || text === "No information available") {
                return text;
            }

            // Approach 1: Check if text contains numbered points using findall-style matching
            const numberedPointsRegex = /(?:^|\s)(\d+\.)\s+([^0-9\.].*?)(?=\s+\d+\.\s+|\s*$)/gs;
            const numberedPoints = [...text.matchAll(numberedPointsRegex)];

            if (numberedPoints.length > 0) {
                let result = '<div class="hypothesis-list">';
                for (const match of numberedPoints) {
                    result += `<div class="hypothesis-item"><span class="hypothesis-number">${match[1]}</span> ${match[2].trim()}</div>`;
                }
                result += '</div>';
                return result;
            }

            // Approach 2 (fallback): Look for numbers at beginning of paragraphs/lines
            const paragraphs = text.split('\n');
            const hasNumberedParagraphs = paragraphs.some(p => p.trim() && /^\s*\d+[\.\)\-]/.test(p));

            if (hasNumberedParagraphs) {
                let result = '<div class="hypothesis-list">';
                for (const para of paragraphs) {
                    if (!para.trim()) continue;

                    const lineMatch = para.match(/^\s*(\d+[\.\)\-])\s*(.*)/);
                    if (lineMatch) {
                        result += `<div class="hypothesis-item"><span class="hypothesis-number">${lineMatch[1]}</span> ${lineMatch[2].trim()}</div>`;
                    } else {
                        result += `<div class="hypothesis-item-continued">${para.trim()}</div>`;
                    }
                }
                result += '</div>';
                return result;
            }

            // Approach 3 (fallback): More aggressive split by numbered items
            const numberedItems = text.split(/(?:^|\s)(\d+[\.\)\-])(?=\s)/);
            if (numberedItems.length > 2) {
                let result = '<div class="hypothesis-list">';
                // Skip the first empty item
                for (let i = 1; i < numberedItems.length; i += 2) {
                    if (i + 1 < numberedItems.length) {
                        const num = numberedItems[i];
                        const content = numberedItems[i + 1].trim();
                        result += `<div class="hypothesis-item"><span class="hypothesis-number">${num}</span> ${content}</div>`;
                    }
                }
                result += '</div>';
                return result;
            }

            // If no patterns match, return original text with line breaks preserved
            return text.replace(/\n/g, '<br>');
        }

        // Helper function to convert markdown to HTML
        function markdownToHtml(text) {
            if (!text || text === "No information available") {
                return text;
            }

            return text
                // Bold: **text** or __text__
                .replace(/\*\*([^*]+)\*\*/g, '<strong>$1</strong>')
                .replace(/__([^_]+)__/g, '<strong>$1</strong>')
                // Italic: *text* or _text_ (but not inside words)
                .replace(/(?<![a-zA-Z])\*([^*]+)\*(?![a-zA-Z])/g, '<em>$1</em>')
                .replace(/(?<![a-zA-Z])_([^_]+)_(?![a-zA-Z])/g, '<em>$1</em>')
                // Line breaks
                .replace(/\n/g, '<br>');
        }

        // Extract sections using regex
        const sections = {};
        
        if (reportStyle.toLowerCase() === "total_summary") {
            // Gene-focused report sections
            const sectionPatterns = {
                'overview': /<OVERVIEW>\s*([\s\S]*?)\s*<\/OVERVIEW>/,
                'initial_hypothesis': /<INITIAL_HYPOTHESIS>\s*([\s\S]*?)\s*<\/INITIAL_HYPOTHESIS>/,
                'genes_analyzed': /<GENES_ANALYZED>\s*([\s\S]*?)\s*<\/GENES_ANALYZED>/,
                'final_conclusion': /<FINAL_CONCLUSION>\s*([\s\S]*?)\s*<\/FINAL_CONCLUSION>/,
                'key_insights': /<KEY_INSIGHTS>\s*([\s\S]*?)\s*<\/KEY_INSIGHTS>/,
                'validation_status': /<VALIDATION_STATUS>\s*([\s\S]*?)\s*<\/VALIDATION_STATUS>/,
            };
            
            for (const [key, pattern] of Object.entries(sectionPatterns)) {
                const match = summaryText.match(pattern);
                sections[key] = match ? match[1].trim() : "No information available";
            }
            
            // Extract gene groups
            const geneGroups = [];
            const genesAnalyzedContent = sections.genes_analyzed || '';
            const geneGroupMatches = genesAnalyzedContent.matchAll(/<GENE_GROUP_\d+>\s*([\s\S]*?)\s*<\/GENE_GROUP_\d+>/g);
            
            for (const match of geneGroupMatches) {
                const groupContent = match[1];
                const title = groupContent.match(/<TITLE>\s*([\s\S]*?)\s*<\/TITLE>/);
                const genes = groupContent.match(/<GENES>\s*([\s\S]*?)\s*<\/GENES>/);
                const findings = groupContent.match(/<FINDINGS>\s*([\s\S]*?)\s*<\/FINDINGS>/);
                
                geneGroups.push({
                    title: title ? title[1].trim() : "Gene Analysis",
                    genes: genes ? genes[1].trim() : "No genes listed",
                    findings: findings ? findings[1].trim() : "No findings available"
                });
            }
            
            // Generate HTML for gene-focused report
            return generateGeneFocusedHTML(sections, geneGroups, searchStrategy);
            
        } else {
            // Original per-iteration sections
            const sectionPatterns = {
                'overview': /<OVERVIEW>\s*([\s\S]*?)\s*<\/OVERVIEW>/,
                'initial_assessment': /<INITIAL_ASSESSMENT>\s*([\s\S]*?)\s*<\/INITIAL_ASSESSMENT>/,
                'final_annotation': /<FINAL_ANNOTATION>\s*([\s\S]*?)\s*<\/FINAL_ANNOTATION>/,
                'marker_summary': /<MARKER_SUMMARY>\s*([\s\S]*?)\s*<\/MARKER_SUMMARY>/,
            };

            for (const [key, pattern] of Object.entries(sectionPatterns)) {
                let match = summaryText.match(pattern);

                // Fallback: try case-insensitive if exact match fails
                if (!match) {
                    const caseInsensitivePattern = new RegExp(pattern.source, 'i');
                    match = summaryText.match(caseInsensitivePattern);
                    if (match) {
                        console.log(`Extracted ${key} using case-insensitive match`);
                    }
                }

                sections[key] = match ? match[1].trim() : "No information available";
            }

            // Warn if critical section missing
            if (sections.final_annotation === "No information available") {
                console.warn('⚠️ FINAL_ANNOTATION not found in LLM response');
            }
            
            // Extract iterations
            const iterations = [];
            const iterationMatches = summaryText.matchAll(/<ITERATION_(\d+)>\s*([\s\S]*?)\s*<\/ITERATION_\1>/g);
            
            for (const match of iterationMatches) {
                const iterNum = match[1];
                const iterContent = match[2];
                
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
            
            // Generate HTML for per-iteration report
            return generatePerIterationHTML(sections, iterations, searchStrategy, formatHypotheses, markdownToHtml);
        }
        
    } catch (error) {
        console.error('Error formatting summary to HTML:', error);
        return `<html><body><h1>Error</h1><p>Error formatting summary: ${error.message}</p></body></html>`;
    }
}

/**
 * Generate HTML for gene-focused report
 */
function generateGeneFocusedHTML(sections, geneGroups, searchStrategy) {
    const strategyDisplay = searchStrategy.charAt(0).toUpperCase() + searchStrategy.slice(1) + "-First Analysis";
    const now = new Date().toLocaleString();
    
    let geneGroupsHTML = '';
    for (const group of geneGroups) {
        geneGroupsHTML += `
            <div class="gene-group">
                <h4 class="gene-group-title">${group.title}</h4>
                <div class="gene-list"><strong>Genes:</strong> ${group.genes}</div>
                <div class="gene-findings">${group.findings}</div>
            </div>
        `;
    }
    
    return `
    <!DOCTYPE html>
    <html>
    <head>
        <title>CASSIA Annotation Boost - Gene-Focused Summary (${strategyDisplay})</title>
        <meta charset="UTF-8">
        <style>
            body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; background: #f8f9fa; }
            .container { max-width: 1000px; margin: 20px auto; background: white; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.1); overflow: hidden; }
            .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; text-align: center; }
            .header h1 { margin: 0 0 10px 0; font-size: 28px; font-weight: 600; }
            .header .subtitle { font-size: 16px; opacity: 0.9; }
            .content { padding: 30px; }
            .section { margin-bottom: 35px; }
            .section h2 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 8px; margin-bottom: 20px; font-size: 22px; }
            .section h3 { color: #34495e; margin-bottom: 15px; font-size: 18px; }
            .overview-box { background: #f8f9fa; border-left: 5px solid #3498db; padding: 20px; border-radius: 5px; white-space: pre-line; }
            .gene-group { background: #fff; border: 1px solid #e1e8ed; border-radius: 8px; padding: 20px; margin-bottom: 20px; }
            .gene-group-title { color: #2c3e50; margin: 0 0 15px 0; font-size: 16px; font-weight: 600; }
            .gene-list { background: #f1f8ff; padding: 12px; border-radius: 5px; margin-bottom: 12px; font-family: 'Courier New', monospace; }
            .gene-findings { line-height: 1.6; color: #2c3e50; white-space: pre-line; }
            .conclusion-box { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 25px; border-radius: 8px; white-space: pre-line; }
            .insights-box { background: #fff3cd; border: 1px solid #ffeaa7; padding: 20px; border-radius: 8px; white-space: pre-line; }
            .validation-box { background: #d1ecf1; border: 1px solid #bee5eb; padding: 20px; border-radius: 8px; white-space: pre-line; }
            .meta-info { text-align: center; color: #6c757d; font-size: 14px; padding: 20px; background: #f8f9fa; }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 CASSIA Annotation Boost Report</h1>
                <div class="subtitle">Gene-Focused Analysis Summary (${strategyDisplay})</div>
            </div>
            
            <div class="content">
                <div class="section">
                    <h2>📋 Overview</h2>
                    <div class="overview-box">${sections.overview}</div>
                </div>
                
                <div class="section">
                    <h2>🎯 Initial Hypothesis</h2>
                    <div>${sections.initial_hypothesis}</div>
                </div>
                
                <div class="section">
                    <h2>🧪 Genes Analyzed</h2>
                    ${geneGroupsHTML}
                </div>
                
                <div class="section">
                    <h2>✅ Final Conclusion</h2>
                    <div class="conclusion-box">${sections.final_conclusion}</div>
                </div>
                
                <div class="section">
                    <h2>💡 Key Insights</h2>
                    <div class="insights-box">${sections.key_insights}</div>
                </div>
                
                <div class="section">
                    <h2>🔬 Validation Status</h2>
                    <div class="validation-box">${sections.validation_status}</div>
                </div>
            </div>
            
            <div class="meta-info">
                Report generated on ${now}
            </div>
        </div>
    </body>
    </html>
    `;
}

/**
 * Generate HTML for per-iteration report
 */
function generatePerIterationHTML(sections, iterations, searchStrategy, formatHypotheses, markdownToHtml) {
    const strategyDisplay = searchStrategy.charAt(0).toUpperCase() + searchStrategy.slice(1) + "-First Analysis";
    const now = new Date().toLocaleString();

    let iterationsHTML = '';
    for (const iter of iterations) {
        iterationsHTML += `
            <div class="iteration-section">
                <h3>🔄 Iteration ${iter.number}</h3>
                <div class="iteration-content">
                    <div class="hypothesis-section">
                        <h4>Hypotheses Tested:</h4>
                        <div class="hypothesis-content">${formatHypotheses(iter.hypotheses)}</div>
                    </div>

                    <div class="genes-section">
                        <h4>Genes Checked:</h4>
                        <div class="genes-content">${markdownToHtml(iter.genes_checked)}</div>
                    </div>
                    
                    <div class="findings-section">
                        <h4>Key Findings:</h4>
                        <div class="findings-content">${markdownToHtml(iter.key_findings)}</div>
                    </div>
                </div>
            </div>
        `;
    }
    
    return `
    <!DOCTYPE html>
    <html>
    <head>
        <title>CASSIA Annotation Boost - Iteration Summary (${strategyDisplay})</title>
        <meta charset="UTF-8">
        <style>
            body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; background: #f8f9fa; }
            .container { max-width: 1000px; margin: 20px auto; background: white; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.1); overflow: hidden; }
            .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; text-align: center; }
            .header h1 { margin: 0 0 10px 0; font-size: 28px; font-weight: 600; }
            .header .subtitle { font-size: 16px; opacity: 0.9; }
            .content { padding: 30px; }
            .section { margin-bottom: 35px; }
            .section h2 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 8px; margin-bottom: 20px; font-size: 22px; }
            .overview-box { background: #f8f9fa; border-left: 5px solid #3498db; padding: 20px; border-radius: 5px; white-space: pre-line; }
            .iteration-section { background: #fff; border: 1px solid #e1e8ed; border-radius: 8px; padding: 25px; margin-bottom: 25px; }
            .iteration-section h3 { color: #2c3e50; margin: 0 0 20px 0; padding-bottom: 10px; border-bottom: 2px solid #ecf0f1; }
            .iteration-content h4 { color: #34495e; margin: 15px 0 8px 0; font-size: 16px; }
            .hypothesis-content, .genes-content, .findings-content { background: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px; white-space: pre-line; }
            .hypothesis-list { margin: 0; white-space: normal; }
            .hypothesis-item { margin-bottom: 12px; padding: 10px; background: white; border-radius: 5px; border-left: 4px solid #3498db; }
            .hypothesis-item-continued { padding-left: 1.5rem; margin-top: -0.3rem; color: #6c757d; }
            .hypothesis-number { font-weight: bold; color: #3498db; margin-right: 8px; }
            .genes-content { font-family: 'Courier New', monospace; background: #f1f8ff; }
            .conclusion-box, .marker-box, .recommendations-box { white-space: pre-line; }
            .conclusion-box { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 25px; border-radius: 8px; }
            .marker-box { background: #e8f5e9; border: 1px solid #c8e6c9; padding: 20px; border-radius: 8px; }
            .recommendations-box { background: #fff3e0; border: 1px solid #ffcc02; padding: 20px; border-radius: 8px; }
            .meta-info { text-align: center; color: #6c757d; font-size: 14px; padding: 20px; background: #f8f9fa; }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🧬 CASSIA Annotation Boost Report</h1>
                <div class="subtitle">Per-Iteration Analysis Summary (${strategyDisplay})</div>
            </div>
            
            <div class="content">
                <div class="section">
                    <h2>📋 Overview</h2>
                    <div class="overview-box">${markdownToHtml(sections.overview)}</div>
                </div>

                <div class="section">
                    <h2>🎯 Initial Assessment</h2>
                    <div>${markdownToHtml(sections.initial_assessment)}</div>
                </div>

                <div class="section">
                    <h2>🔄 Iterative Analysis</h2>
                    ${iterationsHTML}
                </div>

                <div class="section">
                    <h2>✅ Final Annotation</h2>
                    <div class="conclusion-box">${markdownToHtml(sections.final_annotation)}</div>
                </div>

                <div class="section">
                    <h2>🧬 Marker Summary</h2>
                    <div class="marker-box">${markdownToHtml(sections.marker_summary)}</div>
                </div>
            </div>
            
            <div class="meta-info">
                Report generated on ${now}
            </div>
        </div>
    </body>
    </html>
    `;
}