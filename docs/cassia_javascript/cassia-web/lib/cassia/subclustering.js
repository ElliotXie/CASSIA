/**
 * CASSIA Subclustering Agent - JavaScript Implementation for Browser
 * 100% Python-compatible subclustering functionality
 *
 * Replicates all functionality from cassia/subclustering.py:
 * - LLM-based subcluster annotation
 * - Prompt construction for subclusters
 * - Result extraction and CSV writing
 * - Single and batch processing
 * - Report generation integration (TODO)
 */

import { callLLM } from './llm_utils.js';

// ----------------- Helper Functions from runCASSIA_batch.js -----------------

/**
 * Process marker data to extract top genes per cluster (equivalent to dataUtils.getTopMarkers)
 */
function getTopMarkers(df, nGenes = 10, rankingMethod = "avg_log2FC", ascending = null, formatType = null) {
    // Auto-detect format if not specified
    const columns = Object.keys(df[0] || {});

    let isSeurat;
    if (formatType === "seurat") {
        isSeurat = true;
    } else if (formatType === "scanpy") {
        isSeurat = false;
        // TODO: Add Scanpy format support
        console.warn("Scanpy format not yet implemented, treating as pre-processed");
        return df;
    } else {
        // Auto-detect
        isSeurat = columns.includes('cluster') && columns.includes('gene') && columns.includes('avg_log2FC');
    }

    if (isSeurat) {
        return processSeurat(df, nGenes, rankingMethod, ascending);
    } else {
        // Assume pre-processed format with 2 columns
        return df;
    }
}

/**
 * Process Seurat-format differential expression data
 */
function processSeurat(df, nGenes, rankingMethod, ascending) {
    // Filter genes
    const filtered = df.filter(row => {
        const pVal = parseFloat(row.p_val_adj);
        const logFC = parseFloat(row.avg_log2FC);
        const pct1 = parseFloat(row['pct.1']);
        const pct2 = parseFloat(row['pct.2']);

        return !isNaN(pVal) && pVal < 0.05 &&
               !isNaN(logFC) && logFC > 0.25 &&
               (!isNaN(pct1) && pct1 >= 0.1 || !isNaN(pct2) && pct2 >= 0.1);
    });

    // Group by cluster
    const clusters = {};
    filtered.forEach(row => {
        const cluster = row.cluster;
        if (!clusters[cluster]) {
            clusters[cluster] = [];
        }
        clusters[cluster].push(row);
    });

    // Sort and get top genes for each cluster
    const result = [];
    for (const [cluster, genes] of Object.entries(clusters)) {
        // Sort by ranking method
        genes.sort((a, b) => {
            let aVal, bVal;
            switch (rankingMethod) {
                case "avg_log2FC":
                    aVal = parseFloat(a.avg_log2FC);
                    bVal = parseFloat(b.avg_log2FC);
                    return ascending === true ? aVal - bVal : bVal - aVal;
                case "p_val_adj":
                    aVal = parseFloat(a.p_val_adj);
                    bVal = parseFloat(b.p_val_adj);
                    return ascending === false ? bVal - aVal : aVal - bVal;
                case "pct_diff":
                    aVal = parseFloat(a['pct.1']) - parseFloat(a['pct.2']);
                    bVal = parseFloat(b['pct.1']) - parseFloat(b['pct.2']);
                    return ascending === true ? aVal - bVal : bVal - aVal;
                case "Score":
                    aVal = parseFloat(a.Score || a.score || 0);
                    bVal = parseFloat(b.Score || b.score || 0);
                    return ascending === true ? aVal - bVal : bVal - aVal;
                default:
                    aVal = parseFloat(a.avg_log2FC);
                    bVal = parseFloat(b.avg_log2FC);
                    return bVal - aVal;
            }
        });

        // Get top n genes
        const topGenes = genes.slice(0, nGenes).map(g => g.gene);
        result.push({
            cluster: cluster,
            markers: topGenes.join(', ')
        });
    }

    return result;
}

class Semaphore {
    constructor(max) {
        this.max = max;
        this.current = 0;
        this.queue = [];
    }

    async acquire() {
        if (this.current < this.max) {
            this.current++;
            return Promise.resolve();
        }

        return new Promise(resolve => {
            this.queue.push(resolve);
        });
    }

    release() {
        this.current--;
        if (this.queue.length > 0) {
            this.queue.shift()();
        }
    }
}


/**
 * Unified function to call LLM for subcluster annotation.
 * Replicates Python subcluster_agent_annotate_subcluster function.
 *
 * @param {string} userMessage - The prompt message for subcluster annotation
 * @param {string} apiKey - API key for the provider
 * @param {string|null} model - Model to use (defaults to provider's default if null)
 * @param {number} temperature - Temperature for generation (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @returns {Promise<string|Array>} The generated annotation as a string or structured data
 */
export async function subclusterAgentAnnotateSubcluster(userMessage, apiKey, model = null, temperature = 0, provider = "anthropic") {
    // Set default model based on provider if not specified
    if (model === null) {
        if (provider === "openai") {
            model = "gpt-4o";
        } else if (provider === "anthropic") {
            model = "claude-3-5-sonnet-20241022";
        } else if (provider === "openrouter") {
            model = "anthropic/claude-3.5-sonnet";
        } else if (provider.startsWith("http")) {
            // For custom API endpoints, use a default model if none specified
            model = model || "deepseek-chat";
            console.log(`Using model: ${model} with custom provider: ${provider}`);
        }
    }

    // Add JSON tags for providers that need them
    let modifiedMessage = userMessage;
    if (provider === "anthropic" || provider === "openrouter" || provider.startsWith("http")) {
        // Add JSON tags to help with structured output for these providers
        if (!modifiedMessage.toLowerCase().includes("<json>")) {
            modifiedMessage = `${modifiedMessage}\n\nPlease format your response as JSON:\n<json>\n{"response": "Your detailed analysis here"}\n</json>`;
        }
    }

    // Use the unified call_llm function
    const result = await callLLM(
        modifiedMessage,
        provider,
        model,
        apiKey,
        temperature,
        7000   // maxTokens
    );

    // Process the result for providers that use JSON tags
    if (provider === "anthropic" || provider === "openrouter" || provider.startsWith("http")) {
        // Extract content from JSON tags if present

        // First check for JSON format with ```json or <json> tags
        const jsonMatch = result.match(/(?:```json|<json>)(.*?)(?:```|<\/json>)/s);
        if (jsonMatch) {
            try {
                const jsonContent = jsonMatch[1].trim();
                // Try to parse as JSON
                const parsed = JSON.parse(jsonContent);
                if (typeof parsed === 'object' && parsed.response) {
                    // If it's an object with a response key, return the response
                    if (Array.isArray(parsed.response)) {
                        return parsed.response;
                    }
                    return parsed.response;
                }
                return parsed;
            } catch (jsonError) {
                // If JSON parsing fails, continue with other methods
            }
        }

        // Check if the entire response is a JSON string (common with DeepSeek API)
        if (result.trim().startsWith('{') && result.trim().endsWith('}')) {
            try {
                const parsed = JSON.parse(result);
                if (typeof parsed === 'object' && parsed.response) {
                    // If it's an object with a response key, return the response
                    if (Array.isArray(parsed.response)) {
                        return parsed.response;
                    }
                    return parsed.response;
                }
                return parsed;
            } catch (jsonError) {
                // If JSON parsing fails, return the original result
            }
        }
    }

    return result || '';
}

/**
 * Construct prompt from CSV for subcluster analysis.
 * Replicates Python construct_prompt_from_csv_subcluster function.
 *
 * @param {Array|Object} marker - DataFrame-like object or array containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {number} nGenes - Number of top genes to use (default: 50)
 * @returns {string} Generated prompt for subcluster analysis
 */
export function constructPromptFromCsvSubcluster(marker, majorClusterInfo, nGenes = 50) {
    let processedMarker;

    // Process data if it has more than 2 columns
    if (Array.isArray(marker) && marker.length > 0) {
        // Handle array of objects
        const firstRow = marker[0];
        const numColumns = Object.keys(firstRow).length;

        if (numColumns > 2) {
            console.log(`Processing input dataframe to get top ${nGenes} markers`);
            processedMarker = getTopMarkers(marker, nGenes);
        } else {
            console.log("Using input dataframe directly as it appears to be pre-processed (2 columns)");
            processedMarker = [...marker]; // Create a copy
        }
    } else if (marker && typeof marker === 'object' && marker.columns && marker.columns.length > 2) {
        // Handle DataFrame-like object
        console.log(`Processing input dataframe to get top ${nGenes} markers`);
        processedMarker = getTopMarkers(marker, nGenes);
    } else {
        console.log("Using input dataframe directly as it appears to be pre-processed (2 columns)");
        processedMarker = marker;
    }

    // Initialize the prompt with the major cluster information
    let prompt = `

You are an expert biologist specializing in cell type annotation, with deep expertise in immunology, cancer biology, and developmental biology.You will be given sets of highly expressed markers ranked by significance for some subclusters from the ${majorClusterInfo} cluster, identify what is the most likely top2 cell type each marker set implies.

Take a deep breath and work step by step. You'd better do a really good job or 1000 grandma are going to be in danger.
You will be tipped $10,000 if you do a good job.

For each output, provide:
1.Key marker:
2.Explanation:
3.Most likely top2 cell types:

Remember these subclusters are from a ${majorClusterInfo} big cluster. You must include all clusters mentioned in the analysis.
`;

    // Iterate over each row in the data
    if (Array.isArray(processedMarker)) {
        processedMarker.forEach((row, index) => {
            const clusterName = Object.values(row)[0];  // First column
            const markers = Object.values(row)[1];      // Second column
            prompt += `${index + 1}.${markers}\n`;
        });
    } else if (processedMarker && typeof processedMarker === 'object') {
        // Handle DataFrame-like structure
        const keys = Object.keys(processedMarker);
        const firstColKey = keys[0];
        const secondColKey = keys[1];

        for (let i = 0; i < processedMarker[firstColKey].length; i++) {
            const clusterName = processedMarker[firstColKey][i];
            const markers = processedMarker[secondColKey][i];
            prompt += `${i + 1}.${markers}\n`;
        }
    }

    return prompt;
}

/**
 * Annotate subclusters using an LLM.
 * Replicates Python annotate_subclusters function.
 *
 * @param {Array|Object} marker - DataFrame containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {string} apiKey - API key for the provider
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} nGenes - Number of top genes to use
 * @returns {Promise<string|Array>} The generated annotation as a string or structured data
 */
export async function annotateSubclusters(marker, majorClusterInfo, apiKey, model = "claude-3-5-sonnet-20241022", temperature = 0, provider = "anthropic", nGenes = 50) {
    const prompt = constructPromptFromCsvSubcluster(marker, majorClusterInfo, nGenes);
    const outputText = await subclusterAgentAnnotateSubcluster(prompt, apiKey, model, temperature, provider);
    return outputText;
}

/**
 * Extract multiple output results from subcluster analysis text.
 * Replicates Python extract_subcluster_results_with_llm_multiple_output function.
 *
 * @param {string|Array} analysisText - Text containing the analysis results or structured data
 * @param {string} apiKey - API key for the provider
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @returns {Promise<string|Array>} Extracted results in the format: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.
 */
export async function extractSubclusterResultsWithLlmMultipleOutput(analysisText, apiKey, provider = "anthropic", model = "claude-3-5-sonnet-20241022", temperature = 0) {
    // If the analysis_text is already in a structured format (array of objects),
    // we can return it directly for processing
    if (Array.isArray(analysisText) && analysisText.every(item => typeof item === 'object')) {
        console.log(`Analysis text is already in structured format (array of ${analysisText.length} dictionaries)`);
        return analysisText;
    }

    // Define the prompt to instruct the LLM
    const prompt = `You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

${analysisText}`;

    // Use the subcluster_agent_annotate function to get the extraction
    const result = await subclusterAgentAnnotateSubcluster(prompt, apiKey, model, temperature, provider);

    // Check if the result is already in structured format (e.g., from DeepSeek API)
    if (Array.isArray(result) && result.every(item => typeof item === 'object')) {
        return result;
    }

    return result;
}

/**
 * Extract results with reasons from subcluster analysis text.
 * Replicates Python extract_subcluster_results_with_llm function.
 *
 * @param {string|Array} analysisText - Text containing the analysis results or structured data
 * @param {string} apiKey - API key for the provider
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @returns {Promise<string|Array>} Extracted results in the format: results1(celltype1, celltype2, reason), results2(celltype1, celltype2, reason), etc.
 */
export async function extractSubclusterResultsWithLlm(analysisText, apiKey, provider = "anthropic", model = "claude-3-5-sonnet-20241022", temperature = 0) {
    // If the analysis_text is already in a structured format (array of objects),
    // we can return it directly for processing by write_results_to_csv
    if (Array.isArray(analysisText) && analysisText.every(item => typeof item === 'object')) {
        console.log(`Analysis text is already in structured format (array of ${analysisText.length} dictionaries)`);
        return analysisText;
    }

    // Define the prompt to instruct the LLM
    const prompt = `You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2,reason), results2(celltype1, celltype2,reason), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

${analysisText}`;

    // Use the subcluster_agent_annotate function to get the extraction
    const result = await subclusterAgentAnnotateSubcluster(prompt, apiKey, model, temperature, provider);

    // Check if the result is already in structured format (e.g., from DeepSeek API)
    if (Array.isArray(result) && result.every(item => typeof item === 'object')) {
        return result;
    }

    return result;
}


/**
 * Extract cell type results from LLM output and return as CSV string.
 *
 * @param {string|Array} results - LLM analysis results as string or structured data
 * @returns {Promise<Object>} Object containing csvContent and dataFrame
 */
export async function getResultsAsCsv(results) {
    // Handle different result formats
    if (Array.isArray(results)) {
        try {
            const rows = [];
            for (let i = 0; i < results.length; i++) {
                const item = results[i];
                if (typeof item === 'object') {
                    // Always check both capitalized and lowercase keys
                    const clusterId = String(item.cluster || (i + 1));

                    // Main type
                    let mainType = 'Unknown';
                    if (item.most_likely_top2_cell_types && Array.isArray(item.most_likely_top2_cell_types) && item.most_likely_top2_cell_types.length > 0) {
                        mainType = item.most_likely_top2_cell_types[0];
                    } else if (item.main_cell_type) {
                        mainType = item.main_cell_type;
                    } else if (item['Most likely top2 cell types'] && Array.isArray(item['Most likely top2 cell types']) && item['Most likely top2 cell types'].length > 0) {
                        mainType = item['Most likely top2 cell types'][0];
                    } else if (item['Main cell type']) {
                        mainType = item['Main cell type'];
                    }

                    // Sub type
                    let subType = 'Unknown';
                    if (item.most_likely_top2_cell_types && Array.isArray(item.most_likely_top2_cell_types) && item.most_likely_top2_cell_types.length > 1) {
                        subType = item.most_likely_top2_cell_types[1];
                    } else if (item.sub_cell_type) {
                        subType = item.sub_cell_type;
                    } else if (item['Most likely top2 cell types'] && Array.isArray(item['Most likely top2 cell types']) && item['Most likely top2 cell types'].length > 1) {
                        subType = item['Most likely top2 cell types'][1];
                    } else if (item['Sub cell type']) {
                        subType = item['Sub cell type'];
                    }

                    // Key markers
                    const keyMarkers = item.key_markers || item['Key marker'] || item['Key Marker'] || '';

                    // Reason/explanation
                    const reason = item.explanation || item.Explanation || item.reason || '';

                    rows.push([clusterId, mainType, subType, keyMarkers, reason]);
                }
            }

            if (rows.length > 0) {
                const columns = ['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'];
                const csvContent = [columns.join(','), ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))].join('\n');

                // Return DataFrame-like object
                const df = { columns };
                columns.forEach((col, idx) => {
                    df[col] = rows.map(row => row[idx]);
                });
                return { csvContent, dataFrame: df };
            }
        } catch (error) {
            console.log(`Error processing list results: ${error.message}`);
            console.log("Attempting to convert to string and process...");
            // Fall back to string processing
            results = String(results);
        }
    }

    // Process as string (original method)
    if (typeof results === 'string') {
        // Updated regex pattern to capture the reason
        const pattern = /results(\d+)\(([^,]+),\s*([^)]+)\)/g;
        const matches = [];
        let match;

        while ((match = pattern.exec(results)) !== null) {
            matches.push([match[1], match[2], match[3]]);
        }

        if (matches.length > 0) {
            // Convert matches to a DataFrame with the reason column
            const columns = ['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'];
            const csvContent = [columns.join(','), ...matches.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))].join('\n');

            // Return DataFrame-like object
            const df = { columns };
            columns.forEach((col, idx) => {
                df[col] = matches.map(row => row[idx]);
            });
            return { csvContent, dataFrame: df };
        } else {
            // Try alternative pattern without reason
            const altPattern = /results(\d+)\(([^,]+),\s*([^)]+)\)/g;
            const altMatches = [];
            let altMatch;

            while ((altMatch = altPattern.exec(results)) !== null) {
                altMatches.push([altMatch[1], altMatch[2], altMatch[3]]);
            }

            if (altMatches.length > 0) {
                // Convert matches to a DataFrame without reason column
                const columns = ['Result ID', 'main_cell_type', 'sub_cell_type'];
                const rows = altMatches.map(row => [...row, '']); // Add empty reason column
                const csvContent = [columns.concat(['reason']).join(','), ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))].join('\n');

                // Return DataFrame-like object
                const df = { columns: columns.concat(['reason']) };
                df.columns.forEach((col, idx) => {
                    df[col] = rows.map(row => row[idx]);
                });
                return { csvContent, dataFrame: df };
            }
        }
    }

    // If we get here, we couldn't process the results
    console.log(`Warning: Could not extract results in expected format. Results type: ${typeof results}`);
    const rawResults = String(results);

    // Create a minimal DataFrame to avoid errors
    const columns = ['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'];
    const rows = [['1', 'Unknown', 'Unknown', 'Could not parse results']];
    const csvContent = [columns.join(','), ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))].join('\n');

    const df = { columns };
    columns.forEach((col, idx) => {
        df[col] = rows.map(row => row[idx]);
    });
    return { csvContent, dataFrame: df, rawResults };
}


/**
 * Process subclusters from marker data and generate annotated results.
 * Replicates Python runCASSIA_subclusters function.
 *
 * @param {Array|Object} marker - DataFrame containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {string} apiKey - API key for the provider
 * @param {string} model - Model name to use
 * @param {number} temperature - Temperature parameter for API calls (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} nGenes - Number of top genes to use for analysis
 * @returns {Promise<Object>} Object with csvContent and dataFrame.
 */
export async function runCASSIASubclusters(marker, majorClusterInfo, apiKey, model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter", nGenes = 50) {
    // Construct prompt and get analysis from LLM
    const prompt = constructPromptFromCsvSubcluster(marker, majorClusterInfo, nGenes);
    const outputText = await subclusterAgentAnnotateSubcluster(prompt, apiKey, model, temperature, provider);

    // Extract structured results from the analysis text
    const results = await extractSubclusterResultsWithLlm(outputText, apiKey, provider, model, temperature);

    // Get results as CSV
    const { csvContent, dataFrame } = await getResultsAsCsv(results);

    // --- Generate HTML report for the single run CSV ---
    // TODO: Implement report generation in the browser
    // try {
    //     const { processEvaluationCsv } = await import('./generate_reports.js');
    //     // This part needs to be adapted for browser environment
    // } catch (error) {
    //     console.warn('Warning: Could not generate HTML report:', error.message);
    // }

    return { csvContent, dataFrame };
}

/**
 * Run multiple subcluster analyses in parallel and save results.
 * Replicates Python runCASSIA_n_subcluster function.
 *
 * @param {number} n - Number of analyses to run
 * @param {Array|Object} marker - DataFrame containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {string} apiKey - API key for the provider
 * @param {string} model - Model name to use
 * @param {number} temperature - Temperature parameter for API calls (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} maxWorkers - Maximum number of parallel workers
 * @param {number} nGenes - Number of top genes to use for analysis
 * @returns {Promise<Array>} Array of objects, each containing csvContent and dataFrame for a run.
 */
export async function runCASSIANSubcluster(n, marker, majorClusterInfo, apiKey, model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter", maxWorkers = 5, nGenes = 50) {
    // Force temperature to 0.3 for all runs
    temperature = 0.3;
    const semaphore = new Semaphore(maxWorkers);

    /**
     * Run a single analysis iteration
     * @param {number} i - Iteration index
     * @returns {Promise<Object>} Path to the generated CSV file
     */
    async function runSingleAnalysis(i) {
        await semaphore.acquire();
        try {
            console.log(`Starting iteration ${i + 1}...`);
            // Run the annotation process
            const outputText = await annotateSubclusters(marker, majorClusterInfo, apiKey, model, temperature, provider, nGenes);

            // Extract results
            const results = await extractSubclusterResultsWithLlmMultipleOutput(outputText, apiKey, provider, model, temperature);

            // Create DataFrame based on result type
            let rows = [];

            if (Array.isArray(results) && results.every(item => typeof item === 'object')) {
                // Handle structured data (array of objects)
                for (let idx = 0; idx < results.length; idx++) {
                    const item = results[idx];
                    if (typeof item === 'object') {
                        // Format 1: Dictionary with 'cluster' key and 'most_likely_top2_cell_types'
                        if (item.cluster && item.most_likely_top2_cell_types) {
                            const clusterId = String(item.cluster || (idx + 1));
                            const cellTypes = item.most_likely_top2_cell_types || ['Unknown', 'Unknown'];
                            const mainType = Array.isArray(cellTypes) && cellTypes.length > 0 ? cellTypes[0] : "Unknown";
                            const subType = Array.isArray(cellTypes) && cellTypes.length > 1 ? cellTypes[1] : "Unknown";
                            rows.push([clusterId, mainType, subType]);
                        }
                        // Format 2: Dictionary with key like 'results1'
                        else if (Object.keys(item).some(key => key.startsWith('result'))) {
                            const key = Object.keys(item).find(k => k.startsWith('result'));
                            const value = item[key];

                            if (Array.isArray(value) && value.length >= 2) {
                                // Format: {'results1': ['celltype1', 'celltype2']}
                                const mainType = value.length > 0 ? value[0] : "";
                                const subType = value.length > 1 ? value[1] : "";
                                const resultId = key.replace(/\D/g, '') || String(idx + 1);
                                rows.push([resultId, mainType, subType]);
                            } else if (typeof value === 'object' && value.celltype1 && value.celltype2) {
                                // Format: {'results1': {'celltype1': 'type1', 'celltype2': 'type2'}}
                                const mainType = value.celltype1 || "";
                                const subType = value.celltype2 || "";
                                const resultId = key.replace(/\D/g, '') || String(idx + 1);
                                rows.push([resultId, mainType, subType]);
                            }
                        }
                        // Format 3: Dictionary with 'key_markers', 'explanation', etc. (with underscores)
                        else if (item.key_markers || item.explanation) {
                            const clusterId = String(item.cluster || (idx + 1));
                            const cellTypes = item.most_likely_top2_cell_types || ['Unknown', 'Unknown'];
                            const mainType = Array.isArray(cellTypes) && cellTypes.length > 0 ? cellTypes[0] : "Unknown";
                            const subType = Array.isArray(cellTypes) && cellTypes.length > 1 ? cellTypes[1] : "Unknown";
                            rows.push([clusterId, mainType, subType]);
                        }
                        // Format 4: Dictionary with 'Key marker', 'Explanation', etc. (with spaces and capital letters)
                        else if (item['Key marker'] || item['Explanation'] || item['Most likely top2 cell types']) {
                            const clusterId = String(idx + 1);  // Use index as cluster ID since no cluster field
                            const cellTypes = item['Most likely top2 cell types'] || ['Unknown', 'Unknown'];
                            const mainType = Array.isArray(cellTypes) && cellTypes.length > 0 ? cellTypes[0] : "Unknown";
                            const subType = Array.isArray(cellTypes) && cellTypes.length > 1 ? cellTypes[1] : "Unknown";
                            rows.push([clusterId, mainType, subType]);
                        }
                    }
                }

                if (rows.length === 0) {
                    // Fallback if we couldn't extract structured data
                    rows = results.map((_, i) => [String(i + 1), "Unknown", "Unknown"]);
                }
            } else {
                // Use regex to extract the results (original method)
                const pattern = /results(\d+)\(([^,]+),\s*([^)]+)\)/g;
                const matches = [];
                let match;

                while ((match = pattern.exec(results)) !== null) {
                    matches.push([match[1], match[2], match[3]]);
                }

                if (matches.length > 0) {
                    rows = matches;
                } else {
                    // Fallback if regex didn't match
                    rows = [[String(i + 1), "Unknown", "Unknown"]];
                }
            }

            let markerDf;
            try {
                markerDf = getTopMarkers(marker, nGenes);
            } catch (error) {
                console.log(`Warning: ${error.message}. Using original marker dataframe.`);
                markerDf = Array.isArray(marker) ? [...marker] : marker;
            }

            const dataFrameRows = [];
            const minRows = Math.min(rows.length, markerDf.length);
            for (let j = 0; j < minRows; j++) {
                dataFrameRows.push({
                    'True Cell Type': Object.values(markerDf[j])[0],
                    'main_cell_type': rows[j][1],
                    'sub_cell_type': rows[j][2]
                });
            }

            const columns = ['True Cell Type', 'main_cell_type', 'sub_cell_type'];
            const csvContent = [
                columns.join(','),
                ...dataFrameRows.map(row => columns.map(col => `"${String(row[col]).replace(/"/g, '""')}"`).join(','))
            ].join('\n');
            
            console.log(`Finished iteration ${i + 1}.`);
            return { csvContent, dataFrame: dataFrameRows, iteration: i + 1 };
        } catch (error) {
            console.error(`Iteration ${i + 1} generated an exception: ${error.message}`);
            return { error: error.message, iteration: i + 1 };
        } finally {
            semaphore.release();
        }
    }

    const analysisPromises = [];
    for (let i = 0; i < n; i++) {
        analysisPromises.push(runSingleAnalysis(i));
    }

    const results = await Promise.all(analysisPromises);
    
    // TODO: implement batch report generation
    console.log("Batch analysis complete. Report generation is not yet implemented in the browser.");

    return results;
} 