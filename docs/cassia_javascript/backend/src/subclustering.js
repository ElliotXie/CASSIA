/**
 * CASSIA Subclustering Agent - JavaScript Implementation
 * 100% Python-compatible subclustering functionality
 * 
 * Replicates all functionality from cassia/subclustering.py:
 * - LLM-based subcluster annotation
 * - Prompt construction for subclusters  
 * - Result extraction and CSV writing
 * - Single and batch processing
 * - Report generation integration
 */

import fs from 'fs';
import path from 'path';
import { Worker } from 'worker_threads';

// Import utility functions
try {
    // Handle both package and direct imports
    const { callLLM } = await import('./llm_utils.js');
    global.callLLM = callLLM;
} catch (error) {
    console.warn('Warning: Could not import llm_utils.js:', error.message);
}

try {
    const { getTopMarkers } = await import('./utils/dataUtils.js'); 
    global.getTopMarkers = getTopMarkers;
} catch (error) {
    console.warn('Warning: Could not import dataUtils.js:', error.message);
}

/**
 * Unified function to call LLM for subcluster annotation.
 * Replicates Python subcluster_agent_annotate_subcluster function.
 * 
 * @param {string} userMessage - The prompt message for subcluster annotation
 * @param {string|null} model - Model to use (defaults to provider's default if null)
 * @param {number} temperature - Temperature for generation (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @returns {Promise<string|Array>} The generated annotation as a string or structured data
 */
export async function subclusterAgentAnnotateSubcluster(userMessage, model = null, temperature = 0, provider = "anthropic") {
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
    const result = await global.callLLM(
        modifiedMessage,
        provider,
        model,
        null,  // apiKey will be taken from environment
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
            processedMarker = global.getTopMarkers ? global.getTopMarkers(marker, nGenes) : marker;
        } else {
            console.log("Using input dataframe directly as it appears to be pre-processed (2 columns)");
            processedMarker = [...marker]; // Create a copy
        }
    } else if (marker && typeof marker === 'object' && marker.columns && marker.columns.length > 2) {
        // Handle DataFrame-like object
        console.log(`Processing input dataframe to get top ${nGenes} markers`);
        processedMarker = global.getTopMarkers ? global.getTopMarkers(marker, nGenes) : marker;
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
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} nGenes - Number of top genes to use
 * @returns {Promise<string|Array>} The generated annotation as a string or structured data
 */
export async function annotateSubclusters(marker, majorClusterInfo, model = "claude-3-5-sonnet-20241022", temperature = 0, provider = "anthropic", nGenes = 50) {
    const prompt = constructPromptFromCsvSubcluster(marker, majorClusterInfo, nGenes);
    const outputText = await subclusterAgentAnnotateSubcluster(prompt, model, temperature, provider);
    return outputText;
}

/**
 * Extract multiple output results from subcluster analysis text.
 * Replicates Python extract_subcluster_results_with_llm_multiple_output function.
 * 
 * @param {string|Array} analysisText - Text containing the analysis results or structured data
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @returns {Promise<string|Array>} Extracted results in the format: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.
 */
export async function extractSubclusterResultsWithLlmMultipleOutput(analysisText, provider = "anthropic", model = "claude-3-5-sonnet-20241022", temperature = 0) {
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
    const result = await subclusterAgentAnnotateSubcluster(prompt, model, temperature, provider);

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
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)  
 * @param {string} model - Model to use (defaults to Claude 3.5 Sonnet)
 * @param {number} temperature - Temperature for generation (0-1)
 * @returns {Promise<string|Array>} Extracted results in the format: results1(celltype1, celltype2, reason), results2(celltype1, celltype2, reason), etc.
 */
export async function extractSubclusterResultsWithLlm(analysisText, provider = "anthropic", model = "claude-3-5-sonnet-20241022", temperature = 0) {
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
    const result = await subclusterAgentAnnotateSubcluster(prompt, model, temperature, provider);
    
    // Check if the result is already in structured format (e.g., from DeepSeek API)
    if (Array.isArray(result) && result.every(item => typeof item === 'object')) {
        return result;
    }
    
    return result;
}

/**
 * Extract cell type results from LLM output and write to CSV file.
 * Replicates Python write_results_to_csv function.
 * 
 * @param {string|Array} results - LLM analysis results as string or structured data
 * @param {string} outputName - Base name for output file (will add .csv if not present)
 * @returns {Promise<Object>} DataFrame-like object containing the extracted results
 */
export async function writeResultsToCsv(results, outputName = 'subcluster_results') {
    // Add .csv suffix if not present
    if (!outputName.toLowerCase().endsWith('.csv')) {
        outputName = outputName + '.csv';
    }
    
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
                
                fs.writeFileSync(outputName, csvContent);
                console.log(`Results have been written to ${outputName}`);
                
                // Return DataFrame-like object
                const df = { columns };
                columns.forEach((col, idx) => {
                    df[col] = rows.map(row => row[idx]);
                });
                return df;
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
            
            fs.writeFileSync(outputName, csvContent);
            console.log(`Results have been written to ${outputName}`);
            
            // Return DataFrame-like object
            const df = { columns };
            columns.forEach((col, idx) => {
                df[col] = matches.map(row => row[idx]);
            });
            return df;
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
                
                fs.writeFileSync(outputName, csvContent);
                console.log(`Results have been written to ${outputName} (without reasons)`);
                
                // Return DataFrame-like object
                const df = { columns: columns.concat(['reason']) };
                df.columns.forEach((col, idx) => {
                    df[col] = rows.map(row => row[idx]);
                });
                return df;
            }
        }
    }
    
    // If we get here, we couldn't process the results
    console.log(`Warning: Could not extract results in expected format. Results type: ${typeof results}`);
    console.log(`Saving raw results to ${outputName}.txt for inspection`);
    
    // Save raw results for debugging
    fs.writeFileSync(`${outputName}.txt`, String(results));
    
    // Create a minimal DataFrame to avoid errors
    const columns = ['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'];
    const rows = [['1', 'Unknown', 'Unknown', 'Could not parse results']];
    const csvContent = [columns.join(','), ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))].join('\n');
    
    fs.writeFileSync(outputName, csvContent);
    
    // Return DataFrame-like object
    const df = { columns };
    columns.forEach((col, idx) => {
        df[col] = rows.map(row => row[idx]);
    });
    return df;
}

/**
 * Process subclusters from marker data and generate annotated results.
 * Replicates Python runCASSIA_subclusters function.
 * 
 * @param {Array|Object} marker - DataFrame containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {string} outputName - Base name for output file (will add .csv if not present)
 * @param {string} model - Model name to use
 * @param {number} temperature - Temperature parameter for API calls (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} nGenes - Number of top genes to use for analysis
 * @returns {Promise<null>} Results are saved to a CSV file
 */
export async function runCASSIASubclusters(marker, majorClusterInfo, outputName, model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter", nGenes = 50) {
    // Construct prompt and get analysis from LLM
    const prompt = constructPromptFromCsvSubcluster(marker, majorClusterInfo, nGenes);
    const outputText = await subclusterAgentAnnotateSubcluster(prompt, model, temperature, provider);

    // Extract structured results from the analysis text
    const results = await extractSubclusterResultsWithLlm(outputText, provider, model, temperature);
    
    // Save results to CSV
    await writeResultsToCsv(results, outputName);

    // --- Generate HTML report for the single run CSV ---
    try {
        const { processEvaluationCsv } = await import('./generate_reports.js');
        const csvFile = outputName.toLowerCase().endsWith('.csv') ? outputName : outputName + '.csv';
        if (fs.existsSync(csvFile)) {
            await processEvaluationCsv(csvFile, true);
        }
    } catch (error) {
        console.warn('Warning: Could not generate HTML report:', error.message);
    }
    
    return null;
}

/**
 * Run multiple subcluster analyses in parallel and save results.
 * Replicates Python runCASSIA_n_subcluster function.
 * 
 * @param {number} n - Number of analyses to run
 * @param {Array|Object} marker - DataFrame containing marker data
 * @param {string} majorClusterInfo - Description of the major cluster type
 * @param {string} baseOutputName - Base name for output files
 * @param {string} model - Model name to use
 * @param {number} temperature - Temperature parameter for API calls (0-1)
 * @param {string} provider - LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
 * @param {number} maxWorkers - Maximum number of parallel workers
 * @param {number} nGenes - Number of top genes to use for analysis
 * @returns {Promise<null>} Results are saved to CSV files
 */
export async function runCASSIANSubcluster(n, marker, majorClusterInfo, baseOutputName, model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter", maxWorkers = 5, nGenes = 50) {
    // Force temperature to 0.3 for all runs
    temperature = 0.3;
    
    /**
     * Run a single analysis iteration
     * @param {number} i - Iteration index
     * @returns {Promise<string>} Path to the generated CSV file
     */
    async function runSingleAnalysis(i) {
        // Run the annotation process
        const outputText = await annotateSubclusters(marker, majorClusterInfo, model, temperature, provider, nGenes);
        
        // Extract results
        const results = await extractSubclusterResultsWithLlmMultipleOutput(outputText, provider, model, temperature);
        
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

        try {
            // Try to get top markers, but handle the case where required columns are missing
            let markerDf;
            try {
                markerDf = global.getTopMarkers ? global.getTopMarkers(marker, nGenes) : marker;
            } catch (error) {
                // If get_top_markers fails due to missing columns, use the original marker dataframe
                console.log(`Warning: ${error.message}. Using original marker dataframe.`);
                markerDf = Array.isArray(marker) ? [...marker] : marker;
            }
            
            // Convert DataFrame to format we can work with
            const df = {
                'True Cell Type': rows.map(row => String(row[0])),
                'main_cell_type': rows.map(row => row[1]),
                'sub_cell_type': rows.map(row => row[2])
            };
            
            // Make a copy of the original values before swapping
            const originalTrueCellTypes = [...df['True Cell Type']];
            
            // Check if marker_df has data and perform the swap safely
            if (Array.isArray(markerDf) && markerDf.length > 0) {
                const minRows = Math.min(df['True Cell Type'].length, markerDf.length);
                if (minRows > 0) {
                    for (let j = 0; j < minRows; j++) {
                        df['True Cell Type'][j] = Object.values(markerDf[j])[0]; // First column of marker data
                    }
                }
            } else if (markerDf && typeof markerDf === 'object' && Object.keys(markerDf).length > 0) {
                const firstColKey = Object.keys(markerDf)[0];
                const firstColData = markerDf[firstColKey];
                if (Array.isArray(firstColData)) {
                    const minRows = Math.min(df['True Cell Type'].length, firstColData.length);
                    for (let j = 0; j < minRows; j++) {
                        df['True Cell Type'][j] = firstColData[j];
                    }
                }
            } else {
                console.log("Warning: Marker dataframe is empty or has no columns. Skipping column swap.");
            }
        } catch (error) {
            console.log(`Warning: Error during column swap: ${error.message}. Continuing without swapping.`);
        }

        // Write the DataFrame to a CSV file with an index
        const indexedCsvFilePath = `${baseOutputName}_${i + 1}.csv`;
        const csvContent = [
            'True Cell Type,main_cell_type,sub_cell_type',
            ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))
        ].join('\n');
        
        fs.writeFileSync(indexedCsvFilePath, csvContent);
        
        return indexedCsvFilePath;
    }

    // Run analyses with limited concurrency
    const resultFiles = [];
    const runningPromises = [];
    
    for (let i = 0; i < n; i++) {
        const promise = runSingleAnalysis(i).then(resultFile => {
            console.log(`Results for iteration ${i + 1} have been written to ${resultFile}`);
            resultFiles.push(resultFile);
            return resultFile;
        }).catch(error => {
            console.log(`Iteration ${i + 1} generated an exception: ${error.message}`);
            throw error;
        });
        
        runningPromises.push(promise);
        
        // Limit concurrency
        if (runningPromises.length >= maxWorkers) {
            await Promise.race(runningPromises);
            // Remove completed promises
            for (let j = runningPromises.length - 1; j >= 0; j--) {
                const p = runningPromises[j];
                const isResolved = await Promise.race([p.then(() => true, () => true), Promise.resolve(false)]);
                if (isResolved) {
                    runningPromises.splice(j, 1);
                }
            }
        }
    }
    
    // Wait for all remaining promises
    await Promise.allSettled(runningPromises);

    // --- Generate HTML reports for all batch CSVs ---
    try {
        const { processEvaluationCsv, createIndexHtml } = await import('./generate_reports.js');

        // Generate HTML report for each CSV
        for (const csvFile of resultFiles) {
            if (fs.existsSync(csvFile)) {
                await processEvaluationCsv(csvFile, true);
            }
        }

        // Create an index.html summary in the same directory as the first result file
        if (resultFiles.length > 0) {
            const outputDir = path.dirname(resultFiles[0]) || '.';
            await createIndexHtml(resultFiles, outputDir);
            console.log(`Batch HTML reports and index generated in ${outputDir}`);
        }
    } catch (error) {
        console.warn('Warning: Could not generate HTML reports:', error.message);
    }
}

/**
 * Test function to simulate a response from a custom API provider and test the parsing functionality.
 * This is useful for debugging the parsing logic without making actual API calls.
 * Replicates Python test_custom_api_parsing function.
 */
export async function testCustomApiParsing() {
    // Sample structured response that mimics what DeepSeek or other custom APIs might return
    const sampleResponse = [
        {
            'cluster': 1,
            'key_markers': 'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'explanation': 'The presence of IL7R, CD8A, and CD8B suggests a CD8+ T cell identity. CCL4 is associated with effector functions, while KLRB1 (CD161) and ITK indicate a memory-like or tissue-resident phenotype.',
            'most_likely_top2_cell_types': ['CD8+ memory T cells', 'Tissue-resident memory CD8+ T cells (TRM)']
        },
        {
            'cluster': 2,
            'key_markers': 'LAYN, HAVCR2 (TIM-3), TIGIT, IKZF2, KLRC2, KLRC3',
            'explanation': 'LAYN (Lag-3) and HAVCR2 (TIM-3) are markers of exhausted or chronically stimulated CD8+ T cells.',
            'most_likely_top2_cell_types': ['Exhausted CD8+ T cells', 'NK-like CD8+ T cells']
        },
        {
            'cluster': 3,
            'key_markers': 'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'explanation': 'GZMK, GZMH, PRF1, and NKG7 are markers of cytotoxic activity, typical of effector CD8+ T cells.',
            'most_likely_top2_cell_types': ['Effector CD8+ T cells', 'Central memory CD8+ T cells']
        },
        {
            'cluster': 4,
            'key_markers': 'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'explanation': 'WFDC2, CEACAM7, and CLDN8 are markers associated with epithelial or secretory cells.',
            'most_likely_top2_cell_types': ['Regulatory CD8+ T cells', 'Epithelial-like CD8+ T cells (rare subset)']
        }
    ];
    
    // Test the write_results_to_csv function
    console.log("Testing writeResultsToCsv with structured data...");
    const df = await writeResultsToCsv(sampleResponse, 'test_parsing_result');
    console.log(`Generated DataFrame:`, df);
    
    // Test the extract_subcluster_results_with_llm function
    console.log("\nTesting extractSubclusterResultsWithLlm with structured data...");
    const result = await extractSubclusterResultsWithLlm(sampleResponse);
    console.log(`Result type: ${typeof result}`);
    
    // Test the extract_subcluster_results_with_llm_multiple_output function
    console.log("\nTesting extractSubclusterResultsWithLlmMultipleOutput with structured data...");
    const result2 = await extractSubclusterResultsWithLlmMultipleOutput(sampleResponse);
    console.log(`Result type: ${typeof result2}`);
    
    // Create a simple marker dataframe for testing
    console.log("\nCreating test marker dataframe...");
    const testMarkerDf = [
        { cluster: 'cluster1', markers: 'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK' },
        { cluster: 'cluster2', markers: 'LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3' },
        { cluster: 'cluster3', markers: 'GZMK, GZMH, PRF1, NKG7, CCR7, CD27' },
        { cluster: 'cluster4', markers: 'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13' },
        { cluster: 'cluster5', markers: 'GNLY, KLRF1, FCER1G, TYROBP, CD38, KIR2DL4' },
        { cluster: 'cluster6', markers: 'LPL, SNAI2, HAND2, SOX2, NES, PDGFRA' }
    ];
    console.log(`Created test marker dataframe with shape: ${testMarkerDf.length} rows, ${Object.keys(testMarkerDf[0]).length} columns`);
    
    // Save the test dataframe to a temporary file
    const tempFile = 'test_markers_temp.csv';
    const csvContent = [
        'cluster,markers',
        ...testMarkerDf.map(row => `"${row.cluster}","${row.markers}"`)
    ].join('\n');
    fs.writeFileSync(tempFile, csvContent);
    console.log(`Saved test marker dataframe to ${tempFile}`);
    
    try {
        // Create a mock function to simulate annotate_subclusters
        const originalAnnotate = global.annotateSubclusters;
        
        // Replace with mock function  
        global.annotateSubclusters = async () => sampleResponse;
        
        console.log("\nTesting runCASSIANSubcluster with mock data...");
        // Run with n=1 to test a single iteration
        await runCASSIANSubcluster(
            1,
            testMarkerDf,
            "cd8 t cell",
            "test_n_subcluster",
            "dummy-model",
            0,
            "dummy-provider"
        );
        
        // Check if the output file was created
        const outputFile = "test_n_subcluster_1.csv";
        if (fs.existsSync(outputFile)) {
            console.log(`Successfully created output file: ${outputFile}`);
            const resultContent = fs.readFileSync(outputFile, 'utf-8');
            console.log(`Output file contents:\n${resultContent}`);
        } else {
            console.log(`Error: Output file ${outputFile} was not created`);
        }
        
        // Restore the original function
        if (originalAnnotate) {
            global.annotateSubclusters = originalAnnotate;
        }
        
        console.log("\nAll tests completed.");
        
    } catch (error) {
        console.log(`Error during testing: ${error.message}`);
        console.error(error.stack);
    } finally {
        // Clean up temporary files
        if (fs.existsSync(tempFile)) {
            fs.unlinkSync(tempFile);
            console.log(`Removed temporary file: ${tempFile}`);
        }
    }
}

// Note: All functions are already exported inline with their definitions