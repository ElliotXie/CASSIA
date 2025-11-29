import fs from 'fs';
import path from 'path';
import { parse as parseCSV } from 'csv-parse/sync';
import { stringify as stringifyCSV } from 'csv-stringify/sync';
import { callLLM } from './llm_utils.js';

/**
 * Agent function that reads a CSV file with cell cluster annotations and merges/groups them.
 * 
 * @param {string} csvPath - Path to the CSV file containing cluster annotations
 * @param {string|null} outputPath - Path to save the results (if null, returns data without saving)
 * @param {string} provider - LLM provider to use ("openai", "anthropic", or "openrouter")
 * @param {string|null} model - Specific model to use (if null, uses default for provider)
 * @param {string|null} apiKey - API key for the provider (if null, gets from environment)
 * @param {string|null} additionalContext - Optional domain-specific context to help with annotation
 * @param {number} batchSize - Number of clusters to process in each LLM call (for efficiency)
 * @param {string} detailLevel - Level of detail for the groupings:
 *                               - "broad": More general cell categories (e.g., "Myeloid cells" for macrophages and dendritic cells)
 *                               - "detailed": More specific groupings that still consolidate very specific clusters
 *                               - "very_detailed": Most specific groupings with normalized and consistent naming
 * 
 * @returns {Promise<Array>} Array of objects with original annotations and suggested cell groupings
 */
export async function mergeAnnotations({
    csvPath,
    outputPath = null,
    provider = "openai",
    model = null,
    apiKey = null,
    additionalContext = null,
    batchSize = 20,
    detailLevel = "broad"
}) {
    // Validate detail_level parameter
    if (!["broad", "detailed", "very_detailed"].includes(detailLevel)) {
        throw new Error("detailLevel must be one of: 'broad', 'detailed', or 'very_detailed'");
    }
    
    // Set column name for results based on detail level
    const resultColumnMap = {
        "broad": "Merged_Grouping_1",
        "detailed": "Merged_Grouping_2",
        "very_detailed": "Merged_Grouping_3"
    };
    const resultColumn = resultColumnMap[detailLevel];
    
    // Read the CSV file
    let csvData;
    try {
        const csvContent = fs.readFileSync(csvPath, 'utf-8');
        csvData = parseCSV(csvContent, { 
            columns: true, 
            skip_empty_lines: true,
            cast: true
        });
        console.log(`Successfully read CSV file with ${csvData.length} rows.`);
    } catch (error) {
        throw new Error(`Error reading CSV file: ${error.message}`);
    }
    
    // Map expected column names to actual column names
    // Check for the expected column names based on the new information
    const columnMapping = {
        "cluster": "True Cell Type",  // "True Cell Type" is actually the cluster ID column
        "general_annotation": "Predicted Main Cell Type",
        "subtype_annotation": "Predicted Sub Cell Types"
    };
    
    // Verify that we found the necessary columns
    const availableColumns = Object.keys(csvData[0] || {});
    const missingColumns = Object.values(columnMapping).filter(col => !availableColumns.includes(col));
    if (missingColumns.length > 0) {
        throw new Error(`Required columns not found: ${missingColumns.join(', ')}. Available columns: ${availableColumns.join(', ')}`);
    }
    
    // Create a working copy of the data
    const workingData = csvData.map(row => ({ ...row }));
    
    // Extract first subtype if the subtype column contains comma-separated values
    const subtypeCol = columnMapping["subtype_annotation"];
    workingData.forEach(row => {
        const subtypeValue = row[subtypeCol];
        if (typeof subtypeValue === 'string' && subtypeValue.includes(',')) {
            row["processed_subtype"] = subtypeValue.split(',')[0].trim();
        } else {
            row["processed_subtype"] = subtypeValue;
        }
    });
    
    // Initialize the result column
    workingData.forEach(row => {
        row[resultColumn] = "";
    });
    
    // Process in batches for efficiency
    const totalRows = workingData.length;
    for (let i = 0; i < totalRows; i += batchSize) {
        const batchEnd = Math.min(i + batchSize, totalRows);
        const batch = workingData.slice(i, batchEnd);
        
        // Prepare prompt for LLM based on detail level
        const prompt = _createAnnotationPrompt(
            batch, 
            columnMapping["cluster"],
            columnMapping["general_annotation"],
            "processed_subtype",
            additionalContext,
            detailLevel
        );
        
        // Call LLM to get suggested groupings
        try {
            const response = await callLLM({
                prompt: prompt,
                provider: provider,
                model: model,
                apiKey: apiKey,
                temperature: 0.3,  // Lower temperature for more consistent results
                maxTokens: 2000,
                systemPrompt: "You are an expert cell biologist specializing in single-cell analysis. Your task is to analyze cluster annotations and suggest general cell groupings."
            });
            
            // Parse LLM response and update data
            const groupings = _parseLLMResponse(response, batch);
            groupings.forEach((grouping, idx) => {
                workingData[i + idx][resultColumn] = grouping;
            });
                
            console.log(`Processed clusters ${i+1}-${batchEnd} out of ${totalRows}`);
        } catch (error) {
            console.error(`Error processing batch ${i+1}-${batchEnd}: ${error.message}`);
        }
    }
    
    // Add the result column to the original data
    csvData.forEach((row, idx) => {
        row[resultColumn] = workingData[idx][resultColumn];
    });
    
    // Save results if output path is provided
    if (outputPath) {
        const csvOutput = stringifyCSV(csvData, { header: true });
        fs.writeFileSync(outputPath, csvOutput);
        console.log(`Results saved to ${outputPath}`);
    }
    
    return csvData;
}

/**
 * Create a prompt for the LLM to suggest groupings based on cluster annotations.
 * 
 * @param {Array} batch - Array of objects containing clusters to process
 * @param {string} clusterCol - Name of the cluster ID column
 * @param {string} generalCol - Name of the general annotation column
 * @param {string} subtypeCol - Name of the subtype annotation column
 * @param {string|null} additionalContext - Optional domain-specific context
 * @param {string} detailLevel - Level of detail for the groupings ("broad", "detailed", or "very_detailed")
 * 
 * @returns {string} Formatted prompt string
 */
function _createAnnotationPrompt(
    batch, 
    clusterCol,
    generalCol, 
    subtypeCol,
    additionalContext,
    detailLevel = "broad"
) {
    // Format clusters data
    const clustersText = batch.map(row => 
        `Cluster ${row[clusterCol]}: General annotation: ${row[generalCol]}, Subtype: ${row[subtypeCol]}`
    ).join('\n');
    
    let prompt;
    
    // Create the prompt based on detail level
    if (detailLevel === "broad") {
        // Broad groupings prompt (original)
        prompt = `I have single-cell RNA-seq cluster annotations and need to suggest broader cell groupings.
For each cluster, I'll provide the general annotation and subtype annotation.
Based on these annotations, suggest an appropriate broader cell grouping category.

For example:
- "macrophage, inflammatory macrophage" → "Myeloid cells"
- "CD4 T cell, naive CD4 T cell" → "T cells"
- "B cell, memory B cell" → "B cells"

Use general cell lineage categories when possible, combining related cell types into a single group.
Prioritize creating broader categories that span multiple specific cell types.

Annotations to process:
${clustersText}

Please respond with a JSON object where keys are cluster identifiers and values are the suggested groupings. 
For example:
{
  "1": "Myeloid cells",
  "2": "T cells"
}`;
    } else if (detailLevel === "detailed") {
        // Detailed groupings prompt
        prompt = `I have single-cell RNA-seq cluster annotations and need to suggest intermediate-level cell groupings.
For each cluster, I'll provide the general annotation and subtype annotation.
Based on these annotations, suggest a moderately specific cell grouping that balances detail and generality.

For example:
- "macrophage, inflammatory macrophage" → "Macrophages" (not as broad as "Myeloid cells")
- "CD4 T cell, naive CD4 T cell" → "CD4 T cells" (more specific than just "T cells")
- "CD8 T cell, cytotoxic CD8 T cell" → "CD8 T cells" (more specific than just "T cells")
- "B cell, memory B cell" → "B cells" (specific cell type)

Maintain biological specificity when important, but still group very similar subtypes together.
Aim for a middle ground - not too general, but also not too specific.
The grouping should be more detailed than broad categories like "Myeloid cells" or "Lymphoid cells", 
but should still consolidate highly specific annotations.

Annotations to process:
${clustersText}

Please respond with a JSON object where keys are cluster identifiers and values are the suggested groupings. 
For example:
{
  "1": "Macrophages",
  "2": "CD4 T cells",
  "3": "CD8 T cells"
}`;
    } else {  // very_detailed
        // Very detailed groupings prompt
        prompt = `I have single-cell RNA-seq cluster annotations and need to normalize and standardize cell type names 
while preserving the most specific and detailed biological information.
For each cluster, I'll provide the general annotation and subtype annotation.

Your task is to create a consistent and standardized cell type label that:
1. Maintains the highest level of biological specificity from the annotations
2. Uses consistent nomenclature across similar cell types
3. Follows standard cell type naming conventions
4. Preserves functional or activation state information when present
5. Normalizes naming variations (e.g., "inflammatory macrophage" vs "M1 macrophage" should use one consistent term)

Examples:
- "macrophage, inflammatory macrophage" → "Inflammatory macrophages" (preserve activation state)
- "CD4 T cell, naive CD4 T cell" → "Naive CD4+ T cells" (preserve naive state, standardize CD4+)
- "CD8 T cell, cytotoxic CD8 T cell" → "Cytotoxic CD8+ T cells" (preserve function, standardize CD8+)
- "dendritic cell, plasmacytoid dendritic cell" → "Plasmacytoid dendritic cells" (preserve specific subtype)
- "B cell, memory B cell" → "Memory B cells" (preserve memory state)
- "NK cell, CD56bright NK cell" → "CD56bright NK cells" (preserve specific marker)

Annotations to process:
${clustersText}

Please respond with a JSON object where keys are cluster identifiers and values are the normalized, specific cell type labels.
For example:
{
  "1": "Inflammatory macrophages",
  "2": "Naive CD4+ T cells",
  "3": "Memory B cells"
}`;
    }
    
    // Add additional context if provided
    if (additionalContext) {
        prompt += `\n\nAdditional context that may help with the analysis:\n${additionalContext}`;
    }
    
    return prompt;
}

/**
 * Parse the LLM response to extract suggested groupings.
 * 
 * @param {string} response - LLM response text
 * @param {Array} batch - Array of objects for the batch
 * 
 * @returns {Array<string>} Array of suggested groupings corresponding to batch order
 */
function _parseLLMResponse(response, batch) {
    const groupings = new Array(batch.length).fill("Error parsing response");
    
    // Try to find and parse JSON in the response
    try {
        // Extract JSON if it's embedded in text
        const jsonMatch = response.match(/({[\s\S]*})/);
        if (jsonMatch) {
            const jsonStr = jsonMatch[1];
            const parsed = JSON.parse(jsonStr);
            
            // Map the parsed results to batch indices
            let i = 0;
            for (const [clusterId, grouping] of Object.entries(parsed)) {
                if (i < batch.length) {
                    groupings[i] = grouping;
                    i++;
                }
            }
        } else {
            // Fallback: Try parsing line by line
            const lines = response.split('\n').filter(line => line.trim());
            for (let i = 0; i < Math.min(lines.length, batch.length); i++) {
                const line = lines[i];
                if (line.includes(':')) {
                    groupings[i] = line.split(':', 2)[1].trim();
                }
            }
        }
    } catch (error) {
        console.error(`Error parsing LLM response: ${error.message}`);
        // Fallback: Use the raw response
        for (let i = 0; i < batch.length; i++) {
            groupings[i] = "Error parsing response";
        }
    }
    
    return groupings;
}

/**
 * Process all three detail levels in parallel and return a combined data array.
 * 
 * @param {string} csvPath - Path to the CSV file containing cluster annotations
 * @param {string|null} outputPath - Path to save the results (if null, returns data without saving)
 * @param {string} provider - LLM provider to use ("openai", "anthropic", or "openrouter")
 * @param {string|null} model - Specific model to use (if null, uses default for provider)
 * @param {string|null} apiKey - API key for the provider (if null, gets from environment)
 * @param {string|null} additionalContext - Optional domain-specific context to help with annotation
 * @param {number} batchSize - Number of clusters to process in each LLM call (for efficiency)
 * 
 * @returns {Promise<Array>} Array of objects with original annotations and all three levels of suggested cell groupings
 */
export async function mergeAnnotationsAll({
    csvPath,
    outputPath = null,
    provider = "openai",
    model = null,
    apiKey = null,
    additionalContext = null,
    batchSize = 20
}) {
    console.log("Processing all three detail levels in parallel...");
    
    // Define the detail levels to process
    const detailLevels = ["broad", "detailed", "very_detailed"];
    
    // Run all three detail levels in parallel
    const promises = detailLevels.map(detailLevel =>
        mergeAnnotations({
            csvPath,
            outputPath: null,  // We'll save the combined result at the end
            provider,
            model,
            apiKey,
            additionalContext,
            batchSize,
            detailLevel
        }).then(result => ({ detailLevel, result }))
    );
    
    try {
        const results = await Promise.all(promises);
        console.log("All detail levels completed successfully");
        
        // Combine results
        let combinedData = null;
        const resultColumnMap = {
            "broad": "Merged_Grouping_1",
            "detailed": "Merged_Grouping_2",
            "very_detailed": "Merged_Grouping_3"
        };
        
        for (const { detailLevel, result } of results) {
            const resultColumn = resultColumnMap[detailLevel];
            
            if (combinedData === null) {
                combinedData = result.map(row => ({ ...row }));
            } else {
                // Add this level's results to the combined data
                result.forEach((row, idx) => {
                    if (idx < combinedData.length) {
                        combinedData[idx][resultColumn] = row[resultColumn];
                    }
                });
            }
            
            console.log(`Completed processing for ${detailLevel} detail level`);
        }
        
        // Save combined results if output path is provided
        if (outputPath && combinedData) {
            const csvOutput = stringifyCSV(combinedData, { header: true });
            fs.writeFileSync(outputPath, csvOutput);
            console.log(`Combined results saved to ${outputPath}`);
        }
        
        return combinedData;
    } catch (error) {
        throw new Error(`Parallel processing failed: ${error.message}`);
    }
}