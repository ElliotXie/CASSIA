import { callLLM } from './llm_utils.js';

/**
 * Browser-compatible version of the annotation merging agent.
 * Processes CSV data (array of objects) and merges/groups cell cluster annotations.
 */

/**
 * Merge cell cluster annotations using LLM analysis.
 * 
 * @param {Array} csvData - Array of objects containing cluster annotations
 * @param {string} provider - LLM provider to use ("openai", "anthropic", or "openrouter")
 * @param {string|null} model - Specific model to use (if null, uses default for provider)
 * @param {string} apiKey - API key for the provider
 * @param {string|null} additionalContext - Optional domain-specific context
 * @param {number} batchSize - Number of clusters to process in each LLM call
 * @param {string} detailLevel - Level of detail: "broad", "detailed", or "very_detailed"
 * @param {function} onProgress - Optional callback for progress updates
 * @returns {Promise<Array>} Array of objects with original annotations and suggested groupings
 */
export async function mergeAnnotations({
    csvData,
    provider = "openrouter",
    model = null,
    apiKey,
    additionalContext = null,
    batchSize = 20,
    detailLevel = "broad",
    onProgress = null
}) {
    // Validate parameters
    if (!csvData || !Array.isArray(csvData) || csvData.length === 0) {
        throw new Error("csvData must be a non-empty array of objects");
    }

    if (!apiKey) {
        throw new Error("API key is required");
    }

    if (!["broad", "detailed", "very_detailed"].includes(detailLevel)) {
        throw new Error("detailLevel must be one of: 'broad', 'detailed', or 'very_detailed'");
    }

    if (onProgress) onProgress(`Starting annotation merging with ${detailLevel} detail level...`);

    // Map result column names based on detail level
    const resultColumnMap = {
        "broad": "Merged_Grouping_1",
        "detailed": "Merged_Grouping_2", 
        "very_detailed": "Merged_Grouping_3"
    };
    const resultColumn = resultColumnMap[detailLevel];

    // Column mapping for the expected CSV structure
    const columnMapping = {
        "cluster": "True Cell Type",
        "general_annotation": "Predicted Main Cell Type",
        "subtype_annotation": "Predicted Sub Cell Types"
    };

    // Validate required columns exist
    const availableColumns = Object.keys(csvData[0] || {});
    const missingColumns = Object.values(columnMapping).filter(col => !availableColumns.includes(col));
    if (missingColumns.length > 0) {
        throw new Error(`Required columns not found: ${missingColumns.join(', ')}. Available columns: ${availableColumns.join(', ')}`);
    }

    // Create working copy of data
    const workingData = csvData.map(row => ({ ...row }));

    // Process subtype column (extract first subtype if comma-separated)
    const subtypeCol = columnMapping["subtype_annotation"];
    workingData.forEach(row => {
        const subtypeValue = row[subtypeCol];
        if (typeof subtypeValue === 'string' && subtypeValue.includes(',')) {
            row["processed_subtype"] = subtypeValue.split(',')[0].trim();
        } else {
            row["processed_subtype"] = subtypeValue;
        }
    });

    // Initialize result column
    workingData.forEach(row => {
        row[resultColumn] = "";
    });

    // Process in batches
    const totalRows = workingData.length;
    for (let i = 0; i < totalRows; i += batchSize) {
        const batchEnd = Math.min(i + batchSize, totalRows);
        const batch = workingData.slice(i, batchEnd);

        if (onProgress) onProgress(`Processing clusters ${i + 1}-${batchEnd} of ${totalRows}...`);

        // Create prompt for this batch
        const prompt = _createAnnotationPrompt(
            batch,
            columnMapping["cluster"],
            columnMapping["general_annotation"],
            "processed_subtype",
            additionalContext,
            detailLevel
        );

        try {
            // Call LLM
            const response = await callLLM(
                prompt,
                provider,
                model,
                apiKey,
                0.3, // Lower temperature for consistency
                2000,
                "You are an expert cell biologist specializing in single-cell analysis. Your task is to analyze cluster annotations and suggest general cell groupings."
            );

            // Parse response and update data
            const groupings = _parseLLMResponse(response, batch);
            groupings.forEach((grouping, idx) => {
                workingData[i + idx][resultColumn] = grouping;
            });

        } catch (error) {
            console.error(`Error processing batch ${i + 1}-${batchEnd}: ${error.message}`);
            // Fill with error message for failed batch
            batch.forEach((_, idx) => {
                workingData[i + idx][resultColumn] = "Error processing";
            });
        }
    }

    // Add result column to original data
    csvData.forEach((row, idx) => {
        row[resultColumn] = workingData[idx][resultColumn];
    });

    if (onProgress) onProgress(`Completed ${detailLevel} annotation merging for ${totalRows} clusters`);

    return csvData;
}

/**
 * Process all three detail levels in parallel.
 * 
 * @param {Array} csvData - Array of objects containing cluster annotations
 * @param {string} provider - LLM provider to use
 * @param {string|null} model - Specific model to use
 * @param {string} apiKey - API key for the provider
 * @param {string|null} additionalContext - Optional domain-specific context
 * @param {number} batchSize - Number of clusters to process in each LLM call
 * @param {function} onProgress - Optional callback for progress updates
 * @returns {Promise<Array>} Array of objects with all three levels of groupings
 */
export async function mergeAnnotationsAll({
    csvData,
    provider = "openrouter",
    model = null,
    apiKey,
    additionalContext = null,
    batchSize = 20,
    onProgress = null
}) {
    if (onProgress) onProgress("Processing all three detail levels in parallel...");

    const detailLevels = ["broad", "detailed", "very_detailed"];
    
    // Run all three levels in parallel
    const promises = detailLevels.map(detailLevel =>
        mergeAnnotations({
            csvData: csvData.map(row => ({ ...row })), // Deep copy for each process
            provider,
            model,
            apiKey,
            additionalContext,
            batchSize,
            detailLevel,
            onProgress: onProgress ? (msg) => onProgress(`[${detailLevel}] ${msg}`) : null
        }).then(result => ({ detailLevel, result }))
    );

    try {
        const results = await Promise.all(promises);
        
        // Combine results
        let combinedData = csvData.map(row => ({ ...row }));
        const resultColumnMap = {
            "broad": "Merged_Grouping_1",
            "detailed": "Merged_Grouping_2",
            "very_detailed": "Merged_Grouping_3"
        };

        for (const { detailLevel, result } of results) {
            const resultColumn = resultColumnMap[detailLevel];
            result.forEach((row, idx) => {
                if (idx < combinedData.length) {
                    combinedData[idx][resultColumn] = row[resultColumn];
                }
            });
        }

        if (onProgress) onProgress("All detail levels completed successfully!");
        return combinedData;

    } catch (error) {
        throw new Error(`Parallel processing failed: ${error.message}`);
    }
}

/**
 * Create a prompt for the LLM based on the detail level and batch data.
 */
function _createAnnotationPrompt(
    batch,
    clusterCol,
    generalCol,
    subtypeCol,
    additionalContext,
    detailLevel = "broad"
) {
    // Format cluster data
    const clustersText = batch.map(row =>
        `Cluster ${row[clusterCol]}: General annotation: ${row[generalCol]}, Subtype: ${row[subtypeCol]}`
    ).join('\n');

    let prompt;

    if (detailLevel === "broad") {
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

Annotations to process:
${clustersText}

Please respond with a JSON object where keys are cluster identifiers and values are the suggested groupings.
For example:
{
  "1": "Macrophages",
  "2": "CD4 T cells",
  "3": "CD8 T cells"
}`;

    } else { // very_detailed
        prompt = `I have single-cell RNA-seq cluster annotations and need to normalize and standardize cell type names 
while preserving the most specific and detailed biological information.
For each cluster, I'll provide the general annotation and subtype annotation.

Your task is to create a consistent and standardized cell type label that:
1. Maintains the highest level of biological specificity from the annotations
2. Uses consistent nomenclature across similar cell types
3. Follows standard cell type naming conventions
4. Preserves functional or activation state information when present
5. Normalizes naming variations

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
 */
function _parseLLMResponse(response, batch) {
    const groupings = new Array(batch.length).fill("Error parsing response");

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
        // Keep default "Error parsing response" values
    }

    return groupings;
}

export default { mergeAnnotations, mergeAnnotationsAll };