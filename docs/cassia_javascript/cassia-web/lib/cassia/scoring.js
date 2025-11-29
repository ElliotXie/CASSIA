import { callLLM } from './llm_utils.js';
import { parseCSV, formatAsCSV } from '../utils/csv-parser.js';

// ----------------- Core Scoring Functions -----------------

/**
 * Create the scoring prompt exactly like Python version
 * @param {string} majorClusterInfo - Information about species and tissue
 * @param {string} marker - Comma-separated list of marker genes
 * @param {string} annotationHistory - Complete conversation history
 * @returns {string} Scoring prompt
 */
function promptCreatorScore(majorClusterInfo, marker, annotationHistory) {
    const prompt = `
        You are an expert in single-cell annotation analysis. Your task is to evaluate and rate single-cell annotation results, focusing on their correctness and ability to capture the overall picture of the data. You will provide a score from 0 to 100 and justify your rating.

Here are the single-cell annotation results to evaluate:


<marker>
${marker}
</marker>

<Cluster Origin>
${majorClusterInfo}
</Cluster Origin>

<annotation_history>
${annotationHistory}
</annotation_history>

Carefully analyze these results, paying particular attention to the following aspects:
1. Correctness of the annotations
2. Balanced consideration of multiple markers rather than over-focusing on a specific one
3. Ability to capture the general picture of the cell populations

When evaluating, consider:
- Are the annotations scientifically accurate?
- Is there a good balance in the use of different markers?
- Does the annotation provide a comprehensive view of the cell types present?
- Are there any obvious misclassifications or oversights?
- Did it consider the rank of the marker? marker appear first is more important.

Provide your analysis in the following format:
1. Start with a <reasoning> tag, where you explain your evaluation of the annotation results. Discuss the strengths and weaknesses you've identified, referring to specific examples from the results where possible.
2. After your reasoning, use a <score> tag to provide a numerical score from 0 to 100, where 0 represents completely incorrect or unusable results, and 100 represents perfect annotation that captures all aspects of the data correctly.

Your response should look like this:

<reasoning>
[Your detailed analysis and justification here]
</reasoning>

<score>[Your numerical score between 0 and 100]</score>

Remember, the focus is on correctness and the ability to see the general picture, rather than the structure of the results. Be critical but fair in your assessment.
    `;
    return prompt;
}

/**
 * Extract score and reasoning from LLM response with multiple fallback patterns
 * @param {string} text - LLM response text
 * @returns {Object} {score: number|null, reasoning: string|null}
 */
function extractScoreAndReasoning(text) {
    try {
        let score = null;
        let reasoning = null;
        
        // Extract score - try multiple patterns (exactly matching Python)
        const scorePatterns = [
            /\<score\>(\d+)\<\/score\>/i,         // <score>85</score>
            /Score:\s*(\d+)/i,                    // "Score: 85"
            /score:\s*(\d+)/i,                    // "score: 85"
            /(\d+)\/100/i,                        // "85/100"
            /(\d+)\s*out\s*of\s*100/i,           // "85 out of 100"
            /rating.*?(\d+)/i,                    // "rating of 85"
            /(\d+)%/i                             // "85%"
        ];
        
        for (const pattern of scorePatterns) {
            const scoreMatch = text.match(pattern);
            if (scoreMatch) {
                score = parseInt(scoreMatch[1]);
                break;
            }
        }
        
        // Extract reasoning - try multiple patterns (exactly matching Python)
        const reasoningPatterns = [
            /\<reasoning\>(.*?)\<\/reasoning\>/is,                    // <reasoning>...</reasoning>
            /Reasoning:\s*(.*?)(?=Score:|$)/is,                      // "Reasoning: ..." until "Score:" or end
            /reasoning:\s*(.*?)(?=score:|$)/is,                      // lowercase version
            /Analysis:\s*(.*?)(?=Score:|$)/is,                       // "Analysis: ..."
            /Evaluation:\s*(.*?)(?=Score:|$)/is                      // "Evaluation: ..."
        ];
        
        for (const pattern of reasoningPatterns) {
            const reasoningMatch = text.match(pattern);
            if (reasoningMatch) {
                reasoning = reasoningMatch[1].trim();
                break;
            }
        }
        
        // If no specific reasoning found, use the entire text as reasoning
        if (reasoning === null && text.trim()) {
            reasoning = text.trim();
        }
        
        return { score, reasoning };
        
    } catch (error) {
        console.error(`Error extracting data: ${error.message}`);
        return { score: null, reasoning: null };
    }
}

/**
 * Score a single cell type annotation analysis
 * @param {string} majorClusterInfo - Information about species and tissue
 * @param {string} marker - Comma-separated list of marker genes
 * @param {string} annotationHistory - History of annotation conversation
 * @param {string} model - Model to use (default: "deepseek/deepseek-chat-v3-0324")
 * @param {string} provider - AI provider to use ('openai', 'anthropic', or 'openrouter')
 * @param {string} apiKey - API key for the provider
 * @returns {Promise<Object>} {score: number, reasoning: string}
 */
export async function scoreSingleAnalysis(
    majorClusterInfo, 
    marker, 
    annotationHistory, 
    model = "deepseek/deepseek-chat-v3-0324", 
    provider = "openrouter",
    apiKey
) {
    const prompt = promptCreatorScore(majorClusterInfo, marker, annotationHistory);
    
    // Add explicit max_tokens to ensure responses aren't truncated
    const response = await callLLM(
        prompt,
        provider,
        model,
        apiKey,
        0.7,  // temperature
        2000, // max_tokens - ensure enough tokens for reasoning + score
        null  // system_prompt
    );
    
    const { score, reasoning } = extractScoreAndReasoning(response);
    return { score, reasoning };
}

// ----------------- CSV Helper Functions -----------------
// Note: CSV parsing functions have been moved to ../utils/csv-parser.js for reusability

// ----------------- Batch Processing Functions -----------------

/**
 * Process a single row of scoring data with retry logic
 * @param {Object} row - Row data containing columns like Species, Tissue, etc.
 * @param {number} idx - Row index for logging
 * @param {string} model - Model to use
 * @param {string} provider - AI provider to use
 * @param {string} apiKey - API key
 * @param {number} maxRetriesForNone - Max retries for None scores
 * @returns {Promise<Object>} {score, reasoning}
 */
async function processSingleRow(row, idx, model = "deepseek/deepseek-chat-v3-0324", provider = "openrouter", apiKey, maxRetriesForNone = 3) {
    try {
        // Robust column detection for Species and Tissue
        const findColumn = (options) => {
            for (const col of options) {
                if (col in row && row[col]) return row[col];
            }
            // Try case-insensitive match
            const availableColumns = Object.keys(row);
            for (const availCol of availableColumns) {
                const normalizedAvail = availCol.toLowerCase().replace(/[\s._-]+/g, '');
                for (const option of options) {
                    const normalizedOption = option.toLowerCase().replace(/[\s._-]+/g, '');
                    if (normalizedAvail === normalizedOption) {
                        return row[availCol];
                    }
                }
            }
            return null;
        };
        
        const species = findColumn(['Species', 'species', 'SPECIES', 'organism', 'Organism']) || 'Unknown';
        const tissue = findColumn(['Tissue', 'tissue', 'TISSUE', 'organ', 'Organ', 'sample', 'Sample']) || 'Unknown';
        const majorClusterInfo = `${species} ${tissue}`;
        
        console.log(`Processing row ${idx + 1}: ${majorClusterInfo}`);
        
        // Handle marker column with comprehensive naming options and fuzzy matching
        const markerColumnOptions = [
            'Marker List', 'Marker.List', 'marker_list', 'Marker_List', 'Marker-List',
            'MarkerList', 'markerlist', 'marker list', 'MARKER LIST',
            'markers', 'Markers', 'marker', 'Marker', 'genes', 'Genes',
            'Gene List', 'Gene.List', 'gene_list', 'Gene_List', 'gene-list',
            'GeneList', 'genelist', 'gene list', 'GENE LIST'
        ];
        
        let marker = null;
        let markerColumn = null;
        
        // First try exact match
        for (const col of markerColumnOptions) {
            if (col in row && row[col] && row[col].trim() !== '') {
                marker = row[col].trim();
                markerColumn = col;
                break;
            }
        }
        
        // If no exact match, try case-insensitive and normalized match
        if (!marker) {
            const availableColumns = Object.keys(row);
            for (const availCol of availableColumns) {
                const normalizedAvail = availCol.toLowerCase().replace(/[\s._-]+/g, '');
                for (const option of markerColumnOptions) {
                    const normalizedOption = option.toLowerCase().replace(/[\s._-]+/g, '');
                    if ((normalizedAvail === normalizedOption || normalizedAvail.includes('marker') || normalizedAvail.includes('gene')) 
                        && row[availCol] && row[availCol].trim() !== '') {
                        marker = row[availCol].trim();
                        markerColumn = availCol;
                        break;
                    }
                }
                if (marker) break;
            }
        }
        
        if (!marker || marker.trim() === '') {
            const availableColumns = Object.keys(row);
            console.error(`Row ${idx + 1} data for debugging:`, {
                availableColumns,
                'Marker List': row['Marker List'],
                'marker_list': row['marker_list'],
                'Markers': row['Markers'],
                markerColumnFound: markerColumn,
                markerValue: marker
            });
            
            // Try to find any column that might contain marker data
            const possibleMarkerColumns = availableColumns.filter(col => 
                col.toLowerCase().includes('marker') || 
                col.toLowerCase().includes('gene') ||
                col.toLowerCase().includes('list')
            );
            
            console.error(`Possible marker columns found:`, possibleMarkerColumns.map(col => ({
                column: col,
                value: row[col] ? row[col].substring(0, 100) : 'EMPTY'
            })));
            
            throw new Error(`Could not find marker column with data. Available columns: ${availableColumns.join(', ')}. Looking for columns like: ${markerColumnOptions.slice(0, 5).join(', ')}`);
        }
        
        console.log(`Found marker column: "${markerColumn}" with value length: ${marker.length}, preview: ${marker.substring(0, 100)}...`);
        
        // Handle conversation history column with comprehensive naming options
        // This is the key column that contains the annotation conversation from runCASSIA_batch
        const historyColumnOptions = [
            'Conversation History', 'Conversation.History', 'conversation_history', 
            'Conversation_History', 'Conversation-History', 'ConversationHistory',
            'conversationhistory', 'conversation history', 'CONVERSATION HISTORY',
            'annotation_history', 'Annotation History', 'Annotation.History',
            'annotation_conversation', 'Annotation Conversation',
            'conversation', 'Conversation', 'history', 'History',
            'chat_history', 'Chat History', 'chat-history', 'ChatHistory'
        ];
        
        let annotationHistory = null;
        let historyColumn = null;
        
        // First try exact match
        for (const col of historyColumnOptions) {
            if (col in row && row[col] && row[col].trim() !== '') {
                annotationHistory = row[col].trim();
                historyColumn = col;
                break;
            }
        }
        
        // If no exact match, try case-insensitive and normalized match
        if (!annotationHistory) {
            const availableColumns = Object.keys(row);
            for (const availCol of availableColumns) {
                const normalizedAvail = availCol.toLowerCase().replace(/[\s._-]+/g, '');
                for (const option of historyColumnOptions) {
                    const normalizedOption = option.toLowerCase().replace(/[\s._-]+/g, '');
                    if ((normalizedAvail === normalizedOption || 
                        (normalizedAvail.includes('conversation') && normalizedAvail.includes('history')) ||
                        normalizedAvail.includes('conversationhistory') ||
                        normalizedAvail.includes('annotationhistory'))
                        && row[availCol] && row[availCol].trim() !== '') {
                        annotationHistory = row[availCol].trim();
                        historyColumn = availCol;
                        break;
                    }
                }
                if (annotationHistory) break;
            }
        }
        
        if (!annotationHistory || annotationHistory.trim() === '') {
            const availableColumns = Object.keys(row);
            const hasConversationColumn = availableColumns.some(col => 
                col.toLowerCase().includes('conversation') || 
                col.toLowerCase().includes('history')
            );
            
            if (hasConversationColumn) {
                const conversationColumns = availableColumns.filter(col => 
                    col.toLowerCase().includes('conversation') || 
                    col.toLowerCase().includes('history')
                );
                console.error(`Row data:`, row);
                throw new Error(`Found conversation-related columns [${conversationColumns.join(', ')}] but they appear to be empty. Please ensure your CSV contains the conversation history from runCASSIA_batch analysis.`);
            } else {
                throw new Error(`Could not find conversation history column. This should be the output from runCASSIA_batch with 'Conversation History' column. Available columns: ${availableColumns.join(', ')}`);
            }
        }
        
        console.log(`Found history column: "${historyColumn}" with value length: ${annotationHistory.length} chars`);
        
        // Try up to 3 times for a valid score if we get null
        let score = null;
        let reasoning = null;
        let retryCount = 0;
        
        while (score === null && retryCount < maxRetriesForNone) {
            if (retryCount > 0) {
                console.log(`Retry ${retryCount}/${maxRetriesForNone} for row ${idx + 1} due to null score`);
            }
            
            const result = await scoreSingleAnalysis(
                majorClusterInfo,
                marker,
                annotationHistory,
                model,
                provider,
                apiKey
            );
            
            score = result.score;
            reasoning = result.reasoning;
            
            if (score !== null) {
                break;
            }
            
            retryCount++;
        }
        
        console.log(`Processed row ${idx + 1}: Score = ${score}`);
        return { score, reasoning };
        
    } catch (error) {
        console.error(`Error processing row ${idx + 1}: ${error.message}`);
        return { score: null, reasoning: `Error: ${error.message}` };
    }
}

/**
 * Main scoring function for batch processing
 * @param {Array|string} csvData - Array of objects or CSV string content
 * @param {Object} options - Configuration options
 * @returns {Promise<Object>} Results with scored data and CSV content
 */
export async function scoreAnnotationBatch({
    csvData,
    apiKey,
    maxWorkers = 4,
    model = "deepseek/deepseek-chat-v3-0324",
    provider = "openrouter",
    maxRetries = 1,
    onProgress = null,
    onLog = null
} = {}) {
    const logMessage = `🎯 Starting scoring process with ${maxWorkers} workers using ${provider} (${model})...`;
    console.log(logMessage);
    if (onLog) onLog(logMessage);
    
    try {
        // Parse CSV data if it's a string
        let results;
        if (typeof csvData === 'string') {
            results = parseCSV(csvData, { debug: true });
            console.log(`Parsed CSV: ${results.length} rows found`);
            if (results.length > 0) {
                console.log(`Available columns: ${Object.keys(results[0]).join(', ')}`);
                console.log(`First row sample:`, Object.entries(results[0]).slice(0, 5).map(([k, v]) => `${k}: "${String(v).substring(0, 30)}..."`).join(', '));
                
                // Debug ALL columns for first row to see misalignment
                console.log(`🔍 FULL FIRST ROW DEBUG:`);
                Object.entries(results[0]).forEach(([k, v]) => {
                    console.log(`  "${k}": "${String(v).substring(0, 100)}..."`);
                });
                
                // Debug marker column specifically
                const markerCol = results[0]['Marker List'];
                if (markerCol) {
                    console.log(`✅ Marker List column found with ${markerCol.length} characters`);
                } else {
                    console.log(`❌ Marker List column not found`);
                    console.log(`Available marker-like columns:`, Object.keys(results[0]).filter(k => k.toLowerCase().includes('marker') || k.toLowerCase().includes('gene')));
                }
                
                // Debug specific expected columns
                console.log(`🔍 EXPECTED COLUMNS DEBUG:`);
                console.log(`  "Species": "${results[0]['Species'] || 'NOT FOUND'}"`);
                console.log(`  "Tissue": "${results[0]['Tissue'] || 'NOT FOUND'}"`);
                console.log(`  "Marker Number": "${results[0]['Marker Number'] || 'NOT FOUND'}"`);
                console.log(`  "Marker List": "${(results[0]['Marker List'] || 'NOT FOUND').substring(0, 100)}..."`);
                console.log(`  "Conversation History": "${(results[0]['Conversation History'] || 'NOT FOUND').substring(0, 100)}..."`);
            }
        } else if (Array.isArray(csvData)) {
            results = [...csvData]; // Create a copy
        } else {
            throw new Error('csvData must be either a string or array');
        }
        
        if (results.length === 0) {
            throw new Error('No data found in the input');
        }
        
        // Initialize new columns if they don't exist
        results.forEach(row => {
            if (!('Score' in row)) {
                row.Score = null;
            }
            if (!('Scoring_Reasoning' in row)) {
                row.Scoring_Reasoning = null;
            }
        });
        
        // Create a list of unscored rows to process
        const rowsToProcess = results
            .map((row, idx) => ({ idx, row }))
            .filter(({ row }) => !row.Score || row.Score === null || row.Score === '');
        
        if (rowsToProcess.length === 0) {
            console.log("All rows already scored!");
            
            // Calculate statistics for UI display
            const validScores = results.filter(row => row.Score !== null && row.Score !== undefined && row.Score !== '').map(row => parseFloat(row.Score));
            const averageScore = validScores.length > 0 ? validScores.reduce((sum, score) => sum + score, 0) / validScores.length : 0;
            const highQuality = validScores.filter(score => score >= 80).length;
            
            return {
                results,
                csvContent: formatAsCSV(results),
                totalRows: results.length,
                processedRows: 0,
                alreadyScored: results.length,
                totalProcessed: validScores.length,
                averageScore: averageScore,
                highQuality: highQuality
            };
        }
        
        const processMessage = `📊 Found ${rowsToProcess.length} unscored rows to process...`;
        console.log(processMessage);
        if (onLog) onLog(processMessage);
        
        // Create semaphore for controlled concurrency
        const semaphore = new Semaphore(maxWorkers);
        let completedCount = 0;
        
        // Define a function that includes retry logic
        async function processWithRetry(rowData) {
            const { idx, row } = rowData;
            
            for (let attempt = 0; attempt <= maxRetries; attempt++) {
                try {
                    const result = await processSingleRow(row, idx, model, provider, apiKey);
                    
                    // Update the results array
                    results[idx].Score = result.score;
                    results[idx].Scoring_Reasoning = result.reasoning;
                    
                    completedCount++;
                    if (onProgress) {
                        onProgress({
                            completed: completedCount,
                            total: rowsToProcess.length,
                            percentage: Math.round((completedCount / rowsToProcess.length) * 100)
                        });
                    }
                    
                    return result;
                } catch (error) {
                    // Don't retry authentication errors
                    const errorStr = error.message.toLowerCase();
                    if (errorStr.includes("401") || errorStr.includes("api key") || errorStr.includes("authentication")) {
                        console.error(`⚠️  Row ${idx + 1} API ERROR: ${error.message}`);
                        console.error(`⚠️  This appears to be an API authentication error. Please check your API key.`);
                        // Return error info instead of throwing
                        results[idx].Score = null;
                        results[idx].Scoring_Reasoning = `API error: ${error.message}`;
                        return { score: null, reasoning: `API error: ${error.message}` };
                    }
                    
                    if (attempt < maxRetries) {
                        console.log(`Retrying row ${idx + 1} (attempt ${attempt + 2}/${maxRetries + 1})...`);
                        await new Promise(resolve => setTimeout(resolve, 1000)); // Wait 1 second before retry
                    } else {
                        console.error(`Row ${idx + 1} failed after ${maxRetries + 1} attempts: ${error.message}`);
                        results[idx].Score = null;
                        results[idx].Scoring_Reasoning = `Failed after ${maxRetries + 1} attempts: ${error.message}`;
                        return { score: null, reasoning: `Failed after ${maxRetries + 1} attempts: ${error.message}` };
                    }
                }
            }
        }
        
        // Process rows with controlled parallelism
        const promises = rowsToProcess.map(async (rowData) => {
            await semaphore.acquire();
            try {
                return await processWithRetry(rowData);
            } finally {
                semaphore.release();
            }
        });
        
        // Wait for all tasks to complete
        await Promise.all(promises);
        
        const completeMessage = `✅ Scoring complete! Processed ${rowsToProcess.length} rows.`;
        console.log(completeMessage);
        if (onLog) onLog(completeMessage);
        
        // Calculate statistics for UI display
        const validScores = results.filter(row => row.Score !== null && row.Score !== undefined && row.Score !== '').map(row => parseFloat(row.Score));
        const averageScore = validScores.length > 0 ? validScores.reduce((sum, score) => sum + score, 0) / validScores.length : 0;
        const highQuality = validScores.filter(score => score >= 80).length;
        
        return {
            results,
            csvContent: formatAsCSV(results),
            totalRows: results.length,
            processedRows: rowsToProcess.length,
            alreadyScored: results.length - rowsToProcess.length,
            totalProcessed: validScores.length,
            averageScore: averageScore,
            highQuality: highQuality
        };
        
    } catch (error) {
        console.error(`Error in scoring process: ${error.message}`);
        throw error;
    }
}

// ----------------- Semaphore for Concurrency Control -----------------

class Semaphore {
    constructor(max) {
        this.max = max;
        this.current = 0;
        this.queue = [];
    }
    
    async acquire() {
        return new Promise((resolve) => {
            if (this.current < this.max) {
                this.current++;
                resolve();
            } else {
                this.queue.push(resolve);
            }
        });
    }
    
    release() {
        this.current--;
        if (this.queue.length > 0) {
            const next = this.queue.shift();
            this.current++;
            next();
        }
    }
}

// ----------------- Exports -----------------

export default {
    scoreSingleAnalysis,
    scoreAnnotationBatch,
    promptCreatorScore,
    extractScoreAndReasoning
};