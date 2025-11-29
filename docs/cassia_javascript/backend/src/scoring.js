import { callLLM } from './llm_utils.js';
import { promises as fs, createReadStream } from 'fs';
import path from 'path';
import csvParser from 'csv-parser';
import { createObjectCsvWriter } from 'csv-writer';

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
 * @returns {Promise<Object>} {score: number, reasoning: string}
 */
export async function scoreSingleAnalysis(
    majorClusterInfo, 
    marker, 
    annotationHistory, 
    model = "deepseek/deepseek-chat-v3-0324", 
    provider = "openrouter"
) {
    const prompt = promptCreatorScore(majorClusterInfo, marker, annotationHistory);
    
    // Add explicit max_tokens to ensure responses aren't truncated
    const response = await callLLM(
        prompt,
        provider,
        model,
        null, // api_key - let it get from environment
        0.7,  // temperature
        2000, // max_tokens - ensure enough tokens for reasoning + score
        null  // system_prompt
    );
    
    const { score, reasoning } = extractScoreAndReasoning(response);
    return { score, reasoning };
}

// ----------------- CSV Helper Functions -----------------

/**
 * Read CSV file and return as array of objects
 */
async function readCSV(filePath) {
    return new Promise((resolve, reject) => {
        const results = [];
        createReadStream(filePath)
            .pipe(csvParser())
            .on('data', (data) => results.push(data))
            .on('end', () => resolve(results))
            .on('error', reject);
    });
}

/**
 * Write CSV file with headers and data
 */
async function writeCSV(filename, data) {
    if (!data || data.length === 0) {
        throw new Error("No data to write");
    }
    
    const outputDir = path.dirname(filename);
    
    // Create directory if it doesn't exist
    try {
        await fs.mkdir(outputDir, { recursive: true });
    } catch (error) {
        // Directory might already exist
    }
    
    // Get headers from first object
    const headers = Object.keys(data[0]);
    
    const csvWriter = createObjectCsvWriter({
        path: filename,
        header: headers.map(header => ({ id: header, title: header }))
    });
    
    await csvWriter.writeRecords(data);
}

// ----------------- Batch Processing Functions -----------------

/**
 * Process a single row of scoring data with retry logic
 * @param {Object} rowData - {idx, row} containing index and row data
 * @param {string} model - Model to use
 * @param {string} provider - AI provider to use
 * @param {number} maxRetriesForNone - Max retries for None scores
 * @returns {Promise<Object>} {idx, score, reasoning}
 */
async function processSingleRow(rowData, model = "deepseek/deepseek-chat-v3-0324", provider = "openrouter", maxRetriesForNone = 3) {
    const { idx, row } = rowData;
    
    try {
        const majorClusterInfo = `${row['Species']} ${row['Tissue']}`;
        
        // Handle both 'Marker List' and 'Marker.List' column names (flexible naming)
        const markerColumnOptions = ['Marker List', 'Marker.List', 'marker_list', 'Marker_List'];
        let marker = null;
        for (const col of markerColumnOptions) {
            if (col in row) {
                marker = row[col];
                break;
            }
        }
        if (marker === null) {
            throw new Error(`Could not find marker column. Available columns: ${Object.keys(row).join(', ')}`);
        }
        
        // Handle both 'Conversation History' and 'Conversation.History' column names
        const historyColumnOptions = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History'];
        let annotationHistory = null;
        for (const col of historyColumnOptions) {
            if (col in row) {
                annotationHistory = row[col];
                break;
            }
        }
        if (annotationHistory === null) {
            throw new Error(`Could not find conversation history column. Available columns: ${Object.keys(row).join(', ')}`);
        }
        
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
                provider
            );
            
            score = result.score;
            reasoning = result.reasoning;
            
            if (score !== null) {
                break;
            }
            
            retryCount++;
        }
        
        console.log(`Processed row ${idx + 1}: Score = ${score}`);
        return { idx, score, reasoning };
        
    } catch (error) {
        console.error(`Error processing row ${idx + 1}: ${error.message}`);
        return { idx, score: null, reasoning: `Error: ${error.message}` };
    }
}

/**
 * Process and score all rows in a results CSV file in parallel
 * @param {string} resultsFilePath - Path to the results CSV file
 * @param {string} outputFilePath - Path to save the updated results (optional)
 * @param {number} maxWorkers - Maximum number of parallel threads
 * @param {string} model - Model to use
 * @param {string} provider - AI provider to use
 * @returns {Promise<Array>} Original results with added score and reasoning columns
 */
export async function scoreAnnotationBatch(
    resultsFilePath, 
    outputFilePath = null, 
    maxWorkers = 4, 
    model = "deepseek/deepseek-chat-v3-0324", 
    provider = "openrouter"
) {
    // Read results file
    const results = await readCSV(resultsFilePath);
    
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
        return results;
    }
    
    console.log(`Processing ${rowsToProcess.length} unscored rows with up to ${maxWorkers} workers...`);
    
    // Create a semaphore for controlled concurrency
    const semaphore = new Semaphore(maxWorkers);
    
    // Process rows with controlled parallelism
    const promises = rowsToProcess.map(async (rowData) => {
        await semaphore.acquire();
        try {
            const result = await processSingleRow(rowData, model, provider);
            
            // Update the results array
            results[result.idx].Score = result.score;
            results[result.idx].Scoring_Reasoning = result.reasoning;
            
            // Save intermediate results
            if (outputFilePath === null) {
                outputFilePath = resultsFilePath.replace('.csv', '_scored.csv');
            }
            await writeCSV(outputFilePath, results);
            
            return result;
        } finally {
            semaphore.release();
        }
    });
    
    // Wait for all tasks to complete
    await Promise.all(promises);
    
    return results;
}

/**
 * Main scoring function with progress updates and retry logic
 * @param {string} inputFile - Path to input CSV file (with or without .csv extension)
 * @param {string} outputFile - Path to output CSV file (with or without .csv extension)
 * @param {number} maxWorkers - Maximum number of parallel workers
 * @param {string} model - Model to use
 * @param {string} provider - AI provider to use
 * @param {number} maxRetries - Maximum number of retries for failed analyses
 * @returns {Promise<Array>} Results array with scores
 */
export async function runCASSIAScoreBatch({
    inputFile,
    outputFile = null,
    maxWorkers = 4,
    model = "deepseek/deepseek-chat-v3-0324",
    provider = "openrouter",
    maxRetries = 1
} = {}) {
    // Add .csv extension if not present
    if (!inputFile.toLowerCase().endsWith('.csv')) {
        inputFile = inputFile + '.csv';
    }
    
    if (outputFile && !outputFile.toLowerCase().endsWith('.csv')) {
        outputFile = outputFile + '.csv';
    }
    
    console.log(`Starting scoring process with ${maxWorkers} workers using ${provider} (${model})...`);
    
    try {
        // Read the input file
        const results = await readCSV(inputFile);
        
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
            return results;
        }
        
        console.log(`Found ${rowsToProcess.length} unscored rows to process...`);
        
        // Create semaphore for controlled concurrency
        const semaphore = new Semaphore(maxWorkers);
        
        // Define a function that includes retry logic
        async function processWithRetry(rowData) {
            const { idx, row } = rowData;
            
            for (let attempt = 0; attempt <= maxRetries; attempt++) {
                try {
                    return await processSingleRow(rowData, model, provider);
                } catch (error) {
                    // Don't retry authentication errors
                    const errorStr = error.message.toLowerCase();
                    if (errorStr.includes("401") || errorStr.includes("api key") || errorStr.includes("authentication")) {
                        console.error(`⚠️  Row ${idx + 1} API ERROR: ${error.message}`);
                        console.error(`⚠️  This appears to be an API authentication error. Please check your API key.`);
                        // Return error info instead of throwing
                        return { idx, score: null, reasoning: `API error: ${error.message}` };
                    }
                    
                    if (attempt < maxRetries) {
                        console.log(`Retrying row ${idx + 1} (attempt ${attempt + 2}/${maxRetries + 1})...`);
                        await new Promise(resolve => setTimeout(resolve, 1000)); // Wait 1 second before retry
                    } else {
                        console.error(`Row ${idx + 1} failed after ${maxRetries + 1} attempts: ${error.message}`);
                        return { idx, score: null, reasoning: `Failed after ${maxRetries + 1} attempts: ${error.message}` };
                    }
                }
            }
        }
        
        // Process rows with controlled parallelism
        const promises = rowsToProcess.map(async (rowData) => {
            await semaphore.acquire();
            try {
                const result = await processWithRetry(rowData);
                
                // Update the results array
                results[result.idx].Score = result.score;
                results[result.idx].Scoring_Reasoning = result.reasoning;
                
                // Save intermediate results
                if (outputFile === null) {
                    outputFile = inputFile.replace('.csv', '_scored.csv');
                }
                await writeCSV(outputFile, results);
                
                return result;
            } finally {
                semaphore.release();
            }
        });
        
        // Wait for all tasks to complete
        await Promise.all(promises);
        
        console.log(`Scoring complete! Results saved to: ${outputFile}`);
        return results;
        
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
    runCASSIAScoreBatch,
    promptCreatorScore,
    extractScoreAndReasoning
};