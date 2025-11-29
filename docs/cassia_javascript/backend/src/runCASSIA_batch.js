import { runCASSIA } from './runCASSIA.js';
import { promises as fs, createReadStream } from 'fs';
import path from 'path';
import csvParser from 'csv-parser';
import { createObjectCsvWriter } from 'csv-writer';

// ----------------- Helper Functions -----------------

/**
 * Safely navigate nested objects with fallback to null
 */
function safeGet(obj, ...keys) {
    for (const key of keys) {
        if (obj && typeof obj === 'object' && key in obj) {
            obj = obj[key];
        } else {
            return null;
        }
    }
    return obj;
}

/**
 * Intelligent marker string splitting with multiple fallback strategies
 */
function splitMarkers(markerString) {
    if (!markerString) return [];
    
    let markers;
    
    // Strategy 1: Split by comma + space
    if (markerString.includes(', ')) {
        markers = markerString.split(/,\s*/);
    }
    // Strategy 2: Split by comma only
    else if (markerString.includes(',')) {
        markers = markerString.split(',');
    }
    // Strategy 3: Split by space
    else {
        markers = markerString.split(' ');
    }
    
    // Remove empty strings and trim whitespace
    return markers.map(m => m.trim()).filter(m => m.length > 0);
}

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
 * Process marker data to extract top genes per cluster
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

/**
 * Auto-detect column names for celltype and gene columns
 */
function detectColumns(df) {
    const columns = Object.keys(df[0] || {});
    
    let celltypeColumn = null;
    let geneColumn = null;
    
    // Look for common celltype column names
    const celltypeNames = ['cluster', 'celltype', 'cell_type', 'cell type', 'cell_cluster'];
    for (const name of celltypeNames) {
        if (columns.find(col => col.toLowerCase() === name.toLowerCase())) {
            celltypeColumn = columns.find(col => col.toLowerCase() === name.toLowerCase());
            break;
        }
    }
    
    // Look for common gene/marker column names
    const geneNames = ['markers', 'marker', 'genes', 'gene', 'marker_genes'];
    for (const name of geneNames) {
        if (columns.find(col => col.toLowerCase() === name.toLowerCase())) {
            geneColumn = columns.find(col => col.toLowerCase() === name.toLowerCase());
            break;
        }
    }
    
    // Fallback: use first two columns
    if (!celltypeColumn) celltypeColumn = columns[0];
    if (!geneColumn) geneColumn = columns[1];
    
    return { celltypeColumn, geneColumn };
}

/**
 * Write CSV file with headers and data
 */
async function writeCSV(filename, headers, rowData) {
    const outputDir = path.dirname(filename);
    
    // Create directory if it doesn't exist
    try {
        await fs.mkdir(outputDir, { recursive: true });
    } catch (error) {
        // Directory might already exist
    }
    
    const csvWriter = createObjectCsvWriter({
        path: filename,
        header: headers.map((header, index) => ({ id: `col${index}`, title: header }))
    });
    
    // Convert row data to objects
    const records = rowData.map(row => {
        const record = {};
        row.forEach((value, index) => {
            record[`col${index}`] = value;
        });
        return record;
    });
    
    await csvWriter.writeRecords(records);
}

/**
 * Sleep function for delays
 */
function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}

// ----------------- Main Batch Function -----------------

/**
 * Run cell type analysis on multiple clusters in parallel.
 * 
 * @param {Array|string} marker - Input DataFrame (array of objects) or CSV file path
 * @param {string} outputName - Base output filename
 * @param {number} nGenes - Number of top genes per cluster
 * @param {string} model - LLM model to use
 * @param {number} temperature - Model temperature
 * @param {string} tissue - Tissue type
 * @param {string} species - Species
 * @param {string} additionalInfo - Additional context information
 * @param {string} celltypeColumn - Column name containing cell types
 * @param {string} geneColumnName - Column name containing marker genes
 * @param {number} maxWorkers - Maximum number of parallel workers
 * @param {string} provider - AI provider to use
 * @param {number} maxRetries - Maximum retry attempts per cluster
 * @param {string} rankingMethod - Gene ranking method
 * @param {boolean} ascending - Sort direction
 * @param {string} validatorInvolvement - Validator version
 * @returns {Promise<object>} Results summary
 */
export async function runCASSIABatch({
    marker,
    outputName = "cell_type_analysis_results",
    nGenes = 50,
    model = "google/gemini-2.5-flash-preview",
    temperature = 0,
    tissue = "lung",
    species = "human",
    additionalInfo = null,
    celltypeColumn = null,
    geneColumnName = null,
    maxWorkers = 10,
    provider = "openrouter",
    maxRetries = 1,
    rankingMethod = "avg_log2FC",
    ascending = null,
    validatorInvolvement = "v1",
    formatType = null
} = {}) {
    console.log(`Starting CASSIA batch analysis with ${maxWorkers} workers...`);
    
    // Input validation and loading
    let df;
    if (typeof marker === 'string') {
        // Load from CSV file
        console.log(`Loading marker data from ${marker}...`);
        df = await readCSV(marker);
    } else if (Array.isArray(marker)) {
        // Use provided DataFrame
        df = marker;
    } else {
        throw new Error("Marker must be either a file path (string) or an array of objects");
    }
    
    if (!df || df.length === 0) {
        throw new Error("No data found in the marker input");
    }
    
    console.log(`Loaded ${df.length} rows of marker data`);
    
    // Data processing
    if (Object.keys(df[0]).length > 2) {
        console.log("Processing differential expression data...");
        df = getTopMarkers(df, nGenes, rankingMethod, ascending, formatType);
        console.log(`Extracted top ${nGenes} markers for ${df.length} clusters`);
    }
    
    // Auto-detect columns if not specified
    if (!celltypeColumn || !geneColumnName) {
        const detected = detectColumns(df);
        celltypeColumn = celltypeColumn || detected.celltypeColumn;
        geneColumnName = geneColumnName || detected.geneColumn;
        
        console.log(`Auto-detected columns: celltype='${celltypeColumn}', markers='${geneColumnName}'`);
    }
    
    // Validate columns exist
    if (!df[0][celltypeColumn] || !df[0][geneColumnName]) {
        throw new Error(`Specified columns not found. Available columns: ${Object.keys(df[0]).join(', ')}`);
    }
    
    // Worker function for individual cell type analysis
    async function analyzeCellType(cellType, markerList) {
        for (let attempt = 0; attempt <= maxRetries; attempt++) {
            try {
                const [result, conversationHistory] = await runCASSIA(
                    model,
                    temperature,
                    markerList,
                    tissue,
                    species,
                    additionalInfo,
                    provider,
                    validatorInvolvement
                );
                
                // Add metadata
                result.num_markers = markerList.length;
                result.marker_list = markerList;
                
                return [cellType, result, conversationHistory];
            } catch (error) {
                // Special handling for authentication errors (no retry)
                const errorStr = error.message.toLowerCase();
                if (errorStr.includes("401") || errorStr.includes("api key") || errorStr.includes("authentication")) {
                    console.error(`Authentication error for ${cellType}. Please check your API key.`);
                    throw error;
                }
                
                // Retry logic for other errors
                if (attempt < maxRetries) {
                    console.log(`Retrying analysis for ${cellType} (attempt ${attempt + 2}/${maxRetries + 1})...`);
                    await sleep(1000); // Wait 1 second before retry
                } else {
                    console.error(`${cellType} failed after ${maxRetries + 1} attempts with error: ${error.message}`);
                    throw error;
                }
            }
        }
    }
    
    // Prepare tasks
    const tasks = df.map(row => ({
        cellType: row[celltypeColumn],
        markerList: splitMarkers(row[geneColumnName])
    }));
    
    console.log(`Processing ${tasks.length} cell types with up to ${maxWorkers} parallel workers...`);
    
    // Execute tasks with controlled parallelism
    const results = {};
    const semaphore = new Semaphore(maxWorkers);
    
    const promises = tasks.map(async ({ cellType, markerList }) => {
        await semaphore.acquire();
        try {
            console.log(`Starting analysis for: ${cellType}`);
            const [resultCellType, result, conversationHistory] = await analyzeCellType(cellType, markerList);
            
            if (result) {
                results[resultCellType] = {
                    analysis_result: result,
                    conversation_history: conversationHistory,
                    iterations: result.iterations || 1
                };
                console.log(`Completed analysis for: ${cellType}`);
            }
        } catch (error) {
            console.error(`Failed analysis for ${cellType}: ${error.message}`);
        } finally {
            semaphore.release();
        }
    });
    
    // Wait for all tasks to complete
    await Promise.all(promises);
    
    console.log(`Batch analysis completed. ${Object.keys(results).length} successful analyses.`);
    
    // Process results for CSV output
    const fullData = [];
    const summaryData = [];
    
    for (const [trueCellType, details] of Object.entries(results)) {
        const mainCellType = safeGet(details, 'analysis_result', 'main_cell_type') || '';
        const subCellTypes = (safeGet(details, 'analysis_result', 'sub_cell_types') || []).join(', ');
        const possibleMixedCellTypes = (safeGet(details, 'analysis_result', 'possible_mixed_cell_types') || []).join(', ');
        const markerNumber = safeGet(details, 'analysis_result', 'num_markers') || '';
        const markerList = (safeGet(details, 'analysis_result', 'marker_list') || []).join(', ');
        const iterations = safeGet(details, 'analysis_result', 'iterations') || '';
        
        // Process conversation history
        const conversationHistory = details.conversation_history;
        let rawConversationHistory = '';
        
        if (conversationHistory && conversationHistory.all_iterations) {
            rawConversationHistory = conversationHistory.all_iterations
                .map(iter => 
                    iter.annotation.map(entry => `${entry[0]}: ${entry[1]}`).join(' | ')
                ).join(' | ');
        }
        
        // Full data row
        const fullRow = [
            trueCellType,
            mainCellType,
            subCellTypes,
            possibleMixedCellTypes,
            markerNumber,
            markerList,
            iterations,
            model,
            provider,
            tissue,
            species,
            additionalInfo || '',
            rawConversationHistory
        ];
        
        // Summary data row (without conversation history)
        const summaryRow = [
            trueCellType,
            mainCellType,
            subCellTypes,
            possibleMixedCellTypes,
            markerList,
            iterations,
            model,
            provider,
            tissue,
            species
        ];
        
        fullData.push(fullRow);
        summaryData.push(summaryRow);
    }
    
    // Sort data by cell type name
    fullData.sort((a, b) => a[0].localeCompare(b[0]));
    summaryData.sort((a, b) => a[0].localeCompare(b[0]));
    
    // Generate output filenames
    const baseName = path.parse(outputName).name;
    const outputDir = path.dirname(outputName);
    const fullCsvName = path.join(outputDir, `${baseName}_full.csv`);
    const summaryCsvName = path.join(outputDir, `${baseName}_summary.csv`);
    
    // Write CSV files
    const fullHeaders = [
        'True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types',
        'Possible Mixed Cell Types', 'Marker Number', 'Marker List', 'Iterations',
        'Model', 'Provider', 'Tissue', 'Species', 'Additional Info', 'Conversation History'
    ];
    
    const summaryHeaders = [
        'True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types',
        'Possible Mixed Cell Types', 'Marker List', 'Iterations', 'Model', 'Provider',
        'Tissue', 'Species'
    ];
    
    console.log(`Writing results to ${fullCsvName} and ${summaryCsvName}...`);
    
    await writeCSV(fullCsvName, fullHeaders, fullData);
    await writeCSV(summaryCsvName, summaryHeaders, summaryData);
    
    console.log("Batch analysis complete!");
    
    return {
        total_clusters: tasks.length,
        successful_analyses: Object.keys(results).length,
        failed_analyses: tasks.length - Object.keys(results).length,
        output_files: [fullCsvName, summaryCsvName],
        results: results
    };
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

export default { runCASSIABatch };