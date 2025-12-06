import { runCASSIA } from './runCASSIA.js';
import { generateBatchHtmlReportFromData } from './generateBatchReport.js';

// ----------------- Model Presets -----------------

/**
 * Model presets for runCASSIA_batch
 */
const MODEL_PRESETS = {
    performance: {
        name: 'Performance',
        description: 'Best quality results',
        provider: 'openrouter',
        model: 'anthropic/claude-sonnet-4'
    },
    balanced: {
        name: 'Balanced',
        description: 'Good quality with reasonable speed',
        provider: 'openrouter',
        model: 'google/gemini-2.5-flash'
    }
};

/**
 * Apply model preset to get provider and model configuration
 * @param {string} preset - Preset name ('performance', 'balanced', or null for manual)
 * @param {string} manualProvider - Manual provider if preset is null
 * @param {string} manualModel - Manual model if preset is null
 * @returns {Object} {provider, model} configuration
 */
function applyModelPreset(preset, manualProvider = 'openrouter', manualModel = 'google/gemini-2.5-flash') {
    if (preset && MODEL_PRESETS[preset]) {
        const presetConfig = MODEL_PRESETS[preset];
        return {
            provider: presetConfig.provider,
            model: presetConfig.model
        };
    }
    
    // Manual configuration
    return {
        provider: manualProvider,
        model: manualModel
    };
}

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

// CSV reading handled by frontend file processing

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
 * Convert data to CSV format for download with proper escaping
 */
function formatAsCSV(headers, rowData) {
    // Helper function to properly escape CSV cells
    const escapeCsvCell = (cell) => {
        if (cell === null || cell === undefined) {
            return '""';
        }
        
        const cellStr = String(cell);
        
        // If cell contains special characters, wrap in quotes and escape internal quotes
        if (cellStr.includes(',') || cellStr.includes('"') || cellStr.includes('\n') || cellStr.includes('\r')) {
            return `"${cellStr.replace(/"/g, '""')}"`;
        }
        
        // For simple cells, still wrap in quotes for consistency
        return `"${cellStr}"`;
    };
    
    const csvContent = [
        headers.map(escapeCsvCell).join(','),
        ...rowData.map(row => row.map(escapeCsvCell).join(','))
    ].join('\n');
    return csvContent;
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
    apiKey,
    outputName = "cell_type_analysis_results",
    nGenes = 50,
    model = "google/gemini-2.5-flash",
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
    formatType = null,
    preset = null, // New parameter for model presets
    onLog = null
} = {}) {
    // Apply model preset if specified
    const modelConfig = applyModelPreset(preset, provider, model);
    const finalProvider = modelConfig.provider;
    const finalModel = modelConfig.model;
    const startMessage = `🚀 ===== CASSIA BATCH ANALYSIS STARTED =====`;
    console.log(`\n${startMessage}`);
    if (onLog) onLog(startMessage);
    
    // Input validation and loading
    let df;
    if (Array.isArray(marker)) {
        // Use provided DataFrame
        df = marker;
    } else {
        throw new Error("Marker must be an array of objects in browser environment");
    }
    
    if (!df || df.length === 0) {
        throw new Error("No data found in the marker input");
    }
    
    // Now that we have df, we can show the full configuration
    const configMessage = `📋 Configuration: ${maxWorkers} workers, ${df.length} clusters`;
    const presetMessage = preset ? `🎛️ Using preset: ${MODEL_PRESETS[preset].name} (${MODEL_PRESETS[preset].description})` : '🔧 Manual configuration';
    const targetMessage = `🎯 Target: ${tissue} ${species} using ${finalProvider}/${finalModel}`;
    const settingsMessage = `⚙️ Settings: ${nGenes} genes, ${validatorInvolvement} validator`;
    const loadedMessage = `📊 Loaded ${df.length} rows of marker data`;
    
    console.log(configMessage);
    console.log(presetMessage);
    console.log(targetMessage);
    console.log(settingsMessage);
    console.log(`===============================================\n`);
    console.log(loadedMessage);
    
    if (onLog) {
        onLog(configMessage);
        onLog(presetMessage);
        onLog(targetMessage);
        onLog(settingsMessage);
        onLog(loadedMessage);
    }
    
    // Data processing
    if (Object.keys(df[0]).length > 2) {
        const processingDataMessage = "🔬 Processing differential expression data...";
        console.log(processingDataMessage);
        if (onLog) onLog(processingDataMessage);
        
        df = getTopMarkers(df, nGenes, rankingMethod, ascending, formatType);
        
        const extractedMessage = `🧬 Extracted top ${nGenes} markers for ${df.length} clusters`;
        console.log(extractedMessage);
        if (onLog) onLog(extractedMessage);
    }
    
    // Auto-detect columns if not specified
    if (!celltypeColumn || !geneColumnName) {
        const detected = detectColumns(df);
        celltypeColumn = celltypeColumn || detected.celltypeColumn;
        geneColumnName = geneColumnName || detected.geneColumn;
        
        const detectionMessage = `🔍 Auto-detected columns: celltype='${celltypeColumn}', markers='${geneColumnName}'`;
        console.log(detectionMessage);
        if (onLog) onLog(detectionMessage);
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
                    finalModel,
                    temperature,
                    markerList,
                    tissue,
                    species,
                    additionalInfo,
                    finalProvider,
                    validatorInvolvement,
                    apiKey
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
    
    const processingMessage = `⚙️ Processing ${tasks.length} cell types with up to ${maxWorkers} parallel workers...`;
    console.log(processingMessage);
    if (onLog) onLog(processingMessage);
    
    // Execute tasks with controlled parallelism
    const results = {};
    const semaphore = new Semaphore(maxWorkers);
    
    const promises = tasks.map(async ({ cellType, markerList }, index) => {
        await semaphore.acquire();
        try {
            const startAnalysisMessage = `[${index + 1}/${tasks.length}] 🧬 Starting analysis for: ${cellType} (${markerList.length} markers)`;
            console.log(startAnalysisMessage);
            if (onLog) onLog(startAnalysisMessage);
            const startTime = Date.now();
            
            const [resultCellType, result, conversationHistory] = await analyzeCellType(cellType, markerList);
            
            const duration = ((Date.now() - startTime) / 1000).toFixed(1);
            const completed = Object.keys(results).length + 1;
            
            if (result) {
                results[resultCellType] = {
                    analysis_result: result,
                    conversation_history: conversationHistory,
                    iterations: result.iterations || 1
                };
                
                const mainCellType = result.main_cell_type || 'Unknown';
                const completionMessage = `[${completed}/${tasks.length}] ✅ Completed ${cellType} → ${mainCellType} (${duration}s, ${result.iterations || 1} iterations)`;
                console.log(completionMessage);
                if (onLog) onLog(completionMessage);
                
                // Progress indicator
                const progress = Math.round((completed / tasks.length) * 100);
                const progressMessage = `📊 Progress: ${progress}% (${completed}/${tasks.length} clusters completed)`;
                console.log(progressMessage);
                if (onLog) onLog(progressMessage);
            }
        } catch (error) {
            console.error(`[${index + 1}/${tasks.length}] ❌ Failed analysis for ${cellType}: ${error.message}`);
            if (onLog) onLog(`❌ Failed analysis for ${cellType}: ${error.message}`);
            
            // Log more details about the error for debugging
            console.error('Error details:', error);
            if (error.stack) {
                console.error('Stack trace:', error.stack);
            }
        } finally {
            semaphore.release();
        }
    });
    
    // Wait for all tasks to complete
    await Promise.all(promises);
    
    const totalAnalyses = Object.keys(results).length;
    const failedAnalyses = tasks.length - totalAnalyses;
    
    const completionMessage = `🎉 ===== BATCH ANALYSIS COMPLETED =====`;
    const resultsMessage = `📊 Final results: ${totalAnalyses}/${tasks.length} successful, ${failedAnalyses} failed`;
    
    console.log(`\n${completionMessage}`);
    console.log(resultsMessage);
    
    if (onLog) {
        onLog(completionMessage);
        onLog(resultsMessage);
    }
    
    if (totalAnalyses > 0) {
        const avgIterations = Object.values(results)
            .reduce((sum, r) => sum + (r.iterations || 1), 0) / totalAnalyses;
        const avgMessage = `⚡ Average iterations per cluster: ${avgIterations.toFixed(1)}`;
        console.log(avgMessage);
        if (onLog) onLog(avgMessage);
    }
    console.log(`==========================================\n`);
    
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
        
        // Process conversation history with proper CSV escaping
        const conversationHistory = details.conversation_history;
        let rawConversationHistory = '';
        
        if (conversationHistory && conversationHistory.all_iterations) {
            rawConversationHistory = conversationHistory.all_iterations
                .map(iter => {
                    // Add robust type checking for iter.annotation
                    if (!iter.annotation) {
                        console.warn('⚠️ iter.annotation is missing for iteration');
                        return '';
                    }
                    
                    if (!Array.isArray(iter.annotation)) {
                        console.warn('⚠️ iter.annotation is not an array:', typeof iter.annotation);
                        return String(iter.annotation);
                    }
                    
                    return iter.annotation.map(entry => {
                        // Handle different entry formats
                        if (Array.isArray(entry) && entry.length >= 2) {
                            // Standard [role, content] format
                            const role = String(entry[0] || '').trim();
                            const content = String(entry[1] || '')
                                .replace(/\r?\n/g, ' ')  // Replace newlines with spaces
                                .replace(/\s+/g, ' ')    // Collapse multiple spaces
                                .replace(/"/g, '""')     // Escape quotes for CSV
                                .trim();
                            return `${role}: ${content}`;
                        } else if (typeof entry === 'object' && entry.role && entry.content) {
                            // Object format {role, content}
                            const role = String(entry.role || '').trim();
                            const content = String(entry.content || '')
                                .replace(/\r?\n/g, ' ')
                                .replace(/\s+/g, ' ')
                                .replace(/"/g, '""')
                                .trim();
                            return `${role}: ${content}`;
                        } else {
                            // Fallback for other formats
                            return String(entry).replace(/\r?\n/g, ' ').replace(/\s+/g, ' ').trim();
                        }
                    }).join(' | ');
                })
                .filter(item => item.length > 0) // Remove empty items
                .join(' | ');
            
            // Additional sanitization for the entire conversation history
            rawConversationHistory = rawConversationHistory
                .replace(/\r?\n/g, ' ')     // Remove any remaining newlines
                .replace(/\s+/g, ' ')       // Collapse spaces again
                .replace(/\|+/g, '|')       // Remove duplicate separators
                .trim();
            
            // Truncate extremely long conversation histories to prevent CSV bloat
            const maxConversationLength = 5000; // 5KB limit per conversation
            if (rawConversationHistory.length > maxConversationLength) {
                rawConversationHistory = rawConversationHistory.substring(0, maxConversationLength) + '... [truncated]';
            }
        }
        
        // Helper function to sanitize cell values
        const sanitizeCell = (value) => {
            if (value === null || value === undefined) return '';
            return String(value)
                .replace(/\r?\n/g, ' ')    // Replace newlines with spaces
                .replace(/\s+/g, ' ')      // Collapse multiple spaces
                .trim();
        };
        
        // Full data row with sanitized values
        const fullRow = [
            sanitizeCell(trueCellType),
            sanitizeCell(mainCellType),
            sanitizeCell(subCellTypes),
            sanitizeCell(possibleMixedCellTypes),
            sanitizeCell(markerNumber),
            sanitizeCell(markerList),
            sanitizeCell(iterations),
            sanitizeCell(model),
            sanitizeCell(provider),
            sanitizeCell(tissue),
            sanitizeCell(species),
            sanitizeCell(additionalInfo || ''),
            sanitizeCell(rawConversationHistory)
        ];
        
        // Summary data row (without conversation history) with sanitized values
        const summaryRow = [
            sanitizeCell(trueCellType),
            sanitizeCell(mainCellType),
            sanitizeCell(subCellTypes),
            sanitizeCell(possibleMixedCellTypes),
            sanitizeCell(markerList),
            sanitizeCell(iterations),
            sanitizeCell(model),
            sanitizeCell(provider),
            sanitizeCell(tissue),
            sanitizeCell(species)
        ];
        
        fullData.push(fullRow);
        summaryData.push(summaryRow);
    }
    
    // Sort data by cell type name
    fullData.sort((a, b) => a[0].localeCompare(b[0]));
    summaryData.sort((a, b) => a[0].localeCompare(b[0]));
    
    // Prepare CSV data for download
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
    
    console.log("📄 Preparing CSV data for download...");
    
    const fullCsv = formatAsCSV(fullHeaders, fullData);
    const summaryCsv = formatAsCSV(summaryHeaders, summaryData);
    
    const csvReadyMessage = `📄 CSV files ready! Full: ${(fullCsv.length / 1024).toFixed(1)}KB, Summary: ${(summaryCsv.length / 1024).toFixed(1)}KB`;

    console.log(csvReadyMessage);
    if (onLog) onLog(csvReadyMessage);

    // Generate HTML report (same as Python version)
    const htmlReportMessage = "📊 Generating HTML report...";
    console.log(htmlReportMessage);
    if (onLog) onLog(htmlReportMessage);

    const htmlRows = [];
    for (const [trueCellType, details] of Object.entries(results)) {
        const mainCellType = safeGet(details, 'analysis_result', 'main_cell_type') || '';
        const subCellTypes = (safeGet(details, 'analysis_result', 'sub_cell_types') || []).join(', ');
        const possibleMixedCellTypes = (safeGet(details, 'analysis_result', 'possible_mixed_cell_types') || []).join(', ');
        const markerNumber = safeGet(details, 'analysis_result', 'num_markers') || '';
        const markerList = (safeGet(details, 'analysis_result', 'marker_list') || []).join(', ');
        const iterations = safeGet(details, 'analysis_result', 'iterations') || '';

        // Process conversation history for HTML - pass as JSON string
        // The parseConversationHistory() function in generateBatchReport.js will handle parsing
        const conversationHistory = details.conversation_history;
        let rawConversationHistory = '';

        if (conversationHistory) {
            // Pass the full conversation history as JSON string
            // This preserves all data including validation_result and Formatting Agent
            rawConversationHistory = JSON.stringify(conversationHistory);
        }

        htmlRows.push({
            'Cluster ID': trueCellType,
            'Predicted General Cell Type': mainCellType,
            'Predicted Detailed Cell Type': subCellTypes,
            'Possible Mixed Cell Types': possibleMixedCellTypes,
            'Marker Number': markerNumber,
            'Marker List': markerList,
            'Iterations': iterations,
            'Model': finalModel,
            'Provider': finalProvider,
            'Tissue': tissue,
            'Species': species,
            'Additional Info': additionalInfo || 'N/A',
            'Conversation History': rawConversationHistory
        });
    }

    // Sort HTML rows by cluster ID
    htmlRows.sort((a, b) => String(a['Cluster ID']).localeCompare(String(b['Cluster ID'])));

    const htmlContent = generateBatchHtmlReportFromData(
        htmlRows,
        `CASSIA Batch Analysis - ${tissue} (${species})`
    );

    const htmlReadyMessage = `📊 HTML report ready! Size: ${(htmlContent.length / 1024).toFixed(1)}KB`;
    const batchCompleteMessage = "✅ Batch analysis complete!";

    console.log(htmlReadyMessage);
    console.log(batchCompleteMessage);

    if (onLog) {
        onLog(htmlReadyMessage);
        onLog(batchCompleteMessage);
    }

    return {
        total_clusters: tasks.length,
        successful_analyses: Object.keys(results).length,
        failed_analyses: tasks.length - Object.keys(results).length,
        config: {
            preset: preset,
            provider: finalProvider,
            model: finalModel,
            tissue: tissue,
            species: species,
            additionalInfo: additionalInfo,
            nGenes: nGenes,
            maxWorkers: maxWorkers,
            temperature: temperature
        },
        csv_data: {
            full: {
                filename: `${outputName}_full.csv`,
                content: fullCsv,
                headers: fullHeaders,
                data: fullData
            },
            summary: {
                filename: `${outputName}_summary.csv`,
                content: summaryCsv,
                headers: summaryHeaders,
                data: summaryData
            }
        },
        html_report: {
            filename: `${outputName}_report.html`,
            content: htmlContent,
            type: 'text/html',
            size: htmlContent.length
        },
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

export { MODEL_PRESETS };
export default { runCASSIABatch };