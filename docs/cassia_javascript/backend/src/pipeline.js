import fs from 'fs';
import path from 'path';
import { parse as parseCSV } from 'csv-parse/sync';
import { stringify as stringifyCSV } from 'csv-stringify/sync';
import { runCASSIABatch } from './runCASSIA_batch.js';
import { runCASSIAScoreBatch } from './scoring.js';
import { runCASSIAGenerateScoreReport } from './report_generation.js';
import { runCASSIAAnnotationboost } from './annotationBoost.js';
import { mergeAnnotationsAll } from './mergingAnnotation.js';

/**
 * Run the complete cell analysis pipeline including annotation, scoring, and report generation.
 * This is a 100% JavaScript replication of the Python runCASSIA_pipeline function.
 * 
 * @param {string} outputFileName - Base name for output files
 * @param {string} tissue - Tissue type being analyzed
 * @param {string} species - Species being analyzed
 * @param {Array|string} marker - Marker data (array of objects or path to CSV file)
 * @param {number} maxWorkers - Maximum number of concurrent workers
 * @param {string} annotationModel - Model to use for initial annotation
 * @param {string} annotationProvider - Provider for initial annotation
 * @param {string} scoreModel - Model to use for scoring
 * @param {string} scoreProvider - Provider for scoring
 * @param {string} annotationboostModel - Model to use for boosting low-scoring annotations
 * @param {string} annotationboostProvider - Provider for boosting low-scoring annotations
 * @param {number} scoreThreshold - Threshold for identifying low-scoring clusters
 * @param {string} additionalInfo - Additional information for analysis
 * @param {number} maxRetries - Maximum number of retries for failed analyses
 * @param {boolean} mergeAnnotations - Whether to merge annotations from LLM
 * @param {string} mergeModel - Model to use for merging annotations
 * @param {string} mergeProvider - Provider to use for merging annotations
 * @param {string} conversationHistoryMode - Mode for extracting conversation history ("full", "final", or "none")
 * @param {string} rankingMethod - Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
 * @param {boolean|null} ascending - Sort direction (null uses default for each method)
 * @param {string} reportStyle - Style of report generation ("per_iteration" or "total_summary")
 * @param {string} validatorInvolvement - Validator involvement level
 * 
 * @returns {Promise<string>} Path to the final combined results file
 */
export async function runCASSIAPipeline({
    outputFileName,
    tissue,
    species,
    marker,
    maxWorkers = 4,
    annotationModel = "meta-llama/llama-4-maverick",
    annotationProvider = "openrouter",
    scoreModel = "google/gemini-2.5-pro-preview-03-25",
    scoreProvider = "openrouter",
    annotationboostModel = "google/gemini-2.5-flash-preview",
    annotationboostProvider = "openrouter",
    scoreThreshold = 75,
    additionalInfo = "None",
    maxRetries = 1,
    mergeAnnotations = true,
    mergeModel = "deepseek/deepseek-chat-v3-0324",
    mergeProvider = "openrouter",
    conversationHistoryMode = "final",
    rankingMethod = "avg_log2FC",
    ascending = null,
    reportStyle = "per_iteration",
    validatorInvolvement = "v1"
}) {
    // Create a main folder based on tissue and species for organizing reports
    let mainFolderName = `CASSIA_${tissue}_${species}`;
    mainFolderName = mainFolderName.replace(/[^a-zA-Z0-9\s\-_]/g, '').trim();
    mainFolderName = mainFolderName.replace(/\s+/g, '_');
    
    // Remove .csv extension if present
    if (outputFileName.toLowerCase().endsWith('.csv')) {
        outputFileName = outputFileName.slice(0, -4); // Remove last 4 characters (.csv)
    }

    // Add timestamp to prevent overwriting existing folders with the same name
    const timestamp = new Date().toISOString().replace(/[:.]/g, '-').slice(0, 19);
    mainFolderName = `${mainFolderName}_${timestamp}`;
    
    // Create the main folder if it doesn't exist
    if (!fs.existsSync(mainFolderName)) {
        fs.mkdirSync(mainFolderName, { recursive: true });
        console.log(`Created main folder: ${mainFolderName}`);
    }
        
    // Create organized subfolders according to user's specifications
    const annotationResultsFolder = path.join(mainFolderName, "01_annotation_results");  // All CSV files
    const reportsFolder = path.join(mainFolderName, "02_reports");  // All HTML reports except annotation boost
    const boostFolder = path.join(mainFolderName, "03_boost_analysis");   // All annotation boost related results
    
    // Create all subfolders
    const folders = [annotationResultsFolder, reportsFolder, boostFolder];
    for (const folder of folders) {
        if (!fs.existsSync(folder)) {
            fs.mkdirSync(folder, { recursive: true });
            console.log(`Created subfolder: ${folder}`);
        }
    }
    
    // Define derived file names with folder paths
    // All CSV files go to the annotation_results_folder
    const rawFullCsv = path.join(annotationResultsFolder, `${outputFileName}_full.csv`);
    const rawSummaryCsv = path.join(annotationResultsFolder, `${outputFileName}_summary.csv`);
    const rawSortedCsv = path.join(annotationResultsFolder, `${outputFileName}_sorted_full.csv`);
    const scoreFileName = path.join(annotationResultsFolder, `${outputFileName}_scored.csv`);
    const mergedAnnotationFile = path.join(annotationResultsFolder, `${outputFileName}_merged.csv`);
    
    // Reports go to the reports_folder - ALL HTML reports should be in this folder
    const reportBaseName = path.join(reportsFolder, outputFileName);
    
    // First annotation output is in the current directory but will be moved later
    const annotationOutput = outputFileName;

    console.log("\n=== Starting cell type analysis ===");
    // Run initial cell type analysis
    await runCASSIABatch({
        marker: marker,
        outputName: annotationOutput,
        model: annotationModel,
        tissue: tissue,
        species: species,
        additionalInfo: additionalInfo,
        provider: annotationProvider,
        maxWorkers: maxWorkers,
        maxRetries: maxRetries,
        rankingMethod: rankingMethod,
        ascending: ascending,
        validatorInvolvement: validatorInvolvement
    });
    console.log("✓ Cell type analysis completed");
    
    // Copy the generated files to the organized folders
    const originalFullCsv = annotationOutput + "_full.csv";
    const originalSummaryCsv = annotationOutput + "_summary.csv";
    
    // Copy the files if they exist
    if (fs.existsSync(originalFullCsv)) {
        // Read and write instead of just copying to ensure compatibility
        const csvContent = fs.readFileSync(originalFullCsv, 'utf-8');
        const dfFull = parseCSV(csvContent, { columns: true, skip_empty_lines: true });
        const csvOutput = stringifyCSV(dfFull, { header: true });
        fs.writeFileSync(rawFullCsv, csvOutput);
        console.log(`Copied full results to ${rawFullCsv}`);
    }
    if (fs.existsSync(originalSummaryCsv)) {
        const csvContent = fs.readFileSync(originalSummaryCsv, 'utf-8');
        const dfSummary = parseCSV(csvContent, { columns: true, skip_empty_lines: true });
        const csvOutput = stringifyCSV(dfSummary, { header: true });
        fs.writeFileSync(rawSummaryCsv, csvOutput);
        console.log(`Copied summary results to ${rawSummaryCsv}`);
    }

    // Merge annotations if requested
    if (mergeAnnotations) {
        console.log("\n=== Starting annotation merging ===");
        
        try {
            // Sort the CSV file by True Cell Type before merging to ensure consistent order
            console.log("Sorting CSV by True Cell Type before merging...");
            const csvContent = fs.readFileSync(rawFullCsv, 'utf-8');
            let df = parseCSV(csvContent, { columns: true, skip_empty_lines: true });
            
            // Sort by True Cell Type
            df = df.sort((a, b) => {
                const aVal = parseInt(a['True Cell Type']) || 0;
                const bVal = parseInt(b['True Cell Type']) || 0;
                return aVal - bVal;
            });
            
            const sortedCsvOutput = stringifyCSV(df, { header: true });
            fs.writeFileSync(rawSortedCsv, sortedCsvOutput);
            
            // Run the merging process on the sorted CSV
            await mergeAnnotationsAll({
                csvPath: rawSortedCsv,
                outputPath: mergedAnnotationFile,
                provider: mergeProvider,
                model: mergeModel,
                additionalContext: `These are cell clusters from ${species} ${tissue}. ${additionalInfo}`
            });
            console.log(`✓ Annotations merged and saved to ${mergedAnnotationFile}`);
        } catch (error) {
            console.log(`! Error during annotation merging: ${error.message}`);
        }
    }
    
    console.log("\n=== Starting scoring process ===");
    // Run scoring
    try {
        if (fs.existsSync(rawFullCsv)) {
            await runCASSIAScoreBatch({
                inputFile: rawFullCsv,
                outputFile: scoreFileName,
                maxWorkers: maxWorkers,
                model: scoreModel,
                provider: scoreProvider,
                maxRetries: maxRetries
            });
            console.log("✓ Scoring process completed");
        } else {
            console.log(`⚠️  Input file for scoring not found: ${rawFullCsv}`);
            // Create empty scored file to prevent further errors
            const emptyCsv = 'True Cell Type,Predicted Main Cell Type,Score\n';
            fs.writeFileSync(scoreFileName, emptyCsv);
        }
    } catch (error) {
        console.log(`⚠️  Scoring process failed: ${error.message}`);
        // Create empty scored file to prevent further errors
        const emptyCsv = 'True Cell Type,Predicted Main Cell Type,Score\n';
        fs.writeFileSync(scoreFileName, emptyCsv);
    }

    console.log("\n=== Creating final combined results ===");
    let finalCombinedFile;
    // Create final combined CSV with all results
    try {
        // Read the scored file (which has all the original data plus scores)
        const scoredContent = fs.readFileSync(scoreFileName, 'utf-8');
        let finalDf = parseCSV(scoredContent, { columns: true, skip_empty_lines: true });
        
        // If merged annotations exist, add merged columns
        if (fs.existsSync(mergedAnnotationFile)) {
            const mergedContent = fs.readFileSync(mergedAnnotationFile, 'utf-8');
            const mergedDf = parseCSV(mergedContent, { columns: true, skip_empty_lines: true });
            
            // Merge on 'True Cell Type' to add merged annotation columns
            if (mergedDf.length > 0 && 'True Cell Type' in mergedDf[0]) {
                // Create a map for merging
                const mergedMap = new Map();
                mergedDf.forEach(row => {
                    mergedMap.set(row['True Cell Type'], row);
                });
                
                // Add merged columns to final data
                finalDf = finalDf.map(row => {
                    const clusterType = row['True Cell Type'];
                    const mergedRow = mergedMap.get(clusterType);
                    if (mergedRow) {
                        // Add only the merged columns (not duplicating existing ones)
                        const mergeColumns = Object.keys(mergedRow).filter(col => 
                            !Object.keys(row).includes(col) || col === 'True Cell Type'
                        );
                        mergeColumns.forEach(col => {
                            if (col !== 'True Cell Type') {
                                row[col] = mergedRow[col];
                            }
                        });
                    }
                    return row;
                });
            }
        }
        
        // Sort the final results by True Cell Type
        finalDf = finalDf.sort((a, b) => {
            const aVal = parseInt(a['True Cell Type']) || 0;
            const bVal = parseInt(b['True Cell Type']) || 0;
            return aVal - bVal;
        });
        
        // Save the final combined results
        finalCombinedFile = path.join(annotationResultsFolder, `${outputFileName}_FINAL_RESULTS.csv`);
        const finalCsvOutput = stringifyCSV(finalDf, { header: true });
        fs.writeFileSync(finalCombinedFile, finalCsvOutput);
        console.log(`✓ Final combined results saved to ${finalCombinedFile}`);
        
    } catch (error) {
        console.log(`Warning: Could not create final combined results: ${error.message}`);
        finalCombinedFile = scoreFileName;  // Fallback to scored file
    }

    console.log("\n=== Generating main reports ===");
    // Process reports - ensure they go to reports_folder
    try {
        if (fs.existsSync(scoreFileName)) {
            await runCASSIAGenerateScoreReport(scoreFileName, reportBaseName);
        } else {
            console.log(`⚠️  Score file not found, skipping report generation: ${scoreFileName}`);
        }
    } catch (error) {
        console.log(`⚠️  Report generation failed: ${error.message}`);
    }
    
    // Move any HTML files from annotation_results_folder to reports_folder
    const annotationFiles = fs.readdirSync(annotationResultsFolder);
    for (const file of annotationFiles) {
        if (file.endsWith('.html')) {
            const srcPath = path.join(annotationResultsFolder, file);
            const dstPath = path.join(reportsFolder, file);
            try {
                fs.copyFileSync(srcPath, dstPath);
                fs.unlinkSync(srcPath);  // Remove from original location after copying
                console.log(`Moved HTML report ${file} to reports folder`);
            } catch (error) {
                console.log(`Error moving HTML file ${file}: ${error.message}`);
            }
        }
    }
    
    console.log("✓ Main reports generated");

    console.log("\n=== Analyzing low-scoring clusters ===");
    // Handle low-scoring clusters
    let lowScoreClusters = [];
    
    try {
        if (fs.existsSync(scoreFileName)) {
            const scoredContent = fs.readFileSync(scoreFileName, 'utf-8');
            const df = parseCSV(scoredContent, { columns: true, skip_empty_lines: true });
            lowScoreClusters = df
                .filter(row => row.Score && parseFloat(row.Score) < scoreThreshold)
                .map(row => row['True Cell Type']);
        } else {
            console.log(`⚠️  Score file not found: ${scoreFileName}`);
        }
    } catch (error) {
        console.log(`⚠️  Error reading score file: ${error.message}`);
    }

    console.log(`Found ${lowScoreClusters.length} clusters with scores below ${scoreThreshold}:`);
    console.log(lowScoreClusters);
    
    if (lowScoreClusters.length > 0) {
        console.log("\n=== Starting boost annotation for low-scoring clusters ===");
        
        // Create boosted reports list - we will NOT generate a combined report
        for (const cluster of lowScoreClusters) {
            console.log(`Processing low score cluster: ${cluster}`);
            
            // Keep the original cluster name for data lookup
            const originalClusterName = cluster;
            
            // Sanitize the cluster name only for file naming purposes
            const sanitizedClusterName = String(cluster).replace(/[^a-zA-Z0-9\s\-_]/g, '').trim();
            
            // Create individual folder for this cluster's boost analysis
            const clusterBoostFolder = path.join(boostFolder, sanitizedClusterName);
            if (!fs.existsSync(clusterBoostFolder)) {
                fs.mkdirSync(clusterBoostFolder, { recursive: true });
            }
                
            // Define output name for the cluster boost report
            const clusterOutputName = path.join(clusterBoostFolder, `${outputFileName}_${sanitizedClusterName}_boosted`);
            
            // Use the original name for data lookup
            try {
                // major_cluster_info should be simple user-provided information like "human large intestine"
                // NOT complex data extracted from CSV
                const majorClusterInfo = `${species} ${tissue}`;
                
                // Run annotation boost - use original cluster name for data lookup, but sanitized name for output file
                // NOTE: Using the rawFullCsv path to ensure the CSV can be found
                await runCASSIAAnnotationboost({
                    fullResultPath: rawFullCsv,  // This is in the annotation_results_folder
                    marker: marker,
                    clusterName: originalClusterName,
                    majorClusterInfo: majorClusterInfo,
                    outputName: clusterOutputName,
                    numIterations: 5,
                    model: annotationboostModel,
                    provider: annotationboostProvider,
                    temperature: 0,
                    conversationHistoryMode: conversationHistoryMode,
                    reportStyle: reportStyle
                });
            } catch (error) {
                if (error.message.includes('No data found')) {
                    console.log(`Error in pipeline: No data found for cluster: ${originalClusterName}`);
                } else {
                    console.log(`Error in pipeline processing cluster ${originalClusterName}: ${error.message}`);
                }
            }
        }
        
        console.log("✓ Boost annotation completed");
    }
    
    console.log("\n=== Organizing intermediate files ===");
    // Create intermediate files folder
    const intermediateFolder = path.join(annotationResultsFolder, "intermediate_files");
    if (!fs.existsSync(intermediateFolder)) {
        fs.mkdirSync(intermediateFolder, { recursive: true });
    }
    
    // List of intermediate files to move
    const intermediateFiles = [
        rawFullCsv,
        rawSummaryCsv, 
        rawSortedCsv,
        scoreFileName,
        mergedAnnotationFile
    ];
    
    // Move intermediate files to intermediate folder
    for (const filePath of intermediateFiles) {
        if (fs.existsSync(filePath)) {
            try {
                const filename = path.basename(filePath);
                const destination = path.join(intermediateFolder, filename);
                fs.renameSync(filePath, destination);
                console.log(`Moved ${filename} to intermediate_files folder`);
            } catch (error) {
                console.log(`Warning: Could not move ${path.basename(filePath)}: ${error.message}`);
            }
        }
    }
    
    console.log("✓ Intermediate files organized");
    
    // Try to clean up the original files in the root directory
    try {
        const filesToRemove = [originalFullCsv, originalSummaryCsv, `${annotationOutput}_sorted_full.csv`];
        for (const fileToRemove of filesToRemove) {
            if (fs.existsSync(fileToRemove)) {
                fs.unlinkSync(fileToRemove);
                console.log(`Removed original file: ${fileToRemove}`);
            }
        }
    } catch (error) {
        console.log(`Warning: Could not remove some temporary files: ${error.message}`);
    }
    
    console.log("\n=== Cell type analysis pipeline completed ===");
    console.log(`All results have been organized in the '${mainFolderName}' folder:`);
    console.log(`  📊 MAIN RESULTS: ${finalCombinedFile}`);
    console.log(`  📁 HTML Reports: ${reportsFolder}`);
    console.log(`  🔍 Annotation Boost Results: ${boostFolder}`);
    console.log(`  📂 Intermediate Files: ${intermediateFolder}`);
    console.log(`\n✅ Your final results are in: ${path.basename(finalCombinedFile)}`);

    return finalCombinedFile;
}

/**
 * Load built-in marker files.
 * 
 * @param {string} markerType - Type of markers to load. Options:
 *   - "processed": For processed marker data
 *   - "unprocessed": For raw unprocessed marker data
 *   - "subcluster_results": For subcluster analysis results
 * 
 * @returns {Promise<Array>} Marker data as array of objects
 * 
 * @throws {Error} If marker_type is not recognized
 */
export async function loadMarker(markerType = "processed") {
    const markerFiles = {
        "processed": "processed.csv",
        "unprocessed": "unprocessed.csv",
        "subcluster_results": "subcluster_results.csv"
    };
    
    if (!markerFiles[markerType]) {
        throw new Error(`Unknown marker type: ${markerType}. Available types: ${Object.keys(markerFiles).join(', ')}`);
    }
    
    const filename = markerFiles[markerType];
    
    try {
        // For JavaScript implementation, we'll look for the data files in a relative path
        const dataPath = path.join(__dirname, '..', '..', 'CASSIA', 'data', filename);
        
        if (!fs.existsSync(dataPath)) {
            throw new Error(`Marker file not found: ${dataPath}`);
        }
        
        const csvContent = fs.readFileSync(dataPath, 'utf-8');
        const data = parseCSV(csvContent, { columns: true, skip_empty_lines: true });
        
        console.log(`Loaded ${data.length} rows from ${filename}`);
        return data;
        
    } catch (error) {
        throw new Error(`Error loading marker file ${filename}: ${error.message}`);
    }
}