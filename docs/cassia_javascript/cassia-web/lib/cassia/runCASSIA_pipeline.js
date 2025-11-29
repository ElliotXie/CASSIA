/**
 * CASSIA Pipeline Orchestrator - Complete End-to-End Analysis
 * 
 * This module orchestrates the complete CASSIA pipeline:
 * 1. Initial annotation (runCASSIA_batch)
 * 2. Quality scoring (scoreAnnotationBatch)
 * 3. Annotation boost for low-scoring clusters (iterativeMarkerAnalysis)
 * 4. Results organization and download preparation
 */

import { runCASSIABatch } from './runCASSIA_batch.js';
import { scoreAnnotationBatch } from './scoring.js';
import { iterativeMarkerAnalysis, generateSummaryReport } from './annotationBoost.js';
import { parseCSV, formatAsCSV } from '../utils/csv-parser.js';

/**
 * Pipeline execution state for progress tracking
 */
class PipelineState {
    constructor() {
        this.currentStep = '';
        this.totalSteps = 0;
        this.completedSteps = 0;
        this.currentProgress = 0;
        this.logs = [];
        this.results = {};
        this.lowScoreClusters = [];
        this.isRunning = false;
        this.error = null;
    }

    updateStep(step, progress = 0) {
        this.currentStep = step;
        this.currentProgress = progress;
        this.addLog(`🔄 ${step} - ${progress}%`);
    }

    completeStep(step) {
        this.completedSteps++;
        this.currentProgress = Math.round((this.completedSteps / this.totalSteps) * 100);
        this.addLog(`✅ ${step} completed`);
    }

    addLog(message) {
        this.logs.push({
            timestamp: new Date().toISOString(),
            message
        });
        console.log(message);
    }

    setError(error) {
        this.error = error;
        this.addLog(`❌ Error: ${error.message}`);
        this.isRunning = false;
    }
}

/**
 * Main pipeline orchestrator function
 * 
 * @param {Object} config - Pipeline configuration
 * @param {Array|string} config.marker - Marker data (array or file path)
 * @param {string} config.outputName - Base name for output files
 * @param {string} config.tissue - Tissue type
 * @param {string} config.species - Species
 * @param {Object} config.models - Model configuration for each step
 * @param {Object} config.models.annotation - {provider, model} for annotation
 * @param {Object} config.models.scoring - {provider, model} for scoring
 * @param {Object} config.models.annotationBoost - {provider, model} for boost
 * @param {number} config.scoreThreshold - Threshold for low-scoring clusters (default: 75)
 * @param {number} config.maxWorkers - Maximum parallel workers (default: 4)
 * @param {number} config.maxRetries - Maximum retry attempts (default: 1)
 * @param {string} config.additionalInfo - Additional context information
 * @param {string} config.apiKey - API key for authentication
 * @param {Function} config.onProgress - Progress callback function
 * @param {Function} config.onLog - Log callback function
 * @param {Function} config.onComplete - Completion callback function
 * @param {Function} config.onError - Error callback function
 * 
 * @returns {Promise<Object>} Pipeline results and file organization
 */
export async function runCASSIAPipeline(config) {
    const {
        marker,
        outputName = 'cassia_analysis',
        tissue = 'peripheral_blood',
        species = 'human',
        models = {
            annotation: { provider: 'openrouter', model: 'meta-llama/llama-3.1-8b-instruct' },
            scoring: { provider: 'openrouter', model: 'google/gemini-2.5-flash' },
            annotationBoost: { provider: 'openrouter', model: 'anthropic/claude-3.5-sonnet' }
        },
        scoreThreshold = 75,
        maxWorkers = 4,
        maxRetries = 1,
        additionalInfo = '',
        apiKey,
        onProgress = null,
        onLog = null,
        onComplete = null,
        onError = null
    } = config;

    const state = new PipelineState();
    state.totalSteps = 5; // annotation, scoring, boost analysis, report generation, file organization
    state.isRunning = true;

    // Helper function to call progress callback
    const notifyProgress = () => {
        if (onProgress) {
            onProgress({
                step: state.currentStep,
                progress: state.currentProgress,
                totalSteps: state.totalSteps,
                completedSteps: state.completedSteps,
                logs: state.logs,
                results: state.results,
                lowScoreClusters: state.lowScoreClusters,
                isRunning: state.isRunning
            });
        }
    };

    // Helper function to call log callback
    const notifyLog = (message) => {
        state.addLog(message);
        if (onLog) onLog(message);
        notifyProgress();
    };

    try {
        notifyLog('🚀 Starting CASSIA Pipeline Analysis');
        notifyLog(`📋 Configuration: ${tissue} ${species}, Threshold: ${scoreThreshold}`);
        notifyLog(`🤖 Models: Annotation=${models.annotation.model}, Scoring=${models.scoring.model}, Boost=${models.annotationBoost.model}`);

        // Step 1: Initial Annotation
        state.updateStep('Running initial cell type annotation', 0);
        notifyProgress();

        const annotationResults = await runCASSIABatch({
            marker,
            apiKey,
            outputName,
            model: models.annotation.model,
            provider: models.annotation.provider,
            tissue,
            species,
            additionalInfo,
            maxWorkers,
            maxRetries,
            onLog: notifyLog // Pass logging callback to show progress in UI
        });

        state.results.annotation = annotationResults;
        state.completeStep('Initial annotation');
        notifyLog(`📊 Annotation completed: ${Object.keys(annotationResults.results).length} clusters analyzed`);
        notifyProgress();

        // Step 2: Quality Scoring
        state.updateStep('Calculating quality scores', 0);
        notifyProgress();

        // Convert annotation results to CSV format for scoring
        // Ensure config is properly set with pipeline values
        if (!annotationResults.config) {
            annotationResults.config = {};
        }
        annotationResults.config.tissue = tissue;
        annotationResults.config.species = species;
        annotationResults.config.model = models.annotation.model;
        annotationResults.config.provider = models.annotation.provider;
        annotationResults.config.additionalInfo = additionalInfo;
        
        const csvData = convertAnnotationResultsToCSV(annotationResults);
        
        const scoringResults = await scoreAnnotationBatch({
            csvData,
            apiKey,
            model: models.scoring.model,
            provider: models.scoring.provider,
            maxWorkers,
            maxRetries,
            onProgress: (progress) => {
                // Extract percentage from progress object
                const percentage = typeof progress === 'object' ? progress.percentage : progress;
                state.updateStep('Calculating quality scores', percentage);
                notifyProgress();
            },
            onLog: notifyLog // Pass logging callback to show progress in UI
        });

        state.results.scoring = scoringResults;
        state.completeStep('Quality scoring');
        
        // Debug scoring results format
        console.log(`🔍 Scoring results type: ${typeof scoringResults}, isArray: ${Array.isArray(scoringResults)}`);
        console.log(`🔍 Scoring results structure:`, Object.keys(scoringResults));
        
        // Extract the actual results array from the scoring response
        const actualScoringResults = scoringResults.results || scoringResults;
        console.log(`🔍 Actual scoring results array length: ${Array.isArray(actualScoringResults) ? actualScoringResults.length : 'not array'}`);
        
        // Identify low-scoring clusters
        state.lowScoreClusters = identifyLowScoreClusters(actualScoringResults, scoreThreshold);
        notifyLog(`🎯 Scoring completed: ${state.lowScoreClusters.length} clusters below threshold (${scoreThreshold})`);
        notifyProgress();

        // Step 3: Annotation Boost for Low-Scoring Clusters
        state.updateStep('Running annotation boost for low-scoring clusters', 0);
        notifyProgress();

        const boostResults = {};
        if (state.lowScoreClusters.length > 0) {
            notifyLog(`🚀 Starting boost analysis for ${state.lowScoreClusters.length} clusters`);
            
            for (let i = 0; i < state.lowScoreClusters.length; i++) {
                const cluster = state.lowScoreClusters[i];
                const progress = Math.round((i / state.lowScoreClusters.length) * 100);
                
                state.updateStep(`Boosting cluster ${cluster.cellType} (${i + 1}/${state.lowScoreClusters.length})`, progress);
                notifyProgress();

                try {
                    const boostResult = await runAnnotationBoostForCluster(
                        cluster,
                        marker,
                        annotationResults,
                        models.annotationBoost,
                        apiKey,
                        additionalInfo
                    );
                    
                    boostResults[cluster.cellType] = boostResult;
                    notifyLog(`✅ Boost completed for cluster ${cluster.cellType}`);
                } catch (error) {
                    notifyLog(`⚠️ Boost failed for cluster ${cluster.cellType}: ${error.message}`);
                    boostResults[cluster.cellType] = { error: error.message };
                }
            }
        } else {
            notifyLog('ℹ️ No clusters below threshold - skipping boost analysis');
        }

        state.results.boost = boostResults;
        state.completeStep('Annotation boost');
        notifyProgress();

        // Skip HTML report generation - only keep CSV files
        notifyLog('📋 Skipping HTML report generation - keeping CSV files only');

        // Step 4: Organize Final Results
        state.updateStep('Organizing final results', 0);
        notifyProgress();

        const finalResults = await organizeFinalResults(state.results, config);
        
        state.results.final = finalResults;
        state.completeStep('File organization');
        state.isRunning = false;
        
        notifyLog('🎉 Pipeline analysis completed successfully!');
        notifyLog(`📁 Results organized in: ${finalResults.summary.totalFiles} files`);
        notifyProgress();

        if (onComplete) {
            onComplete(finalResults);
        }

        return finalResults;

    } catch (error) {
        state.setError(error);
        notifyProgress();
        
        if (onError) {
            onError(error);
        }
        
        throw error;
    }
}

/**
 * Convert annotation results to CSV format for scoring
 */
function convertAnnotationResultsToCSV(annotationResults) {
    const csvRows = [];
    
    // Add debugging info
    console.log(`🔍 Converting annotation results for ${Object.keys(annotationResults.results).length} clusters`);
    
    for (const [cellType, details] of Object.entries(annotationResults.results)) {
        const analysisResult = details.analysis_result || {};
        const conversationHistory = details.conversation_history || [];
        
        // Add debugging for conversation history format
        console.log(`🔍 Processing ${cellType}: conversationHistory type = ${typeof conversationHistory}, isArray = ${Array.isArray(conversationHistory)}`);
        
        // Format conversation history with robust type checking
        let formattedHistory = '';
        
        if (Array.isArray(conversationHistory)) {
            // Handle array format (expected)
            formattedHistory = conversationHistory
                .map(entry => {
                    if (Array.isArray(entry) && entry.length >= 2) {
                        return `${entry[0]}: ${entry[1]}`;
                    } else if (typeof entry === 'object' && entry.role && entry.content) {
                        return `${entry.role}: ${entry.content}`;
                    } else if (typeof entry === 'string') {
                        return entry;
                    }
                    return String(entry);
                })
                .join(' | ');
        } else if (typeof conversationHistory === 'string') {
            // Handle string format (already formatted)
            formattedHistory = conversationHistory;
        } else if (conversationHistory && typeof conversationHistory === 'object') {
            // Handle object format
            formattedHistory = JSON.stringify(conversationHistory);
        } else {
            // Handle other types or null/undefined
            formattedHistory = conversationHistory ? String(conversationHistory) : 'N/A';
        }
        
        csvRows.push({
            'True Cell Type': cellType,
            'Predicted Main Cell Type': analysisResult.main_cell_type || '',
            'Predicted Sub Cell Types': (analysisResult.sub_cell_types || []).join(', '),
            'Possible Mixed Cell Types': (analysisResult.possible_mixed_cell_types || []).join(', '),
            'Marker Number': analysisResult.num_markers || 0,
            'Marker List': (analysisResult.marker_list || []).join(', '),
            'Iterations': analysisResult.iterations || 1,
            'Model': annotationResults.config?.model || 'unknown',
            'Provider': annotationResults.config?.provider || 'unknown',
            'Tissue': annotationResults.config?.tissue || 'unknown',
            'Species': annotationResults.config?.species || 'unknown',
            'Additional Info': annotationResults.config?.additionalInfo || 'N/A',
            'Conversation History': formattedHistory
        });
    }
    
    return csvRows;
}

/**
 * Identify clusters with scores below threshold
 */
function identifyLowScoreClusters(scoringResults, threshold) {
    const lowScoreClusters = [];
    
    // Add robust type checking for scoringResults
    console.log(`🔍 identifyLowScoreClusters: scoringResults type = ${typeof scoringResults}, isArray = ${Array.isArray(scoringResults)}`);
    
    let resultsArray = [];
    
    if (Array.isArray(scoringResults)) {
        resultsArray = scoringResults;
    } else if (scoringResults && typeof scoringResults === 'object') {
        // Handle case where scoringResults might be an object with results property
        if (scoringResults.results && Array.isArray(scoringResults.results)) {
            resultsArray = scoringResults.results;
        } else if (scoringResults.data && Array.isArray(scoringResults.data)) {
            resultsArray = scoringResults.data;
        } else {
            // Try to convert object values to array
            resultsArray = Object.values(scoringResults);
        }
    } else {
        console.error(`❌ scoringResults is not iterable:`, scoringResults);
        return lowScoreClusters; // Return empty array
    }
    
    console.log(`🔍 Processing ${resultsArray.length} scoring results`);
    
    for (const result of resultsArray) {
        if (result.Score && parseFloat(result.Score) < threshold) {
            lowScoreClusters.push({
                cellType: result['True Cell Type'],
                score: parseFloat(result.Score),
                predictedType: result['Predicted Main Cell Type'],
                markers: result['Marker List'],
                conversationHistory: result['Conversation History']
            });
        }
    }
    
    return lowScoreClusters.sort((a, b) => a.score - b.score); // Sort by score (lowest first)
}

/**
 * Run annotation boost for a specific cluster
 */
async function runAnnotationBoostForCluster(cluster, marker, annotationResults, model, apiKey, additionalInfo) {
    const majorClusterInfo = `${annotationResults.config?.species || 'unknown'} ${annotationResults.config?.tissue || 'unknown'}`;
    const commaSeparatedGenes = cluster.markers;
    const annotationHistory = cluster.conversationHistory;
    
    // Run iterative marker analysis
    const analysisResult = await iterativeMarkerAnalysis(
        majorClusterInfo,
        marker,
        commaSeparatedGenes,
        annotationHistory,
        5, // numIterations
        model.provider,
        model.model,
        null, // additionalTask
        0, // temperature
        'breadth', // searchStrategy
        apiKey
    );
    
    // Generate summary report
    const summaryReport = await generateSummaryReport(
        analysisResult.messages,
        'breadth', // searchStrategy
        'per_iteration', // reportStyle
        model.provider,
        model.model,
        apiKey
    );
    
    return {
        conversation: analysisResult.conversation,
        messages: analysisResult.messages,
        summaryReport,
        cluster: cluster.cellType,
        originalScore: cluster.score
    };
}

/**
 * Generate comprehensive analysis reports
 */
async function generateAnalysisReports(results, lowScoreClusters, config) {
    const reports = {};
    
    // Main scoring report
    if (results.scoring) {
        reports.scoring = {
            type: 'scoring_summary',
            title: 'Quality Scoring Summary',
            data: results.scoring,
            metadata: {
                totalClusters: results.scoring.results ? results.scoring.results.length : (Array.isArray(results.scoring) ? results.scoring.length : 0),
                lowScoreClusters: lowScoreClusters.length,
                averageScore: calculateAverageScore(results.scoring),
                threshold: config.scoreThreshold
            }
        };
    }
    
    // Annotation boost reports
    if (results.boost && Object.keys(results.boost).length > 0) {
        reports.boostSummary = {
            type: 'boost_summary',
            title: 'Annotation Boost Summary',
            data: results.boost,
            metadata: {
                clustersProcessed: Object.keys(results.boost).length,
                successfulBoosts: Object.values(results.boost).filter(r => !r.error).length
            }
        };
        
        // Individual boost reports
        for (const [clusterName, boostResult] of Object.entries(results.boost)) {
            if (!boostResult.error && boostResult.summaryReport) {
                reports[`boost_${clusterName}`] = {
                    type: 'boost_detail',
                    title: `Annotation Boost: ${clusterName}`,
                    data: boostResult.summaryReport,
                    cluster: clusterName,
                    originalScore: boostResult.originalScore
                };
            }
        }
    }
    
    return reports;
}

/**
 * Calculate total unique genes from annotation results
 */
function calculateTotalGenes(annotationResults) {
    const allGenes = new Set();
    
    if (annotationResults.results) {
        for (const details of Object.values(annotationResults.results)) {
            const analysisResult = details.analysis_result || {};
            const markerList = analysisResult.marker_list || [];
            markerList.forEach(gene => {
                if (gene && typeof gene === 'string') {
                    allGenes.add(gene.trim());
                }
            });
        }
    }
    
    return allGenes.size;
}

/**
 * Calculate average score from scoring results
 */
function calculateAverageScore(scoringResults) {
    // Handle the scoring results format (object with results property or direct array)
    let resultsArray = [];
    
    if (Array.isArray(scoringResults)) {
        resultsArray = scoringResults;
    } else if (scoringResults && scoringResults.results && Array.isArray(scoringResults.results)) {
        resultsArray = scoringResults.results;
    } else {
        console.warn('⚠️ calculateAverageScore: scoringResults is not in expected format');
        return 0;
    }
    
    const validScores = resultsArray
        .map(r => parseFloat(r.Score))
        .filter(score => !isNaN(score));
    
    if (validScores.length === 0) return 0;
    
    return Math.round((validScores.reduce((a, b) => a + b, 0) / validScores.length) * 100) / 100;
}

/**
 * Organize final results into downloadable format
 */
async function organizeFinalResults(results, config) {
    const finalResults = {
        summary: {
            outputName: config.outputName,
            tissue: config.tissue,
            species: config.species,
            timestamp: new Date().toISOString(),
            totalFiles: 0,
            totalClusters: 0,
            lowScoreClusters: results.boost ? Object.keys(results.boost).length : 0
        },
        files: {
            annotation: null,
            scoring: null,
            reports: {},
            boost: {}
        },
        downloadUrls: {},
        metadata: {
            models: config.models,
            scoreThreshold: config.scoreThreshold,
            averageScore: results.scoring ? calculateAverageScore(results.scoring) : 0,
            totalGenes: results.annotation ? calculateTotalGenes(results.annotation) : 0,
            configuration: {
                maxWorkers: config.maxWorkers,
                maxRetries: config.maxRetries,
                additionalInfo: config.additionalInfo
            }
        }
    };
    
    // Prepare annotation results file
    if (results.annotation) {
        const annotationCSV = formatAsCSV(convertAnnotationResultsToCSV(results.annotation));
        finalResults.files.annotation = {
            filename: `${config.outputName}_annotation_results.csv`,
            content: annotationCSV,
            type: 'text/csv',
            size: annotationCSV.length
        };
        finalResults.summary.totalFiles++;
        finalResults.summary.totalClusters = Object.keys(results.annotation.results).length;
    }
    
    // Prepare scoring results file
    if (results.scoring) {
        // Extract the actual results array from the scoring response
        const actualScoringData = results.scoring.results || results.scoring;
        const scoringCSV = formatAsCSV(actualScoringData);
        finalResults.files.scoring = {
            filename: `${config.outputName}_scoring_results.csv`,
            content: scoringCSV,
            type: 'text/csv',
            size: scoringCSV.length
        };
        // Store the raw scoring data for UI display (same data that creates perfect CSV)
        finalResults.rawScoringData = actualScoringData;
        finalResults.summary.totalFiles++;
    }
    
    // Skip HTML report generation - only keep CSV files
    
    // Prepare boost analysis files
    if (results.boost) {
        for (const [clusterName, boostResult] of Object.entries(results.boost)) {
            if (!boostResult.error) {
                // Conversation file
                const conversationFilename = `${config.outputName}_boost_${clusterName}_conversation.txt`;
                finalResults.files.boost[`${clusterName}_conversation`] = {
                    filename: conversationFilename,
                    content: boostResult.conversation || '',
                    type: 'text/plain',
                    size: (boostResult.conversation || '').length
                };
                
                // Summary report file
                if (boostResult.summaryReport) {
                    const reportFilename = `${config.outputName}_boost_${clusterName}_report.html`;
                    finalResults.files.boost[`${clusterName}_report`] = {
                        filename: reportFilename,
                        content: boostResult.summaryReport,
                        type: 'text/html',
                        size: boostResult.summaryReport.length
                    };
                }
                
                finalResults.totalFiles += 2;
            }
        }
    }
    
    // Create download URLs (would be implemented by the calling component)
    finalResults.downloadUrls = {
        individual: {}, // Individual file downloads
        archive: null   // Complete ZIP download
    };
    
    return finalResults;
}

/**
 * Create downloadable blob URLs for results
 * This function should be called by the frontend component
 */
export function createDownloadUrls(finalResults) {
    const urls = {
        individual: {},
        blobs: [] // Keep track of blobs for cleanup
    };
    
    // Create URLs for individual files
    const allFiles = {
        ...finalResults.files.reports,
        ...finalResults.files.boost
    };
    
    if (finalResults.files.annotation) {
        allFiles.annotation = finalResults.files.annotation;
    }
    
    if (finalResults.files.scoring) {
        allFiles.scoring = finalResults.files.scoring;
    }
    
    for (const [key, file] of Object.entries(allFiles)) {
        const blob = new Blob([file.content], { type: file.type });
        const url = URL.createObjectURL(blob);
        
        urls.individual[key] = {
            url,
            filename: file.filename,
            size: file.size,
            type: file.type
        };
        
        urls.blobs.push(blob);
    }
    
    return urls;
}

/**
 * Cleanup download URLs to prevent memory leaks
 */
export function cleanupDownloadUrls(urls) {
    if (urls && urls.individual) {
        for (const fileInfo of Object.values(urls.individual)) {
            if (fileInfo.url) {
                URL.revokeObjectURL(fileInfo.url);
            }
        }
    }
}