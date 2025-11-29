// CASSIA JavaScript Package - Complete Implementation
// 100% Replication of Python CASSIA Package

// Main Pipeline Function
export { runCASSIAPipeline, loadMarker } from './src/pipeline.js';

// Individual Agent Functions
export { runCASSIAAnnotationboost } from './src/annotationBoost.js';
export { mergeAnnotations, mergeAnnotationsAll } from './src/mergingAnnotation.js';
export { symphonyCompare } from './src/symphonyCompare.js';
export { runCASSIASubclusters, runCASSIANSubcluster } from './src/subclustering.js';

// Core Processing Functions  
export { runCASSIA } from './src/runCASSIA.js';
export { runCASSIABatch } from './src/runCASSIA_batch.js';
export { runCASSIAScoreBatch, scoreSingleAnalysis, scoreAnnotationBatch } from './src/scoring.js';

// Report Generation Functions
export { runCASSIAGenerateScoreReport } from './src/report_generation.js';
export { generateSubclusteringReport, processEvaluationCsv, createIndexHtml } from './src/generate_reports.js';

// LLM Utilities
export { callLLM } from './src/llm_utils.js';

// Main Function Code (core annotation logic)
export * from './src/main_function_code.js';

// Version and Package Info
export const CASSIA_VERSION = "2.0.0-js";
export const PACKAGE_NAME = "cassia-javascript";

/**
 * CASSIA (Cell Annotation using Single-cell Sequencing Intelligence Assistant)
 * 
 * This is a complete JavaScript replication of the Python CASSIA package for
 * automated cell type annotation in single-cell RNA sequencing data.
 * 
 * Key Features:
 * - runCASSIAPipeline: Complete end-to-end analysis pipeline
 * - Annotation Boost: Iterative marker analysis for challenging clusters
 * - Symphony Compare: Multi-model consensus cell type comparison
 * - Merging Annotations: LLM-powered annotation grouping and standardization
 * - Subclustering: Automated subcluster identification and annotation  
 * - Comprehensive Scoring: Quality assessment of annotations
 * - Rich Reporting: Interactive HTML reports with visualizations
 * 
 * Usage Example:
 * ```javascript
 * import { runCASSIAPipeline } from 'cassia-javascript';
 * 
 * const results = await runCASSIAPipeline({
 *   outputFileName: 'my_analysis',
 *   tissue: 'peripheral_blood', 
 *   species: 'human',
 *   marker: markerData,
 *   annotationModel: 'claude-3-5-sonnet-20241022',
 *   annotationProvider: 'anthropic'
 * });
 * ```
 * 
 * Supported LLM Providers:
 * - OpenAI (GPT-4, GPT-3.5-turbo)
 * - Anthropic (Claude models)  
 * - OpenRouter (Multiple models)
 * - Custom endpoints
 * 
 * For detailed documentation and examples, see the README.md file.
 */

// Create default export object
const CASSIA = {
    runCASSIAPipeline,
    loadMarker,
    runCASSIAAnnotationboost,
    mergeAnnotations,
    mergeAnnotationsAll,
    symphonyCompare,
    runCASSIASubclusters,
    runCASSIANSubcluster,
    runCASSIA,
    runCASSIABatch,
    runCASSIAScoreBatch,
    scoreSingleAnalysis,
    scoreAnnotationBatch,
    runCASSIAGenerateScoreReport,
    generateSubclusteringReport,
    processEvaluationCsv,
    createIndexHtml,
    callLLM,
    CASSIA_VERSION,
    PACKAGE_NAME
};

console.log(`
🧬 CASSIA JavaScript Package v${CASSIA_VERSION}
   Cell Annotation using Single-cell Sequencing Intelligence Assistant
   
   Available Functions:
   📊 runCASSIAPipeline()     - Complete analysis pipeline
   🚀 runCASSIAAnnotationboost() - Iterative marker analysis  
   🎼 symphonyCompare()       - Multi-model consensus analysis
   🔗 mergeAnnotationsAll()   - Annotation grouping & standardization
   🧩 runCASSIASubclusters()  - Subcluster identification
   📈 runCASSIAScoreBatch()   - Annotation quality scoring
   📋 runCASSIAGenerateScoreReport() - Interactive HTML reports
   
   Ready for single-cell analysis! 🎯
`);

export default CASSIA;