/**
 * CASSIA Annotation Boost Example
 * Demonstrates iterative marker analysis for enhanced cell type annotation
 */

import { runCASSIAAnnotationboost, runCASSIAAnnotationboostAdditionalTask } from '../index.js';
import path from 'path';

// Configuration
const API_KEY = "your-api-key-here";
const OUTPUT_DIR = "./annotation_boost_results";

/**
 * Example 1: Basic Annotation Boost
 * Performs iterative marker analysis to refine cell type annotation
 */
async function basicAnnotationBoostExample() {
    console.log("🔬 Running Basic Annotation Boost Example...");
    
    try {
        const result = await runCASSIAAnnotationboost({
            fullResultPath: "../test/test_results/mock_full_results.csv",
            marker: "../test/test_results/mock_markers.csv",
            clusterName: "T cell",
            majorClusterInfo: "Human lung tissue scRNA-seq dataset",
            outputName: path.join(OUTPUT_DIR, "t_cell_annotation_boost"),
            numIterations: 3,
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            temperature: 0,
            conversationHistoryMode: "final", // Uses summarization
            searchStrategy: "breadth", // Test multiple hypotheses
            reportStyle: "per_iteration"
        });
        
        if (result.status === 'success') {
            console.log("✅ Annotation boost completed successfully!");
            console.log(`   Execution time: ${result.execution_time}s`);
            console.log(`   Raw conversation: ${result.raw_text_path}`);
            console.log(`   Summary report: ${result.summary_report_path}`);
        } else {
            console.log(`❌ Annotation boost failed: ${result.error_message}`);
        }
        
    } catch (error) {
        console.log(`❌ Error: ${error.message}`);
    }
}

/**
 * Example 2: Depth-First Analysis with Additional Task
 * Uses focused analysis approach with a specific biological question
 */
async function depthFirstAnalysisExample() {
    console.log("\n🎯 Running Depth-First Analysis with Additional Task...");
    
    try {
        const result = await runCASSIAAnnotationboostAdditionalTask({
            fullResultPath: "../test/test_results/mock_full_results.csv",
            marker: "../test/test_results/mock_markers.csv",
            clusterName: "T cell",
            majorClusterInfo: "Human lung tissue scRNA-seq from cancer patients",
            outputName: path.join(OUTPUT_DIR, "t_cell_regulatory_analysis"),
            numIterations: 2,
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            additionalTask: "determine if these are regulatory T cells (Tregs) and assess their activation state",
            temperature: 0,
            conversationHistoryMode: "none", // Skip history for speed
            searchStrategy: "depth", // Focus on one hypothesis at a time
            reportStyle: "total_summary" // Gene-focused report
        });
        
        if (result.status === 'success') {
            console.log("✅ Depth-first analysis completed successfully!");
            console.log(`   Execution time: ${result.execution_time}s`);
            console.log(`   Analysis focused on regulatory T cell identification`);
            console.log(`   Summary report: ${result.summary_report_path}`);
        } else {
            console.log(`❌ Analysis failed: ${result.error_message}`);
        }
        
    } catch (error) {
        console.log(`❌ Error: ${error.message}`);
    }
}

/**
 * Example 3: Custom Marker Data Analysis
 * Shows how to use custom marker data instead of CSV files
 */
async function customMarkerDataExample() {
    console.log("\n🧬 Running Custom Marker Data Analysis...");
    
    // Custom marker data as JavaScript array
    const customMarkerData = [
        {
            gene: "CD3D",
            avg_log2FC: 2.5,
            p_val_adj: 1.2e-15,
            pct_1: 0.85,
            pct_2: 0.12
        },
        {
            gene: "CD3E",
            avg_log2FC: 2.1,
            p_val_adj: 3.4e-12,
            pct_1: 0.82,
            pct_2: 0.15
        },
        {
            gene: "CD4",
            avg_log2FC: 1.2,
            p_val_adj: 2.3e-05,
            pct_1: 0.65,
            pct_2: 0.35
        },
        {
            gene: "FOXP3",
            avg_log2FC: 0.3,
            p_val_adj: 0.12,
            pct_1: 0.15,
            pct_2: 0.85
        },
        {
            gene: "IL2RA",
            avg_log2FC: 1.8,
            p_val_adj: 4.5e-08,
            pct_1: 0.45,
            pct_2: 0.25
        }
    ];
    
    try {
        const result = await runCASSIAAnnotationboost({
            fullResultPath: "../test/test_results/mock_full_results.csv",
            marker: customMarkerData, // Using custom data instead of CSV
            clusterName: "T cell", 
            majorClusterInfo: "Human PBMC scRNA-seq data",
            outputName: path.join(OUTPUT_DIR, "custom_marker_analysis"),
            numIterations: 2,
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            temperature: 0.1, // Slightly higher for creativity
            conversationHistoryMode: "full", // Use full history
            searchStrategy: "breadth",
            reportStyle: "per_iteration"
        });
        
        if (result.status === 'success') {
            console.log("✅ Custom marker analysis completed!");
            console.log(`   Used ${customMarkerData.length} custom marker genes`);
            console.log(`   Generated reports in: ${OUTPUT_DIR}`);
        } else {
            console.log(`❌ Custom analysis failed: ${result.error_message}`);
        }
        
    } catch (error) {
        console.log(`❌ Error: ${error.message}`);
    }
}

/**
 * Main execution function
 */
async function runExamples() {
    console.log("🚀 CASSIA Annotation Boost Examples");
    console.log("=" .repeat(50));
    
    // Note: You would need to set your actual API key here
    console.log("⚠️  Note: Set your API key in the configuration section before running");
    
    // Create output directory if it doesn't exist
    try {
        const fs = await import('fs');
        if (!fs.existsSync(OUTPUT_DIR)) {
            fs.mkdirSync(OUTPUT_DIR, { recursive: true });
        }
    } catch (error) {
        console.log("Could not create output directory");
    }
    
    // Run examples (uncomment to test with real API)
    // await basicAnnotationBoostExample();
    // await depthFirstAnalysisExample();
    // await customMarkerDataExample();
    
    console.log("\n📚 Example Usage Summary:");
    console.log("1. Basic Annotation Boost: Iterative analysis with breadth-first search");
    console.log("2. Depth-First + Additional Task: Focused analysis with specific questions");
    console.log("3. Custom Marker Data: Use JavaScript arrays instead of CSV files");
    console.log("\n🔧 Key Parameters:");
    console.log("- searchStrategy: 'breadth' (multiple hypotheses) or 'depth' (focused)");
    console.log("- reportStyle: 'per_iteration' or 'total_summary'");
    console.log("- conversationHistoryMode: 'full', 'final' (summarized), or 'none'");
    console.log("- numIterations: Number of analysis rounds (1-10 recommended)");
}

// Run examples
runExamples().catch(console.error);

// Export for use in other modules
export {
    basicAnnotationBoostExample,
    depthFirstAnalysisExample,
    customMarkerDataExample
};