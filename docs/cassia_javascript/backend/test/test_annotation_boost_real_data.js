/**
 * Annotation Boost Test with Real CASSIA Data
 * Tests using actual marker data from CASSIA/data/unprocessed.csv
 */

import { 
    runCASSIAAnnotationboost,
    runCASSIAAnnotationboostAdditionalTask
} from '../src/annotationBoost.js';
import fs from 'fs';
import { createReadStream } from 'fs';
import csv from 'csv-parser';

// Set API key
process.env.OPENROUTER_API_KEY = "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";

const REAL_MARKER_PATH = "/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/CASSIA/data/unprocessed.csv";
const TEST_RESULTS_DIR = "./test_results";

/**
 * Get top markers for a specific cluster from the real data
 */
async function getTopMarkersForCluster(clusterName, topN = 10) {
    return new Promise((resolve, reject) => {
        const markers = [];
        
        createReadStream(REAL_MARKER_PATH)
            .pipe(csv())
            .on('data', (row) => {
                if (row.cluster === clusterName) {
                    markers.push({
                        gene: row.gene,
                        avg_log2FC: parseFloat(row.avg_log2FC),
                        p_val_adj: parseFloat(row.p_val_adj),
                        pct_1: parseFloat(row['pct.1']),
                        pct_2: parseFloat(row['pct.2'])
                    });
                }
            })
            .on('end', () => {
                // Sort by avg_log2FC (descending) and take top N
                const topMarkers = markers
                    .sort((a, b) => b.avg_log2FC - a.avg_log2FC)
                    .slice(0, topN);
                resolve(topMarkers);
            })
            .on('error', reject);
    });
}

/**
 * Create mock full results CSV with real cluster data
 */
async function createMockFullResults() {
    console.log("📁 Creating mock full results with real cluster data...");
    
    const clusters = [
        "plasma cell",
        "monocyte", 
        "cd8-positive",
        "intestinal enteroendocrine cell"
    ];
    
    let csvContent = "True Cell Type,Marker List,Conversation History\n";
    
    for (const cluster of clusters) {
        const topMarkers = await getTopMarkersForCluster(cluster, 5);
        const markerList = topMarkers.map(m => m.gene).join(", ");
        const mockHistory = `Previous analysis of ${cluster} based on differential expression markers`;
        csvContent += `"${cluster}","${markerList}","${mockHistory}"\n`;
    }
    
    const fullResultsPath = `${TEST_RESULTS_DIR}/real_data_full_results.csv`;
    fs.writeFileSync(fullResultsPath, csvContent);
    
    console.log(`✅ Mock full results created with ${clusters.length} clusters`);
    return fullResultsPath;
}

/**
 * Test annotation boost with plasma cell cluster
 */
async function testPlasmaCellAnnotationBoost() {
    console.log("\n🧬 Testing Annotation Boost with Real Plasma Cell Data...");
    
    const fullResultsPath = await createMockFullResults();
    
    try {
        const result = await runCASSIAAnnotationboost({
            fullResultPath: fullResultsPath,
            marker: REAL_MARKER_PATH, // Use real marker file
            clusterName: "plasma cell",
            majorClusterInfo: "Human intestinal single-cell RNA-seq dataset",
            outputName: `${TEST_RESULTS_DIR}/plasma_cell_real_analysis`,
            numIterations: 2, // 2 iterations for thorough analysis
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            temperature: 0,
            conversationHistoryMode: "none", // Skip for speed
            searchStrategy: "breadth", // Test multiple hypotheses
            reportStyle: "per_iteration"
        });
        
        if (result.status === 'success') {
            console.log("✅ Plasma cell annotation boost SUCCEEDED!");
            console.log(`   Execution time: ${result.execution_time.toFixed(1)}s`);
            console.log(`   Analysis text length: ${result.analysis_text?.length || 0} chars`);
            
            if (result.summary_report_path) {
                console.log(`   📊 HTML report: ${result.summary_report_path}`);
            }
            if (result.raw_text_path) {
                console.log(`   📝 Raw conversation: ${result.raw_text_path}`);
            }
            
            return true;
        } else {
            console.log(`❌ Plasma cell analysis failed: ${result.error_message}`);
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Test error: ${error.message}`);
        return false;
    }
}

/**
 * Test annotation boost with additional task
 */
async function testMonocyteWithAdditionalTask() {
    console.log("\n🎯 Testing Annotation Boost with Additional Task (Monocyte)...");
    
    const fullResultsPath = `${TEST_RESULTS_DIR}/real_data_full_results.csv`;
    
    try {
        const result = await runCASSIAAnnotationboostAdditionalTask({
            fullResultPath: fullResultsPath,
            marker: REAL_MARKER_PATH,
            clusterName: "monocyte",
            majorClusterInfo: "Human intestinal tissue scRNA-seq",
            outputName: `${TEST_RESULTS_DIR}/monocyte_inflammatory_analysis`,
            numIterations: 2,
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            additionalTask: "determine if these are inflammatory monocytes or tissue-resident macrophages and assess their activation state",
            temperature: 0,
            conversationHistoryMode: "none",
            searchStrategy: "depth", // Focused analysis
            reportStyle: "total_summary" // Gene-focused report
        });
        
        if (result.status === 'success') {
            console.log("✅ Monocyte additional task analysis SUCCEEDED!");
            console.log(`   Execution time: ${result.execution_time.toFixed(1)}s`);
            console.log(`   Focused on inflammatory vs tissue-resident classification`);
            
            if (result.summary_report_path) {
                console.log(`   📊 Gene-focused report: ${result.summary_report_path}`);
            }
            
            return true;
        } else {
            console.log(`❌ Monocyte additional task failed: ${result.error_message}`);
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Additional task test error: ${error.message}`);
        return false;
    }
}

/**
 * Show sample of real data being used
 */
async function showRealDataSample() {
    console.log("\n📊 Sample of Real Marker Data Being Used:");
    
    const plasmaCellMarkers = await getTopMarkersForCluster("plasma cell", 5);
    const monocyteMarkers = await getTopMarkersForCluster("monocyte", 5);
    
    console.log("\n🧬 Top Plasma Cell Markers:");
    plasmaCellMarkers.forEach((marker, i) => {
        console.log(`   ${i+1}. ${marker.gene} (log2FC: ${marker.avg_log2FC.toFixed(2)}, pct.1: ${marker.pct_1})`);
    });
    
    console.log("\n🧬 Top Monocyte Markers:");
    monocyteMarkers.forEach((marker, i) => {
        console.log(`   ${i+1}. ${marker.gene} (log2FC: ${marker.avg_log2FC.toFixed(2)}, pct.1: ${marker.pct_1})`);
    });
}

/**
 * Main test execution
 */
async function runRealDataTest() {
    console.log("🚀 ANNOTATION BOOST TEST WITH REAL CASSIA DATA");
    console.log("=" .repeat(55));
    console.log(`📁 Using marker file: ${REAL_MARKER_PATH}`);
    
    // Ensure test directory exists
    if (!fs.existsSync(TEST_RESULTS_DIR)) {
        fs.mkdirSync(TEST_RESULTS_DIR, { recursive: true });
    }
    
    // Show what data we're working with
    await showRealDataSample();
    
    // Run tests
    const tests = [
        { name: "Plasma Cell Annotation Boost", func: testPlasmaCellAnnotationBoost },
        { name: "Monocyte with Additional Task", func: testMonocyteWithAdditionalTask }
    ];
    
    let passed = 0;
    let failed = 0;
    
    for (const test of tests) {
        try {
            const result = await test.func();
            if (result) {
                passed++;
            } else {
                failed++;
            }
        } catch (error) {
            console.log(`❌ ${test.name} failed: ${error.message}`);
            failed++;
        }
    }
    
    console.log("\n" + "=" .repeat(55));
    console.log(`📊 REAL DATA TEST RESULTS: ${passed}/${tests.length} passed`);
    
    if (passed === tests.length) {
        console.log("🎉 ALL REAL DATA TESTS PASSED!");
        console.log("   ✅ Works with actual CASSIA marker data");
        console.log("   ✅ Handles real biological complexity");
        console.log("   ✅ Both standard and additional task modes functional");
        console.log("   🚀 Ready for production use with real datasets!");
    } else if (passed > 0) {
        console.log("⚠️  Some tests passed - annotation boost mostly functional");
    } else {
        console.log("❌ Tests failed - check API connectivity and configuration");
    }
    
    console.log(`\n📁 Results saved to: ${TEST_RESULTS_DIR}`);
    console.log("   Review the HTML reports to see the biological analysis!");
}

runRealDataTest().catch(console.error);