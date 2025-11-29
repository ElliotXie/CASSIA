/**
 * Quick API Test for Annotation Boost
 * Tests real LLM integration with minimal iterations
 */

import { 
    runCASSIAAnnotationboost,
    summarizeConversationHistory
} from '../src/annotationBoost.js';
import fs from 'fs';
import path from 'path';

const API_KEY = "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";
const TEST_RESULTS_DIR = "./test_results";

// Set API key in environment
process.env.OPENROUTER_API_KEY = API_KEY;

// Create mock CSV files
async function createMockFiles() {
    console.log("📁 Creating mock test files...");
    
    // Mock full results CSV
    const fullResultsContent = `True Cell Type,Marker List,Conversation History
T cell,"CD3D, CD3E, CD3G","Mock conversation about T cell analysis with CD3 markers"
Monocyte,"CD14, CD68, CD163","Mock conversation about monocyte markers"`;
    
    const fullResultsPath = path.join(TEST_RESULTS_DIR, "mock_full_results.csv");
    fs.writeFileSync(fullResultsPath, fullResultsContent);
    
    // Mock marker CSV
    const markerContent = `gene,avg_log2FC,p_val_adj,pct_1,pct_2
CD3D,2.5,1.2e-15,0.85,0.12
CD3E,2.1,3.4e-12,0.82,0.15
CD3G,1.8,5.6e-10,0.78,0.18
CD4,1.2,2.3e-05,0.65,0.35
CD8A,0.8,1.1e-03,0.45,0.55
FOXP3,0.3,0.12,0.15,0.85
CCR7,1.8,4.5e-08,0.75,0.25`;
    
    const markerPath = path.join(TEST_RESULTS_DIR, "mock_markers.csv");
    fs.writeFileSync(markerPath, markerContent);
    
    console.log("✅ Mock files created");
    return { fullResultsPath, markerPath };
}

async function testConversationSummarization() {
    console.log("\n📝 Testing Conversation Summarization...");
    
    const mockHistory = `
    ## ASSISTANT
    Based on CD3D, CD3E, CD3G markers, this cluster represents T lymphocytes.
    
    **Analysis:**
    - CD3D, CD3E, CD3G are core components of the CD3 complex
    - Essential for T-cell receptor signaling
    - High confidence for T cell identification
    
    **Conclusion:** T cells with high confidence based on CD3 complex expression.
    `;
    
    try {
        const summary = await summarizeConversationHistory(
            mockHistory,
            "openrouter",
            "google/gemini-2.5-flash-preview",
            0.1
        );
        
        if (summary && summary.length < mockHistory.length && summary.toLowerCase().includes("t cell")) {
            console.log("✅ Conversation summarization working");
            console.log(`   Original: ${mockHistory.length} chars → Summary: ${summary.length} chars`);
            return true;
        } else {
            console.log("⚠️  Summarization may have issues");
            return false;
        }
    } catch (error) {
        console.log(`❌ Summarization failed: ${error.message}`);
        return false;
    }
}

async function testQuickAnnotationBoost() {
    console.log("\n🚀 Testing Quick Annotation Boost...");
    
    const mockFiles = await createMockFiles();
    const fullResultsPath = mockFiles.fullResultsPath;
    const markerPath = mockFiles.markerPath;
    
    try {
        const result = await runCASSIAAnnotationboost({
            fullResultPath,
            marker: markerPath,
            clusterName: "T cell",
            majorClusterInfo: "Human lung tissue scRNA-seq",
            outputName: path.join(TEST_RESULTS_DIR, "quick_annotation_test"),
            numIterations: 1, // Minimal for speed
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            temperature: 0,
            conversationHistoryMode: "none", // Skip for speed
            searchStrategy: "breadth",
            reportStyle: "per_iteration"
        });
        
        if (result.status === 'success') {
            console.log("✅ Annotation boost API test PASSED");
            console.log(`   Execution time: ${result.execution_time}s`);
            
            if (result.summary_report_path && fs.existsSync(result.summary_report_path)) {
                console.log(`   ✅ HTML report generated: ${result.summary_report_path}`);
            }
            
            if (result.raw_text_path && fs.existsSync(result.raw_text_path)) {
                console.log(`   ✅ Raw conversation saved: ${result.raw_text_path}`);
            }
            
            return true;
        } else {
            console.log(`❌ Annotation boost failed: ${result.error_message}`);
            return false;
        }
        
    } catch (error) {
        console.log(`❌ API test error: ${error.message}`);
        return false;
    }
}

async function runQuickAPITest() {
    console.log("🧪 ANNOTATION BOOST API TEST");
    console.log("=" .repeat(40));
    console.log("Running minimal tests with real API calls");
    
    // Ensure test directory exists
    if (!fs.existsSync(TEST_RESULTS_DIR)) {
        fs.mkdirSync(TEST_RESULTS_DIR, { recursive: true });
    }
    
    const tests = [
        { name: "Conversation Summarization", func: testConversationSummarization },
        { name: "Quick Annotation Boost", func: testQuickAnnotationBoost }
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
    
    console.log("\n" + "=" .repeat(40));
    console.log(`📊 API Test Results: ${passed} passed, ${failed} failed`);
    
    if (passed === 2) {
        console.log("🎉 ALL API TESTS PASSED!");
        console.log("   ✅ Real LLM integration working");
        console.log("   ✅ End-to-end workflow functional");
        console.log("   🚀 Annotation Boost is fully operational!");
    } else if (passed > 0) {
        console.log("⚠️  Some API tests passed - check network/API key");
    } else {
        console.log("❌ API tests failed - check configuration");
    }
    
    console.log(`\n📁 Test results in: ${TEST_RESULTS_DIR}`);
}

runQuickAPITest().catch(console.error);