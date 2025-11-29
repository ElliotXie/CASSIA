/**
 * Test suite for CASSIA Annotation Boost functionality
 * Tests all functions with mock data and real API calls
 */

import { 
    summarizeConversationHistory,
    promptHypothesisGenerator2,
    promptHypothesisGeneratorDepthFirst,
    promptHypothesisGenerator,
    promptHypothesisGeneratorAdditionalTaskDepthFirst,
    promptHypothesisGeneratorAdditionalTask,
    getMarkerInfo,
    extractGenesFromConversation,
    iterativeMarkerAnalysis,
    prepareAnalysisData,
    saveRawConversationText,
    generateSummaryReport,
    formatSummaryToHtml,
    runCASSIAAnnotationboost,
    runCASSIAAnnotationboostAdditionalTask
} from '../src/annotationBoost.js';

import fs from 'fs';
import path from 'path';

// Test configuration
const API_KEY = "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";
const TEST_RESULTS_DIR = "./test_results";

// Ensure test results directory exists
if (!fs.existsSync(TEST_RESULTS_DIR)) {
    fs.mkdirSync(TEST_RESULTS_DIR, { recursive: true });
}

/**
 * Create mock marker data for testing
 */
function createMockMarkerData() {
    return [
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
            gene: "CD3G",
            avg_log2FC: 1.8,
            p_val_adj: 5.6e-10,
            pct_1: 0.78,
            pct_2: 0.18
        },
        {
            gene: "CD14",
            avg_log2FC: 3.2,
            p_val_adj: 2.1e-20,
            pct_1: 0.95,
            pct_2: 0.08
        },
        {
            gene: "CD68",
            avg_log2FC: 2.8,
            p_val_adj: 4.5e-18,
            pct_1: 0.88,
            pct_2: 0.11
        },
        {
            gene: "CD79A",
            avg_log2FC: 3.0,
            p_val_adj: 1.8e-16,
            pct_1: 0.92,
            pct_2: 0.09
        }
    ];
}

/**
 * Create mock conversation history for testing
 */
function createMockConversationHistory() {
    return `## ASSISTANT
    Based on the provided marker genes CD3D, CD3E, and CD3G, I can analyze this cluster from human lung tissue.

    **Initial Assessment:**
    All three markers (CD3D, CD3E, CD3G) are core components of the CD3 complex, which is essential for T-cell receptor signaling. This is a strong indication of T lymphocytes.

    **Hypothesis 1: T cells (T lymphocytes)**
    - Strong evidence from CD3 complex markers
    - High confidence for general T cell identification

    To validate this hypothesis and determine T cell subtypes, I need to check additional markers:

    <check_genes>CD4,CD8A,CD8B,FOXP3,CCR7</check_genes>

    <reasoning>
    These genes will help distinguish between major T cell subtypes:
    - CD4: T helper cells
    - CD8A/CD8B: Cytotoxic T lymphocytes  
    - FOXP3: Regulatory T cells
    - CCR7: Naive vs effector/memory status
    </reasoning>

    ## USER
    gene	avg_log2FC	p_val_adj	pct_1	pct_2
    CD4	1.2	2.3e-05	0.65	0.35
    CD8A	0.8	1.1e-03	0.45	0.55
    CD8B	0.7	2.4e-03	0.42	0.58
    FOXP3	0.3	0.12	0.15	0.85
    CCR7	1.8	4.5e-08	0.75	0.25

    ## ASSISTANT
    **Analysis of Additional Markers:**

    The results show a mixed T cell population with both CD4+ and CD8+ cells present:
    - CD4 is moderately expressed (avg_log2FC: 1.2, 65% of cells)
    - CD8A/CD8B show lower expression (45% and 42% respectively)
    - FOXP3 is low, indicating few regulatory T cells
    - CCR7 is highly expressed, suggesting naive or central memory T cells

    **Final Assessment:**
    This cluster represents a mixed population of T lymphocytes, predominantly CD4+ T helper cells with some CD8+ cytotoxic T lymphocytes. The high CCR7 expression suggests these are primarily naive or central memory T cells rather than effector cells.

    FINAL ANNOTATION COMPLETED

    **Cell Type: T lymphocytes (mixed CD4+/CD8+ population)**
    **Confidence: High**
    **Key Supporting Markers:** CD3D, CD3E, CD3G (T cell identity), CD4 (helper T cells), CCR7 (naive/central memory)
    **Alternative Possibilities:** None - the CD3 complex expression definitively identifies these as T cells.`;
}

/**
 * Test prompt generation functions
 */
async function testPromptGenerators() {
    console.log("\n🔬 Testing Prompt Generators...");
    
    const majorClusterInfo = "Human lung tissue scRNA-seq";
    const genes = "CD3D, CD3E, CD3G";
    const history = "Previous analysis suggested T cells based on CD3 expression";
    const additionalTask = "determine if these are regulatory T cells";
    
    // Test basic prompt generator
    const prompt1 = promptHypothesisGenerator2(majorClusterInfo, genes, history);
    console.log("✅ Basic prompt generator working");
    
    // Test depth-first prompt generator
    const prompt2 = promptHypothesisGeneratorDepthFirst(majorClusterInfo, genes, history);
    console.log("✅ Depth-first prompt generator working");
    
    // Test breadth-first prompt generator
    const prompt3 = promptHypothesisGenerator(majorClusterInfo, genes, history);
    console.log("✅ Breadth-first prompt generator working");
    
    // Test additional task prompt generators
    const prompt4 = promptHypothesisGeneratorAdditionalTaskDepthFirst(majorClusterInfo, genes, history, additionalTask);
    const prompt5 = promptHypothesisGeneratorAdditionalTask(majorClusterInfo, genes, history, additionalTask);
    console.log("✅ Additional task prompt generators working");
    
    // Verify prompts contain expected elements
    const prompts = [prompt1, prompt2, prompt3, prompt4, prompt5];
    for (const prompt of prompts) {
        if (!prompt.includes(majorClusterInfo) || !prompt.includes(genes)) {
            throw new Error("Prompt missing required information");
        }
    }
    
    console.log("✅ All prompt generators test passed");
    return true;
}

/**
 * Test marker information extraction
 */
async function testMarkerInfo() {
    console.log("\n🧬 Testing Marker Info Extraction...");
    
    const mockData = createMockMarkerData();
    const geneList = ["CD3D", "CD3E", "UNKNOWN_GENE"];
    
    // Test marker info extraction
    const markerInfo = await getMarkerInfo(geneList, mockData);
    
    // Verify results
    if (!markerInfo.includes("CD3D") || !markerInfo.includes("CD3E")) {
        throw new Error("Marker info missing expected genes");
    }
    
    if (!markerInfo.includes("UNKNOWN_GENE")) {
        throw new Error("Marker info should mention unknown genes");
    }
    
    console.log("✅ Marker info extraction test passed");
    return true;
}

/**
 * Test gene extraction from conversation
 */
async function testGeneExtraction() {
    console.log("\n🔍 Testing Gene Extraction...");
    
    const conversation = `
    Based on the analysis, I need to check these markers:
    
    <check_genes>CD4,CD8A,CD8B</check_genes>
    
    And also these additional genes:
    
    <check_genes>FOXP3,CCR7,GZMB</check_genes>
    `;
    
    const genes = extractGenesFromConversation(conversation);
    
    // Verify extraction
    const expectedGenes = ["CD4", "CD8A", "CD8B", "FOXP3", "CCR7", "GZMB"];
    for (const gene of expectedGenes) {
        if (!genes.includes(gene)) {
            throw new Error(`Gene ${gene} not extracted`);
        }
    }
    
    console.log("✅ Gene extraction test passed");
    return true;
}

/**
 * Test conversation summarization
 */
async function testConversationSummarization() {
    console.log("\n📝 Testing Conversation Summarization...");
    
    const fullHistory = createMockConversationHistory();
    
    try {
        // Test summarization with API
        const summary = await summarizeConversationHistory(
            fullHistory,
            "openrouter",
            "google/gemini-2.5-flash-preview",
            0.1
        );
        
        // Verify summary is shorter and contains key information
        if (summary.length >= fullHistory.length) {
            throw new Error("Summary should be shorter than original");
        }
        
        if (!summary.toLowerCase().includes("t cell") || !summary.toLowerCase().includes("cd3")) {
            throw new Error("Summary missing key biological information");
        }
        
        console.log("✅ Conversation summarization test passed");
        return true;
        
    } catch (error) {
        console.log(`⚠️  Summarization test failed: ${error.message}`);
        console.log("   This may be due to API issues - continuing with other tests");
        return false;
    }
}

/**
 * Test iterative marker analysis (mock version)
 */
async function testIterativeAnalysisQuick() {
    console.log("\n🔄 Testing Iterative Analysis (Quick Mock)...");
    
    const mockData = createMockMarkerData();
    
    try {
        // Test with minimal iterations for speed
        const [result, messages] = await iterativeMarkerAnalysis({
            majorClusterInfo: "Human lung tissue",
            marker: mockData,
            commaSeparatedGenes: "CD3D, CD3E, CD3G",
            annotationHistory: "No previous analysis",
            numIterations: 1, // Minimal for testing
            provider: "openrouter",
            model: "google/gemini-2.5-flash-preview",
            temperature: 0,
            searchStrategy: "breadth"
        });
        
        // Verify basic structure
        if (!result || !Array.isArray(messages)) {
            throw new Error("Invalid analysis result structure");
        }
        
        if (messages.length === 0) {
            throw new Error("No messages generated");
        }
        
        console.log("✅ Iterative analysis test passed");
        return true;
        
    } catch (error) {
        console.log(`⚠️  Iterative analysis test failed: ${error.message}`);
        console.log("   This may be due to API issues - continuing with other tests");
        return false;
    }
}

/**
 * Test report generation functions
 */
async function testReportGeneration() {
    console.log("\n📊 Testing Report Generation...");
    
    // Create mock conversation messages
    const mockMessages = [
        {
            role: "user",
            content: "Analyze these T cell markers: CD3D, CD3E, CD3G"
        },
        {
            role: "assistant", 
            content: createMockConversationHistory()
        }
    ];
    
    // Test raw conversation text saving
    const rawTextPath = path.join(TEST_RESULTS_DIR, "test_raw_conversation.txt");
    await saveRawConversationText(mockMessages, rawTextPath);
    
    if (!fs.existsSync(rawTextPath)) {
        throw new Error("Raw conversation text file not created");
    }
    
    // Create mock summary text
    const mockSummary = `
    <OVERVIEW>
    This analysis examined T cell markers in human lung tissue, focusing on the CD3 complex components.
    </OVERVIEW>
    
    <INITIAL_ASSESSMENT>
    Initial markers CD3D, CD3E, CD3G strongly suggested T lymphocytes based on CD3 complex expression.
    </INITIAL_ASSESSMENT>
    
    <ITERATION_1>
    <HYPOTHESES>
    1. CD4+ T helper cells based on CD3 expression
    2. CD8+ cytotoxic T lymphocytes as alternative
    3. Mixed T cell population possibility
    </HYPOTHESES>
    
    <GENES_CHECKED>
    CD4, CD8A, CD8B, FOXP3, CCR7
    </GENES_CHECKED>
    
    <KEY_FINDINGS>
    Results confirmed mixed T cell population with predominant CD4+ cells and some CD8+ cells. High CCR7 expression indicates naive/central memory phenotype.
    </KEY_FINDINGS>
    </ITERATION_1>
    
    <FINAL_ANNOTATION>
    T lymphocytes (mixed CD4+/CD8+ population) with high confidence. Key markers: CD3D, CD3E, CD3G, CD4, CCR7.
    </FINAL_ANNOTATION>
    
    <MARKER_SUMMARY>
    Core T cell markers: CD3D, CD3E, CD3G
    Subtype markers: CD4 (helper), CD8A/CD8B (cytotoxic)
    Functional state: CCR7 (naive/central memory)
    </MARKER_SUMMARY>
    
    <RECOMMENDATIONS>
    Further analysis could examine specific T helper subtypes (Th1, Th2, Th17) or activation markers.
    </RECOMMENDATIONS>
    `;
    
    // Test HTML formatting
    const htmlPath = path.join(TEST_RESULTS_DIR, "test_summary_report.html");
    await formatSummaryToHtml(mockSummary, htmlPath, "breadth", "per_iteration");
    
    if (!fs.existsSync(htmlPath)) {
        throw new Error("HTML report not created");
    }
    
    // Verify HTML content
    const htmlContent = fs.readFileSync(htmlPath, 'utf-8');
    if (!htmlContent.includes("T lymphocytes") || !htmlContent.includes("CD3D")) {
        throw new Error("HTML report missing expected content");
    }
    
    console.log("✅ Report generation test passed");
    return true;
}

/**
 * Create mock CSV files for testing main functions
 */
async function createMockCSVFiles() {
    console.log("\n📁 Creating Mock CSV Files...");
    
    // Create mock full results CSV
    const fullResultsPath = path.join(TEST_RESULTS_DIR, "mock_full_results.csv");
    const fullResultsContent = `True Cell Type,Marker List,Conversation History
T cell,"CD3D, CD3E, CD3G","${createMockConversationHistory().replace(/"/g, '""')}"
Monocyte,"CD14, CD68, CD163","Previous analysis of monocyte markers"
B cell,"CD79A, CD79B, CD19","Previous analysis of B cell markers"`;
    
    fs.writeFileSync(fullResultsPath, fullResultsContent);
    
    // Create mock marker CSV
    const markerPath = path.join(TEST_RESULTS_DIR, "mock_markers.csv");
    const markerContent = `gene,avg_log2FC,p_val_adj,pct_1,pct_2
CD3D,2.5,1.2e-15,0.85,0.12
CD3E,2.1,3.4e-12,0.82,0.15
CD3G,1.8,5.6e-10,0.78,0.18
CD4,1.2,2.3e-05,0.65,0.35
CD8A,0.8,1.1e-03,0.45,0.55
CD8B,0.7,2.4e-03,0.42,0.58
FOXP3,0.3,0.12,0.15,0.85
CCR7,1.8,4.5e-08,0.75,0.25
CD14,3.2,2.1e-20,0.95,0.08
CD68,2.8,4.5e-18,0.88,0.11
CD163,2.0,3.2e-14,0.76,0.19
CD79A,3.0,1.8e-16,0.92,0.09
CD79B,2.7,2.5e-15,0.89,0.12
CD19,2.4,5.1e-13,0.84,0.16`;
    
    fs.writeFileSync(markerPath, markerContent);
    
    console.log(`✅ Mock CSV files created:
    - ${fullResultsPath}
    - ${markerPath}`);
    
    return { fullResultsPath, markerPath };
}

/**
 * Test main annotation boost function (quick version)
 */
async function testMainAnnotationBoostQuick() {
    console.log("\n🚀 Testing Main Annotation Boost Function (Quick)...");
    
    const { fullResultsPath, markerPath } = await createMockCSVFiles();
    
    try {
        // Test with minimal parameters for speed
        const result = await runCASSIAAnnotationboost({
            fullResultPath: fullResultsPath,
            marker: markerPath,
            clusterName: "T cell",
            majorClusterInfo: "Human lung tissue scRNA-seq",
            outputName: path.join(TEST_RESULTS_DIR, "test_annotation_boost"),
            numIterations: 1, // Minimal for testing
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            temperature: 0,
            conversationHistoryMode: "none", // Skip summarization for speed
            searchStrategy: "breadth",
            reportStyle: "per_iteration"
        });
        
        // Verify result structure
        if (!result || result.status !== 'success') {
            throw new Error(`Analysis failed: ${result?.error_message || 'Unknown error'}`);
        }
        
        // Check if files were created
        if (result.raw_text_path && !fs.existsSync(result.raw_text_path)) {
            throw new Error("Raw text file not created");
        }
        
        if (result.summary_report_path && !fs.existsSync(result.summary_report_path)) {
            throw new Error("Summary report file not created");
        }
        
        console.log("✅ Main annotation boost test passed");
        console.log(`   - Execution time: ${result.execution_time}s`);
        
        if (result.raw_text_path) {
            console.log(`   - Raw conversation: ${result.raw_text_path}`);
        }
        if (result.summary_report_path) {
            console.log(`   - Summary report: ${result.summary_report_path}`);
        }
        
        return true;
        
    } catch (error) {
        console.log(`⚠️  Main annotation boost test failed: ${error.message}`);
        console.log("   This may be due to API issues - check network connection");
        return false;
    }
}

/**
 * Test additional task function
 */
async function testAdditionalTaskFunction() {
    console.log("\n🎯 Testing Additional Task Function...");
    
    const { fullResultsPath, markerPath } = await createMockCSVFiles();
    
    try {
        const result = await runCASSIAAnnotationboostAdditionalTask({
            fullResultPath: fullResultsPath,
            marker: markerPath,
            clusterName: "T cell",
            majorClusterInfo: "Human lung tissue scRNA-seq",
            outputName: path.join(TEST_RESULTS_DIR, "test_additional_task"),
            numIterations: 1,
            model: "google/gemini-2.5-flash-preview",
            provider: "openrouter",
            additionalTask: "determine if these are regulatory T cells",
            temperature: 0,
            conversationHistoryMode: "none",
            searchStrategy: "depth",
            reportStyle: "total_summary"
        });
        
        // Verify result structure
        if (!result || result.status !== 'success') {
            throw new Error(`Analysis failed: ${result?.error_message || 'Unknown error'}`);
        }
        
        console.log("✅ Additional task function test passed");
        console.log(`   - Execution time: ${result.execution_time}s`);
        
        return true;
        
    } catch (error) {
        console.log(`⚠️  Additional task test failed: ${error.message}`);
        console.log("   This may be due to API issues - continuing");
        return false;
    }
}

/**
 * Run all tests
 */
async function runAllTests() {
    console.log("🧪 Starting CASSIA Annotation Boost Test Suite");
    console.log("=" .repeat(60));
    
    const tests = [
        { name: "Prompt Generators", func: testPromptGenerators },
        { name: "Marker Info Extraction", func: testMarkerInfo },
        { name: "Gene Extraction", func: testGeneExtraction },
        { name: "Conversation Summarization", func: testConversationSummarization },
        { name: "Iterative Analysis (Quick)", func: testIterativeAnalysisQuick },
        { name: "Report Generation", func: testReportGeneration },
        { name: "Main Annotation Boost (Quick)", func: testMainAnnotationBoostQuick },
        { name: "Additional Task Function", func: testAdditionalTaskFunction }
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
    
    console.log("\n" + "=" .repeat(60));
    console.log(`📊 Test Results: ${passed} passed, ${failed} failed`);
    
    if (failed === 0) {
        console.log("🎉 All tests passed! Annotation Boost implementation is working correctly.");
    } else if (passed > failed) {
        console.log("⚠️  Most tests passed. Some failures may be due to API connectivity.");
    } else {
        console.log("❌ Multiple test failures detected. Please check implementation.");
    }
    
    console.log(`\n📁 Test results saved to: ${TEST_RESULTS_DIR}`);
    console.log("   You can examine the generated files to verify functionality.");
    
    return { passed, failed };
}

// Export for module use
export { runAllTests };

// Run tests if called directly
if (import.meta.url === `file://${process.argv[1]}`) {
    runAllTests().catch(console.error);
}