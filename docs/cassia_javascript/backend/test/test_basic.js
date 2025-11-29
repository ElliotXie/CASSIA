import { callLLM, runCASSIA, runCASSIABatch, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Set API key for testing
const API_KEY = "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";
setApiKey(API_KEY, "openrouter");

// Basic test configuration
const testConfig = {
    model: "google/gemini-2.5-flash-preview",
    temperature: 0,
    tissue: "blood",
    species: "human",
    provider: "openrouter",
    maxWorkers: 2  // Reduced for testing
};

async function testCallLLM() {
    console.log("\n=== Testing callLLM function ===");
    
    try {
        const prompt = "What is the primary function of T cells in the immune system? Answer in one sentence.";
        
        console.log("Testing OpenRouter provider...");
        const response = await callLLM(
            prompt,
            "openrouter",
            "google/gemini-2.5-flash-preview",
            null, // Let it get API key from environment
            0.1,
            100
        );
        
        console.log("Response:", response.substring(0, 200) + "...");
        console.log("✅ callLLM test passed");
        return true;
    } catch (error) {
        console.error("❌ callLLM test failed:", error.message);
        return false;
    }
}

async function testRunCASSIA() {
    console.log("\n=== Testing runCASSIA function ===");
    
    try {
        const testMarkers = ["CD3D", "CD3E", "CD3G", "CD8A", "CD4"];
        
        console.log(`Testing single cluster analysis with markers: ${testMarkers.join(', ')}`);
        
        const [result, conversationHistory] = await runCASSIA(
            testConfig.model,
            testConfig.temperature,
            testMarkers,
            testConfig.tissue,
            testConfig.species,
            "These are highly expressed markers from a cell cluster.",
            testConfig.provider,
            "v1"
        );
        
        console.log("Analysis Result:");
        console.log("- Main Cell Type:", result.main_cell_type);
        console.log("- Sub Cell Types:", result.sub_cell_types);
        console.log("- Iterations:", result.iterations);
        console.log("- Validation Passed:", result.validation_passed);
        
        console.log("✅ runCASSIA test passed");
        return true;
    } catch (error) {
        console.error("❌ runCASSIA test failed:", error.message);
        return false;
    }
}

async function testRunCASSIABatch() {
    console.log("\n=== Testing runCASSIA_batch function ===");
    
    try {
        const csvPath = path.join(__dirname, 'test_sample_data.csv');
        
        console.log(`Testing batch analysis with sample data from: ${csvPath}`);
        
        const results = await runCASSIABatch({
            marker: csvPath,
            outputName: path.join(__dirname, 'test_results', 'basic_test_results'),
            nGenes: 5,
            model: testConfig.model,
            temperature: testConfig.temperature,
            tissue: testConfig.tissue,
            species: testConfig.species,
            provider: testConfig.provider,
            maxWorkers: testConfig.maxWorkers,
            maxRetries: 0  // No retries for testing
        });
        
        console.log("Batch Analysis Results:");
        console.log("- Total Clusters:", results.total_clusters);
        console.log("- Successful Analyses:", results.successful_analyses);
        console.log("- Failed Analyses:", results.failed_analyses);
        console.log("- Output Files:", results.output_files);
        
        if (results.successful_analyses > 0) {
            console.log("\nSample Results:");
            for (const [cellType, details] of Object.entries(results.results)) {
                console.log(`- ${cellType}: ${details.analysis_result.main_cell_type}`);
                break; // Just show first result
            }
        }
        
        console.log("✅ runCASSIA_batch test passed");
        return true;
    } catch (error) {
        console.error("❌ runCASSIA_batch test failed:", error.message);
        return false;
    }
}

async function runAllTests() {
    console.log("🧪 Starting CASSIA JavaScript Tests");
    console.log("====================================");
    
    console.log("✅ API key set for testing");
    
    const tests = [
        { name: "callLLM", func: testCallLLM },
        { name: "runCASSIA", func: testRunCASSIA },
        { name: "runCASSIA_batch", func: testRunCASSIABatch }
    ];
    
    let passed = 0;
    let total = tests.length;
    
    for (const test of tests) {
        console.log(`\n📋 Running ${test.name} test...`);
        try {
            const result = await test.func();
            if (result) passed++;
        } catch (error) {
            console.error(`❌ ${test.name} test failed with error:`, error.message);
        }
    }
    
    console.log("\n" + "=".repeat(50));
    console.log(`📊 Test Results: ${passed}/${total} tests passed`);
    
    if (passed === total) {
        console.log("🎉 All tests passed! CASSIA JavaScript implementation is working correctly.");
    } else {
        console.log("⚠️  Some tests failed. Please check the error messages above.");
    }
}

// Run tests automatically
runAllTests().catch(console.error);

export { runAllTests };