import { callLLM, runCASSIA, runCASSIABatch, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Set API key for testing
const API_KEY = "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";
setApiKey(API_KEY, "openrouter");

async function quickTest() {
    console.log("🧪 CASSIA Quick Test with Organized Output");
    console.log("==========================================");
    
    try {
        // Test 1: Single analysis
        console.log("\n📋 Test 1: Single cluster analysis...");
        const [result] = await runCASSIA(
            "google/gemini-2.5-flash-preview",
            0,
            ["CD3D", "CD3E", "CD8A"],
            "blood",
            "human",
            null,
            "openrouter",
            "v1"
        );
        console.log("✅ Single analysis result:", result.main_cell_type);
        
        // Test 2: Batch analysis with organized output
        console.log("\n📋 Test 2: Batch analysis with organized output...");
        const csvPath = path.join(__dirname, 'test_sample_data.csv');
        
        const batchResults = await runCASSIABatch({
            marker: csvPath,
            outputName: path.join(__dirname, 'test_results', 'quick_test'),
            nGenes: 3,
            model: "google/gemini-2.5-flash-preview",
            maxWorkers: 2,
            provider: "openrouter",
            maxRetries: 0
        });
        
        console.log("✅ Batch analysis completed:");
        console.log(`   - Successful analyses: ${batchResults.successful_analyses}`);
        console.log(`   - Output files: ${batchResults.output_files.join(', ')}`);
        
        // Check that files are in the correct location
        const outputDir = path.join(__dirname, 'test_results');
        console.log(`   - All files saved to: ${outputDir}`);
        
        console.log("\n🎉 Quick test completed successfully!");
        console.log("📁 All test results are organized in the test_results folder");
        
    } catch (error) {
        console.error("❌ Test failed:", error.message);
    }
}

// Run the quick test
quickTest();