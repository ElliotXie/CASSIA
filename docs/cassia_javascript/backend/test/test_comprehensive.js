import { callLLM, runCASSIA, runCASSIABatch, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Set API key for testing
const API_KEY = "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";
setApiKey(API_KEY, "openrouter");

async function testAllParameters() {
    console.log("🧪 CASSIA Comprehensive Parameter Test");
    console.log("=====================================");
    
    console.log("\n📋 Test 1: Testing all validator involvement versions...");
    
    const testMarkers = ["CD3D", "CD3E", "CD8A"];
    
    // Test v0 validator
    try {
        console.log("Testing validator_involvement='v0'...");
        const [result_v0] = await runCASSIA(
            "google/gemini-2.5-flash-preview",
            0,
            testMarkers,
            "blood",
            "human",
            null,
            "openrouter",
            "v0"  // Should use coupling_validator_system_v0
        );
        console.log("✅ v0 validator works:", result_v0.main_cell_type);
    } catch (error) {
        console.error("❌ v0 validator failed:", error.message);
    }
    
    // Test v1 validator
    try {
        console.log("Testing validator_involvement='v1'...");
        const [result_v1] = await runCASSIA(
            "google/gemini-2.5-flash-preview",
            0,
            testMarkers,
            "blood",
            "human",
            null,
            "openrouter",
            "v1"  // Should use coupling_validator_system_v1
        );
        console.log("✅ v1 validator works:", result_v1.main_cell_type);
    } catch (error) {
        console.error("❌ v1 validator failed:", error.message);
    }
    
    // Test tissue-blind (should use v2 validator)
    try {
        console.log("Testing tissue-blind analysis...");
        const [result_blind] = await runCASSIA(
            "google/gemini-2.5-flash-preview",
            0,
            testMarkers,
            "tissue blind",  // Should trigger tissue-blind mode
            "human",
            null,
            "openrouter",
            "v1"
        );
        console.log("✅ Tissue-blind analysis works:", result_blind.main_cell_type);
    } catch (error) {
        console.error("❌ Tissue-blind analysis failed:", error.message);
    }
    
    console.log("\n📋 Test 2: Testing ranking methods...");
    
    const csvPath = path.join(__dirname, 'test_sample_data.csv');
    
    // Test different ranking methods
    const rankingMethods = ["avg_log2FC", "p_val_adj", "pct_diff", "Score"];
    
    for (const method of rankingMethods) {
        try {
            console.log(`Testing ranking method: ${method}...`);
            const results = await runCASSIABatch({
                marker: csvPath,
                outputName: path.join(__dirname, 'test_results', `ranking_${method}`),
                nGenes: 3,
                rankingMethod: method,
                ascending: method === "p_val_adj" ? true : false,  // p_val_adj should be ascending
                model: "google/gemini-2.5-flash-preview",
                maxWorkers: 1,  // Sequential for testing
                maxRetries: 0,
                provider: "openrouter"
            });
            console.log(`✅ ${method} ranking works: ${results.successful_analyses} analyses completed`);
        } catch (error) {
            console.error(`❌ ${method} ranking failed:`, error.message);
        }
    }
    
    console.log("\n📋 Test 3: Testing all batch parameters...");
    
    try {
        const fullParamResults = await runCASSIABatch({
            marker: csvPath,
            outputName: path.join(__dirname, 'test_results', 'full_params_test'),
            nGenes: 4,
            model: "google/gemini-2.5-flash-preview",
            temperature: 0.1,  // Non-zero temperature
            tissue: "immune system",
            species: "human",
            additionalInfo: "Test dataset with comprehensive parameters",
            celltypeColumn: "cluster",  // Explicitly specified
            geneColumnName: "gene",     // Explicitly specified
            maxWorkers: 2,
            provider: "openrouter",
            maxRetries: 1,
            rankingMethod: "avg_log2FC",
            ascending: false,
            validatorInvolvement: "v1",
            formatType: "seurat"  // Explicitly specified format
        });
        
        console.log("✅ All parameters test passed:");
        console.log(`   - Total clusters: ${fullParamResults.total_clusters}`);
        console.log(`   - Successful: ${fullParamResults.successful_analyses}`);
        console.log(`   - Output files: ${fullParamResults.output_files.join(', ')}`);
        
    } catch (error) {
        console.error("❌ Full parameters test failed:", error.message);
    }
    
    console.log("\n🎉 Comprehensive parameter testing complete!");
}

// Run the comprehensive test
testAllParameters().catch(console.error);