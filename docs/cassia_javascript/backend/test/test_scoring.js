import { scoreSingleAnalysis, scoreAnnotationBatch, runCASSIAScoreBatch, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Set API key for testing
const API_KEY = "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";
setApiKey(API_KEY, "openrouter");

async function testScoring() {
    console.log("🧪 CASSIA Scoring System Test");
    console.log("=============================");
    
    try {
        // Test 1: Single annotation scoring
        console.log("\n📋 Test 1: Single annotation scoring...");
        
        const testMarkers = "CD3D, CD3E, CD3G, CD8A, CD4";
        const testHistory = `Final Annotation Agent: Based on the provided markers CD3D, CD3E, CD3G, CD8A, and CD4, I can identify this as a T cell population. The presence of CD3 complex markers (CD3D, CD3E, CD3G) clearly indicates T cells, while CD8A and CD4 suggest both cytotoxic and helper T cell subpopulations. FINAL ANNOTATION COMPLETED: This is a mixed T lymphocyte population containing both CD8+ cytotoxic T cells and CD4+ helper T cells.`;
        
        const singleResult = await scoreSingleAnalysis(
            "human blood",
            testMarkers,
            testHistory,
            "google/gemini-2.5-flash-preview",
            "openrouter"
        );
        
        console.log("✅ Single scoring result:");
        console.log(`   - Score: ${singleResult.score}/100`);
        console.log(`   - Reasoning: ${singleResult.reasoning ? singleResult.reasoning.substring(0, 100) + "..." : "None"}`);
        
        // Test 2: Check if we have any existing batch results to score
        console.log("\n📋 Test 2: Looking for existing batch results to score...");
        
        const testResultsDir = path.join(__dirname, 'test_results');
        const fs = await import('fs');
        
        // Look for any existing _full.csv files
        let testFile = null;
        try {
            const files = await fs.promises.readdir(testResultsDir);
            const fullFiles = files.filter(f => f.includes('_full.csv'));
            if (fullFiles.length > 0) {
                testFile = path.join(testResultsDir, fullFiles[0]);
                console.log(`Found existing batch results: ${fullFiles[0]}`);
            }
        } catch (error) {
            console.log("No test_results directory found");
        }
        
        if (testFile) {
            console.log("\n📋 Test 3: Batch scoring on existing results...");
            
            try {
                const scoringResults = await runCASSIAScoreBatch({
                    inputFile: testFile,
                    outputFile: path.join(testResultsDir, 'scored_results.csv'),
                    maxWorkers: 2,
                    model: "google/gemini-2.5-flash-preview",
                    provider: "openrouter",
                    maxRetries: 1
                });
                
                console.log("✅ Batch scoring completed:");
                console.log(`   - Total rows processed: ${scoringResults.length}`);
                
                // Show sample scores
                const scoredRows = scoringResults.filter(row => row.Score !== null && row.Score !== '');
                if (scoredRows.length > 0) {
                    console.log(`   - Scored rows: ${scoredRows.length}`);
                    console.log(`   - Sample scores: ${scoredRows.slice(0, 3).map(row => `${row['True Cell Type']}: ${row.Score}`).join(', ')}`);
                }
                
            } catch (error) {
                console.error("❌ Batch scoring failed:", error.message);
            }
        } else {
            console.log("No existing batch results found. Run a batch analysis first to test scoring.");
            console.log("You can run: npm run test-basic");
        }
        
        console.log("\n🎉 Scoring system test completed!");
        
    } catch (error) {
        console.error("❌ Scoring test failed:", error.message);
    }
}

// Run the scoring test
testScoring();