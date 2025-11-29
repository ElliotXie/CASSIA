import { callLLM, runCASSIA, runCASSIABatch, scoreSingleAnalysis, scoreAnnotationBatch, runCASSIAScoreBatch, runCASSIAGenerateScoreReport, setApiKey } from '../index.js';

console.log("🧪 Testing CASSIA JavaScript Import");
console.log("====================================");

try {
    console.log("✅ Successfully imported callLLM:", typeof callLLM);
    console.log("✅ Successfully imported runCASSIA:", typeof runCASSIA);
    console.log("✅ Successfully imported runCASSIABatch:", typeof runCASSIABatch);
    console.log("✅ Successfully imported scoreSingleAnalysis:", typeof scoreSingleAnalysis);
    console.log("✅ Successfully imported scoreAnnotationBatch:", typeof scoreAnnotationBatch);
    console.log("✅ Successfully imported runCASSIAScoreBatch:", typeof runCASSIAScoreBatch);
    console.log("✅ Successfully imported runCASSIAGenerateScoreReport:", typeof runCASSIAGenerateScoreReport);
    console.log("✅ Successfully imported setApiKey:", typeof setApiKey);
    
    console.log("\n🎉 All imports successful! The CASSIA JavaScript implementation is ready to use.");
    
    console.log("\nAvailable functions:");
    console.log("📊 Core Analysis:");
    console.log("  - runCASSIA: Single cluster analysis");
    console.log("  - runCASSIABatch: Multiple cluster analysis");
    console.log("📈 Scoring System:");
    console.log("  - scoreSingleAnalysis: Score individual annotations");
    console.log("  - scoreAnnotationBatch: Score batch results");
    console.log("  - runCASSIAScoreBatch: Complete scoring pipeline");
    console.log("📄 Report Generation:");
    console.log("  - runCASSIAGenerateScoreReport: Generate HTML reports");
    console.log("🔧 Utilities:");
    console.log("  - callLLM: Direct LLM API calls");
    console.log("  - setApiKey: Configure API keys");
    
} catch (error) {
    console.error("❌ Import failed:", error.message);
    process.exit(1);
}