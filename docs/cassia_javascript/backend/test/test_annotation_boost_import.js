/**
 * Quick import test for annotation boost functionality
 * Verifies all functions can be imported without API calls
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

async function testImports() {
    console.log("🧪 Testing Annotation Boost Imports...");
    
    // Test prompt generators (no API calls)
    const prompt1 = promptHypothesisGenerator2("Human PBMC", "CD3D,CD3E", "Previous analysis");
    const prompt2 = promptHypothesisGeneratorDepthFirst("Human PBMC", "CD3D,CD3E", "Previous analysis");
    const prompt3 = promptHypothesisGenerator("Human PBMC", "CD3D,CD3E", "Previous analysis");
    const prompt4 = promptHypothesisGeneratorAdditionalTaskDepthFirst("Human PBMC", "CD3D,CD3E", "Previous analysis", "additional task");
    const prompt5 = promptHypothesisGeneratorAdditionalTask("Human PBMC", "CD3D,CD3E", "Previous analysis", "additional task");
    
    console.log("✅ Prompt generators imported successfully");
    
    // Test gene extraction (no API calls)
    const conversation = "<check_genes>CD4,CD8A</check_genes>";
    const genes = extractGenesFromConversation(conversation);
    
    if (genes.includes("CD4") && genes.includes("CD8A")) {
        console.log("✅ Gene extraction working");
    } else {
        throw new Error("Gene extraction failed");
    }
    
    // Test marker info (no API calls)
    const mockData = [
        { gene: "CD3D", avg_log2FC: 2.5, p_val_adj: 1e-15 },
        { gene: "CD3E", avg_log2FC: 2.1, p_val_adj: 1e-12 }
    ];
    
    const markerInfo = await getMarkerInfo(["CD3D", "UNKNOWN"], mockData);
    if (markerInfo.includes("CD3D")) {
        console.log("✅ Marker info extraction working");
    } else {
        throw new Error("Marker info extraction failed");
    }
    
    console.log("✅ All annotation boost functions imported successfully!");
    console.log("🎉 Annotation Boost module is ready to use!");
    
    return true;
}

testImports().catch(console.error);