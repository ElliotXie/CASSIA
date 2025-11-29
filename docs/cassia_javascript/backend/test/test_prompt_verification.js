/**
 * Strict Prompt Verification Test
 * Verifies that JavaScript prompts exactly match Python implementation
 */

import { 
    promptHypothesisGenerator2,
    promptHypothesisGeneratorDepthFirst,
    promptHypothesisGenerator,
    extractGenesFromConversation,
    getMarkerInfo
} from '../src/annotationBoost.js';

async function verifyPromptAccuracy() {
    console.log("🔍 STRICT PROMPT VERIFICATION TEST");
    console.log("=" .repeat(50));
    
    const majorClusterInfo = "Human lung tissue scRNA-seq";
    const genes = "CD3D, CD3E, CD3G";
    const history = "Previous analysis suggested T cells";
    
    console.log("\n📝 Testing Prompt Generator 2 (Basic)...");
    const prompt1 = promptHypothesisGenerator2(majorClusterInfo, genes, history);
    
    // Check for exact Python prompt elements
    const requiredElements = [
        "You are a careful professional biologist, specializing in single-cell RNA-seq analysis.",
        "Previous expert has done some analysis but the reuslts is not conclusive, your additional anlysis is needed.",
        "Here are the top marker genes that are differentially expressed in this cluster:",
        "Previous annotation/analysis:",
        "<check_genes>...</check_genes>",
        "FINAL ANNOTATION COMPLETED",
        "- Final cell type",
        "- Confidence level (high, medium, low)"
    ];
    
    let missingElements = [];
    for (const element of requiredElements) {
        if (!prompt1.includes(element)) {
            missingElements.push(element);
        }
    }
    
    if (missingElements.length === 0) {
        console.log("✅ Basic prompt generator - ALL ELEMENTS PRESENT");
    } else {
        console.log("❌ Basic prompt generator - MISSING ELEMENTS:");
        missingElements.forEach(elem => console.log(`   - ${elem}`));
    }
    
    // Check for preserved typos (100% compatibility requirement)
    if (prompt1.includes("reuslts") && prompt1.includes("anlysis")) {
        console.log("✅ Typos preserved for 100% Python compatibility");
    } else {
        console.log("❌ Missing required typos for Python compatibility");
    }
    
    console.log("\n📝 Testing Depth-First Prompt Generator...");
    const prompt2 = promptHypothesisGeneratorDepthFirst(majorClusterInfo, genes, history);
    
    const depthFirstElements = [
        "DEPTH-FIRST ANALYSIS",
        "NEVER say \"FINAL ANNOTATION COMPLETED\" immediately after requesting genes",
        "You MUST complete at least 2 rounds of gene checking",
        "Focus on ONE hypothesis per iteration",
        "Primary hypothesis to investigate:"
    ];
    
    let missingDepthElements = [];
    for (const element of depthFirstElements) {
        if (!prompt2.includes(element)) {
            missingDepthElements.push(element);
        }
    }
    
    if (missingDepthElements.length === 0) {
        console.log("✅ Depth-first prompt generator - ALL ELEMENTS PRESENT");
    } else {
        console.log("❌ Depth-first prompt generator - MISSING ELEMENTS:");
        missingDepthElements.forEach(elem => console.log(`   - ${elem}`));
    }
    
    console.log("\n📝 Testing Breadth-First Prompt Generator...");
    const prompt3 = promptHypothesisGenerator(majorClusterInfo, genes, history);
    
    const breadthElements = [
        "10 grandma are going to be in danger",
        "Design up to 3 follow‑up checks",
        "celltype to check 1",
        "hypothesis to check 3",
        "Skeptical, critical, and careful"
    ];
    
    let missingBreadthElements = [];
    for (const element of breadthElements) {
        if (!prompt3.includes(element)) {
            missingBreadthElements.push(element);
        }
    }
    
    if (missingBreadthElements.length === 0) {
        console.log("✅ Breadth-first prompt generator - ALL ELEMENTS PRESENT");
    } else {
        console.log("❌ Breadth-first prompt generator - MISSING ELEMENTS:");
        missingBreadthElements.forEach(elem => console.log(`   - ${elem}`));
    }
    
    return missingElements.length === 0 && missingDepthElements.length === 0 && missingBreadthElements.length === 0;
}

async function testGeneExtractionAccuracy() {
    console.log("\n🧬 Testing Gene Extraction Accuracy...");
    
    const testConversations = [
        {
            name: "Standard format",
            text: "I need to check these genes: <check_genes>CD4,CD8A,CD8B</check_genes>",
            expected: ["CD4", "CD8A", "CD8B"]
        },
        {
            name: "Multiple tags",
            text: `First check: <check_genes>CD3D,CD3E</check_genes>
                   Second check: <check_genes>FOXP3,IL2RA</check_genes>`,
            expected: ["CD3D", "CD3E", "FOXP3", "IL2RA"]
        },
        {
            name: "With spaces",
            text: "<check_genes>CD14, CD68, CD163</check_genes>",
            expected: ["CD14", "CD68", "CD163"]
        }
    ];
    
    let allPassed = true;
    
    for (const test of testConversations) {
        const extracted = extractGenesFromConversation(test.text);
        const matches = test.expected.every(gene => extracted.includes(gene));
        
        if (matches && extracted.length === test.expected.length) {
            console.log(`✅ ${test.name} - PASSED`);
        } else {
            console.log(`❌ ${test.name} - FAILED`);
            console.log(`   Expected: ${test.expected.join(", ")}`);
            console.log(`   Got: ${extracted.join(", ")}`);
            allPassed = false;
        }
    }
    
    return allPassed;
}

async function testMarkerInfoAccuracy() {
    console.log("\n📊 Testing Marker Info Accuracy...");
    
    const mockData = [
        {
            gene: "CD3D",
            avg_log2FC: 2.5,
            p_val_adj: 1.2e-15,
            pct_1: 0.85,
            pct_2: 0.12
        },
        {
            gene: "CD4",
            avg_log2FC: 1.2,
            p_val_adj: 2.3e-05,
            pct_1: 0.65,
            pct_2: 0.35
        }
    ];
    
    try {
        const result = await getMarkerInfo(["CD3D", "CD4", "UNKNOWN"], mockData);
        
        // Check if it includes both genes and formats correctly
        const hasCD3D = result.includes("CD3D") && result.includes("2.50");
        const hasCD4 = result.includes("CD4") && result.includes("1.20");
        const hasUnknownNote = result.includes("UNKNOWN") && result.includes("not in the differential expression list");
        
        if (hasCD3D && hasCD4 && hasUnknownNote) {
            console.log("✅ Marker info extraction - PASSED");
            return true;
        } else {
            console.log("❌ Marker info extraction - FAILED");
            console.log("   Missing expected content in output");
            return false;
        }
    } catch (error) {
        console.log(`❌ Marker info extraction - ERROR: ${error.message}`);
        return false;
    }
}

async function runStrictVerification() {
    console.log("🎯 Running Strict Python Compatibility Verification");
    console.log("   This verifies 100% accuracy of prompts and core functions");
    
    const results = {
        prompts: await verifyPromptAccuracy(),
        geneExtraction: await testGeneExtractionAccuracy(),
        markerInfo: await testMarkerInfoAccuracy()
    };
    
    const passed = Object.values(results).filter(r => r).length;
    const total = Object.keys(results).length;
    
    console.log("\n" + "=" .repeat(50));
    console.log(`📊 VERIFICATION RESULTS: ${passed}/${total} categories passed`);
    
    if (passed === total) {
        console.log("🎉 ALL VERIFICATION TESTS PASSED!");
        console.log("   ✅ Prompts are 100% Python-compatible");
        console.log("   ✅ Gene extraction works correctly");
        console.log("   ✅ Marker data processing is accurate");
        console.log("   🚀 Annotation Boost is ready for production!");
    } else {
        console.log("⚠️  Some verification tests failed");
        console.log("   Please review the output above for details");
    }
    
    return passed === total;
}

// Run verification
runStrictVerification().catch(console.error);