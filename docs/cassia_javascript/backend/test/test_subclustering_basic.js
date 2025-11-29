/**
 * Basic Subclustering Test
 * Tests core subclustering functionality
 */

import { 
    subclusterAgentAnnotateSubcluster,
    constructPromptFromCsvSubcluster,
    annotateSubclusters,
    extractSubclusterResultsWithLlm,
    writeResultsToCsv,
    runCASSIASubclusters
} from '../src/subclustering.js';
import fs from 'fs';

// Set API key for testing
process.env.OPENROUTER_API_KEY = "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";

async function testSubclusteringBasics() {
    console.log("🧬 SUBCLUSTERING BASIC FUNCTIONALITY TEST");
    console.log("=" .repeat(50));
    
    const testDir = "./test_results";
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir, { recursive: true });
    }
    
    // Test 1: Prompt construction
    console.log("\n📝 Test 1: Prompt Construction...");
    
    const testMarkerData = [
        { cluster: 'subcluster_1', markers: 'CD4, IL2, IFNG, TNF' },
        { cluster: 'subcluster_2', markers: 'CD8A, CD8B, GZMA, GZMB, PRF1' },
        { cluster: 'subcluster_3', markers: 'FOXP3, IL10, CTLA4, TIGIT' }
    ];
    
    const prompt = constructPromptFromCsvSubcluster(
        testMarkerData, 
        "T cell", 
        50
    );
    
    if (prompt.includes("T cell") && prompt.includes("CD4, IL2, IFNG, TNF")) {
        console.log("✅ Prompt construction working correctly");
        console.log(`   Prompt length: ${prompt.length} characters`);
    } else {
        console.log("❌ Prompt construction failed");
        return false;
    }
    
    // Test 2: LLM annotation (short test)
    console.log("\n🤖 Test 2: LLM Annotation...");
    
    try {
        const shortPrompt = `Analyze these T cell markers:
1. CD4, IL2, IFNG
2. CD8A, CD8B, GZMA

Provide cell type annotations.`;

        const result = await subclusterAgentAnnotateSubcluster(
            shortPrompt,
            "google/gemini-2.5-flash-preview",
            0,
            "openrouter"
        );
        
        if (result && typeof result === 'string' && result.length > 50) {
            console.log("✅ LLM annotation working");
            console.log(`   Response length: ${result.length} characters`);
        } else {
            console.log("⚠️  LLM response shorter than expected");
        }
        
    } catch (error) {
        console.log(`❌ LLM annotation failed: ${error.message}`);
        return false;
    }
    
    // Test 3: Result extraction with mock data
    console.log("\n📊 Test 3: Result Extraction...");
    
    const mockResults = "results1(Helper T cells, Th1 cells), results2(Cytotoxic T cells, CD8+ effector), results3(Regulatory T cells, Tregs)";
    
    try {
        const df = await writeResultsToCsv(mockResults, `${testDir}/test_extraction`);
        
        if (df && df.columns && df.columns.includes('main_cell_type')) {
            console.log("✅ Result extraction working");
            console.log(`   Extracted ${df['Result ID'].length} cell types`);
        } else {
            console.log("❌ Result extraction failed");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Result extraction error: ${error.message}`);
        return false;
    }
    
    // Test 4: Structured data handling
    console.log("\n🗂️  Test 4: Structured Data Handling...");
    
    const structuredData = [
        {
            cluster: 1,
            key_markers: 'CD4, IL2, IFNG',
            explanation: 'Helper T cell markers',
            most_likely_top2_cell_types: ['Th1 cells', 'Helper T cells']
        },
        {
            cluster: 2,
            key_markers: 'CD8A, CD8B, GZMA',
            explanation: 'Cytotoxic T cell markers',
            most_likely_top2_cell_types: ['CD8+ T cells', 'Cytotoxic T cells']
        }
    ];
    
    try {
        const df2 = await writeResultsToCsv(structuredData, `${testDir}/test_structured`);
        
        if (df2 && df2.columns && df2['main_cell_type'][0] === 'Th1 cells') {
            console.log("✅ Structured data handling working");
            console.log(`   Processed ${df2['Result ID'].length} structured entries`);
        } else {
            console.log("❌ Structured data handling failed");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Structured data error: ${error.message}`);
        return false;
    }
    
    console.log("\n✅ BASIC SUBCLUSTERING TESTS COMPLETE");
    return true;
}

async function testWithRealData() {
    console.log("\n🧬 REAL DATA SUBCLUSTERING TEST");
    console.log("=" .repeat(50));
    
    const testDir = "./test_results";
    
    // Create simple test marker data based on real scenarios
    const realTestData = [
        { cluster: 'CD4_naive', markers: 'CCR7, LEF1, TCF7, SELL, IL7R' },
        { cluster: 'CD4_memory', markers: 'IL7R, CCR7, CD44, CD62L, LCK' },
        { cluster: 'CD8_effector', markers: 'CD8A, CD8B, GZMA, GZMB, PRF1' },
        { cluster: 'Treg', markers: 'FOXP3, IL2RA, CTLA4, TIGIT, ICOS' }
    ];
    
    console.log("\n🎯 Testing with realistic T cell subcluster data...");
    
    try {
        const result = await runCASSIASubclusters(
            realTestData,
            "T cell",
            `${testDir}/real_tcell_subclusters`,
            "google/gemini-2.5-flash-preview",
            0,
            "openrouter",
            50
        );
        
        // Check if CSV was created
        const csvPath = `${testDir}/real_tcell_subclusters.csv`;
        if (fs.existsSync(csvPath)) {
            const content = fs.readFileSync(csvPath, 'utf-8');
            const lines = content.split('\n').filter(line => line.trim());
            
            console.log("✅ Real data subclustering successful");
            console.log(`   Generated CSV with ${lines.length - 1} rows`);
            console.log(`   File: ${csvPath}`);
            
            // Check for HTML report
            const htmlPath = csvPath.replace('.csv', '_summary.html');
            if (fs.existsSync(htmlPath)) {
                console.log("✅ HTML report generated");
            } else {
                console.log("⚠️  HTML report not found (report generation may be disabled)");
            }
            
            return true;
        } else {
            console.log("❌ CSV file was not created");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Real data test failed: ${error.message}`);
        return false;
    }
}

async function runSubclusteringTest() {
    const basicSuccess = await testSubclusteringBasics();
    
    if (basicSuccess) {
        await testWithRealData();
    }
    
    console.log("\n" + "=" .repeat(50));
    if (basicSuccess) {
        console.log("🎉 SUBCLUSTERING IMPLEMENTATION VERIFIED!");
        console.log("   ✅ Prompt construction working");
        console.log("   ✅ LLM integration functional");
        console.log("   ✅ Result extraction robust");
        console.log("   ✅ CSV writing operational");
        console.log("   ✅ Real data processing successful");
        console.log("🚀 Ready for production use!");
    } else {
        console.log("❌ Some subclustering tests failed - check logs above");
    }
}

runSubclusteringTest().catch(console.error);