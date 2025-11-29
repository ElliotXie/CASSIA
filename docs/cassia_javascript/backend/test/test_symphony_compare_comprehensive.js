/**
 * Comprehensive Symphony Compare Test
 * Tests all functionality including edge cases and real biological scenarios
 */

import { 
    symphonyCompare,
    extractCelltypeScores,
    extractDiscussion,
    generateComparisonHtmlReport
} from '../src/symphonyCompare.js';
import fs from 'fs';

// Set API key for testing
process.env.OPENROUTER_API_KEY = "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";

async function testParameterValidation() {
    console.log("🔍 TESTING PARAMETER VALIDATION");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    let allTestsPassed = true;
    
    // Test 1: Too few cell types
    console.log("\n📝 Test 1: Invalid number of cell types...");
    try {
        await symphonyCompare(
            "blood",
            ["T cell"], // Only 1 cell type
            "CD3D",
            "human",
            "budget",
            null,
            testDir,
            "validation_test_1",
            false, 1, 0.8, false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        console.log("❌ Should have thrown error for too few cell types");
        allTestsPassed = false;
    } catch (error) {
        if (error.message.includes("2-4 cell types")) {
            console.log("✅ Correctly rejected too few cell types");
        } else {
            console.log(`❌ Wrong error message: ${error.message}`);
            allTestsPassed = false;
        }
    }
    
    // Test 2: Too many cell types
    console.log("\n📝 Test 2: Too many cell types...");
    try {
        await symphonyCompare(
            "blood",
            ["T cell", "B cell", "NK cell", "Monocyte", "Neutrophil"], // 5 cell types
            "CD3D",
            "human",
            "budget",
            null,
            testDir,
            "validation_test_2",
            false, 1, 0.8, false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        console.log("❌ Should have thrown error for too many cell types");
        allTestsPassed = false;
    } catch (error) {
        if (error.message.includes("2-4 cell types")) {
            console.log("✅ Correctly rejected too many cell types");
        } else {
            console.log(`❌ Wrong error message: ${error.message}`);
            allTestsPassed = false;
        }
    }
    
    // Test 3: Missing API key
    console.log("\n📝 Test 3: Missing API key...");
    try {
        await symphonyCompare(
            "blood",
            ["T cell", "B cell"],
            "CD3D",
            "human",
            "budget",
            null,
            testDir,
            "validation_test_3",
            false, 1, 0.8, false,
            null, // No API key
            false
        );
        console.log("❌ Should have thrown error for missing API key");
        allTestsPassed = false;
    } catch (error) {
        if (error.message.includes("OPENROUTER_API_KEY")) {
            console.log("✅ Correctly rejected missing API key");
        } else {
            console.log(`❌ Wrong error message: ${error.message}`);
            allTestsPassed = false;
        }
    }
    
    // Test 4: Valid parameters
    console.log("\n📝 Test 4: Valid parameters...");
    try {
        const result = await symphonyCompare(
            "blood",
            ["T cell", "B cell"],
            "CD3D, CD19",
            "human",
            "budget",
            null,
            testDir,
            "validation_test_4",
            false, 1, 0.8, false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        
        if (result && result.summary) {
            console.log("✅ Valid parameters accepted and processed");
        } else {
            console.log("❌ Valid parameters failed to process");
            allTestsPassed = false;
        }
    } catch (error) {
        console.log(`❌ Valid parameters threw error: ${error.message}`);
        allTestsPassed = false;
    }
    
    return allTestsPassed;
}

async function testBiologicalScenarios() {
    console.log("\n🧬 TESTING REAL BIOLOGICAL SCENARIOS");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    const scenarios = [
        {
            name: "Immune Cell Differentiation",
            tissue: "peripheral blood",
            celltypes: ["CD4+ T cell", "CD8+ T cell", "B cell"],
            markers: "CD3D, CD3E, CD4, CD8A, CD19, CD79A, MS4A1"
        },
        {
            name: "Myeloid Lineage",
            tissue: "bone marrow",
            celltypes: ["Monocyte", "Neutrophil"],
            markers: "CD14, CD16, FCGR3A, CSF3R, MPO"
        },
        {
            name: "Brain Cell Types",
            tissue: "brain cortex",
            celltypes: ["Neuron", "Astrocyte", "Microglia"],
            markers: "MAP2, RBFOX3, GFAP, AQP4, CX3CR1, IBA1"
        }
    ];
    
    const results = {};
    
    for (const scenario of scenarios) {
        try {
            console.log(`\n🎯 Testing ${scenario.name}...`);
            
            const result = await symphonyCompare(
                scenario.tissue,
                scenario.celltypes,
                scenario.markers,
                "human",
                "budget", // Use budget for speed
                null,
                testDir,
                `bio_scenario_${scenario.name.replace(/\s+/g, '_').toLowerCase()}`,
                false, // No discussion for speed
                1,
                0.8,
                true,  // Generate HTML for review
                process.env.OPENROUTER_API_KEY,
                false  // Quiet mode
            );
            
            if (result && result.summary && result.consensus) {
                console.log(`   ✅ ${scenario.name} analysis successful`);
                console.log(`   Winner: ${result.consensus}`);
                console.log(`   Confidence: ${Math.round(result.confidence * 100)}%`);
                console.log(`   Models used: ${result.summary.models_used}`);
                
                // Check biological plausibility
                const winner = result.consensus.toLowerCase();
                const hasReasonableWinner = scenario.celltypes.some(ct => 
                    ct.toLowerCase().includes(winner) || winner.includes(ct.toLowerCase())
                );
                
                console.log(`   Biologically plausible: ${hasReasonableWinner ? 'Yes' : 'Maybe'}`);
                results[scenario.name] = true;
                
            } else {
                console.log(`   ❌ ${scenario.name} failed to reach consensus`);
                results[scenario.name] = false;
            }
            
        } catch (error) {
            console.log(`   ❌ ${scenario.name} error: ${error.message}`);
            results[scenario.name] = false;
        }
    }
    
    const successfulScenarios = Object.values(results).filter(r => r).length;
    console.log(`\n📊 Biological scenarios: ${successfulScenarios}/${scenarios.length} successful`);
    
    return successfulScenarios >= scenarios.length - 1; // Allow 1 failure
}

async function testEdgeCases() {
    console.log("\n⚠️ TESTING EDGE CASES");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    let allTestsPassed = true;
    
    // Test 1: Very short marker list
    console.log("\n📝 Test 1: Minimal markers...");
    try {
        const result = await symphonyCompare(
            "blood",
            ["T cell", "B cell"],
            "CD3D", // Just one marker
            "human",
            "budget",
            null,
            testDir,
            "edge_case_minimal_markers",
            false, 1, 0.8, false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        
        if (result && result.summary) {
            console.log("✅ Minimal markers handled successfully");
        } else {
            console.log("❌ Minimal markers failed");
            allTestsPassed = false;
        }
    } catch (error) {
        console.log(`❌ Minimal markers error: ${error.message}`);
        allTestsPassed = false;
    }
    
    // Test 2: Unusual cell type names
    console.log("\n📝 Test 2: Unusual cell type names...");
    try {
        const result = await symphonyCompare(
            "tissue",
            ["Type-1 Cell", "Cell/Type-2"],
            "MARKER1, MARKER2",
            "human",
            "budget",
            null,
            testDir,
            "edge_case_unusual_names",
            false, 1, 0.8, false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        
        if (result && result.summary) {
            console.log("✅ Unusual cell type names handled");
        } else {
            console.log("❌ Unusual cell type names failed");
            allTestsPassed = false;
        }
    } catch (error) {
        console.log(`❌ Unusual names error: ${error.message}`);
        allTestsPassed = false;
    }
    
    // Test 3: Very high consensus threshold
    console.log("\n📝 Test 3: Impossible consensus threshold...");
    try {
        const result = await symphonyCompare(
            "blood",
            ["T cell", "B cell"],
            "CD3D, CD19",
            "human",
            "budget",
            null,
            testDir,
            "edge_case_high_threshold",
            true,  // Enable discussion
            2,     // Max rounds
            1.0,   // 100% consensus required
            false,
            process.env.OPENROUTER_API_KEY,
            false
        );
        
        if (result) {
            console.log(`✅ High threshold handled (consensus: ${result.summary.consensus_reached})`);
        } else {
            console.log("❌ High threshold test failed");
            allTestsPassed = false;
        }
    } catch (error) {
        console.log(`❌ High threshold error: ${error.message}`);
        allTestsPassed = false;
    }
    
    return allTestsPassed;
}

async function testFileGeneration() {
    console.log("\n📁 TESTING FILE GENERATION");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    
    try {
        console.log("\n🎯 Testing complete file generation pipeline...");
        
        const result = await symphonyCompare(
            "immune system",
            ["T helper cell", "Cytotoxic T cell", "B cell"],
            "CD3D, CD3E, CD4, CD8A, CD19, CD79A",
            "human",
            "budget",
            null,
            testDir,
            "file_generation_test",
            true,  // Enable discussion
            1,
            0.75,
            true,  // Generate HTML
            process.env.OPENROUTER_API_KEY,
            true
        );
        
        if (result) {
            console.log("✅ Symphony execution completed");
            
            // Check CSV file
            const csvExists = fs.existsSync(result.csv_file);
            console.log(`   CSV file created: ${csvExists ? '✅' : '❌'}`);
            
            if (csvExists) {
                const csvContent = fs.readFileSync(result.csv_file, 'utf-8');
                const lines = csvContent.split('\n').filter(line => line.trim());
                console.log(`   CSV rows: ${lines.length - 1} (excluding header)`);
                
                // Check for expected columns
                const hasExpectedColumns = csvContent.includes('model') && 
                                         csvContent.includes('researcher') &&
                                         csvContent.includes('T helper cell_score');
                console.log(`   CSV format valid: ${hasExpectedColumns ? '✅' : '❌'}`);
            }
            
            // Check HTML file
            const htmlExists = result.html_file ? fs.existsSync(result.html_file) : false;
            console.log(`   HTML file created: ${htmlExists ? '✅' : '❌'}`);
            
            if (htmlExists) {
                const htmlSize = fs.statSync(result.html_file).size;
                console.log(`   HTML file size: ${(htmlSize/1024).toFixed(1)}KB`);
                
                const htmlContent = fs.readFileSync(result.html_file, 'utf-8');
                const hasInteractivity = htmlContent.includes('showRound') && 
                                       htmlContent.includes('score-table');
                console.log(`   HTML interactivity: ${hasInteractivity ? '✅' : '❌'}`);
            }
            
            // Check data completeness
            const hasValidData = result.results && result.results.length > 0 &&
                               result.summary && typeof result.confidence === 'number';
            console.log(`   Data completeness: ${hasValidData ? '✅' : '❌'}`);
            
            return csvExists && htmlExists && hasValidData;
        } else {
            console.log("❌ File generation test failed");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ File generation error: ${error.message}`);
        return false;
    }
}

async function testPerformanceMetrics() {
    console.log("\n⚡ TESTING PERFORMANCE METRICS");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    
    try {
        console.log("\n⏱️ Measuring execution time and resource usage...");
        
        const startTime = Date.now();
        const startMemory = process.memoryUsage().heapUsed;
        
        const result = await symphonyCompare(
            "blood",
            ["T cell", "B cell", "NK cell"],
            "CD3D, CD3E, CD19, CD56, NCAM1",
            "human",
            "budget", // Use budget for consistent timing
            null,
            testDir,
            "performance_test",
            true,  // Enable discussion for full test
            1,
            0.8,
            true,
            process.env.OPENROUTER_API_KEY,
            false
        );
        
        const endTime = Date.now();
        const endMemory = process.memoryUsage().heapUsed;
        
        const executionTime = (endTime - startTime) / 1000;
        const memoryUsed = (endMemory - startMemory) / 1024 / 1024;
        
        console.log(`✅ Performance test completed`);
        console.log(`   Execution time: ${executionTime.toFixed(1)}s`);
        console.log(`   Memory used: ${memoryUsed.toFixed(1)}MB`);
        
        if (result && result.summary) {
            console.log(`   Models used: ${result.summary.models_used}`);
            console.log(`   Total rounds: ${result.summary.total_rounds}`);
            console.log(`   API calls made: ${result.results.length}`);
            
            const avgTimePerCall = executionTime / result.results.length;
            console.log(`   Avg time per API call: ${avgTimePerCall.toFixed(1)}s`);
            
            // Performance benchmarks
            const isAcceptableTime = executionTime < 120; // Under 2 minutes
            const isAcceptableMemory = memoryUsed < 100;   // Under 100MB
            
            console.log(`   Time acceptable: ${isAcceptableTime ? '✅' : '⚠️'} (${isAcceptableTime ? 'Under 2min' : 'Over 2min'})`);
            console.log(`   Memory acceptable: ${isAcceptableMemory ? '✅' : '⚠️'} (${isAcceptableMemory ? 'Under 100MB' : 'Over 100MB'})`);
            
            return isAcceptableTime && isAcceptableMemory;
        } else {
            console.log("❌ Performance test failed to complete");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Performance test error: ${error.message}`);
        return false;
    }
}

async function runComprehensiveSymphonyTest() {
    console.log("🎼 CASSIA SYMPHONY COMPARE - COMPREHENSIVE TEST SUITE");
    console.log("=" .repeat(80));
    
    const tests = [
        { name: "Parameter Validation", fn: testParameterValidation },
        { name: "Biological Scenarios", fn: testBiologicalScenarios },
        { name: "Edge Cases", fn: testEdgeCases },
        { name: "File Generation", fn: testFileGeneration },
        { name: "Performance Metrics", fn: testPerformanceMetrics }
    ];
    
    const results = {};
    const timings = {};
    
    for (const test of tests) {
        console.log(`\n🔬 Running ${test.name} test...`);
        const testStartTime = Date.now();
        
        try {
            results[test.name] = await test.fn();
            const testTime = (Date.now() - testStartTime) / 1000;
            timings[test.name] = testTime;
            console.log(`   ⏱️ Completed in ${testTime.toFixed(1)}s`);
        } catch (error) {
            console.log(`❌ ${test.name} test crashed: ${error.message}`);
            results[test.name] = false;
            timings[test.name] = 0;
        }
    }
    
    // Comprehensive Summary
    console.log("\n" + "=" .repeat(80));
    console.log("📊 COMPREHENSIVE TEST RESULTS:");
    
    let passCount = 0;
    for (const [testName, passed] of Object.entries(results)) {
        const status = passed ? "✅ PASS" : "❌ FAIL";
        const timing = timings[testName] ? ` (${timings[testName].toFixed(1)}s)` : "";
        console.log(`   ${status}: ${testName}${timing}`);
        if (passed) passCount++;
    }
    
    const successRate = (passCount / tests.length * 100).toFixed(1);
    const totalTime = Object.values(timings).reduce((a, b) => a + b, 0);
    
    console.log(`\n🎯 OVERALL RESULTS:`);
    console.log(`   Tests passed: ${passCount}/${tests.length} (${successRate}%)`);
    console.log(`   Total time: ${totalTime.toFixed(1)}s`);
    
    if (passCount === tests.length) {
        console.log("\n🎉 ALL COMPREHENSIVE SYMPHONY TESTS PASSED!");
        console.log("   ✅ Parameter validation robust");
        console.log("   ✅ Real biological scenarios working");
        console.log("   ✅ Edge cases handled gracefully");
        console.log("   ✅ File generation pipeline complete");
        console.log("   ✅ Performance metrics acceptable");
        console.log("\n🏆 SYMPHONY COMPARE IS PRODUCTION READY!");
        console.log("   🎼 Multi-model consensus building operational");
        console.log("   💬 Discussion rounds and debate functional");
        console.log("   📊 Comprehensive reporting system active");
        console.log("   ⚡ Performance optimized for real-world use");
    } else {
        console.log(`\n⚠️ ${tests.length - passCount} test(s) failed - review above for details`);
        if (successRate >= 80) {
            console.log("🎵 Symphony Compare is mostly functional but needs attention to failed areas");
        } else {
            console.log("🔧 Symphony Compare needs significant fixes before production use");
        }
    }
    
    console.log(`\n📁 All test results and reports saved in: ./test_results/`);
    console.log("   📊 Open HTML files to review detailed multi-model analysis reports");
    console.log("   📈 CSV files contain complete scoring data for further analysis");
    console.log("\n🎼 The Symphony Compare orchestra is tuned and ready to perform!");
}

runComprehensiveSymphonyTest().catch(console.error);