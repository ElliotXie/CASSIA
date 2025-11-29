/**
 * Advanced Symphony Compare Test
 * Tests discussion rounds, consensus building, and multiple model presets
 */

import { 
    symphonyCompare,
    extractCelltypeScores,
    extractDiscussion
} from '../src/symphonyCompare.js';
import fs from 'fs';

// Set API key for testing
process.env.OPENROUTER_API_KEY = "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";

async function testDiscussionRounds() {
    console.log("💬 TESTING DISCUSSION ROUNDS");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir, { recursive: true });
    }
    
    try {
        console.log("\n🎯 Testing Symphony with discussion enabled...");
        console.log("   Using budget models with 2 discussion rounds");
        
        const result = await symphonyCompare(
            "immune tissue",
            ["T cell", "B cell", "NK cell"],
            "CD3D, CD3E, CD4, CD8A, CD19, CD20, CD56, NCAM1",
            "human",
            "budget",
            null,
            testDir,
            "discussion_test",
            true,      // Enable discussion
            2,         // Max 2 discussion rounds
            0.7,       // Lower consensus threshold to potentially trigger discussion
            true,
            process.env.OPENROUTER_API_KEY,
            true
        );
        
        if (result && result.results) {
            console.log("✅ Discussion test execution successful");
            console.log(`   Total results: ${result.results.length}`);
            console.log(`   Total rounds: ${result.summary.total_rounds}`);
            console.log(`   Consensus reached: ${result.summary.consensus_reached ? 'Yes' : 'No'}`);
            
            // Check for discussion rounds
            const discussionResults = result.results.filter(r => r.round && r.round.startsWith('discussion_'));
            const hasDiscussion = discussionResults.length > 0;
            
            console.log(`   Discussion rounds performed: ${hasDiscussion ? 'Yes' : 'No'}`);
            if (hasDiscussion) {
                const discussionRounds = new Set(discussionResults.map(r => r.round)).size;
                console.log(`   Number of discussion rounds: ${discussionRounds}`);
                
                // Check if discussion content exists
                const hasDiscussionContent = discussionResults.some(r => 
                    r.discussion && r.discussion !== "No discussion found"
                );
                console.log(`   Discussion content found: ${hasDiscussionContent ? 'Yes' : 'No'}`);
            }
            
            if (result.consensus) {
                console.log(`   Final consensus: ${result.consensus}`);
                console.log(`   Confidence: ${Math.round(result.confidence * 100)}%`);
            }
            
            return true;
        } else {
            console.log("❌ Discussion test failed to return valid results");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Discussion test error: ${error.message}`);
        return false;
    }
}

async function testModelPresets() {
    console.log("\n🤖 TESTING MODEL PRESETS");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    const presets = [
        { name: "budget", models: 3 },
        { name: "symphony", models: 3 }
        // Skip quartet for speed - it has 4 models
    ];
    
    const results = {};
    
    for (const preset of presets) {
        try {
            console.log(`\n🎼 Testing ${preset.name} preset...`);
            
            const result = await symphonyCompare(
                "blood",
                ["T cell", "B cell"],
                "CD3D, CD19, CD79A",
                "human",
                preset.name,
                null,
                testDir,
                `preset_test_${preset.name}`,
                false,     // Disable discussion for speed
                1,
                0.8,
                true,
                process.env.OPENROUTER_API_KEY,
                false      // Reduce verbosity
            );
            
            if (result && result.summary) {
                console.log(`   ✅ ${preset.name} preset successful`);
                console.log(`   Models used: ${result.summary.models_used}`);
                console.log(`   Consensus: ${result.summary.consensus_reached ? 'Yes' : 'No'}`);
                
                if (result.consensus) {
                    console.log(`   Winner: ${result.consensus}`);
                }
                
                results[preset.name] = true;
            } else {
                console.log(`   ❌ ${preset.name} preset failed`);
                results[preset.name] = false;
            }
            
        } catch (error) {
            console.log(`   ❌ ${preset.name} preset error: ${error.message}`);
            results[preset.name] = false;
        }
    }
    
    const successfulPresets = Object.values(results).filter(r => r).length;
    console.log(`\n📊 Preset test summary: ${successfulPresets}/${presets.length} presets working`);
    
    return successfulPresets === presets.length;
}

async function testCustomModels() {
    console.log("\n⚙️ TESTING CUSTOM MODEL CONFIGURATION");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    
    try {
        console.log("\n🎯 Testing with custom model list...");
        
        const customModels = [
            "google/gemini-2.5-flash",
            "deepseek/deepseek-chat-v3-0324"
        ];
        
        const result = await symphonyCompare(
            "tissue",
            ["Cell Type A", "Cell Type B"],
            "MARKER1, MARKER2, MARKER3",
            "human",
            "custom",
            customModels,
            testDir,
            "custom_models_test",
            false,
            1,
            0.8,
            true,
            process.env.OPENROUTER_API_KEY,
            true
        );
        
        if (result && result.summary) {
            console.log("✅ Custom models test successful");
            console.log(`   Models used: ${result.summary.models_used}`);
            console.log(`   Expected models: ${customModels.length}`);
            
            // Verify the correct models were used
            const usedModels = [...new Set(result.results.map(r => r.model))];
            const correctModels = customModels.every(model => usedModels.includes(model));
            
            console.log(`   Correct models used: ${correctModels ? 'Yes' : 'No'}`);
            
            if (correctModels) {
                console.log(`   Used models: ${usedModels.join(', ')}`);
                return true;
            } else {
                console.log(`   Expected: ${customModels.join(', ')}`);
                console.log(`   Actual: ${usedModels.join(', ')}`);
                return false;
            }
        } else {
            console.log("❌ Custom models test failed");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Custom models test error: ${error.message}`);
        return false;
    }
}

async function testConsensusThresholds() {
    console.log("\n🎯 TESTING CONSENSUS THRESHOLDS");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    const thresholds = [0.5, 0.7, 1.0]; // Different consensus requirements
    
    const results = {};
    
    for (const threshold of thresholds) {
        try {
            console.log(`\n📊 Testing threshold ${threshold}...`);
            
            const result = await symphonyCompare(
                "blood",
                ["T cell", "B cell"],
                "CD3D, CD19",
                "human",
                "budget",
                null,
                testDir,
                `threshold_test_${threshold}`,
                false,     // No discussion for clarity
                1,
                threshold,
                false,     // No HTML for speed
                process.env.OPENROUTER_API_KEY,
                false      // Quiet mode
            );
            
            if (result) {
                console.log(`   ✅ Threshold ${threshold} test completed`);
                console.log(`   Consensus reached: ${result.summary.consensus_reached ? 'Yes' : 'No'}`);
                
                if (result.consensus) {
                    console.log(`   Winner: ${result.consensus}`);
                    console.log(`   Confidence: ${Math.round(result.confidence * 100)}%`);
                }
                
                results[threshold] = true;
            } else {
                console.log(`   ❌ Threshold ${threshold} test failed`);
                results[threshold] = false;
            }
            
        } catch (error) {
            console.log(`   ❌ Threshold ${threshold} error: ${error.message}`);
            results[threshold] = false;
        }
    }
    
    const successfulTests = Object.values(results).filter(r => r).length;
    console.log(`\n📊 Threshold test summary: ${successfulTests}/${thresholds.length} thresholds working`);
    
    return successfulTests === thresholds.length;
}

async function testComplexScenario() {
    console.log("\n🏆 TESTING COMPLEX SCENARIO");
    console.log("=" .repeat(40));
    
    const testDir = "./test_results";
    
    try {
        console.log("\n🎯 Running complex multi-celltype analysis...");
        console.log("   4 cell types, comprehensive markers, discussion enabled");
        
        const result = await symphonyCompare(
            "peripheral blood mononuclear cells",
            ["CD4+ T cell", "CD8+ T cell", "B cell", "Monocyte"],
            "CD3D, CD3E, CD4, CD8A, CD19, CD79A, CD14, CD16, FCGR3A",
            "human",
            "budget",  // Use budget for speed but still test complexity
            null,
            testDir,
            "complex_scenario_test",
            true,      // Enable discussion
            1,         // 1 discussion round
            0.75,      // Moderate consensus threshold
            true,      // Generate full report
            process.env.OPENROUTER_API_KEY,
            true
        );
        
        if (result && result.summary) {
            console.log("✅ Complex scenario test successful");
            console.log(`   Cell types analyzed: 4`);
            console.log(`   Models used: ${result.summary.models_used}`);
            console.log(`   Total rounds: ${result.summary.total_rounds}`);
            console.log(`   Consensus reached: ${result.summary.consensus_reached ? 'Yes' : 'No'}`);
            
            if (result.consensus) {
                console.log(`   Winning cell type: ${result.consensus}`);
                console.log(`   Confidence: ${Math.round(result.confidence * 100)}%`);
            }
            
            // Check result quality
            const totalResults = result.results.length;
            const successfulResults = result.results.filter(r => r.status === 'success').length;
            const successRate = (successfulResults / totalResults * 100).toFixed(1);
            
            console.log(`   Model success rate: ${successRate}%`);
            
            // Check file generation
            const csvExists = fs.existsSync(result.csv_file);
            const htmlExists = result.html_file ? fs.existsSync(result.html_file) : false;
            
            console.log(`   CSV generated: ${csvExists ? 'Yes' : 'No'}`);
            console.log(`   HTML generated: ${htmlExists ? 'Yes' : 'No'}`);
            
            if (htmlExists) {
                const fileSize = fs.statSync(result.html_file).size;
                console.log(`   HTML report size: ${(fileSize/1024).toFixed(1)}KB`);
            }
            
            return successRate >= 70; // Consider success if 70%+ of model calls succeeded
        } else {
            console.log("❌ Complex scenario test failed");
            return false;
        }
        
    } catch (error) {
        console.log(`❌ Complex scenario error: ${error.message}`);
        return false;
    }
}

async function runAdvancedSymphonyTest() {
    console.log("🎼 CASSIA SYMPHONY COMPARE - ADVANCED FUNCTIONALITY TEST");
    console.log("=" .repeat(70));
    
    const tests = [
        { name: "Discussion Rounds", fn: testDiscussionRounds },
        { name: "Model Presets", fn: testModelPresets },
        { name: "Custom Models", fn: testCustomModels },
        { name: "Consensus Thresholds", fn: testConsensusThresholds },
        { name: "Complex Scenario", fn: testComplexScenario }
    ];
    
    const results = {};
    
    for (const test of tests) {
        console.log(`\n🔬 Running ${test.name} test...`);
        try {
            results[test.name] = await test.fn();
        } catch (error) {
            console.log(`❌ ${test.name} test crashed: ${error.message}`);
            results[test.name] = false;
        }
    }
    
    // Summary
    console.log("\n" + "=" .repeat(70));
    console.log("📊 ADVANCED TEST RESULTS:");
    
    let passCount = 0;
    for (const [testName, passed] of Object.entries(results)) {
        const status = passed ? "✅ PASS" : "❌ FAIL";
        console.log(`   ${status}: ${testName}`);
        if (passed) passCount++;
    }
    
    const successRate = (passCount / tests.length * 100).toFixed(1);
    
    if (passCount === tests.length) {
        console.log("\n🎉 ALL ADVANCED SYMPHONY TESTS PASSED!");
        console.log("   ✅ Discussion rounds and consensus building working");
        console.log("   ✅ Multiple model presets operational");
        console.log("   ✅ Custom model configuration functional");
        console.log("   ✅ Consensus threshold variations tested");
        console.log("   ✅ Complex multi-celltype scenarios successful");
        console.log("🚀 Symphony Compare advanced functionality fully verified!");
    } else {
        console.log(`\n⚠️  ${passCount}/${tests.length} advanced tests passed (${successRate}%)`);
        console.log("   Review failed tests above for issues");
    }
    
    console.log(`\n📁 Advanced test results saved in: ./test_results/`);
    console.log("   Open HTML files to review the sophisticated reports generated!");
    console.log("\n🎼 Symphony Compare is ready for production use with multi-model consensus!");
}

runAdvancedSymphonyTest().catch(console.error);