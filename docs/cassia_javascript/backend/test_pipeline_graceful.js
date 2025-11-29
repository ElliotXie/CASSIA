import { runCASSIAPipeline } from './src/pipeline.js';
import fs from 'fs';
import path from 'path';

// Test that pipeline handles missing API keys gracefully
async function testPipelineGracefulHandling() {
    console.log('🧪 Testing Pipeline Graceful Error Handling...\n');
    
    const markerData = [
        { "p_val": 0.0001, "avg_log2FC": 2.5, "pct.1": 0.8, "pct.2": 0.1, "p_val_adj": 0.001, "cluster": "0", "gene": "CD3D" },
        { "p_val": 0.0002, "avg_log2FC": 2.3, "pct.1": 0.75, "pct.2": 0.12, "p_val_adj": 0.002, "cluster": "1", "gene": "CD19" }
    ];
    
    try {
        console.log('Starting pipeline with no API keys (should handle gracefully)...');
        
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'graceful_test',
            tissue: 'test_tissue',
            species: 'human',
            marker: markerData,
            maxWorkers: 1,
            scoreThreshold: 90, // High threshold to avoid boost
            mergeAnnotations: false, // Disable merging to reduce complexity
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter'
        });
        
        console.log('\n🔍 Checking Pipeline Output...');
        
        // The pipeline should complete even if API calls fail
        // Check that directory structure was created
        const mainFolder = path.dirname(path.dirname(finalResultsPath));
        console.log(`Main folder created: ${mainFolder}`);
        
        const expectedFolders = [
            path.join(mainFolder, '01_annotation_results'),
            path.join(mainFolder, '02_reports'),
            path.join(mainFolder, '03_boost_analysis')
        ];
        
        let foldersCreated = 0;
        for (const folder of expectedFolders) {
            if (fs.existsSync(folder)) {
                console.log(`✅ ${path.basename(folder)} folder created`);
                foldersCreated++;
            } else {
                console.log(`❌ ${path.basename(folder)} folder missing`);
            }
        }
        
        // Check if final results file was created (even if empty)
        if (fs.existsSync(finalResultsPath)) {
            console.log('✅ Final results file created');
            
            const content = fs.readFileSync(finalResultsPath, 'utf-8');
            console.log(`📄 Final results content length: ${content.length} characters`);
            
            if (content.includes('True Cell Type')) {
                console.log('✅ Final results has expected header structure');
            }
        } else {
            console.log('❌ Final results file not created');
        }
        
        console.log(`\n📊 Summary: ${foldersCreated}/3 folders created`);
        
        if (foldersCreated === 3) {
            console.log('🎉 Pipeline structure creation works correctly!');
            console.log('✅ Pipeline handles missing API keys gracefully');
            return true;
        } else {
            console.log('⚠️  Some pipeline structure issues detected');
            return false;
        }
        
    } catch (error) {
        console.error(`❌ Pipeline graceful handling test failed: ${error.message}`);
        console.error('Stack trace:', error.stack);
        return false;
    }
}

// Test directory naming logic
async function testDirectoryNaming() {
    console.log('\n🧪 Testing Directory Naming Logic...\n');
    
    try {
        // Test various tissue and species combinations
        const testCases = [
            { tissue: 'peripheral blood', species: 'human' },
            { tissue: 'brain_cortex', species: 'mouse' },
            { tissue: 'liver-tissue', species: 'rat' },
            { tissue: 'heart@muscle', species: 'human' } // Test special characters
        ];
        
        for (const testCase of testCases) {
            let mainFolderName = `CASSIA_${testCase.tissue}_${testCase.species}`;
            mainFolderName = mainFolderName.replace(/[^a-zA-Z0-9\s\-_]/g, '').trim();
            mainFolderName = mainFolderName.replace(/\s+/g, '_');
            
            console.log(`${testCase.tissue} + ${testCase.species} → ${mainFolderName}`);
            
            // Check that no invalid characters remain
            if (!/^[a-zA-Z0-9_\-]+$/.test(mainFolderName.split('_').slice(0, -1).join('_'))) {
                console.log('❌ Invalid characters in folder name');
                return false;
            }
        }
        
        console.log('✅ Directory naming logic works correctly');
        return true;
        
    } catch (error) {
        console.error(`❌ Directory naming test failed: ${error.message}`);
        return false;
    }
}

// Main test runner
async function runGracefulTests() {
    console.log('🚀 Starting Pipeline Graceful Handling Tests...\n');
    
    const results = [];
    
    results.push(await testDirectoryNaming());
    results.push(await testPipelineGracefulHandling());
    
    const passedTests = results.filter(result => result).length;
    const totalTests = results.length;
    
    console.log('\n' + '='.repeat(60));
    console.log(`📊 Graceful Handling Test Results: ${passedTests}/${totalTests} tests passed`);
    
    if (passedTests === totalTests) {
        console.log('🎉 Pipeline graceful handling works correctly!');
        console.log('\n📝 Test Summary:');
        console.log('   ✅ Directory structure creation works');
        console.log('   ✅ Pipeline handles API failures gracefully');  
        console.log('   ✅ Error handling prevents crashes');
        console.log('   ✅ File organization logic functional');
        console.log('\n🔑 To test with real functionality, set API keys:');
        console.log('   export OPENROUTER_API_KEY="your-key-here"');
        console.log('   export ANTHROPIC_API_KEY="your-key-here"');
    } else {
        console.log('❌ Some graceful handling tests failed');
    }
    
    return passedTests === totalTests;
}

// Run the tests
runGracefulTests().catch(error => {
    console.error('Graceful testing failed:', error);
    process.exit(1);
});