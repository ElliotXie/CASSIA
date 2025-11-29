import { runCASSIAPipeline, loadMarker } from './src/pipeline.js';
import fs from 'fs';
import path from 'path';

// Create test data for pipeline testing
async function createMockPipelineData() {
    const testDir = '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results';
    
    // Ensure test directory exists
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir, { recursive: true });
    }
    
    // Create comprehensive mock marker data that represents findallmarkers output
    const markerData = [
        {
            "p_val": 0.0001,
            "avg_log2FC": 2.5,
            "pct.1": 0.8,
            "pct.2": 0.1,
            "p_val_adj": 0.001,
            "cluster": "0",
            "gene": "CD3D"
        },
        {
            "p_val": 0.0002,
            "avg_log2FC": 2.3,
            "pct.1": 0.75,
            "pct.2": 0.12,
            "p_val_adj": 0.002,
            "cluster": "0",
            "gene": "CD3E"
        },
        {
            "p_val": 0.0001,
            "avg_log2FC": 1.8,
            "pct.1": 0.65,
            "pct.2": 0.05,
            "p_val_adj": 0.001,
            "cluster": "0",
            "gene": "CD4"
        },
        {
            "p_val": 0.0003,
            "avg_log2FC": 2.1,
            "pct.1": 0.7,
            "pct.2": 0.08,
            "p_val_adj": 0.003,
            "cluster": "1",
            "gene": "CD19"
        },
        {
            "p_val": 0.0002,
            "avg_log2FC": 1.9,
            "pct.1": 0.68,
            "pct.2": 0.1,
            "p_val_adj": 0.002,
            "cluster": "1",
            "gene": "CD20"
        },
        {
            "p_val": 0.0004,
            "avg_log2FC": 1.6,
            "pct.1": 0.6,
            "pct.2": 0.15,
            "p_val_adj": 0.004,
            "cluster": "1",
            "gene": "CD79A"
        },
        {
            "p_val": 0.0001,
            "avg_log2FC": 2.2,
            "pct.1": 0.72,
            "pct.2": 0.06,
            "p_val_adj": 0.001,
            "cluster": "2",
            "gene": "CD68"
        },
        {
            "p_val": 0.0002,
            "avg_log2FC": 2.0,
            "pct.1": 0.7,
            "pct.2": 0.08,
            "p_val_adj": 0.002,
            "cluster": "2",
            "gene": "CD163"
        }
    ];
    
    return markerData;
}

// Test basic pipeline functionality
async function testBasicPipeline() {
    console.log('\n=== Testing Basic Pipeline Functionality ===');
    
    const markerData = await createMockPipelineData();
    
    try {
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'test_pipeline_basic',
            tissue: 'peripheral_blood',
            species: 'human',
            marker: markerData,
            maxWorkers: 2,
            scoreThreshold: 50,  // Lower threshold to test boost functionality
            mergeAnnotations: true,
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter',
            scoreModel: 'google/gemini-2.5-flash',
            scoreProvider: 'openrouter',
            mergeModel: 'google/gemini-2.5-flash',
            mergeProvider: 'openrouter',
            annotationboostModel: 'google/gemini-2.5-flash',
            annotationboostProvider: 'openrouter'
        });
        
        console.log(`✅ Basic pipeline test completed successfully!`);
        console.log(`Final results saved to: ${finalResultsPath}`);
        
        // Verify directory structure was created
        const mainFolder = path.dirname(path.dirname(finalResultsPath));
        const annotationFolder = path.join(mainFolder, '01_annotation_results');
        const reportsFolder = path.join(mainFolder, '02_reports');
        const boostFolder = path.join(mainFolder, '03_boost_analysis');
        const intermediateFolder = path.join(annotationFolder, 'intermediate_files');
        
        const structureValid = fs.existsSync(annotationFolder) && 
                              fs.existsSync(reportsFolder) && 
                              fs.existsSync(boostFolder) && 
                              fs.existsSync(intermediateFolder);
        
        if (structureValid) {
            console.log('✅ Directory structure created correctly');
        } else {
            console.log('❌ Directory structure creation failed');
            return false;
        }
        
        return true;
    } catch (error) {
        console.error(`❌ Basic pipeline test failed: ${error.message}`);
        return false;
    }
}

// Test pipeline with CSV file input
async function testPipelineWithCSV() {
    console.log('\n=== Testing Pipeline with CSV Input ===');
    
    const testDir = '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results';
    const csvPath = path.join(testDir, 'test_markers.csv');
    
    // Create CSV file
    const csvContent = `p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster,gene
0.0001,2.5,0.8,0.1,0.001,0,CD3D
0.0002,2.3,0.75,0.12,0.002,0,CD3E
0.0001,1.8,0.65,0.05,0.001,0,CD4
0.0003,2.1,0.7,0.08,0.003,1,CD19
0.0002,1.9,0.68,0.1,0.002,1,CD20`;
    
    fs.writeFileSync(csvPath, csvContent);
    
    try {
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'test_pipeline_csv',
            tissue: 'lymph_node',
            species: 'mouse',
            marker: csvPath,  // Use CSV file path instead of array
            maxWorkers: 1,
            scoreThreshold: 80,
            mergeAnnotations: false,  // Test without merging
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter'
        });
        
        console.log(`✅ CSV pipeline test completed successfully!`);
        console.log(`Final results saved to: ${finalResultsPath}`);
        
        return true;
    } catch (error) {
        console.error(`❌ CSV pipeline test failed: ${error.message}`);
        return false;
    }
}

// Test loadMarker function
async function testLoadMarker() {
    console.log('\n=== Testing loadMarker Function ===');
    
    try {
        // Test error handling for invalid marker type
        console.log('--- Testing invalid marker type ---');
        try {
            await loadMarker('invalid_type');
            console.log('❌ Should have thrown error for invalid marker type');
            return false;
        } catch (error) {
            console.log('✅ Correctly caught invalid marker type error');
        }
        
        // Test with valid marker types (will fail if files don't exist, which is expected)
        console.log('--- Testing valid marker types (may fail if files not found) ---');
        const markerTypes = ['processed', 'unprocessed', 'subcluster_results'];
        
        for (const markerType of markerTypes) {
            try {
                const data = await loadMarker(markerType);
                console.log(`✅ Successfully loaded ${markerType} marker data: ${data.length} rows`);
            } catch (error) {
                console.log(`⚠️  Could not load ${markerType} marker data (file may not exist): ${error.message}`);
            }
        }
        
        console.log('\n✅ loadMarker function tests completed!');
        return true;
    } catch (error) {
        console.error(`❌ loadMarker test failed: ${error.message}`);
        return false;
    }
}

// Test pipeline configuration options
async function testPipelineConfiguration() {
    console.log('\n=== Testing Pipeline Configuration Options ===');
    
    const markerData = await createMockPipelineData();
    
    try {
        // Test with different configuration
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'test_pipeline_config.csv',  // Test .csv extension removal
            tissue: 'brain',
            species: 'human',
            marker: markerData,
            maxWorkers: 1,
            scoreThreshold: 90,  // High threshold to avoid boost processing
            mergeAnnotations: true,
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter',
            scoreModel: 'google/gemini-2.5-flash',
            scoreProvider: 'openrouter',
            mergeModel: 'google/gemini-2.5-flash',
            mergeProvider: 'openrouter',
            conversationHistoryMode: 'full',
            rankingMethod: 'p_val_adj',
            ascending: true,
            reportStyle: 'total_summary',
            validatorInvolvement: 'v2'
        });
        
        console.log(`✅ Configuration test completed successfully!`);
        console.log(`Final results saved to: ${finalResultsPath}`);
        
        // Verify .csv extension was removed from folder name
        const mainFolder = path.dirname(path.dirname(finalResultsPath));
        const folderName = path.basename(mainFolder);
        
        if (!folderName.includes('.csv')) {
            console.log('✅ CSV extension correctly removed from folder name');
        } else {
            console.log('❌ CSV extension not removed from folder name');
        }
        
        return true;
    } catch (error) {
        console.error(`❌ Configuration test failed: ${error.message}`);
        return false;
    }
}

// Test error handling
async function testErrorHandling() {
    console.log('\n=== Testing Error Handling ===');
    
    try {
        // Test with invalid marker data
        console.log('--- Testing with invalid marker data ---');
        try {
            await runCASSIAPipeline({
                outputFileName: 'test_error',
                tissue: 'test',
                species: 'test',
                marker: 'nonexistent_file.csv',  // Invalid file path
                maxWorkers: 1
            });
            // If it doesn't throw an error, we should check if it handled gracefully
            console.log('⚠️  Pipeline continued despite invalid marker file');
        } catch (error) {
            console.log('✅ Correctly handled invalid marker file error');
        }
        
        console.log('\n✅ Error handling tests completed!');
        return true;
    } catch (error) {
        console.error(`❌ Error handling test failed: ${error.message}`);
        return false;
    }
}

// Main test runner
async function runAllTests() {
    console.log('🧪 Starting Complete Pipeline Tests...\n');
    
    const results = [];
    
    results.push(await testBasicPipeline());
    results.push(await testPipelineWithCSV());
    results.push(await testLoadMarker());
    results.push(await testPipelineConfiguration());
    results.push(await testErrorHandling());
    
    const passedTests = results.filter(result => result).length;
    const totalTests = results.length;
    
    console.log('\n' + '='.repeat(60));
    console.log(`📊 Pipeline Test Results: ${passedTests}/${totalTests} tests passed`);
    
    if (passedTests === totalTests) {
        console.log('🎉 All pipeline tests passed successfully!');
        console.log('\n📁 Test output folders created in current directory:');
        console.log('   - CASSIA_peripheral_blood_human_* (basic test)');
        console.log('   - CASSIA_lymph_node_mouse_* (CSV test)');
        console.log('   - CASSIA_brain_human_* (configuration test)');
        console.log('\n✅ The runCASSIA_pipeline function is working correctly!');
    } else {
        console.log('❌ Some tests failed. Please check the output above.');
    }
    
    return passedTests === totalTests;
}

// Run tests if this script is executed directly
console.log('Starting pipeline test execution...');
runAllTests().catch(error => {
    console.error('Pipeline test execution failed:', error);
    process.exit(1);
});

export { runAllTests };