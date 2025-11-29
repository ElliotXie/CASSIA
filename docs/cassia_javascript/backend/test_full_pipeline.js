import { runCASSIAPipeline } from './src/pipeline.js';
import fs from 'fs';
import path from 'path';

// Create realistic test data for complete pipeline testing
async function createRealisticTestData() {
    // Create realistic marker data that mimics findallmarkers output
    const markerData = [
        // T cell cluster
        { "p_val": 0.0001, "avg_log2FC": 2.5, "pct.1": 0.8, "pct.2": 0.1, "p_val_adj": 0.001, "cluster": "0", "gene": "CD3D" },
        { "p_val": 0.0002, "avg_log2FC": 2.3, "pct.1": 0.75, "pct.2": 0.12, "p_val_adj": 0.002, "cluster": "0", "gene": "CD3E" },
        { "p_val": 0.0001, "avg_log2FC": 1.8, "pct.1": 0.65, "pct.2": 0.05, "p_val_adj": 0.001, "cluster": "0", "gene": "CD4" },
        
        // B cell cluster
        { "p_val": 0.0003, "avg_log2FC": 2.1, "pct.1": 0.7, "pct.2": 0.08, "p_val_adj": 0.003, "cluster": "1", "gene": "CD19" },
        { "p_val": 0.0002, "avg_log2FC": 1.9, "pct.1": 0.68, "pct.2": 0.1, "p_val_adj": 0.002, "cluster": "1", "gene": "CD20" },
        { "p_val": 0.0004, "avg_log2FC": 1.6, "pct.1": 0.6, "pct.2": 0.15, "p_val_adj": 0.004, "cluster": "1", "gene": "CD79A" },
        
        // Macrophage cluster
        { "p_val": 0.0001, "avg_log2FC": 2.2, "pct.1": 0.72, "pct.2": 0.06, "p_val_adj": 0.001, "cluster": "2", "gene": "CD68" },
        { "p_val": 0.0002, "avg_log2FC": 2.0, "pct.1": 0.7, "pct.2": 0.08, "p_val_adj": 0.002, "cluster": "2", "gene": "CD163" },
        { "p_val": 0.0003, "avg_log2FC": 1.7, "pct.1": 0.65, "pct.2": 0.1, "p_val_adj": 0.003, "cluster": "2", "gene": "CD14" }
    ];
    
    return markerData;
}

// Test the complete pipeline with actual LLM calls
async function testCompletePipeline() {
    console.log('🧪 Testing Complete Pipeline with Real LLM Calls...\n');
    
    const markerData = await createRealisticTestData();
    
    try {
        console.log('Starting complete pipeline test...');
        const startTime = Date.now();
        
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'full_pipeline_test',
            tissue: 'peripheral_blood',
            species: 'human',
            marker: markerData,
            maxWorkers: 2,
            scoreThreshold: 60,  // Lower threshold to potentially trigger boost analysis
            mergeAnnotations: true,
            
            // Use fast, reliable models for testing
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter',
            scoreModel: 'google/gemini-2.5-flash', 
            scoreProvider: 'openrouter',
            mergeModel: 'google/gemini-2.5-flash',
            mergeProvider: 'openrouter',
            annotationboostModel: 'google/gemini-2.5-flash',
            annotationboostProvider: 'openrouter',
            
            // Other settings
            conversationHistoryMode: 'final',
            rankingMethod: 'avg_log2FC',
            reportStyle: 'per_iteration',
            validatorInvolvement: 'v1',
            additionalInfo: 'Test dataset for pipeline validation'
        });
        
        const endTime = Date.now();
        const duration = (endTime - startTime) / 1000;
        
        console.log(`\n✅ Complete pipeline test finished in ${duration.toFixed(1)} seconds!`);
        console.log(`📊 Final results: ${finalResultsPath}`);
        
        // Verify the output structure
        await verifyPipelineOutput(finalResultsPath);
        
        return true;
        
    } catch (error) {
        console.error(`❌ Complete pipeline test failed: ${error.message}`);
        console.error('Stack trace:', error.stack);
        return false;
    }
}

// Verify the pipeline output structure and content
async function verifyPipelineOutput(finalResultsPath) {
    console.log('\n🔍 Verifying Pipeline Output...');
    
    try {
        // Check if final results file exists
        if (!fs.existsSync(finalResultsPath)) {
            throw new Error(`Final results file not found: ${finalResultsPath}`);
        }
        console.log('✅ Final results file exists');
        
        // Verify directory structure
        const mainFolder = path.dirname(path.dirname(finalResultsPath));
        const expectedFolders = [
            path.join(mainFolder, '01_annotation_results'),
            path.join(mainFolder, '02_reports'),
            path.join(mainFolder, '03_boost_analysis'),
            path.join(mainFolder, '01_annotation_results', 'intermediate_files')
        ];
        
        for (const folder of expectedFolders) {
            if (!fs.existsSync(folder)) {
                throw new Error(`Expected folder not found: ${folder}`);
            }
            console.log(`✅ Found folder: ${path.basename(folder)}`);
        }
        
        // Check final results content
        const csvContent = fs.readFileSync(finalResultsPath, 'utf-8');
        const lines = csvContent.split('\n').filter(line => line.trim());
        
        if (lines.length < 2) {
            throw new Error('Final results file appears empty or invalid');
        }
        console.log(`✅ Final results has ${lines.length - 1} data rows`);
        
        // Check for expected columns
        const header = lines[0];
        const expectedColumns = ['True Cell Type', 'Predicted Main Cell Type', 'Score'];
        
        for (const col of expectedColumns) {
            if (!header.includes(col)) {
                throw new Error(`Missing expected column: ${col}`);
            }
        }
        console.log('✅ All expected columns present');
        
        // Check reports folder
        const reportsFolder = path.join(mainFolder, '02_reports');
        const reportFiles = fs.readdirSync(reportsFolder);
        const htmlFiles = reportFiles.filter(f => f.endsWith('.html'));
        
        if (htmlFiles.length === 0) {
            console.log('⚠️  No HTML reports found (may be expected if report generation failed)');
        } else {
            console.log(`✅ Found ${htmlFiles.length} HTML report(s)`);
        }
        
        // Check intermediate files
        const intermediateFolder = path.join(mainFolder, '01_annotation_results', 'intermediate_files');
        const intermediateFiles = fs.readdirSync(intermediateFolder);
        
        if (intermediateFiles.length === 0) {
            console.log('⚠️  No intermediate files found');
        } else {
            console.log(`✅ Found ${intermediateFiles.length} intermediate file(s)`);
        }
        
        console.log('\n🎉 Pipeline output verification completed successfully!');
        
    } catch (error) {
        console.error(`❌ Output verification failed: ${error.message}`);
        throw error;
    }
}

// Test with minimal configuration
async function testMinimalPipeline() {
    console.log('\n🧪 Testing Minimal Pipeline Configuration...\n');
    
    const markerData = await createRealisticTestData();
    
    try {
        const finalResultsPath = await runCASSIAPipeline({
            outputFileName: 'minimal_pipeline_test',
            tissue: 'test_tissue',
            species: 'human',
            marker: markerData.slice(0, 3), // Only first 3 markers for speed
            maxWorkers: 1,
            scoreThreshold: 90, // High threshold to avoid boost analysis
            mergeAnnotations: false, // Disable merging for speed
            annotationModel: 'google/gemini-2.5-flash',
            annotationProvider: 'openrouter'
        });
        
        console.log(`✅ Minimal pipeline test completed!`);
        console.log(`📊 Results: ${finalResultsPath}`);
        
        return true;
        
    } catch (error) {
        console.error(`❌ Minimal pipeline test failed: ${error.message}`);
        return false;
    }
}

// Main test runner
async function runPipelineTests() {
    console.log('🚀 Starting Full Pipeline Testing...\n');
    
    const results = [];
    
    // Test minimal configuration first
    results.push(await testMinimalPipeline());
    
    // Test complete pipeline if minimal works
    if (results[0]) {
        results.push(await testCompletePipeline());
    } else {
        console.log('⚠️  Skipping complete pipeline test due to minimal test failure');
        results.push(false);
    }
    
    const passedTests = results.filter(result => result).length;
    const totalTests = results.length;
    
    console.log('\n' + '='.repeat(60));
    console.log(`📊 Pipeline Test Results: ${passedTests}/${totalTests} tests passed`);
    
    if (passedTests === totalTests) {
        console.log('🎉 All pipeline tests passed! The pipeline is working correctly!');
        console.log('\n📁 Check the generated CASSIA_* folders for output');
    } else {
        console.log('❌ Some pipeline tests failed. Check the output above for details.');
    }
    
    return passedTests === totalTests;
}

// Run the tests
runPipelineTests().catch(error => {
    console.error('Pipeline testing failed:', error);
    process.exit(1);
});