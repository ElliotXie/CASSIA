import { mergeAnnotations, mergeAnnotationsAll } from './src/mergingAnnotation.js';
import fs from 'fs';
import path from 'path';

// Create test data
async function createMockMergingData() {
    const testDir = '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results';
    
    // Ensure test directory exists
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir, { recursive: true });
    }
    
    const csvPath = path.join(testDir, 'test_cluster_annotations.csv');
    
    // Create test CSV data with cluster annotations
    const csvContent = `True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types,Additional Info
1,T cell,naive CD4 T cell,Helper T cells
2,T cell,cytotoxic CD8 T cell,Killer T cells
3,B cell,memory B cell,Antibody producing
4,macrophage,inflammatory macrophage,M1 type
5,macrophage,tissue-resident macrophage,M2 type
6,dendritic cell,plasmacytoid dendritic cell,Antigen presenting
7,NK cell,CD56bright NK cell,Natural killer
8,B cell,plasma cell,Antibody secreting
9,T cell,regulatory T cell,Immune suppressor
10,monocyte,classical monocyte,Blood derived`;
    
    fs.writeFileSync(csvPath, csvContent);
    console.log(`Created test data at: ${csvPath}`);
    
    return { csvPath };
}

// Test individual mergeAnnotations function
async function testMergeAnnotations() {
    console.log('\n=== Testing mergeAnnotations function ===');
    
    const { csvPath } = await createMockMergingData();
    
    try {
        console.log('\n--- Testing broad detail level ---');
        const broadResult = await mergeAnnotations({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/broad_merge_test.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            detailLevel: 'broad',
            batchSize: 5
        });
        
        console.log('Broad grouping results:');
        broadResult.slice(0, 5).forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']} → ${row['Merged_Grouping_1']}`);
        });
        
        console.log('\n--- Testing detailed detail level ---');
        const detailedResult = await mergeAnnotations({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/detailed_merge_test.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            detailLevel: 'detailed',
            batchSize: 5
        });
        
        console.log('Detailed grouping results:');
        detailedResult.slice(0, 5).forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']} → ${row['Merged_Grouping_2']}`);
        });
        
        console.log('\n--- Testing very_detailed detail level ---');
        const veryDetailedResult = await mergeAnnotations({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/very_detailed_merge_test.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            detailLevel: 'very_detailed',
            batchSize: 5
        });
        
        console.log('Very detailed grouping results:');
        veryDetailedResult.slice(0, 5).forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']} → ${row['Merged_Grouping_3']}`);
        });
        
        console.log('\n✅ mergeAnnotations function test completed successfully!');
        return true;
    } catch (error) {
        console.error(`❌ mergeAnnotations test failed: ${error.message}`);
        return false;
    }
}

// Test parallel mergeAnnotationsAll function
async function testMergeAnnotationsAll() {
    console.log('\n=== Testing mergeAnnotationsAll function ===');
    
    const { csvPath } = await createMockMergingData();
    
    try {
        const allResults = await mergeAnnotationsAll({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/all_merge_test.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            batchSize: 5
        });
        
        console.log('\nParallel processing results (all three detail levels):');
        allResults.slice(0, 5).forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']}`);
            console.log(`  → Broad: ${row['Merged_Grouping_1']}`);
            console.log(`  → Detailed: ${row['Merged_Grouping_2']}`);
            console.log(`  → Very Detailed: ${row['Merged_Grouping_3']}`);
            console.log('');
        });
        
        // Verify all three columns exist
        const hasAllColumns = allResults.every(row => 
            'Merged_Grouping_1' in row && 
            'Merged_Grouping_2' in row && 
            'Merged_Grouping_3' in row
        );
        
        if (hasAllColumns) {
            console.log('✅ All three detail level columns present in results');
        } else {
            console.log('❌ Missing some detail level columns in results');
            return false;
        }
        
        console.log('\n✅ mergeAnnotationsAll function test completed successfully!');
        return true;
    } catch (error) {
        console.error(`❌ mergeAnnotationsAll test failed: ${error.message}`);
        return false;
    }
}

// Test error handling
async function testErrorHandling() {
    console.log('\n=== Testing Error Handling ===');
    
    try {
        // Test invalid detail level
        console.log('--- Testing invalid detail level ---');
        try {
            await mergeAnnotations({
                csvPath: '/fake/path.csv',
                detailLevel: 'invalid_level'
            });
            console.log('❌ Should have thrown error for invalid detail level');
            return false;
        } catch (error) {
            console.log('✅ Correctly caught invalid detail level error');
        }
        
        // Test missing file
        console.log('--- Testing missing file ---');
        try {
            await mergeAnnotations({
                csvPath: '/fake/nonexistent/path.csv',
                detailLevel: 'broad'
            });
            console.log('❌ Should have thrown error for missing file');
            return false;
        } catch (error) {
            console.log('✅ Correctly caught missing file error');
        }
        
        console.log('\n✅ Error handling tests completed successfully!');
        return true;
    } catch (error) {
        console.error(`❌ Error handling test failed: ${error.message}`);
        return false;
    }
}

// Test with additional context
async function testAdditionalContext() {
    console.log('\n=== Testing Additional Context ===');
    
    const { csvPath } = await createMockMergingData();
    
    try {
        const contextResult = await mergeAnnotations({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/context_merge_test.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            detailLevel: 'broad',
            additionalContext: 'This is peripheral blood sample from a healthy human adult. Focus on immune cell classifications.',
            batchSize: 5
        });
        
        console.log('Grouping results with additional context:');
        contextResult.slice(0, 3).forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']} → ${row['Merged_Grouping_1']}`);
        });
        
        console.log('\n✅ Additional context test completed successfully!');
        return true;
    } catch (error) {
        console.error(`❌ Additional context test failed: ${error.message}`);
        return false;
    }
}

// Main test runner
async function runAllTests() {
    console.log('🧪 Starting Merging Annotation Agent Tests...\n');
    
    const results = [];
    
    results.push(await testMergeAnnotations());
    results.push(await testMergeAnnotationsAll());
    results.push(await testErrorHandling());
    results.push(await testAdditionalContext());
    
    const passedTests = results.filter(result => result).length;
    const totalTests = results.length;
    
    console.log('\n' + '='.repeat(50));
    console.log(`📊 Test Results: ${passedTests}/${totalTests} tests passed`);
    
    if (passedTests === totalTests) {
        console.log('🎉 All merging annotation tests passed successfully!');
        console.log('\n📁 Test output files saved to:');
        console.log('   - test_results/broad_merge_test.csv');
        console.log('   - test_results/detailed_merge_test.csv');
        console.log('   - test_results/very_detailed_merge_test.csv');
        console.log('   - test_results/all_merge_test.csv');
        console.log('   - test_results/context_merge_test.csv');
    } else {
        console.log('❌ Some tests failed. Please check the output above.');
    }
    
    return passedTests === totalTests;
}

// Run tests if this script is executed directly
console.log('Starting test execution...');
runAllTests().catch(error => {
    console.error('Test execution failed:', error);
    process.exit(1);
});

export { runAllTests };