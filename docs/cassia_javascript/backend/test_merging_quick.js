import { mergeAnnotations } from './src/mergingAnnotation.js';
import fs from 'fs';
import path from 'path';

// Create small test data
async function createQuickTestData() {
    const testDir = '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results';
    
    // Ensure test directory exists
    if (!fs.existsSync(testDir)) {
        fs.mkdirSync(testDir, { recursive: true });
    }
    
    const csvPath = path.join(testDir, 'quick_test_clusters.csv');
    
    // Create minimal test CSV data
    const csvContent = `True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types
1,T cell,naive CD4 T cell
2,macrophage,inflammatory macrophage
3,B cell,memory B cell`;
    
    fs.writeFileSync(csvPath, csvContent);
    console.log(`Created quick test data at: ${csvPath}`);
    
    return { csvPath };
}

// Quick test with OpenRouter
async function quickTest() {
    console.log('🧪 Running Quick Merging Test with LLM...\n');
    
    const { csvPath } = await createQuickTestData();
    
    try {
        // Test broad detail level only
        const result = await mergeAnnotations({
            csvPath: csvPath,
            outputPath: '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/quick_merge_result.csv',
            provider: 'openrouter',
            model: 'google/gemini-2.5-flash',
            detailLevel: 'broad',
            batchSize: 3  // Process all in one batch
        });
        
        console.log('✅ Quick test results:');
        result.forEach(row => {
            console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']} → ${row['Merged_Grouping_1']}`);
        });
        
        console.log('\n🎉 Merging agent working correctly!');
        return true;
        
    } catch (error) {
        console.error(`❌ Quick test failed: ${error.message}`);
        return false;
    }
}

// Run quick test
quickTest().catch(console.error);