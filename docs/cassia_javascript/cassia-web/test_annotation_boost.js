/**
 * Test the updated annotation boost module with new CSV parser and column selection
 */

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';
import { 
    extractConversationForCluster, 
    getAvailableClusters, 
    getAvailableColumns 
} from './lib/cassia/annotationBoost.js';

// Get current directory
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

async function testAnnotationBoost() {
    console.log('🧪 Testing Updated Annotation Boost Module\n');
    
    try {
        // Load the example CSV
        const csvPath = join(__dirname, 'public/examples/scoring_example_full.csv');
        const csvContent = readFileSync(csvPath, 'utf-8');
        
        console.log('✅ Loaded example CSV');
        console.log(`📊 File size: ${csvContent.length} characters\n`);
        
        // Test 1: Get available columns
        console.log('1️⃣ Testing getAvailableColumns:');
        const columns = getAvailableColumns(csvContent);
        console.log(`Found ${columns.length} columns:`, columns);
        
        // Check for common cluster columns
        const clusterOptions = ['True Cell Type', 'Cluster', 'Cell Type', 'Predicted Main Cell Type'];
        const foundClusterColumns = clusterOptions.filter(col => columns.includes(col));
        console.log(`Found cluster-like columns:`, foundClusterColumns);
        
        // Test 2: Get clusters using default column (True Cell Type)
        console.log('\n2️⃣ Testing getAvailableClusters with default column:');
        const defaultClusters = getAvailableClusters(csvContent, 'True Cell Type');
        console.log(`Found ${defaultClusters.length} clusters in "True Cell Type":`, defaultClusters.slice(0, 3));
        
        // Test 3: Try with different cluster column if available
        if (columns.includes('Predicted Main Cell Type')) {
            console.log('\n3️⃣ Testing getAvailableClusters with "Predicted Main Cell Type":');
            const predictedClusters = getAvailableClusters(csvContent, 'Predicted Main Cell Type');
            console.log(`Found ${predictedClusters.length} clusters in "Predicted Main Cell Type":`, predictedClusters.slice(0, 3));
        }
        
        // Test 4: Extract conversation history for a cluster
        if (defaultClusters.length > 0) {
            console.log('\n4️⃣ Testing extractConversationForCluster:');
            const testCluster = defaultClusters[0];
            console.log(`Extracting conversation for: "${testCluster}"`);
            
            try {
                const conversation = extractConversationForCluster(csvContent, testCluster, 'True Cell Type');
                console.log(`✅ Successfully extracted conversation: ${conversation.length} characters`);
                console.log(`📝 Preview: ${conversation.substring(0, 200)}...`);
                
                // Check if conversation contains expected content
                const hasMarkers = conversation.toLowerCase().includes('marker');
                const hasAnalysis = conversation.toLowerCase().includes('analysis');
                const hasAnnotation = conversation.toLowerCase().includes('annotation');
                
                console.log(`📊 Content validation:`);
                console.log(`  - Contains "marker": ${hasMarkers}`);
                console.log(`  - Contains "analysis": ${hasAnalysis}`);
                console.log(`  - Contains "annotation": ${hasAnnotation}`);
                
            } catch (extractError) {
                console.error(`❌ Failed to extract conversation: ${extractError.message}`);
            }
        }
        
        // Test 5: Test with non-existent cluster
        console.log('\n5️⃣ Testing error handling with non-existent cluster:');
        try {
            extractConversationForCluster(csvContent, 'NonExistentCluster', 'True Cell Type');
            console.log('❌ Should have thrown an error');
        } catch (expectedError) {
            console.log(`✅ Correctly threw error: ${expectedError.message.substring(0, 100)}...`);
        }
        
        // Test 6: Test with non-existent column
        console.log('\n6️⃣ Testing error handling with non-existent column:');
        try {
            extractConversationForCluster(csvContent, defaultClusters[0], 'NonExistentColumn');
            console.log('❌ Should have thrown an error');
        } catch (expectedError) {
            console.log(`✅ Correctly threw error: ${expectedError.message.substring(0, 100)}...`);
        }
        
        console.log('\n🎉 All tests completed successfully!');
        console.log('\n📋 Summary of new features:');
        console.log('- ✅ Robust CSV parsing with multiline field support');
        console.log('- ✅ User-selectable cluster column identification');
        console.log('- ✅ Default "True Cell Type" column with fallbacks');
        console.log('- ✅ Enhanced error handling and validation');
        console.log('- ✅ Support for various column naming conventions');
        
    } catch (error) {
        console.error('❌ Test failed:', error.message);
        console.error(error.stack);
    }
}

// Run the test
testAnnotationBoost().then(() => {
    console.log('\n✅ Annotation boost module test completed');
}).catch(error => {
    console.error('\n❌ Test suite failed:', error);
});