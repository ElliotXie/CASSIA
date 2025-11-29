#!/usr/bin/env node

// Test script for CASSIA Batch Integration
import { runCASSIABatch } from './backend/src/runCASSIA_batch.js';
import fs from 'fs/promises';
import path from 'path';

console.log('🧪 Testing CASSIA Batch Integration...\n');

async function testBatchFunctionality() {
    try {
        // Test data - simple marker format
        const testData = [
            { cluster: '0', markers: 'CD3D, CD3E, CD4' },
            { cluster: '1', markers: 'CD19, CD20, CD79A' },
            { cluster: '2', markers: 'CD68, CD163, CD14' }
        ];

        console.log('📊 Test Data:');
        console.table(testData);

        console.log('\n🚀 Running Batch Analysis...');
        
        const results = await runCASSIABatch({
            marker: testData,
            outputName: 'test_batch_output',
            nGenes: 10,
            model: 'google/gemini-2.5-flash-preview',
            temperature: 0,
            tissue: 'peripheral_blood',
            species: 'human',
            maxWorkers: 2,
            provider: 'openrouter',
            maxRetries: 0 // Skip retries for test
        });

        console.log('\n✅ Batch Analysis Results:');
        console.log(`📈 Total Clusters: ${results.total_clusters}`);
        console.log(`✅ Successful Analyses: ${results.successful_analyses}`);
        console.log(`❌ Failed Analyses: ${results.failed_analyses}`);
        console.log(`📁 Output Files: ${results.output_files.join(', ')}`);

        // Check if output files were created
        for (const file of results.output_files) {
            try {
                const stats = await fs.stat(file);
                console.log(`📄 ${path.basename(file)}: ${(stats.size / 1024).toFixed(2)} KB`);
            } catch (error) {
                console.log(`❌ ${path.basename(file)}: File not found`);
            }
        }

        console.log('\n🎯 Batch functionality test completed successfully!');
        return true;

    } catch (error) {
        console.error('\n❌ Batch test failed:', error.message);
        console.error(error.stack);
        return false;
    }
}

// Run the test
testBatchFunctionality()
    .then(success => {
        process.exit(success ? 0 : 1);
    })
    .catch(error => {
        console.error('Unexpected error:', error);
        process.exit(1);
    });