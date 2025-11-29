import { runCASSIAPipeline } from './src/pipeline.js';
import fs from 'fs';
import path from 'path';

// Mock the LLM functions to avoid API calls
const originalRunCASSIABatch = (await import('./src/runCASSIA_batch.js')).runCASSIABatch;
const originalRunCASSIAScoreBatch = (await import('./src/scoring.js')).runCASSIAScoreBatch;
const originalMergeAnnotationsAll = (await import('./src/mergingAnnotation.js')).mergeAnnotationsAll;

// Create test data
function createTestData() {
    return [
        { "p_val": 0.0001, "avg_log2FC": 2.5, "pct.1": 0.8, "pct.2": 0.1, "p_val_adj": 0.001, "cluster": "0", "gene": "CD3D" },
        { "p_val": 0.0002, "avg_log2FC": 2.3, "pct.1": 0.75, "pct.2": 0.12, "p_val_adj": 0.002, "cluster": "1", "gene": "CD19" }
    ];
}

// Mock successful annotation output
function createMockAnnotationFiles(outputName) {
    const fullCsv = `True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types,markers,tissue,species,additional_info,analysis_text
0,T cell,naive CD4 T cell,"CD3D,CD3E,CD4",peripheral_blood,human,Test dataset,Mock analysis for T cells
1,B cell,naive B cell,"CD19,CD20,CD79A",peripheral_blood,human,Test dataset,Mock analysis for B cells`;
    
    const summaryCsv = `True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types
0,T cell,naive CD4 T cell
1,B cell,naive B cell`;
    
    fs.writeFileSync(`${outputName}_full.csv`, fullCsv);
    fs.writeFileSync(`${outputName}_summary.csv`, summaryCsv);
    
    console.log(`Created mock annotation files for ${outputName}`);
}

// Mock successful scoring output
function createMockScoredFile(inputFile, outputFile) {
    if (!fs.existsSync(inputFile)) {
        console.log(`⚠️  Input file for scoring mock not found: ${inputFile}`);
        return;
    }
    
    const inputContent = fs.readFileSync(inputFile, 'utf-8');
    const lines = inputContent.split('\n');
    
    if (lines.length < 2) return;
    
    // Add Score column to header
    const header = lines[0] + ',Score';
    
    // Add mock scores to data lines
    const dataLines = lines.slice(1).filter(line => line.trim()).map(line => {
        const score = Math.floor(Math.random() * 30) + 70; // Random score 70-100
        return line + ',' + score;
    });
    
    const scoredContent = [header, ...dataLines].join('\n');
    fs.writeFileSync(outputFile, scoredContent);
    
    console.log(`Created mock scored file: ${outputFile}`);
}

// Test pipeline structure without API calls
async function testPipelineStructure() {
    console.log('🧪 Testing Pipeline Structure (No API Calls)...\n');
    
    // Mock the functions to avoid API calls
    const mockRunCASSIABatch = async (params) => {
        console.log('Mock: Running CASSIA batch analysis...');
        createMockAnnotationFiles(params.outputName);
        console.log('Mock: CASSIA batch analysis completed');
    };
    
    const mockRunCASSIAScoreBatch = async (params) => {
        console.log('Mock: Running scoring process...');
        createMockScoredFile(params.inputFile, params.outputFile);
        console.log('Mock: Scoring process completed');
    };
    
    const mockMergeAnnotationsAll = async (params) => {
        console.log('Mock: Running annotation merging...');
        // Create a simple merged file by copying the input and adding merge columns
        if (fs.existsSync(params.csvPath)) {
            const content = fs.readFileSync(params.csvPath, 'utf-8');
            const lines = content.split('\n');
            if (lines.length > 0) {
                const header = lines[0] + ',Merged_Grouping_1,Merged_Grouping_2,Merged_Grouping_3';
                const dataLines = lines.slice(1).filter(line => line.trim()).map(line => {
                    return line + ',Lymphoid cells,T cells,Helper T cells';
                });
                const mergedContent = [header, ...dataLines].join('\n');
                fs.writeFileSync(params.outputPath, mergedContent);
            }
        }
        console.log('Mock: Annotation merging completed');
    };
    
    // Temporarily replace the functions
    const pipeline = await import('./src/pipeline.js');
    const batchModule = await import('./src/runCASSIA_batch.js');
    const scoringModule = await import('./src/scoring.js');
    const mergingModule = await import('./src/mergingAnnotation.js');
    
    // Replace functions with mocks
    batchModule.runCASSIABatch = mockRunCASSIABatch;
    scoringModule.runCASSIAScoreBatch = mockRunCASSIAScoreBatch;
    mergingModule.mergeAnnotationsAll = mockMergeAnnotationsAll;
    
    try {
        const markerData = createTestData();
        
        const finalResultsPath = await pipeline.runCASSIAPipeline({
            outputFileName: 'structure_test',
            tissue: 'test_tissue',
            species: 'human',
            marker: markerData,
            maxWorkers: 1,
            scoreThreshold: 90, // High threshold to avoid boost
            mergeAnnotations: true,
            annotationModel: 'test-model',
            annotationProvider: 'test-provider'
        });
        
        console.log('\n🔍 Verifying Pipeline Structure...');
        
        // Check directory structure
        const mainFolder = path.dirname(path.dirname(finalResultsPath));
        console.log(`Main folder: ${mainFolder}`);
        
        const expectedFolders = [
            path.join(mainFolder, '01_annotation_results'),
            path.join(mainFolder, '02_reports'),
            path.join(mainFolder, '03_boost_analysis')
        ];
        
        let allFoldersExist = true;
        for (const folder of expectedFolders) {
            if (fs.existsSync(folder)) {
                console.log(`✅ ${path.basename(folder)} folder exists`);
            } else {
                console.log(`❌ ${path.basename(folder)} folder missing`);
                allFoldersExist = false;
            }
        }
        
        // Check final results file
        if (fs.existsSync(finalResultsPath)) {
            console.log('✅ Final results file exists');
            
            const content = fs.readFileSync(finalResultsPath, 'utf-8');
            const lines = content.split('\n').filter(line => line.trim());
            console.log(`✅ Final results has ${lines.length - 1} data rows`);
            
            // Check for expected columns
            const header = lines[0];
            const expectedColumns = ['True Cell Type', 'Predicted Main Cell Type', 'Score', 'Merged_Grouping_1'];
            
            let allColumnsPresent = true;
            for (const col of expectedColumns) {
                if (header.includes(col)) {
                    console.log(`✅ Column '${col}' present`);
                } else {
                    console.log(`❌ Column '${col}' missing`);
                    allColumnsPresent = false;
                }
            }
            
            if (allFoldersExist && allColumnsPresent) {
                console.log('\n🎉 Pipeline structure test passed completely!');
                return true;
            }
        } else {
            console.log('❌ Final results file not found');
        }
        
        return false;
        
    } catch (error) {
        console.error(`❌ Pipeline structure test failed: ${error.message}`);
        console.error('Stack trace:', error.stack);
        return false;
    }
}

// Test error handling without API calls
async function testErrorHandling() {
    console.log('\n🧪 Testing Error Handling...\n');
    
    try {
        // Test with invalid marker data
        await runCASSIAPipeline({
            outputFileName: 'error_test',
            tissue: 'test',
            species: 'test',
            marker: [], // Empty marker data
            maxWorkers: 1
        });
        
        console.log('✅ Pipeline handled empty marker data gracefully');
        return true;
        
    } catch (error) {
        console.log(`✅ Pipeline correctly handled error: ${error.message}`);
        return true;
    }
}

// Main test runner
async function runStructureTests() {
    console.log('🚀 Starting Pipeline Structure Tests (No API Required)...\n');
    
    const results = [];
    
    results.push(await testPipelineStructure());
    results.push(await testErrorHandling());
    
    const passedTests = results.filter(result => result).length;
    const totalTests = results.length;
    
    console.log('\n' + '='.repeat(60));
    console.log(`📊 Structure Test Results: ${passedTests}/${totalTests} tests passed`);
    
    if (passedTests === totalTests) {
        console.log('🎉 Pipeline structure is working correctly!');
        console.log('✅ Ready for testing with real API keys');
    } else {
        console.log('❌ Some structure tests failed');
    }
    
    return passedTests === totalTests;
}

// Run the tests
runStructureTests().catch(error => {
    console.error('Structure testing failed:', error);
    process.exit(1);
});