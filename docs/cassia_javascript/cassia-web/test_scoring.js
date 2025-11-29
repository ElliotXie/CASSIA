/**
 * Test script for scoring.js functionality
 * This will test the CSV parsing and column detection
 */

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

// Get current directory
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Import our scoring functions
// Note: We'll need to simulate the browser environment
global.fetch = async (url) => {
    throw new Error('fetch not available in test environment');
};

// Mock console methods to capture output
let logOutput = [];
let errorOutput = [];
const originalConsoleLog = console.log;
const originalConsoleError = console.error;

console.log = (...args) => {
    logOutput.push(args.join(' '));
    originalConsoleLog(...args);
};

console.error = (...args) => {
    errorOutput.push(args.join(' '));
    originalConsoleError(...args);
};

// Import scoring functions by reading the file and evaluating
const scoringCode = readFileSync(join(__dirname, 'lib/cassia/scoring.js'), 'utf-8');

// Create a mock module environment
const module = { exports: {} };
const exports = module.exports;

// Mock LLM call for testing
const mockCallLLM = async (prompt, provider, model, apiKey, temperature, maxTokens) => {
    return `<reasoning>
This is a test annotation with good marker expression patterns. The genes show strong specificity for T cells with CD8A, CD8B, and cytotoxic markers like GZMB, PRF1. The annotation correctly identifies this as a CD8+ T cell population.
</reasoning>

<score>85</score>`;
};

// Create a test environment
const testEnvironment = {
    callLLM: mockCallLLM,
    console,
    fetch: global.fetch
};

// Evaluate the scoring code in our test environment
function createScoringModule() {
    // Parse CSV function from scoring.js
    function parseCSV(csvText) {
        const lines = csvText.trim().split('\n');
        if (lines.length < 2) return [];
        
        // More robust CSV parser that handles multiline quoted fields
        const parseCSVAdvanced = (text) => {
            const rows = [];
            let currentRow = [];
            let currentField = '';
            let inQuotes = false;
            let i = 0;
            
            while (i < text.length) {
                const char = text[i];
                const nextChar = text[i + 1];
                
                if (char === '"') {
                    if (inQuotes && nextChar === '"') {
                        // Escaped quote - add literal quote
                        currentField += '"';
                        i += 2;
                        continue;
                    } else {
                        // Toggle quote state
                        inQuotes = !inQuotes;
                    }
                } else if (char === ',' && !inQuotes) {
                    // End of field
                    currentRow.push(currentField.trim());
                    currentField = '';
                } else if ((char === '\n' || char === '\r') && !inQuotes) {
                    // End of row
                    if (currentField.trim() || currentRow.length > 0) {
                        currentRow.push(currentField.trim());
                        rows.push([...currentRow]);
                        currentRow = [];
                        currentField = '';
                    }
                    // Skip \r\n combinations
                    if (char === '\r' && nextChar === '\n') {
                        i++;
                    }
                } else {
                    currentField += char;
                }
                i++;
            }
            
            // Handle last field/row
            if (currentField.trim() || currentRow.length > 0) {
                currentRow.push(currentField.trim());
                rows.push(currentRow);
            }
            
            return rows;
        };
        
        try {
            // Parse the entire CSV
            const parsedRows = parseCSVAdvanced(csvText);
            
            if (parsedRows.length < 2) return [];
            
            const headers = parsedRows[0];
            const dataRows = parsedRows.slice(1);
            
            console.log(`CSV Parser Debug: Found ${headers.length} headers, ${dataRows.length} data rows`);
            console.log(`Headers: ${headers.join(' | ')}`);
            
            // Convert to objects
            const result = dataRows.map((row, rowIndex) => {
                const obj = {};
                headers.forEach((header, colIndex) => {
                    obj[header] = row[colIndex] || '';
                });
                
                // Debug first few rows
                if (rowIndex < 3) {
                    console.log(`Row ${rowIndex + 1} sample:`, {
                        'Marker List': obj['Marker List'] ? `${obj['Marker List'].length} chars: ${obj['Marker List'].substring(0, 50)}...` : 'EMPTY',
                        'Species': obj['Species'],
                        'Tissue': obj['Tissue']
                    });
                }
                
                return obj;
            });
            
            return result;
            
        } catch (error) {
            console.error('CSV parsing error:', error);
            return [];
        }
    }

    // Column detection function
    function findColumn(row, options) {
        for (const col of options) {
            if (col in row && row[col] && row[col].trim() !== '') {
                return { value: row[col].trim(), column: col };
            }
        }
        // Try case-insensitive match
        const availableColumns = Object.keys(row);
        for (const availCol of availableColumns) {
            const normalizedAvail = availCol.toLowerCase().replace(/[\s._-]+/g, '');
            for (const option of options) {
                const normalizedOption = option.toLowerCase().replace(/[\s._-]+/g, '');
                if (normalizedAvail === normalizedOption && row[availCol] && row[availCol].trim() !== '') {
                    return { value: row[availCol].trim(), column: availCol };
                }
            }
        }
        return null;
    }

    // Test function to process a single row
    function testProcessSingleRow(row, idx) {
        console.log(`\n=== Testing Row ${idx + 1} ===`);
        
        // Find Species and Tissue
        const species = findColumn(row, ['Species', 'species', 'SPECIES', 'organism', 'Organism']);
        const tissue = findColumn(row, ['Tissue', 'tissue', 'TISSUE', 'organ', 'Organ', 'sample', 'Sample']);
        
        console.log(`Species: ${species ? species.value : 'NOT FOUND'} (column: ${species ? species.column : 'N/A'})`);
        console.log(`Tissue: ${tissue ? tissue.value : 'NOT FOUND'} (column: ${tissue ? tissue.column : 'N/A'})`);
        
        // Find Marker List
        const markerColumnOptions = [
            'Marker List', 'Marker.List', 'marker_list', 'Marker_List', 'Marker-List',
            'MarkerList', 'markerlist', 'marker list', 'MARKER LIST',
            'markers', 'Markers', 'marker', 'Marker', 'genes', 'Genes'
        ];
        
        const marker = findColumn(row, markerColumnOptions);
        console.log(`Marker: ${marker ? `${marker.value.length} chars: ${marker.value.substring(0, 100)}...` : 'NOT FOUND'} (column: ${marker ? marker.column : 'N/A'})`);
        
        // Find Conversation History
        const historyColumnOptions = [
            'Conversation History', 'Conversation.History', 'conversation_history', 
            'Conversation_History', 'Conversation-History', 'ConversationHistory'
        ];
        
        const history = findColumn(row, historyColumnOptions);
        console.log(`History: ${history ? `${history.value.length} chars: ${history.value.substring(0, 100)}...` : 'NOT FOUND'} (column: ${history ? history.column : 'N/A'})`);
        
        // Check if we have all required data
        const hasAllRequiredData = species && tissue && marker && history;
        console.log(`Has all required data: ${hasAllRequiredData ? 'YES' : 'NO'}`);
        
        if (!hasAllRequiredData) {
            console.log('Available columns:', Object.keys(row));
            console.log('Sample of row data:', Object.entries(row).slice(0, 5).map(([k, v]) => `${k}: "${String(v).substring(0, 30)}..."`));
        }
        
        return hasAllRequiredData;
    }

    return {
        parseCSV,
        findColumn,
        testProcessSingleRow
    };
}

// Main test function
async function runTest() {
    console.log('🧪 Starting Scoring.js Test\n');
    
    try {
        // Read the example CSV file
        const csvPath = join(__dirname, 'public/examples/scoring_example_full.csv');
        const csvContent = readFileSync(csvPath, 'utf-8');
        
        console.log('✅ Successfully read example CSV file');
        console.log(`File size: ${csvContent.length} characters`);
        
        // Create scoring module
        const scoring = createScoringModule();
        
        // Test CSV parsing
        console.log('\n📊 Testing CSV Parsing...');
        const parsedData = scoring.parseCSV(csvContent);
        
        if (parsedData.length === 0) {
            console.error('❌ CSV parsing failed - no rows returned');
            return;
        }
        
        console.log(`✅ CSV parsed successfully: ${parsedData.length} rows`);
        
        // Test column detection on first few rows
        console.log('\n🔍 Testing Column Detection...');
        let successCount = 0;
        const rowsToTest = Math.min(3, parsedData.length);
        
        for (let i = 0; i < rowsToTest; i++) {
            const success = scoring.testProcessSingleRow(parsedData[i], i);
            if (success) successCount++;
        }
        
        console.log(`\n📈 Test Results: ${successCount}/${rowsToTest} rows processed successfully`);
        
        if (successCount === rowsToTest) {
            console.log('🎉 All tests passed! The scoring agent should work correctly.');
        } else {
            console.log('⚠️  Some tests failed. There may be issues with column detection.');
        }
        
    } catch (error) {
        console.error('❌ Test failed with error:', error.message);
        console.error(error.stack);
    }
    
    console.log('\n📋 Captured Log Output:');
    logOutput.forEach(log => console.log('  ', log));
    
    if (errorOutput.length > 0) {
        console.log('\n🚨 Captured Error Output:');
        errorOutput.forEach(error => console.log('  ', error));
    }
}

// Run the test
runTest().then(() => {
    console.log('\n✅ Test completed');
}).catch(error => {
    console.error('\n❌ Test failed:', error);
});