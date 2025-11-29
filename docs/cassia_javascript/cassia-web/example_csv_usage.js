/**
 * Example usage of the standalone CSV parser module
 * Demonstrates how other modules can easily import and use the CSV parser
 */

import { parseCSV, formatAsCSV, validateCSV, getCSVStats } from './lib/utils/csv-parser.js';

// Example CSV data with complex content
const sampleCSV = `Name,Age,Notes,Markers
"John Doe",25,"Has multiple
lines and ""quotes""","Gene1, Gene2, Gene3"
"Jane Smith",30,"Simple note","Gene4, Gene5"
"Bob Johnson",35,"Contains, commas and
newlines in text","Gene6"`;

console.log('📋 CSV Parser Usage Examples\n');

// 1. Basic parsing
console.log('1️⃣ Basic CSV Parsing:');
const parsedData = parseCSV(sampleCSV, { debug: true });
console.log('Parsed data:', parsedData);

// 2. Parsing without debug output
console.log('\n2️⃣ Quiet Parsing (debug: false):');
const quietData = parseCSV(sampleCSV, { debug: false });
console.log(`Quietly parsed ${quietData.length} rows`);

// 3. Validation
console.log('\n3️⃣ CSV Validation:');
const validationResult = validateCSV(parsedData, {
    requiredColumns: ['Name', 'Age', 'Markers'],
    debug: true
});
console.log('Validation result:', validationResult);

// 4. Statistics
console.log('\n4️⃣ CSV Statistics:');
const stats = getCSVStats(parsedData);
console.log('CSV Stats:', stats);

// 5. Converting back to CSV
console.log('\n5️⃣ Converting back to CSV:');
const csvOutput = formatAsCSV(parsedData);
console.log('CSV Output:');
console.log(csvOutput);

// 6. Usage in other modules
console.log('\n6️⃣ How to use in other modules:');
console.log(`
// Import the parser functions
import { parseCSV, formatAsCSV, validateCSV } from './lib/utils/csv-parser.js';

// Parse CSV from file or string
const data = parseCSV(csvString, { debug: false });

// Validate required columns exist
const validation = validateCSV(data, { 
    requiredColumns: ['required_column1', 'required_column2'] 
});

if (!validation.valid) {
    console.error('CSV validation failed:', validation.errors);
    return;
}

// Process your data
data.forEach(row => {
    // Your processing logic here
    console.log(row.column_name);
});

// Export back to CSV if needed
const outputCSV = formatAsCSV(processedData);
`);

console.log('\n✅ CSV Parser module is ready for use across your codebase!');
console.log('Key features:');
console.log('- ✅ Handles multiline quoted fields');
console.log('- ✅ Properly processes escaped quotes');
console.log('- ✅ Supports mixed line endings');
console.log('- ✅ Includes validation and statistics');
console.log('- ✅ Can export back to CSV format');
console.log('- ✅ Debug mode for troubleshooting');