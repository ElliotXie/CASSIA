/**
 * Industrial-Strength CSV Parser
 * 
 * A robust CSV parser that handles all edge cases including:
 * - Multiline quoted fields (conversation history with newlines)
 * - Escaped quotes within fields
 * - Commas inside quoted content
 * - Mixed line endings (\r\n, \n, \r)
 * - Empty fields and rows
 * - Special characters and unicode
 * 
 * @author CASSIA Team
 * @version 2.0.0
 */

/**
 * Parse CSV text into array of objects
 * @param {string} csvText - CSV content as string
 * @param {Object} options - Parsing options
 * @param {boolean} options.header - Whether first row contains headers (default: true)
 * @param {boolean} options.skipEmptyLines - Skip empty lines (default: true)
 * @param {boolean} options.debug - Enable debug logging (default: false)
 * @returns {Array} Array of objects representing CSV rows
 */
export function parseCSV(csvText, options = {}) {
    const {
        header = true,
        skipEmptyLines = true,
        debug = false
    } = options;

    if (!csvText || typeof csvText !== 'string') {
        if (debug) console.error('Invalid CSV input:', typeof csvText);
        return [];
    }
    
    const text = csvText.trim();
    if (text.length === 0) return [];
    
    if (debug) console.log(`🔍 Starting robust CSV parsing for ${text.length} character input`);
    
    try {
        const result = parseCSVRobust(text, debug);
        if (debug) console.log(`✅ CSV parsing successful: ${result.length} rows`);
        
        if (!header) {
            return result;
        }
        
        // Convert to objects with headers
        if (result.length === 0) return [];
        
        const headers = result[0].map(h => h.trim());
        const dataRows = result.slice(1);
        
        return dataRows.map((row, rowIndex) => {
            const obj = {};
            headers.forEach((headerName, colIndex) => {
                let value = '';
                if (colIndex < row.length) {
                    value = row[colIndex] || '';
                    value = value.trim();
                }
                obj[headerName] = value;
            });
            
            // Debug first few rows
            if (debug && rowIndex < 2) {
                const debugInfo = {
                    rowNumber: rowIndex + 1,
                    columnCount: row.length,
                    expectedColumns: headers.length
                };
                
                // Check specific important columns
                if ('Marker List' in obj) {
                    debugInfo['Marker List'] = obj['Marker List'] ? `${obj['Marker List'].length} chars` : 'EMPTY';
                }
                if ('Conversation History' in obj) {
                    debugInfo['Conversation History'] = obj['Conversation History'] ? `${obj['Conversation History'].length} chars` : 'EMPTY';
                }
                
                console.log(`🔍 Row ${rowIndex + 1} debug:`, debugInfo);
            }
            
            return obj;
        });
        
    } catch (error) {
        if (debug) console.error('❌ Robust CSV parsing failed:', error.message);
        if (debug) console.log('🔄 Attempting fallback parsing...');
        
        try {
            const result = parseCSVFallback(text, debug);
            if (debug) console.log(`⚠️ Fallback parsing succeeded: ${result.length} rows`);
            return result;
        } catch (fallbackError) {
            if (debug) console.error('❌ All CSV parsing methods failed:', fallbackError.message);
            return [];
        }
    }
}

/**
 * Primary robust CSV parser - handles all edge cases
 * @param {string} text - CSV text to parse
 * @param {boolean} debug - Enable debug logging
 * @returns {Array} Array of arrays representing CSV rows
 */
function parseCSVRobust(text, debug = false) {
    const rows = [];
    let pos = 0;
    let currentRow = [];
    let currentField = '';
    let inQuotes = false;
    let quoteCount = 0;
    
    const length = text.length;
    let lineNumber = 1;
    let fieldNumber = 1;
    
    while (pos < length) {
        const char = text[pos];
        const nextChar = pos + 1 < length ? text[pos + 1] : null;
        
        try {
            if (char === '"') {
                quoteCount++;
                
                if (!inQuotes) {
                    // Starting quoted field
                    if (currentField.length === 0) {
                        inQuotes = true;
                    } else {
                        // Quote in middle of unquoted field - treat as literal
                        currentField += char;
                    }
                } else {
                    // Inside quoted field
                    if (nextChar === '"') {
                        // Escaped quote - add literal quote and skip next
                        currentField += '"';
                        pos += 2;
                        continue;
                    } else {
                        // End quote
                        inQuotes = false;
                    }
                }
            } else if (char === ',' && !inQuotes) {
                // Field separator
                currentRow.push(currentField);
                currentField = '';
                fieldNumber++;
            } else if ((char === '\n' || (char === '\r' && nextChar !== '\n') || (char === '\r' && nextChar === '\n')) && !inQuotes) {
                // Row separator (handle \n, \r, or \r\n)
                currentRow.push(currentField);
                
                if (currentRow.length > 0 || currentField.length > 0) {
                    rows.push([...currentRow]);
                }
                
                currentRow = [];
                currentField = '';
                fieldNumber = 1;
                lineNumber++;
                
                // Skip \n if we just processed \r\n
                if (char === '\r' && nextChar === '\n') {
                    pos++;
                }
            } else {
                // Regular character
                currentField += char;
            }
            
            pos++;
            
        } catch (charError) {
            if (debug) console.error(`Error processing character at position ${pos}, line ${lineNumber}, field ${fieldNumber}:`, charError);
            // Skip problematic character and continue
            pos++;
        }
    }
    
    // Handle last field/row
    if (currentField.length > 0 || currentRow.length > 0) {
        currentRow.push(currentField);
        rows.push(currentRow);
    }
    
    // Validate and clean up
    if (rows.length === 0) {
        throw new Error('No rows found in CSV');
    }
    
    if (debug) console.log(`📊 Parsed ${rows.length} rows, ${quoteCount} quotes processed`);
    
    return rows;
}

/**
 * Fallback CSV parser for when robust parsing fails
 * @param {string} text - CSV text to parse
 * @param {boolean} debug - Enable debug logging
 * @returns {Array} Array of objects representing CSV rows
 */
function parseCSVFallback(text, debug = false) {
    if (debug) console.log('🔄 Using fallback CSV parser');
    
    // Split by lines first
    const lines = text.split(/\r?\n/);
    if (lines.length < 2) {
        throw new Error('Not enough lines in CSV');
    }
    
    // Simple split for headers (assuming no commas in header names)
    const headers = lines[0].split(',').map(h => h.trim().replace(/^"|"$/g, ''));
    
    const rows = [];
    for (let i = 1; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.length === 0) continue;
        
        // Very simple parsing - just split by comma and remove outer quotes
        const values = line.split(',').map(v => {
            v = v.trim();
            // Remove outer quotes if present
            if (v.startsWith('"') && v.endsWith('"')) {
                v = v.slice(1, -1);
            }
            // Handle escaped quotes
            v = v.replace(/""/g, '"');
            return v;
        });
        
        const obj = {};
        headers.forEach((header, index) => {
            obj[header] = values[index] || '';
        });
        
        rows.push(obj);
    }
    
    if (debug) console.log(`📊 Fallback parser: ${rows.length} rows`);
    return rows;
}

/**
 * Convert array of objects to CSV format
 * @param {Array} data - Array of objects
 * @param {Object} options - Formatting options
 * @param {boolean} options.includeHeaders - Include header row (default: true)
 * @returns {string} CSV formatted string
 */
export function formatAsCSV(data, options = {}) {
    const { includeHeaders = true } = options;
    
    if (!data || data.length === 0) return '';
    
    const headers = Object.keys(data[0]);
    
    // Helper function to properly escape CSV cells (matching runCASSIA_batch approach)
    const escapeCsvCell = (cell) => {
        if (cell === null || cell === undefined) {
            return ''; // Return empty string for null/undefined, don't quote it
        }
        
        const cellStr = String(cell);
        
        // Only wrap in quotes if it contains special characters or is already quoted
        const needsQuotes = cellStr.includes(',') || cellStr.includes('"') || cellStr.includes('\n') || cellStr.includes('\r');
        
        if (needsQuotes) {
            // Escape any existing double quotes by doubling them up
            const escapedStr = cellStr.replace(/"/g, '""');
            return `"${escapedStr}"`;
        }
        
        // If no special characters, return the string as is
        return cellStr;
    };
    
    const csvRows = [];
    
    if (includeHeaders) {
        csvRows.push(headers.map(escapeCsvCell).join(','));
    }
    
    // Process data rows
    data.forEach(row => {
        const rowValues = headers.map(header => escapeCsvCell(row[header] || ''));
        csvRows.push(rowValues.join(','));
    });
    
    return csvRows.join('\n');
}

/**
 * Validate CSV structure and content
 * @param {Array} data - Parsed CSV data
 * @param {Object} options - Validation options
 * @param {Array} options.requiredColumns - Array of required column names
 * @param {boolean} options.debug - Enable debug logging
 * @returns {Object} Validation result with { valid: boolean, errors: Array, warnings: Array }
 */
export function validateCSV(data, options = {}) {
    const { requiredColumns = [], debug = false } = options;
    
    const result = {
        valid: true,
        errors: [],
        warnings: []
    };
    
    if (!Array.isArray(data)) {
        result.valid = false;
        result.errors.push('Data must be an array');
        return result;
    }
    
    if (data.length === 0) {
        result.valid = false;
        result.errors.push('No data rows found');
        return result;
    }
    
    // Check for required columns
    const firstRow = data[0];
    const availableColumns = Object.keys(firstRow);
    
    for (const requiredCol of requiredColumns) {
        if (!availableColumns.includes(requiredCol)) {
            result.valid = false;
            result.errors.push(`Required column missing: ${requiredCol}`);
        }
    }
    
    // Check for consistent column structure
    const expectedColumnCount = availableColumns.length;
    for (let i = 0; i < data.length; i++) {
        const rowColumnCount = Object.keys(data[i]).length;
        if (rowColumnCount !== expectedColumnCount) {
            result.warnings.push(`Row ${i + 1} has ${rowColumnCount} columns, expected ${expectedColumnCount}`);
        }
    }
    
    if (debug) {
        console.log('CSV Validation Result:', result);
        console.log(`Available columns: ${availableColumns.join(', ')}`);
    }
    
    return result;
}

/**
 * Get basic statistics about CSV data
 * @param {Array} data - Parsed CSV data
 * @returns {Object} Statistics object
 */
export function getCSVStats(data) {
    if (!Array.isArray(data) || data.length === 0) {
        return {
            rowCount: 0,
            columnCount: 0,
            columns: [],
            emptyFields: 0,
            totalFields: 0
        };
    }
    
    const columns = Object.keys(data[0]);
    let emptyFields = 0;
    let totalFields = 0;
    
    for (const row of data) {
        for (const col of columns) {
            totalFields++;
            if (!row[col] || row[col].trim() === '') {
                emptyFields++;
            }
        }
    }
    
    return {
        rowCount: data.length,
        columnCount: columns.length,
        columns,
        emptyFields,
        totalFields,
        completeness: ((totalFields - emptyFields) / totalFields * 100).toFixed(1) + '%'
    };
}

// Default export with all functions
export default {
    parseCSV,
    formatAsCSV,
    validateCSV,
    getCSVStats
};