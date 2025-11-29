/**
 * Generate HTML reports for CASSIA subclustering results
 * Simplified JavaScript implementation of Python generate_reports.py
 */

import fs from 'fs';
import path from 'path';

/**
 * Generate a beautiful HTML report for subclustering results
 * @param {string} csvPath - Path to the subclustering results CSV
 * @param {string} htmlReportPath - Path for output HTML (optional)
 * @param {string} modelName - Model name (optional)
 */
export function generateSubclusteringReport(csvPath, htmlReportPath = null, modelName = null) {
    // Read CSV content
    const csvContent = fs.readFileSync(csvPath, 'utf-8');
    const lines = csvContent.split('\n').filter(line => line.trim());
    
    if (lines.length < 2) {
        console.log(`Warning: CSV file ${csvPath} appears to be empty`);
        return;
    }
    
    // Parse CSV
    const headers = lines[0].split(',').map(h => h.trim().replace(/"/g, ''));
    const rows = [];
    
    for (let i = 1; i < lines.length; i++) {
        const values = lines[i].match(/("([^"]|"")*"|[^,]+)/g) || [];
        const row = {};
        headers.forEach((header, idx) => {
            let value = values[idx] || '';
            // Remove surrounding quotes and unescape internal quotes
            value = value.trim();
            if (value.startsWith('"') && value.endsWith('"')) {
                value = value.slice(1, -1).replace(/""/g, '"');
            }
            row[header] = value;
        });
        if (Object.keys(row).length > 0) {
            rows.push(row);
        }
    }
    
    // Set output path and model name
    if (!htmlReportPath) {
        htmlReportPath = csvPath.replace('.csv', '_summary.html');
    }
    if (!modelName) {
        modelName = path.basename(csvPath).replace('.csv', '');
    }
    
    // Build HTML table rows
    const tableRows = rows.map((row, idx) => {
        const clusterId = row['Result ID'] || row['True Cell Type'] || '';
        const mainType = row['main_cell_type'] || '';
        const subType = row['sub_cell_type'] || '';
        const keyMarkers = row['key_markers'] || '';
        const reason = row['reason'] || row['explanation'] || '';
        const markerId = `marker_${idx + 1}`;
        
        return `
        <tr>
            <td class="cluster-col">${clusterId}</td>
            <td class="annotation-col">
                <div class="main-type">${mainType}</div>
                <div class="sub-type">${subType}</div>
            </td>
            <td class="marker-col">
                ${keyMarkers ? `
                <a href="javascript:void(0)" onclick="toggleMarkers('${markerId}')" class="toggle-link">Show markers</a>
                <div id="${markerId}" class="collapsible">${keyMarkers}</div>
                ` : '<span class="no-data">No markers</span>'}
            </td>
            <td class="reason-col">${reason || '<span class="no-data">No explanation provided</span>'}</td>
        </tr>`;
    }).join('');
    
    // Generate HTML
    const html = `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Subclustering Report - ${modelName}</title>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px; 
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        h1 { 
            color: #2c3e50; 
            text-align: center;
            margin-bottom: 10px;
        }
        .meta {
            text-align: center;
            color: #666;
            margin-bottom: 30px;
            font-size: 14px;
        }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin-top: 20px;
            background-color: white;
        }
        th { 
            background-color: #3498db; 
            color: white; 
            padding: 12px; 
            text-align: left;
            font-weight: bold;
        }
        td { 
            padding: 10px; 
            border-bottom: 1px solid #ddd;
            vertical-align: top;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        .cluster-col { 
            width: 10%; 
            font-weight: bold;
            color: #2c3e50;
        }
        .annotation-col { 
            width: 25%; 
        }
        .main-type {
            font-weight: bold;
            color: #27ae60;
            margin-bottom: 3px;
        }
        .sub-type {
            color: #7f8c8d;
            font-size: 0.9em;
            font-style: italic;
        }
        .marker-col { 
            width: 25%; 
        }
        .reason-col { 
            width: 40%; 
            font-size: 0.95em;
            line-height: 1.4;
        }
        .toggle-link {
            color: #3498db;
            text-decoration: none;
            font-size: 0.9em;
            cursor: pointer;
        }
        .toggle-link:hover {
            text-decoration: underline;
        }
        .collapsible {
            display: none;
            margin-top: 5px;
            padding: 8px;
            background-color: #ecf0f1;
            border-radius: 4px;
            font-size: 0.9em;
            line-height: 1.6;
        }
        .no-data {
            color: #bdc3c7;
            font-style: italic;
        }
    </style>
    <script>
    function toggleMarkers(id) {
        var el = document.getElementById(id);
        var link = event.target;
        if (el.style.display === 'none' || el.style.display === '') {
            el.style.display = 'block';
            link.textContent = 'Hide markers';
        } else {
            el.style.display = 'none';
            link.textContent = 'Show markers';
        }
    }
    </script>
</head>
<body>
    <div class="container">
        <h1>Subclustering Annotation Report</h1>
        <div class="meta">
            <b>Model:</b> ${modelName} &nbsp; | &nbsp; 
            <b>Generated:</b> ${new Date().toLocaleString()} &nbsp; | &nbsp;
            <b>Total Subclusters:</b> ${rows.length}
        </div>
        <table>
            <tr>
                <th>Cluster</th>
                <th>Annotation<br><span style="font-weight:normal;font-size:12px">(Main / Subtype)</span></th>
                <th>Top Markers</th>
                <th>Reasoning</th>
            </tr>
            ${tableRows}
        </table>
    </div>
</body>
</html>`;
    
    // Write HTML file
    fs.writeFileSync(htmlReportPath, html);
    console.log(`Subclustering HTML report saved to ${htmlReportPath}`);
}

/**
 * Process evaluation CSV and generate appropriate HTML report
 * @param {string} csvPath - Path to CSV file
 * @param {boolean} overwrite - Whether to overwrite existing HTML
 */
export async function processEvaluationCsv(csvPath, overwrite = false) {
    try {
        if (!fs.existsSync(csvPath)) {
            console.log(`File not found: ${csvPath}`);
            return;
        }
        
        const modelName = path.basename(csvPath).replace('.csv', '');
        const htmlPath = csvPath.replace('.csv', '_summary.html');
        
        if (fs.existsSync(htmlPath) && !overwrite) {
            console.log(`HTML report already exists for ${modelName}. Skipping.`);
            return;
        }
        
        // Read CSV to check format
        const content = fs.readFileSync(csvPath, 'utf-8');
        const firstLine = content.split('\n')[0];
        
        // Check if it's a subclustering format
        // Accept both 'Result ID' and 'True Cell Type' as valid cluster column names
        const hasMainCellType = firstLine.includes('main_cell_type');
        const hasSubCellType = firstLine.includes('sub_cell_type');
        const hasClusterCol = firstLine.includes('Result ID') || firstLine.includes('True Cell Type');
        const hasSubclusteringCols = hasMainCellType && hasSubCellType && hasClusterCol;
        
        if (hasSubclusteringCols) {
            generateSubclusteringReport(csvPath, htmlPath, modelName);
        } else {
            console.log(`Warning: CSV format not recognized as subclustering results`);
            // Could implement other report types here in the future
        }
        
    } catch (error) {
        console.log(`Error processing ${csvPath}: ${error.message}`);
    }
}

/**
 * Create an index.html file that links to all the reports
 * @param {Array<string>} csvFiles - List of CSV file paths
 * @param {string} outputDir - Output directory for index.html
 */
export async function createIndexHtml(csvFiles, outputDir) {
    const reports = [];
    
    // Process each CSV file
    for (const csvFile of csvFiles) {
        try {
            if (!fs.existsSync(csvFile)) continue;
            
            const basename = path.basename(csvFile);
            const htmlFile = csvFile.replace('.csv', '_summary.html');
            const htmlBasename = path.basename(htmlFile);
            
            // Count rows in CSV
            const content = fs.readFileSync(csvFile, 'utf-8');
            const lines = content.split('\n').filter(line => line.trim());
            const rowCount = Math.max(0, lines.length - 1); // Subtract header
            
            reports.push({
                name: basename.replace('.csv', ''),
                csvFile: basename,
                htmlFile: htmlBasename,
                rowCount: rowCount,
                hasHtml: fs.existsSync(htmlFile)
            });
            
        } catch (error) {
            console.log(`Error processing ${csvFile} for index: ${error.message}`);
        }
    }
    
    // Sort reports by name
    reports.sort((a, b) => a.name.localeCompare(b.name));
    
    // Generate index HTML
    const indexHtml = `<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Subclustering Batch Results</title>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px; 
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1000px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        h1 { 
            color: #2c3e50; 
            text-align: center;
        }
        .summary {
            text-align: center;
            color: #666;
            margin-bottom: 30px;
        }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin-top: 20px;
        }
        th { 
            background-color: #3498db; 
            color: white; 
            padding: 12px; 
            text-align: left;
        }
        td { 
            padding: 10px; 
            border-bottom: 1px solid #ddd;
        }
        tr:hover {
            background-color: #f8f9fa;
        }
        a {
            color: #3498db;
            text-decoration: none;
        }
        a:hover {
            text-decoration: underline;
        }
        .no-report {
            color: #bdc3c7;
            font-style: italic;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Subclustering Batch Results</h1>
        <div class="summary">
            <b>Total Runs:</b> ${reports.length} &nbsp; | &nbsp;
            <b>Generated:</b> ${new Date().toLocaleString()}
        </div>
        <table>
            <tr>
                <th>Analysis Name</th>
                <th>Subclusters</th>
                <th>CSV File</th>
                <th>HTML Report</th>
            </tr>
            ${reports.map(report => `
            <tr>
                <td><strong>${report.name}</strong></td>
                <td>${report.rowCount}</td>
                <td><a href="${report.csvFile}">${report.csvFile}</a></td>
                <td>${report.hasHtml ? 
                    `<a href="${report.htmlFile}">View Report</a>` : 
                    '<span class="no-report">Not generated</span>'
                }</td>
            </tr>`).join('')}
        </table>
    </div>
</body>
</html>`;
    
    // Write index.html
    const indexPath = path.join(outputDir, 'index.html');
    fs.writeFileSync(indexPath, indexHtml);
    console.log(`Batch index HTML created at ${indexPath}`);
}

// Export all functions
export default {
    generateSubclusteringReport,
    processEvaluationCsv,
    createIndexHtml
};