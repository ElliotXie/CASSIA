import { promises as fs, createReadStream } from 'fs';
import path from 'path';
import csvParser from 'csv-parser';

// ----------------- HTML Template Functions -----------------

/**
 * Generate complete HTML report from analysis text (100% Python compatible)
 * @param {string} analysisText - Combined analysis text with sections
 * @returns {string} Complete HTML report
 */
function generateHTMLReport(analysisText) {
    // Split the text into sections based on agents
    const sections = analysisText.split(" | ");
    
    // HTML template with CSS styling - exact Python replica
    const htmlTemplate = `
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body { 
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif; 
                max-width: 1200px; 
                margin: 0 auto; 
                padding: 20px; 
                background-color: #f0f2f5;
                line-height: 1.6;
            }
            .container { 
                background-color: white; 
                padding: 40px; 
                border-radius: 16px; 
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }
            .agent-section { 
                margin-bottom: 35px; 
                padding: 25px; 
                border-radius: 12px; 
                transition: all 0.3s ease;
            }
            .agent-section:hover {
                transform: translateY(-2px);
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }
            .final-annotation { 
                background-color: #f0f7ff; 
                border-left: 5px solid #2196f3; 
            }
            .validator { 
                background-color: #f0fdf4; 
                border-left: 5px solid #22c55e; 
            }
            .formatting { 
                background: linear-gradient(145deg, #fff7ed, #ffe4c4);
                border-left: 5px solid #f97316; 
                box-shadow: 0 4px 15px rgba(249, 115, 22, 0.1);
            }
            h2 { 
                color: #1a2b3c; 
                margin-top: 0; 
                font-size: 1.5rem;
                font-weight: 600;
                display: flex;
                align-items: center;
                gap: 10px;
            }
            ul { 
                margin: 15px 0; 
                padding-left: 20px; 
            }
            pre { 
                background-color: #f8fafc; 
                padding: 20px; 
                border-radius: 8px; 
                overflow-x: auto;
                font-family: 'Consolas', 'Monaco', monospace;
                font-size: 0.9rem;
                line-height: 1.5;
            }
            .validation-result { 
                font-weight: 600; 
                color: #16a34a; 
                padding: 12px 20px;
                background-color: #dcfce7; 
                border-radius: 8px; 
                display: inline-block;
                margin: 10px 0;
            }
            br { 
                margin-bottom: 8px; 
            }
            p {
                margin: 12px 0;
                color: #374151;
            }
            .summary-content {
                display: flex;
                flex-direction: column;
                gap: 24px;
            }
            .summary-item {
                display: flex;
                flex-direction: column;
                gap: 8px;
                background: rgba(255, 255, 255, 0.7);
                padding: 16px;
                border-radius: 12px;
                backdrop-filter: blur(8px);
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
            }
            .summary-label {
                font-weight: 600;
                color: #c2410c;
                font-size: 0.95rem;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }
            .summary-value {
                color: #1f2937;
                font-size: 1.1rem;
                padding: 8px 16px;
                background-color: rgba(255, 255, 255, 0.9);
                border-radius: 8px;
                display: inline-block;
                box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            }
            .summary-list {
                margin: 0;
                padding-left: 24px;
                list-style-type: none;
            }
            .summary-list li {
                color: #1f2937;
                padding: 8px 0;
                position: relative;
            }
            .summary-list li:before {
                content: "•";
                color: #f97316;
                font-weight: bold;
                position: absolute;
                left: -20px;
            }
            .report-header {
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }
            
            .report-title {
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }
            
            .report-subtitle {
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }
            .scoring { 
                background: linear-gradient(145deg, #f0fdf4, #dcfce7);
                border-left: 5px solid #22c55e;
                box-shadow: 0 4px 15px rgba(34, 197, 94, 0.1);
            }
            .scoring-content {
                display: flex;
                flex-direction: column;
                gap: 16px;
                color: #1f2937;
                line-height: 1.8;
            }
            .scoring-content br + br {
                content: "";
                display: block;
                margin: 12px 0;
            }
            .empty-list {
                color: #6b7280;
                font-style: italic;
            }
            .error-message {
                color: #dc2626;
                padding: 12px;
                background-color: #fef2f2;
                border-radius: 6px;
                border-left: 4px solid #dc2626;
            }
            .score-badge {
                background: linear-gradient(135deg, #22c55e, #16a34a);
                color: white;
                padding: 8px 16px;
                border-radius: 12px;
                font-size: 1.5rem;
                font-weight: 700;
                display: inline-block;
                margin: 12px 0;
                box-shadow: 0 4px 12px rgba(34, 197, 94, 0.2);
                position: relative;
                top: -10px;
            }
            .score-badge::before {
                content: "Score:";
                font-size: 0.9rem;
                font-weight: 500;
                margin-right: 8px;
                opacity: 0.9;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <h1 class="report-title">CASSIA Analysis Report</h1>
                <p class="report-subtitle">Comprehensive Cell Type Analysis and Annotation</p>
            </div>
            {content}
        </div>
    </body>
    </html>
    `;
    
    const content = [];
    
    // Process each section
    for (const section of sections) {
        if (section.startsWith("Final Annotation Agent:")) {
            const annotationContent = section.replace("Final Annotation Agent:", "").trim();
            content.push(`
                <div class="agent-section final-annotation">
                    <h2>🔍 Final Annotation Analysis</h2>
                    ${annotationContent.replace(/\n/g, '<br>')}
                </div>
            `);
            
        } else if (section.startsWith("Coupling Validator:")) {
            const validatorContent = section.replace("Coupling Validator:", "").trim();
            const validationResult = validatorContent.includes("VALIDATION PASSED") 
                ? '<div class="validation-result">✅ VALIDATION PASSED</div>' 
                : "";
            
            content.push(`
                <div class="agent-section validator">
                    <h2>✓ Validation Check</h2>
                    ${validationResult}
                    ${validatorContent.replace(/\n/g, '<br>')}
                </div>
            `);
            
        } else if (section.startsWith("Formatting Agent:")) {
            try {
                // Get the content after "Formatting Agent:"
                const jsonText = section.replace("Formatting Agent:", "").trim();
                
                // Since the JSON is consistently formatted with newlines,
                // we can find where it ends (the last '}' followed by a newline or end of string)
                const jsonEnd = jsonText.lastIndexOf('}');
                if (jsonEnd !== -1) {
                    const jsonContent = jsonText.substring(0, jsonEnd + 1);
                    const data = JSON.parse(jsonContent);
                    
                    // Process the data...
                    const mainCellType = data.main_cell_type || 'Not specified';
                    const subCellTypes = data.sub_cell_types || [];
                    const mixedTypes = data.possible_mixed_cell_types || [];
                    const numMarkers = data.num_markers || 'Not specified';
                    
                    // Format the content...
                    const formattedContent = `
                        <div class="summary-content">
                            <div class="summary-item">
                                <span class="summary-label">Main Cell Type:</span>
                                <span class="summary-value">${mainCellType}</span>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Sub Cell Types:</span>
                                <ul class="summary-list">
                                    ${subCellTypes.length > 0 
                                        ? subCellTypes.map(item => `<li>${item}</li>`).join('')
                                        : '<li class="empty-list">No sub cell types identified</li>'}
                                </ul>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Possible Mixed Cell Types:</span>
                                <ul class="summary-list">
                                    ${mixedTypes.length > 0 
                                        ? mixedTypes.map(item => `<li>${item}</li>`).join('')
                                        : '<li class="empty-list">No mixed cell types identified</li>'}
                                </ul>
                            </div>
                            
                            <div class="summary-item">
                                <span class="summary-label">Number of Markers:</span>
                                <span class="summary-value">${numMarkers}</span>
                            </div>
                        </div>
                    `;
                    
                    content.push(`
                        <div class="agent-section formatting">
                            <h2>📋 Summary</h2>
                            ${formattedContent}
                        </div>
                    `);
                } else {
                    throw new Error("Could not find JSON content");
                }
                    
            } catch (error) {
                content.push(`
                    <div class="agent-section formatting">
                        <h2>📋 Summary</h2>
                        <p class="error-message">Error formatting data: ${error.message}</p>
                    </div>
                `);
            }
        } else if (section.startsWith("Scoring Agent:")) {
            try {
                // Get the content after "Scoring Agent:"
                const scoringText = section.split("Scoring Agent:")[1].trim();
                
                // Split the score from the main text
                const lastScoreIndex = scoringText.lastIndexOf("Score:");
                const mainText = scoringText.substring(0, lastScoreIndex).trim();
                const score = scoringText.substring(lastScoreIndex + 6).trim();
                
                content.push(`
                    <div class="agent-section scoring">
                        <h2>🎯 Quality Assessment</h2>
                        <div class="score-badge">${score}</div>
                        <div class="scoring-content">
                            ${mainText.replace(/\n/g, '<br>')}
                        </div>
                    </div>
                `);
            } catch (error) {
                content.push(`
                    <div class="agent-section scoring">
                        <h2>🎯 Quality Assessment</h2>
                        <p class="error-message">Error formatting scoring data: ${error.message}</p>
                    </div>
                `);
            }
        }
    }
    
    // Combine all sections
    const finalHTML = htmlTemplate.replace('{content}', content.join(''));
    return finalHTML;
}

/**
 * Process single report data into HTML (100% Python compatible)
 * @param {string} text - Conversation history text
 * @param {string} scoreReasoning - Scoring reasoning text
 * @param {number} score - Numerical score
 * @returns {string} Complete HTML report
 */
function processSingleReport(text, scoreReasoning, score) {
    const combined = `${text}\n | Scoring Agent: ${scoreReasoning}\nScore: ${score}`;
    return generateHTMLReport(combined);
}

/**
 * Generate index page HTML (100% Python compatible)
 * @param {Array<string>} reportFiles - Array of report filenames
 * @returns {string} Complete index page HTML
 */
function generateIndexPage(reportFiles) {
    const indexTemplate = `
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body { 
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif; 
                max-width: 1200px; 
                margin: 0 auto; 
                padding: 20px; 
                background-color: #f0f2f5;
                line-height: 1.6;
            }
            .container { 
                background-color: white; 
                padding: 40px; 
                border-radius: 16px; 
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }
            .report-list {
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
                gap: 20px;
                margin-top: 30px;
            }
            .report-card {
                background: linear-gradient(145deg, #fff, #f8fafc);
                border-radius: 12px;
                padding: 24px;
                text-decoration: none;
                color: #1a2b3c;
                transition: all 0.3s ease;
                border: 1px solid #e2e8f0;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }
            .report-card:hover {
                transform: translateY(-4px);
                box-shadow: 0 8px 25px rgba(0,0,0,0.15);
                background: linear-gradient(145deg, #fff, #f1f5f9);
                text-decoration: none;
                color: #1a2b3c;
            }
            .report-title {
                font-size: 1.2rem;
                font-weight: 600;
                margin-bottom: 8px;
                color: #f97316;
            }
            .report-description {
                color: #64748b;
                font-size: 0.95rem;
                line-height: 1.5;
            }
            .index-header {
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }
            .index-title {
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }
            .index-subtitle {
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }
            .report-count {
                background: linear-gradient(135deg, #f97316, #c2410c);
                color: white;
                padding: 8px 16px;
                border-radius: 20px;
                font-size: 0.9rem;
                font-weight: 600;
                display: inline-block;
                margin-top: 16px;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div class="index-header">
                <h1 class="index-title">CASSIA Analysis Reports</h1>
                <p class="index-subtitle">Comprehensive Cell Type Analysis and Annotation Results</p>
                <div class="report-count">${reportFiles.length} Reports Generated</div>
            </div>
            
            <div class="report-list">
                ${reportFiles.map(file => {
                    const displayName = file.replace('report_', '').replace('.html', '');
                    const formattedName = displayName.charAt(0).toUpperCase() + displayName.slice(1).replace(/_/g, ' ');
                    return `
                        <a href="${file}" class="report-card">
                            <div class="report-title">${formattedName}</div>
                            <div class="report-description">
                                Detailed analysis report for ${formattedName.toLowerCase()} cell type annotation
                            </div>
                        </a>
                    `;
                }).join('')}
            </div>
        </div>
    </body>
    </html>
    `;
    
    return indexTemplate;
}

// ----------------- CSV Helper Functions -----------------

/**
 * Read CSV file and return as array of objects
 */
async function readCSV(filePath) {
    return new Promise((resolve, reject) => {
        const results = [];
        createReadStream(filePath)
            .pipe(csvParser())
            .on('data', (data) => results.push(data))
            .on('end', () => resolve(results))
            .on('error', reject);
    });
}

/**
 * Sanitize filename for safe file system usage
 * @param {string} filename - Original filename
 * @returns {string} Sanitized filename
 */
function sanitizeFilename(filename) {
    return filename.toString().trim()
        .replace(/[^a-zA-Z0-9\s\-_]/g, '')  // Only alphanumeric, spaces, hyphens, underscores
        .trim();
}

// ----------------- Main Report Generation Function -----------------

/**
 * Generate HTML reports from a scored CSV file and create an index page (100% Python compatible)
 * @param {string} csvPath - Path to the CSV file containing the score results
 * @param {string} indexName - Base name for the index file (without .html extension)
 * @returns {Promise<void>}
 */
export async function runCASSIAGenerateScoreReport(csvPath, indexName = "CASSIA_reports_summary") {
    console.log(`Generating HTML reports from: ${csvPath}`);
    
    // Read the CSV file
    const report = await readCSV(csvPath);
    const reportFiles = [];
    
    // Determine output folder (same folder as the CSV file)
    let outputFolder = path.dirname(csvPath);
    if (!outputFolder || outputFolder === '.') {
        outputFolder = '.';
    }
    
    console.log(`Output folder: ${outputFolder}`);
    
    // Process each row
    for (let index = 0; index < report.length; index++) {
        const row = report[index];
        
        // Get the first column value for the filename
        const firstColumnKey = Object.keys(row)[0];
        const filename = sanitizeFilename(row[firstColumnKey]);
        
        // Handle multiple possible column name variations for conversation history
        const historyColumnOptions = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History'];
        let text = null;
        for (const col of historyColumnOptions) {
            if (col in row) {
                text = row[col];
                break;
            }
        }
        if (text === null) {
            throw new Error(`Could not find conversation history column. Available columns: ${Object.keys(row).join(', ')}`);
        }
        
        // Handle multiple possible column name variations for scoring reasoning
        const reasoningColumnOptions = ['Scoring_Reasoning', 'Scoring.Reasoning', 'scoring_reasoning', 'Scoring_reasoning'];
        let scoreReasoning = null;
        for (const col of reasoningColumnOptions) {
            if (col in row) {
                scoreReasoning = row[col];
                break;
            }
        }
        if (scoreReasoning === null) {
            throw new Error(`Could not find scoring reasoning column. Available columns: ${Object.keys(row).join(', ')}`);
        }
        
        const score = row["Score"];
        
        // Generate HTML for this row
        const htmlContent = processSingleReport(text, scoreReasoning, score);
        
        // Save using the first column value as filename in the output folder
        const outputPath = path.join(outputFolder, `report_${filename}.html`);
        await fs.writeFile(outputPath, htmlContent, 'utf8');
        
        // Store just the filename for the index (not the full path)
        reportFiles.push(path.basename(outputPath));
        console.log(`Report saved to ${outputPath}`);
    }
    
    // Generate and save index page in the same folder
    const indexHTML = generateIndexPage(reportFiles);
    const indexFilename = path.join(outputFolder, `${path.basename(indexName)}.html`);
    await fs.writeFile(indexFilename, indexHTML, 'utf8');
    console.log(`Index page saved to ${indexFilename}`);
    
    console.log(`Report generation complete! Generated ${reportFiles.length} reports with index page.`);
    
    return {
        indexFile: indexFilename,
        reportFiles: reportFiles.map(file => path.join(outputFolder, file)),
        totalReports: reportFiles.length
    };
}

// ----------------- Exports -----------------

export default {
    runCASSIAGenerateScoreReport,
    generateHTMLReport,
    processSingleReport,
    generateIndexPage
};