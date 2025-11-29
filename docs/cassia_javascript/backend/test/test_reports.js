import { runCASSIAGenerateScoreReport, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';
import { promises as fs } from 'fs';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Set API key for testing
const API_KEY = "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";
setApiKey(API_KEY, "openrouter");

async function testReportGeneration() {
    console.log("🧪 CASSIA Report Generation Test");
    console.log("=================================");
    
    try {
        // Test 1: Check for existing scored results
        console.log("\n📋 Test 1: Looking for scored results to generate reports from...");
        
        const testResultsDir = path.join(__dirname, 'test_results');
        let scoredFile = null;
        
        try {
            const files = await fs.readdir(testResultsDir);
            const scoredFiles = files.filter(f => f.includes('scored') && f.endsWith('.csv'));
            if (scoredFiles.length > 0) {
                scoredFile = path.join(testResultsDir, scoredFiles[0]);
                console.log(`Found scored results: ${scoredFiles[0]}`);
            } else {
                console.log("No scored results found. Looking for any _full.csv files...");
                const fullFiles = files.filter(f => f.includes('_full.csv'));
                if (fullFiles.length > 0) {
                    scoredFile = path.join(testResultsDir, fullFiles[0]);
                    console.log(`Using batch results: ${fullFiles[0]}`);
                    console.log("⚠️  Note: This file may not have scores, but we'll try to generate reports anyway");
                }
            }
        } catch (error) {
            console.log("No test_results directory found");
        }
        
        if (!scoredFile) {
            console.log("❌ No suitable CSV files found for report generation");
            console.log("Please run the scoring test first: npm run test-scoring");
            return;
        }
        
        // Test 2: Generate reports from scored data
        console.log("\n📋 Test 2: Generating HTML reports...");
        
        try {
            const results = await runCASSIAGenerateScoreReport(
                scoredFile,
                "test_reports_index"
            );
            
            console.log("✅ Report generation completed:");
            console.log(`   - Index file: ${results.indexFile}`);
            console.log(`   - Total reports: ${results.totalReports}`);
            console.log(`   - Report files:`);
            results.reportFiles.forEach((file, index) => {
                console.log(`     ${index + 1}. ${path.basename(file)}`);
            });
            
            // Test 3: Verify files were created
            console.log("\n📋 Test 3: Verifying generated files...");
            
            // Check index file
            try {
                await fs.access(results.indexFile);
                const indexStats = await fs.stat(results.indexFile);
                console.log(`✅ Index file created: ${Math.round(indexStats.size / 1024)}KB`);
            } catch (error) {
                console.error(`❌ Index file not found: ${error.message}`);
            }
            
            // Check individual report files
            let validReports = 0;
            for (const reportFile of results.reportFiles) {
                try {
                    await fs.access(reportFile);
                    validReports++;
                } catch (error) {
                    console.error(`❌ Report file not found: ${path.basename(reportFile)}`);
                }
            }
            
            console.log(`✅ ${validReports}/${results.reportFiles.length} report files verified`);
            
            // Test 4: Check file contents
            console.log("\n📋 Test 4: Checking report content...");
            
            try {
                const indexContent = await fs.readFile(results.indexFile, 'utf8');
                const hasValidHTML = indexContent.includes('<!DOCTYPE html>') && 
                                  indexContent.includes('CASSIA Analysis Reports') &&
                                  indexContent.includes('report-card');
                console.log(`✅ Index HTML structure: ${hasValidHTML ? 'Valid' : 'Invalid'}`);
                
                if (results.reportFiles.length > 0) {
                    const firstReportContent = await fs.readFile(results.reportFiles[0], 'utf8');
                    const hasReportStructure = firstReportContent.includes('<!DOCTYPE html>') && 
                                             firstReportContent.includes('CASSIA Analysis Report') &&
                                             firstReportContent.includes('agent-section');
                    console.log(`✅ Report HTML structure: ${hasReportStructure ? 'Valid' : 'Invalid'}`);
                }
                
            } catch (error) {
                console.error(`❌ Error checking file contents: ${error.message}`);
            }
            
            console.log("\n🎉 Report generation test completed successfully!");
            console.log(`📁 Open ${results.indexFile} in your browser to view the reports`);
            
        } catch (error) {
            console.error("❌ Report generation failed:", error.message);
        }
        
    } catch (error) {
        console.error("❌ Report generation test failed:", error.message);
    }
}

// Run the report generation test
testReportGeneration();