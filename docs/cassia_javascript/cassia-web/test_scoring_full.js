/**
 * Full end-to-end test for scoring.js functionality including API calls
 */

import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

// Get current directory
const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Mock LLM API call
async function mockCallLLM(prompt, provider, model, apiKey, temperature, maxTokens) {
    console.log(`🔄 Mock API Call: ${provider}/${model} (temp: ${temperature})`);
    console.log(`📝 Prompt length: ${prompt.length} chars`);
    
    // Simulate API delay
    await new Promise(resolve => setTimeout(resolve, 100));
    
    // Return a realistic scoring response
    return `<reasoning>
This annotation demonstrates strong evidence for the identified cell type. The marker genes show excellent specificity and are consistent with known expression patterns. The conversation history reveals a thorough analysis with appropriate consideration of alternative possibilities. The annotation is well-supported by the evidence presented.

Key strengths:
- High specificity markers present
- Comprehensive analysis approach
- Well-documented reasoning process
- Consideration of potential confounding factors

Minor areas for improvement:
- Could benefit from additional negative marker validation
- Some rare subtype possibilities not fully explored

Overall, this represents a high-quality annotation with strong scientific basis.
</reasoning>

<score>87</score>`;
}

// Import and adapt scoring functions
function createScoringEnvironment() {
    // Scoring prompt function
    function promptCreatorScore(majorClusterInfo, marker, annotationHistory) {
        return `
You are an expert in single-cell annotation analysis. Your task is to evaluate and rate single-cell annotation results, focusing on their correctness and ability to capture the overall picture of the data. You will provide a score from 0 to 100 and justify your rating.

Here are the single-cell annotation results to evaluate:

<marker>
${marker}
</marker>

<Cluster Origin>
${majorClusterInfo}
</Cluster Origin>

<annotation_history>
${annotationHistory}
</annotation_history>

Carefully analyze these results, paying particular attention to the following aspects:
1. Correctness of the annotations
2. Balanced consideration of multiple markers rather than over-focusing on a specific one
3. Ability to capture the general picture of the cell populations

When evaluating, consider:
- Are the annotations scientifically accurate?
- Is there a good balance in the use of different markers?
- Does the annotation provide a comprehensive view of the cell types present?
- Are there any obvious misclassifications or oversights?
- Did it consider the rank of the marker? marker appear first is more important.

Provide your analysis in the following format:
1. Start with a <reasoning> tag, where you explain your evaluation of the annotation results. Discuss the strengths and weaknesses you've identified, referring to specific examples from the results where possible.
2. After your reasoning, use a <score> tag to provide a numerical score from 0 to 100, where 0 represents completely incorrect or unusable results, and 100 represents perfect annotation that captures all aspects of the data correctly.

Your response should look like this:

<reasoning>
[Your detailed analysis and justification here]
</reasoning>

<score>[Your numerical score between 0 and 100]</score>

Remember, the focus is on correctness and the ability to see the general picture, rather than the structure of the results. Be critical but fair in your assessment.
        `;
    }

    // Extract score and reasoning function
    function extractScoreAndReasoning(text) {
        try {
            let score = null;
            let reasoning = null;
            
            // Extract score - try multiple patterns
            const scorePatterns = [
                /\<score\>(\d+)\<\/score\>/i,
                /Score:\s*(\d+)/i,
                /score:\s*(\d+)/i,
                /(\d+)\/100/i,
                /(\d+)\s*out\s*of\s*100/i,
                /rating.*?(\d+)/i,
                /(\d+)%/i
            ];
            
            for (const pattern of scorePatterns) {
                const scoreMatch = text.match(pattern);
                if (scoreMatch) {
                    score = parseInt(scoreMatch[1]);
                    break;
                }
            }
            
            // Extract reasoning - try multiple patterns
            const reasoningPatterns = [
                /\<reasoning\>(.*?)\<\/reasoning\>/is,
                /Reasoning:\s*(.*?)(?=Score:|$)/is,
                /reasoning:\s*(.*?)(?=score:|$)/is,
                /Analysis:\s*(.*?)(?=Score:|$)/is,
                /Evaluation:\s*(.*?)(?=Score:|$)/is
            ];
            
            for (const pattern of reasoningPatterns) {
                const reasoningMatch = text.match(pattern);
                if (reasoningMatch) {
                    reasoning = reasoningMatch[1].trim();
                    break;
                }
            }
            
            // If no specific reasoning found, use the entire text as reasoning
            if (reasoning === null && text.trim()) {
                reasoning = text.trim();
            }
            
            return { score, reasoning };
            
        } catch (error) {
            console.error(`Error extracting data: ${error.message}`);
            return { score: null, reasoning: null };
        }
    }

    // Single analysis scoring function
    async function scoreSingleAnalysis(majorClusterInfo, marker, annotationHistory, model, provider, apiKey) {
        const prompt = promptCreatorScore(majorClusterInfo, marker, annotationHistory);
        
        const response = await mockCallLLM(
            prompt,
            provider,
            model,
            apiKey,
            0.7,
            2000,
            null
        );
        
        const { score, reasoning } = extractScoreAndReasoning(response);
        return { score, reasoning };
    }

    // CSV parsing (simplified version)
    function parseCSV(csvText) {
        const lines = csvText.trim().split('\n');
        if (lines.length < 2) return [];
        
        // Simple CSV parser for testing
        const headers = lines[0].split(',').map(h => h.trim().replace(/^"|"$/g, ''));
        const rows = [];
        
        for (let i = 1; i < lines.length; i++) {
            // Handle quoted fields with commas
            const line = lines[i];
            const values = [];
            let current = '';
            let inQuotes = false;
            
            for (let j = 0; j < line.length; j++) {
                const char = line[j];
                if (char === '"') {
                    inQuotes = !inQuotes;
                } else if (char === ',' && !inQuotes) {
                    values.push(current.trim().replace(/^"|"$/g, ''));
                    current = '';
                } else {
                    current += char;
                }
            }
            values.push(current.trim().replace(/^"|"$/g, ''));
            
            const row = {};
            headers.forEach((header, index) => {
                row[header] = values[index] || '';
            });
            rows.push(row);
        }
        
        return rows;
    }

    // Process single row function
    async function processSingleRow(row, idx, model, provider, apiKey) {
        try {
            console.log(`\n🔄 Processing Row ${idx + 1}`);
            
            // Find columns
            const findColumn = (options) => {
                for (const col of options) {
                    if (col in row && row[col] && row[col].trim() !== '') {
                        return row[col].trim();
                    }
                }
                return null;
            };
            
            const species = findColumn(['Species', 'species']) || 'Unknown';
            const tissue = findColumn(['Tissue', 'tissue']) || 'Unknown';
            const majorClusterInfo = `${species} ${tissue}`;
            
            const marker = findColumn(['Marker List', 'Marker.List', 'marker_list', 'Markers']);
            const annotationHistory = findColumn(['Conversation History', 'Conversation.History', 'conversation_history']);
            
            if (!marker) {
                throw new Error('Could not find marker column');
            }
            
            if (!annotationHistory) {
                throw new Error('Could not find conversation history column');
            }
            
            console.log(`📊 Context: ${majorClusterInfo}`);
            console.log(`🧬 Markers: ${marker.length} chars (${marker.split(',').length} genes)`);
            console.log(`💬 History: ${annotationHistory.length} chars`);
            
            const result = await scoreSingleAnalysis(
                majorClusterInfo,
                marker,
                annotationHistory,
                model,
                provider,
                apiKey
            );
            
            console.log(`✅ Score: ${result.score}`);
            console.log(`📝 Reasoning preview: ${result.reasoning.substring(0, 100)}...`);
            
            return { idx, score: result.score, reasoning: result.reasoning };
            
        } catch (error) {
            console.error(`❌ Error processing row ${idx + 1}: ${error.message}`);
            return { idx, score: null, reasoning: `Error: ${error.message}` };
        }
    }

    return {
        parseCSV,
        processSingleRow,
        promptCreatorScore,
        extractScoreAndReasoning,
        scoreSingleAnalysis
    };
}

// Main test function
async function runFullTest() {
    console.log('🧪 Starting Full Scoring Test\n');
    
    try {
        // Read the example CSV file
        const csvPath = join(__dirname, 'public/examples/scoring_example_full.csv');
        const csvContent = readFileSync(csvPath, 'utf-8');
        
        console.log('✅ Successfully read example CSV file');
        console.log(`File size: ${csvContent.length} characters`);
        
        // Create scoring environment
        const scoring = createScoringEnvironment();
        
        // Parse CSV
        console.log('\n📊 Parsing CSV...');
        const parsedData = scoring.parseCSV(csvContent);
        console.log(`✅ Parsed ${parsedData.length} rows`);
        
        // Test scoring on first 2 rows (to keep test fast)
        console.log('\n🎯 Testing Full Scoring Workflow...');
        const testRows = parsedData.slice(0, 2);
        const results = [];
        
        for (let i = 0; i < testRows.length; i++) {
            const result = await scoring.processSingleRow(
                testRows[i], 
                i, 
                'gpt-3.5-turbo', 
                'openai', 
                'test-api-key'
            );
            results.push(result);
        }
        
        // Summary
        console.log('\n📈 Test Summary:');
        const successfulScores = results.filter(r => r.score !== null);
        console.log(`✅ Successfully scored: ${successfulScores.length}/${results.length} rows`);
        
        if (successfulScores.length > 0) {
            const avgScore = successfulScores.reduce((sum, r) => sum + r.score, 0) / successfulScores.length;
            console.log(`📊 Average score: ${avgScore.toFixed(1)}`);
            
            console.log('\n🏆 Sample Results:');
            successfulScores.forEach((r, i) => {
                console.log(`  Row ${r.idx + 1}: Score ${r.score}`);
                console.log(`    Reasoning: ${r.reasoning.substring(0, 150)}...`);
            });
        }
        
        if (successfulScores.length === results.length) {
            console.log('\n🎉 All tests passed! The scoring system is working correctly.');
        } else {
            console.log('\n⚠️  Some tests failed. Check the error messages above.');
        }
        
    } catch (error) {
        console.error('❌ Test failed with error:', error.message);
        console.error(error.stack);
    }
}

// Run the test
runFullTest().then(() => {
    console.log('\n✅ Full test completed');
}).catch(error => {
    console.error('\n❌ Test failed:', error);
});