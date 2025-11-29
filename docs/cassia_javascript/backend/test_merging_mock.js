import fs from 'fs';
import path from 'path';
import { parse as parseCSV } from 'csv-parse/sync';
import { stringify as stringifyCSV } from 'csv-stringify/sync';

// Mock the LLM utils to test parsing logic
const mockLLMResponse = {
    broad: `Here are the suggested broader cell groupings:

{
  "1": "T cells",
  "2": "Myeloid cells", 
  "3": "B cells"
}`,
    detailed: `{
  "1": "CD4 T cells",
  "2": "Macrophages",
  "3": "B cells"
}`,
    very_detailed: `{
  "1": "Naive CD4+ T cells",
  "2": "Inflammatory macrophages", 
  "3": "Memory B cells"
}`
};

// Import and modify the response parsing function
function _parseLLMResponse(response, batch) {
    const groupings = new Array(batch.length).fill("Error parsing response");
    
    try {
        const jsonMatch = response.match(/({[\s\S]*})/);
        if (jsonMatch) {
            const jsonStr = jsonMatch[1];
            const parsed = JSON.parse(jsonStr);
            
            let i = 0;
            for (const [clusterId, grouping] of Object.entries(parsed)) {
                if (i < batch.length) {
                    groupings[i] = grouping;
                    i++;
                }
            }
        }
    } catch (error) {
        console.error(`Error parsing LLM response: ${error.message}`);
    }
    
    return groupings;
}

async function testMockMerging() {
    console.log('🧪 Testing Merging with Mock LLM Responses...\n');
    
    // Create test data
    const testData = [
        { 'True Cell Type': '1', 'Predicted Main Cell Type': 'T cell', 'Predicted Sub Cell Types': 'naive CD4 T cell', 'processed_subtype': 'naive CD4 T cell' },
        { 'True Cell Type': '2', 'Predicted Main Cell Type': 'macrophage', 'Predicted Sub Cell Types': 'inflammatory macrophage', 'processed_subtype': 'inflammatory macrophage' },
        { 'True Cell Type': '3', 'Predicted Main Cell Type': 'B cell', 'Predicted Sub Cell Types': 'memory B cell', 'processed_subtype': 'memory B cell' }
    ];
    
    // Test each detail level
    const detailLevels = ['broad', 'detailed', 'very_detailed'];
    
    for (const detailLevel of detailLevels) {
        console.log(`--- Testing ${detailLevel} detail level ---`);
        
        const response = mockLLMResponse[detailLevel];
        const groupings = _parseLLMResponse(response, testData);
        
        console.log('Mock LLM Response:');
        console.log(response);
        console.log('\nParsed Groupings:');
        groupings.forEach((grouping, idx) => {
            console.log(`Cluster ${testData[idx]['True Cell Type']}: ${testData[idx]['Predicted Main Cell Type']} → ${grouping}`);
        });
        console.log('');
    }
    
    // Test combining all results
    console.log('--- Testing Combined Results ---');
    const combinedData = testData.map(row => ({ ...row }));
    
    // Add all three grouping columns
    const resultColumnMap = {
        "broad": "Merged_Grouping_1",
        "detailed": "Merged_Grouping_2", 
        "very_detailed": "Merged_Grouping_3"
    };
    
    for (const detailLevel of detailLevels) {
        const response = mockLLMResponse[detailLevel];
        const groupings = _parseLLMResponse(response, testData);
        const resultColumn = resultColumnMap[detailLevel];
        
        groupings.forEach((grouping, idx) => {
            combinedData[idx][resultColumn] = grouping;
        });
    }
    
    console.log('Combined results for all detail levels:');
    combinedData.forEach(row => {
        console.log(`Cluster ${row['True Cell Type']}: ${row['Predicted Main Cell Type']}`);
        console.log(`  → Broad: ${row['Merged_Grouping_1']}`);
        console.log(`  → Detailed: ${row['Merged_Grouping_2']}`);
        console.log(`  → Very Detailed: ${row['Merged_Grouping_3']}`);
        console.log('');
    });
    
    // Save to CSV to verify format
    const csvOutput = stringifyCSV(combinedData, { header: true });
    const outputPath = '/mnt/c/Users/ellio/OneDrive - UW-Madison/Revision_cassia/cassia_javascript/backend/test_results/mock_merge_result.csv';
    fs.writeFileSync(outputPath, csvOutput);
    console.log(`✅ Mock test results saved to: ${outputPath}`);
    
    console.log('\n🎉 All mock merging tests completed successfully!');
}

testMockMerging().catch(console.error);