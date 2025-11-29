import { readFileSync, writeFileSync } from 'fs';
import { join } from 'path';
import { scoreAnnotationBatch } from './lib/cassia/scoring.js';
import { formatAsCSV } from './lib/utils/csv-parser.js';

// This is a special script to run the scoring batch process on a specific file
// and save the results to a new CSV file.

async function main() {
    console.log('🚀 Starting scoring process...');

    const inputFile = 'my_annotationnew1_full.csv';
    const outputFile = 'scored_output.csv';

    try {
        // Read the input CSV file
        const csvContent = readFileSync(join(process.cwd(), 'cassia-web', inputFile), 'utf-8');
        console.log(`✅ Successfully read ${inputFile}`);

        // Define a progress callback to see what's happening
        const onProgress = ({ completed, total, percentage }) => {
            process.stdout.write(`🔄 Progress: ${completed}/${total} rows processed (${percentage}%)\r`);
        };
        
        const onLog = (message) => {
            console.log(message);
        };

        // Run the scoring batch process
        // NOTE: This uses the actual `scoreAnnotationBatch` function.
        // It will make mock LLM calls as defined in the test setup for scoring.js if run in a test env,
        // but here we expect it to run fully. We need to provide a dummy API key.
        const { results, csvContent: scoredCsvContent } = await scoreAnnotationBatch({
            csvData: csvContent,
            apiKey: process.env.OPENROUTER_API_KEY || 'dummy-key-for-testing', // Use env var or a dummy key
            maxWorkers: 4,
            onProgress,
            onLog
        });

        console.log('\n✅ Scoring complete.');

        // Save the results to the output CSV file
        writeFileSync(join(process.cwd(), 'cassia-web', outputFile), scoredCsvContent, 'utf-8');
        console.log(`💾 Results saved to ${outputFile}`);

        console.log('🎉 Script finished successfully.');

    } catch (error) {
        console.error('❌ An error occurred during the scoring process:');
        console.error(error);
        process.exit(1);
    }
}

main(); 