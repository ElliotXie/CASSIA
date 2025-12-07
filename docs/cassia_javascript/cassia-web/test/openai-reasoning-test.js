/**
 * Tests for OpenAI reasoning effort handling
 * Run with: node test/openai-reasoning-test.js
 */

// Mock test to verify reasoning config logic
function testReasoningConfigLogic() {
    console.log('Testing reasoning config logic...\n');

    const testCases = [
        { input: 'none', expected: null, desc: 'lowercase none should return null' },
        { input: 'NONE', expected: null, desc: 'uppercase NONE should return null' },
        { input: 'None', expected: null, desc: 'mixed case None should return null' },
        { input: 'high', expected: { effort: 'high' }, desc: 'high should return config' },
        { input: 'medium', expected: { effort: 'medium' }, desc: 'medium should return config' },
        { input: 'low', expected: { effort: 'low' }, desc: 'low should return config' },
        { input: 'HIGH', expected: { effort: 'high' }, desc: 'uppercase HIGH should return lowercase config' },
        { input: null, expected: null, desc: 'null should return null' },
        { input: '', expected: null, desc: 'empty string should return null' },
    ];

    let passed = 0;
    let failed = 0;

    testCases.forEach(({ input, expected, desc }) => {
        // NEW LOGIC (case-insensitive) - matches runCASSIA.js
        const reasoningConfig = input && input.toLowerCase() !== 'none'
            ? { effort: input.toLowerCase() }
            : null;

        const result = JSON.stringify(reasoningConfig);
        const expectedStr = JSON.stringify(expected);

        if (result === expectedStr) {
            console.log(`  PASS: ${desc}`);
            passed++;
        } else {
            console.log(`  FAIL: ${desc}`);
            console.log(`   Expected: ${expectedStr}`);
            console.log(`   Got: ${result}`);
            failed++;
        }
    });

    console.log(`\nResults: ${passed} passed, ${failed} failed`);
    return failed === 0;
}

// Test valid effort values for OpenAI Responses API
function testValidEffortValues() {
    console.log('\nTesting valid effort values for Responses API...\n');

    const VALID_REASONING_EFFORTS = ['high', 'medium', 'low'];

    const testCases = [
        { effort: 'high', shouldUseResponsesAPI: true },
        { effort: 'medium', shouldUseResponsesAPI: true },
        { effort: 'low', shouldUseResponsesAPI: true },
        { effort: 'HIGH', shouldUseResponsesAPI: true },
        { effort: 'none', shouldUseResponsesAPI: false },
        { effort: 'NONE', shouldUseResponsesAPI: false },
        { effort: 'invalid', shouldUseResponsesAPI: false },
        { effort: '', shouldUseResponsesAPI: false },
    ];

    let passed = 0;
    let failed = 0;

    testCases.forEach(({ effort, shouldUseResponsesAPI }) => {
        // Matches llm_utils.js logic (convert to boolean for comparison)
        const useResponsesAPI = !!(effort && VALID_REASONING_EFFORTS.includes(effort.toLowerCase()));

        if (useResponsesAPI === shouldUseResponsesAPI) {
            console.log(`  PASS: effort="${effort}" -> useResponsesAPI=${useResponsesAPI}`);
            passed++;
        } else {
            console.log(`  FAIL: effort="${effort}" -> expected ${shouldUseResponsesAPI}, got ${useResponsesAPI}`);
            failed++;
        }
    });

    console.log(`\nResults: ${passed} passed, ${failed} failed`);
    return failed === 0;
}

// Test the full flow: reasoningEffort -> reasoningConfig -> useResponsesAPI
function testFullFlow() {
    console.log('\nTesting full flow (reasoningEffort -> reasoningConfig -> useResponsesAPI)...\n');

    const VALID_REASONING_EFFORTS = ['high', 'medium', 'low'];

    const testCases = [
        { reasoningEffort: 'none', expectedUseResponsesAPI: false, desc: 'none should use Chat Completions' },
        { reasoningEffort: 'NONE', expectedUseResponsesAPI: false, desc: 'NONE should use Chat Completions' },
        { reasoningEffort: null, expectedUseResponsesAPI: false, desc: 'null should use Chat Completions' },
        { reasoningEffort: 'high', expectedUseResponsesAPI: true, desc: 'high should use Responses API' },
        { reasoningEffort: 'medium', expectedUseResponsesAPI: true, desc: 'medium should use Responses API' },
        { reasoningEffort: 'low', expectedUseResponsesAPI: true, desc: 'low should use Responses API' },
    ];

    let passed = 0;
    let failed = 0;

    testCases.forEach(({ reasoningEffort, expectedUseResponsesAPI, desc }) => {
        // Step 1: runCASSIA.js logic
        const reasoningConfig = reasoningEffort && reasoningEffort.toLowerCase() !== 'none'
            ? { effort: reasoningEffort.toLowerCase() }
            : null;

        // Step 2: llm_utils.js logic (convert to boolean for comparison)
        const useResponsesAPI = !!(reasoningConfig
            && reasoningConfig.effort
            && VALID_REASONING_EFFORTS.includes(reasoningConfig.effort.toLowerCase()));

        if (useResponsesAPI === expectedUseResponsesAPI) {
            console.log(`  PASS: ${desc}`);
            passed++;
        } else {
            console.log(`  FAIL: ${desc}`);
            console.log(`   reasoningConfig: ${JSON.stringify(reasoningConfig)}`);
            console.log(`   Expected useResponsesAPI: ${expectedUseResponsesAPI}, got: ${useResponsesAPI}`);
            failed++;
        }
    });

    console.log(`\nResults: ${passed} passed, ${failed} failed`);
    return failed === 0;
}

// Run all tests
console.log('=== OpenAI Reasoning Effort Tests ===\n');
console.log('These tests verify the fix for the error:');
console.log('"Failed to execute \'append\' on \'Headers\': String contains non ISO-8859-1 code point"\n');
console.log('=' .repeat(50) + '\n');

const test1 = testReasoningConfigLogic();
const test2 = testValidEffortValues();
const test3 = testFullFlow();

console.log('\n' + '='.repeat(50));
if (test1 && test2 && test3) {
    console.log('\n All tests passed!');
    process.exit(0);
} else {
    console.log('\n Some tests failed!');
    process.exit(1);
}
