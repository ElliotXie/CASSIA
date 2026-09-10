/**
 * Test DeepSeek API connection
 * Run with: node test/test-deepseek.js
 */

const OpenAI = require('openai').default;

const API_KEY = process.env.DEEPSEEK_API_KEY || 'YOUR_API_KEY_HERE';

// Test different base URL formats
const baseUrls = [
    'https://api.deepseek.com',
    'https://api.deepseek.com/v1',
];

async function testConnection(baseUrl) {
    console.log(`\n🧪 Testing: ${baseUrl}`);

    try {
        const client = new OpenAI({
            apiKey: API_KEY,
            baseURL: baseUrl,
        });

        const response = await client.chat.completions.create({
            model: 'deepseek-v4-flash',
            messages: [
                { role: 'user', content: 'Say hello in one word.' }
            ],
            thinking: { type: 'enabled' },
            reasoning_effort: 'low',
            max_tokens: 128,
        });

        console.log(`✅ SUCCESS with ${baseUrl}`);
        console.log(`   Response: ${response.choices[0].message.content}`);
        console.log(`   Thinking: ${response.choices[0].message.reasoning_content ? 'enabled' : 'missing'}`);
        return true;
    } catch (error) {
        console.log(`❌ FAILED with ${baseUrl}`);
        console.log(`   Error: ${error.message}`);
        if (error.status) {
            console.log(`   Status: ${error.status}`);
        }
        return false;
    }
}

async function main() {
    console.log('🔍 DeepSeek API Connection Test\n');

    for (const url of baseUrls) {
        const success = await testConnection(url);
        if (success) {
            console.log(`\n✅ Working base URL: ${url}`);
            break;
        }
    }
}

main();
