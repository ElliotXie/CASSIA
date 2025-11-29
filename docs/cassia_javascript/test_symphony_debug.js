#!/usr/bin/env node

/**
 * Test script to debug Symphony Compare issues
 */

import { callLLM } from './cassia-web/lib/cassia/llm_utils.js';
import { symphonyCompare, extractCelltypeScores } from './cassia-web/lib/cassia/symphonyCompare.js';

async function testBasicAPI() {
    console.log('🧪 Testing basic API call...');
    
    const API_KEY = process.env.OPENROUTER_API_KEY;
    if (!API_KEY) {
        console.error('❌ OPENROUTER_API_KEY not set in environment');
        return false;
    }
    
    try {
        const response = await callLLM(
            'Say "Hello World" and nothing else.',
            'openrouter',
            'google/gemini-2.5-flash',
            API_KEY,
            0,
            10
        );
        
        console.log('✅ Basic API call successful');
        console.log('📝 Response:', response);
        return true;
    } catch (error) {
        console.error('❌ Basic API call failed:', error.message);
        return false;
    }
}

async function testScoreExtraction() {
    console.log('\n🧪 Testing score extraction...');
    
    const mockResponse = `
<celltype>T cells</celltype>
<reasoning>
T cells express CD3, CD4, and CD8 markers which are present in this marker set.
</reasoning>
<score>85</score>

<celltype>B cells</celltype>
<reasoning>
B cells typically express CD19 and CD20, but I don't see these markers clearly.
</reasoning>
<score>35</score>
`;

    const celltypes = ['T cells', 'B cells'];
    const extracted = extractCelltypeScores(mockResponse, celltypes);
    
    console.log('📊 Extracted scores:', extracted);
    
    // Validate extraction worked
    if (extracted['T cells'].score === '85' && extracted['B cells'].score === '35') {
        console.log('✅ Score extraction working correctly');
        return true;
    } else {
        console.log('❌ Score extraction failed');
        return false;
    }
}

async function testSymphonyCompare() {
    console.log('\n🧪 Testing Symphony Compare...');
    
    const API_KEY = process.env.OPENROUTER_API_KEY;
    if (!API_KEY) {
        console.error('❌ OPENROUTER_API_KEY not set in environment');
        return false;
    }
    
    try {
        const result = await symphonyCompare({
            tissue: 'blood',
            celltypes: ['T cells', 'B cells'],
            markerSet: 'CD3D, CD3E, CD19, CD20',
            species: 'human',
            modelPreset: 'budget',  // Use budget to test fewer models
            enableDiscussion: false,  // Disable discussion for simpler test
            apiKey: API_KEY,
            provider: 'openrouter'
        });
        
        console.log('✅ Symphony Compare completed');
        console.log('📊 Results summary:', {
            consensus: result.consensus,
            confidence: result.confidence,
            resultsCount: result.results?.length || 0
        });
        
        if (result.results && result.results.length > 0) {
            console.log('📝 First result sample:', {
                model: result.results[0].model,
                status: result.results[0].status,
                extracted_scores: result.results[0].extracted_scores
            });
        }
        
        return true;
    } catch (error) {
        console.error('❌ Symphony Compare failed:', error.message);
        console.error('🔍 Error details:', error);
        return false;
    }
}

async function runTests() {
    console.log('🎼 Symphony Compare Debug Tests\n');
    
    const basicAPI = await testBasicAPI();
    if (!basicAPI) {
        console.log('\n❌ Basic API test failed - stopping here');
        return;
    }
    
    const scoreExtraction = await testScoreExtraction();
    if (!scoreExtraction) {
        console.log('\n❌ Score extraction test failed - stopping here');
        return;
    }
    
    const symphonyTest = await testSymphonyCompare();
    if (!symphonyTest) {
        console.log('\n❌ Symphony Compare test failed');
        return;
    }
    
    console.log('\n✅ All tests passed! Symphony Compare should work in frontend.');
}

// Run tests
runTests().catch(console.error);