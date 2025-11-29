"use client";

import { useState } from 'react';
import { callLLM } from '@/lib/cassia/llm_utils';
import { symphonyCompare, extractCelltypeScores } from '@/lib/cassia/symphonyCompare';

export default function TestSymphonyPage() {
    const [apiKey, setApiKey] = useState('');
    const [testResults, setTestResults] = useState<any[]>([]);
    const [isRunning, setIsRunning] = useState(false);

    const addResult = (test: string, status: 'pass' | 'fail', details: any) => {
        setTestResults(prev => [...prev, { test, status, details, timestamp: new Date().toISOString() }]);
    };

    const testBasicAPI = async () => {
        if (!apiKey) {
            addResult('Basic API', 'fail', 'No API key provided');
            return false;
        }

        try {
            const response = await callLLM(
                'Say "Hello World" and nothing else.',
                'openrouter',
                'google/gemini-2.5-flash',
                apiKey,
                0,
                10
            );
            
            addResult('Basic API', 'pass', { response: response.substring(0, 100) });
            return true;
        } catch (error: any) {
            addResult('Basic API', 'fail', { error: error.message });
            return false;
        }
    };

    const testScoreExtraction = async () => {
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
        
        if (extracted['T cells'].score === '85' && extracted['B cells'].score === '35') {
            addResult('Score Extraction', 'pass', extracted);
            return true;
        } else {
            addResult('Score Extraction', 'fail', extracted);
            return false;
        }
    };

    const testPromptGeneration = async () => {
        if (!apiKey) {
            addResult('Prompt Generation', 'fail', 'No API key provided');
            return false;
        }

        try {
            const response = await callLLM(
                `You are a professional biologist. Your task is to analyze how well a given marker set matches a list of cell types from human blood.

For EACH of the following cell types, you must provide your analysis in a specific structured format.
The cell types to analyze are:
- T cells
- B cells

The required output format for EACH cell type is:
<celltype>cell type name</celltype>
<reasoning>
Your detailed reasoning for the match, considering each marker's relevance.
</reasoning>
<score>A score from 0-100 indicating the match quality.</score>

Please provide a complete block of <celltype>, <reasoning>, and <score> for every cell type listed above.

Ranked marker set: CD3D, CD3E, CD19, CD20`,
                'openrouter',
                'google/gemini-2.5-flash',
                apiKey,
                0.5,
                4000
            );
            
            const extracted = extractCelltypeScores(response, ['T cells', 'B cells']);
            const hasValidScores = Object.values(extracted).some((result: any) => result.score !== "No score found");
            
            if (hasValidScores) {
                addResult('Prompt Generation', 'pass', { response: response.substring(0, 200), extracted });
                return true;
            } else {
                addResult('Prompt Generation', 'fail', { response: response.substring(0, 200), extracted });
                return false;
            }
        } catch (error: any) {
            addResult('Prompt Generation', 'fail', { error: error.message });
            return false;
        }
    };

    const testSymphonyCompare = async () => {
        if (!apiKey) {
            addResult('Symphony Compare', 'fail', 'No API key provided');
            return false;
        }

        try {
            const result = await symphonyCompare({
                tissue: 'blood',
                celltypes: ['T cells', 'B cells'],
                markerSet: 'CD3D, CD3E, CD19, CD20',
                species: 'human',
                modelPreset: 'budget',
                enableDiscussion: false,
                apiKey: apiKey,
                provider: 'openrouter'
            });
            
            const hasResults = result.results && result.results.length > 0;
            const hasValidScores = result.results?.some((r: any) => 
                r.extracted_scores && Object.keys(r.extracted_scores).length > 0
            );
            
            if (hasResults && hasValidScores) {
                addResult('Symphony Compare', 'pass', {
                    consensus: result.consensus,
                    confidence: result.confidence,
                    resultsCount: result.results?.length || 0,
                    sampleResult: result.results[0]
                });
                return true;
            } else {
                addResult('Symphony Compare', 'fail', {
                    resultsCount: result.results?.length || 0,
                    results: result.results
                });
                return false;
            }
        } catch (error: any) {
            addResult('Symphony Compare', 'fail', { error: error.message, stack: error.stack });
            return false;
        }
    };

    const runAllTests = async () => {
        setIsRunning(true);
        setTestResults([]);

        console.log('🧪 Running Symphony Compare tests...');

        // Test 1: Score extraction (no API needed)
        await testScoreExtraction();
        
        // Test 2: Basic API call
        const basicAPIWorks = await testBasicAPI();
        if (!basicAPIWorks) {
            setIsRunning(false);
            return;
        }

        // Test 3: Prompt generation and score extraction
        const promptWorks = await testPromptGeneration();
        if (!promptWorks) {
            setIsRunning(false);
            return;
        }

        // Test 4: Full symphony compare
        await testSymphonyCompare();

        setIsRunning(false);
    };

    return (
        <div className="min-h-screen bg-gradient-to-br from-purple-50 via-blue-50 to-indigo-50 p-8">
            <div className="max-w-4xl mx-auto">
                <h1 className="text-3xl font-bold mb-8">🧪 Symphony Compare Test Suite</h1>
                
                <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
                    <h2 className="text-xl font-semibold mb-4">Test Configuration</h2>
                    <div className="space-y-4">
                        <div>
                            <label className="block text-sm font-medium mb-2">OpenRouter API Key</label>
                            <input
                                type="password"
                                value={apiKey}
                                onChange={(e) => setApiKey(e.target.value)}
                                placeholder="Enter your OpenRouter API key"
                                className="w-full px-3 py-2 border border-gray-300 rounded-md"
                            />
                        </div>
                        <button
                            onClick={runAllTests}
                            disabled={isRunning || !apiKey}
                            className="bg-blue-600 text-white px-6 py-2 rounded-md hover:bg-blue-700 disabled:opacity-50"
                        >
                            {isRunning ? 'Running Tests...' : 'Run All Tests'}
                        </button>
                    </div>
                </div>

                <div className="bg-white rounded-lg shadow-lg p-6">
                    <h2 className="text-xl font-semibold mb-4">Test Results</h2>
                    {testResults.length === 0 ? (
                        <p className="text-gray-500">No tests run yet</p>
                    ) : (
                        <div className="space-y-4">
                            {testResults.map((result, index) => (
                                <div
                                    key={index}
                                    className={`p-4 rounded-md border-l-4 ${
                                        result.status === 'pass'
                                            ? 'bg-green-50 border-green-500'
                                            : 'bg-red-50 border-red-500'
                                    }`}
                                >
                                    <div className="flex items-center justify-between mb-2">
                                        <h3 className="font-medium">
                                            {result.status === 'pass' ? '✅' : '❌'} {result.test}
                                        </h3>
                                        <span className="text-sm text-gray-500">
                                            {new Date(result.timestamp).toLocaleTimeString()}
                                        </span>
                                    </div>
                                    <pre className="text-sm bg-gray-100 p-2 rounded overflow-auto max-h-40">
                                        {JSON.stringify(result.details, null, 2)}
                                    </pre>
                                </div>
                            ))}
                        </div>
                    )}
                </div>
            </div>
        </div>
    );
}