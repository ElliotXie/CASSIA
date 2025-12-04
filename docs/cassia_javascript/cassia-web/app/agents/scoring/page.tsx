"use client";

import { useState, ChangeEvent, useEffect } from 'react';
import Link from 'next/link';
import { scoreAnnotationBatch } from '@/lib/cassia/scoring';
import { useApiKeyStore } from '@/lib/stores/api-key-store';
import { parseCSV } from '@/lib/utils/csv-parser';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { ArrowLeft, Play, HelpCircle, Target, Upload, Download } from 'lucide-react';
import modelSettings from '../../../public/examples/model_settings.json';

export default function ScoringAgentPage() {
    const globalApiKey = useApiKeyStore((state) => state.getApiKey());
    const globalProvider = useApiKeyStore((state) => state.provider);
    
    // State management
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState('openrouter');
    const [model, setModel] = useState('google/gemini-2.5-flash');
    const [maxWorkers, setMaxWorkers] = useState(4);
    const [maxRetries, setMaxRetries] = useState(1);
    const [csvFile, setCsvFile] = useState<File | null>(null);
    const [csvData, setCsvData] = useState<any[]>([]);
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState<any>(null);
    const [error, setError] = useState<string | null>(null);
    const [progress, setProgress] = useState({ completed: 0, total: 0, percentage: 0 });

    // Initialize with global settings
    useEffect(() => {
        if (globalApiKey && !apiKey) {
            setApiKey(globalApiKey);
        }
        if (globalProvider && !provider) {
            setProvider(globalProvider);
        }
    }, [globalApiKey, globalProvider]);

    // Update model when provider changes
    useEffect(() => {
        const providerData = modelSettings.providers[provider as keyof typeof modelSettings.providers];
        if (providerData) {
            setModel(providerData.default_model);
        }
    }, [provider]);

    // Get available models for current provider
    const getAvailableModels = () => {
        const providerData = modelSettings.providers[provider as keyof typeof modelSettings.providers];
        if (!providerData) return [];
        return Object.entries(providerData.models || {}).map(([key, model]) => ({
            id: model.actual_name,
            name: model.description.split(' - ')[0] || model.actual_name
        }));
    };

    const handleFileUpload = (e: ChangeEvent<HTMLInputElement>) => {
        const file = e.target.files?.[0];
        if (file) {
            setCsvFile(file);
            const reader = new FileReader();
            reader.onload = (event) => {
                try {
                    const text = event.target?.result as string;
                    
                    // Use proper CSV parser instead of naive split
                    const data = parseCSV(text, { debug: true });
                    
                    setCsvData(data);
                    console.log('CSV data loaded:', data.length, 'rows');
                } catch (err) {
                    console.error('Error parsing CSV:', err);
                    setError('Failed to parse CSV file');
                }
            };
            reader.readAsText(file);
        }
    };

    const loadExampleData = async () => {
        try {
            setError(null);
            const response = await fetch('/examples/cassia_analysis_full.csv');
            const text = await response.text();
            
            // Use proper CSV parser instead of naive split
            const data = parseCSV(text, { debug: true });
            
            setCsvData(data);
            const mockFile = new File([text], 'cassia_analysis_full.csv', { type: 'text/csv' });
            setCsvFile(mockFile);
            console.log('Example data loaded successfully');
        } catch (error) {
            setError('Failed to load example data.');
        }
    };

    const handleRunScoring = async () => {
        if (!csvData.length || !apiKey) {
            setError('Please upload data and provide API key');
            return;
        }

        setIsLoading(true);
        setError(null);
        setResults(null);
        setProgress({ completed: 0, total: 0, percentage: 0 });

        try {
            const result = await scoreAnnotationBatch({
                csvData,
                apiKey,
                maxWorkers: Number(maxWorkers),
                model,
                provider,
                maxRetries: Number(maxRetries),
                onProgress: (progressData: any) => {
                    setProgress(progressData);
                }
            } as any);
            
            setResults(result);
            console.log('Scoring completed:', result);
        } catch (err: any) {
            setError(`Scoring failed: ${err.message}`);
        } finally {
            setIsLoading(false);
        }
    };

    return (
        <div className="min-h-screen">
            {/* Modern Header */}
            <header className="glass border-b border-white/20 sticky top-0 z-50">
                <div className="container mx-auto px-6 py-6">
                    <div className="flex items-center justify-between">
                        <div className="flex items-center space-x-4">
                            <Button variant="ghost" size="sm" asChild className="glass border-white/30 hover:bg-white/20">
                                <Link href="/">
                                    <ArrowLeft className="h-4 w-4 mr-2" />
                                    Back
                                </Link>
                            </Button>
                            <div className="flex items-center space-x-4">
                                <div className="relative">
                                    <div className="w-12 h-12 bg-gradient-to-br from-purple-600 to-pink-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                                        <Target className="h-6 w-6 text-white" />
                                    </div>
                                    <div className="absolute -top-1 -right-1 w-4 h-4 bg-green-400 rounded-full animate-pulse"></div>
                                </div>
                                <div>
                                    <h1 className="text-3xl font-bold gradient-text">Scoring Agent</h1>
                                    <p className="text-gray-600 dark:text-gray-300">Quality evaluation for annotations</p>
                                </div>
                            </div>
                        </div>
                        <Button variant="outline" size="sm" className="glass border-white/30 hover:bg-white/20 btn-modern">
                            <HelpCircle className="h-4 w-4 mr-2" />
                            Help
                        </Button>
                    </div>
                </div>
            </header>

            <main className="container mx-auto px-6 py-12">
                <div className="grid grid-cols-1 xl:grid-cols-5 gap-8">
                    {/* Left Panel - Configuration */}
                    <div className="xl:col-span-2">
                        <div className="glass rounded-2xl p-6 border border-white/20 sticky top-24 space-y-6">
                            <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                                ⚙️ <span className="ml-2">Configuration</span>
                            </h2>

                            {/* API Settings */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🔑 API Configuration</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">API Key</label>
                                        <Input
                                            type="password"
                                            value={apiKey}
                                            onChange={(e) => setApiKey(e.target.value)}
                                            placeholder="Enter your API key"
                                        />
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Provider</label>
                                        <select
                                            value={provider}
                                            onChange={(e) => setProvider(e.target.value)}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        >
                                            <option value="openrouter">OpenRouter</option>
                                            <option value="openai">OpenAI</option>
                                            <option value="anthropic">Anthropic</option>
                                        </select>
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Model</label>
                                        <select
                                            value={model}
                                            onChange={(e) => setModel(e.target.value)}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        >
                                            {getAvailableModels().map((m) => (
                                                <option key={m.id} value={m.id}>
                                                    {m.name}
                                                </option>
                                            ))}
                                        </select>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Analysis Settings */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🎯 Analysis Settings</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">
                                            Max Workers: {maxWorkers}
                                        </label>
                                        <input
                                            type="range"
                                            min="1"
                                            max="10"
                                            step="1"
                                            value={maxWorkers}
                                            onChange={(e) => setMaxWorkers(Number(e.target.value))}
                                            className="w-full"
                                        />
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">
                                            Max Retries: {maxRetries}
                                        </label>
                                        <input
                                            type="range"
                                            min="0"
                                            max="5"
                                            step="1"
                                            value={maxRetries}
                                            onChange={(e) => setMaxRetries(Number(e.target.value))}
                                            className="w-full"
                                        />
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Data Upload */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">📁 Data Upload</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Upload CSV File</label>
                                        <input
                                            type="file"
                                            accept=".csv"
                                            onChange={handleFileUpload}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        />
                                    </div>
                                    <Button
                                        onClick={loadExampleData}
                                        variant="outline"
                                        size="sm"
                                        className="w-full"
                                    >
                                        <Upload className="h-4 w-4 mr-2" />
                                        Load Example Data
                                    </Button>
                                    {csvFile && (
                                        <div className="glass rounded-lg p-3 border border-white/20">
                                            <p className="text-sm text-gray-700 dark:text-gray-300">
                                                <strong>File:</strong> {csvFile.name}
                                            </p>
                                            <p className="text-sm text-gray-600 dark:text-gray-400">
                                                Size: {(csvFile.size / 1024).toFixed(1)}KB | Rows: {csvData.length}
                                            </p>
                                        </div>
                                    )}
                                </CardContent>
                            </Card>

                            {/* Run Analysis Button */}
                            <Button
                                onClick={handleRunScoring}
                                disabled={isLoading || !csvData.length || !apiKey}
                                className="w-full bg-gradient-to-r from-purple-600 to-pink-600 hover:from-purple-700 hover:to-pink-700 text-white btn-modern"
                                size="lg"
                            >
                                {isLoading ? (
                                    <>
                                        <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2"></div>
                                        Running Analysis...
                                    </>
                                ) : (
                                    <>
                                        <Play className="h-4 w-4 mr-2" />
                                        Start Scoring Analysis
                                    </>
                                )}
                            </Button>
                        </div>
                    </div>

                    {/* Right Panel - Results */}
                    <div className="xl:col-span-3">
                        <div className="glass rounded-2xl p-6 border border-white/20">
                            <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                                📊 <span className="ml-2">Results</span>
                            </h2>

                            {isLoading && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 bg-gradient-to-r from-purple-600 to-pink-600 rounded-full flex items-center justify-center mx-auto mb-4 animate-pulse">
                                        <Target className="h-8 w-8 text-white" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Running Scoring Analysis
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300 mb-4">
                                        Evaluating annotation quality...
                                    </p>
                                    {progress.percentage > 0 && (
                                        <div className="glass rounded-lg p-4 border border-white/20">
                                            <div className="w-full bg-gray-200 rounded-full h-2 mb-2">
                                                <div
                                                    className="bg-gradient-to-r from-purple-600 to-pink-600 h-2 rounded-full transition-all duration-300"
                                                    style={{ width: `${progress.percentage}%` }}
                                                ></div>
                                            </div>
                                            <p className="text-sm text-gray-600 dark:text-gray-300">
                                                {progress.completed} / {progress.total} completed ({progress.percentage.toFixed(1)}%)
                                            </p>
                                        </div>
                                    )}
                                </div>
                            )}

                            {error && (
                                <div className="glass rounded-lg p-4 border border-red-400/30 bg-red-50/20">
                                    <h3 className="text-lg font-semibold text-red-800 dark:text-red-400 mb-2">Error</h3>
                                    <p className="text-red-700 dark:text-red-300">{error}</p>
                                </div>
                            )}

                            {results && (
                                <div className="space-y-6">
                                    <div className="glass rounded-lg p-4 border border-green-400/30 bg-green-50/20">
                                        <h3 className="text-lg font-semibold text-green-800 dark:text-green-400 mb-2">
                                            Analysis Complete!
                                        </h3>
                                        <p className="text-green-700 dark:text-green-300">
                                            Scoring analysis has been completed successfully.
                                        </p>
                                    </div>

                                    {/* Results Display */}
                                    <Card>
                                        <CardHeader>
                                            <CardTitle>📋 Scoring Results</CardTitle>
                                        </CardHeader>
                                        <CardContent>
                                            <div className="space-y-4">
                                                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.totalProcessed || 0}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">Total Processed</div>
                                                    </div>
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.averageScore?.toFixed(1) || 'N/A'}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">Average Score</div>
                                                    </div>
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.highQuality || 0}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">High Quality</div>
                                                    </div>
                                                </div>

                                                {results.csvContent && (
                                                    <Button
                                                        onClick={() => {
                                                            const blob = new Blob([results.csvContent], { type: 'text/csv' });
                                                            const url = URL.createObjectURL(blob);
                                                            const a = document.createElement('a');
                                                            a.href = url;
                                                            a.download = 'scoring_results.csv';
                                                            document.body.appendChild(a);
                                                            a.click();
                                                            document.body.removeChild(a);
                                                            URL.revokeObjectURL(url);
                                                        }}
                                                        variant="outline"
                                                        className="w-full"
                                                    >
                                                        <Download className="h-4 w-4 mr-2" />
                                                        Download Results CSV
                                                    </Button>
                                                )}
                                            </div>
                                        </CardContent>
                                    </Card>
                                </div>
                            )}

                            {!isLoading && !results && !error && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 glass rounded-full flex items-center justify-center mx-auto mb-4 border border-white/20">
                                        <Target className="h-8 w-8 text-gray-400" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Ready for Analysis
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300">
                                        Upload your annotation data and configure settings to start scoring analysis.
                                    </p>
                                </div>
                            )}
                        </div>
                    </div>
                </div>
            </main>
        </div>
    );
}