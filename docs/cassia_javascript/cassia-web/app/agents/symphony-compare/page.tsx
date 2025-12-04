"use client";

import { useState, ChangeEvent, useEffect } from 'react';
import Link from 'next/link';
import { symphonyCompare } from '@/lib/cassia/symphonyCompare';
import { useApiKeyStore } from '@/lib/stores/api-key-store';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { ArrowLeft, Play, HelpCircle, Music, Upload, Download, ChevronDown, ChevronUp } from 'lucide-react';
import { MODEL_PRESETS, DEFAULT_MODEL_PRESET, getModelPreset, getModelPersona } from '@/lib/config/model-presets';

// Component to parse and display model responses
const ParsedModelResponse = ({ response }: { response: string }) => {
    const parseModelResponse = (text: string) => {
        const cellTypes = [];
        const cellTypeRegex = /<celltype>(.*?)<\/celltype>\s*<reasoning>(.*?)<\/reasoning>\s*<score>(.*?)<\/score>/gs;
        let match;
        
        while ((match = cellTypeRegex.exec(text)) !== null) {
            cellTypes.push({
                cellType: match[1].trim(),
                reasoning: match[2].trim(),
                score: match[3].trim()
            });
        }
        
        return cellTypes;
    };

    const parsedData = parseModelResponse(response);
    
    if (parsedData.length === 0) {
        // Fallback to raw display if parsing fails
        return (
            <div className="bg-gray-100 dark:bg-gray-800 p-3 rounded text-xs font-mono max-h-32 overflow-y-auto">
                {response}
            </div>
        );
    }

    return (
        <div className="space-y-3">
            {parsedData.map((data, index) => (
                <div key={index} className="border border-gray-200 dark:border-gray-700 rounded-lg p-4 bg-white/50 dark:bg-gray-800/50">
                    <div className="flex items-center justify-between mb-2">
                        <h4 className="font-semibold text-gray-900 dark:text-white text-sm">
                            {data.cellType}
                        </h4>
                        <span className="text-lg font-bold text-indigo-600 dark:text-indigo-400">
                            {data.score}
                        </span>
                    </div>
                    <p className="text-sm text-gray-600 dark:text-gray-300 italic leading-relaxed">
                        {data.reasoning}
                    </p>
                </div>
            ))}
        </div>
    );
};

export default function SymphonyComparePage() {
    const globalApiKey = useApiKeyStore((state) => state.getApiKey());
    const globalProvider = useApiKeyStore((state) => state.provider);
    
    // State management
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState('openrouter');
    const [model, setModel] = useState('google/gemini-2.5-flash');
    const [markerGenes, setMarkerGenes] = useState('');
    const [candidateCellTypes, setCandidateCellTypes] = useState('');
    const [majorClusterInfo, setMajorClusterInfo] = useState('Human PBMC');
    const [temperature, setTemperature] = useState(0);
    const [additionalTask, setAdditionalTask] = useState('');
    const [modelPreset, setModelPreset] = useState(DEFAULT_MODEL_PRESET);
    const [selectedModels, setSelectedModels] = useState<string[]>([]);
    const [customModels, setCustomModels] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const [results, setResults] = useState<any>(null);
    const [error, setError] = useState<string | null>(null);
    const [showDetailedResults, setShowDetailedResults] = useState(false);
    const [showDebugPanel, setShowDebugPanel] = useState(false);

    // Initialize with global settings
    useEffect(() => {
        if (globalApiKey && !apiKey) {
            setApiKey(globalApiKey);
        }
        // Symphony Compare only supports OpenRouter for multi-model access
        setProvider('openrouter');
    }, [globalApiKey, globalProvider]);

    // Initialize selected models when preset changes
    useEffect(() => {
        if (modelPreset !== 'custom') {
            const preset = getModelPreset(modelPreset);
            if (preset) {
                setSelectedModels(preset.models.slice(0, 3)); // Select first 3 models by default
            }
        }
    }, [modelPreset]);

    const loadExampleData = () => {
        setMarkerGenes('CD3D, CD3E, CD3G, IL7R, CCR7, CD8A, CD8B, PRF1, GZMB, GNLY, NKG7, FCGR3A, MS4A1, CD79A, CD79B, IGHM, CD14, LYZ, S100A9, FCGR1A, CST3');
        setCandidateCellTypes('T cells\nB cells\nNK cells\nMonocytes');
        setMajorClusterInfo('Human PBMC');
        console.log('Example data loaded successfully');
    };

    const handleModelToggle = (model: string) => {
        setSelectedModels(prev => {
            if (prev.includes(model)) {
                return prev.filter(m => m !== model);
            } else {
                return [...prev, model];
            }
        });
    };

    const handleRunSymphony = async () => {
        if (!markerGenes.trim() || !candidateCellTypes.trim() || !apiKey) {
            setError('Please provide marker genes, candidate cell types, and API key');
            return;
        }

        // Parse candidate cell types
        const cellTypes = candidateCellTypes
            .split('\n')
            .map(type => type.trim())
            .filter(type => type.length > 0);

        if (cellTypes.length < 2 || cellTypes.length > 4) {
            setError('Please provide 2-4 candidate cell types');
            return;
        }

        // Parse custom models if using custom preset, otherwise use selected models
        let customModelsList = null;
        if (modelPreset === 'custom') {
            customModelsList = customModels
                .split('\n')
                .map(model => model.trim())
                .filter(model => model.length > 0);
            
            if (customModelsList.length === 0) {
                setError('Please provide custom models when using custom preset');
                return;
            }
        } else {
            // Use selected models from the preset
            customModelsList = selectedModels;
            
            if (customModelsList.length < 2) {
                setError('Please select at least 2 models from the preset');
                return;
            }
        }

        setIsLoading(true);
        setError(null);
        setResults(null);

        try {
            const result = await symphonyCompare({
                apiKey: apiKey,
                provider: provider,
                tissue: majorClusterInfo,
                species: 'human',
                celltypes: cellTypes,
                markerSet: markerGenes,
                modelPreset: "custom", // Always use custom since we're providing selected models
                customModels: customModelsList,
                enableDiscussion: true,
                maxDiscussionRounds: 2,
                consensusThreshold: 0.8
            });
            
            setResults(result);
            console.log('Symphony Compare completed:', result);
        } catch (err: any) {
            setError(`Symphony Compare failed: ${err.message}`);
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
                                    <div className="w-12 h-12 bg-gradient-to-br from-indigo-600 to-cyan-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                                        <Music className="h-6 w-6 text-white" />
                                    </div>
                                    <div className="absolute -top-1 -right-1 w-4 h-4 bg-purple-400 rounded-full animate-pulse"></div>
                                </div>
                                <div>
                                    <h1 className="text-3xl font-bold gradient-text">Symphony Compare</h1>
                                    <p className="text-gray-600 dark:text-gray-300">Multi-model consensus analysis</p>
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
                {/* Info Panel */}
                <div className="glass rounded-2xl p-8 mb-12 border border-indigo-400/30 bg-indigo-50/20">
                    <div className="flex items-start space-x-4">
                        <div className="w-12 h-12 bg-gradient-to-br from-indigo-600 to-cyan-600 rounded-xl flex items-center justify-center shadow-lg">
                            <Music className="h-6 w-6 text-white" />
                        </div>
                        <div>
                            <h3 className="text-xl font-semibold text-gray-900 dark:text-white mb-2">Symphony Compare</h3>
                            <p className="text-gray-600 dark:text-gray-300 leading-relaxed">
                                Symphony Compare orchestrates multiple AI models to compare candidate cell types against your marker genes. 
                                Models engage in discussion rounds to reach consensus, providing confidence scores and detailed analysis 
                                for each candidate cell type to determine the best match.
                            </p>
                        </div>
                    </div>
                </div>

                <div className="grid grid-cols-1 xl:grid-cols-5 gap-8">
                    {/* Left Panel - Configuration */}
                    <div className="xl:col-span-2">
                        <div className="glass rounded-2xl p-6 border border-white/20 sticky top-24 space-y-6">
                            <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                                🎛️ <span className="ml-2">Configuration</span>
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
                                            disabled
                                        >
                                            <option value="openrouter">OpenRouter (Required for Symphony Compare)</option>
                                        </select>
                                        <div className="text-xs text-blue-600 dark:text-blue-400 mt-1">
                                            💡 Symphony Compare requires OpenRouter to access multiple models for consensus analysis
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Model Configuration */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🤖 Model Configuration</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Model Preset</label>
                                        <select
                                            value={modelPreset}
                                            onChange={(e) => setModelPreset(e.target.value)}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        >
                                            {Object.entries(MODEL_PRESETS).map(([key, preset]) => (
                                                <option key={key} value={key}>
                                                    {preset.name} ({preset.description})
                                                </option>
                                            ))}
                                            <option value="custom">Custom Models</option>
                                        </select>
                                        <div className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                            Choose a preset or use custom models
                                        </div>
                                    </div>
                                    {modelPreset === 'custom' && (
                                        <div>
                                            <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Custom Models</label>
                                            <textarea
                                                value={customModels}
                                                onChange={(e) => setCustomModels(e.target.value)}
                                                placeholder="Enter model names, one per line..."
                                                className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20 min-h-[80px] resize-y"
                                            />
                                            <div className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                                Example: anthropic/claude-3-sonnet
                                            </div>
                                        </div>
                                    )}
                                    
                                    {/* Model Selection */}
                                    {modelPreset !== 'custom' && (
                                        <div className="glass rounded-lg p-3 border border-white/20 bg-blue-50/20">
                                            <h5 className="font-medium text-gray-900 dark:text-white mb-3">
                                                🎭 Select Models ({selectedModels.length} selected):
                                            </h5>
                                            <div className="space-y-2 max-h-48 overflow-y-auto">
                                                {(() => {
                                                    const preset = getModelPreset(modelPreset);
                                                    const availableModels = preset?.models || [];
                                                    return availableModels.map((model, index) => (
                                                        <label key={index} className="flex items-center space-x-2 cursor-pointer hover:bg-white/20 dark:hover:bg-gray-700/20 rounded p-2 transition-colors">
                                                            <input
                                                                type="checkbox"
                                                                checked={selectedModels.includes(model)}
                                                                onChange={() => handleModelToggle(model)}
                                                                className="rounded border-gray-300 text-indigo-600 focus:ring-indigo-500"
                                                            />
                                                            <span className="text-sm text-gray-700 dark:text-gray-300 flex-1">
                                                                <div className="font-medium">{getModelPersona(model)}</div>
                                                                <div className="text-xs text-gray-500 dark:text-gray-400">{model}</div>
                                                            </span>
                                                        </label>
                                                    ));
                                                })()}
                                            </div>
                                            {selectedModels.length < 2 && (
                                                <div className="text-xs text-amber-600 dark:text-amber-400 mt-2">
                                                    ⚠️ Please select at least 2 models
                                                </div>
                                            )}
                                        </div>
                                    )}

                                    {/* Show custom models input */}
                                    {modelPreset === 'custom' && (
                                        <div className="glass rounded-lg p-3 border border-white/20 bg-blue-50/20">
                                            <h5 className="font-medium text-gray-900 dark:text-white mb-2">
                                                🎭 Custom Models:
                                            </h5>
                                            <div className="space-y-1">
                                                {customModels.split('\n').filter(m => m.trim()).length > 0 ? (
                                                    customModels.split('\n').filter(m => m.trim()).map((model, index) => (
                                                        <div key={index} className="text-sm text-gray-700 dark:text-gray-300">
                                                            • {model.trim()}
                                                        </div>
                                                    ))
                                                ) : (
                                                    <div className="text-sm text-gray-500 dark:text-gray-400 italic">
                                                        No custom models specified
                                                    </div>
                                                )}
                                            </div>
                                        </div>
                                    )}
                                </CardContent>
                            </Card>

                            {/* Analysis Settings */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🎯 Analysis Settings</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Sample Context</label>
                                        <Input
                                            type="text"
                                            value={majorClusterInfo}
                                            onChange={(e) => setMajorClusterInfo(e.target.value)}
                                            placeholder="e.g., Human PBMC, Mouse Brain"
                                        />
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">
                                            Temperature: {temperature}
                                        </label>
                                        <input
                                            type="range"
                                            min="0"
                                            max="1"
                                            step="0.1"
                                            value={temperature}
                                            onChange={(e) => setTemperature(Number(e.target.value))}
                                            className="w-full"
                                        />
                                        <div className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                            0 = Focused, 1 = Creative
                                        </div>
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Additional Task (Optional)</label>
                                        <textarea
                                            value={additionalTask}
                                            onChange={(e) => setAdditionalTask(e.target.value)}
                                            placeholder="Any additional analysis instructions..."
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20 min-h-[80px] resize-y"
                                        />
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Candidate Cell Types Input */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🔬 Candidate Cell Types</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Cell Types to Compare (2-4 types)</label>
                                        <textarea
                                            value={candidateCellTypes}
                                            onChange={(e) => setCandidateCellTypes(e.target.value)}
                                            placeholder="Enter cell types, one per line..."
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20 min-h-[100px] resize-y"
                                        />
                                        <div className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                            Example: T cells, B cells, NK cells, Monocytes
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Marker Genes Input */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🧬 Marker Genes</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Marker Genes</label>
                                        <textarea
                                            value={markerGenes}
                                            onChange={(e) => setMarkerGenes(e.target.value)}
                                            placeholder="Enter marker genes separated by commas..."
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20 min-h-[120px] resize-y"
                                        />
                                        <div className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                            Example: CD3D, CD3E, CD3G, IL7R, CCR7
                                        </div>
                                    </div>
                                    <Button
                                        onClick={loadExampleData}
                                        variant="outline"
                                        size="sm"
                                        className="w-full"
                                    >
                                        <Upload className="h-4 w-4 mr-2" />
                                        Load Example Genes
                                    </Button>
                                </CardContent>
                            </Card>

                            {/* Run Analysis Button */}
                            <Button
                                onClick={handleRunSymphony}
                                disabled={isLoading || !markerGenes.trim() || !candidateCellTypes.trim() || !apiKey || 
                                         (modelPreset !== 'custom' && selectedModels.length < 2) ||
                                         (modelPreset === 'custom' && !customModels.trim())}
                                className="w-full bg-gradient-to-r from-indigo-600 to-cyan-600 hover:from-indigo-700 hover:to-cyan-700 text-white btn-modern"
                                size="lg"
                            >
                                {isLoading ? (
                                    <>
                                        <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-white mr-2"></div>
                                        Running Symphony...
                                    </>
                                ) : (
                                    <>
                                        <Play className="h-4 w-4 mr-2" />
                                        Start Symphony Compare
                                    </>
                                )}
                            </Button>
                        </div>
                    </div>

                    {/* Right Panel - Results */}
                    <div className="xl:col-span-3">
                        <div className="glass rounded-2xl p-6 border border-white/20">
                            <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                                🎼 <span className="ml-2">Symphony Results</span>
                            </h2>

                            {isLoading && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 bg-gradient-to-r from-indigo-600 to-cyan-600 rounded-full flex items-center justify-center mx-auto mb-4 animate-pulse">
                                        <Music className="h-8 w-8 text-white" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Running Symphony Compare
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300 mb-4">
                                        Orchestrating multiple AI models for consensus analysis...
                                    </p>
                                    <div className="glass rounded-lg p-4 border border-white/20">
                                        <div className="flex justify-center space-x-2">
                                            <div className="w-2 h-2 bg-indigo-600 rounded-full animate-bounce"></div>
                                            <div className="w-2 h-2 bg-cyan-600 rounded-full animate-bounce" style={{animationDelay: '0.1s'}}></div>
                                            <div className="w-2 h-2 bg-purple-600 rounded-full animate-bounce" style={{animationDelay: '0.2s'}}></div>
                                        </div>
                                        <p className="text-sm text-gray-600 dark:text-gray-300 mt-2">
                                            Comparing model outputs...
                                        </p>
                                    </div>
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
                                            Symphony Complete! 🎵
                                        </h3>
                                        <p className="text-green-700 dark:text-green-300">
                                            Multi-model consensus analysis has been completed successfully.
                                        </p>
                                    </div>

                                    {/* Consensus/Winner Display - PROMINENT */}
                                    {(results.consensus || results.summary?.final_consensus) && (
                                        <Card className="border-2 border-indigo-500 shadow-lg bg-gradient-to-br from-indigo-50/50 to-cyan-50/50 dark:from-indigo-950/30 dark:to-cyan-950/30">
                                            <CardHeader className="pb-3">
                                                <CardTitle className="text-2xl flex items-center justify-center">
                                                    🏆 <span className="ml-2 gradient-text">Consensus Cell Type</span>
                                                </CardTitle>
                                            </CardHeader>
                                            <CardContent>
                                                <div className="text-center mb-6">
                                                    <h2 className="text-4xl font-bold text-indigo-600 dark:text-indigo-400 mb-2">
                                                        {results.consensus || results.summary?.final_consensus}
                                                    </h2>
                                                    {results.confidence !== undefined && (
                                                        <div className="flex items-center justify-center gap-4">
                                                            <div className="text-sm text-gray-600 dark:text-gray-400">
                                                                Confidence: {Math.round((results.confidence || 0) * 100)}%
                                                            </div>
                                                            <div className="w-32 bg-gray-200 rounded-full h-2">
                                                                <div
                                                                    className="bg-gradient-to-r from-indigo-600 to-cyan-600 h-2 rounded-full transition-all duration-300"
                                                                    style={{ width: `${(results.confidence || 0) * 100}%` }}
                                                                ></div>
                                                            </div>
                                                        </div>
                                                    )}
                                                </div>

                                                {/* Model Votes and Reasoning for Winner */}
                                                <div className="space-y-4">
                                                    <h4 className="font-semibold text-gray-900 dark:text-white text-center">
                                                        Model Analysis for {results.consensus || results.summary?.final_consensus}
                                                    </h4>
                                                    {results.raw_responses && results.raw_responses
                                                        .filter(resp => resp.round === 'initial')
                                                        .map((resp, index) => {
                                                            const winnerCellType = results.consensus || results.summary?.final_consensus || '';
                                                            const cellData = resp.extracted[winnerCellType];
                                                            if (!cellData) return null;
                                                            
                                                            return (
                                                                <div key={index} className="glass rounded-lg p-4 border border-white/20">
                                                                    <div className="flex items-center justify-between mb-2">
                                                                        <h5 className="font-semibold text-gray-900 dark:text-white">
                                                                            {resp.researcher}
                                                                        </h5>
                                                                        <span className="text-2xl font-bold gradient-text">
                                                                            {cellData.score}
                                                                        </span>
                                                                    </div>
                                                                    <p className="text-sm text-gray-600 dark:text-gray-300 italic">
                                                                        "{cellData.reasoning}"
                                                                    </p>
                                                                </div>
                                                            );
                                                        })}
                                                </div>
                                            </CardContent>
                                        </Card>
                                    )}

                                    {/* No Consensus - Show Top Voted */}
                                    {!results.consensus && !results.summary?.final_consensus && results.summary?.celltype_scores && (
                                        <Card className="border-2 border-amber-500 shadow-lg bg-gradient-to-br from-amber-50/50 to-orange-50/50 dark:from-amber-950/30 dark:to-orange-950/30">
                                            <CardHeader className="pb-3">
                                                <CardTitle className="text-2xl flex items-center justify-center">
                                                    ⚖️ <span className="ml-2 gradient-text">No Consensus Reached</span>
                                                </CardTitle>
                                            </CardHeader>
                                            <CardContent>
                                                <div className="text-center mb-4">
                                                    <p className="text-gray-600 dark:text-gray-300">
                                                        Models did not reach consensus. Here are the votes:
                                                    </p>
                                                </div>
                                                <div className="space-y-3">
                                                    {Object.entries(results.summary.celltype_scores)
                                                        .sort(([,a]: any, [,b]: any) => (b.mean || 0) - (a.mean || 0))
                                                        .map(([celltype, scores]: [string, any]) => (
                                                            <div key={celltype} className="glass rounded-lg p-3 border border-white/20">
                                                                <div className="flex items-center justify-between">
                                                                    <h5 className="font-semibold text-gray-900 dark:text-white">
                                                                        {celltype}
                                                                    </h5>
                                                                    <div className="text-right">
                                                                        <span className="text-lg font-bold gradient-text">
                                                                            {scores.mean?.toFixed(1) || 'N/A'}
                                                                        </span>
                                                                        <span className="text-sm text-gray-600 dark:text-gray-400 ml-2">
                                                                            avg score
                                                                        </span>
                                                                    </div>
                                                                </div>
                                                                {scores.votes !== undefined && (
                                                                    <p className="text-sm text-gray-600 dark:text-gray-300 mt-1">
                                                                        Votes: {scores.votes} / {results.summary.models_used || 0} models
                                                                    </p>
                                                                )}
                                                            </div>
                                                        ))}
                                                </div>
                                            </CardContent>
                                        </Card>
                                    )}

                                    {/* Results Display - Collapsible */}
                                    <Card>
                                        <CardHeader className="cursor-pointer hover:bg-gray-50 dark:hover:bg-gray-800 transition-colors" onClick={() => setShowDetailedResults(!showDetailedResults)}>
                                            <CardTitle className="flex items-center justify-between">
                                                <span className="flex items-center">
                                                    📊 <span className="ml-2">Detailed Analysis</span>
                                                </span>
                                                <span className="flex items-center text-sm font-normal text-gray-600 dark:text-gray-400">
                                                    {showDetailedResults ? 'Hide' : 'Show'}
                                                    {showDetailedResults ? <ChevronUp className="h-4 w-4 ml-1" /> : <ChevronDown className="h-4 w-4 ml-1" />}
                                                </span>
                                            </CardTitle>
                                        </CardHeader>
                                        {showDetailedResults && (
                                            <CardContent>
                                                <div className="space-y-4">
                                                {/* Summary Stats */}
                                                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.summary?.models_used || 0}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">Models Consulted</div>
                                                    </div>
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.summary?.total_rounds || 0}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">Total Rounds</div>
                                                    </div>
                                                    <div className="glass rounded-lg p-4 border border-white/20 text-center">
                                                        <div className="text-2xl font-bold gradient-text">
                                                            {results.summary?.consensus_reached ? '✅' : '❌'}
                                                        </div>
                                                        <div className="text-sm text-gray-600 dark:text-gray-300">Consensus Reached</div>
                                                    </div>
                                                </div>

                                                {results.consensus && (
                                                    <div className="glass rounded-lg p-4 border border-white/20">
                                                        <h4 className="font-semibold text-gray-900 dark:text-white mb-2">
                                                            🏆 Consensus Cell Type
                                                        </h4>
                                                        <p className="text-lg gradient-text font-semibold">
                                                            {results.consensus}
                                                        </p>
                                                    </div>
                                                )}

                                                {results.confidence !== undefined && (
                                                    <div className="glass rounded-lg p-4 border border-white/20">
                                                        <h4 className="font-semibold text-gray-900 dark:text-white mb-2">
                                                            📊 Confidence Score
                                                        </h4>
                                                        <div className="w-full bg-gray-200 rounded-full h-3 mb-2">
                                                            <div
                                                                className="bg-gradient-to-r from-indigo-600 to-cyan-600 h-3 rounded-full transition-all duration-300"
                                                                style={{ width: `${results.confidence * 100}%` }}
                                                            ></div>
                                                        </div>
                                                        <p className="text-sm text-gray-600 dark:text-gray-300">
                                                            {Math.round(results.confidence * 100)}% confidence
                                                        </p>
                                                    </div>
                                                )}

                                                {results.summary?.celltype_scores && (
                                                    <div className="glass rounded-lg p-4 border border-white/20">
                                                        <h4 className="font-semibold text-gray-900 dark:text-white mb-3">
                                                            📈 Cell Type Scores
                                                        </h4>
                                                        <div className="space-y-3">
                                                            {Object.entries(results.summary.celltype_scores).map(([celltype, scores]: [string, any]) => (
                                                                <div key={celltype} className="border-l-4 border-indigo-500 pl-4">
                                                                    <h5 className="font-medium text-gray-900 dark:text-white">
                                                                        {celltype}
                                                                    </h5>
                                                                    <p className="text-sm text-gray-600 dark:text-gray-300">
                                                                        Average: {scores.mean?.toFixed(1) || 'N/A'} 
                                                                        {scores.std !== undefined && ` (±${scores.std.toFixed(1)})`}
                                                                        {scores.min !== undefined && scores.max !== undefined && ` [${scores.min}-${scores.max}]`}
                                                                    </p>
                                                                </div>
                                                            ))}
                                                        </div>
                                                    </div>
                                                )}

                                                {/* Debug Panel - Raw Responses */}
                                                {results.raw_responses && (
                                                    <Card>
                                                        <CardHeader className="cursor-pointer hover:bg-gray-50 dark:hover:bg-gray-800 transition-colors" onClick={() => setShowDebugPanel(!showDebugPanel)}>
                                                            <CardTitle className="flex items-center justify-between">
                                                                <span className="flex items-center">
                                                                    🔍 <span className="ml-2">Model Responses</span>
                                                                </span>
                                                                <span className="flex items-center text-sm font-normal text-gray-600 dark:text-gray-400">
                                                                    {showDebugPanel ? 'Hide' : 'Show'}
                                                                    {showDebugPanel ? <ChevronUp className="h-4 w-4 ml-1" /> : <ChevronDown className="h-4 w-4 ml-1" />}
                                                                </span>
                                                            </CardTitle>
                                                        </CardHeader>
                                                        {showDebugPanel && (
                                                            <CardContent>
                                                                <div className="space-y-4 max-h-96 overflow-y-auto">
                                                                {results.raw_responses.map((resp, index) => (
                                                                    <div key={index} className="glass rounded-lg p-4 border border-white/20">
                                                                        <h5 className="font-semibold text-gray-900 dark:text-white mb-2">
                                                                            {getModelPersona(resp.model)}
                                                                        </h5>
                                                                        <div className="text-xs text-gray-600 dark:text-gray-300 mb-2">
                                                                            Model: {resp.model} | Round: {resp.round}
                                                                        </div>
                                                                        <div className="space-y-3">
                                                                            <ParsedModelResponse response={resp.response} />
                                                                        </div>
                                                                    </div>
                                                                ))}
                                                                </div>
                                                            </CardContent>
                                                        )}
                                                    </Card>
                                                )}

                                                {results.reportHtml && (
                                                    <Button
                                                        onClick={() => {
                                                            const blob = new Blob([results.reportHtml], { type: 'text/html' });
                                                            const url = URL.createObjectURL(blob);
                                                            const a = document.createElement('a');
                                                            a.href = url;
                                                            a.download = 'symphony_analysis.html';
                                                            document.body.appendChild(a);
                                                            a.click();
                                                            document.body.removeChild(a);
                                                            URL.revokeObjectURL(url);
                                                        }}
                                                        variant="outline"
                                                        className="w-full"
                                                    >
                                                        <Download className="h-4 w-4 mr-2" />
                                                        Download Analysis Report
                                                    </Button>
                                                )}
                                            </div>
                                        </CardContent>
                                    )}
                                </Card>
                                </div>
                            )}

                            {!isLoading && !results && !error && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 glass rounded-full flex items-center justify-center mx-auto mb-4 border border-white/20">
                                        <Music className="h-8 w-8 text-gray-400" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Ready for Symphony
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300">
                                        Enter your marker genes, candidate cell types (2-4), and configuration to start the multi-model consensus analysis.
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