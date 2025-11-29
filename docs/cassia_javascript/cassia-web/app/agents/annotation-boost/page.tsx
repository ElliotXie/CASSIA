"use client";

import { useState, useEffect } from 'react';
import Link from 'next/link';
import { useDropzone } from 'react-dropzone';
import Papa from 'papaparse';
import { iterativeMarkerAnalysis, generateSummaryReport, extractConversationForCluster, getAvailableClusters, extractTopMarkerGenes } from '@/lib/cassia/annotationBoost';
import { useApiKeyStore } from '@/lib/stores/api-key-store';
import { useConfigStore } from '@/lib/stores/config-store';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Badge } from '@/components/ui/badge';
import { cn } from '@/lib/utils';
import { ArrowLeft, Play, HelpCircle, Zap, Upload, Download, FileText, History, RefreshCcw, BookOpen, Database, AlertCircle, CheckCircle, Loader2, Sparkles, Settings } from 'lucide-react';
import modelSettings from '../../../public/examples/model_settings.json';

export default function AnnotationBoostPage() {
    const [isInitialized, setIsInitialized] = useState(false);
    
    // Use safer store access with error handling
    const globalApiKey = useApiKeyStore((state) => {
        try {
            return state.getApiKey();
        } catch (error) {
            console.warn('Error accessing global API key:', error);
            return '';
        }
    });
    const globalProvider = useApiKeyStore((state) => {
        try {
            return state.provider;
        } catch (error) {
            console.warn('Error accessing global provider:', error);
            return 'openrouter';
        }
    });
    const globalModel = useConfigStore((state) => {
        try {
            return state.model;
        } catch (error) {
            console.warn('Error accessing global model:', error);
            return modelSettings.use_case_recommendations.annotation_boost.best;
        }
    });
    
    // State management
    const [conversationFile, setConversationFile] = useState<File | null>(null);
    const [conversationData, setConversationData] = useState<string>('');
    const [markerFile, setMarkerFile] = useState<File | null>(null);
    const [markerData, setMarkerData] = useState<any[]>([]);
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState('openrouter');
    const [model, setModel] = useState(modelSettings.use_case_recommendations.annotation_boost.best);
    const [majorClusterInfo, setMajorClusterInfo] = useState('Human PBMC');
    const [numIterations, setNumIterations] = useState(3);
    const [searchStrategy, setSearchStrategy] = useState('breadth');
    const [selectedMode, setSelectedMode] = useState<'performance' | 'balanced'>('balanced');
    const [customModelId, setCustomModelId] = useState('');
    
    // Provider-specific models from model_settings.json
    const getModelsByProvider = (provider: string) => {
        const providerData = modelSettings.providers[provider as keyof typeof modelSettings.providers];
        if (!providerData || !providerData.models) {
            return [];
        }
        
        return Object.entries(providerData.models).map(([key, model]) => ({
            value: model.actual_name,
            label: model.description.split(' - ')[0] || model.actual_name
        }));
    };

    const getDefaultModelForProvider = (provider: string) => {
        const providerData = modelSettings.providers[provider as keyof typeof modelSettings.providers];
        return providerData?.default_model || modelSettings.use_case_recommendations.annotation_boost.best;
    };
    const [historySource, setHistorySource] = useState<'csv' | 'manual'>('csv');
    const [manualHistory, setManualHistory] = useState('');
    const [selectedCluster, setSelectedCluster] = useState('');
    const [availableClusters, setAvailableClusters] = useState<string[]>([]);
    const [extractedHistory, setExtractedHistory] = useState('');
    const [availableColumns, setAvailableColumns] = useState<string[]>([]);
    const [selectedClusterColumn, setSelectedClusterColumn] = useState('True Cell Type');
    const [isLoading, setIsLoading] = useState(false);
    const [resultsHtml, setResultsHtml] = useState('');
    const [error, setError] = useState<string | null>(null);
    const [consoleOutput, setConsoleOutput] = useState<string[]>([]);
    const [showInstructions, setShowInstructions] = useState(true);

    // Console output capture
    useEffect(() => {
        try {
            const originalLog = console.log;
            console.log = function(...args) {
                try {
                    originalLog.apply(console, args);
                    const message = args.map(arg => 
                        typeof arg === 'string' ? arg : JSON.stringify(arg)
                    ).join(' ');
                    setConsoleOutput(prev => [...prev.slice(-19), `${new Date().toLocaleTimeString()}: ${message}`]);
                } catch (error) {
                    originalLog('Error in console capture:', error);
                }
            };
            
            return () => {
                console.log = originalLog;
            };
        } catch (error) {
            console.error('Error setting up console capture:', error);
        }
    }, []);

    // Initialize with global settings
    useEffect(() => {
        try {
            if (!isInitialized) {
                if (globalApiKey && !apiKey) {
                    setApiKey(globalApiKey);
                }
                if (globalProvider && provider !== globalProvider) {
                    setProvider(globalProvider);
                }
                // For annotation boost, we prefer Gemini 2.5 Flash as default
                // Only override if there's a different saved preference
                const savedModel = localStorage.getItem('annotation-boost-model');
                if (savedModel && savedModel !== model) {
                    setModel(savedModel);
                } else if (!savedModel) {
                    // Ensure Gemini 2.5 Flash is the default for new users
                    setModel('google/gemini-2.5-flash');
                }
                setIsInitialized(true);
            }
        } catch (error) {
            console.error('Error initializing global settings:', error);
            setIsInitialized(true); // Set to true even on error to prevent infinite loops
        }
    }, [globalApiKey, globalProvider, globalModel, apiKey, provider, model, isInitialized]);

    // File upload handlers
    const onConversationDrop = (acceptedFiles: File[]) => {
        const file = acceptedFiles[0];
        if (file) {
            setConversationFile(file);
            const reader = new FileReader();
            reader.onload = (e) => {
                const text = e.target?.result as string;
                setConversationData(text);
                
                // Parse and get available columns
                Papa.parse(text, {
                    header: true,
                    skipEmptyLines: true,
                    complete: (result) => {
                        if (result.data.length > 0) {
                            const columns = Object.keys(result.data[0] as any);
                            setAvailableColumns(columns);
                        }
                    }
                });
            };
            reader.readAsText(file);
        }
    };

    const onMarkerDrop = (acceptedFiles: File[]) => {
        const file = acceptedFiles[0];
        if (file) {
            setMarkerFile(file);
            Papa.parse(file, {
                header: true,
                skipEmptyLines: true,
                complete: (result) => {
                    setMarkerData(result.data as any[]);
                },
                error: (err) => {
                    setError(`Error parsing marker CSV: ${err.message}`);
                }
            });
        }
    };

    const { 
        getRootProps: getConversationRootProps, 
        getInputProps: getConversationInputProps, 
        isDragActive: isConversationDragActive 
    } = useDropzone({ 
        onDrop: onConversationDrop, 
        accept: { 'text/csv': ['.csv'] },
        multiple: false
    });

    const { 
        getRootProps: getMarkerRootProps, 
        getInputProps: getMarkerInputProps, 
        isDragActive: isMarkerDragActive 
    } = useDropzone({ 
        onDrop: onMarkerDrop, 
        accept: { 'text/csv': ['.csv'] },
        multiple: false
    });

    // Load available clusters
    const loadClusters = async () => {
        if (!conversationData || !selectedClusterColumn) return;
        
        try {
            const clusters = await getAvailableClusters(conversationData, selectedClusterColumn);
            setAvailableClusters(clusters);
        } catch (err: any) {
            setError(`Failed to load clusters: ${err.message}`);
        }
    };

    // Extract conversation for selected cluster
    const extractConversation = async () => {
        if (!conversationData || !selectedCluster || !selectedClusterColumn) return;
        
        try {
            const extracted = await extractConversationForCluster(conversationData, selectedCluster, selectedClusterColumn);
            setExtractedHistory(extracted);
        } catch (err: any) {
            setError(`Failed to extract conversation: ${err.message}`);
        }
    };

    // Load example data by simulating file upload
    const loadExampleData = async () => {
        try {
            setError(null);
            console.log('🔄 Loading example conversation history...');
            
            // Try to load from public directory (standard Next.js static serving)
            let conversationText, markerText;
            
            try {
                const conversationResponse = await fetch('/examples/cassia_analysis_full.csv');
                if (conversationResponse.ok) {
                    conversationText = await conversationResponse.text();
                } else {
                    throw new Error('Static file fetch failed');
                }
            } catch {
                console.log('📝 Static file loading failed, using embedded fallback...');
                // Fallback: Use a small sample of the conversation data
                conversationText = `"True Cell Type","Predicted Main Cell Type","Predicted Sub Cell Types","Possible Mixed Cell Types","Marker Number","Marker List","Iterations","Model","Provider","Tissue","Species","Additional Info","Conversation History"
"monocyte","Monocytes","Classical Monocytes, Intermediate Monocytes","","50","CDH19, NRXN1, PLP1, LGI4, GPM6B, SCN7A, SPARCL1, CRYAB, PMP22, IGFBP7, SPP1, CLU, MYOT, TIMP3, GPX3, CHL1, NRXN3, NCAM1, GPR155, SAMHD1","2","anthropic/claude-3.5-sonnet","openrouter","peripheral_blood","human","Example analysis","Example conversation history for monocyte analysis demonstrating iterative marker analysis workflow."`;
            }
            
            console.log('🔄 Loading example marker data...');
            try {
                const markerResponse = await fetch('/examples/batch_raw_seurat_example.csv');
                if (markerResponse.ok) {
                    markerText = await markerResponse.text();
                } else {
                    throw new Error('Static file fetch failed');
                }
            } catch {
                console.log('📝 Static marker file loading failed, using embedded fallback...');
                // Fallback: Use a small sample of marker data
                markerText = `"","p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene"
"CDH19",0,8.86571097085819,0.996,0.014,0,"monocyte","CDH19"
"NRXN1",0,8.41654303243543,0.974,0.013,0,"monocyte","NRXN1"
"PLP1",0,8.58194567864456,0.917,0.01,0,"monocyte","PLP1"
"LGI4",0,6.89879043296571,0.958,0.059,0,"monocyte","LGI4"
"GPM6B",0,7.23456789012345,0.892,0.023,0,"monocyte","GPM6B"`;
            }
            
            // Create File objects (simulating real file upload)
            const conversationBlob = new Blob([conversationText], { type: 'text/csv' });
            const markerBlob = new Blob([markerText], { type: 'text/csv' });
            
            const conversationFile = new File([conversationBlob], 'cassia_analysis_full.csv', { type: 'text/csv' });
            const markerFile = new File([markerBlob], 'batch_raw_seurat_example.csv', { type: 'text/csv' });
            
            console.log('📁 Simulating file upload...');
            
            // Simulate the exact file drop process
            onConversationDrop([conversationFile]);
            onMarkerDrop([markerFile]);
            
            console.log('🎉 Example data loaded successfully!');
            
        } catch (err: any) {
            setError(`Failed to load example data: ${err.message}`);
            console.error('❌ Error loading example data:', err);
        }
    };

    // Clear console output
    const clearConsole = () => {
        setConsoleOutput([]);
        console.log('🧹 Console cleared');
    };

    // Run iterative analysis
    const handleRunAnalysis = async () => {
        if (!markerData.length || !apiKey) {
            setError('Please upload marker data and provide API key');
            return;
        }

        const historyToUse = historySource === 'csv' ? extractedHistory : manualHistory;
        if (!historyToUse.trim()) {
            setError('Please provide conversation history');
            return;
        }

        setIsLoading(true);
        setError(null);
        setResultsHtml('');
        console.log('🚀 Starting Annotation Boost analysis...');

        try {
            const results = await iterativeMarkerAnalysis(
                majorClusterInfo,
                markerData,
                extractTopMarkerGenes(markerData, 20), // Extract top genes as comma-separated string
                historyToUse,
                numIterations,
                provider,
                model,
                null, // additionalTask
                0, // temperature
                searchStrategy as 'breadth' | 'depth',
                apiKey
            );

            const htmlReport = await generateSummaryReport(
                results.messages,
                searchStrategy as 'breadth' | 'depth',
                'per_iteration',
                provider,
                model,
                apiKey
            );
            setResultsHtml(htmlReport);
        } catch (err: any) {
            setError(`Analysis failed: ${err.message}`);
        } finally {
            setIsLoading(false);
        }
    };

    // Show loading state if stores are not properly initialized
    if (!isInitialized) {
        return (
            <div className="min-h-screen flex items-center justify-center">
                <div className="text-center">
                    <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-orange-600 mx-auto mb-4"></div>
                    <p className="text-gray-600 dark:text-gray-300">Initializing Annotation Boost...</p>
                </div>
            </div>
        );
    }

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
                                    <div className="w-12 h-12 bg-gradient-to-br from-orange-600 to-red-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                                        <Zap className="h-6 w-6 text-white" />
                                    </div>
                                    <div className="absolute -top-1 -right-1 w-4 h-4 bg-yellow-400 rounded-full animate-pulse"></div>
                                </div>
                                <div>
                                    <h1 className="text-3xl font-bold gradient-text">Annotation Boost</h1>
                                    <p className="text-gray-600 dark:text-gray-300">Iterative marker analysis for enhanced annotation confidence</p>
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
                <div className="glass rounded-2xl p-8 mb-8 border border-orange-400/30 bg-orange-50/20">
                    <div className="flex items-start justify-between">
                        <div className="flex items-start space-x-4">
                            <div className="w-12 h-12 bg-gradient-to-br from-orange-600 to-red-600 rounded-xl flex items-center justify-center shadow-lg">
                                <Zap className="h-6 w-6 text-white" />
                            </div>
                            <div>
                                <h3 className="text-xl font-semibold text-gray-900 dark:text-white mb-2">Annotation Boost</h3>
                                <p className="text-gray-600 dark:text-gray-300 leading-relaxed">
                                    Enhance your cell annotation confidence through iterative marker analysis. Upload conversation 
                                    history and marker genes to perform deep analysis with multiple search strategies and iterations.
                                </p>
                            </div>
                        </div>
                        <div className="flex space-x-2">
                            {!showInstructions && (
                                <Button
                                    onClick={() => setShowInstructions(true)}
                                    variant="outline"
                                    size="sm"
                                    className="glass border-blue-400/30 hover:bg-blue-50/20"
                                >
                                    <BookOpen className="h-4 w-4 mr-2" />
                                    Guide
                                </Button>
                            )}
                        </div>
                    </div>
                </div>

                {/* Instructions Panel - Now under Annotation Boost */}
                {showInstructions && (
                    <div className="glass rounded-2xl p-6 mb-8 border border-blue-400/30 bg-blue-50/20">
                        <div className="flex items-start justify-between">
                            <div className="flex items-start space-x-4">
                                <div className="w-12 h-12 bg-gradient-to-br from-blue-600 to-indigo-600 rounded-xl flex items-center justify-center shadow-lg">
                                    <BookOpen className="h-6 w-6 text-white" />
                                </div>
                                <div className="flex-1">
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-3 flex items-center">
                                        <span>Quick Start Guide</span>
                                        <span className="ml-2 px-2 py-1 bg-blue-100 dark:bg-blue-900 text-blue-700 dark:text-blue-300 text-xs rounded-full">New</span>
                                    </h3>
                                    <div className="space-y-2 text-sm text-gray-700 dark:text-gray-300">
                                        <p><strong>What is Annotation Boost?</strong> An iterative AI-powered analysis tool that examines marker genes and conversation history to provide detailed cell type annotations with high confidence.</p>
                                        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mt-4">
                                            <div className="flex items-start space-x-2">
                                                <div className="w-6 h-6 bg-green-100 dark:bg-green-900 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                                                    <span className="text-green-600 dark:text-green-400 text-xs font-bold">1</span>
                                                </div>
                                                <div>
                                                    <p className="font-medium">Load Data</p>
                                                    <p className="text-xs text-gray-600 dark:text-gray-400">Upload conversation history CSV + marker genes CSV, or try example data</p>
                                                </div>
                                            </div>
                                            <div className="flex items-start space-x-2">
                                                <div className="w-6 h-6 bg-yellow-100 dark:bg-yellow-900 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                                                    <span className="text-yellow-600 dark:text-yellow-400 text-xs font-bold">2</span>
                                                </div>
                                                <div>
                                                    <p className="font-medium">Configure</p>
                                                    <p className="text-xs text-gray-600 dark:text-gray-400">Set API key, select cluster, choose analysis parameters</p>
                                                </div>
                                            </div>
                                            <div className="flex items-start space-x-2">
                                                <div className="w-6 h-6 bg-purple-100 dark:bg-purple-900 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                                                    <span className="text-purple-600 dark:text-purple-400 text-xs font-bold">3</span>
                                                </div>
                                                <div>
                                                    <p className="font-medium">Analyze</p>
                                                    <p className="text-xs text-gray-600 dark:text-gray-400">Watch real-time progress and download detailed HTML report</p>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <Button 
                                variant="ghost" 
                                size="sm" 
                                onClick={() => setShowInstructions(false)}
                                className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
                            >
                                ×
                            </Button>
                        </div>
                    </div>
                )}

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
                                            onChange={(e) => {
                                                const newProvider = e.target.value;
                                                setProvider(newProvider);
                                                // Auto-switch to appropriate model for the provider
                                                const defaultModel = getDefaultModelForProvider(newProvider);
                                                setModel(defaultModel);
                                                localStorage.setItem('annotation-boost-model', defaultModel);
                                            }}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        >
                                            <option value="openrouter">OpenRouter</option>
                                            <option value="openai">OpenAI</option>
                                            <option value="anthropic">Anthropic</option>
                                        </select>
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Model</label>
                                        <div className="space-y-2">
                                            <select
                                                value={getModelsByProvider(provider).some(m => m.value === model) ? model : 'custom'}
                                                onChange={(e) => {
                                                    const newModel = e.target.value;
                                                    if (newModel === 'custom') {
                                                        setCustomModelId(model);
                                                    } else {
                                                        setModel(newModel);
                                                        setCustomModelId('');
                                                        // Save model preference for annotation boost
                                                        localStorage.setItem('annotation-boost-model', newModel);
                                                    }
                                                }}
                                                className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                            >
                                                {getModelsByProvider(provider).map(modelOption => (
                                                    <option key={modelOption.value} value={modelOption.value}>
                                                        {modelOption.label}
                                                    </option>
                                                ))}
                                                <option value="custom">Custom Model ID</option>
                                            </select>
                                            {!getModelsByProvider(provider).some(m => m.value === model) && (
                                                <Input
                                                    type="text"
                                                    value={model}
                                                    onChange={(e) => {
                                                        const newModel = e.target.value;
                                                        setModel(newModel);
                                                        // Save custom model preference
                                                        localStorage.setItem('annotation-boost-model', newModel);
                                                    }}
                                                    placeholder="Enter custom model ID"
                                                    className="mt-2"
                                                />
                                            )}
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Model Selection Guide - Only show for OpenRouter */}
                            {provider === 'openrouter' && (
                                <Card>
                                    <CardHeader>
                                        <CardTitle className="flex items-center gap-2">
                                            <Sparkles className="h-5 w-5" />
                                            Quick Model Selection
                                        </CardTitle>
                                    </CardHeader>
                                    <CardContent>
                                        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                                            {/* Performance Mode */}
                                            <button
                                                onClick={() => {
                                                    setSelectedMode('performance')
                                                    setModel('anthropic/claude-sonnet-4')
                                                }}
                                                className={cn(
                                                    "relative p-4 rounded-lg border-2 transition-all text-left",
                                                    selectedMode === 'performance' && model === 'anthropic/claude-sonnet-4'
                                                        ? "border-blue-500 bg-blue-50 dark:bg-blue-950/20"
                                                        : "border-gray-200 dark:border-gray-700 hover:border-gray-300 dark:hover:border-gray-600"
                                                )}
                                            >
                                                <div className="flex items-center gap-3 mb-3">
                                                    <div className="w-10 h-10 bg-gradient-to-br from-blue-500 to-purple-600 rounded-lg flex items-center justify-center">
                                                        <Zap className="h-5 w-5 text-white" />
                                                    </div>
                                                    <div>
                                                        <h3 className="font-semibold text-lg">Performance</h3>
                                                        <p className="text-sm text-muted-foreground">Best quality results</p>
                                                    </div>
                                                </div>
                                                <div className="space-y-2">
                                                    <div className="flex items-center gap-2">
                                                        <span className="text-sm font-medium">Model:</span>
                                                        <span className="text-sm">Claude Sonnet 4</span>
                                                    </div>
                                                    <div className="flex gap-1 flex-wrap">
                                                        <Badge variant="default" className="text-xs">
                                                            high quality
                                                        </Badge>
                                                        <Badge variant="secondary" className="text-xs">
                                                            medium speed
                                                        </Badge>
                                                        <Badge variant="outline" className="text-xs">
                                                            $high cost
                                                        </Badge>
                                                    </div>
                                                </div>
                                                {selectedMode === 'performance' && model === 'anthropic/claude-sonnet-4' && (
                                                    <div className="absolute top-2 right-2">
                                                        <Badge className="bg-blue-500">Selected</Badge>
                                                    </div>
                                                )}
                                            </button>

                                            {/* Balanced Mode */}
                                            <button
                                                onClick={() => {
                                                    setSelectedMode('balanced')
                                                    setModel('google/gemini-2.5-flash')
                                                }}
                                                className={cn(
                                                    "relative p-4 rounded-lg border-2 transition-all text-left",
                                                    selectedMode === 'balanced' && model === 'google/gemini-2.5-flash'
                                                        ? "border-green-500 bg-green-50 dark:bg-green-950/20"
                                                        : "border-gray-200 dark:border-gray-700 hover:border-gray-300 dark:hover:border-gray-600"
                                                )}
                                            >
                                                <div className="flex items-center gap-3 mb-3">
                                                    <div className="w-10 h-10 bg-gradient-to-br from-green-500 to-emerald-600 rounded-lg flex items-center justify-center">
                                                        <Settings className="h-5 w-5 text-white" />
                                                    </div>
                                                    <div>
                                                        <h3 className="font-semibold text-lg">Balanced</h3>
                                                        <p className="text-sm text-muted-foreground">Good quality, fast speed</p>
                                                    </div>
                                                </div>
                                                <div className="space-y-2">
                                                    <div className="flex items-center gap-2">
                                                        <span className="text-sm font-medium">Model:</span>
                                                        <span className="text-sm">Gemini 2.5 Flash</span>
                                                    </div>
                                                    <div className="flex gap-1 flex-wrap">
                                                        <Badge variant="secondary" className="text-xs">
                                                            medium quality
                                                        </Badge>
                                                        <Badge variant="default" className="text-xs">
                                                            fast speed
                                                        </Badge>
                                                        <Badge variant="default" className="text-xs">
                                                            $low cost
                                                        </Badge>
                                                    </div>
                                                </div>
                                                {selectedMode === 'balanced' && model === 'google/gemini-2.5-flash' && (
                                                    <div className="absolute top-2 right-2">
                                                        <Badge className="bg-green-500">Selected</Badge>
                                                    </div>
                                                )}
                                            </button>
                                        </div>
                                        <div className="mt-4 p-3 bg-muted/50 rounded-lg">
                                            <p className="text-sm text-muted-foreground">
                                                <strong>Tip:</strong> Use Performance mode for complex or critical analyses. 
                                                Use Balanced mode for routine processing or large datasets where speed matters.
                                            </p>
                                        </div>
                                    </CardContent>
                                </Card>
                            )}

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
                                            Iterations: {numIterations}
                                        </label>
                                        <input
                                            type="range"
                                            min="1"
                                            max="5"
                                            step="1"
                                            value={numIterations}
                                            onChange={(e) => setNumIterations(Number(e.target.value))}
                                            className="w-full"
                                        />
                                    </div>
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Search Strategy</label>
                                        <select
                                            value={searchStrategy}
                                            onChange={(e) => setSearchStrategy(e.target.value)}
                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                        >
                                            <option value="breadth">Breadth-First</option>
                                            <option value="depth">Depth-First</option>
                                        </select>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* Data Upload */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg flex items-center justify-between">
                                        <span>📁 Data Upload</span>
                                        <Button
                                            onClick={loadExampleData}
                                            variant="outline"
                                            size="sm"
                                            className="text-xs bg-green-50 border-green-200 text-green-700 hover:bg-green-100 dark:bg-green-900/20 dark:border-green-700 dark:text-green-300"
                                        >
                                            <Database className="h-3 w-3 mr-1" />
                                            Load Example
                                        </Button>
                                    </CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-4">
                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Conversation History CSV</label>
                                        <div 
                                            {...getConversationRootProps()} 
                                            className={`glass border-2 border-dashed rounded-xl p-6 text-center cursor-pointer transition-all duration-300 ${
                                                isConversationDragActive 
                                                    ? 'border-orange-400 bg-orange-50/20' 
                                                    : 'border-white/30 hover:border-orange-300 hover:bg-white/10'
                                            }`}
                                        >
                                            <input {...getConversationInputProps()} />
                                            <History className="h-8 w-8 text-gray-400 mx-auto mb-2" />
                                            {conversationFile ? (
                                                <div>
                                                    <p className="text-sm font-medium text-gray-900 dark:text-white">{conversationFile.name}</p>
                                                    <p className="text-xs text-gray-500 dark:text-gray-400">
                                                        {availableColumns.length} columns detected
                                                    </p>
                                                </div>
                                            ) : (
                                                <div>
                                                    <p className="text-sm text-gray-600 dark:text-gray-300">
                                                        {isConversationDragActive ? "Drop conversation CSV..." : "Upload conversation history"}
                                                    </p>
                                                    <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                                        CSV files only
                                                    </p>
                                                </div>
                                            )}
                                        </div>
                                    </div>

                                    <div>
                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Marker Genes CSV</label>
                                        <div 
                                            {...getMarkerRootProps()} 
                                            className={`glass border-2 border-dashed rounded-xl p-6 text-center cursor-pointer transition-all duration-300 ${
                                                isMarkerDragActive 
                                                    ? 'border-orange-400 bg-orange-50/20' 
                                                    : 'border-white/30 hover:border-orange-300 hover:bg-white/10'
                                            }`}
                                        >
                                            <input {...getMarkerInputProps()} />
                                            <FileText className="h-8 w-8 text-gray-400 mx-auto mb-2" />
                                            {markerFile ? (
                                                <div>
                                                    <p className="text-sm font-medium text-gray-900 dark:text-white">{markerFile.name}</p>
                                                    <p className="text-xs text-gray-500 dark:text-gray-400">
                                                        {markerData.length} rows loaded
                                                    </p>
                                                </div>
                                            ) : (
                                                <div>
                                                    <p className="text-sm text-gray-600 dark:text-gray-300">
                                                        {isMarkerDragActive ? "Drop marker CSV..." : "Upload marker genes"}
                                                    </p>
                                                    <p className="text-xs text-gray-500 dark:text-gray-400 mt-1">
                                                        CSV files only
                                                    </p>
                                                </div>
                                            )}
                                        </div>
                                    </div>
                                </CardContent>
                            </Card>

                            {/* History Configuration */}
                            {conversationFile && (
                                <Card>
                                    <CardHeader>
                                        <CardTitle className="text-lg">📜 History Configuration</CardTitle>
                                    </CardHeader>
                                    <CardContent className="space-y-4">
                                        <div>
                                            <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">History Source</label>
                                            <select
                                                value={historySource}
                                                onChange={(e) => setHistorySource(e.target.value as 'csv' | 'manual')}
                                                className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                            >
                                                <option value="csv">Extract from CSV</option>
                                                <option value="manual">Manual Input</option>
                                            </select>
                                        </div>

                                        {historySource === 'csv' && (
                                            <>
                                                <div>
                                                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Cluster Column</label>
                                                    <select
                                                        value={selectedClusterColumn}
                                                        onChange={(e) => setSelectedClusterColumn(e.target.value)}
                                                        className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                                    >
                                                        {availableColumns.map(col => (
                                                            <option key={col} value={col}>{col}</option>
                                                        ))}
                                                    </select>
                                                </div>
                                                
                                                <Button
                                                    onClick={loadClusters}
                                                    variant="outline"
                                                    size="sm"
                                                    className="w-full"
                                                >
                                                    <RefreshCcw className="h-4 w-4 mr-2" />
                                                    Load Clusters
                                                </Button>

                                                {availableClusters.length > 0 && (
                                                    <div>
                                                        <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Select Cluster</label>
                                                        <select
                                                            value={selectedCluster}
                                                            onChange={(e) => setSelectedCluster(e.target.value)}
                                                            className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20"
                                                        >
                                                            <option value="">Choose cluster...</option>
                                                            {availableClusters.map(cluster => (
                                                                <option key={cluster} value={cluster}>{cluster}</option>
                                                            ))}
                                                        </select>
                                                    </div>
                                                )}

                                                {selectedCluster && (
                                                    <Button
                                                        onClick={extractConversation}
                                                        variant="outline"
                                                        size="sm"
                                                        className="w-full"
                                                    >
                                                        <Upload className="h-4 w-4 mr-2" />
                                                        Extract Conversation
                                                    </Button>
                                                )}
                                            </>
                                        )}

                                        {historySource === 'manual' && (
                                            <div>
                                                <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Conversation History</label>
                                                <textarea
                                                    value={manualHistory}
                                                    onChange={(e) => setManualHistory(e.target.value)}
                                                    placeholder="Enter conversation history manually..."
                                                    className="w-full px-3 py-2 glass border border-white/30 rounded-xl form-modern text-gray-900 dark:text-white bg-white/20 dark:bg-black/20 min-h-[120px] resize-y"
                                                />
                                            </div>
                                        )}

                                        {extractedHistory && historySource === 'csv' && (
                                            <div className="glass rounded-lg p-3 border border-white/20">
                                                <p className="text-sm font-medium text-gray-900 dark:text-white mb-2">Extracted History Preview:</p>
                                                <p className="text-xs text-gray-600 dark:text-gray-300 max-h-24 overflow-y-auto">
                                                    {extractedHistory.substring(0, 200)}...
                                                </p>
                                            </div>
                                        )}
                                    </CardContent>
                                </Card>
                            )}

                            {/* Run Analysis Button */}
                            <Button
                                onClick={handleRunAnalysis}
                                disabled={isLoading || !markerData.length || !apiKey || !(extractedHistory || manualHistory)}
                                className="w-full bg-gradient-to-r from-orange-600 to-red-600 hover:from-orange-700 hover:to-red-700 text-white btn-modern"
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
                                        Start Annotation Boost
                                    </>
                                )}
                            </Button>
                        </div>
                    </div>

                    {/* Right Panel - Results & Console */}
                    <div className="xl:col-span-3 space-y-6">
                        {/* Console Output */}
                        <div className="glass rounded-2xl p-6 border border-white/20">
                            <div className="flex items-center justify-between mb-4">
                                <h2 className="text-lg font-bold gradient-text flex items-center">
                                    🔍 <span className="ml-2">Analysis Progress</span>
                                </h2>
                                <Button
                                    onClick={clearConsole}
                                    variant="outline"
                                    size="sm"
                                    className="text-xs"
                                    disabled={consoleOutput.length === 0}
                                >
                                    Clear
                                </Button>
                            </div>
                            <div 
                                className="bg-gray-900 dark:bg-black rounded-lg p-4 min-h-[200px] max-h-[300px] overflow-y-auto font-mono text-sm"
                                ref={(el) => {
                                    if (el && consoleOutput.length > 0) {
                                        el.scrollTop = el.scrollHeight;
                                    }
                                }}
                            >
                                {consoleOutput.length === 0 ? (
                                    <div className="text-gray-500 dark:text-gray-400 flex items-center">
                                        <AlertCircle className="h-4 w-4 mr-2" />
                                        Console output will appear here during analysis...
                                    </div>
                                ) : (
                                    <div className="space-y-1">
                                        {consoleOutput.map((line, index) => (
                                            <div key={index} className="text-green-400 dark:text-green-300">
                                                {line}
                                            </div>
                                        ))}
                                        {isLoading && (
                                            <div className="flex items-center text-yellow-400 dark:text-yellow-300 animate-pulse">
                                                <Loader2 className="h-3 w-3 mr-2 animate-spin" />
                                                Analysis in progress...
                                            </div>
                                        )}
                                    </div>
                                )}
                            </div>
                        </div>

                        {/* Results */}
                        <div className="glass rounded-2xl p-6 border border-white/20">
                            <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                                ⚡ <span className="ml-2">Boost Results</span>
                            </h2>

                            {isLoading && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 bg-gradient-to-r from-orange-600 to-red-600 rounded-full flex items-center justify-center mx-auto mb-4 animate-pulse">
                                        <Zap className="h-8 w-8 text-white" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Running Annotation Boost
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300 mb-4">
                                        Performing iterative marker analysis...
                                    </p>
                                    <div className="glass rounded-lg p-4 border border-white/20">
                                        <div className="flex justify-center space-x-2">
                                            <div className="w-2 h-2 bg-orange-600 rounded-full animate-bounce"></div>
                                            <div className="w-2 h-2 bg-red-600 rounded-full animate-bounce" style={{animationDelay: '0.1s'}}></div>
                                            <div className="w-2 h-2 bg-yellow-600 rounded-full animate-bounce" style={{animationDelay: '0.2s'}}></div>
                                        </div>
                                        <p className="text-sm text-gray-600 dark:text-gray-300 mt-2">
                                            Iteration {Math.min(3, numIterations)} of {numIterations}...
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

                            {resultsHtml && (
                                <div className="space-y-6">
                                    <div className="glass rounded-lg p-4 border border-green-400/30 bg-green-50/20">
                                        <h3 className="text-lg font-semibold text-green-800 dark:text-green-400 mb-2">
                                            Annotation Boost Complete! ⚡
                                        </h3>
                                        <p className="text-green-700 dark:text-green-300">
                                            Iterative analysis completed with {numIterations} iterations using {searchStrategy}-first strategy.
                                        </p>
                                    </div>

                                    <Card>
                                        <CardHeader>
                                            <CardTitle className="flex items-center justify-between">
                                                <span>🔬 Analysis Results</span>
                                                <Button
                                                    onClick={() => {
                                                        const blob = new Blob([resultsHtml], { type: 'text/html' });
                                                        const url = URL.createObjectURL(blob);
                                                        const a = document.createElement('a');
                                                        a.href = url;
                                                        a.download = 'annotation_boost_results.html';
                                                        document.body.appendChild(a);
                                                        a.click();
                                                        document.body.removeChild(a);
                                                        URL.revokeObjectURL(url);
                                                    }}
                                                    variant="outline"
                                                    size="sm"
                                                >
                                                    <Download className="h-4 w-4 mr-2" />
                                                    Download Report
                                                </Button>
                                            </CardTitle>
                                        </CardHeader>
                                        <CardContent>
                                            <div 
                                                className="glass rounded-lg p-4 border border-white/20 min-h-[500px] max-h-[600px] overflow-y-auto prose prose-sm dark:prose-invert max-w-none"
                                                dangerouslySetInnerHTML={{ __html: resultsHtml }}
                                            />
                                        </CardContent>
                                    </Card>
                                </div>
                            )}

                            {!isLoading && !resultsHtml && !error && (
                                <div className="text-center py-12">
                                    <div className="w-16 h-16 glass rounded-full flex items-center justify-center mx-auto mb-4 border border-white/20">
                                        <Zap className="h-8 w-8 text-gray-400" />
                                    </div>
                                    <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                                        Ready for Annotation Boost
                                    </h3>
                                    <p className="text-gray-600 dark:text-gray-300">
                                        Upload your conversation history and marker genes to start iterative annotation analysis.
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