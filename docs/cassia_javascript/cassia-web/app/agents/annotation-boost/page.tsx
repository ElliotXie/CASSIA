"use client";

import { useState, useEffect } from 'react';
import Link from 'next/link';
import { useDropzone } from 'react-dropzone';
import Papa from 'papaparse';
import { iterativeMarkerAnalysis, generateSummaryReport, extractConversationForCluster, getAvailableClusters, extractTopMarkerGenes } from '@/lib/cassia/annotationBoost';
import { useApiKeyStore, Provider } from '@/lib/stores/api-key-store';
import { useAuthStore } from '@/lib/stores/auth-store';
import { ReasoningEffort } from '@/lib/config/model-presets';
import { useConfigStore } from '@/lib/stores/config-store';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { ArrowLeft, Play, HelpCircle, Zap, Upload, Download, FileText, History, BookOpen, Database, AlertCircle, Loader2, CheckCircle } from 'lucide-react';
import { AgentModelSelector } from '@/components/AgentModelSelector';
import { testApiKey } from '@/lib/cassia/llm_utils';

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
            return 'google/gemini-2.5-flash';
        }
    });
    const loadApiKeys = useApiKeyStore((state) => state.loadApiKeys);

    // Auth state
    const { isAuthenticated, user } = useAuthStore();

    // State management
    const [conversationFile, setConversationFile] = useState<File | null>(null);
    const [conversationData, setConversationData] = useState<string>('');
    const [markerFile, setMarkerFile] = useState<File | null>(null);
    const [markerData, setMarkerData] = useState<any[]>([]);
    const [apiKey, setApiKey] = useState('');
    const [provider, setProvider] = useState<Provider>('openrouter');
    const [model, setModel] = useState('google/gemini-2.5-flash');
    const [customBaseUrl, setCustomBaseUrl] = useState('');
    const [reasoningEffort, setReasoningEffort] = useState<ReasoningEffort | null>(null);
    const [majorClusterInfo, setMajorClusterInfo] = useState('Human PBMC');
    const [numIterations, setNumIterations] = useState(3);
    const [searchStrategy, setSearchStrategy] = useState('breadth');
    
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
    const [isTestingApi, setIsTestingApi] = useState(false);
    const [testStatus, setTestStatus] = useState<'idle' | 'success' | 'error'>('idle');
    const [testErrorMessage, setTestErrorMessage] = useState<string>('');

    // Load API key from account state
    const [isLoadingKeys, setIsLoadingKeys] = useState(false);
    const [loadKeyStatus, setLoadKeyStatus] = useState<'idle' | 'success' | 'error'>('idle');
    const [loadKeyError, setLoadKeyError] = useState('');

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
                    setProvider(globalProvider as Provider);
                }
                // For annotation boost, we prefer Gemini 2.5 Flash as default
                const savedModel = localStorage.getItem('annotation-boost-model');
                if (savedModel && savedModel !== model) {
                    setModel(savedModel);
                }
                setIsInitialized(true);
            }
        } catch (error) {
            console.error('Error initializing global settings:', error);
            setIsInitialized(true);
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

    // Auto-load clusters when conversation data or cluster column changes
    useEffect(() => {
        if (conversationData && selectedClusterColumn) {
            loadClusters();
        }
    }, [conversationData, selectedClusterColumn]);

    // Auto-extract conversation when cluster is selected (no manual button needed)
    useEffect(() => {
        if (conversationData && selectedCluster && selectedClusterColumn) {
            extractConversation();
        }
    }, [selectedCluster, selectedClusterColumn, conversationData]);

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
"plasma cell","Plasma Cell","Plasmablast, Circulating Plasma Cell, Short-lived Plasma Cell","","50","IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3","1","anthropic/claude-3.5-sonnet","openrouter","large_intestine","human","Example analysis","Example conversation history for plasma cell analysis demonstrating iterative marker analysis workflow."`;
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

    // Load API keys from Supabase account
    const handleLoadApiKeys = async () => {
        if (!isAuthenticated || !user) {
            setLoadKeyError('Please sign in to load API keys');
            setLoadKeyStatus('error');
            setTimeout(() => setLoadKeyStatus('idle'), 5000);
            return;
        }

        setIsLoadingKeys(true);
        setLoadKeyStatus('idle');
        setLoadKeyError('');

        try {
            await loadApiKeys();
            // Get the loaded key for current provider
            const loadedKey = useApiKeyStore.getState().apiKeys[provider];
            if (loadedKey) {
                setApiKey(loadedKey);
                setLoadKeyStatus('success');
                setTimeout(() => setLoadKeyStatus('idle'), 3000);
            } else {
                setLoadKeyError('No API key found for this provider');
                setLoadKeyStatus('error');
                setTimeout(() => setLoadKeyStatus('idle'), 5000);
            }
        } catch (err: any) {
            setLoadKeyError(err.message || 'Failed to load API keys');
            setLoadKeyStatus('error');
            setTimeout(() => setLoadKeyStatus('idle'), 5000);
        } finally {
            setIsLoadingKeys(false);
        }
    };

    // Test API key
    const handleTestApiKey = async () => {
        if (!apiKey.trim()) {
            setTestStatus('error');
            setTestErrorMessage('Please enter an API key first');
            setTimeout(() => setTestStatus('idle'), 3000);
            return;
        }

        setIsTestingApi(true);
        setTestStatus('idle');
        setTestErrorMessage('');

        try {
            const baseUrl = provider === 'custom' ? customBaseUrl : null;
            const result = await testApiKey(provider, apiKey, baseUrl);

            if (result.success) {
                setTestStatus('success');
                setTimeout(() => setTestStatus('idle'), 5000);
            } else {
                setTestStatus('error');
                setTestErrorMessage(result.error || 'API test failed');
                setTimeout(() => setTestStatus('idle'), 5000);
            }
        } catch (err: any) {
            setTestStatus('error');
            setTestErrorMessage(err.message || 'API test failed');
            setTimeout(() => setTestStatus('idle'), 5000);
        } finally {
            setIsTestingApi(false);
        }
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
            // Use customBaseUrl as provider for custom endpoints
            const effectiveProvider = provider === 'custom' ? customBaseUrl : provider;

            const results = await iterativeMarkerAnalysis(
                majorClusterInfo,
                markerData,
                extractTopMarkerGenes(markerData, 20), // Extract top genes as comma-separated string
                historyToUse,
                numIterations,
                effectiveProvider,
                model,
                null, // additionalTask
                0, // temperature
                searchStrategy as 'breadth' | 'depth',
                apiKey,
                reasoningEffort
            );

            // Skip the first message which contains the prompt (matching Python implementation)
            const conversationWithoutPrompt = results.messages.length > 1 ? results.messages.slice(1) : results.messages;

            const htmlReport = await generateSummaryReport(
                conversationWithoutPrompt,
                searchStrategy as 'breadth' | 'depth',
                'per_iteration',
                effectiveProvider,
                model,
                apiKey,
                reasoningEffort
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

                            {/* Model Selection */}
                            <AgentModelSelector
                                provider={provider}
                                model={model}
                                onProviderChange={setProvider}
                                onModelChange={setModel}
                                customBaseUrl={customBaseUrl}
                                onCustomBaseUrlChange={setCustomBaseUrl}
                                reasoningEffort={reasoningEffort}
                                onReasoningEffortChange={setReasoningEffort}
                            />

                            {/* API Key */}
                            <Card>
                                <CardHeader>
                                    <CardTitle className="text-lg">🔑 API Key</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-3">
                                    <Input
                                        type="password"
                                        value={apiKey}
                                        onChange={(e) => setApiKey(e.target.value)}
                                        placeholder="Enter your API key"
                                    />
                                    {isAuthenticated && (
                                        <Button
                                            onClick={handleLoadApiKeys}
                                            disabled={isLoadingKeys}
                                            variant={loadKeyStatus === 'success' ? 'default' : loadKeyStatus === 'error' ? 'destructive' : 'outline'}
                                            size="sm"
                                            className="w-full"
                                        >
                                            {isLoadingKeys ? (
                                                <>
                                                    <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                                                    Loading Keys...
                                                </>
                                            ) : loadKeyStatus === 'success' ? (
                                                <>
                                                    <CheckCircle className="mr-2 h-4 w-4" />
                                                    API Key Loaded
                                                </>
                                            ) : loadKeyStatus === 'error' ? (
                                                <>
                                                    <AlertCircle className="mr-2 h-4 w-4" />
                                                    Load Failed
                                                </>
                                            ) : (
                                                <>
                                                    <Download className="mr-2 h-4 w-4" />
                                                    Load from Account
                                                </>
                                            )}
                                        </Button>
                                    )}
                                    {loadKeyStatus === 'error' && loadKeyError && (
                                        <p className="text-xs text-red-600">{loadKeyError}</p>
                                    )}
                                    {loadKeyStatus === 'success' && (
                                        <p className="text-xs text-green-600">API key loaded from your account</p>
                                    )}
                                    <Button
                                        onClick={handleTestApiKey}
                                        disabled={isTestingApi || !apiKey}
                                        variant={testStatus === 'success' ? 'default' : testStatus === 'error' ? 'destructive' : 'outline'}
                                        size="sm"
                                        className="w-full"
                                    >
                                        {isTestingApi ? (
                                            <>
                                                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                                                Testing API...
                                            </>
                                        ) : testStatus === 'success' ? (
                                            <>
                                                <CheckCircle className="mr-2 h-4 w-4" />
                                                API Key Valid
                                            </>
                                        ) : testStatus === 'error' ? (
                                            <>
                                                <AlertCircle className="mr-2 h-4 w-4" />
                                                Test Failed
                                            </>
                                        ) : (
                                            <>
                                                <Zap className="mr-2 h-4 w-4" />
                                                Test API Key
                                            </>
                                        )}
                                    </Button>
                                    {testStatus === 'error' && testErrorMessage && (
                                        <p className="text-xs text-red-600">{testErrorMessage}</p>
                                    )}
                                    {testStatus === 'success' && (
                                        <p className="text-xs text-green-600">API key is valid and working</p>
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

                                                {/* Clusters load automatically when conversation file is uploaded */}
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

                                                {/* Conversation is extracted automatically when cluster is selected */}
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
                                            <div className="glass rounded-lg p-3 border border-green-500/30 bg-green-500/5">
                                                <p className="text-sm font-medium text-green-700 dark:text-green-400 mb-2 flex items-center">
                                                    <CheckCircle className="h-4 w-4 mr-2" />
                                                    Conversation Extracted
                                                </p>
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
                                            <iframe
                                                srcDoc={resultsHtml}
                                                className="w-full min-h-[500px] h-[600px] rounded-lg border border-white/20"
                                                title="Annotation Boost Results"
                                                sandbox="allow-same-origin"
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