"use client";

import { useState, ChangeEvent, useEffect } from 'react';
import Link from 'next/link';
import { useDropzone } from 'react-dropzone';
import Papa from 'papaparse';
import { runCASSIASubclusters, runCASSIANSubcluster } from '@/lib/cassia/subclustering';
import { useApiKeyStore } from '@/lib/stores/api-key-store';
import { Button } from '@/components/ui/button';
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { ArrowLeft, Play, HelpCircle, Layers, Upload, Download, FileText } from 'lucide-react';
import { ContactDialog } from '@/components/ContactDialog';
import modelSettings from '../public/examples/model_settings.json';

export default function SubclusteringPage() {
  const [showContactModal, setShowContactModal] = useState(false);
  const globalApiKey = useApiKeyStore((state) => state.getApiKey());
  const globalProvider = useApiKeyStore((state) => state.provider);
  
  // State management
  const [markerFile, setMarkerFile] = useState<File | null>(null);
  const [markerData, setMarkerData] = useState<any[]>([]);
  const [majorClusterInfo, setMajorClusterInfo] = useState('CD8 T cells');
  const [outputName, setOutputName] = useState('subcluster_results');
  const [nGenes, setNGenes] = useState(50);
  const [model, setModel] = useState('google/gemini-2.5-flash');
  const [temperature, setTemperature] = useState(0);
  const [provider, setProvider] = useState('openrouter');
  const [apiKey, setApiKey] = useState('');
  
  // Batch settings
  const [numRuns, setNumRuns] = useState(3);
  const [maxWorkers, setMaxWorkers] = useState(5);

  const [isLoading, setIsLoading] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [batchResults, setBatchResults] = useState<any[]>([]);
  const [error, setError] = useState<string | null>(null);

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

  const onDrop = (acceptedFiles: File[]) => {
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
          setError(`Error parsing CSV: ${err.message}`);
        }
      });
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ 
    onDrop, 
    accept: { 'text/csv': ['.csv'] },
    multiple: false
  });

  const handleRunSingle = async () => {
    if (!markerData.length || !apiKey) {
      setError('Please upload a marker file and provide an API key.');
      return;
    }
    setIsLoading(true);
    setError(null);
    setResults(null);
    setBatchResults([]);

    try {
      const result = await runCASSIASubclusters(
        markerData,
        majorClusterInfo,
        apiKey,
        model,
        Number(temperature),
        provider,
        Number(nGenes)
      );
      setResults(result);
    } catch (err: any) {
      setError(`Analysis failed: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };

  const handleRunBatch = async () => {
    if (!markerData.length || !apiKey) {
      setError('Please upload a marker file and provide an API key.');
      return;
    }
    setIsLoading(true);
    setError(null);
    setResults(null);
    setBatchResults([]);

    try {
      const batchResult = await runCASSIANSubcluster(
        Number(numRuns),
        markerData,
        majorClusterInfo,
        apiKey,
        model,
        Number(temperature),
        provider,
        Number(maxWorkers),
        Number(nGenes)
      );
      setBatchResults(batchResult);
    } catch (err: any) {
      setError(`Batch analysis failed: ${err.message}`);
    } finally {
      setIsLoading(false);
    }
  };
  
  const downloadCSV = (csvContent: string, fileName: string) => {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', fileName);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
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
                  <div className="w-12 h-12 bg-gradient-to-br from-teal-600 to-green-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                    <Layers className="h-6 w-6 text-white" />
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-orange-400 rounded-full animate-pulse"></div>
                </div>
                <div>
                  <h1 className="text-3xl font-bold gradient-text">CASSIA Subclustering</h1>
                  <p className="text-gray-600 dark:text-gray-300">Subtype analysis for detailed cell classification</p>
                </div>
              </div>
            </div>
            <Button 
              variant="outline" 
              size="sm" 
              onClick={() => setShowContactModal(true)}
              className="glass border-white/30 hover:bg-white/20 btn-modern"
            >
              <HelpCircle className="h-4 w-4 mr-2" />
              Help
            </Button>
          </div>
        </div>
      </header>

      <main className="container mx-auto px-6 py-12">
        {/* Info Panel */}
        <div className="glass rounded-2xl p-8 mb-12 border border-teal-400/30 bg-teal-50/20">
          <div className="flex items-start space-x-4">
            <div className="w-12 h-12 bg-gradient-to-br from-teal-600 to-green-600 rounded-xl flex items-center justify-center shadow-lg">
              <Layers className="h-6 w-6 text-white" />
            </div>
            <div>
              <h3 className="text-xl font-semibold text-gray-900 dark:text-white mb-2">Subclustering Analysis</h3>
              <p className="text-gray-600 dark:text-gray-300 leading-relaxed">
                Perform detailed subtype analysis within major cell clusters. Upload your marker gene data
                to identify subclusters and generate comprehensive annotation reports with confidence scores.
              </p>
            </div>
          </div>
        </div>

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
                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Major Cluster Info</label>
                    <Input
                      type="text"
                      value={majorClusterInfo}
                      onChange={(e) => setMajorClusterInfo(e.target.value)}
                      placeholder="e.g., CD8 T cells, B cells"
                    />
                  </div>
                  <div>
                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Output Name</label>
                    <Input
                      type="text"
                      value={outputName}
                      onChange={(e) => setOutputName(e.target.value)}
                      placeholder="Output file name base"
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
                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">
                      Number of Genes: {nGenes}
                    </label>
                    <input
                      type="range"
                      min="10"
                      max="200"
                      step="10"
                      value={nGenes}
                      onChange={(e) => setNGenes(Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
                </CardContent>
              </Card>

              {/* Batch Settings */}
              <Card>
                <CardHeader>
                  <CardTitle className="text-lg">🔄 Batch Settings</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">
                      Number of Runs: {numRuns}
                    </label>
                    <input
                      type="range"
                      min="1"
                      max="10"
                      step="1"
                      value={numRuns}
                      onChange={(e) => setNumRuns(Number(e.target.value))}
                      className="w-full"
                    />
                  </div>
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
                </CardContent>
              </Card>

              {/* File Upload */}
              <Card>
                <CardHeader>
                  <CardTitle className="text-lg">📁 Data Upload</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                  <div>
                    <label className="block text-sm font-medium text-gray-900 dark:text-white mb-2">Marker CSV File</label>
                    <div 
                      {...getRootProps()} 
                      className={`glass border-2 border-dashed rounded-xl p-6 text-center cursor-pointer transition-all duration-300 ${
                        isDragActive 
                          ? 'border-teal-400 bg-teal-50/20' 
                          : 'border-white/30 hover:border-teal-300 hover:bg-white/10'
                      }`}
                    >
                      <input {...getInputProps()} />
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
                            {isDragActive ? "Drop the file here..." : "Drag & drop CSV file or click to select"}
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

              {/* Action Buttons */}
              <div className="space-y-3">
                <Button
                  onClick={handleRunSingle}
                  disabled={isLoading || !markerData.length || !apiKey}
                  className="w-full bg-gradient-to-r from-teal-600 to-green-600 hover:from-teal-700 hover:to-green-700 text-white btn-modern"
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
                      Run Single Analysis
                    </>
                  )}
                </Button>

                <Button
                  onClick={handleRunBatch}
                  disabled={isLoading || !markerData.length || !apiKey}
                  variant="outline"
                  className="w-full btn-modern"
                  size="lg"
                >
                  {isLoading ? (
                    <>
                      <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-current mr-2"></div>
                      Running Batch...
                    </>
                  ) : (
                    <>
                      <Layers className="h-4 w-4 mr-2" />
                      Run Batch Analysis
                    </>
                  )}
                </Button>
              </div>
            </div>
          </div>

          {/* Right Panel - Results */}
          <div className="xl:col-span-3">
            <div className="glass rounded-2xl p-6 border border-white/20">
              <h2 className="text-xl font-bold gradient-text mb-6 flex items-center">
                🔬 <span className="ml-2">Analysis Results</span>
              </h2>

              {isLoading && (
                <div className="text-center py-12">
                  <div className="w-16 h-16 bg-gradient-to-r from-teal-600 to-green-600 rounded-full flex items-center justify-center mx-auto mb-4 animate-pulse">
                    <Layers className="h-8 w-8 text-white" />
                  </div>
                  <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                    Running Subclustering Analysis
                  </h3>
                  <p className="text-gray-600 dark:text-gray-300 mb-4">
                    Analyzing marker genes for subtype classification...
                  </p>
                  <div className="glass rounded-lg p-4 border border-white/20">
                    <div className="flex justify-center space-x-2">
                      <div className="w-2 h-2 bg-teal-600 rounded-full animate-bounce"></div>
                      <div className="w-2 h-2 bg-green-600 rounded-full animate-bounce" style={{animationDelay: '0.1s'}}></div>
                      <div className="w-2 h-2 bg-emerald-600 rounded-full animate-bounce" style={{animationDelay: '0.2s'}}></div>
                    </div>
                    <p className="text-sm text-gray-600 dark:text-gray-300 mt-2">
                      Processing data...
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
                      Single Analysis Complete! 🎯
                    </h3>
                    <p className="text-green-700 dark:text-green-300">
                      Subclustering analysis has been completed successfully.
                    </p>
                  </div>

                  <Card>
                    <CardHeader>
                      <CardTitle>📊 Single Analysis Results</CardTitle>
                    </CardHeader>
                    <CardContent>
                      <Button
                        onClick={() => downloadCSV(results.csvContent, `${outputName}.csv`)}
                        className="w-full mb-4"
                        variant="outline"
                      >
                        <Download className="h-4 w-4 mr-2" />
                        Download Results CSV
                      </Button>
                      
                      {results.dataFrame && (
                        <div className="glass rounded-lg p-4 border border-white/20 overflow-auto max-h-96">
                          <table className="w-full text-sm">
                            <thead>
                              <tr className="border-b border-white/20">
                                {results.dataFrame.columns?.map((header: string) => (
                                  <th key={header} className="text-left p-2 font-medium">
                                    {header}
                                  </th>
                                ))}
                              </tr>
                            </thead>
                            <tbody>
                              {results.dataFrame.columns?.[0] && results.dataFrame[results.dataFrame.columns[0]]?.map((_: any, i: number) => (
                                <tr key={i} className="border-b border-white/10">
                                  {results.dataFrame.columns.map((col: string) => (
                                    <td key={`${col}-${i}`} className="p-2">
                                      {results.dataFrame[col][i]}
                                    </td>
                                  ))}
                                </tr>
                              ))}
                            </tbody>
                          </table>
                        </div>
                      )}
                    </CardContent>
                  </Card>
                </div>
              )}

              {batchResults.length > 0 && (
                <div className="space-y-6">
                  <div className="glass rounded-lg p-4 border border-green-400/30 bg-green-50/20">
                    <h3 className="text-lg font-semibold text-green-800 dark:text-green-400 mb-2">
                      Batch Analysis Complete! 🚀
                    </h3>
                    <p className="text-green-700 dark:text-green-300">
                      All {batchResults.length} batch runs have been completed.
                    </p>
                  </div>

                  <Card>
                    <CardHeader>
                      <CardTitle>📊 Batch Analysis Results</CardTitle>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-4">
                        {batchResults.map((result, index) => (
                          <div key={index} className="glass rounded-lg p-4 border border-white/20">
                            <div className="flex items-center justify-between mb-3">
                              <h4 className="font-semibold text-gray-900 dark:text-white">
                                Run {result.iteration}
                              </h4>
                              {result.error ? (
                                <span className="text-red-500 text-sm">Error</span>
                              ) : (
                                <span className="text-green-500 text-sm">Success</span>
                              )}
                            </div>
                            
                            {result.error ? (
                              <p className="text-red-600 dark:text-red-400 text-sm">
                                Error: {result.error}
                              </p>
                            ) : (
                              <div className="space-y-3">
                                <Button
                                  onClick={() => downloadCSV(result.csvContent, `${outputName}_${result.iteration}.csv`)}
                                  variant="outline"
                                  size="sm"
                                  className="w-full"
                                >
                                  <Download className="h-4 w-4 mr-2" />
                                  Download Run {result.iteration} CSV
                                </Button>
                                
                                {result.dataFrame && (
                                  <div className="glass rounded-lg p-3 border border-white/20 overflow-auto max-h-48">
                                    <table className="w-full text-xs">
                                      <thead>
                                        <tr className="border-b border-white/20">
                                          {Object.keys(result.dataFrame[0] || {}).map((header) => (
                                            <th key={header} className="text-left p-1 font-medium">
                                              {header}
                                            </th>
                                          ))}
                                        </tr>
                                      </thead>
                                      <tbody>
                                        {result.dataFrame.slice(0, 5).map((row: any, i: number) => (
                                          <tr key={i} className="border-b border-white/10">
                                            {Object.keys(result.dataFrame[0] || {}).map((col) => (
                                              <td key={`${col}-${i}`} className="p-1">
                                                {row[col]}
                                              </td>
                                            ))}
                                          </tr>
                                        ))}
                                      </tbody>
                                    </table>
                                    {result.dataFrame.length > 5 && (
                                      <p className="text-xs text-gray-500 mt-2 text-center">
                                        Showing 5 of {result.dataFrame.length} rows
                                      </p>
                                    )}
                                  </div>
                                )}
                              </div>
                            )}
                          </div>
                        ))}
                      </div>
                    </CardContent>
                  </Card>
                </div>
              )}

              {!isLoading && !results && !error && batchResults.length === 0 && (
                <div className="text-center py-12">
                  <div className="w-16 h-16 glass rounded-full flex items-center justify-center mx-auto mb-4 border border-white/20">
                    <Layers className="h-8 w-8 text-gray-400" />
                  </div>
                  <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                    Ready for Subclustering
                  </h3>
                  <p className="text-gray-600 dark:text-gray-300">
                    Upload your marker gene data and configure settings to start subclustering analysis.
                  </p>
                </div>
              )}
            </div>
          </div>
        </div>
      </main>

      {/* Contact Dialog Modal */}
      {showContactModal && (
        <div className="fixed inset-0 bg-black/80 backdrop-blur-md flex items-center justify-center z-50 p-4">
          <div className="bg-white dark:bg-gray-900 rounded-2xl p-8 max-w-lg w-full border-2 border-blue-200 dark:border-blue-800 shadow-2xl">
            <ContactDialog />
            <div className="flex justify-end pt-4">
              <Button 
                variant="outline" 
                onClick={() => setShowContactModal(false)}
                className="border-2 border-gray-300 dark:border-gray-600 hover:bg-gray-50 dark:hover:bg-gray-800 font-semibold py-3"
              >
                Close
              </Button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}