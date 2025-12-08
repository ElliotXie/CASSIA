'use client';

import { useState, useEffect } from 'react';
import Link from 'next/link';
import { ProgressTracker } from '@/components/ProgressTracker';
import { Button } from '@/components/ui/button';
import { Card } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'; // Still used for detail-level
import { Textarea } from '@/components/ui/textarea';
import { Label } from '@/components/ui/label';
import { useApiKeyStore, Provider } from '@/lib/stores/api-key-store';
import { useAuthStore } from '@/lib/stores/auth-store';
import { ReasoningEffort } from '@/lib/config/model-presets';
import { mergeAnnotations, mergeAnnotationsAll } from '@/lib/cassia/mergingAnnotation';
import { parseCSV } from '@/lib/utils/csv-parser';
import { Upload, File, CheckCircle, AlertCircle, X, Zap, Loader2, Download, ArrowLeft } from 'lucide-react';
import { useCallback } from 'react';
import { useDropzone } from 'react-dropzone';
import { AgentModelSelector } from '@/components/AgentModelSelector';
import { CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { testApiKey } from '@/lib/cassia/llm_utils';

export default function AnnotationMergingPage() {
  const [csvData, setCsvData] = useState(null);
  const [isProcessing, setIsProcessing] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState('');
  const [progress, setProgress] = useState('');

  // Form state
  const [detailLevel, setDetailLevel] = useState('broad');
  const [processAllLevels, setProcessAllLevels] = useState(false);
  const [batchSize, setBatchSize] = useState(20);
  const [additionalContext, setAdditionalContext] = useState('');
  const [provider, setProvider] = useState<Provider>('openrouter');
  const [model, setModel] = useState('google/gemini-2.5-flash');
  const [customBaseUrl, setCustomBaseUrl] = useState('');
  const [reasoningEffort, setReasoningEffort] = useState<ReasoningEffort | null>(null);
  const [apiKey, setApiKey] = useState('');
  const [isTestingApi, setIsTestingApi] = useState(false);
  const [testStatus, setTestStatus] = useState<'idle' | 'success' | 'error'>('idle');
  const [testErrorMessage, setTestErrorMessage] = useState<string>('');

  // Load API key from account state
  const [isLoadingKeys, setIsLoadingKeys] = useState(false);
  const [loadKeyStatus, setLoadKeyStatus] = useState<'idle' | 'success' | 'error'>('idle');
  const [loadKeyError, setLoadKeyError] = useState('');

  // Get global API key for initial value
  const globalApiKey = useApiKeyStore((state) => state.getApiKey());
  const loadApiKeys = useApiKeyStore((state) => state.loadApiKeys);

  // Auth state
  const { isAuthenticated, user } = useAuthStore();

  // Initialize with global API key
  useEffect(() => {
    if (globalApiKey && !apiKey) {
      setApiKey(globalApiKey);
    }
  }, [globalApiKey]);

  // Custom file upload for CASSIA results files
  const [uploadedFile, setUploadedFile] = useState(null);
  const [isUploading, setIsUploading] = useState(false);

  const onDrop = useCallback(async (acceptedFiles) => {
    const file = acceptedFiles[0];
    if (!file) return;

    setIsUploading(true);
    setError('');
    
    try {
      const text = await file.text();
      const parsed = parseCSV(text);
      
      // Basic validation for CASSIA results format
      if (parsed.length === 0) {
        throw new Error('CSV file is empty');
      }

      const firstRow = parsed[0];
      const requiredColumns = ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types'];
      const missingColumns = requiredColumns.filter(col => !firstRow.hasOwnProperty(col));
      
      if (missingColumns.length > 0) {
        throw new Error(`Missing required columns for CASSIA results: ${missingColumns.join(', ')}`);
      }

      setCsvData(parsed);
      setUploadedFile(file);
      setResults(null);
      setProgress(`Loaded ${parsed.length} rows from CASSIA results file`);
      
    } catch (err) {
      setError(`Error parsing CSV file: ${err.message}`);
      setCsvData(null);
      setUploadedFile(null);
    } finally {
      setIsUploading(false);
    }
  }, []);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      'text/csv': ['.csv']
    },
    multiple: false,
    disabled: isUploading || isProcessing
  });

  const clearFile = () => {
    setCsvData(null);
    setUploadedFile(null);
    setResults(null);
    setError('');
    setProgress('');
  };

  // Load example CASSIA results data
  const loadExampleData = async () => {
    setIsUploading(true);
    setError('');

    try {
      const response = await fetch('/examples/cassia_analysis_full.csv');
      if (!response.ok) {
        throw new Error('Failed to fetch example file');
      }
      const text = await response.text();
      const parsed = parseCSV(text);

      if (parsed.length === 0) {
        throw new Error('Example CSV file is empty');
      }

      setCsvData(parsed);
      const mockFile = new File([text], 'cassia_analysis_full.csv', { type: 'text/csv' });
      setUploadedFile(mockFile);
      setResults(null);
      setProgress(`Loaded ${parsed.length} rows from example CASSIA results`);
      console.log('Example data loaded successfully');
    } catch (err: any) {
      setError(`Error loading example data: ${err.message}`);
      setCsvData(null);
      setUploadedFile(null);
    } finally {
      setIsUploading(false);
    }
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

  const handleProcess = async () => {
    if (!csvData || csvData.length === 0) {
      setError('Please upload a CSV file first');
      return;
    }

    if (!apiKey) {
      setError('Please set your API key first');
      return;
    }

    setIsProcessing(true);
    setError('');
    setResults(null);
    setProgress('Starting annotation merging process...');

    try {
      let result;

      // Use customBaseUrl as provider for custom endpoints
      const effectiveProvider = provider === 'custom' ? customBaseUrl : provider;

      if (processAllLevels) {
        setProgress('Processing all detail levels in parallel...');
        result = await mergeAnnotationsAll({
          csvData,
          provider: effectiveProvider,
          model,
          apiKey,
          additionalContext: additionalContext || null,
          batchSize: parseInt(batchSize),
          reasoningEffort,
          onProgress: (msg) => setProgress(msg)
        });
      } else {
        setProgress(`Processing with ${detailLevel} detail level...`);
        result = await mergeAnnotations({
          csvData,
          provider: effectiveProvider,
          model,
          apiKey,
          additionalContext: additionalContext || null,
          batchSize: parseInt(batchSize),
          detailLevel,
          reasoningEffort,
          onProgress: (msg) => setProgress(msg)
        });
      }

      setResults(result);
      setProgress('Annotation merging completed successfully!');
    } catch (err) {
      setError(`Error during processing: ${err.message}`);
      setProgress('');
    } finally {
      setIsProcessing(false);
    }
  };

  const downloadResults = () => {
    if (!results) return;

    const csvContent = convertToCSV(results);
    const blob = new Blob([csvContent], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `annotation_merging_results_${new Date().toISOString().split('T')[0]}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  const convertToCSV = (data) => {
    if (!data || data.length === 0) return '';
    
    const headers = Object.keys(data[0]);
    const csvRows = [
      headers.join(','),
      ...data.map(row => 
        headers.map(header => {
          const value = row[header];
          if (value === null || value === undefined) return '';
          return `"${String(value).replace(/"/g, '""')}"`;
        }).join(',')
      )
    ];
    
    return csvRows.join('\n');
  };

  return (
    <div className="container mx-auto p-6 max-w-4xl">
      <div className="mb-8">
        <Link
          href="/"
          className="group inline-flex items-center gap-2 px-4 py-2 mb-4 rounded-xl bg-gradient-to-r from-gray-100 to-gray-50 dark:from-gray-800 dark:to-gray-700 border border-gray-200 dark:border-gray-600 shadow-sm hover:shadow-md hover:from-blue-50 hover:to-purple-50 dark:hover:from-blue-900/30 dark:hover:to-purple-900/30 hover:border-blue-300 dark:hover:border-blue-600 transition-all duration-300"
        >
          <ArrowLeft className="h-4 w-4 text-gray-600 dark:text-gray-300 group-hover:text-blue-600 dark:group-hover:text-blue-400 group-hover:-translate-x-1 transition-all duration-300" />
          <span className="text-sm font-medium text-gray-700 dark:text-gray-200 group-hover:text-blue-600 dark:group-hover:text-blue-400">Back</span>
        </Link>
        <h1 className="text-3xl font-bold mb-2">Annotation Merging Agent</h1>
        <p className="text-gray-600">
          Merge and group cell cluster annotations using AI to create broader cell type categories.
          Upload CASSIA results files (not raw marker data) with existing annotations to group them into hierarchical categories.
        </p>
      </div>

      <div className="space-y-6">
        {/* File Upload */}
        <Card className="p-6">
          <h2 className="text-xl font-semibold mb-4">Upload CASSIA Results File</h2>
          
          {!uploadedFile ? (
            <div
              {...getRootProps()}
              className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${
                isDragActive
                  ? 'border-blue-500 bg-blue-50'
                  : 'border-gray-300 hover:border-blue-400'
              } ${isUploading ? 'opacity-50 cursor-not-allowed' : ''}`}
            >
              <input {...getInputProps()} />
              <div className="space-y-4">
                <div className="mx-auto w-12 h-12 bg-blue-100 rounded-lg flex items-center justify-center">
                  {isUploading ? (
                    <div className="w-6 h-6 border-2 border-blue-500 border-t-transparent rounded-full animate-spin" />
                  ) : (
                    <Upload className="h-6 w-6 text-blue-500" />
                  )}
                </div>
                
                <div className="space-y-2">
                  <h3 className="font-medium">
                    {isUploading
                      ? 'Processing file...'
                      : isDragActive
                      ? 'Drop your file here'
                      : 'Upload CASSIA Results File'
                    }
                  </h3>
                  <p className="text-sm text-gray-600">
                    {isUploading
                      ? 'Please wait while we process your file'
                      : 'Drag & drop a CSV file, or click to browse'
                    }
                  </p>
                </div>
                
                <div className="text-xs text-gray-500 space-y-1">
                  <div>📄 Supported format: CSV</div>
                  <div>📊 Expected columns: "True Cell Type", "Predicted Main Cell Type", "Predicted Sub Cell Types"</div>
                </div>
              </div>
            </div>
          ) : (
            <div className="p-4 bg-green-50 border border-green-200 rounded-lg">
              <div className="flex items-center justify-between">
                <div className="flex items-center space-x-3">
                  <CheckCircle className="h-5 w-5 text-green-500" />
                  <div>
                    <h3 className="font-medium text-green-800">{uploadedFile.name}</h3>
                    <p className="text-sm text-green-600">
                      {csvData.length} rows loaded successfully
                    </p>
                  </div>
                </div>
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={clearFile}
                  className="text-green-600 hover:text-green-700"
                >
                  <X className="h-4 w-4" />
                </Button>
              </div>
            </div>
          )}
        </Card>

        {/* Configuration */}
        <Card className="p-6">
          <h2 className="text-xl font-semibold mb-4">Configuration</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <Label htmlFor="detail-level">Detail Level</Label>
              <Select value={detailLevel} onValueChange={setDetailLevel} disabled={processAllLevels}>
                <SelectTrigger>
                  <SelectValue />
                </SelectTrigger>
                <SelectContent>
                  <SelectItem value="broad">Broad (General categories)</SelectItem>
                  <SelectItem value="detailed">Detailed (Intermediate specificity)</SelectItem>
                  <SelectItem value="very_detailed">Very Detailed (High specificity)</SelectItem>
                </SelectContent>
              </Select>
            </div>

            <div>
              <Label htmlFor="batch-size">Batch Size</Label>
              <Input
                id="batch-size"
                type="number"
                value={batchSize}
                onChange={(e) => setBatchSize(e.target.value)}
                min="1"
                max="50"
              />
            </div>

          </div>

          {/* Model Configuration */}
          <div className="mt-4">
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
          </div>

          {/* API Key */}
          <Card className="mt-4">
            <CardHeader>
              <CardTitle className="text-lg">API Key</CardTitle>
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

          <div className="mt-4">
            <Label htmlFor="process-all">
              <input
                id="process-all"
                type="checkbox"
                checked={processAllLevels}
                onChange={(e) => setProcessAllLevels(e.target.checked)}
                className="mr-2"
              />
              Process all detail levels in parallel
            </Label>
          </div>

          <div className="mt-4">
            <Label htmlFor="additional-context">Additional Context (Optional)</Label>
            <Textarea
              id="additional-context"
              value={additionalContext}
              onChange={(e) => setAdditionalContext(e.target.value)}
              placeholder="Provide any additional context about the tissue type, species, or experimental conditions that might help with annotation merging..."
              rows={3}
            />
          </div>
        </Card>

        {/* Process Button */}
        <Card className="p-6">
          <Button 
            onClick={handleProcess}
            disabled={isProcessing || !csvData || !apiKey}
            className="w-full"
          >
            {isProcessing ? 'Processing...' : 'Start Annotation Merging'}
          </Button>
        </Card>

        {/* Progress */}
        {(progress || isProcessing) && (
          <Card className="p-6">
            <ProgressTracker 
              currentStep={progress || 'Processing...'}
              isProcessing={isProcessing}
            />
          </Card>
        )}

        {/* Error Display */}
        {error && (
          <Card className="p-6 border-red-200 bg-red-50">
            <div className="text-red-800">
              <h3 className="font-semibold mb-2">Error</h3>
              <p>{error}</p>
            </div>
          </Card>
        )}

        {/* Results */}
        {results && (
          <Card className="p-6">
            <div className="flex justify-between items-center mb-4">
              <h3 className="text-lg font-semibold">Results</h3>
              <Button onClick={downloadResults}>
                Download CSV
              </Button>
            </div>
            <div className="space-y-2">
              <p>Successfully processed {results.length} clusters</p>
              {processAllLevels && (
                <p className="text-sm text-gray-600">
                  Results include all three detail levels: Merged_Grouping_1 (broad), Merged_Grouping_2 (detailed), Merged_Grouping_3 (very detailed)
                </p>
              )}
            </div>
            
            {/* Preview of results */}
            <div className="mt-4 max-h-64 overflow-auto">
              <table className="w-full text-sm">
                <thead className="bg-gray-50">
                  <tr>
                    {Object.keys(results[0] || {}).slice(0, 5).map(key => (
                      <th key={key} className="px-2 py-1 text-left">{key}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {results.slice(0, 10).map((row, idx) => (
                    <tr key={idx} className="border-t">
                      {Object.values(row).slice(0, 5).map((value, i) => (
                        <td key={i} className="px-2 py-1 truncate max-w-32">
                          {String(value || '')}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
              {results.length > 10 && (
                <p className="text-sm text-gray-500 mt-2">
                  Showing first 10 of {results.length} results
                </p>
              )}
            </div>
          </Card>
        )}
      </div>
    </div>
  );
}