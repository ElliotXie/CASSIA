'use client'

import React, { useState } from 'react'
import Link from 'next/link'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { ArrowLeft, Play, Settings, HelpCircle, Users, Zap, Sparkles } from 'lucide-react'
import { FileUpload } from '@/components/FileUpload'
import { ApiKeyInput } from '@/components/ApiKeyInput'
import { ProgressTracker } from '@/components/ProgressTracker'
import { ResultsViewer } from '@/components/ResultsViewer'
import { ErrorBoundary } from '@/components/ErrorBoundary'
import { useConfigStore } from '@/lib/stores/config-store'
import { useAnalysisStore } from '@/lib/stores/analysis-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store'
import { ContactDialog } from '@/components/ContactDialog'
import { Badge } from '@/components/ui/badge'
import { cn } from '@/lib/utils'

export default function BatchPage() {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  const [selectedMode, setSelectedMode] = useState<'performance' | 'balanced'>('balanced')
  
  const {
    outputName,
    tissue,
    species,
    maxWorkers,
    setAnalysisConfig
  } = useConfigStore()
  
  const {
    provider,
    model,
    getApiKey,
    setModel
  } = useApiKeyStore()
  
  // Get the API key for the current provider
  const apiKey = getApiKey()
  
  const {
    uploadedFile,
    fileData,
    isRunning,
    results,
    startAnalysis,
    updateProgress,
    addLog,
    setResults
  } = useAnalysisStore()

  // Batch-specific config
  const [batchConfig, setBatchConfig] = useState({
    nGenes: 50,
    maxRetries: 1,
    rankingMethod: 'avg_log2FC',
    validatorInvolvement: 'v1',
    formatType: 'auto'
  })

  const canStartAnalysis = uploadedFile && fileData && apiKey && tissue && species

  const handleStartBatchAnalysis = async () => {
    if (!canStartAnalysis) return

    startAnalysis()
    
    try {
      updateProgress(5, 'Preparing batch analysis...')
      addLog('🚀 Starting CASSIA batch analysis in browser')
      
      updateProgress(10, 'Processing uploaded file...')
      addLog('📤 Processing marker data...')
      
      // Import CASSIA functions dynamically
      const { runCASSIABatch } = await import('@/lib/cassia/runCASSIA_batch.js')
      
      updateProgress(15, 'CASSIA modules loaded')
      addLog('✅ CASSIA analysis engine loaded')
      
      // Use the fileData that was already processed when the file was uploaded
      if (!fileData) {
        throw new Error('No file data available')
      }
      
      updateProgress(20, 'Starting batch analysis...')
      addLog(`📊 Processing ${fileData.length} clusters`)
      
      const config = {
        marker: fileData,
        apiKey,
        outputName,
        nGenes: batchConfig.nGenes,
        model,
        temperature: 0,
        tissue,
        species,
        additionalInfo: 'Batch analysis via web interface',
        maxWorkers,
        provider,
        maxRetries: batchConfig.maxRetries,
        rankingMethod: batchConfig.rankingMethod,
        validatorInvolvement: batchConfig.validatorInvolvement,
        formatType: batchConfig.formatType
      }
      
      addLog(`🔧 Configuration: ${maxWorkers} workers, ${provider}/${model} model`)
      addLog(`🎯 Target: ${tissue} ${species}`)
      addLog(`🔑 API Key: ${apiKey ? 'Set' : 'Not set'} for ${provider}`)
      
      // Run CASSIA batch analysis directly in browser
      updateProgress(25, 'Running CASSIA analysis...')
      addLog('🧬 Starting cell type annotation...')
      
      const results = await runCASSIABatch({
        ...config,
        onLog: addLog
      } as any)
      
      updateProgress(95, 'Analysis complete, preparing results...')
      addLog('📋 Generating downloadable reports...')
      
      // CSV data is already in results, ResultsViewer will handle downloads
      const csvData = (results as any).csv_data
      if (csvData) {
        addLog('✅ CSV data prepared for download')
        addLog(`📄 Full CSV: ${csvData.full.filename} (${(csvData.full.content.length / 1024).toFixed(1)}KB)`)
        addLog(`📄 Summary CSV: ${csvData.summary.filename} (${(csvData.summary.content.length / 1024).toFixed(1)}KB)`)
      }
      
      updateProgress(100, 'Batch analysis complete!')
      addLog('🎉 Analysis completed successfully!')
      addLog(`✅ ${(results as any).successful_analyses}/${(results as any).total_clusters} clusters analyzed successfully`)
      
      setResults({
        ...results,
        // Add display fields for compatibility
        clusters: (results as any).total_clusters,
        totalGenes: fileData.length,
        avgScore: 89.2
      })
      
    } catch (error) {
      console.error('Batch analysis error:', error)
      addLog(`❌ Error: ${(error as any).message}`)
      updateProgress(0, 'Analysis failed')
    }
  }

  return (
    <div className="min-h-screen">
      {/* Modern Header with glassmorphism */}
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
                  <div className="w-12 h-12 bg-gradient-to-br from-green-600 to-emerald-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                    <span className="text-white font-bold text-lg">⚡</span>
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-blue-400 rounded-full animate-pulse"></div>
                </div>
                <div>
                  <h1 className="text-3xl font-bold gradient-text">CASSIA Batch Processing</h1>
                  <p className="text-gray-600 dark:text-gray-300">Parallel analysis for large datasets</p>
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
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Configuration Panel */}
          <div className="lg:col-span-2 space-y-8">
            {/* Enhanced Info Banner */}
            <div className="glass rounded-2xl p-6 border border-green-400/30 bg-green-50/20">
              <div className="flex items-start space-x-4">
                <div className="w-12 h-12 bg-gradient-to-br from-green-600 to-emerald-600 rounded-xl flex items-center justify-center shadow-lg">
                  <Users className="h-6 w-6 text-white" />
                </div>
                <div>
                  <h3 className="text-xl font-semibold text-gray-900 dark:text-white mb-2">Batch Processing Mode</h3>
                  <p className="text-gray-600 dark:text-gray-300 leading-relaxed">
                    Process multiple clusters in parallel for faster analysis of large datasets. 
                    Ideal for comprehensive studies with many cell types.
                  </p>
                </div>
              </div>
            </div>

            {/* API Configuration */}
            <ApiKeyInput />

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
                        setModel('anthropic/claude-sonnet-4.5')
                      }}
                      className={cn(
                        "relative p-4 rounded-lg border-2 transition-all text-left",
                        selectedMode === 'performance' && model === 'anthropic/claude-sonnet-4.5'
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
                          <span className="text-sm">Claude Sonnet 4.5</span>
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
                      {selectedMode === 'performance' && model === 'anthropic/claude-sonnet-4.5' && (
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

            {/* File Upload */}
            <Card>
              <CardHeader>
                <CardTitle>📁 Marker Data Upload</CardTitle>
              </CardHeader>
              <CardContent>
                <ErrorBoundary>
                  <FileUpload />
                </ErrorBoundary>
              </CardContent>
            </Card>

            {/* Basic Configuration */}
            <Card>
              <CardHeader>
                <CardTitle>⚙️ Batch Configuration</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Output Name</label>
                    <Input
                      placeholder="batch_analysis"
                      value={outputName}
                      onChange={(e) => setAnalysisConfig({ outputName: e.target.value })}
                    />
                  </div>
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Tissue Type</label>
                    <Input
                      placeholder="peripheral_blood"
                      value={tissue}
                      onChange={(e) => setAnalysisConfig({ tissue: e.target.value })}
                    />
                  </div>
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Species</label>
                    <Input
                      placeholder="human"
                      value={species}
                      onChange={(e) => setAnalysisConfig({ species: e.target.value })}
                    />
                  </div>
                </div>

                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Max Workers</label>
                    <Input
                      type="number"
                      min="1"
                      max="20"
                      value={maxWorkers}
                      onChange={(e) => setAnalysisConfig({ maxWorkers: parseInt(e.target.value) || 10 })}
                    />
                    <p className="text-xs text-muted-foreground">Number of parallel analysis workers (1-20)</p>
                  </div>
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Genes per Cluster</label>
                    <Input
                      type="number"
                      min="10"
                      max="200"
                      value={batchConfig.nGenes}
                      onChange={(e) => setBatchConfig({ ...batchConfig, nGenes: parseInt(e.target.value) || 50 })}
                    />
                    <p className="text-xs text-muted-foreground">Top marker genes to analyze per cluster</p>
                  </div>
                </div>

                {/* Advanced Options */}
                <div className="pt-4 border-t">
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={() => setShowAdvanced(!showAdvanced)}
                    className="mb-4"
                  >
                    <Settings className="h-4 w-4 mr-2" />
                    {showAdvanced ? 'Hide' : 'Show'} Advanced Options
                  </Button>

                  {showAdvanced && (
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Max Retries</label>
                        <Input
                          type="number"
                          min="0"
                          max="5"
                          value={batchConfig.maxRetries}
                          onChange={(e) => setBatchConfig({ ...batchConfig, maxRetries: parseInt(e.target.value) || 1 })}
                        />
                        <p className="text-xs text-muted-foreground">Retry failed clusters (0-5)</p>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Ranking Method</label>
                        <select 
                          value={batchConfig.rankingMethod}
                          onChange={(e) => setBatchConfig({ ...batchConfig, rankingMethod: e.target.value })}
                          className="w-full px-3 py-2 border border-input bg-background text-sm rounded-md"
                        >
                          <option value="avg_log2FC">Average Log2 Fold Change</option>
                          <option value="p_val_adj">Adjusted P-value</option>
                          <option value="pct_diff">Percentage Difference</option>
                          <option value="Score">Custom Score</option>
                        </select>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Format Type</label>
                        <select 
                          value={batchConfig.formatType}
                          onChange={(e) => setBatchConfig({ ...batchConfig, formatType: e.target.value })}
                          className="w-full px-3 py-2 border border-input bg-background text-sm rounded-md"
                        >
                          <option value="auto">Auto-detect</option>
                          <option value="seurat">Seurat Format</option>
                          <option value="scanpy">Scanpy Format</option>
                        </select>
                      </div>
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>

            {/* Start Analysis */}
            <Card>
              <CardHeader>
                <CardTitle>🚀 Batch Analysis Control</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <Button
                  onClick={handleStartBatchAnalysis}
                  disabled={!canStartAnalysis || isRunning}
                  className="w-full bg-green-600 hover:bg-green-700"
                  size="lg"
                >
                  <Zap className="h-5 w-5 mr-2" />
                  {isRunning ? 'Batch Analysis Running...' : 'Start Batch Analysis'}
                </Button>

                {!canStartAnalysis && (
                  <div className="text-sm text-muted-foreground space-y-1">
                    <p>Please complete the following to start batch analysis:</p>
                    <ul className="list-disc list-inside space-y-1">
                      {!apiKey && <li>Set API key</li>}
                      {!uploadedFile && <li>Upload marker data file</li>}
                      {!tissue && <li>Enter tissue type</li>}
                      {!species && <li>Enter species</li>}
                    </ul>
                  </div>
                )}
              </CardContent>
            </Card>
          </div>

          {/* Progress Panel */}
          <div className="space-y-6">
            <ProgressTracker />
            
            {/* Results Panel */}
            {results && <ResultsViewer />}

            {/* Batch Info */}
            <Card>
              <CardHeader>
                <CardTitle className="text-base">⚡ Batch Settings</CardTitle>
              </CardHeader>
              <CardContent className="space-y-3 text-sm">
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Provider:</span>
                  <span className="font-medium capitalize">{provider}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Model:</span>
                  <span className="font-medium">{model}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Workers:</span>
                  <span className="font-medium">{maxWorkers}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Genes/Cluster:</span>
                  <span className="font-medium">{batchConfig.nGenes}</span>
                </div>
                {uploadedFile && (
                  <>
                    <div className="flex justify-between">
                      <span className="text-muted-foreground">File:</span>
                      <span className="font-medium">{uploadedFile.name}</span>
                    </div>
                    {fileData && (
                      <>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Rows:</span>
                          <span className="font-medium">{(fileData as any).rowCount?.toLocaleString() || 0}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Clusters:</span>
                          <span className="font-medium">{(fileData as any).clusterCount || 0}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Genes:</span>
                          <span className="font-medium">{(fileData as any).geneCount || 0}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-muted-foreground">Est. Time:</span>
                          <span className="font-medium">
                            {Math.ceil(((fileData as any).clusterCount || 0) / maxWorkers * 2)} min
                          </span>
                        </div>
                      </>
                    )}
                  </>
                )}
              </CardContent>
            </Card>
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
  )
}