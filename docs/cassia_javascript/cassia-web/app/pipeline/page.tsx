'use client'

import React, { useState } from 'react'
import Link from 'next/link'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { ArrowLeft, Play, Settings, HelpCircle } from 'lucide-react'
import { FileUpload } from '@/components/FileUpload'
import { ApiKeyInput } from '@/components/ApiKeyInput'
import { ProgressTracker } from '@/components/ProgressTracker'
import { ResultsViewer } from '@/components/ResultsViewer'
import { ErrorBoundary } from '@/components/ErrorBoundary'
import { ModelSelectionMatrix } from '@/components/ModelSelectionMatrix'
import { ResultsDownloader } from '@/components/ResultsDownloader'
import { useConfigStore } from '@/lib/stores/config-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store'
import { useAnalysisStore } from '@/lib/stores/analysis-store'
import { ContactDialog } from '@/components/ContactDialog'

export default function PipelinePage() {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  
  const {
    pipelineModels,
    outputName,
    tissue,
    species,
    scoreThreshold,
    maxWorkers,
    mergeAnnotations,
    additionalInfo,
    maxRetries,
    setAnalysisConfig
  } = useConfigStore()
  
  const { provider, model, getApiKey } = useApiKeyStore()
  const apiKey = getApiKey()
  
  const {
    uploadedFile,
    fileData,
    fileMetadata,
    isRunning,
    results,
    startAnalysis,
    updateProgress,
    addLog,
    setResults
  } = useAnalysisStore()

  const canStartAnalysis = uploadedFile && fileData && apiKey && tissue && species

  const handleStartAnalysis = async () => {
    if (!canStartAnalysis) return

    startAnalysis()
    
    try {
      // Import the pipeline orchestrator
      const { runCASSIAPipeline } = await import('@/lib/cassia/runCASSIA_pipeline')
      
      // Configure the pipeline
      const pipelineConfig = {
        marker: fileData,
        outputName: outputName || 'cassia_analysis',
        tissue,
        species,
        models: pipelineModels,
        scoreThreshold,
        maxWorkers,
        maxRetries,
        additionalInfo,
        apiKey,
        onProgress: (progressData) => {
          updateProgress(progressData.progress, progressData.step)
        },
        onLog: (message) => {
          addLog(message)
        },
        onComplete: (finalResults) => {
          // Set final results
          setResults({
            clusters: finalResults.summary.totalClusters,
            avgScore: finalResults.metadata?.averageScore || 0,
            totalFiles: finalResults.summary.totalFiles,
            lowScoreClusters: finalResults.summary.lowScoreClusters,
            downloadUrls: finalResults.downloadUrls,
            finalResults
          })
        },
        onError: (error) => {
          console.error('Pipeline error:', error)
          addLog(`❌ Pipeline failed: ${error.message}`)
        }
      }
      
      // Run the complete pipeline
      const results = await runCASSIAPipeline(pipelineConfig)
      
      addLog('🎉 Pipeline completed successfully!')
      
    } catch (error) {
      console.error('Pipeline execution failed:', error)
      addLog(`❌ Pipeline failed: ${error.message}`)
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
                  <div className="w-12 h-12 bg-gradient-to-br from-blue-600 to-purple-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                    <span className="text-white font-bold text-lg">🧬</span>
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-yellow-400 rounded-full animate-pulse"></div>
                </div>
                <div>
                  <h1 className="text-3xl font-bold gradient-text">CASSIA Pipeline</h1>
                  <p className="text-gray-600 dark:text-gray-300">Complete end-to-end analysis</p>
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
          <div className="lg:col-span-2 space-y-6">
            {/* API Configuration */}
            <ApiKeyInput />

            {/* Model Selection Matrix - Only show when OpenRouter is selected */}
            {provider === 'openrouter' && <ModelSelectionMatrix />}

            {/* File Upload */}
            <Card>
              <CardHeader>
                <CardTitle>📁 Data Upload</CardTitle>
              </CardHeader>
              <CardContent>
                <ErrorBoundary>
                  <FileUpload />
                </ErrorBoundary>
              </CardContent>
            </Card>

            {/* Analysis Configuration */}
            <Card>
              <CardHeader>
                <CardTitle>⚙️ Analysis Configuration</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Output Name</label>
                    <Input
                      placeholder="my_analysis"
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
                        <label className="text-sm font-medium">Score Threshold</label>
                        <Input
                          type="number"
                          min="0"
                          max="100"
                          value={scoreThreshold}
                          onChange={(e) => setAnalysisConfig({ scoreThreshold: parseInt(e.target.value) || 75 })}
                        />
                        <p className="text-xs text-muted-foreground">Clusters below this score will trigger boost analysis</p>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Max Workers</label>
                        <Input
                          type="number"
                          min="1"
                          max="8"
                          value={maxWorkers}
                          onChange={(e) => setAnalysisConfig({ maxWorkers: parseInt(e.target.value) || 4 })}
                        />
                        <p className="text-xs text-muted-foreground">Parallel processing workers</p>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Merge Annotations</label>
                        <div className="flex items-center space-x-2">
                          <input
                            type="checkbox"
                            checked={mergeAnnotations}
                            onChange={(e) => setAnalysisConfig({ mergeAnnotations: e.target.checked })}
                            className="rounded"
                          />
                          <span className="text-sm">Enable LLM-powered annotation grouping</span>
                        </div>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Additional Info</label>
                        <Input
                          placeholder="Additional context for analysis..."
                          value={additionalInfo}
                          onChange={(e) => setAnalysisConfig({ additionalInfo: e.target.value })}
                        />
                        <p className="text-xs text-muted-foreground">Extra information to help with analysis</p>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Max Retries</label>
                        <Input
                          type="number"
                          min="1"
                          max="5"
                          value={maxRetries}
                          onChange={(e) => setAnalysisConfig({ maxRetries: parseInt(e.target.value) || 1 })}
                        />
                        <p className="text-xs text-muted-foreground">Number of retry attempts for failed requests</p>
                      </div>
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>

            {/* Start Analysis */}
            <Card>
              <CardHeader>
                <CardTitle>🚀 Analysis Control</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <Button
                  onClick={handleStartAnalysis}
                  disabled={!canStartAnalysis || isRunning}
                  className="w-full"
                  size="lg"
                >
                  <Play className="h-5 w-5 mr-2" />
                  {isRunning ? 'Analysis Running...' : 'Start Pipeline Analysis'}
                </Button>

                {!canStartAnalysis && (
                  <div className="text-sm text-muted-foreground space-y-1">
                    <p>Please complete the following to start analysis:</p>
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

            {/* Results Download Panel */}
            {results && results.finalResults && (
              <ErrorBoundary>
                <ResultsDownloader 
                  finalResults={results.finalResults}
                  isVisible={!isRunning && results.finalResults}
                />
              </ErrorBoundary>
            )}

            {/* Quick Info */}
            <Card>
              <CardHeader>
                <CardTitle className="text-base">📊 Analysis Info</CardTitle>
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
                
                {/* Agent Model Configuration */}
                <div className="pt-2 border-t">
                  <div className="mb-2">
                    <span className="text-muted-foreground text-xs font-medium">Agent Models:</span>
                  </div>
                  <div className="space-y-2">
                    <div className="flex justify-between items-center">
                      <span className="text-muted-foreground text-xs">🔍 Annotation:</span>
                      <span className="font-medium text-xs">{pipelineModels.annotation.model}</span>
                    </div>
                    <div className="flex justify-between items-center">
                      <span className="text-muted-foreground text-xs">📊 Scoring:</span>
                      <span className="font-medium text-xs">{pipelineModels.scoring.model}</span>
                    </div>
                    <div className="flex justify-between items-center">
                      <span className="text-muted-foreground text-xs">🚀 Boost:</span>
                      <span className="font-medium text-xs">{pipelineModels.annotationBoost.model}</span>
                    </div>
                  </div>
                </div>

                {uploadedFile && (
                  <>
                    <div className="pt-2 border-t">
                      <div className="flex justify-between">
                        <span className="text-muted-foreground">File:</span>
                        <span className="font-medium">{uploadedFile.name}</span>
                      </div>
                      {fileMetadata && (
                        <>
                          <div className="flex justify-between">
                            <span className="text-muted-foreground">Rows:</span>
                            <span className="font-medium">{fileMetadata.rowCount?.toLocaleString() || 0}</span>
                          </div>
                          <div className="flex justify-between">
                            <span className="text-muted-foreground">Clusters:</span>
                            <span className="font-medium">{fileMetadata.clusterCount || 0}</span>
                          </div>
                          <div className="flex justify-between">
                            <span className="text-muted-foreground">Genes:</span>
                            <span className="font-medium">{fileMetadata.geneCount || 0}</span>
                          </div>
                        </>
                      )}
                    </div>
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