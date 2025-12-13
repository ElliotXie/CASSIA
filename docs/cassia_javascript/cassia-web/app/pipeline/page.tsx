'use client'

import React, { useState } from 'react'
import Link from 'next/link'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { ArrowLeft, Play, Settings, HelpCircle, FileText } from 'lucide-react'
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
import { useResultsStore } from '@/lib/stores/results-store'
import { useAuthStore } from '@/lib/stores/auth-store'
import { ContactDialog } from '@/components/ContactDialog'
import { ReportViewerModal, type BatchReportData } from '@/components/reports'

export default function PipelinePage() {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  const [showReportModal, setShowReportModal] = useState(false)
  
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
    boostIterations,
    boostSearchStrategy,
    setAnalysisConfig
  } = useConfigStore()
  
  const { provider, model, getApiKey, getCustomBaseUrl } = useApiKeyStore()
  const apiKey = getApiKey()
  const customBaseUrl = getCustomBaseUrl()

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

  const { saveResult } = useResultsStore()
  const { isAuthenticated } = useAuthStore()

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
        boostIterations,
        boostSearchStrategy,
        apiKey,
        provider,
        customBaseUrl,
        onProgress: (progressData) => {
          updateProgress(progressData.progress, progressData.step)
        },
        onLog: (message) => {
          addLog(message)
        },
        onComplete: async (finalResults) => {
          // Set final results
          setResults({
            clusters: finalResults.summary.totalClusters,
            avgScore: finalResults.metadata?.averageScore || 0,
            totalFiles: finalResults.summary.totalFiles,
            lowScoreClusters: finalResults.summary.lowScoreClusters,
            downloadUrls: finalResults.downloadUrls,
            finalResults
          })

          // Save to database for authenticated users
          if (isAuthenticated) {
            try {
              await saveResult({
                analysis_type: 'batch',
                title: `${tissue} - ${species} (${finalResults.summary.totalClusters} clusters)`,
                description: outputName || 'Pipeline analysis',
                results: finalResults,
                settings: { provider, model, tissue, species, scoreThreshold }
              })
              addLog('💾 Results saved to your dashboard')
            } catch (err) {
              console.error('Failed to save results:', err)
            }
          }
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
              <Link
                href="/"
                className="group flex items-center gap-2 px-4 py-2 rounded-xl bg-gradient-to-r from-gray-100 to-gray-50 dark:from-gray-800 dark:to-gray-700 border border-gray-200 dark:border-gray-600 shadow-sm hover:shadow-md hover:from-blue-50 hover:to-purple-50 dark:hover:from-blue-900/30 dark:hover:to-purple-900/30 hover:border-blue-300 dark:hover:border-blue-600 transition-all duration-300"
              >
                <ArrowLeft className="h-4 w-4 text-gray-600 dark:text-gray-300 group-hover:text-blue-600 dark:group-hover:text-blue-400 group-hover:-translate-x-1 transition-all duration-300" />
                <span className="text-sm font-medium text-gray-700 dark:text-gray-200 group-hover:text-blue-600 dark:group-hover:text-blue-400">Back</span>
              </Link>
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

            {/* Model Selection Matrix - Show for OpenRouter, OpenAI, and Anthropic */}
            {['openrouter', 'openai', 'anthropic'].includes(provider) && (
              <ModelSelectionMatrix
                lockedProvider={provider !== 'openrouter' ? (provider as 'openai' | 'anthropic') : null}
              />
            )}

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
                      placeholder="large_intestine"
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
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Boost Iterations</label>
                        <Input
                          type="number"
                          min="1"
                          max="10"
                          value={boostIterations}
                          onChange={(e) => setAnalysisConfig({ boostIterations: parseInt(e.target.value) || 5 })}
                        />
                        <p className="text-xs text-muted-foreground">Max iterations for annotation boost analysis</p>
                      </div>
                      <div className="space-y-2">
                        <label className="text-sm font-medium">Boost Search Strategy</label>
                        <select
                          value={boostSearchStrategy}
                          onChange={(e) => setAnalysisConfig({ boostSearchStrategy: e.target.value as 'breadth' | 'depth' })}
                          className="w-full h-10 px-3 rounded-md border border-input bg-background text-sm"
                        >
                          <option value="breadth">Breadth-First</option>
                          <option value="depth">Depth-First</option>
                        </select>
                        <p className="text-xs text-muted-foreground">Search strategy for marker gene analysis</p>
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

      {/* Interactive Report Viewer Modal */}
      {results && results.finalResults && (
        <ReportViewerModal
          open={showReportModal}
          onOpenChange={setShowReportModal}
          data={{
            batch: results.finalResults.batchResults?.map((row: any) => ({
              clusterId: row.cluster_id || row.Cluster || row['Cluster ID'] || `Cluster ${row.index || 0}`,
              mainType: row.predicted_cell_type || row.main_cell_type || row['Predicted Cell Type'] || 'Unknown',
              subTypes: row.sub_cell_types || row['Sub Cell Types'] || '',
              mixedTypes: row.possible_mixed_cell_types || row['Mixed Cell Types'] || '',
              markers: row.markers || row.marker_list || row.Markers || '',
              markerCount: row.marker_count || row['Marker Count'] || 0,
              score: row.score || row.quality_score || row.Score,
              tissue: row.tissue || tissue,
              species: row.species || species,
              model: row.model || pipelineModels.annotation.model,
              provider: row.provider || 'openrouter',
              iterations: row.iterations || 1,
              conversationHistory: row.conversation_history || row.conversationHistory || '',
              mergedGrouping1: row.merged_grouping_1 || row['Merged Grouping 1'] || '',
              mergedGrouping2: row.merged_grouping_2 || row['Merged Grouping 2'] || '',
              mergedGrouping3: row.merged_grouping_3 || row['Merged Grouping 3'] || ''
            })) || []
          }}
          title="CASSIA Analysis Report"
        />
      )}
    </div>
  )
}