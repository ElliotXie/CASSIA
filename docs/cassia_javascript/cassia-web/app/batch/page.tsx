'use client'

import React, { useState } from 'react'
import Link from 'next/link'
import dynamic from 'next/dynamic'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { ArrowLeft, Play, Settings, HelpCircle, Users, Zap, Sparkles, FileText, Brain, Square } from 'lucide-react'
import { modelSupportsReasoning, getReasoningEffortOptions, ReasoningEffort } from '@/lib/config/model-presets'
import { ErrorBoundary } from '@/components/ErrorBoundary'

// Lazy-load heavy components to reduce initial page load time
const FileUpload = dynamic(() => import('@/components/FileUpload').then(mod => ({ default: mod.FileUpload })), { ssr: false })
const ApiKeyInput = dynamic(() => import('@/components/ApiKeyInput').then(mod => ({ default: mod.ApiKeyInput })), { ssr: false })
const ProgressTracker = dynamic(() => import('@/components/ProgressTracker').then(mod => ({ default: mod.ProgressTracker })), { ssr: false })
const ResultsViewer = dynamic(() => import('@/components/ResultsViewer').then(mod => ({ default: mod.ResultsViewer })), { ssr: false })
import { useConfigStore } from '@/lib/stores/config-store'
import { useAnalysisStore } from '@/lib/stores/analysis-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store'
import { Badge } from '@/components/ui/badge'
import { cn } from '@/lib/utils'

// Lazy-load components that are only shown after user interaction
const ContactDialog = dynamic(() => import('@/components/ContactDialog').then(mod => mod.ContactDialog), { ssr: false })
const ReportViewerModal = dynamic(() => import('@/components/reports').then(mod => mod.ReportViewerModal), { ssr: false })

export default function BatchPage() {
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  const [showReportModal, setShowReportModal] = useState(false)
  const [selectedMode, setSelectedMode] = useState<'performance' | 'balanced' | undefined>('balanced')
  const [useCustomModel, setUseCustomModel] = useState(false)
  const [customModel, setCustomModel] = useState('')
  const [openaiIdentityVerified, setOpenaiIdentityVerified] = useState(false)

  const {
    outputName,
    tissue,
    species,
    maxWorkers,
    additionalInfo,
    setAnalysisConfig
  } = useConfigStore()
  
  const {
    provider,
    model,
    reasoningEffort,
    getApiKey,
    setModel,
    setReasoningEffort,
    getCustomBaseUrl
  } = useApiKeyStore()

  const customBaseUrl = getCustomBaseUrl()

  const reasoningOptions = getReasoningEffortOptions()
  
  // Get the API key for the current provider
  const apiKey = getApiKey()
  
  const {
    uploadedFile,
    fileData,
    isRunning,
    results,
    startAnalysis,
    stopAnalysis,
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

  // Model options for each provider (same as CASSIA Pipeline)
  const modelOptions = {
    openrouter: [
      { id: 'anthropic/claude-sonnet-4.5', name: 'Claude Sonnet 4.5', cost: 'high', speed: 'medium' },
      { id: 'anthropic/claude-haiku-4.5', name: 'Claude Haiku 4.5', cost: 'low', speed: 'fast' },
      { id: 'anthropic/claude-opus-4.5', name: 'Claude Opus 4.5', cost: 'high', speed: 'slow' },
      { id: 'openai/gpt-5.2', name: 'GPT-5.2', cost: 'high', speed: 'medium' },
      { id: 'openai/gpt-5-mini', name: 'GPT-5 Mini', cost: 'low', speed: 'fast' },
      { id: 'openai/gpt-4o', name: 'GPT-4o', cost: 'medium', speed: 'fast' },
      { id: 'google/gemini-3-flash-preview', name: 'Gemini 3 Flash', cost: 'low', speed: 'fast' },
      { id: 'google/gemini-2.5-flash', name: 'Gemini 2.5 Flash', cost: 'low', speed: 'fast' },
      { id: 'google/gemini-2.5-pro', name: 'Gemini 2.5 Pro', cost: 'high', speed: 'medium' },
      { id: 'deepseek/deepseek-chat-v3.1', name: 'DeepSeek Chat V3.1', cost: 'very low', speed: 'fast' },
      { id: 'meta-llama/llama-4-maverick', name: 'Llama 4 Maverick', cost: 'low', speed: 'fast' },
      { id: 'x-ai/grok-4.1-fast', name: 'Grok 4.1 Fast', cost: 'medium', speed: 'fast' }
    ],
    openai: [
      { id: 'gpt-5.2', name: 'GPT-5.2', cost: 'high', speed: 'medium' },
      { id: 'gpt-5-mini', name: 'GPT-5 Mini', cost: 'low', speed: 'fast' },
      { id: 'gpt-4o', name: 'GPT-4o', cost: 'medium', speed: 'fast' }
    ],
    anthropic: [
      { id: 'claude-sonnet-4.5', name: 'Claude Sonnet 4.5', cost: 'high', speed: 'medium' },
      { id: 'claude-haiku-4.5', name: 'Claude Haiku 4.5', cost: 'low', speed: 'fast' },
      { id: 'claude-opus-4.5', name: 'Claude Opus 4.5', cost: 'high', speed: 'slow' }
    ],
    custom: []  // User enters model name manually
  }

  const canStartAnalysis = uploadedFile && fileData && apiKey && tissue && species

  // Preload analysis modules on button hover for instant start
  const preloadAnalysisModules = () => {
    import('@/lib/cassia/runCASSIA_batch.js')
  }

  const handleStartBatchAnalysis = async () => {
    if (!canStartAnalysis) return

    const controller = startAnalysis()
    const signal = controller.signal

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
        additionalInfo: additionalInfo || '',
        maxWorkers,
        provider,
        customBaseUrl,  // Custom provider base URL
        maxRetries: batchConfig.maxRetries,
        rankingMethod: batchConfig.rankingMethod,
        validatorInvolvement: batchConfig.validatorInvolvement,
        formatType: batchConfig.formatType,
        reasoningEffort: reasoningEffort
      }
      
      addLog(`🔧 Configuration: ${maxWorkers} workers, ${provider}/${model} model`)
      addLog(`🧠 Reasoning effort: ${reasoningEffort || 'none'}`)
      addLog(`🎯 Target: ${tissue} ${species}`)
      addLog(`🔑 API Key: ${apiKey ? 'Set' : 'Not set'} for ${provider}`)
      
      // Run CASSIA batch analysis directly in browser
      updateProgress(25, 'Running CASSIA analysis...')
      addLog('🧬 Starting cell type annotation...')
      
      const results = await runCASSIABatch({
        ...config,
        onLog: addLog,
        signal
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
      if ((error as any).name === 'AbortError' || signal.aborted) {
        addLog('🛑 Analysis was stopped by user')
        updateProgress(0, 'Analysis stopped')
      } else {
        console.error('Batch analysis error:', error)
        addLog(`❌ Error: ${(error as any).message}`)
        updateProgress(0, 'Analysis failed')
      }
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
                  <div className="w-12 h-12 bg-gradient-to-br from-green-600 to-emerald-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                    <span className="text-white font-bold text-lg">⚡</span>
                  </div>
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
                        setModel('google/gemini-3-flash-preview')
                      }}
                      className={cn(
                        "relative p-4 rounded-lg border-2 transition-all text-left",
                        selectedMode === 'balanced' && model === 'google/gemini-3-flash-preview'
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
                          <span className="text-sm">Gemini 3 Flash</span>
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
                      {selectedMode === 'balanced' && model === 'google/gemini-3-flash-preview' && (
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

            {/* Model Selection Dropdown */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Settings className="h-5 w-5" />
                  Model Selection
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Select Model</label>
                    {provider === 'custom' ? (
                      // Custom provider: always show model input
                      <div className="space-y-2">
                        <Input
                          placeholder="Enter model name (e.g., meta-llama/Llama-3-70b-chat-hf)"
                          value={model}
                          onChange={(e) => setModel(e.target.value)}
                        />
                        <p className="text-xs text-muted-foreground">
                          Enter the model name as required by your custom endpoint
                        </p>
                      </div>
                    ) : (
                      // Standard providers: show dropdown
                      <>
                        <Select
                          value={model}
                          onValueChange={(value) => {
                            setModel(value)
                            setSelectedMode(undefined)
                          }}
                        >
                          <SelectTrigger>
                            <SelectValue placeholder="Choose a model" />
                          </SelectTrigger>
                          <SelectContent>
                            {modelOptions[provider as keyof typeof modelOptions]?.map((m) => (
                              <SelectItem key={m.id} value={m.id}>
                                <div className="flex items-center gap-2">
                                  <span>{m.name}</span>
                                  <span className="text-xs text-muted-foreground">
                                    (${m.cost} cost, {m.speed})
                                  </span>
                                </div>
                              </SelectItem>
                            ))}
                          </SelectContent>
                        </Select>
                        <p className="text-xs text-muted-foreground">
                          Currently selected: <span className="font-medium">{useCustomModel ? customModel || '(enter model name)' : model}</span>
                        </p>

                        {/* Custom Model Option */}
                        <div className="flex items-center gap-2 mt-3">
                          <input
                            type="checkbox"
                            id="useCustomModel"
                            checked={useCustomModel}
                            onChange={(e) => {
                              setUseCustomModel(e.target.checked)
                              if (!e.target.checked) {
                                setCustomModel('')
                              }
                            }}
                            className="rounded border-gray-300"
                          />
                          <label htmlFor="useCustomModel" className="text-sm cursor-pointer">
                            Use custom model name
                          </label>
                        </div>
                        {useCustomModel && (
                          <Input
                            placeholder="Enter model name (e.g., gpt-4-turbo, claude-3-opus)"
                            value={customModel}
                            onChange={(e) => {
                              setCustomModel(e.target.value)
                              setModel(e.target.value)
                            }}
                            className="mt-2"
                          />
                        )}
                      </>
                    )}
                  </div>

                  {/* Reasoning Effort Section - with OpenAI verification requirement */}
                  {provider === 'openai' ? (
                    // OpenAI-specific: requires identity verification for reasoning
                    <div className="space-y-3">
                      {/* Info note */}
                      <div className="p-3 bg-amber-50 border border-amber-200 rounded-md text-sm">
                        <p className="text-amber-800">
                          <strong>Note:</strong> OpenAI requires identity verification to use reasoning models (GPT-5 with high/medium/low effort).
                          If you haven&apos;t verified, use <strong>OpenRouter</strong> provider instead for GPT-5 reasoning.
                        </p>
                      </div>

                      {/* Identity verification checkbox */}
                      <div className="flex items-center gap-2">
                        <input
                          type="checkbox"
                          id="openaiIdentityVerified"
                          checked={openaiIdentityVerified}
                          onChange={(e) => {
                            setOpenaiIdentityVerified(e.target.checked)
                            if (!e.target.checked) {
                              // Reset reasoning to none when unchecked
                              setReasoningEffort(null)
                            }
                          }}
                          className="rounded border-gray-300"
                        />
                        <label htmlFor="openaiIdentityVerified" className="text-sm cursor-pointer">
                          I have verified my identity with OpenAI
                        </label>
                      </div>

                      {/* Show reasoning dropdown only if verified AND using GPT-5 model */}
                      {openaiIdentityVerified && (model.toLowerCase().includes('gpt-5') || model.toLowerCase().includes('gpt5')) && (
                        <div className="space-y-2">
                          <label className="text-sm font-medium flex items-center gap-1">
                            <Brain className="h-4 w-4" />
                            Reasoning Effort
                          </label>
                          <Select
                            value={reasoningEffort || 'none'}
                            onValueChange={(value) => {
                              setReasoningEffort(value === 'none' ? null : value as ReasoningEffort)
                            }}
                          >
                            <SelectTrigger>
                              <SelectValue placeholder="Choose reasoning effort" />
                            </SelectTrigger>
                            <SelectContent>
                              {reasoningOptions.map((option) => (
                                <SelectItem key={option.value} value={option.value}>
                                  <div className="flex items-center gap-2">
                                    <span>{option.label}</span>
                                    <span className="text-xs text-muted-foreground">
                                      ({option.description})
                                    </span>
                                  </div>
                                </SelectItem>
                              ))}
                            </SelectContent>
                          </Select>
                          <p className="text-xs text-muted-foreground">
                            Higher effort = more thorough reasoning, slower response
                          </p>
                        </div>
                      )}
                    </div>
                  ) : (
                    // Non-OpenAI providers: show reasoning if model supports it
                    modelSupportsReasoning(provider, model) && (
                      <div className="space-y-2">
                        <label className="text-sm font-medium flex items-center gap-1">
                          <Brain className="h-4 w-4" />
                          Reasoning Effort
                        </label>
                        <Select
                          value={reasoningEffort || 'none'}
                          onValueChange={(value) => {
                            setReasoningEffort(value === 'none' ? null : value as ReasoningEffort)
                          }}
                        >
                          <SelectTrigger>
                            <SelectValue placeholder="Choose reasoning effort" />
                          </SelectTrigger>
                          <SelectContent>
                            {reasoningOptions.map((option) => (
                              <SelectItem key={option.value} value={option.value}>
                                <div className="flex items-center gap-2">
                                  <span>{option.label}</span>
                                  <span className="text-xs text-muted-foreground">
                                    ({option.description})
                                  </span>
                                </div>
                              </SelectItem>
                            ))}
                          </SelectContent>
                        </Select>
                        <p className="text-xs text-muted-foreground">
                          Higher effort = more thorough reasoning, slower response
                        </p>
                      </div>
                    )
                  )}
                </div>
              </CardContent>
            </Card>

            {/* File Upload */}
            <Card>
              <CardHeader>
                <CardTitle>📁 Marker Data Upload</CardTitle>
              </CardHeader>
              <CardContent>
                <ErrorBoundary>
                  <FileUpload showSimpleFormat />
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
                    <>
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
                    <div className="space-y-2 mt-4">
                      <label className="text-sm font-medium">Additional Information</label>
                      <textarea
                        placeholder="e.g., This dataset contains cells from a tumor microenvironment with expected immune infiltration..."
                        value={additionalInfo}
                        onChange={(e) => setAnalysisConfig({ additionalInfo: e.target.value })}
                        className="w-full px-3 py-2 border border-input bg-background text-sm rounded-md min-h-[80px] resize-y"
                        rows={3}
                      />
                      <p className="text-xs text-muted-foreground">
                        Provide any additional context about your dataset to help improve annotation accuracy (optional)
                      </p>
                    </div>
                    </>
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
                <div className="flex gap-3">
                  <Button
                    onClick={handleStartBatchAnalysis}
                    onMouseEnter={preloadAnalysisModules}
                    onFocus={preloadAnalysisModules}
                    disabled={!canStartAnalysis || isRunning}
                    className="flex-1 bg-green-600 hover:bg-green-700"
                    size="lg"
                  >
                    <Zap className="h-5 w-5 mr-2" />
                    {isRunning ? 'Batch Analysis Running...' : 'Start Batch Analysis'}
                  </Button>
                  {isRunning && (
                    <Button
                      onClick={stopAnalysis}
                      variant="destructive"
                      size="lg"
                    >
                      <Square className="h-5 w-5 mr-2" />
                      Stop
                    </Button>
                  )}
                </div>
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
                {provider === 'openai' ? (
                  openaiIdentityVerified && (model.toLowerCase().includes('gpt-5') || model.toLowerCase().includes('gpt5')) && (
                    <div className="flex justify-between">
                      <span className="text-muted-foreground">Reasoning:</span>
                      <span className="font-medium capitalize">{reasoningEffort || 'none'}</span>
                    </div>
                  )
                ) : (
                  modelSupportsReasoning(provider, model) && (
                    <div className="flex justify-between">
                      <span className="text-muted-foreground">Reasoning:</span>
                      <span className="font-medium capitalize">{reasoningEffort || 'none'}</span>
                    </div>
                  )
                )}
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