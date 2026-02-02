'use client'

import React, { useState, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { Settings, Zap, Clock, Brain } from 'lucide-react'
import { useConfigStore } from '@/lib/stores/config-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store'
import { modelSupportsReasoning, getReasoningEffortOptions, ReasoningEffort } from '@/lib/config/model-presets'
import modelSettings from '../public/examples/model_settings.json'

interface ModelOption {
  id: string
  name: string
  provider: string
  performance: 'high' | 'medium' | 'low'
  speed: 'fast' | 'medium' | 'slow'
  cost: 'low' | 'medium' | 'high'
}

// Generate model options from model_settings.json
function generateModelOptions(): ModelOption[] {
  const options: ModelOption[] = []
  
  Object.entries(modelSettings.providers).forEach(([providerKey, provider]) => {
    Object.entries(provider.models || {}).forEach(([modelKey, model]) => {
      const costTier = model.cost_tier
      const performance = costTier === 'high' ? 'high' : costTier === 'medium' ? 'medium' : 'low'
      const speed = costTier === 'very_low' || costTier === 'low' ? 'fast' : 'medium'
      const cost = costTier === 'high' ? 'high' : costTier === 'medium' ? 'medium' : 'low'
      
      options.push({
        id: model.actual_name,
        name: model.description.split(' - ')[0] || model.actual_name,
        provider: providerKey,
        performance,
        speed,
        cost
      })
    })
  })
  
  return options
}

const modelOptions: ModelOption[] = generateModelOptions()

const providers = [
  { id: 'openrouter', name: 'OpenRouter', description: 'Access to multiple models' },
  { id: 'openai', name: 'OpenAI', description: 'Direct OpenAI API' },
  { id: 'anthropic', name: 'Anthropic', description: 'Direct Anthropic API' },
]

const pipelineSteps = [
  { 
    id: 'annotation' as const, 
    name: 'Initial Annotation', 
    description: 'Identify cell types from marker genes',
    icon: '🔍',
    recommendedModels: [modelSettings.use_case_recommendations.annotation.best, ...modelSettings.use_case_recommendations.annotation.alternatives.slice(0, 1)]
  },
  { 
    id: 'scoring' as const, 
    name: 'Quality Scoring', 
    description: 'Evaluate annotation quality and confidence',
    icon: '📊',
    recommendedModels: [modelSettings.use_case_recommendations.scoring.best, ...modelSettings.use_case_recommendations.scoring.alternatives.slice(0, 1)]
  },
  { 
    id: 'annotationBoost' as const, 
    name: 'Annotation Boost', 
    description: 'Improve low-scoring annotations with advanced analysis',
    icon: '🚀',
    recommendedModels: [modelSettings.use_case_recommendations.annotation_boost.best, ...modelSettings.use_case_recommendations.annotation_boost.alternatives.slice(0, 1)]
  },
]

interface ModelSelectionMatrixProps {
  lockedProvider?: 'openai' | 'anthropic' | null  // null = OpenRouter mode (full control)
}

// Generate presets dynamically based on whether provider is locked
function getPresets(lockedProvider: 'openai' | 'anthropic' | null) {
  if (lockedProvider === 'openai') {
    return {
      performance: {
        name: 'Performance',
        description: 'Best quality results',
        icon: <Zap className="h-4 w-4" />,
        models: {
          annotation: { provider: 'openai', model: 'gpt-5.2' },
          scoring: { provider: 'openai', model: 'gpt-5.2' },
          annotationBoost: { provider: 'openai', model: 'gpt-5.2' }
        }
      },
      balanced: {
        name: 'Balanced',
        description: 'Good quality with reasonable speed',
        icon: <Settings className="h-4 w-4" />,
        models: {
          annotation: { provider: 'openai', model: 'gpt-4o' },
          scoring: { provider: 'openai', model: 'gpt-5-mini' },
          annotationBoost: { provider: 'openai', model: 'gpt-4o' }
        }
      }
    }
  } else if (lockedProvider === 'anthropic') {
    return {
      performance: {
        name: 'Performance',
        description: 'Best quality results',
        icon: <Zap className="h-4 w-4" />,
        models: {
          annotation: { provider: 'anthropic', model: 'claude-sonnet-4.5' },
          scoring: { provider: 'anthropic', model: 'claude-sonnet-4.5' },
          annotationBoost: { provider: 'anthropic', model: 'claude-sonnet-4.5' }
        }
      },
      balanced: {
        name: 'Balanced',
        description: 'Good quality with reasonable speed',
        icon: <Settings className="h-4 w-4" />,
        models: {
          annotation: { provider: 'anthropic', model: 'claude-haiku-4.5' },
          scoring: { provider: 'anthropic', model: 'claude-haiku-4.5' },
          annotationBoost: { provider: 'anthropic', model: 'claude-haiku-4.5' }
        }
      }
    }
  } else {
    // OpenRouter mode - full control
    return {
      performance: {
        name: 'Performance',
        description: 'Best quality results',
        icon: <Zap className="h-4 w-4" />,
        models: {
          annotation: { provider: 'openrouter', model: 'anthropic/claude-sonnet-4.5' },
          scoring: { provider: 'openrouter', model: 'openai/gpt-5.2' },
          annotationBoost: { provider: 'openrouter', model: 'anthropic/claude-sonnet-4.5' }
        }
      },
      balanced: {
        name: 'Balanced',
        description: 'Good quality with reasonable speed',
        icon: <Settings className="h-4 w-4" />,
        models: {
          annotation: { provider: 'openrouter', model: 'anthropic/claude-haiku-4.5' },
          scoring: { provider: 'openrouter', model: 'google/gemini-3-flash-preview' },
          annotationBoost: { provider: 'openrouter', model: 'anthropic/claude-haiku-4.5' }
        }
      }
    }
  }
}

export function ModelSelectionMatrix({ lockedProvider = null }: ModelSelectionMatrixProps) {
  const { pipelineModels, setPipelineModel, setReasoningEffort } = useConfigStore()
  const { provider } = useApiKeyStore()
  const reasoningOptions = getReasoningEffortOptions()
  const [selectedPreset, setSelectedPreset] = useState<string>('')
  const [hasAutoSelectedBalanced, setHasAutoSelectedBalanced] = useState(false)
  const [lastLockedProvider, setLastLockedProvider] = useState<string | null>(null)

  // Get presets based on locked provider
  const presets = getPresets(lockedProvider)

  const getModelsByProvider = (providerId: string) => {
    // If provider is locked, only show models from that provider
    const effectiveProvider = lockedProvider || providerId
    return modelOptions.filter(model => model.provider === effectiveProvider)
  }

  const getModelInfo = (modelId: string) => {
    return modelOptions.find(model => model.id === modelId)
  }

  const applyPreset = (presetKey: keyof typeof presets) => {
    const preset = presets[presetKey]
    Object.entries(preset.models).forEach(([step, config]) => {
      setPipelineModel(step as keyof typeof preset.models, config.provider, config.model)
    })
    setSelectedPreset(presetKey)
  }

  // Auto-select balanced preset when OpenRouter is chosen (only once)
  useEffect(() => {
    if (provider === 'openrouter' && !hasAutoSelectedBalanced && selectedPreset === '' && !lockedProvider) {
      applyPreset('performance')
      setHasAutoSelectedBalanced(true)
    }
  }, [provider, hasAutoSelectedBalanced, selectedPreset, lockedProvider])

  // Auto-apply balanced preset when lockedProvider changes
  useEffect(() => {
    if (lockedProvider && lockedProvider !== lastLockedProvider) {
      applyPreset('performance')
      setLastLockedProvider(lockedProvider)
    }
  }, [lockedProvider, lastLockedProvider])

  const PerformanceBadge = ({ performance, speed, cost }: { performance: string, speed: string, cost: string }) => {
    // Map performance levels to better user-friendly labels
    const qualityLabel = performance === 'low' ? 'medium' : performance;
    
    return (
      <div className="flex gap-1 flex-wrap">
        <Badge variant={performance === 'high' ? 'default' : performance === 'medium' ? 'secondary' : 'outline'} className="text-xs">
          {qualityLabel} quality
        </Badge>
        <Badge variant={speed === 'fast' ? 'default' : speed === 'medium' ? 'secondary' : 'outline'} className="text-xs">
          {speed}
        </Badge>
        <Badge variant={cost === 'low' ? 'default' : cost === 'medium' ? 'secondary' : 'outline'} className="text-xs">
          ${cost} cost
        </Badge>
      </div>
    )
  }

  // Get provider display name
  const getProviderDisplayName = (providerId: string) => {
    const providerInfo = providers.find(p => p.id === providerId)
    return providerInfo?.name || providerId
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          ⚙️ Model Configuration
          {lockedProvider && (
            <Badge variant="secondary" className="ml-2">
              Using: {getProviderDisplayName(lockedProvider)}
            </Badge>
          )}
        </CardTitle>
        <div className="flex gap-2 flex-wrap">
          {Object.entries(presets).map(([key, preset]) => (
            <Button
              key={key}
              variant={selectedPreset === key ? 'default' : 'outline'}
              size="sm"
              onClick={() => applyPreset(key as keyof typeof presets)}
              className="flex items-center gap-2"
            >
              {preset.icon}
              {preset.name}
            </Button>
          ))}
        </div>
      </CardHeader>
      <CardContent className="space-y-6">
        {pipelineSteps.map((step) => {
          const currentConfig = pipelineModels[step.id]
          const currentModel = getModelInfo(currentConfig.model)
          
          return (
            <div key={step.id} className="border rounded-lg p-4 space-y-4">
              <div className="flex items-start gap-3">
                <div className="text-2xl">{step.icon}</div>
                <div className="flex-1">
                  <h3 className="font-semibold text-lg">{step.name}</h3>
                  <p className="text-sm text-muted-foreground">{step.description}</p>
                </div>
              </div>

              <div className={`grid grid-cols-1 ${lockedProvider ? 'md:grid-cols-2' : 'md:grid-cols-3'} gap-4`}>
                {/* Provider Selector - only show when provider is not locked */}
                {!lockedProvider && (
                  <div className="space-y-2">
                    <label className="text-sm font-medium">Provider</label>
                    <Select
                      value={currentConfig.provider}
                      onValueChange={(provider) => {
                        // When provider changes, select first available model for that provider
                        const availableModels = getModelsByProvider(provider)
                        if (availableModels.length > 0) {
                          setPipelineModel(step.id, provider, availableModels[0].id)
                        }
                        setSelectedPreset('') // Clear preset selection
                      }}
                    >
                      <SelectTrigger>
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        {providers.map((provider) => (
                          <SelectItem key={provider.id} value={provider.id}>
                            <div>
                              <div className="font-medium">{provider.name}</div>
                              <div className="text-xs text-muted-foreground">{provider.description}</div>
                            </div>
                          </SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                  </div>
                )}

                <div className="space-y-2">
                  <label className="text-sm font-medium">Model</label>
                  <Select
                    value={currentConfig.model}
                    onValueChange={(model) => {
                      // Use locked provider if set, otherwise use current config's provider
                      const effectiveProvider = lockedProvider || currentConfig.provider
                      setPipelineModel(step.id, effectiveProvider, model)
                      setSelectedPreset('') // Clear preset selection
                    }}
                  >
                    <SelectTrigger>
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      {getModelsByProvider(lockedProvider || currentConfig.provider).map((model) => (
                        <SelectItem key={model.id} value={model.id}>
                          <div className="font-medium">{model.name}</div>
                        </SelectItem>
                      ))}
                    </SelectContent>
                  </Select>
                </div>

                {/* Reasoning Effort Selector - shown when model supports reasoning */}
                {modelSupportsReasoning(lockedProvider || currentConfig.provider, currentConfig.model) && (
                  <div className="space-y-2">
                    <label className="text-sm font-medium flex items-center gap-1">
                      <Brain className="h-3 w-3" />
                      Reasoning Effort
                    </label>
                    <Select
                      value={currentConfig.reasoningEffort || 'none'}
                      onValueChange={(effort) => {
                        setReasoningEffort(step.id, effort === 'none' ? null : effort as ReasoningEffort)
                      }}
                    >
                      <SelectTrigger>
                        <SelectValue />
                      </SelectTrigger>
                      <SelectContent>
                        {reasoningOptions.map((option) => (
                          <SelectItem key={option.value} value={option.value}>
                            <div>
                              <div className="font-medium">{option.label}</div>
                              <div className="text-xs text-muted-foreground">{option.description}</div>
                            </div>
                          </SelectItem>
                        ))}
                      </SelectContent>
                    </Select>
                  </div>
                )}
              </div>

              {currentModel && (
                <div className="bg-muted/50 rounded-lg p-3">
                  <div className="flex justify-between items-center">
                    <span className="text-sm font-medium">Selected: {currentModel.name}</span>
                    <PerformanceBadge 
                      performance={currentModel.performance} 
                      speed={currentModel.speed} 
                      cost={currentModel.cost} 
                    />
                  </div>
                </div>
              )}
            </div>
          )
        })}
      </CardContent>
    </Card>
  )
}