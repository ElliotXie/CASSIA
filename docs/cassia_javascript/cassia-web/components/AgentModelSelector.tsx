'use client'

import React, { useState, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { Input } from '@/components/ui/input'
import { Settings, Zap, ChevronDown, Brain } from 'lucide-react'
import { useApiKeyStore, Provider } from '@/lib/stores/api-key-store'
import modelSettings from '../public/examples/model_settings.json'
import { modelSupportsReasoning, getDefaultReasoningEffort, getReasoningEffortOptions, ReasoningEffort } from '@/lib/config/model-presets'

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
    Object.entries(provider.models || {}).forEach(([modelKey, model]: [string, any]) => {
      const costTier = model.cost_tier
      const performance = costTier === 'high' ? 'high' : costTier === 'medium' ? 'medium' : 'low'
      const speed = costTier === 'very_low' || costTier === 'low' ? 'fast' : 'medium'
      const cost = costTier === 'high' ? 'high' : costTier === 'medium' ? 'medium' : 'low'

      options.push({
        id: model.actual_name,
        name: model.description.split(' - ')[0] || model.actual_name,
        provider: providerKey,
        performance: performance as 'high' | 'medium' | 'low',
        speed: speed as 'fast' | 'medium' | 'slow',
        cost: cost as 'low' | 'medium' | 'high'
      })
    })
  })

  return options
}

const modelOptions: ModelOption[] = generateModelOptions()

// Provider configurations
const PROVIDERS = [
  { id: 'openrouter' as const, name: 'OpenRouter', description: 'Access to multiple models' },
  { id: 'openai' as const, name: 'OpenAI', description: 'Direct OpenAI API' },
  { id: 'anthropic' as const, name: 'Anthropic', description: 'Direct Anthropic API' },
  { id: 'custom' as const, name: 'Custom', description: 'Custom API endpoint' },
]

// Custom provider presets
const CUSTOM_PROVIDER_PRESETS = {
  deepseek: {
    name: 'DeepSeek',
    baseUrl: 'https://api.deepseek.com',
    models: ['deepseek-chat', 'deepseek-reasoner'],
  },
  qwen: {
    name: 'Qwen (Alibaba)',
    baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1',
    models: ['qwen-max', 'qwen-plus', 'qwen-turbo'],
  },
  kimi: {
    name: 'Kimi (Moonshot)',
    baseUrl: 'https://api.moonshot.cn/v1',
    models: ['moonshot-v1-8k', 'moonshot-v1-32k', 'moonshot-v1-128k'],
  },
  siliconflow: {
    name: 'SiliconFlow',
    baseUrl: 'https://api.siliconflow.cn/v1',
    models: ['Pro/deepseek-ai/DeepSeek-V3', 'deepseek-ai/DeepSeek-V3', 'Qwen/Qwen2.5-72B-Instruct'],
  },
  minimax: {
    name: 'MiniMax',
    baseUrl: 'https://api.minimax.io/v1',
    models: ['MiniMax-M2'],
  },
  zhipuai: {
    name: 'Zhipu AI (智谱)',
    baseUrl: 'https://open.bigmodel.cn/api/paas/v4',
    models: ['glm-4-plus', 'glm-4-flash', 'glm-4-long'],
  },
  manual: {
    name: 'Manual Entry',
    baseUrl: '',
    models: [],
  }
} as const

type CustomPresetKey = keyof typeof CUSTOM_PROVIDER_PRESETS

// Presets per provider
function getPresets(provider: Provider) {
  if (provider === 'openai') {
    return {
      performance: { model: 'gpt-5.2', name: 'Performance', icon: <Zap className="h-4 w-4" /> },
      balanced: { model: 'gpt-4o', name: 'Balanced', icon: <Settings className="h-4 w-4" /> }
    }
  } else if (provider === 'anthropic') {
    return {
      performance: { model: 'claude-sonnet-4.5', name: 'Performance', icon: <Zap className="h-4 w-4" /> },
      balanced: { model: 'claude-haiku-4.5', name: 'Balanced', icon: <Settings className="h-4 w-4" /> }
    }
  } else if (provider === 'custom') {
    return {
      performance: { model: 'deepseek-chat', name: 'Performance', icon: <Zap className="h-4 w-4" /> },
      balanced: { model: 'deepseek-chat', name: 'Balanced', icon: <Settings className="h-4 w-4" /> }
    }
  } else {
    // OpenRouter
    return {
      performance: { model: 'anthropic/claude-sonnet-4.5', name: 'Performance', icon: <Zap className="h-4 w-4" /> },
      balanced: { model: 'google/gemini-3-flash-preview', name: 'Balanced', icon: <Settings className="h-4 w-4" /> }
    }
  }
}

interface AgentModelSelectorProps {
  provider: Provider
  model: string
  onProviderChange: (provider: Provider) => void
  onModelChange: (model: string) => void
  customBaseUrl?: string
  onCustomBaseUrlChange?: (url: string) => void
  title?: string
  reasoningEffort?: ReasoningEffort | null
  onReasoningEffortChange?: (effort: ReasoningEffort | null) => void
}

export function AgentModelSelector({
  provider,
  model,
  onProviderChange,
  onModelChange,
  customBaseUrl = '',
  onCustomBaseUrlChange,
  title = "Model Configuration",
  reasoningEffort,
  onReasoningEffortChange
}: AgentModelSelectorProps) {
  const [selectedPreset, setSelectedPreset] = useState<string>('')
  const [customPreset, setCustomPreset] = useState<CustomPresetKey>('deepseek')
  const [showPresetDropdown, setShowPresetDropdown] = useState(false)

  const presets = getPresets(provider)

  const getModelsByProvider = (providerId: string) => {
    return modelOptions.filter(m => m.provider === providerId)
  }

  const getModelInfo = (modelId: string) => {
    return modelOptions.find(m => m.id === modelId)
  }

  // Handle model change and auto-set default reasoning effort
  const handleModelChange = (newModel: string) => {
    onModelChange(newModel)
    // Auto-set default reasoning effort for new model
    if (onReasoningEffortChange) {
      const defaultEffort = getDefaultReasoningEffort(provider, newModel)
      onReasoningEffortChange(defaultEffort)
    }
  }

  const applyPreset = (presetKey: 'performance' | 'balanced') => {
    const preset = presets[presetKey]
    handleModelChange(preset.model)
    setSelectedPreset(presetKey)
  }

  // When provider changes, select default model for that provider
  useEffect(() => {
    if (provider === 'custom') {
      // Set default custom preset
      const defaultPreset = CUSTOM_PROVIDER_PRESETS.deepseek
      if (onCustomBaseUrlChange) {
        onCustomBaseUrlChange(defaultPreset.baseUrl)
      }
      const newModel = defaultPreset.models[0]
      onModelChange(newModel)
      setCustomPreset('deepseek')
      // Set default reasoning effort for custom provider
      if (onReasoningEffortChange) {
        onReasoningEffortChange(getDefaultReasoningEffort('custom', newModel))
      }
    } else {
      const availableModels = getModelsByProvider(provider)
      if (availableModels.length > 0 && !availableModels.find(m => m.id === model)) {
        const newModel = availableModels[0].id
        onModelChange(newModel)
        // Set default reasoning effort for new model
        if (onReasoningEffortChange) {
          onReasoningEffortChange(getDefaultReasoningEffort(provider, newModel))
        }
      }
    }
    setSelectedPreset('')
  }, [provider])

  const PerformanceBadge = ({ performance, speed, cost }: { performance: string, speed: string, cost: string }) => {
    const qualityLabel = performance === 'low' ? 'medium' : performance

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

  const currentModel = getModelInfo(model)

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          ⚙️ {title}
        </CardTitle>
        {provider !== 'custom' && (
          <div className="flex gap-2 flex-wrap">
            {Object.entries(presets).map(([key, preset]) => (
              <Button
                key={key}
                variant={selectedPreset === key ? 'default' : 'outline'}
                size="sm"
                onClick={() => applyPreset(key as 'performance' | 'balanced')}
                className="flex items-center gap-2"
              >
                {preset.icon}
                {preset.name}
              </Button>
            ))}
          </div>
        )}
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Provider Selection */}
        <div className="space-y-2">
          <label className="text-sm font-medium">Provider</label>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
            {PROVIDERS.map((p) => (
              <button
                key={p.id}
                onClick={() => onProviderChange(p.id)}
                className={`p-3 text-left border rounded-lg transition-colors min-w-0 ${
                  provider === p.id
                    ? 'border-primary bg-primary/5'
                    : 'border-muted hover:border-primary/50'
                }`}
              >
                <div className="font-medium text-sm truncate">{p.name}</div>
                <div className="text-xs text-muted-foreground line-clamp-2">{p.description}</div>
              </button>
            ))}
          </div>
        </div>

        {/* Custom Provider Configuration */}
        {provider === 'custom' && (
          <div className="space-y-4 p-4 border rounded-lg bg-muted/30">
            {/* Preset Selection */}
            <div className="space-y-2">
              <label className="text-sm font-medium">Select Provider</label>
              <div className="relative">
                <button
                  onClick={() => setShowPresetDropdown(!showPresetDropdown)}
                  className="w-full p-3 text-left border rounded-lg bg-background flex items-center justify-between hover:border-primary/50 transition-colors"
                >
                  <span>{CUSTOM_PROVIDER_PRESETS[customPreset].name}</span>
                  <ChevronDown className={`h-4 w-4 transition-transform ${showPresetDropdown ? 'rotate-180' : ''}`} />
                </button>
                {showPresetDropdown && (
                  <div className="absolute z-10 w-full mt-1 bg-background border rounded-lg shadow-lg">
                    {Object.entries(CUSTOM_PROVIDER_PRESETS).map(([key, preset]) => (
                      <button
                        key={key}
                        onClick={() => {
                          setCustomPreset(key as CustomPresetKey)
                          if (preset.baseUrl && onCustomBaseUrlChange) {
                            onCustomBaseUrlChange(preset.baseUrl)
                            if (preset.models.length > 0) {
                              handleModelChange(preset.models[0])
                            }
                          }
                          setShowPresetDropdown(false)
                        }}
                        className={`w-full p-3 text-left hover:bg-muted transition-colors first:rounded-t-lg last:rounded-b-lg ${
                          customPreset === key ? 'bg-primary/10' : ''
                        }`}
                      >
                        <div className="font-medium text-sm">{preset.name}</div>
                        {preset.baseUrl && (
                          <div className="text-xs text-muted-foreground truncate">{preset.baseUrl}</div>
                        )}
                      </button>
                    ))}
                  </div>
                )}
              </div>
            </div>

            {/* Base URL */}
            <div className="space-y-2">
              <label className="text-sm font-medium">Base URL</label>
              {customPreset === 'manual' ? (
                <Input
                  type="url"
                  placeholder="https://api.example.com/v1"
                  value={customBaseUrl}
                  onChange={(e) => onCustomBaseUrlChange?.(e.target.value)}
                />
              ) : (
                <Input
                  type="url"
                  value={CUSTOM_PROVIDER_PRESETS[customPreset].baseUrl}
                  readOnly
                  className="bg-muted/50"
                />
              )}
            </div>

            {/* Model Selection for Custom */}
            <div className="space-y-2">
              <label className="text-sm font-medium">Model</label>
              {customPreset !== 'manual' && CUSTOM_PROVIDER_PRESETS[customPreset].models.length > 0 ? (
                <div className="space-y-2">
                  <div className="flex flex-wrap gap-2">
                    {CUSTOM_PROVIDER_PRESETS[customPreset].models.map((modelName) => (
                      <button
                        key={modelName}
                        onClick={() => handleModelChange(modelName)}
                        className={`px-3 py-1.5 text-sm border rounded-md transition-colors ${
                          model === modelName
                            ? 'border-primary bg-primary/10 text-primary'
                            : 'border-muted hover:border-primary/50'
                        }`}
                      >
                        {modelName}
                      </button>
                    ))}
                  </div>
                  <Input
                    type="text"
                    placeholder="Or enter a custom model name..."
                    value={CUSTOM_PROVIDER_PRESETS[customPreset].models.includes(model as any) ? '' : model}
                    onChange={(e) => handleModelChange(e.target.value)}
                    className="mt-2"
                  />
                </div>
              ) : (
                <Input
                  type="text"
                  placeholder="Enter model name (e.g., gpt-4, llama-3-70b)"
                  value={model}
                  onChange={(e) => handleModelChange(e.target.value)}
                />
              )}
            </div>
          </div>
        )}

        {/* Model Selection for non-custom providers */}
        {provider !== 'custom' && (
          <div className="space-y-2">
            <label className="text-sm font-medium">Model</label>
            <Select value={model} onValueChange={handleModelChange}>
              <SelectTrigger>
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                {getModelsByProvider(provider).map((m) => (
                  <SelectItem key={m.id} value={m.id}>
                    <div className="font-medium">{m.name}</div>
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>
        )}

        {/* Performance Badge */}
        {currentModel && provider !== 'custom' && (
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

        {/* Reasoning Effort Selector - shown when model supports reasoning */}
        {modelSupportsReasoning(provider, model) && onReasoningEffortChange && (
          <div className="space-y-2">
            <label className="text-sm font-medium flex items-center gap-1">
              <Brain className="h-3 w-3" />
              Reasoning Effort
            </label>
            <Select
              value={reasoningEffort || 'none'}
              onValueChange={(value) => {
                onReasoningEffortChange(value === 'none' ? null : value as ReasoningEffort)
              }}
            >
              <SelectTrigger>
                <SelectValue />
              </SelectTrigger>
              <SelectContent>
                {getReasoningEffortOptions().map((option) => (
                  <SelectItem key={option.value} value={option.value}>
                    <div className="flex flex-col">
                      <span className="font-medium">{option.label}</span>
                      <span className="text-xs text-muted-foreground">{option.description}</span>
                    </div>
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
          </div>
        )}
      </CardContent>
    </Card>
  )
}
