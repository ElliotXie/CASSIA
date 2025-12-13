'use client'

import React, { useState, useEffect } from 'react'
import { Eye, EyeOff, Key, ExternalLink, CheckCircle, AlertCircle, Download, Loader2, ChevronDown, Zap } from 'lucide-react'
import { testApiKey } from '@/lib/cassia/llm_utils'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { useApiKeyStore } from '@/lib/stores/api-key-store'
import { useConfigStore } from '@/lib/stores/config-store'
import { useAuthStore } from '@/lib/stores/auth-store'
import modelSettings from '../public/examples/model_settings.json'

// Custom provider presets
const CUSTOM_PROVIDER_PRESETS = {
  deepseek: {
    name: 'DeepSeek',
    baseUrl: 'https://api.deepseek.com',
    models: ['deepseek-chat', 'deepseek-reasoner'],
    helpUrl: 'https://platform.deepseek.com/api_keys'
  },
  qwen: {
    name: 'Qwen (Alibaba)',
    baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1',
    models: ['qwen-max', 'qwen-plus', 'qwen-turbo'],
    helpUrl: 'https://dashscope.console.aliyun.com/apiKey'
  },
  kimi: {
    name: 'Kimi (Moonshot)',
    baseUrl: 'https://api.moonshot.cn/v1',
    models: ['kimi-k2-turbo-preview'],
    helpUrl: 'https://platform.moonshot.cn/console/api-keys'
  },
  siliconflow: {
    name: 'SiliconFlow',
    baseUrl: 'https://api.siliconflow.cn/v1',
    models: ['deepseek-ai/DeepSeek-V3.2'],
    helpUrl: 'https://cloud.siliconflow.cn/account/ak'
  },
  minimax: {
    name: 'MiniMax',
    baseUrl: 'https://api.minimax.chat/v1',
    models: ['MiniMax-M2'],
    helpUrl: 'https://platform.minimaxi.com/user-center/basic-information/interface-key'
  },
  zhipuai: {
    name: 'Zhipu AI (智谱)',
    baseUrl: 'https://open.bigmodel.cn/api/paas/v4',
    models: ['glm-4.6'],
    helpUrl: 'https://open.bigmodel.cn/usercenter/apikeys'
  },
  manual: {
    name: 'Manual Entry',
    baseUrl: '',
    models: [],
    helpUrl: '#'
  }
} as const

type CustomPresetKey = keyof typeof CUSTOM_PROVIDER_PRESETS

// Generate provider configurations from model_settings.json
function generateProviders() {
  const providerHelp = {
    openrouter: 'https://openrouter.ai/keys',
    anthropic: 'https://console.anthropic.com/settings/keys',
    openai: 'https://platform.openai.com/api-keys',
    custom: '#'
  }

  const providers = Object.entries(modelSettings.providers).map(([providerId, provider]) => ({
    id: providerId as const,
    name: provider.name,
    description: providerId === 'openrouter'
      ? 'Access to multiple models (Recommended)'
      : providerId === 'custom'
        ? 'Any OpenAI-compatible API'
        : `${provider.name} models directly`,
    defaultModel: provider.default_model,
    helpUrl: providerHelp[providerId as keyof typeof providerHelp] || '#',
    models: Object.values(provider.models || {}).map(model => model.actual_name)
  }))

  // Ensure OpenRouter is first, custom is last
  return providers.sort((a, b) => {
    if (a.id === 'openrouter') return -1;
    if (b.id === 'openrouter') return 1;
    if (a.id === 'custom') return 1;
    if (b.id === 'custom') return -1;
    return 0;
  });
}

const PROVIDERS = generateProviders()

export function ApiKeyInput() {
  const [showKey, setShowKey] = useState(false)
  const [isValidating, setIsValidating] = useState(false)
  const [validationStatus, setValidationStatus] = useState<'idle' | 'valid' | 'invalid'>('idle')
  const [loadStatus, setLoadStatus] = useState<'idle' | 'loading' | 'success' | 'error'>('idle')
  const [errorMessage, setErrorMessage] = useState<string>('')
  const [customPreset, setCustomPreset] = useState<CustomPresetKey>('deepseek')
  const [showPresetDropdown, setShowPresetDropdown] = useState(false)
  const [isTestingApi, setIsTestingApi] = useState(false)
  const [testStatus, setTestStatus] = useState<'idle' | 'success' | 'error'>('idle')
  const [testErrorMessage, setTestErrorMessage] = useState<string>('')
  
  // Get API key data directly from the API key store for proper reactivity
  const {
    provider,
    model,
    setApiKey,
    setProvider,
    setModel: setApiModel,
    getApiKey,
    customBaseUrl,
    setCustomBaseUrl
  } = useApiKeyStore()

  // Get auth data for Load API Keys functionality
  const { isAuthenticated, user } = useAuthStore()

  // Get pipeline model configuration functions
  const { setPipelineModel } = useConfigStore()

  // Get the API key for the current provider
  const apiKey = getApiKey(provider)

  const currentProvider = PROVIDERS.find(p => p.id === provider) || PROVIDERS[0]

  // Sync pipeline models with API model when non-OpenRouter provider is used
  useEffect(() => {
    if (provider !== 'openrouter' && model) {
      // Update all pipeline models to use the same model and provider
      setPipelineModel('annotation', provider, model)
      setPipelineModel('scoring', provider, model)
      setPipelineModel('annotationBoost', provider, model)
    }
  }, [provider, model, setPipelineModel])

  const handleProviderChange = (newProvider: typeof provider) => {
    setProvider(newProvider)
    const providerConfig = PROVIDERS.find(p => p.id === newProvider)
    if (providerConfig) {
      setApiModel(providerConfig.defaultModel)
    }
    // Initialize custom provider with default preset
    if (newProvider === 'custom') {
      const defaultPreset = CUSTOM_PROVIDER_PRESETS.deepseek
      setCustomBaseUrl(defaultPreset.baseUrl)
      setApiModel(defaultPreset.models[0])
      setCustomPreset('deepseek')
    }
    setValidationStatus('idle')
  }

  // Get the appropriate help URL for custom provider
  const getHelpUrl = () => {
    if (provider === 'custom') {
      return CUSTOM_PROVIDER_PRESETS[customPreset].helpUrl
    }
    return currentProvider.helpUrl
  }

  const handleKeyChange = (value: string) => {
    setApiKey(value, provider)
    setValidationStatus('idle')
  }

  const validateApiKey = async () => {
    if (!apiKey.trim()) {
      setValidationStatus('invalid')
      return
    }

    setIsValidating(true)
    
    try {
      // Simple validation - just check if the key has a reasonable format
      // Real validation would require making an API call
      const keyPattern: Record<string, RegExp | null> = {
        openrouter: /^sk-or-v1-[a-f0-9]{64}$/i,
        anthropic: /^sk-ant-api03-[a-zA-Z0-9_-]+$/,
        openai: /^sk-[a-zA-Z0-9]{48,}$/,
        custom: null  // Custom provider accepts any key format
      }

      const pattern = keyPattern[provider]
      if (pattern && pattern.test(apiKey)) {
        setValidationStatus('valid')
      } else {
        // For now, just accept any non-empty key
        // In a real app, you'd make a test API call
        setValidationStatus(apiKey.length > 10 ? 'valid' : 'invalid')
      }
    } catch (error) {
      setValidationStatus('invalid')
    } finally {
      setIsValidating(false)
    }
  }

  const loadApiKeys = async () => {
    if (!isAuthenticated || !user) {
      setErrorMessage('Please sign in to load API keys')
      setLoadStatus('error')
      setTimeout(() => setLoadStatus('idle'), 3000)
      return
    }

    setLoadStatus('loading')
    setErrorMessage('')

    try {
      // Import Supabase client
      const { createClient } = await import('@/utils/supabase/client')
      const supabase = createClient()

      // Query user's API keys
      const { data, error } = await supabase
        .from('user_api_keys')
        .select('*')
        .eq('user_id', user.id)

      if (error) {
        throw error
      }

      // Decrypt and populate API keys
      const hasLoadedKeys = data && data.length > 0

      if (hasLoadedKeys) {
        // Process keys and await all setApiKey calls
        const keyPromises = data.map(async (keyData) => {
          const provider = keyData.provider as 'openrouter' | 'anthropic' | 'openai' | 'custom'
          try {
            // Simple base64 decryption (same as in the Supabase store)
            const decryptedKey = atob(keyData.encrypted_key)
            await setApiKey(decryptedKey, provider)
          } catch (decryptError) {
            console.error(`Failed to decrypt key for ${provider}:`, decryptError)
          }
        })
        await Promise.all(keyPromises)

        setLoadStatus('success')
        setTimeout(() => setLoadStatus('idle'), 3000)
      } else {
        setErrorMessage('No API keys found in your account')
        setLoadStatus('error')
        setTimeout(() => setLoadStatus('idle'), 5000)
      }
    } catch (error) {
      console.error('Failed to load API keys:', error)
      setErrorMessage(error instanceof Error ? error.message : 'Failed to load API keys')
      setLoadStatus('error')
      setTimeout(() => setLoadStatus('idle'), 5000)
    }
  }

  const handleTestApiKey = async () => {
    if (!apiKey.trim()) {
      setTestStatus('error')
      setTestErrorMessage('Please enter an API key first')
      setTimeout(() => setTestStatus('idle'), 3000)
      return
    }

    setIsTestingApi(true)
    setTestStatus('idle')
    setTestErrorMessage('')

    try {
      // For custom provider, use the custom base URL
      const baseUrl = provider === 'custom' ? customBaseUrl : null
      const result = await testApiKey(provider, apiKey, baseUrl)

      if (result.success) {
        setTestStatus('success')
        setTimeout(() => setTestStatus('idle'), 5000)
      } else {
        setTestStatus('error')
        setTestErrorMessage(result.error || 'API test failed')
        setTimeout(() => setTestStatus('idle'), 5000)
      }
    } catch (error) {
      setTestStatus('error')
      setTestErrorMessage(error instanceof Error ? error.message : 'API test failed')
      setTimeout(() => setTestStatus('idle'), 5000)
    } finally {
      setIsTestingApi(false)
    }
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center space-x-2">
          <Key className="h-5 w-5" />
          <span>API Configuration</span>
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Provider Selection */}
        <div className="space-y-2">
          <label className="text-sm font-medium">Provider</label>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-2">
            {PROVIDERS.map((p) => (
              <button
                key={p.id}
                onClick={() => handleProviderChange(p.id)}
                className={`p-3 text-left border rounded-lg transition-colors ${
                  provider === p.id
                    ? 'border-primary bg-primary/5'
                    : 'border-muted hover:border-primary/50'
                }`}
              >
                <div className="font-medium text-sm">{p.name}</div>
                <div className="text-xs text-muted-foreground">{p.description}</div>
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
                          if (preset.baseUrl) {
                            setCustomBaseUrl(preset.baseUrl)
                            // Set first model as default
                            if (preset.models.length > 0) {
                              setApiModel(preset.models[0])
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

            {/* Base URL - editable for manual, readonly for presets */}
            <div className="space-y-2">
              <label className="text-sm font-medium">Base URL</label>
              {customPreset === 'manual' ? (
                <Input
                  type="url"
                  placeholder="https://api.example.com/v1"
                  value={customBaseUrl}
                  onChange={(e) => setCustomBaseUrl(e.target.value)}
                />
              ) : (
                <Input
                  type="url"
                  value={CUSTOM_PROVIDER_PRESETS[customPreset].baseUrl}
                  readOnly
                  className="bg-muted/50"
                />
              )}
              <p className="text-xs text-muted-foreground">
                {customPreset === 'manual'
                  ? 'Enter any OpenAI-compatible endpoint URL'
                  : `Pre-configured for ${CUSTOM_PROVIDER_PRESETS[customPreset].name}`
                }
              </p>
            </div>
          </div>
        )}

        {/* API Key Input */}
        <div className="space-y-2">
          <div className="flex items-center justify-between">
            <label className="text-sm font-medium">API Key</label>
            <div className="flex items-center space-x-2">
              {isAuthenticated && (
                <Button
                  onClick={loadApiKeys}
                  disabled={loadStatus === 'loading'}
                  variant={loadStatus === 'success' ? 'default' : loadStatus === 'error' ? 'destructive' : 'outline'}
                  size="sm"
                  className="text-xs"
                >
                  {loadStatus === 'loading' ? (
                    <>
                      <Loader2 className="mr-1 h-3 w-3 animate-spin" />
                      Loading...
                    </>
                  ) : loadStatus === 'success' ? (
                    <>
                      <CheckCircle className="mr-1 h-3 w-3" />
                      Loaded
                    </>
                  ) : loadStatus === 'error' ? (
                    <>
                      <AlertCircle className="mr-1 h-3 w-3" />
                      Failed
                    </>
                  ) : (
                    <>
                      <Download className="mr-1 h-3 w-3" />
                      Load Keys
                    </>
                  )}
                </Button>
              )}
              <Button
                onClick={handleTestApiKey}
                disabled={isTestingApi || !apiKey}
                variant={testStatus === 'success' ? 'default' : testStatus === 'error' ? 'destructive' : 'outline'}
                size="sm"
                className="text-xs"
              >
                {isTestingApi ? (
                  <>
                    <Loader2 className="mr-1 h-3 w-3 animate-spin" />
                    Testing...
                  </>
                ) : testStatus === 'success' ? (
                  <>
                    <CheckCircle className="mr-1 h-3 w-3" />
                    API OK
                  </>
                ) : testStatus === 'error' ? (
                  <>
                    <AlertCircle className="mr-1 h-3 w-3" />
                    Failed
                  </>
                ) : (
                  <>
                    <Zap className="mr-1 h-3 w-3" />
                    Test API
                  </>
                )}
              </Button>
              <a
                href={getHelpUrl()}
                target="_blank"
                rel="noopener noreferrer"
                className="text-xs text-primary hover:underline flex items-center space-x-1"
              >
                <span>Get API key</span>
                <ExternalLink className="h-3 w-3" />
              </a>
            </div>
          </div>
          
          <div className="relative">
            <Input
              type={showKey ? 'text' : 'password'}
              placeholder={`Enter your ${currentProvider.name} API key`}
              value={apiKey}
              onChange={(e) => handleKeyChange(e.target.value)}
              className="pr-20"
            />
            <div className="absolute inset-y-0 right-0 flex items-center space-x-1 pr-2">
              {apiKey && (
                <button
                  onClick={validateApiKey}
                  disabled={isValidating}
                  className="p-1 text-muted-foreground hover:text-foreground"
                >
                  {isValidating ? (
                    <div className="w-4 h-4 border-2 border-primary border-t-transparent rounded-full animate-spin" />
                  ) : validationStatus === 'valid' ? (
                    <CheckCircle className="h-4 w-4 text-green-500" />
                  ) : validationStatus === 'invalid' ? (
                    <AlertCircle className="h-4 w-4 text-red-500" />
                  ) : null}
                </button>
              )}
              <button
                onClick={() => setShowKey(!showKey)}
                className="p-1 text-muted-foreground hover:text-foreground"
              >
                {showKey ? <EyeOff className="h-4 w-4" /> : <Eye className="h-4 w-4" />}
              </button>
            </div>
          </div>
          
          {validationStatus === 'invalid' && (
            <p className="text-xs text-red-600">
              Please enter a valid {currentProvider.name} API key
            </p>
          )}
          
          {loadStatus === 'error' && errorMessage && (
            <p className="text-xs text-red-600">
              {errorMessage}
            </p>
          )}
          
          {loadStatus === 'success' && (
            <p className="text-xs text-green-600">
              API keys loaded successfully from your account
            </p>
          )}

          {testStatus === 'error' && testErrorMessage && (
            <p className="text-xs text-red-600">
              API Test Failed: {testErrorMessage}
            </p>
          )}

          {testStatus === 'success' && (
            <p className="text-xs text-green-600">
              API key is valid and working
            </p>
          )}
        </div>

        {/* Status Indicator */}
        {apiKey && (
          <div className={`p-2 rounded-md text-xs ${
            validationStatus === 'valid'
              ? 'bg-green-50 text-green-700 border border-green-200'
              : validationStatus === 'invalid'
              ? 'bg-red-50 text-red-700 border border-red-200'
              : 'bg-blue-50 text-blue-700 border border-blue-200'
          }`}>
            {validationStatus === 'valid' && '✅ API key is configured'}
            {validationStatus === 'invalid' && '❌ API key validation failed'}
            {validationStatus === 'idle' && '💡 Click the check icon to validate your API key'}
          </div>
        )}
      </CardContent>
    </Card>
  )
}