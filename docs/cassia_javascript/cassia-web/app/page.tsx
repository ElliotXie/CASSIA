'use client'

import Link from 'next/link'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Zap, Settings, Download, HelpCircle, Key, CheckCircle, XCircle, Loader2, BookOpen, FileText, Cpu, Github, FileX, Eye, EyeOff, ChevronDown, ExternalLink } from 'lucide-react'
import { useState, Suspense } from 'react'
import { useApiKeyStore } from '@/lib/stores/api-key-store-simple'
import { useAuthStore } from '@/lib/stores/auth-store'
import { ContactDialog } from '@/components/ContactDialog'
import { AuthButton } from '@/components/auth/AuthButton'
import { UserDashboard } from '@/components/dashboard/UserDashboard'
import { LoadApiKeysButton } from '@/components/LoadApiKeysButton'
import modelSettings from '../public/examples/model_settings.json'
import { EmailConfirmationHandler } from '@/components/EmailConfirmationHandler'
import { testApiKey as testApiKeyFn } from '@/lib/cassia/llm_utils'

// Custom provider presets (same as ApiKeyInput.tsx)
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
    models: ['moonshot-v1-8k', 'moonshot-v1-32k', 'moonshot-v1-128k'],
    helpUrl: 'https://platform.moonshot.cn/console/api-keys'
  },
  siliconflow: {
    name: 'SiliconFlow',
    baseUrl: 'https://api.siliconflow.cn/v1',
    models: ['Pro/deepseek-ai/DeepSeek-V3', 'deepseek-ai/DeepSeek-V3', 'Qwen/Qwen2.5-72B-Instruct'],
    helpUrl: 'https://cloud.siliconflow.cn/account/ak'
  },
  minimax: {
    name: 'MiniMax',
    baseUrl: 'https://api.minimax.io/v1',
    models: ['MiniMax-M2'],
    helpUrl: 'https://platform.minimaxi.com/user-center/basic-information/interface-key'
  },
  zhipuai: {
    name: 'Zhipu AI (智谱)',
    baseUrl: 'https://open.bigmodel.cn/api/paas/v4',
    models: ['glm-4-plus', 'glm-4-flash', 'glm-4-long'],
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

export default function HomePage() {
  const [showApiKeyModal, setShowApiKeyModal] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  const [showDashboard, setShowDashboard] = useState(false)
  const { apiKeys, provider, model, setApiKey, setProvider, setModel, customBaseUrl, setCustomBaseUrl } = useApiKeyStore()
  const { isAuthenticated } = useAuthStore()
  const [tempApiKey, setTempApiKey] = useState(apiKeys[provider])
  const [tempProvider, setTempProvider] = useState(provider)
  const [tempModel, setTempModel] = useState(model)
  const [isTestingApi, setIsTestingApi] = useState(false)
  const [testResult, setTestResult] = useState<'success' | 'error' | null>(null)
  const [testMessage, setTestMessage] = useState('')
  const [showApiKey, setShowApiKey] = useState(false)
  const [isSaving, setIsSaving] = useState(false)
  const [saveSuccess, setSaveSuccess] = useState(false)
  // Custom provider state
  const [customPreset, setCustomPreset] = useState<CustomPresetKey>('deepseek')
  const [showPresetDropdown, setShowPresetDropdown] = useState(false)
  const [tempCustomBaseUrl, setTempCustomBaseUrl] = useState(customBaseUrl || CUSTOM_PROVIDER_PRESETS.deepseek.baseUrl)

  // Get current API key
  const apiKey = apiKeys[provider]

  // Default models for each provider from model_settings.json
  const getDefaultModel = (providerName: string) => {
    const providerData = modelSettings.providers[providerName as keyof typeof modelSettings.providers];
    return providerData?.default_model || modelSettings.providers.openrouter.default_model;
  }

  // Update model when provider changes
  const handleProviderChange = (newProvider: string) => {
    setTempProvider(newProvider as typeof provider)
    setTestResult(null)
    setTestMessage('')

    if (newProvider === 'custom') {
      // Initialize with DeepSeek preset (model will be set in batch/pipeline pages)
      const defaultPreset = CUSTOM_PROVIDER_PRESETS.deepseek
      setTempCustomBaseUrl(defaultPreset.baseUrl)
      setCustomPreset('deepseek')
    } else {
      setTempModel(getDefaultModel(newProvider))
    }
  }

  // Test API key function
  const testApiKey = async () => {
    if (!tempApiKey.trim()) {
      setTestResult('error')
      setTestMessage('Please enter an API key')
      return
    }

    setIsTestingApi(true)
    setTestResult(null)
    setTestMessage('')

    try {
      // Handle custom provider using the imported testApiKeyFn
      if (tempProvider === 'custom') {
        if (!tempCustomBaseUrl) {
          setTestResult('error')
          setTestMessage('Base URL is required for custom provider')
          setIsTestingApi(false)
          return
        }

        // Get test model from the selected preset
        const preset = CUSTOM_PROVIDER_PRESETS[customPreset]
        const testModel = preset?.models?.[0] || tempModel
        if (!testModel) {
          setTestResult('error')
          setTestMessage('Model is required for custom provider testing')
          setIsTestingApi(false)
          return
        }

        const result = await testApiKeyFn('custom', tempApiKey, tempCustomBaseUrl, testModel)

        if (result.success) {
          setTestResult('success')
          setTestMessage('API key is valid and working!')
        } else {
          setTestResult('error')
          setTestMessage(`API test failed: ${result.error}`)
        }
        setIsTestingApi(false)
        return
      }

      let testUrl = ''
      let testHeaders: any = {}
      let testBody: any = {}

      switch (tempProvider) {
        case 'openai':
          testUrl = 'https://api.openai.com/v1/chat/completions'
          testHeaders = {
            'Authorization': `Bearer ${tempApiKey}`,
            'Content-Type': 'application/json'
          }
          testBody = {
            model: tempModel,
            messages: [{ role: 'user', content: 'Hello' }],
            max_tokens: 5
          }
          break
        case 'anthropic':
          testUrl = 'https://api.anthropic.com/v1/messages'
          testHeaders = {
            'x-api-key': tempApiKey,
            'Content-Type': 'application/json',
            'anthropic-version': '2023-06-01'
          }
          testBody = {
            model: tempModel,
            max_tokens: 5,
            messages: [{ role: 'user', content: 'Hello' }]
          }
          break
        case 'openrouter':
          testUrl = 'https://openrouter.ai/api/v1/chat/completions'
          testHeaders = {
            'Authorization': `Bearer ${tempApiKey}`,
            'Content-Type': 'application/json',
            'HTTP-Referer': window.location.origin,
            'X-Title': 'CASSIA'
          }
          testBody = {
            model: tempModel,
            messages: [{ role: 'user', content: 'Hello' }],
            max_tokens: 5
          }
          break
      }

      const response = await fetch(testUrl, {
        method: 'POST',
        headers: testHeaders,
        body: JSON.stringify(testBody)
      })

      if (response.ok) {
        setTestResult('success')
        setTestMessage('API key is valid and working!')
      } else {
        const errorData = await response.json().catch(() => ({ error: { message: 'Unknown error' } }))
        setTestResult('error')
        setTestMessage(`API test failed: ${errorData.error?.message || response.statusText}`)
      }
    } catch (error: any) {
      setTestResult('error')
      setTestMessage(`Connection failed: ${error.message}`)
    } finally {
      setIsTestingApi(false)
    }
  }

  const handleSaveApiKey = async () => {
    setIsSaving(true)
    try {
      await setApiKey(tempApiKey, tempProvider)
      setProvider(tempProvider)
      setModel(tempModel)

      // Save custom base URL if using custom provider
      if (tempProvider === 'custom') {
        setCustomBaseUrl(tempCustomBaseUrl)
      }

      setSaveSuccess(true)
      // Show success message briefly, then close modal
      setTimeout(() => {
        setSaveSuccess(false)
        setShowApiKeyModal(false)
      }, 1500)
    } catch (error) {
      console.error('Failed to save API key:', error)
      setTestResult('error')
      setTestMessage('Failed to save API key. Please try again.')
    } finally {
      setIsSaving(false)
    }
  }
  
  const handleLoadApiKeysSuccess = () => {
    // Get the latest API keys from the store
    const store = useApiKeyStore.getState()
    const loadedKey = store.apiKeys[tempProvider as keyof typeof store.apiKeys]

    if (loadedKey) {
      setTempApiKey(loadedKey)
      setTestResult('success')
      setTestMessage('API key loaded from your account')
    }

    // Also load customBaseUrl if custom provider is selected
    if (tempProvider === 'custom' && store.customBaseUrl) {
      setTempCustomBaseUrl(store.customBaseUrl)
    }
  }

  const downloadExample = (filename: string) => {
    const link = document.createElement('a')
    link.href = `/examples/${filename}`
    link.download = filename
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
  }

  return (
    <div className="min-h-screen">
      {/* Email confirmation notification */}
      <Suspense fallback={null}>
        <EmailConfirmationHandler />
      </Suspense>

      {/* Enhanced Header with glassmorphism */}
      <header className="glass border-b border-white/20 sticky top-0 z-50">
        <div className="container mx-auto px-6 py-6">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-4">
              <div className="relative">
                <div className="w-12 h-12 bg-gradient-to-br from-blue-600 via-purple-600 to-indigo-700 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                  <span className="text-white font-bold text-lg">🧬</span>
                </div>
                <div className="absolute -top-1 -right-1 w-4 h-4 bg-green-400 rounded-full animate-pulse-slow"></div>
              </div>
              <div>
                <h1 className="text-2xl font-bold gradient-text">CASSIA</h1>
                <p className="text-sm text-gray-600 dark:text-gray-300">Collective Agent System for Single-cell Interpretable Annotation</p>
              </div>
            </div>
            <div className="flex items-center space-x-3">
              {isAuthenticated && (
                <Button 
                  variant="outline" 
                  size="sm"
                  onClick={() => setShowDashboard(true)}
                  className="glass border-white/30 hover:bg-white/20 btn-modern"
                >
                  <Settings className="h-4 w-4 mr-2" />
                  Dashboard
                </Button>
              )}
              <Button
                variant="outline"
                size="sm"
                onClick={() => {
                  setTempApiKey(apiKey)
                  setTempProvider(provider)
                  setTempModel(model || getDefaultModel(provider))
                  setTempCustomBaseUrl(customBaseUrl || CUSTOM_PROVIDER_PRESETS.deepseek.baseUrl)
                  setCustomPreset('deepseek')
                  setShowPresetDropdown(false)
                  setTestResult(null)
                  setTestMessage('')
                  setIsSaving(false)
                  setSaveSuccess(false)
                  setShowApiKeyModal(true)
                }}
                className="glass border-white/30 hover:bg-white/20 btn-modern"
              >
                <Key className="h-4 w-4 mr-2" />
                {apiKey ? 'API Key Set' : 'Set API Key'}
              </Button>
              <Button 
                variant="outline" 
                size="sm"
                onClick={() => setShowContactModal(true)}
                className="glass border-white/30 hover:bg-white/20 btn-modern"
              >
                <HelpCircle className="h-4 w-4 mr-2" />
                Help
              </Button>
              <AuthButton />
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="container mx-auto px-6 py-16">
        {/* Enhanced Hero Section */}
        <div className="text-center mb-20">
          
          <h2 className="text-7xl lg:text-8xl font-bold mb-6">
            <span className="gradient-text">CASSIA</span>
          </h2>

          <p className="text-xl text-gray-600 dark:text-gray-300 max-w-4xl mx-auto mb-10 leading-relaxed">
            <span className="font-semibold text-blue-600 dark:text-blue-400">Built for biologists, powered by cutting-edge language models.</span>
          </p>
          
          
          {/* Stats */}
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-6 max-w-4xl mx-auto">
            <div className="text-center">
              <div className="text-3xl font-bold gradient-text mb-2">15-40%</div>
              <div className="text-gray-600 dark:text-gray-300">Better than traditional methods</div>
            </div>
            <div className="text-center">
              <div className="text-3xl font-bold gradient-text mb-2">Interpretable</div>
              <div className="text-gray-600 dark:text-gray-300">Clear biological reasoning</div>
            </div>
            <div className="text-center">
              <div className="text-3xl font-bold gradient-text mb-2">&lt;2min</div>
              <div className="text-gray-600 dark:text-gray-300">Average runtime</div>
            </div>
            <div className="text-center">
              <div className="text-3xl font-bold gradient-text mb-2">1000+</div>
              <div className="text-gray-600 dark:text-gray-300">Cell types benchmarked</div>
            </div>
          </div>
          
          {/* Validation Badge */}
          <div className="mt-8 text-center">
            <p className="text-sm text-gray-500 dark:text-gray-400 italic">
              Validated in real-world single-cell studies and peer-reviewed publications
            </p>
          </div>
        </div>

        {/* Enhanced Quick Start Guide */}
        <div className="glass rounded-2xl p-8 mb-16 border border-white/20">
          <div className="text-center mb-8">
            <h3 className="text-2xl font-bold gradient-text mb-2">Quick Start Guide</h3>
            <p className="text-gray-600 dark:text-gray-300">Get started with CASSIA in just 4 simple steps</p>
          </div>
          
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
            <div className="text-center group">
              <div className="relative mx-auto mb-4 w-16 h-16">
                <div className="w-16 h-16 bg-gradient-to-r from-yellow-500 to-orange-600 rounded-full flex items-center justify-center text-white font-bold text-xl shadow-lg group-hover:scale-110 transition-transform duration-300">
                  1
                </div>
                <div className="absolute -top-2 -right-2 w-6 h-6 bg-red-400 rounded-full animate-pulse"></div>
              </div>
              <h4 className="font-semibold mb-2 text-gray-900 dark:text-white">Configure API</h4>
              <p className="text-sm text-gray-600 dark:text-gray-300">Set up your AI provider key</p>
            </div>
            
            <div className="text-center group">
              <div className="relative mx-auto mb-4 w-16 h-16">
                <div className="w-16 h-16 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full flex items-center justify-center text-white font-bold text-xl shadow-lg group-hover:scale-110 transition-transform duration-300">
                  2
                </div>
                <div className="absolute -top-2 -right-2 w-6 h-6 bg-green-400 rounded-full animate-pulse" style={{animationDelay: '0.5s'}}></div>
              </div>
              <h4 className="font-semibold mb-2 text-gray-900 dark:text-white">Prepare Data</h4>
              <p className="text-sm text-gray-600 dark:text-gray-300">FindAllMarkers CSV/XLSX file</p>
            </div>
            
            <div className="text-center group">
              <div className="relative mx-auto mb-4 w-16 h-16">
                <div className="w-16 h-16 bg-gradient-to-r from-purple-500 to-pink-600 rounded-full flex items-center justify-center text-white font-bold text-xl shadow-lg group-hover:scale-110 transition-transform duration-300">
                  3
                </div>
                <div className="absolute -top-2 -right-2 w-6 h-6 bg-blue-400 rounded-full animate-pulse" style={{animationDelay: '1s'}}></div>
              </div>
              <h4 className="font-semibold mb-2 text-gray-900 dark:text-white">Select Module</h4>
              <p className="text-sm text-gray-600 dark:text-gray-300">Start with Pipeline or Batch</p>
            </div>
            
            <div className="text-center group">
              <div className="relative mx-auto mb-4 w-16 h-16">
                <div className="w-16 h-16 bg-gradient-to-r from-green-500 to-teal-600 rounded-full flex items-center justify-center text-white font-bold text-xl shadow-lg group-hover:scale-110 transition-transform duration-300">
                  4
                </div>
                <div className="absolute -top-2 -right-2 w-6 h-6 bg-purple-400 rounded-full animate-pulse" style={{animationDelay: '1.5s'}}></div>
              </div>
              <h4 className="font-semibold mb-2 text-gray-900 dark:text-white">Optional Agents</h4>
              <p className="text-sm text-gray-600 dark:text-gray-300">Enhance results with specialized tools</p>
            </div>
          </div>
        </div>

        {/* API Key Configuration Section */}
        <div className="glass rounded-2xl p-6 mb-16 border border-white/20">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-4">
              <div className="w-12 h-12 bg-gradient-to-br from-yellow-500 to-orange-600 rounded-xl flex items-center justify-center shadow-lg">
                <Key className="h-6 w-6 text-white" />
              </div>
              <div>
                <h3 className="text-xl font-bold text-gray-900 dark:text-white flex items-center">
                  API Configuration
                  {!apiKey && <span className="ml-3 px-2.5 py-0.5 bg-red-500 text-white text-xs rounded-full">Required</span>}
                  {apiKey && <span className="ml-3 px-2.5 py-0.5 bg-green-500 text-white text-xs rounded-full">Configured</span>}
                </h3>
                <p className="text-gray-600 dark:text-gray-300 text-sm">
                  Supports OpenRouter, OpenAI, Anthropic, or any custom API. At least one API key is required.
                </p>
              </div>
            </div>
            <div className="flex items-center space-x-3">
              {isAuthenticated && <LoadApiKeysButton />}
              <Button
                onClick={() => {
                  setTempApiKey(apiKey)
                  setTempProvider(provider)
                  setTempModel(model || getDefaultModel(provider))
                  setTempCustomBaseUrl(customBaseUrl || CUSTOM_PROVIDER_PRESETS.deepseek.baseUrl)
                  setCustomPreset('deepseek')
                  setShowPresetDropdown(false)
                  setTestResult(null)
                  setTestMessage('')
                  setIsSaving(false)
                  setSaveSuccess(false)
                  setShowApiKeyModal(true)
                }}
                className={`px-6 py-2 font-semibold shadow-md hover:shadow-lg transition-all duration-300 btn-modern ${
                  !apiKey
                    ? 'bg-gradient-to-r from-yellow-500 to-orange-600 hover:from-yellow-600 hover:to-orange-700 text-white'
                    : 'bg-gradient-to-r from-green-600 to-emerald-600 hover:from-green-700 hover:to-emerald-700 text-white'
                }`}
              >
                <Key className="h-4 w-4 mr-2" />
                {!apiKey ? 'Configure API Key' : 'Update'}
              </Button>
            </div>
          </div>
        </div>

        {/* Enhanced Main Options */}
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8 mb-16">
          {/* Full Pipeline */}
          <Card className="relative overflow-hidden glass border border-white/20 card-hover group">
            <div className="absolute inset-0 bg-gradient-to-br from-blue-500/20 to-purple-600/20 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></div>
            <CardHeader className="pb-4 relative z-10">
              <div className="flex items-center space-x-4">
                <div className="relative">
                  <div className="w-14 h-14 bg-gradient-to-br from-blue-600 to-purple-600 rounded-xl flex items-center justify-center shadow-lg">
                    <Zap className="h-7 w-7 text-white" />
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-yellow-400 rounded-full animate-pulse"></div>
                </div>
                <div>
                  <CardTitle className="text-xl text-gray-900 dark:text-white">CASSIA Pipeline</CardTitle>
                  <CardDescription className="text-gray-600 dark:text-gray-300">🚀 Recommended starting point</CardDescription>
                </div>
              </div>
            </CardHeader>
            <CardContent className="space-y-6 relative z-10">
              <ul className="space-y-3 text-sm">
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Automated cell type annotation</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Quality scoring and validation</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Annotation merging and grouping</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Boost analysis for challenging clusters</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-blue-500 to-purple-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Interactive HTML reports</span>
                </li>
              </ul>
              <Button asChild className="w-full bg-gradient-to-r from-blue-600 to-purple-600 hover:from-blue-700 hover:to-purple-700 text-white btn-modern">
                <Link href="/pipeline">
                  Run CASSIA Pipeline
                </Link>
              </Button>
            </CardContent>
          </Card>

          {/* Batch Processing */}
          <Card className="relative overflow-hidden glass border border-white/20 card-hover group">
            <div className="absolute inset-0 bg-gradient-to-br from-green-500/20 to-emerald-600/20 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></div>
            <CardHeader className="pb-4 relative z-10">
              <div className="flex items-center space-x-4">
                <div className="relative">
                  <div className="w-14 h-14 bg-gradient-to-br from-green-600 to-emerald-600 rounded-xl flex items-center justify-center shadow-lg">
                    <Settings className="h-7 w-7 text-white" />
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-blue-400 rounded-full animate-pulse" style={{animationDelay: '0.5s'}}></div>
                </div>
                <div>
                  <CardTitle className="text-xl text-gray-900 dark:text-white">CASSIA Batch</CardTitle>
                  <CardDescription className="text-gray-600 dark:text-gray-300">🔄 Alternative starting point</CardDescription>
                </div>
              </div>
            </CardHeader>
            <CardContent className="space-y-6 relative z-10">
              <ul className="space-y-3 text-sm">
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-green-500 to-emerald-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Process multiple clusters in parallel</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-green-500 to-emerald-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Configurable worker count</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-green-500 to-emerald-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Automatic retry on failures</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-green-500 to-emerald-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Comprehensive CSV reports</span>
                </li>
                <li className="flex items-center space-x-3 group/item">
                  <span className="w-2 h-2 bg-gradient-to-r from-green-500 to-emerald-600 rounded-full group-hover/item:scale-125 transition-transform"></span>
                  <span className="text-gray-700 dark:text-gray-300">Best for large datasets</span>
                </li>
              </ul>
              <Button asChild className="w-full bg-gradient-to-r from-green-600 to-emerald-600 hover:from-green-700 hover:to-emerald-700 text-white btn-modern">
                <Link href="/batch">
                  Run CASSIA Batch
                </Link>
              </Button>
            </CardContent>
          </Card>

          {/* Individual Agents */}
          <Card className="relative overflow-hidden glass border border-white/20 card-hover group">
            <div className="absolute inset-0 bg-gradient-to-br from-purple-500/20 to-pink-600/20 opacity-0 group-hover:opacity-100 transition-opacity duration-300"></div>
            <CardHeader className="pb-4 relative z-10">
              <div className="flex items-center space-x-4">
                <div className="relative">
                  <div className="w-14 h-14 bg-gradient-to-br from-purple-600 to-pink-600 rounded-xl flex items-center justify-center shadow-lg">
                    <Settings className="h-7 w-7 text-white" />
                  </div>
                  <div className="absolute -top-1 -right-1 w-4 h-4 bg-orange-400 rounded-full animate-pulse" style={{animationDelay: '1s'}}></div>
                </div>
                <div>
                  <CardTitle className="text-xl text-gray-900 dark:text-white">Optional Agents</CardTitle>
                  <CardDescription className="text-gray-600 dark:text-gray-300">🔧 Use after Pipeline/Batch</CardDescription>
                </div>
              </div>
            </CardHeader>
            <CardContent className="space-y-6 relative z-10">
              <div className="grid grid-cols-2 gap-3 text-sm">
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/agents/annotation-boost">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Annotation Boost</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Iterative analysis</div>
                    </div>
                  </Link>
                </Button>
{/* Hidden until progress indicator is added
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/agents/symphony-compare">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Symphony Compare</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Multi-model consensus</div>
                    </div>
                  </Link>
                </Button>
                */}
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/agents/scoring">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Scoring Agent</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Quality evaluation</div>
                    </div>
                  </Link>
                </Button>
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/agents/merging">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Annotation Merging</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Group similar types</div>
                    </div>
                  </Link>
                </Button>
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Resources Section */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
          {/* Documentation & Resources */}
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center space-x-2">
                <BookOpen className="h-5 w-5" />
                <span>Documentation & Resources</span>
              </CardTitle>
              <CardDescription>
                Comprehensive guides and examples for CASSIA
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-3">
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://docs.cassia.bio/en/docs/python/introduction" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">🐍</span>
                  <span className="ml-2">Python Documentation/Vignette</span>
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://docs.cassia.bio/en/docs/r/introduction" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">📊</span>
                  <span className="ml-2">R Documentation/Vignette</span>
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://sc-llm-benchmark.com/methods/cassia" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">🤖</span>
                  <span className="ml-2">LLMs Annotation Benchmark</span>
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://doi.org/10.1038/s41467-025-67084-x" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">📄</span>
                  <span className="ml-2">Nature Communications Paper</span>
                </Link>
              </Button>
            </CardContent>
          </Card>

          {/* Installation & Support */}
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center space-x-2">
                <Download className="h-5 w-5" />
                <span>Installation & Support</span>
              </CardTitle>
              <CardDescription>
                Get CASSIA packages and community support
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-3">
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">📦</span>
                  <span className="ml-2">R & Python Package Installation</span>
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-start" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA/issues" target="_blank" rel="noopener noreferrer" className="flex items-center">
                  <span className="w-6 text-center flex-shrink-0">🐛</span>
                  <span className="ml-2">Report Issues</span>
                </Link>
              </Button>
            </CardContent>
          </Card>
        </div>
      </main>

      {/* Enhanced Footer */}
      <footer className="glass border-t border-white/20 mt-20">
        <div className="container mx-auto px-6 py-12">
          <div className="text-center">
            <div className="flex items-center justify-center space-x-2 mb-4">
              <div className="w-8 h-8 bg-gradient-to-br from-blue-600 to-purple-600 rounded-lg flex items-center justify-center">
                <span className="text-white font-bold text-sm">🧬</span>
              </div>
              <h3 className="text-xl font-bold gradient-text">CASSIA v2.0</h3>
            </div>
            <p className="text-gray-600 dark:text-gray-300 mb-6 max-w-2xl mx-auto">
              Collective Agent System for Single-cell Interpretable Annotation
              <br />
              <span className="text-sm">Built with Next.js, powered by Large Language Models</span>
            </p>
            <div className="flex flex-wrap justify-center gap-6 text-sm text-gray-500 dark:text-gray-400">
              <a href="https://docs.cassia.bio/en" target="_blank" rel="noopener noreferrer" className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors">Documentation</a>
              <a href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer" className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors">GitHub</a>
              <button onClick={() => setShowContactModal(true)} className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors">Support</button>
            </div>
            <div className="mt-6 pt-6 border-t border-white/20">
              <p className="text-xs text-gray-500 dark:text-gray-400">
                © 2024 CASSIA. All rights reserved. Built for the scientific community.
              </p>
            </div>
          </div>
        </div>
      </footer>

      {/* Enhanced API Key Modal */}
      {showApiKeyModal && (
        <div className="fixed inset-0 bg-black/80 backdrop-blur-md flex items-center justify-center z-50 p-4">
          <div className="bg-white dark:bg-gray-900 rounded-2xl p-8 max-w-lg w-full border-2 border-blue-200 dark:border-blue-800 shadow-2xl">
            <div className="text-center mb-8">
              <div className="w-20 h-20 bg-gradient-to-br from-blue-500 via-purple-500 to-indigo-600 rounded-2xl flex items-center justify-center mx-auto mb-6 shadow-xl animate-glow">
                <Key className="h-10 w-10 text-white" />
              </div>
              <h2 className="text-3xl font-bold text-gray-900 dark:text-white mb-3">
                <span className="bg-gradient-to-r from-blue-600 via-purple-600 to-indigo-600 bg-clip-text text-transparent">
                  Configure API Settings
                </span>
              </h2>
              <p className="text-gray-600 dark:text-gray-400 text-lg">Set up your AI provider and API key</p>
            </div>
            
            <div className="space-y-6">
              <div>
                <label className="block text-sm font-semibold mb-3 text-gray-900 dark:text-white">AI Provider</label>
                <select
                  value={tempProvider}
                  onChange={(e) => handleProviderChange(e.target.value)}
                  className="w-full px-4 py-3 border-2 border-blue-200 dark:border-blue-700 rounded-xl text-gray-900 dark:text-white bg-blue-50 dark:bg-blue-900/20 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all"
                >
                  <option value="openrouter">OpenRouter (Recommended)</option>
                  <option value="openai">OpenAI</option>
                  <option value="anthropic">Anthropic</option>
                  <option value="custom">Custom Provider</option>
                </select>
              </div>

              {/* Custom Provider Configuration */}
              {tempProvider === 'custom' && (
                <div className="space-y-4 p-4 border-2 border-purple-200 dark:border-purple-700 rounded-xl bg-purple-50/30 dark:bg-purple-900/10">
                  {/* Preset Selection */}
                  <div className="space-y-2">
                    <label className="block text-sm font-semibold text-gray-900 dark:text-white">Provider Preset</label>
                    <div className="relative">
                      <button
                        type="button"
                        onClick={() => setShowPresetDropdown(!showPresetDropdown)}
                        className="w-full p-3 text-left border-2 border-purple-200 dark:border-purple-600 rounded-xl bg-white dark:bg-gray-800 flex items-center justify-between hover:border-purple-400 transition-colors"
                      >
                        <span className="text-gray-900 dark:text-white">{CUSTOM_PROVIDER_PRESETS[customPreset].name}</span>
                        <ChevronDown className={`h-4 w-4 text-gray-500 transition-transform ${showPresetDropdown ? 'rotate-180' : ''}`} />
                      </button>
                      {showPresetDropdown && (
                        <div className="absolute z-20 w-full mt-1 bg-white dark:bg-gray-800 border-2 border-purple-200 dark:border-purple-600 rounded-xl shadow-lg overflow-hidden">
                          {Object.entries(CUSTOM_PROVIDER_PRESETS).map(([key, preset]) => (
                            <button
                              key={key}
                              type="button"
                              onClick={() => {
                                setCustomPreset(key as CustomPresetKey)
                                if (preset.baseUrl) {
                                  setTempCustomBaseUrl(preset.baseUrl)
                                } else {
                                  setTempCustomBaseUrl('')
                                }
                                setShowPresetDropdown(false)
                              }}
                              className={`w-full p-3 text-left hover:bg-purple-50 dark:hover:bg-purple-900/30 transition-colors ${
                                customPreset === key ? 'bg-purple-100 dark:bg-purple-900/50' : ''
                              }`}
                            >
                              <div className="font-medium text-sm text-gray-900 dark:text-white">{preset.name}</div>
                              {preset.baseUrl && (
                                <div className="text-xs text-gray-500 dark:text-gray-400 truncate">{preset.baseUrl}</div>
                              )}
                            </button>
                          ))}
                        </div>
                      )}
                    </div>
                  </div>

                  {/* Base URL */}
                  <div className="space-y-2">
                    <label className="block text-sm font-semibold text-gray-900 dark:text-white">Base URL</label>
                    {customPreset === 'manual' ? (
                      <input
                        type="url"
                        placeholder="https://api.example.com/v1"
                        value={tempCustomBaseUrl}
                        onChange={(e) => setTempCustomBaseUrl(e.target.value)}
                        className="w-full px-4 py-3 border-2 border-purple-200 dark:border-purple-600 rounded-xl text-gray-900 dark:text-white bg-white dark:bg-gray-800 focus:border-purple-500 focus:ring-2 focus:ring-purple-500/20 transition-all"
                      />
                    ) : (
                      <input
                        type="url"
                        value={CUSTOM_PROVIDER_PRESETS[customPreset].baseUrl}
                        readOnly
                        className="w-full px-4 py-3 border-2 border-purple-200 dark:border-purple-600 rounded-xl text-gray-900 dark:text-white bg-gray-100 dark:bg-gray-700 cursor-not-allowed"
                      />
                    )}
                    <p className="text-xs text-gray-500 dark:text-gray-400">
                      {customPreset === 'manual'
                        ? 'Enter any OpenAI-compatible endpoint URL'
                        : `Pre-configured for ${CUSTOM_PROVIDER_PRESETS[customPreset].name}`
                      }
                    </p>
                  </div>

                  {/* Help Link */}
                  {customPreset !== 'manual' && (
                    <a
                      href={CUSTOM_PROVIDER_PRESETS[customPreset].helpUrl}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="inline-flex items-center text-sm text-purple-600 dark:text-purple-400 hover:underline"
                    >
                      Get {CUSTOM_PROVIDER_PRESETS[customPreset].name} API Key
                      <ExternalLink className="h-3 w-3 ml-1" />
                    </a>
                  )}
                </div>
              )}

              {/* Default Model - only show for non-custom providers */}
              {tempProvider !== 'custom' && (
                <div>
                  <label className="block text-sm font-semibold mb-3 text-gray-900 dark:text-white">Default Model</label>
                  <input
                    type="text"
                    value={tempModel}
                    onChange={(e) => setTempModel(e.target.value)}
                    className="w-full px-4 py-3 border-2 border-blue-200 dark:border-blue-700 rounded-xl text-gray-900 dark:text-white bg-blue-50 dark:bg-blue-900/20 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all"
                  />
                  <p className="text-xs text-gray-500 dark:text-gray-400 mt-2">
                    Default: {getDefaultModel(tempProvider)}
                  </p>
                </div>
              )}
              
              <div>
                <label className="block text-sm font-semibold mb-3 text-gray-900 dark:text-white">API Key</label>
                <div className="space-y-3">
                  <div className="relative">
                    <input 
                      type={showApiKey ? 'text' : 'password'} 
                      value={tempApiKey} 
                      onChange={(e) => {
                        setTempApiKey(e.target.value)
                        setTestResult(null)
                        setTestMessage('')
                      }}
                      placeholder="Enter your API key"
                      className="w-full px-4 py-3 pr-12 border-2 border-blue-200 dark:border-blue-700 rounded-xl text-gray-900 dark:text-white bg-blue-50 dark:bg-blue-900/20 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all placeholder-gray-400"
                    />
                    <button
                      type="button"
                      onClick={() => setShowApiKey(!showApiKey)}
                      className="absolute right-3 top-1/2 transform -translate-y-1/2 text-gray-500 dark:text-gray-400 hover:text-gray-700 dark:hover:text-gray-200 transition-colors"
                    >
                      {showApiKey ? (
                        <EyeOff className="h-5 w-5" />
                      ) : (
                        <Eye className="h-5 w-5" />
                      )}
                    </button>
                  </div>
                  
                  {/* Load API Keys Button */}
                  {isAuthenticated && (
                    <LoadApiKeysButton 
                      onSuccess={handleLoadApiKeysSuccess}
                      className="w-full"
                    />
                  )}
                  
                  <Button 
                    onClick={testApiKey}
                    disabled={!tempApiKey.trim() || isTestingApi}
                    variant="outline"
                    className="w-full border-2 border-orange-200 dark:border-orange-700 text-orange-700 dark:text-orange-300 hover:bg-orange-50 dark:hover:bg-orange-900/20"
                  >
                    {isTestingApi ? (
                      <>
                        <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                        Testing API Key...
                      </>
                    ) : (
                      <>
                        <Key className="h-4 w-4 mr-2" />
                        Test API Key
                      </>
                    )}
                  </Button>

                  {testResult && (
                    <div className={`p-4 rounded-xl border-2 ${
                      testResult === 'success' 
                        ? 'bg-green-50 dark:bg-green-900/20 border-green-200 dark:border-green-700' 
                        : 'bg-red-50 dark:bg-red-900/20 border-red-200 dark:border-red-700'
                    }`}>
                      <div className="flex items-start space-x-3">
                        {testResult === 'success' ? (
                          <CheckCircle className="h-5 w-5 text-green-600 dark:text-green-400 flex-shrink-0 mt-0.5" />
                        ) : (
                          <XCircle className="h-5 w-5 text-red-600 dark:text-red-400 flex-shrink-0 mt-0.5" />
                        )}
                        <div>
                          <p className={`text-sm font-semibold ${
                            testResult === 'success' 
                              ? 'text-green-800 dark:text-green-300' 
                              : 'text-red-800 dark:text-red-300'
                          }`}>
                            {testResult === 'success' ? 'Success!' : 'Test Failed'}
                          </p>
                          <p className={`text-xs ${
                            testResult === 'success' 
                              ? 'text-green-700 dark:text-green-400' 
                              : 'text-red-700 dark:text-red-400'
                          }`}>
                            {testMessage}
                          </p>
                        </div>
                      </div>
                    </div>
                  )}
                </div>
              </div>
              
              <div className="bg-blue-50 dark:bg-blue-900/20 border-2 border-blue-200 dark:border-blue-700 rounded-xl p-4">
                <div className="flex items-start space-x-3">
                  <div className="w-6 h-6 bg-blue-500 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                    <span className="text-white text-sm font-bold">🔒</span>
                  </div>
                  <div>
                    <h4 className="font-semibold text-blue-800 dark:text-blue-300 mb-2">Security & Privacy</h4>
                    <p className="text-blue-700 dark:text-blue-400 text-sm leading-relaxed">
                      Your API key and model settings are stored securely in your browser. 
                      {isAuthenticated 
                        ? "When signed in, your API keys are encrypted and saved to your account for access across devices."
                        : "Sign in to save your API keys to your account for secure access from any device."}
                    </p>
                  </div>
                </div>
              </div>
              
              <div className="flex space-x-3 pt-4">
                <Button
                  onClick={handleSaveApiKey}
                  disabled={!tempApiKey.trim() || isSaving || saveSuccess}
                  className={`flex-1 shadow-lg hover:shadow-xl transform hover:scale-105 transition-all font-semibold py-3 ${
                    saveSuccess
                      ? 'bg-green-600 hover:bg-green-600 text-white'
                      : 'bg-gradient-to-r from-blue-600 via-purple-600 to-indigo-600 hover:from-blue-700 hover:via-purple-700 hover:to-indigo-700 text-white'
                  }`}
                >
                  {isSaving ? (
                    <>
                      <Loader2 className="h-5 w-5 mr-2 animate-spin" />
                      Saving...
                    </>
                  ) : saveSuccess ? (
                    <>
                      <CheckCircle className="h-5 w-5 mr-2" />
                      Configuration Saved!
                    </>
                  ) : (
                    <>
                      <CheckCircle className="h-5 w-5 mr-2" />
                      Save Configuration
                    </>
                  )}
                </Button>
                <Button
                  variant="outline"
                  onClick={() => {
                    setTempApiKey(apiKey)
                    setTempProvider(provider)
                    setTempModel(model)
                    setTestResult(null)
                    setTestMessage('')
                    setIsSaving(false)
                    setSaveSuccess(false)
                    setShowApiKeyModal(false)
                  }}
                  className="flex-1 border-2 border-gray-300 dark:border-gray-600 hover:bg-gray-50 dark:hover:bg-gray-800 font-semibold py-3"
                >
                  Cancel
                </Button>
              </div>
            </div>
          </div>
        </div>
      )}

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

      {/* Dashboard Modal */}
      {showDashboard && (
        <div className="fixed inset-0 bg-black/80 backdrop-blur-md flex items-center justify-center z-50 p-4">
          <div className="bg-white dark:bg-gray-900 rounded-2xl max-w-6xl w-full max-h-[90vh] overflow-y-auto border-2 border-blue-200 dark:border-blue-800 shadow-2xl">
            <div className="sticky top-0 bg-white dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 px-6 py-4">
              <div className="flex items-center justify-between">
                <h2 className="text-2xl font-bold gradient-text">User Dashboard</h2>
                <Button 
                  variant="outline" 
                  onClick={() => setShowDashboard(false)}
                  className="border-2 border-gray-300 dark:border-gray-600 hover:bg-gray-50 dark:hover:bg-gray-800"
                >
                  Close
                </Button>
              </div>
            </div>
            <UserDashboard onNavigate={(path) => {
              setShowDashboard(false)
              // Navigate to the path - in a real app, you'd use router.push(path)
              window.location.href = path
            }} />
          </div>
        </div>
      )}
    </div>
  )
}