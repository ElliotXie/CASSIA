'use client'

import Link from 'next/link'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { ArrowRight, Zap, Settings, Download, HelpCircle, Key, CheckCircle, XCircle, Loader2, BookOpen, FileText, Cpu, Github, FileX, Eye, EyeOff } from 'lucide-react'
import { useState } from 'react'
import { useApiKeyStore } from '@/lib/stores/api-key-store-simple'
import { useAuthStore } from '@/lib/stores/auth-store'
import { ContactDialog } from '@/components/ContactDialog'
import { AuthButton } from '@/components/auth/AuthButton'
import { UserDashboard } from '@/components/dashboard/UserDashboard'
import { LoadApiKeysButton } from '@/components/LoadApiKeysButton'
import modelSettings from '../public/examples/model_settings.json'

export default function HomePage() {
  const [showApiKeyModal, setShowApiKeyModal] = useState(false)
  const [showContactModal, setShowContactModal] = useState(false)
  const [showDashboard, setShowDashboard] = useState(false)
  const { apiKeys, provider, model, setApiKey, setProvider, setModel } = useApiKeyStore()
  const { isAuthenticated } = useAuthStore()
  const [tempApiKey, setTempApiKey] = useState(apiKeys[provider])
  const [tempProvider, setTempProvider] = useState(provider)
  const [tempModel, setTempModel] = useState(model)
  const [isTestingApi, setIsTestingApi] = useState(false)
  const [testResult, setTestResult] = useState<'success' | 'error' | null>(null)
  const [testMessage, setTestMessage] = useState('')
  const [showApiKey, setShowApiKey] = useState(false)
  
  // Get current API key
  const apiKey = apiKeys[provider]

  // Default models for each provider from model_settings.json
  const getDefaultModel = (providerName: string) => {
    const providerData = modelSettings.providers[providerName as keyof typeof modelSettings.providers];
    return providerData?.default_model || modelSettings.providers.openrouter.default_model;
  }

  // Update model when provider changes
  const handleProviderChange = (newProvider: string) => {
    setTempProvider(newProvider)
    setTempModel(getDefaultModel(newProvider))
    setTestResult(null)
    setTestMessage('')
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
    await setApiKey(tempApiKey, tempProvider)
    setProvider(tempProvider)
    setModel(tempModel)
    setShowApiKeyModal(false)
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
                <p className="text-sm text-gray-600 dark:text-gray-300">Cell Annotation Intelligence Assistant</p>
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
                  setTestResult(null)
                  setTestMessage('')
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
          
          <h2 className="text-5xl lg:text-6xl font-bold mb-6">
            <span className="gradient-text">Automated Cell Type</span>
            <br />
            <span className="text-gray-900 dark:text-white">Annotation Platform</span>
          </h2>
          
          <p className="text-xl text-gray-600 dark:text-gray-300 max-w-4xl mx-auto mb-10 leading-relaxed">
            Transform your single-cell RNA-seq analysis with AI-powered precision. Upload marker data, 
            configure your analysis, and receive comprehensive cell type annotations with detailed reports. 
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

        {/* Prominent API Key Configuration Section */}
        <div className="glass rounded-2xl p-8 mb-16 border border-yellow-400/30 bg-yellow-50/20">
          <div className="flex items-start justify-between">
            <div className="flex items-start space-x-4">
              <div className="w-16 h-16 bg-gradient-to-br from-yellow-500 to-orange-600 rounded-xl flex items-center justify-center shadow-lg animate-glow">
                <Key className="h-8 w-8 text-white" />
              </div>
              <div className="flex-1">
                <h3 className="text-2xl font-bold text-gray-900 dark:text-white mb-2 flex items-center">
                  🔑 API Configuration Required
                  {!apiKey && <span className="ml-3 px-3 py-1 bg-red-500 text-white text-sm rounded-full animate-pulse">Required</span>}
                  {apiKey && <span className="ml-3 px-3 py-1 bg-green-500 text-white text-sm rounded-full">✓ Configured</span>}
                </h3>
                <p className="text-gray-600 dark:text-gray-300 mb-4 text-lg">
                  {!apiKey 
                    ? "Set up your AI provider API key to start analyzing your data. This is required for all analysis features."
                    : `API key configured for ${provider}. You're ready to run analysis!`
                  }
                </p>
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4 text-sm">
                  <div className="flex items-center space-x-2">
                    <span className="w-2 h-2 bg-green-500 rounded-full"></span>
                    <span className="text-gray-700 dark:text-gray-300">OpenRouter (Recommended)</span>
                  </div>
                  <div className="flex items-center space-x-2">
                    <span className="w-2 h-2 bg-blue-500 rounded-full"></span>
                    <span className="text-gray-700 dark:text-gray-300">OpenAI GPT Models</span>
                  </div>
                  <div className="flex items-center space-x-2">
                    <span className="w-2 h-2 bg-purple-500 rounded-full"></span>
                    <span className="text-gray-700 dark:text-gray-300">Anthropic Claude</span>
                  </div>
                </div>
              </div>
            </div>
            <div className="flex flex-col space-y-3 ml-6">
              <Button 
                onClick={() => {
                  setTempApiKey(apiKey)
                  setTempProvider(provider)
                  setTempModel(model || getDefaultModel(provider))
                  setTestResult(null)
                  setTestMessage('')
                  setShowApiKeyModal(true)
                }}
                size="lg"
                className={`px-8 py-4 text-lg font-semibold shadow-lg hover:shadow-xl transform hover:scale-105 transition-all duration-300 btn-modern ${
                  !apiKey 
                    ? 'bg-gradient-to-r from-yellow-500 to-orange-600 hover:from-yellow-600 hover:to-orange-700 text-white animate-pulse'
                    : 'bg-gradient-to-r from-green-600 to-emerald-600 hover:from-green-700 hover:to-emerald-700 text-white'
                }`}
              >
                <Key className="h-5 w-5 mr-2" />
                {!apiKey ? 'Configure API Key' : 'Update API Key'}
              </Button>
              {apiKey && (
                <div className="text-center">
                  <span className="text-sm text-gray-600 dark:text-gray-400">Provider: {provider}</span>
                </div>
              )}
            </div>
            
            {/* Load API Keys from Supabase */}
            {isAuthenticated && (
              <div className="text-center">
                <LoadApiKeysButton />
              </div>
            )}
          </div>
          
          {!apiKey && (
            <div className="mt-6 p-4 glass border border-orange-400/30 bg-orange-50/20 rounded-xl">
              <div className="flex items-start space-x-3">
                <div className="w-6 h-6 bg-orange-500 rounded-full flex items-center justify-center flex-shrink-0 mt-0.5">
                  <span className="text-white text-sm font-bold">!</span>
                </div>
                <div>
                  <h4 className="font-semibold text-orange-800 dark:text-orange-400 mb-2">Why do I need an API key?</h4>
                  <p className="text-orange-700 dark:text-orange-300 text-sm leading-relaxed">
                    CASSIA uses advanced AI language models to analyze your cell marker data and provide intelligent annotations. 
                    You need an API key from OpenRouter, OpenAI, or Anthropic to access these models. Your key is stored securely 
                    in your browser and never shared with our servers.
                  </p>
                </div>
              </div>
            </div>
          )}
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
                <Link href="/pipeline" className="flex flex-col items-center justify-center">
                  <span>Run CASSIA Pipeline</span>
                  <ArrowRight className="h-4 w-4 mt-1" />
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
                <Link href="/batch" className="flex flex-col items-center justify-center">
                  <span>Run CASSIA Batch</span>
                  <ArrowRight className="h-4 w-4 mt-1" />
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
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/agents/symphony-compare">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Symphony Compare</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Multi-model consensus</div>
                    </div>
                  </Link>
                </Button>
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
                <Button variant="outline" size="sm" asChild className="justify-center h-auto py-3 glass border-white/30 hover:bg-white/20 btn-modern">
                  <Link href="/subclustering">
                    <div className="text-center">
                      <div className="font-medium text-gray-900 dark:text-white">Subclustering</div>
                      <div className="text-xs text-gray-600 dark:text-gray-400">Subtype analysis</div>
                    </div>
                  </Link>
                </Button>
              </div>
              <Button variant="outline" asChild className="w-full glass border-white/30 hover:bg-white/20 btn-modern">
                <Link href="/agents" className="flex flex-col items-center justify-center">
                  <span>View All Agents</span>
                  <ArrowRight className="h-4 w-4 mt-1" />
                </Link>
              </Button>
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
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://cassia-documentation-en-new.vercel.app/" target="_blank" rel="noopener noreferrer">
                  <span>📚 Complete R Documentation/Vignette</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA/blob/main/CASSIA_example/CASSIA_python_tutorial.ipynb" target="_blank" rel="noopener noreferrer">
                  <span>📝 Example Python workflow/Vignette</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://sc-llm-benchmark.com/methods/cassia" target="_blank" rel="noopener noreferrer">
                  <span>🤖 LLMs Annotation Benchmark</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2" target="_blank" rel="noopener noreferrer">
                  <span>📄 BioRxiv Preprint</span>
                  <ArrowRight className="h-4 w-4" />
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
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer">
                  <span>📦 R & Python Package Installation</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://cassia-documentation-en-new.vercel.app/" target="_blank" rel="noopener noreferrer">
                  <span>📚 Complete Documentation</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA/issues" target="_blank" rel="noopener noreferrer">
                  <span>💬 Community Support</span>
                  <ArrowRight className="h-4 w-4" />
                </Link>
              </Button>
              <Button variant="outline" className="w-full justify-between" asChild>
                <Link href="https://github.com/ElliotXie/CASSIA/issues" target="_blank" rel="noopener noreferrer">
                  <span>🐛 Report Issues</span>
                  <ArrowRight className="h-4 w-4" />
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
              Cell Annotation using Single-cell Sequencing Intelligence Assistant
              <br />
              <span className="text-sm">Built with Next.js, powered by Large Language Models</span>
            </p>
            <div className="flex flex-wrap justify-center gap-6 text-sm text-gray-500 dark:text-gray-400">
              <a href="https://cassia-documentation-en-new.vercel.app/" target="_blank" rel="noopener noreferrer" className="hover:text-blue-600 dark:hover:text-blue-400 transition-colors">Documentation</a>
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
              <p className="text-gray-600 dark:text-gray-400 text-lg">Set up your AI provider and default model</p>
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
                </select>
              </div>

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
                  disabled={!tempApiKey.trim()}
                  className="flex-1 bg-gradient-to-r from-blue-600 via-purple-600 to-indigo-600 hover:from-blue-700 hover:via-purple-700 hover:to-indigo-700 text-white shadow-lg hover:shadow-xl transform hover:scale-105 transition-all font-semibold py-3"
                >
                  <CheckCircle className="h-5 w-5 mr-2" />
                  Save Configuration
                </Button>
                <Button 
                  variant="outline" 
                  onClick={() => {
                    setTempApiKey(apiKey)
                    setTempProvider(provider)
                    setTempModel(model)
                    setTestResult(null)
                    setTestMessage('')
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