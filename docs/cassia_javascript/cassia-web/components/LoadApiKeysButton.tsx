'use client'

import { useState } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store-simple'
import { Button } from '@/components/ui/button'
import { Download, Loader2, CheckCircle, AlertCircle } from 'lucide-react'

interface LoadApiKeysButtonProps {
  onSuccess?: () => void
  className?: string
  size?: 'default' | 'sm' | 'lg' | 'icon'
}

export function LoadApiKeysButton({ onSuccess, className, size = 'default' }: LoadApiKeysButtonProps) {
  const { isAuthenticated, user } = useAuthStore()
  const { loadApiKeys, isLoading } = useApiKeyStore()
  const [loadStatus, setLoadStatus] = useState<'idle' | 'loading' | 'success' | 'error'>('idle')
  const [errorMessage, setErrorMessage] = useState<string>('')

  const handleLoadApiKeys = async () => {
    if (!isAuthenticated || !user) {
      setErrorMessage('Please sign in to load API keys')
      setLoadStatus('error')
      setTimeout(() => setLoadStatus('idle'), 3000)
      return
    }

    setLoadStatus('loading')
    setErrorMessage('')

    try {
      await loadApiKeys()
      
      // Check if any API keys were loaded
      const { apiKeys: newApiKeys, customProviders } = useApiKeyStore.getState()
      const hasLoadedStandardKeys = Object.values(newApiKeys).some((key) => key !== '')
      const hasLoadedCustomKeys = Object.values(customProviders).some((providerConfig) => providerConfig.apiKey !== '')
      const hasLoadedKeys = hasLoadedStandardKeys || hasLoadedCustomKeys
      
      if (hasLoadedKeys) {
        setLoadStatus('success')
        onSuccess?.()
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

  if (!isAuthenticated) {
    return null
  }

  return (
    <div className="flex flex-col gap-2">
      <Button
        onClick={handleLoadApiKeys}
        disabled={isLoading || loadStatus === 'loading'}
        variant={loadStatus === 'success' ? 'default' : loadStatus === 'error' ? 'destructive' : 'outline'}
        className={className || "w-full sm:w-auto"}
        size={size}
      >
        {loadStatus === 'loading' ? (
          <>
            <Loader2 className="mr-2 h-4 w-4 animate-spin" />
            Loading API Keys...
          </>
        ) : loadStatus === 'success' ? (
          <>
            <CheckCircle className="mr-2 h-4 w-4" />
            API Keys Loaded
          </>
        ) : loadStatus === 'error' ? (
          <>
            <AlertCircle className="mr-2 h-4 w-4" />
            Failed to Load
          </>
        ) : (
          <>
            <Download className="mr-2 h-4 w-4" />
            Load API Keys from Account
          </>
        )}
      </Button>
      
      {loadStatus === 'error' && errorMessage && (
        <p className="text-sm text-red-600 dark:text-red-400">
          {errorMessage}
        </p>
      )}
      
      {loadStatus === 'success' && (
        <p className="text-sm text-green-600 dark:text-green-400">
          API keys loaded successfully from your account
        </p>
      )}
    </div>
  )
}
