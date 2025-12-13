import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { createClient } from '@/utils/supabase/client'
import { useAuthStore } from './auth-store'
import modelSettings from '../../public/examples/model_settings.json'
import { ReasoningEffort, getDefaultReasoningEffort } from '../config/model-presets'

export type Provider = 'openrouter' | 'anthropic' | 'openai' | 'custom'

interface ApiKeyState {
  apiKeys: {
    openrouter: string
    anthropic: string
    openai: string
    custom: string
  }
  provider: Provider
  model: string
  reasoningEffort: ReasoningEffort | null
  customBaseUrl: string  // Custom provider base URL
  isLoading: boolean
  error: string | null

  // Core API key actions
  setApiKey: (key: string, provider?: Provider) => Promise<void>
  loadApiKeys: () => Promise<void>
  clearApiKey: (provider?: Provider) => Promise<void>

  // UI state actions
  setProvider: (provider: Provider) => void
  setModel: (model: string) => void
  setReasoningEffort: (effort: ReasoningEffort | null) => void
  setCustomBaseUrl: (url: string) => void
  clearError: () => void

  // Helpers
  getApiKey: (provider?: Provider) => string
  apiKey: string
}

// Simple encryption/decryption for demo
const encryptKey = (key: string): string => btoa(key)
const decryptKey = (encryptedKey: string): string => {
  try {
    return atob(encryptedKey)
  } catch {
    return ''
  }
}

export const useApiKeyStore = create<ApiKeyState>()(
  persist(
    (set, get) => ({
      apiKeys: {
        openrouter: '',
        anthropic: '',
        openai: '',
        custom: ''
      },
      provider: 'openrouter',
      model: modelSettings.providers.openrouter.default_model,
      reasoningEffort: getDefaultReasoningEffort('openrouter', modelSettings.providers.openrouter.default_model),
      customBaseUrl: '',
      isLoading: false,
      error: null,

      clearError: () => set({ error: null }),

      setProvider: (provider: Provider) => {
        set({ provider })
      },

      setModel: (model: string) => {
        const provider = get().provider
        set({
          model,
          reasoningEffort: getDefaultReasoningEffort(provider, model)
        })
      },

      setReasoningEffort: (effort: ReasoningEffort | null) => {
        set({ reasoningEffort: effort })
      },

      setCustomBaseUrl: (url: string) => {
        set({ customBaseUrl: url })
      },
      
      setApiKey: async (key: string, provider?: Provider) => {
        const targetProvider = provider || get().provider
        const authStore = useAuthStore.getState()
        
        // Update local state immediately
        set((state) => ({
          apiKeys: {
            ...state.apiKeys,
            [targetProvider]: key
          }
        }))
        
        // Save to Supabase if authenticated
        if (authStore.user) {
          set({ isLoading: true, error: null })
          
          try {
            const supabase = createClient()
            const encryptedKey = encryptKey(key)
            
            const { error } = await supabase
              .from('user_api_keys')
              .upsert({
                user_id: authStore.user.id,
                provider: targetProvider,
                encrypted_key: encryptedKey
              }, {
                onConflict: 'user_id,provider'
              })
            
            if (error) throw error
            console.log('API key saved to Supabase')
          } catch (error) {
            console.error('Error saving API key:', error)
            set({ error: 'Failed to save API key to account' })
          } finally {
            set({ isLoading: false })
          }
        }
      },
      
      loadApiKeys: async () => {
        const authStore = useAuthStore.getState()
        
        if (!authStore.user) {
          console.log('No user authenticated')
          return
        }
        
        set({ isLoading: true, error: null })
        
        try {
          console.log('Loading API keys for user:', authStore.user.id)
          const supabase = createClient()
          
          const { data, error } = await supabase
            .from('user_api_keys')
            .select('provider, encrypted_key')
            .eq('user_id', authStore.user.id)
          
          if (error) {
            console.error('Supabase error:', error)
            throw error
          }
          
          console.log('Loaded', data?.length || 0, 'API keys from Supabase')
          
          // Update local state with loaded keys
          const loadedKeys = {
            openrouter: '',
            anthropic: '',
            openai: '',
            custom: ''
          }
          
          data?.forEach((keyData) => {
            const provider = keyData.provider as Provider
            if (provider in loadedKeys) {
              loadedKeys[provider] = decryptKey(keyData.encrypted_key)
            }
          })
          
          set({ apiKeys: loadedKeys })
          
          const hasKeys = Object.values(loadedKeys).some(k => k !== '')
          if (!hasKeys) {
            set({ error: 'No API keys found in your account' })
          }
          
        } catch (error) {
          const message = error instanceof Error ? error.message : 'Failed to load API keys'
          set({ error: message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      clearApiKey: async (provider?: Provider) => {
        const targetProvider = provider || get().provider
        const authStore = useAuthStore.getState()
        
        // Clear local state
        set((state) => ({
          apiKeys: {
            ...state.apiKeys,
            [targetProvider]: ''
          }
        }))
        
        // Delete from Supabase if authenticated
        if (authStore.user) {
          try {
            const supabase = createClient()
            const { error } = await supabase
              .from('user_api_keys')
              .delete()
              .eq('user_id', authStore.user.id)
              .eq('provider', targetProvider)
            
            if (error) throw error
          } catch (error) {
            console.error('Error deleting API key:', error)
          }
        }
      },
      
      getApiKey: (provider?: Provider) => {
        const state = get()
        return state.apiKeys[provider || state.provider]
      },

      // Note: Can't use getter here as it loses access to Zustand's get()
      // Use getApiKey() instead or access apiKeys[provider] directly
      apiKey: ''  // Placeholder, use getApiKey() for actual value
    }),
    {
      name: 'cassia-api-key-storage',
      partialize: (state) => ({
        provider: state.provider,
        model: state.model,
        reasoningEffort: state.reasoningEffort,
        customBaseUrl: state.customBaseUrl,
        apiKeys: state.apiKeys
      })
    }
  )
)