import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import modelSettings from '../../public/examples/model_settings.json'
import { createClient } from '../../utils/supabase/client'
import { Provider } from '../supabase/types'
import { useAuthStore } from './auth-store'

interface ApiKeyState {
  apiKeys: {
    openrouter: string
    anthropic: string
    openai: string
  }
  provider: Provider
  model: string
  isLoading: boolean
  error: string | null
  
  // Core API key actions
  setApiKey: (key: string, provider?: Provider) => Promise<void>
  setProvider: (provider: Provider) => void
  setModel: (model: string) => void
  clearApiKey: (provider?: Provider) => Promise<void>
  getApiKey: (provider?: Provider) => string
  
  // Supabase integration
  loadApiKeys: () => Promise<void>
  syncApiKeys: () => Promise<void>
  
  // Utility functions
  clearError: () => void
  getCurrentUser: () => Promise<any>
  
  // Compatibility properties
  apiKey?: string
}

// Simple encryption/decryption functions (for demo purposes)
// In production, use proper encryption libraries
const encryptKey = (key: string): string => {
  return btoa(key) // Base64 encoding
}

const decryptKey = (encryptedKey: string): string => {
  try {
    return atob(encryptedKey) // Base64 decoding
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
        openai: ''
      },
      provider: 'openrouter',
      model: modelSettings.providers.openrouter.default_model,
      isLoading: false,
      error: null,
      
      clearError: () => set({ error: null }),
      
      getCurrentUser: async () => {
        // Get user from auth store for consistency
        const authStore = useAuthStore.getState()
        return authStore.user
      },
      
      setApiKey: async (key: string, provider?: Provider) => {
        const targetProvider = provider || get().provider
        const supabase = createClient()
        
        // Update local state immediately for better UX
        set((state) => ({
          apiKeys: {
            ...state.apiKeys,
            [targetProvider]: key
          }
        }))
        
        // Get current user directly from Supabase
        try {
          const user = await get().getCurrentUser()
          
          if (user) {
            set({ isLoading: true, error: null })
            
            try {
              const encryptedKey = encryptKey(key)
              
              // Upsert API key in Supabase
              const { error } = await supabase
                .from('user_api_keys')
                .upsert({
                  user_id: user.id,
                  provider: targetProvider,
                  encrypted_key: encryptedKey
                }, {
                  onConflict: 'user_id,provider'
                })
              
              if (error) {
                console.error('Supabase error saving API key:', error)
                throw error
              }
              
              console.log(`✅ API key for ${targetProvider} saved securely to Supabase`)
            } catch (error) {
              set({ error: (error as Error).message })
              console.error('Error saving API key to Supabase:', error)
            } finally {
              set({ isLoading: false })
            }
          } else {
            console.log('User not authenticated, API key saved locally only')
          }
        } catch (error) {
          console.error('Error getting current user:', error)
          set({ error: 'Failed to authenticate user' })
        }
      },
      
      setProvider: (provider: Provider) => set({ provider }),
      setModel: (model: string) => set({ model }),
      
      clearApiKey: async (provider?: Provider) => {
        const targetProvider = provider || get().provider
        const supabase = createClient()
        
        // Update local state
        if (targetProvider) {
          set((state) => ({
            apiKeys: { ...state.apiKeys, [targetProvider]: '' }
          }))
        } else {
          set({
            apiKeys: { openrouter: '', anthropic: '', openai: '' }
          })
        }
        
        // Get current user directly from Supabase
        try {
          const user = await get().getCurrentUser()
          
          if (user) {
            set({ isLoading: true, error: null })
            
            try {
              if (targetProvider) {
                // Clear specific provider
                const { error } = await supabase
                  .from('user_api_keys')
                  .delete()
                  .eq('user_id', user.id)
                  .eq('provider', targetProvider)
                
                if (error) throw error
              } else {
                // Clear all providers
                const { error } = await supabase
                  .from('user_api_keys')
                  .delete()
                  .eq('user_id', user.id)
                
                if (error) throw error
              }
              
              console.log(`API key(s) cleared for ${targetProvider || 'all providers'}`)
            } catch (error) {
              set({ error: (error as Error).message })
              console.error('Error clearing API key:', error)
            } finally {
              set({ isLoading: false })
            }
          }
        } catch (error) {
          console.error('Error getting current user for clearing keys:', error)
        }
      },
      
      getApiKey: (provider?: Provider) => {
        const state = get()
        return state.apiKeys[provider || state.provider]
      },
      
      loadApiKeys: async () => {
        const supabase = createClient()
        
        try {
          const user = await get().getCurrentUser()
          
          if (!user) {
            console.log('No authenticated user, skipping API key loading')
            return
          }
          
          console.log('Loading API keys for user:', user.id)
          set({ isLoading: true, error: null })
          
          try {
            console.log('Querying user_api_keys table...')
            const startTime = Date.now()
            
            // Add timeout to the query
            const queryPromise = supabase
              .from('user_api_keys')
              .select('*')
              .eq('user_id', user.id)
            
            const timeoutPromise = new Promise((_, reject) => 
              setTimeout(() => reject(new Error('Query timeout after 10 seconds')), 10000)
            )
            
            const { data, error } = await Promise.race([
              queryPromise,
              timeoutPromise
            ]) as any
            
            const queryTime = Date.now() - startTime
            console.log(`Query completed in ${queryTime}ms`)
            
            if (error) {
              console.error('Supabase query error:', error)
              throw error
            }
            
            console.log('Query result:', data?.length || 0, 'keys found')
            
            // Decrypt and populate API keys
            const decryptedKeys = {
              openrouter: '',
              anthropic: '',
              openai: ''
            }
            
            data?.forEach((keyData) => {
              console.log('Processing key for provider:', keyData.provider)
              const provider = keyData.provider as Provider
              if (provider in decryptedKeys) {
                try {
                  decryptedKeys[provider] = decryptKey(keyData.encrypted_key)
                  console.log(`Decrypted key for ${provider}`)
                } catch (decryptError) {
                  console.error(`Failed to decrypt key for ${provider}:`, decryptError)
                }
              }
            })
            
            set({ apiKeys: decryptedKeys })
            const loadedProviders = Object.keys(decryptedKeys).filter(k => decryptedKeys[k as keyof typeof decryptedKeys])
            console.log('✅ API keys loaded from Supabase:', loadedProviders.length > 0 ? loadedProviders : 'No keys found')
            
            if (loadedProviders.length === 0) {
              console.log('ℹ️ No API keys found in your account. Please save an API key first.')
            }
          } catch (error) {
            const errorMessage = error instanceof Error ? error.message : 'Failed to load API keys'
            set({ error: errorMessage })
            console.error('Error loading API keys:', error)
            throw error // Re-throw to let the button component handle it
          } finally {
            set({ isLoading: false })
          }
        } catch (error) {
          console.error('Error getting current user for loading keys:', error)
          throw new Error('Not authenticated. Please sign in first.')
        }
      },
      
      syncApiKeys: async () => {
        console.log('Syncing API keys...')
        try {
          const user = await get().getCurrentUser()
          
          if (!user) {
            // If not authenticated, clear API keys
            set({
              apiKeys: { openrouter: '', anthropic: '', openai: '' }
            })
            console.log('User not authenticated, clearing API keys')
            return
          }
          
          console.log('User authenticated, loading API keys from Supabase')
          // Load API keys from Supabase
          await get().loadApiKeys()
        } catch (error) {
          console.error('Error syncing API keys:', error)
        }
      },
      
      // Compatibility getter for backward compatibility
      get apiKey() {
        const state = get()
        return state.apiKeys[state.provider]
      }
    }),
    {
      name: 'cassia-api-key-storage',
      storage: {
        getItem: (name) => {
          try {
            const item = typeof window !== 'undefined' ? localStorage.getItem(name) : null
            return item ? JSON.parse(item) : null
          } catch (error) {
            console.error('Error reading from localStorage:', error)
            return null
          }
        },
        setItem: (name, value) => {
          try {
            if (typeof window !== 'undefined') {
              localStorage.setItem(name, JSON.stringify(value))
            }
          } catch (error) {
            console.error('Error writing to localStorage:', error)
          }
        },
        removeItem: (name) => {
          try {
            if (typeof window !== 'undefined') {
              localStorage.removeItem(name)
            }
          } catch (error) {
            console.error('Error removing from localStorage:', error)
          }
        }
      },
      partialize: (state) => ({
        // Only persist provider and model settings
        // API keys will be loaded from Supabase
        provider: state.provider,
        model: state.model,
        // Keep API keys for offline fallback
        apiKeys: state.apiKeys
      }) as Partial<ApiKeyState>,
      version: 2, // Increment version for migration
      skipHydration: false,
      migrate: (persistedState: any, version: number) => {
        // Migration from version 1 to 2
        if (version === 1) {
          console.log('Migrating API key store from version 1 to 2')
          // Convert old format to new format
          const oldState = persistedState as any
          return {
            apiKeys: {
              openrouter: oldState.apiKeys?.openrouter || '',
              anthropic: oldState.apiKeys?.anthropic || '',
              openai: oldState.apiKeys?.openai || ''
            },
            provider: oldState.provider || 'openrouter',
            model: oldState.model || modelSettings.providers.openrouter.default_model,
            isLoading: false,
            error: null
          }
        }
        return persistedState
      },
      onRehydrateStorage: () => (state) => {
        if (state) {
          console.log('API key storage hydrated successfully')
          // Sync API keys will be triggered by AuthProvider when authenticated
        }
      }
    }
  )
)