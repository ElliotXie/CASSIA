import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import modelSettings from '../../public/examples/model_settings.json'

type Provider = 'openrouter' | 'anthropic' | 'openai'

interface ApiKeyState {
  apiKeys: {
    openrouter: string
    anthropic: string
    openai: string
  }
  provider: Provider
  model: string
  setApiKey: (key: string, provider?: Provider) => void
  setProvider: (provider: Provider) => void
  setModel: (model: string) => void
  clearApiKey: (provider?: Provider) => void
  getApiKey: (provider?: Provider) => string
  // Compatibility property for backward compatibility
  apiKey?: string
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
      setApiKey: (key, provider) => 
        set((state) => ({
          apiKeys: {
            ...state.apiKeys,
            [provider || state.provider]: key
          }
        })),
      setProvider: (provider) => set({ provider }),
      setModel: (model) => set({ model }),
      clearApiKey: (provider) => 
        set((state) => ({
          apiKeys: provider 
            ? { ...state.apiKeys, [provider]: '' }
            : { openrouter: '', anthropic: '', openai: '' }
        })),
      getApiKey: (provider) => {
        const state = get()
        return state.apiKeys[provider || state.provider]
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
        apiKeys: state.apiKeys,
        provider: state.provider,
        model: state.model
      }) as Partial<ApiKeyState>,
      version: 1,
      skipHydration: false,
      migrate: (persistedState: any, version: number) => {
        // Handle migration from any version to version 1
        if (version === 0 || version === 2) {
          console.log(`Migrating API key store from version ${version} to 1`)
          const oldState = persistedState as any
          return {
            apiKeys: {
              openrouter: oldState.apiKeys?.openrouter || '',
              anthropic: oldState.apiKeys?.anthropic || '',
              openai: oldState.apiKeys?.openai || ''
            },
            provider: oldState.provider || 'openrouter',
            model: oldState.model || modelSettings.providers.openrouter.default_model
          }
        }
        return persistedState
      },
      onRehydrateStorage: () => (state) => {
        if (state) {
          console.log('API key storage hydrated successfully')
        }
      }
    }
  )
)