import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { createClient } from '@/utils/supabase/client'
import { useAuthStore } from './auth-store'
import modelSettings from '../../public/examples/model_settings.json'
import { ReasoningEffort, getDefaultReasoningEffort } from '../config/model-presets'

export type Provider = 'openrouter' | 'anthropic' | 'openai' | 'custom'

// Custom provider preset types
export type CustomPresetKey = 'deepseek' | 'qwen' | 'kimi' | 'siliconflow' | 'minimax' | 'zhipuai' | 'manual'

export interface CustomProviderConfig {
  apiKey: string
  baseUrl: string
}

// Default base URLs for each custom provider preset
const DEFAULT_CUSTOM_PROVIDERS: Record<CustomPresetKey, CustomProviderConfig> = {
  deepseek: { apiKey: '', baseUrl: 'https://api.deepseek.com' },
  qwen: { apiKey: '', baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1' },
  kimi: { apiKey: '', baseUrl: 'https://api.moonshot.cn/v1' },
  siliconflow: { apiKey: '', baseUrl: 'https://api.siliconflow.cn/v1' },
  minimax: { apiKey: '', baseUrl: 'https://api.minimax.io/v1' },
  zhipuai: { apiKey: '', baseUrl: 'https://open.bigmodel.cn/api/paas/v4' },
  manual: { apiKey: '', baseUrl: '' }
}

interface ApiKeyState {
  // Standard provider API keys (openrouter, anthropic, openai)
  apiKeys: {
    openrouter: string
    anthropic: string
    openai: string
  }

  // Per-preset custom provider storage
  customProviders: Record<CustomPresetKey, CustomProviderConfig>
  selectedCustomPreset: CustomPresetKey

  provider: Provider
  model: string
  reasoningEffort: ReasoningEffort | null
  isLoading: boolean
  error: string | null

  // Core API key actions
  setApiKey: (key: string, provider?: Provider) => Promise<void>
  loadApiKeys: () => Promise<void>
  clearApiKey: (provider?: Provider) => Promise<void>

  // Custom provider specific actions
  setCustomProviderKey: (preset: CustomPresetKey, key: string) => Promise<void>
  setCustomProviderBaseUrl: (preset: CustomPresetKey, url: string) => void
  setSelectedCustomPreset: (preset: CustomPresetKey) => void
  clearCustomProviderKey: (preset: CustomPresetKey) => Promise<void>

  // Custom provider getters
  getCustomApiKey: (preset?: CustomPresetKey) => string
  getCustomBaseUrl: (preset?: CustomPresetKey) => string

  // UI state actions
  setProvider: (provider: Provider) => void
  setModel: (model: string) => void
  setReasoningEffort: (effort: ReasoningEffort | null) => void
  setCustomBaseUrl: (url: string) => void  // Legacy - updates selected preset's baseUrl
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
        openai: ''
      },
      customProviders: { ...DEFAULT_CUSTOM_PROVIDERS },
      selectedCustomPreset: 'deepseek',
      provider: 'openrouter',
      model: modelSettings.providers.openrouter.default_model,
      reasoningEffort: getDefaultReasoningEffort('openrouter', modelSettings.providers.openrouter.default_model),
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

      // Legacy method - updates the currently selected custom preset's baseUrl
      setCustomBaseUrl: (url: string) => {
        const preset = get().selectedCustomPreset
        set((state) => ({
          customProviders: {
            ...state.customProviders,
            [preset]: { ...state.customProviders[preset], baseUrl: url }
          }
        }))
      },

      setSelectedCustomPreset: (preset: CustomPresetKey) => {
        set({ selectedCustomPreset: preset })
      },

      setCustomProviderBaseUrl: (preset: CustomPresetKey, url: string) => {
        set((state) => ({
          customProviders: {
            ...state.customProviders,
            [preset]: { ...state.customProviders[preset], baseUrl: url }
          }
        }))
      },

      setCustomProviderKey: async (preset: CustomPresetKey, key: string) => {
        const authStore = useAuthStore.getState()

        // Update local state immediately
        set((state) => ({
          customProviders: {
            ...state.customProviders,
            [preset]: { ...state.customProviders[preset], apiKey: key }
          }
        }))

        // Save to Supabase if authenticated with custom:preset format
        if (authStore.user) {
          set({ isLoading: true, error: null })

          try {
            const supabase = createClient()
            const encryptedKey = encryptKey(key)
            const providerKey = `custom:${preset}`

            const { error } = await supabase
              .from('user_api_keys')
              .upsert({
                user_id: authStore.user.id,
                provider: providerKey,
                encrypted_key: encryptedKey
              }, {
                onConflict: 'user_id,provider'
              })

            if (error) throw error
            console.log(`API key saved to Supabase for ${providerKey}`)
          } catch (error) {
            console.error('Error saving custom provider API key:', error)
            set({ error: 'Failed to save API key to account' })
          } finally {
            set({ isLoading: false })
          }
        }
      },

      clearCustomProviderKey: async (preset: CustomPresetKey) => {
        const authStore = useAuthStore.getState()

        // Clear local state
        set((state) => ({
          customProviders: {
            ...state.customProviders,
            [preset]: { ...state.customProviders[preset], apiKey: '' }
          }
        }))

        // Delete from Supabase if authenticated
        if (authStore.user) {
          try {
            const supabase = createClient()
            const providerKey = `custom:${preset}`

            const { error } = await supabase
              .from('user_api_keys')
              .delete()
              .eq('user_id', authStore.user.id)
              .eq('provider', providerKey)

            if (error) throw error
          } catch (error) {
            console.error('Error deleting custom provider API key:', error)
          }
        }
      },

      getCustomApiKey: (preset?: CustomPresetKey) => {
        const state = get()
        const targetPreset = preset || state.selectedCustomPreset
        return state.customProviders[targetPreset]?.apiKey || ''
      },

      getCustomBaseUrl: (preset?: CustomPresetKey) => {
        const state = get()
        const targetPreset = preset || state.selectedCustomPreset
        return state.customProviders[targetPreset]?.baseUrl || ''
      },

      setApiKey: async (key: string, provider?: Provider) => {
        const targetProvider = provider || get().provider
        const authStore = useAuthStore.getState()

        // For custom provider, delegate to setCustomProviderKey
        if (targetProvider === 'custom') {
          const preset = get().selectedCustomPreset
          await get().setCustomProviderKey(preset, key)
          return
        }

        // Update local state immediately for standard providers
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
            openai: ''
          }

          const loadedCustomProviders = { ...DEFAULT_CUSTOM_PROVIDERS }

          data?.forEach((keyData) => {
            const providerValue = keyData.provider

            // Handle custom:preset format (e.g., custom:deepseek, custom:qwen)
            if (providerValue.startsWith('custom:')) {
              const preset = providerValue.replace('custom:', '') as CustomPresetKey
              if (preset in loadedCustomProviders) {
                loadedCustomProviders[preset] = {
                  ...loadedCustomProviders[preset],
                  apiKey: decryptKey(keyData.encrypted_key)
                }
              }
            }
            // Handle legacy 'custom' provider (migrate to manual)
            else if (providerValue === 'custom') {
              loadedCustomProviders.manual = {
                ...loadedCustomProviders.manual,
                apiKey: decryptKey(keyData.encrypted_key)
              }
            }
            // Handle standard providers
            else if (providerValue in loadedKeys) {
              loadedKeys[providerValue as keyof typeof loadedKeys] = decryptKey(keyData.encrypted_key)
            }
          })

          set({ apiKeys: loadedKeys, customProviders: loadedCustomProviders })

          const hasStandardKeys = Object.values(loadedKeys).some(k => k !== '')
          const hasCustomKeys = Object.values(loadedCustomProviders).some(c => c.apiKey !== '')

          if (!hasStandardKeys && !hasCustomKeys) {
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

        // For custom provider, delegate to clearCustomProviderKey
        if (targetProvider === 'custom') {
          const preset = get().selectedCustomPreset
          await get().clearCustomProviderKey(preset)
          return
        }

        // Clear local state for standard providers
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
        const targetProvider = provider || state.provider

        // For custom provider, return the selected preset's API key
        if (targetProvider === 'custom') {
          return state.customProviders[state.selectedCustomPreset]?.apiKey || ''
        }

        return state.apiKeys[targetProvider as keyof typeof state.apiKeys] || ''
      },

      // Placeholder, use getApiKey() for actual value
      apiKey: ''
    }),
    {
      name: 'cassia-api-key-storage',
      version: 2,
      partialize: (state) => ({
        provider: state.provider,
        model: state.model,
        reasoningEffort: state.reasoningEffort,
        customProviders: state.customProviders,
        selectedCustomPreset: state.selectedCustomPreset,
        apiKeys: state.apiKeys
      }),
      // Migration from old format to new format
      migrate: (persistedState: unknown, version: number) => {
        const state = persistedState as Record<string, unknown>

        if (version < 2) {
          console.log('Migrating API key store from version', version, 'to 2')

          // Migrate legacy custom key to customProviders.manual
          const oldApiKeys = state.apiKeys as Record<string, string> | undefined
          const oldCustomBaseUrl = state.customBaseUrl as string | undefined

          const newCustomProviders = { ...DEFAULT_CUSTOM_PROVIDERS }

          if (oldApiKeys?.custom) {
            newCustomProviders.manual = {
              apiKey: oldApiKeys.custom,
              baseUrl: oldCustomBaseUrl || ''
            }
            console.log('Migrated legacy custom key to manual preset')
          }

          // Build new apiKeys without custom field
          const newApiKeys = {
            openrouter: oldApiKeys?.openrouter || '',
            anthropic: oldApiKeys?.anthropic || '',
            openai: oldApiKeys?.openai || ''
          }

          return {
            ...state,
            apiKeys: newApiKeys,
            customProviders: newCustomProviders,
            selectedCustomPreset: oldApiKeys?.custom ? 'manual' : 'deepseek'
          }
        }

        return state
      }
    }
  )
)
