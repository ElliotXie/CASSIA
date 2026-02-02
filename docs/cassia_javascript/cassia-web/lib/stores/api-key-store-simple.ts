import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { useAuthStore } from './auth-store'
import modelSettings from '../../public/examples/model_settings.json'
import { ReasoningEffort, getDefaultReasoningEffort } from '../config/model-presets'

// ─── Public types (consumed by all pages/components) ───────────────────────

export type Provider = 'openrouter' | 'anthropic' | 'openai' | 'custom'

export type CustomPresetKey = 'deepseek' | 'qwen' | 'kimi' | 'siliconflow' | 'minimax' | 'zhipuai' | 'manual'

export interface CustomProviderConfig {
  apiKey: string
  baseUrl: string
}

// ─── Constants ─────────────────────────────────────────────────────────────

const DEFAULT_CUSTOM_PROVIDERS: Record<CustomPresetKey, CustomProviderConfig> = {
  deepseek: { apiKey: '', baseUrl: 'https://api.deepseek.com' },
  qwen: { apiKey: '', baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1' },
  kimi: { apiKey: '', baseUrl: 'https://api.moonshot.cn/v1' },
  siliconflow: { apiKey: '', baseUrl: 'https://api.siliconflow.cn/v1' },
  minimax: { apiKey: '', baseUrl: 'https://api.minimax.io/v1' },
  zhipuai: { apiKey: '', baseUrl: 'https://open.bigmodel.cn/api/paas/v4' },
  manual: { apiKey: '', baseUrl: '' }
}

// ─── Encryption (base64, demo-level) ──────────────────────────────────────

const encryptKey = (key: string): string => btoa(key)
const decryptKey = (enc: string): string => {
  try { return atob(enc) } catch { return '' }
}

// ─── Supabase REST client ─────────────────────────────────────────────────
//
// Direct fetch against PostgREST — avoids @supabase/ssr client hangs.
// Every function: returns data or throws on failure.

const supabaseUrl = () => process.env.NEXT_PUBLIC_SUPABASE_URL!
const supabaseAnonKey = () => process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY!

function getAuthHeaders(): Record<string, string> {
  const token = useAuthStore.getState().session?.access_token
  const headers: Record<string, string> = {
    'Content-Type': 'application/json',
    'apikey': supabaseAnonKey(),
  }
  if (token) {
    headers['Authorization'] = `Bearer ${token}`
  }
  return headers
}

function getAuthUserId(): string | null {
  const state = useAuthStore.getState()
  return state.userId ?? state.user?.id ?? null
}

/** Upsert a row. Throws on any failure. */
async function dbUpsert(table: string, row: Record<string, unknown>, onConflict: string): Promise<void> {
  const resp = await fetch(`${supabaseUrl()}/rest/v1/${table}?on_conflict=${onConflict}`, {
    method: 'POST',
    headers: { ...getAuthHeaders(), 'Prefer': 'resolution=merge-duplicates' },
    body: JSON.stringify(row),
  })
  if (!resp.ok) {
    const body = await resp.text()
    throw new Error(`Supabase upsert failed (HTTP ${resp.status}): ${body}`)
  }
}

/** Select rows. Throws on failure, returns parsed JSON array. */
async function dbSelect<T = any>(table: string, params: string): Promise<T[]> {
  const resp = await fetch(`${supabaseUrl()}/rest/v1/${table}?${params}`, {
    method: 'GET',
    headers: getAuthHeaders(),
  })
  if (!resp.ok) {
    const body = await resp.text()
    throw new Error(`Supabase select failed (HTTP ${resp.status}): ${body}`)
  }
  return resp.json()
}

/** Delete rows. Throws on failure. */
async function dbDelete(table: string, params: string): Promise<void> {
  const resp = await fetch(`${supabaseUrl()}/rest/v1/${table}?${params}`, {
    method: 'DELETE',
    headers: getAuthHeaders(),
  })
  if (!resp.ok) {
    const body = await resp.text()
    throw new Error(`Supabase delete failed (HTTP ${resp.status}): ${body}`)
  }
}

// ─── Supabase save/delete helpers (single source of truth) ────────────────

/**
 * Save an API key to Supabase.
 * - If not authenticated at all: silently skips (local-only save).
 * - If authenticated but userId/session is missing: throws so UI can show the error.
 * - If Supabase errors: throws.
 */
async function saveKeyToSupabase(providerKey: string, rawKey: string): Promise<void> {
  const authState = useAuthStore.getState()
  const userId = authState.userId ?? authState.user?.id ?? null

  // DEBUG — remove after diagnosing save issue
  console.log('[CASSIA-DEBUG] saveKeyToSupabase:', {
    providerKey,
    userId,
    isAuthenticated: authState.isAuthenticated,
    hasSession: !!authState.session,
    hasToken: !!authState.session?.access_token,
    hasUser: !!authState.user,
    storeUserId: authState.userId,
    userObjId: authState.user?.id,
  })

  if (!userId) {
    if (authState.isAuthenticated) {
      throw new Error('Session error: signed in but user ID is unavailable. Please sign out and sign in again.')
    }
    console.log('[CASSIA-DEBUG] saveKeyToSupabase: SKIPPED (not authenticated)')
    return // Not authenticated — local-only save
  }

  if (!authState.session?.access_token) {
    throw new Error('Session expired. Please sign out and sign in again.')
  }

  console.log('[CASSIA-DEBUG] saveKeyToSupabase: calling dbUpsert...')
  await dbUpsert('user_api_keys', {
    user_id: userId,
    provider: providerKey,
    encrypted_key: encryptKey(rawKey),
  }, 'user_id,provider')
  console.log('[CASSIA-DEBUG] saveKeyToSupabase: SUCCESS')
}

/**
 * Delete an API key from Supabase.
 * Same auth-check logic as saveKeyToSupabase.
 */
async function deleteKeyFromSupabase(providerKey: string): Promise<void> {
  const authState = useAuthStore.getState()
  const userId = authState.userId ?? authState.user?.id ?? null

  if (!userId) {
    if (authState.isAuthenticated) {
      throw new Error('Session error: signed in but user ID is unavailable. Please sign out and sign in again.')
    }
    return
  }

  if (!authState.session?.access_token) {
    throw new Error('Session expired. Please sign out and sign in again.')
  }

  await dbDelete(
    'user_api_keys',
    `user_id=eq.${userId}&provider=eq.${encodeURIComponent(providerKey)}`
  )
}

// ─── Store interface ──────────────────────────────────────────────────────

interface ApiKeyState {
  apiKeys: { openrouter: string; anthropic: string; openai: string }
  customProviders: Record<CustomPresetKey, CustomProviderConfig>
  selectedCustomPreset: CustomPresetKey
  provider: Provider
  model: string
  reasoningEffort: ReasoningEffort | null
  isLoading: boolean
  error: string | null

  // Core actions
  setApiKey: (key: string, provider?: Provider) => Promise<void>
  loadApiKeys: () => Promise<void>
  clearApiKey: (provider?: Provider) => Promise<void>

  // Custom provider actions
  setCustomProviderKey: (preset: CustomPresetKey, key: string) => Promise<void>
  setCustomProviderBaseUrl: (preset: CustomPresetKey, url: string) => void
  setSelectedCustomPreset: (preset: CustomPresetKey) => void
  clearCustomProviderKey: (preset: CustomPresetKey) => Promise<void>

  // Custom provider getters
  getCustomApiKey: (preset?: CustomPresetKey) => string
  getCustomBaseUrl: (preset?: CustomPresetKey) => string

  // UI state
  setProvider: (provider: Provider) => void
  setModel: (model: string) => void
  setReasoningEffort: (effort: ReasoningEffort | null) => void
  setCustomBaseUrl: (url: string) => void
  clearError: () => void

  // Helpers
  getApiKey: (provider?: Provider) => string
  apiKey: string
}

// ─── Store implementation ─────────────────────────────────────────────────

export const useApiKeyStore = create<ApiKeyState>()(
  persist(
    (set, get) => ({
      apiKeys: { openrouter: '', anthropic: '', openai: '' },
      customProviders: { ...DEFAULT_CUSTOM_PROVIDERS },
      selectedCustomPreset: 'deepseek',
      provider: 'openrouter',
      model: modelSettings.providers.openrouter.default_model,
      reasoningEffort: getDefaultReasoningEffort('openrouter', modelSettings.providers.openrouter.default_model),
      isLoading: false,
      error: null,

      // ── UI setters ────────────────────────────────────────────────

      clearError: () => set({ error: null }),

      setProvider: (provider) => set({ provider }),

      setModel: (model) => {
        set({ model, reasoningEffort: getDefaultReasoningEffort(get().provider, model) })
      },

      setReasoningEffort: (effort) => set({ reasoningEffort: effort }),

      setCustomBaseUrl: (url) => {
        const preset = get().selectedCustomPreset
        set((s) => ({
          customProviders: { ...s.customProviders, [preset]: { ...s.customProviders[preset], baseUrl: url } }
        }))
      },

      setSelectedCustomPreset: (preset) => set({ selectedCustomPreset: preset }),

      setCustomProviderBaseUrl: (preset, url) => {
        set((s) => ({
          customProviders: { ...s.customProviders, [preset]: { ...s.customProviders[preset], baseUrl: url } }
        }))
      },

      // ── Custom provider getters ───────────────────────────────────

      getCustomApiKey: (preset?) => {
        const s = get()
        return s.customProviders[preset || s.selectedCustomPreset]?.apiKey || ''
      },

      getCustomBaseUrl: (preset?) => {
        const s = get()
        return s.customProviders[preset || s.selectedCustomPreset]?.baseUrl || ''
      },

      getApiKey: (provider?) => {
        const s = get()
        const p = provider || s.provider
        if (p === 'custom') return s.customProviders[s.selectedCustomPreset]?.apiKey || ''
        return s.apiKeys[p as keyof typeof s.apiKeys] || ''
      },

      apiKey: '',

      // ── Save ──────────────────────────────────────────────────────

      setApiKey: async (key, provider?) => {
        const targetProvider = provider || get().provider

        if (targetProvider === 'custom') {
          await get().setCustomProviderKey(get().selectedCustomPreset, key)
          return
        }

        // Update local state
        set((s) => ({ apiKeys: { ...s.apiKeys, [targetProvider]: key } }))

        // Persist to Supabase (throws on failure)
        await saveKeyToSupabase(targetProvider, key)
      },

      setCustomProviderKey: async (preset, key) => {
        // Update local state
        set((s) => ({
          customProviders: { ...s.customProviders, [preset]: { ...s.customProviders[preset], apiKey: key } }
        }))

        // Persist to Supabase (throws on failure)
        await saveKeyToSupabase(`custom:${preset}`, key)
      },

      // ── Load ──────────────────────────────────────────────────────

      loadApiKeys: async () => {
        const authState = useAuthStore.getState()
        const userId = authState.userId ?? authState.user?.id ?? null

        if (!userId) {
          if (authState.isAuthenticated) {
            set({ error: 'Session error: signed in but user ID is unavailable. Please sign out and sign in again.' })
          }
          return
        }

        set({ isLoading: true, error: null })

        try {
          const rows = await dbSelect<{ provider: string; encrypted_key: string }>(
            'user_api_keys',
            `select=provider,encrypted_key&user_id=eq.${userId}`
          )

          const loadedKeys = { openrouter: '', anthropic: '', openai: '' }
          const loadedCustom = { ...DEFAULT_CUSTOM_PROVIDERS }

          for (const row of rows) {
            const p = row.provider
            const decrypted = decryptKey(row.encrypted_key)

            if (p.startsWith('custom:')) {
              const preset = p.replace('custom:', '') as CustomPresetKey
              if (preset in loadedCustom) {
                loadedCustom[preset] = { ...loadedCustom[preset], apiKey: decrypted }
              }
            } else if (p === 'custom') {
              // Legacy migration: bare 'custom' → manual preset
              loadedCustom.manual = { ...loadedCustom.manual, apiKey: decrypted }
            } else if (p in loadedKeys) {
              loadedKeys[p as keyof typeof loadedKeys] = decrypted
            }
          }

          set({ apiKeys: loadedKeys, customProviders: loadedCustom })

          const hasAny =
            Object.values(loadedKeys).some((k) => k !== '') ||
            Object.values(loadedCustom).some((c) => c.apiKey !== '')

          if (!hasAny) {
            set({ error: 'No API keys found in your account' })
          }
        } catch (err) {
          const msg = err instanceof Error ? err.message : 'Failed to load API keys'
          set({ error: msg })
          throw err
        } finally {
          set({ isLoading: false })
        }
      },

      // ── Delete ────────────────────────────────────────────────────

      clearApiKey: async (provider?) => {
        const targetProvider = provider || get().provider

        if (targetProvider === 'custom') {
          await get().clearCustomProviderKey(get().selectedCustomPreset)
          return
        }

        set((s) => ({ apiKeys: { ...s.apiKeys, [targetProvider]: '' } }))
        await deleteKeyFromSupabase(targetProvider)
      },

      clearCustomProviderKey: async (preset) => {
        set((s) => ({
          customProviders: { ...s.customProviders, [preset]: { ...s.customProviders[preset], apiKey: '' } }
        }))
        await deleteKeyFromSupabase(`custom:${preset}`)
      },
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
        apiKeys: state.apiKeys,
      }),
      migrate: (persistedState: unknown, version: number) => {
        const state = persistedState as Record<string, unknown>

        if (version < 2) {
          const oldApiKeys = state.apiKeys as Record<string, string> | undefined
          const oldCustomBaseUrl = state.customBaseUrl as string | undefined
          const newCustomProviders = { ...DEFAULT_CUSTOM_PROVIDERS }

          if (oldApiKeys?.custom) {
            newCustomProviders.manual = { apiKey: oldApiKeys.custom, baseUrl: oldCustomBaseUrl || '' }
          }

          return {
            ...state,
            apiKeys: {
              openrouter: oldApiKeys?.openrouter || '',
              anthropic: oldApiKeys?.anthropic || '',
              openai: oldApiKeys?.openai || '',
            },
            customProviders: newCustomProviders,
            selectedCustomPreset: oldApiKeys?.custom ? 'manual' : 'deepseek',
          }
        }

        return state
      },
    }
  )
)
