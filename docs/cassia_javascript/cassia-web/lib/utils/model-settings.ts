/**
 * Lazy loader for model settings with caching
 * Improves FCP by deferring JSON parsing until needed
 */

export interface ModelInfo {
  actual_name: string
  aliases: string[]
  description: string
  recommended: boolean
  cost_tier: string
  context_window: number
  use_cases: string[]
}

export interface ProviderConfig {
  name: string
  base_url: string
  default_model: string
  recommended_model: string
  description?: string
  models: Record<string, ModelInfo>
}

export interface ModelSettings {
  version: string
  last_updated: string
  providers: Record<string, ProviderConfig>
  provider_shortcuts: Record<string, Record<string, string>>
  use_case_recommendations: Record<string, {
    best: string
    alternatives: string[]
    description: string
  }>
  cost_tiers: Record<string, {
    description: string
    models: string[]
  }>
}

// In-memory cache
let cachedSettings: ModelSettings | null = null
let loadPromise: Promise<ModelSettings> | null = null

/**
 * Load model settings asynchronously with caching
 * First call fetches from JSON, subsequent calls return cached data
 */
export async function loadModelSettings(): Promise<ModelSettings> {
  // Return cached settings if available
  if (cachedSettings) {
    return cachedSettings
  }

  // Return existing promise if already loading
  if (loadPromise) {
    return loadPromise
  }

  // Start loading
  loadPromise = fetch('/examples/model_settings.json')
    .then(res => {
      if (!res.ok) throw new Error('Failed to load model settings')
      return res.json()
    })
    .then(data => {
      cachedSettings = data as ModelSettings
      return cachedSettings
    })
    .catch(err => {
      loadPromise = null // Allow retry on error
      throw err
    })

  return loadPromise
}

/**
 * Get default model for a provider
 * Falls back to OpenRouter default if provider not found
 */
export function getDefaultModel(settings: ModelSettings, provider: string): string {
  const providerData = settings.providers[provider]
  return providerData?.default_model || settings.providers.openrouter?.default_model || ''
}

/**
 * Check if settings are already cached (for SSR/initial render)
 */
export function isSettingsCached(): boolean {
  return cachedSettings !== null
}

/**
 * Get cached settings synchronously (returns null if not loaded)
 * Useful for initial render before async load completes
 */
export function getCachedSettings(): ModelSettings | null {
  return cachedSettings
}

/**
 * Preload settings - can be called early to start loading
 */
export function preloadModelSettings(): void {
  loadModelSettings().catch(() => {
    // Silently ignore preload errors
  })
}
