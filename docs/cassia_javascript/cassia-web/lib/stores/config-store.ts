import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { useApiKeyStore } from './api-key-store'
import modelSettings from '../../public/examples/model_settings.json'
import { ReasoningEffort, getDefaultReasoningEffort } from '../config/model-presets'

export interface PipelineStepConfig {
  provider: string
  model: string
  reasoningEffort: ReasoningEffort | null
}

export interface PipelineModels {
  annotation: PipelineStepConfig
  scoring: PipelineStepConfig
  annotationBoost: PipelineStepConfig
}

interface ConfigState {
  // Analysis Configuration (API config moved to api-key-store)
  model: string // Legacy single model support
  pipelineModels: PipelineModels
  outputName: string
  tissue: string
  species: string
  scoreThreshold: number
  maxWorkers: number
  mergeAnnotations: boolean
  additionalInfo: string
  maxRetries: number
  
  // Actions
  setModel: (model: string) => void // Legacy single model
  setPipelineModel: (step: keyof PipelineModels, provider: string, model: string) => void
  setReasoningEffort: (step: keyof PipelineModels, effort: ReasoningEffort | null) => void
  setAnalysisConfig: (config: Partial<{
    outputName: string
    tissue: string
    species: string
    scoreThreshold: number
    maxWorkers: number
    mergeAnnotations: boolean
    additionalInfo: string
    maxRetries: number
  }>) => void
  reset: () => void
  
  // API key methods that delegate to api-key-store
  getApiKey: () => string
  getProvider: () => string
  setApiKey: (key: string) => void
  setProvider: (provider: string) => void
  setApiModel: (model: string) => void
}

const defaultConfig = {
  model: modelSettings.use_case_recommendations.annotation.best,
  pipelineModels: {
    annotation: {
      provider: 'openrouter',
      model: 'anthropic/claude-haiku-4.5',
      reasoningEffort: 'high' as ReasoningEffort // Anthropic default
    },
    scoring: {
      provider: 'openrouter',
      model: 'google/gemini-2.5-flash',
      reasoningEffort: 'high' as ReasoningEffort // Gemini default
    },
    annotationBoost: {
      provider: 'openrouter',
      model: 'anthropic/claude-haiku-4.5',
      reasoningEffort: 'high' as ReasoningEffort // Anthropic default
    }
  },
  outputName: 'cassia_analysis',
  tissue: 'large_intestine',
  species: 'human',
  scoreThreshold: 75,
  maxWorkers: 4,
  mergeAnnotations: true,
  additionalInfo: '',
  maxRetries: 1,
}

export const useConfigStore = create<ConfigState>()(
  persist(
    (set, get) => ({
      ...defaultConfig,
      
      setModel: (model: string) => set({ model }),
      
      setPipelineModel: (step: keyof PipelineModels, provider: string, model: string) =>
        set((state) => ({
          pipelineModels: {
            ...state.pipelineModels,
            [step]: {
              provider,
              model,
              reasoningEffort: getDefaultReasoningEffort(provider, model)
            }
          }
        })),

      setReasoningEffort: (step: keyof PipelineModels, effort: ReasoningEffort | null) =>
        set((state) => ({
          pipelineModels: {
            ...state.pipelineModels,
            [step]: {
              ...state.pipelineModels[step],
              reasoningEffort: effort
            }
          }
        })),
      
      setAnalysisConfig: (config) => set((state) => ({ ...state, ...config })),
      
      reset: () => set(defaultConfig),
      
      // API key methods that delegate to api-key-store
      getApiKey: () => {
        try {
          return useApiKeyStore.getState().getApiKey()
        } catch (error) {
          console.warn('Error accessing API key store:', error)
          return ''
        }
      },
      getProvider: () => {
        try {
          return useApiKeyStore.getState().provider
        } catch (error) {
          console.warn('Error accessing provider from API key store:', error)
          return 'openrouter'
        }
      },
      setApiKey: (key: string) => {
        try {
          useApiKeyStore.getState().setApiKey(key)
        } catch (error) {
          console.warn('Error setting API key:', error)
        }
      },
      setProvider: (provider: string) => {
        try {
          useApiKeyStore.getState().setProvider(provider as any)
        } catch (error) {
          console.warn('Error setting provider:', error)
        }
      },
      setApiModel: (model: string) => {
        try {
          useApiKeyStore.getState().setModel(model)
        } catch (error) {
          console.warn('Error setting API model:', error)
        }
      }
    }),
    {
      name: 'cassia-config',
      partialize: (state) => ({
        model: state.model,
        pipelineModels: state.pipelineModels,
        outputName: state.outputName,
        tissue: state.tissue,
        species: state.species,
        scoreThreshold: state.scoreThreshold,
        maxWorkers: state.maxWorkers,
        mergeAnnotations: state.mergeAnnotations,
        additionalInfo: state.additionalInfo,
        maxRetries: state.maxRetries,
      }),
      skipHydration: false,
      onRehydrateStorage: () => (state) => {
        if (state) {
          console.log('Config storage hydrated successfully')
          // Migration: Update old peripheral_blood tissue value to large_intestine
          if (state.tissue === 'peripheral_blood') {
            console.log('Migrating tissue from peripheral_blood to large_intestine')
            state.tissue = 'large_intestine'
            // Force save the updated value
            if (typeof window !== 'undefined') {
              const stored = localStorage.getItem('cassia-config')
              if (stored) {
                try {
                  const parsed = JSON.parse(stored)
                  if (parsed.state && parsed.state.tissue === 'peripheral_blood') {
                    parsed.state.tissue = 'large_intestine'
                    localStorage.setItem('cassia-config', JSON.stringify(parsed))
                  }
                } catch (e) {
                  console.error('Error migrating tissue value:', e)
                }
              }
            }
          }
        }
      },
      storage: {
        getItem: (name) => {
          try {
            const item = typeof window !== 'undefined' ? localStorage.getItem(name) : null
            return item ? JSON.parse(item) : null
          } catch (error) {
            console.error('Error reading config from localStorage:', error)
            return null
          }
        },
        setItem: (name, value) => {
          try {
            if (typeof window !== 'undefined') {
              localStorage.setItem(name, JSON.stringify(value))
            }
          } catch (error) {
            console.error('Error writing config to localStorage:', error)
          }
        },
        removeItem: (name) => {
          try {
            if (typeof window !== 'undefined') {
              localStorage.removeItem(name)
            }
          } catch (error) {
            console.error('Error removing config from localStorage:', error)
          }
        }
      }
    }
  )
)