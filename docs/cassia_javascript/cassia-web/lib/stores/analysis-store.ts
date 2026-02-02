import { create } from 'zustand'

interface AnalysisState {
  // Current analysis state
  isRunning: boolean
  progress: number
  currentStep: string
  logs: string[]
  abortController: AbortController | null

  // File data
  uploadedFile: File | null
  fileData: any[] | null
  fileMetadata: any | null

  // Results
  results: any | null
  downloadLinks: string[]

  // Actions
  setFile: (file: File, data: any[], metadata?: any) => void
  startAnalysis: () => AbortController
  stopAnalysis: () => void
  updateProgress: (progress: number, step: string) => void
  addLog: (log: string) => void
  setResults: (results: any) => void
  addDownloadLink: (link: string) => void
  reset: () => void
}

export const useAnalysisStore = create<AnalysisState>((set, get) => ({
  // Initial state
  isRunning: false,
  progress: 0,
  currentStep: '',
  logs: [],
  abortController: null,
  uploadedFile: null,
  fileData: null,
  fileMetadata: null,
  results: null,
  downloadLinks: [],
  
  // Actions
  setFile: (file: File, data: any[], metadata?: any) => set({
    uploadedFile: file,
    fileData: data,
    fileMetadata: metadata,
    results: null,
    downloadLinks: [],
    logs: [`📁 Loaded file: ${file.name} (${data.length} rows)`],
  }),
  
  startAnalysis: () => {
    const controller = new AbortController()
    set({
      isRunning: true,
      progress: 0,
      currentStep: 'Initializing analysis...',
      results: null,
      downloadLinks: [],
      abortController: controller,
    })
    return controller
  },

  stopAnalysis: () => {
    const { abortController } = get()
    if (abortController) {
      abortController.abort()
    }
    set((state) => ({
      isRunning: false,
      currentStep: 'Analysis stopped by user',
      abortController: null,
      logs: [...state.logs, '🛑 Analysis stopped by user'],
    }))
  },
  
  updateProgress: (progress: number, step: string) => set((state) => ({
    progress,
    currentStep: step,
    logs: [...state.logs, `⏳ ${step}`],
  })),
  
  addLog: (log: string) => set((state) => ({
    logs: [...state.logs, log],
  })),
  
  setResults: (results: any) => set({
    results,
    isRunning: false,
    progress: 100,
    currentStep: 'Analysis complete!',
    abortController: null,
  }),
  
  addDownloadLink: (link: string) => set((state) => ({
    downloadLinks: [...state.downloadLinks, link],
  })),
  
  reset: () => set({
    isRunning: false,
    progress: 0,
    currentStep: '',
    logs: [],
    abortController: null,
    uploadedFile: null,
    fileData: null,
    fileMetadata: null,
    results: null,
    downloadLinks: [],
  }),
}))