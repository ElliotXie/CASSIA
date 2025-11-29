import { create } from 'zustand'

interface AnalysisState {
  // Current analysis state
  isRunning: boolean
  progress: number
  currentStep: string
  logs: string[]
  
  // File data
  uploadedFile: File | null
  fileData: any[] | null
  fileMetadata: any | null
  
  // Results
  results: any | null
  downloadLinks: string[]
  
  // Actions
  setFile: (file: File, data: any[], metadata?: any) => void
  startAnalysis: () => void
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
  
  startAnalysis: () => set({
    isRunning: true,
    progress: 0,
    currentStep: 'Initializing analysis...',
    results: null,
    downloadLinks: [],
  }),
  
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
  }),
  
  addDownloadLink: (link: string) => set((state) => ({
    downloadLinks: [...state.downloadLinks, link],
  })),
  
  reset: () => set({
    isRunning: false,
    progress: 0,
    currentStep: '',
    logs: [],
    uploadedFile: null,
    fileData: null,
    fileMetadata: null,
    results: null,
    downloadLinks: [],
  }),
}))