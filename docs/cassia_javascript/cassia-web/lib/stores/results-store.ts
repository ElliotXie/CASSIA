import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { createClient } from '../../utils/supabase/client'
import { useAuthStore } from './auth-store'
import { AnalysisResult, AnalysisSession, AnalysisType, SessionStatus } from '../supabase/types'

interface ResultsState {
  // State
  results: AnalysisResult[]
  sessions: AnalysisSession[]
  currentSession: AnalysisSession | null
  isLoading: boolean
  error: string | null
  
  // Result actions
  saveResult: (result: Omit<AnalysisResult, 'id' | 'user_id' | 'created_at' | 'updated_at'>) => Promise<string>
  loadResults: (limit?: number) => Promise<void>
  deleteResult: (id: string) => Promise<void>
  updateResult: (id: string, updates: Partial<AnalysisResult>) => Promise<void>
  
  // Session actions
  createSession: (session: Omit<AnalysisSession, 'id' | 'user_id' | 'created_at' | 'updated_at'>) => Promise<string>
  updateSession: (id: string, updates: Partial<AnalysisSession>) => Promise<void>
  loadSessions: (status?: SessionStatus) => Promise<void>
  deleteSession: (id: string) => Promise<void>
  setCurrentSession: (session: AnalysisSession | null) => void
  
  // Utility actions
  clearError: () => void
  getResultsByType: (type: AnalysisType) => AnalysisResult[]
  getSessionsByType: (type: AnalysisType) => AnalysisSession[]
  
  // Export/Import
  exportResults: (resultIds: string[]) => Promise<Blob>
  importResults: (data: any[]) => Promise<void>
}

export const useResultsStore = create<ResultsState>()(
  persist(
    (set, get) => ({
      results: [],
      sessions: [],
      currentSession: null,
      isLoading: false,
      error: null,
      
      clearError: () => set({ error: null }),
      
      saveResult: async (result) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to save results')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('analysis_results')
            .insert({
              ...result,
              user_id: authStore.userId
            })
            .select()
            .single()
          
          if (error) throw error
          
          // Add to local state
          set((state) => ({
            results: [data, ...state.results]
          }))
          
          console.log('Analysis result saved successfully')
          return data.id
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      loadResults: async (limit = 50) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          set({ results: [] })
          return
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('analysis_results')
            .select('*')
            .eq('user_id', authStore.userId)
            .order('created_at', { ascending: false })
            .limit(limit)
          
          if (error) throw error
          
          set({ results: data || [] })
          console.log(`Loaded ${data?.length || 0} analysis results`)
        } catch (error) {
          set({ error: (error as Error).message })
          console.error('Error loading results:', error)
        } finally {
          set({ isLoading: false })
        }
      },
      
      deleteResult: async (id: string) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to delete results')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { error } = await supabase
            .from('analysis_results')
            .delete()
            .eq('id', id)
            .eq('user_id', authStore.userId)
          
          if (error) throw error
          
          // Remove from local state
          set((state) => ({
            results: state.results.filter(r => r.id !== id)
          }))
          
          console.log('Analysis result deleted successfully')
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      updateResult: async (id: string, updates: Partial<AnalysisResult>) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to update results')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('analysis_results')
            .update(updates)
            .eq('id', id)
            .eq('user_id', authStore.userId)
            .select()
            .single()
          
          if (error) throw error
          
          // Update local state
          set((state) => ({
            results: state.results.map(r => r.id === id ? data : r)
          }))
          
          console.log('Analysis result updated successfully')
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      createSession: async (session) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to create sessions')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('analysis_sessions')
            .insert({
              ...session,
              user_id: authStore.userId
            })
            .select()
            .single()
          
          if (error) throw error
          
          // Add to local state
          set((state) => ({
            sessions: [data, ...state.sessions],
            currentSession: data
          }))
          
          console.log('Analysis session created successfully')
          return data.id
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      updateSession: async (id: string, updates: Partial<AnalysisSession>) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to update sessions')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('analysis_sessions')
            .update(updates)
            .eq('id', id)
            .eq('user_id', authStore.userId)
            .select()
            .single()
          
          if (error) throw error
          
          // Update local state
          set((state) => ({
            sessions: state.sessions.map(s => s.id === id ? data : s),
            currentSession: state.currentSession?.id === id ? data : state.currentSession
          }))
          
          console.log('Analysis session updated successfully')
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      loadSessions: async (status?: SessionStatus) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          set({ sessions: [] })
          return
        }
        
        set({ isLoading: true, error: null })
        
        try {
          let query = supabase
            .from('analysis_sessions')
            .select('*')
            .eq('user_id', authStore.userId)
            .order('created_at', { ascending: false })
          
          if (status) {
            query = query.eq('status', status)
          }
          
          const { data, error } = await query
          
          if (error) throw error
          
          set({ sessions: data || [] })
          console.log(`Loaded ${data?.length || 0} analysis sessions`)
        } catch (error) {
          set({ error: (error as Error).message })
          console.error('Error loading sessions:', error)
        } finally {
          set({ isLoading: false })
        }
      },
      
      deleteSession: async (id: string) => {
        const supabase = createClient()
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to delete sessions')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const { error } = await supabase
            .from('analysis_sessions')
            .delete()
            .eq('id', id)
            .eq('user_id', authStore.userId)
          
          if (error) throw error
          
          // Remove from local state
          set((state) => ({
            sessions: state.sessions.filter(s => s.id !== id),
            currentSession: state.currentSession?.id === id ? null : state.currentSession
          }))
          
          console.log('Analysis session deleted successfully')
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      },
      
      setCurrentSession: (session: AnalysisSession | null) => {
        set({ currentSession: session })
      },
      
      getResultsByType: (type: AnalysisType) => {
        const state = get()
        return state.results.filter(r => r.analysis_type === type)
      },
      
      getSessionsByType: (type: AnalysisType) => {
        const state = get()
        return state.sessions.filter(s => s.analysis_type === type)
      },
      
      exportResults: async (resultIds: string[]) => {
        const state = get()
        const resultsToExport = state.results.filter(r => resultIds.includes(r.id))
        
        const exportData = {
          exported_at: new Date().toISOString(),
          cassia_version: '2.0',
          results: resultsToExport
        }
        
        const blob = new Blob([JSON.stringify(exportData, null, 2)], {
          type: 'application/json'
        })
        
        return blob
      },
      
      importResults: async (data: any[]) => {
        const authStore = useAuthStore.getState()
        
        if (!authStore.isAuthenticated || !authStore.userId) {
          throw new Error('User must be authenticated to import results')
        }
        
        set({ isLoading: true, error: null })
        
        try {
          const supabase = createClient()
          
          // Prepare data for import
          const resultsToImport = data.map(result => ({
            ...result,
            user_id: authStore.userId,
            id: undefined, // Let Supabase generate new IDs
            created_at: undefined,
            updated_at: undefined
          }))
          
          const { data: importedResults, error } = await supabase
            .from('analysis_results')
            .insert(resultsToImport)
            .select()
          
          if (error) throw error
          
          // Add to local state
          set((state) => ({
            results: [...(importedResults || []), ...state.results]
          }))
          
          console.log(`Imported ${importedResults?.length || 0} analysis results`)
        } catch (error) {
          set({ error: (error as Error).message })
          throw error
        } finally {
          set({ isLoading: false })
        }
      }
    }),
    {
      name: 'cassia-results-storage',
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
        // Only persist current session for continuity
        currentSession: state.currentSession,
        // Keep minimal results for offline access
        results: state.results.slice(0, 10)
      }) as Partial<ResultsState>,
      version: 1,
      skipHydration: false,
      onRehydrateStorage: () => (state) => {
        if (state) {
          console.log('Results storage hydrated successfully')
          // Automatically load results when hydrated
          const authStore = useAuthStore.getState()
          if (authStore.isAuthenticated) {
            setTimeout(() => {
              useResultsStore.getState().loadResults()
              useResultsStore.getState().loadSessions()
            }, 1000)
          }
        }
      }
    }
  )
)