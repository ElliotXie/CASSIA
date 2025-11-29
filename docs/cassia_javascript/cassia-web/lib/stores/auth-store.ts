import { create } from 'zustand'
import { persist } from 'zustand/middleware'
import { User, Session, AuthError } from '@supabase/supabase-js'
import { createClient } from '../../utils/supabase/client'
import { Profile } from '../supabase/types'

interface AuthState {
  user: User | null
  session: Session | null
  profile: Profile | null
  isLoading: boolean
  error: string | null
  isAuthenticated: boolean
  userId: string | null
  
  // Auth actions
  signIn: (email: string, password: string) => Promise<void>
  signUp: (email: string, password: string, fullName?: string) => Promise<void>
  signOut: () => Promise<void>
  resetPassword: (email: string) => Promise<void>
  updateProfile: (updates: Partial<Profile>) => Promise<void>
  
  // Session management
  initialize: () => Promise<void>
  refreshSession: () => Promise<void>
  loadProfile: () => Promise<void>
  clearError: () => void
  forceResetLoading: () => void
}

export const useAuthStore = create<AuthState>()(
  persist(
    (set, get) => ({
      user: null,
      session: null,
      profile: null,
      isLoading: false,
      error: null,
      isAuthenticated: false,
      userId: null,
      
      clearError: () => set({ error: null }),
      forceResetLoading: () => {
        console.log('Auth store: Force resetting loading state')
        set({ isLoading: false, error: null })
      },
      
      signIn: async (email: string, password: string) => {
        const supabase = createClient()
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase.auth.signInWithPassword({
            email,
            password,
          })
          
          if (error) throw error
          
          set({ 
            user: data.user, 
            session: data.session,
            isAuthenticated: !!data.user,
            userId: data.user?.id || null
          })
          
          // Load profile after successful sign in
          await get().loadProfile()
        } catch (error) {
          const authError = error as AuthError
          let userFriendlyMessage = authError.message
          
          // Make error messages more user-friendly
          if (authError.message.includes('Invalid login credentials')) {
            userFriendlyMessage = 'Invalid email or password. Please check your credentials and try again.'
          } else if (authError.message.includes('Email not confirmed')) {
            userFriendlyMessage = 'Please check your email and click the confirmation link before signing in.'
          } else if (authError.message.includes('User not found')) {
            userFriendlyMessage = 'No account found with this email address.'
          } else if (authError.message.includes('Too many requests')) {
            userFriendlyMessage = 'Too many failed attempts. Please wait a moment before trying again.'
          }
          
          set({ error: userFriendlyMessage })
        } finally {
          set({ isLoading: false })
        }
      },
      
      signUp: async (email: string, password: string, fullName?: string) => {
        const supabase = createClient()
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase.auth.signUp({
            email,
            password,
            options: {
              data: {
                full_name: fullName,
              },
            },
          })
          
          if (error) throw error
          
          set({ 
            user: data.user, 
            session: data.session,
            isAuthenticated: !!data.user,
            userId: data.user?.id || null
          })
        } catch (error) {
          const authError = error as AuthError
          let userFriendlyMessage = authError.message
          
          // Make error messages more user-friendly
          if (authError.message.includes('User already registered')) {
            userFriendlyMessage = 'An account with this email already exists. Please sign in instead.'
          } else if (authError.message.includes('Password should be at least')) {
            userFriendlyMessage = 'Password must be at least 6 characters long.'
          } else if (authError.message.includes('Invalid email')) {
            userFriendlyMessage = 'Please enter a valid email address.'
          } else if (authError.message.includes('Signup is disabled')) {
            userFriendlyMessage = 'Account creation is temporarily disabled. Please contact support.'
          }
          
          set({ error: userFriendlyMessage })
        } finally {
          set({ isLoading: false })
        }
      },
      
      signOut: async () => {
        const supabase = createClient()
        set({ isLoading: true, error: null })
        
        try {
          const { error } = await supabase.auth.signOut()
          if (error) throw error
          
          set({ 
            user: null, 
            session: null, 
            profile: null,
            isAuthenticated: false,
            userId: null
          })
        } catch (error) {
          set({ error: (error as AuthError).message })
        } finally {
          set({ isLoading: false })
        }
      },
      
      resetPassword: async (email: string) => {
        const supabase = createClient()
        set({ isLoading: true, error: null })
        
        try {
          const { error } = await supabase.auth.resetPasswordForEmail(email, {
            redirectTo: `${window.location.origin}/auth/reset-password`,
          })
          
          if (error) throw error
        } catch (error) {
          set({ error: (error as AuthError).message })
        } finally {
          set({ isLoading: false })
        }
      },
      
      updateProfile: async (updates: Partial<Profile>) => {
        const supabase = createClient()
        const { userId } = get()
        
        if (!userId) throw new Error('User not authenticated')
        
        set({ isLoading: true, error: null })
        
        try {
          const { data, error } = await supabase
            .from('profiles')
            .update(updates)
            .eq('id', userId)
            .select()
            .single()
          
          if (error) throw error
          
          set({ profile: data })
        } catch (error) {
          set({ error: (error as Error).message })
        } finally {
          set({ isLoading: false })
        }
      },
      
      loadProfile: async () => {
        const supabase = createClient()
        const { userId } = get()
        
        if (!userId) return
        
        try {
          const { data, error } = await supabase
            .from('profiles')
            .select('*')
            .eq('id', userId)
            .single()
          
          if (error && error.code !== 'PGRST116') throw error
          
          set({ profile: data })
        } catch (error) {
          console.error('Error loading profile:', error)
        }
      },
      
      initialize: async () => {
        const supabase = createClient()
        console.log('Auth store: Starting initialization')
        set({ isLoading: true })
        
        try {
          console.log('Auth store: Getting session from Supabase')
          const { data: { session }, error } = await supabase.auth.getSession()
          
          if (error) {
            console.error('Auth store: Error getting session:', error)
            throw error
          }
          
          console.log('Auth store: Session obtained', session ? 'User logged in' : 'No session')
          
          set({
            user: session?.user || null,
            session: session || null,
            isAuthenticated: !!session?.user,
            userId: session?.user?.id || null
          })
          
          if (session?.user) {
            console.log('Auth store: Loading profile for user:', session.user.id)
            try {
              await get().loadProfile()
              console.log('Auth store: Profile loaded successfully')
            } catch (profileError) {
              console.error('Auth store: Error loading profile:', profileError)
              // Don't throw - profile loading failure shouldn't block initialization
            }
          }
          
          // Listen for auth changes
          supabase.auth.onAuthStateChange(async (event, session) => {
            console.log('Auth store: Auth state changed:', event)
            set({
              user: session?.user || null,
              session: session || null,
              isAuthenticated: !!session?.user,
              userId: session?.user?.id || null
            })
            
            if (event === 'SIGNED_IN' && session?.user) {
              try {
                await get().loadProfile()
              } catch (profileError) {
                console.error('Auth store: Error loading profile on sign in:', profileError)
              }
            } else if (event === 'SIGNED_OUT') {
              set({ profile: null })
            }
          })
          
          console.log('Auth store: Initialization completed successfully')
        } catch (error) {
          console.error('Auth store: Initialization error:', error)
          set({ error: (error as Error).message })
        } finally {
          console.log('Auth store: Setting isLoading to false')
          set({ isLoading: false })
        }
      },
      
      refreshSession: async () => {
        const supabase = createClient()
        
        try {
          const { data: { session }, error } = await supabase.auth.refreshSession()
          
          if (error) throw error
          
          set({
            user: session?.user || null,
            session: session || null,
            isAuthenticated: !!session?.user,
            userId: session?.user?.id || null
          })
        } catch (error) {
          set({ error: (error as Error).message })
        }
      },
    }),
    {
      name: 'cassia-auth-storage',
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
        // Only persist minimal state - session will be handled by Supabase
        // Never persist isLoading or error states
        user: state.user,
        session: state.session,
        profile: state.profile,
        isAuthenticated: state.isAuthenticated,
        userId: state.userId,
      }),
      version: 2, // Bumped to reset persisted state after isAuthenticated changes
      skipHydration: false,
      migrate: (persistedState: any, version: number) => {
        // Handle migration from any version to version 2
        if (version === 0 || version === 1) {
          console.log(`Migrating auth store from version ${version} to 2`)
          const oldState = persistedState as any
          return {
            user: oldState.user || null,
            session: oldState.session || null,
            profile: oldState.profile || null,
            isAuthenticated: oldState.isAuthenticated || false,
            userId: oldState.userId || null,
            isLoading: false,
            error: null
          }
        }
        return persistedState
      },
      onRehydrateStorage: () => (state) => {
        console.log('Auth store: Hydration callback triggered')
        if (state) {
          console.log('Auth store: Hydrated successfully, force resetting loading state')
          // Always reset loading state after hydration
          state.isLoading = false
          state.error = null
        }
        // Also set a timeout to force reset loading state
        setTimeout(() => {
          console.log('Auth store: Timeout force reset loading state')
          useAuthStore.getState().forceResetLoading()
        }, 100)
      },
    }
  )
)