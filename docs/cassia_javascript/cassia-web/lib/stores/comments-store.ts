import { create } from 'zustand'
import { createClient } from '../../utils/supabase/client'
import { useAuthStore } from './auth-store'
import { FeedbackComment, CommentType } from '../supabase/types'

interface CommentsState {
  // State
  comments: FeedbackComment[]
  userComments: FeedbackComment[]
  isLoading: boolean
  isSubmitting: boolean
  error: string | null
  successMessage: string | null

  // Actions
  loadApprovedComments: () => Promise<void>
  loadUserComments: () => Promise<void>
  submitComment: (data: {
    comment_type: CommentType
    name?: string
    email?: string
    message: string
  }) => Promise<{ success: boolean; comment?: FeedbackComment }>
  clearError: () => void
  clearSuccess: () => void
}

export const useCommentsStore = create<CommentsState>()((set, get) => ({
  comments: [],
  userComments: [],
  isLoading: false,
  isSubmitting: false,
  error: null,
  successMessage: null,

  clearError: () => set({ error: null }),
  clearSuccess: () => set({ successMessage: null }),

  loadApprovedComments: async () => {
    const supabase = createClient()
    set({ isLoading: true, error: null })

    try {
      const { data, error } = await supabase
        .from('feedback_comments')
        .select('*')
        .eq('is_approved', true)
        .order('created_at', { ascending: false })
        .limit(50)

      if (error) throw error

      set({ comments: data || [] })
    } catch (error) {
      set({ error: (error as Error).message })
      console.error('Error loading comments:', error)
    } finally {
      set({ isLoading: false })
    }
  },

  loadUserComments: async () => {
    const supabase = createClient()
    const authStore = useAuthStore.getState()

    if (!authStore.isAuthenticated || !authStore.userId) {
      set({ userComments: [] })
      return
    }

    try {
      const { data, error } = await supabase
        .from('feedback_comments')
        .select('*')
        .eq('user_id', authStore.userId)
        .order('created_at', { ascending: false })

      if (error) throw error

      set({ userComments: data || [] })
    } catch (error) {
      console.error('Error loading user comments:', error)
    }
  },

  submitComment: async (data) => {
    const supabase = createClient()
    const authStore = useAuthStore.getState()

    if (!authStore.isAuthenticated || !authStore.userId) {
      set({ error: 'You must be logged in to submit feedback' })
      return { success: false }
    }

    set({ isSubmitting: true, error: null, successMessage: null })

    try {
      const { data: comment, error } = await supabase
        .from('feedback_comments')
        .insert({
          user_id: authStore.userId,
          comment_type: data.comment_type,
          name: data.name || null,
          email: data.email || null,
          message: data.message
        })
        .select()
        .single()

      if (error) throw error

      // Add to user comments
      set((state) => ({
        userComments: [comment, ...state.userComments],
        successMessage: 'Feedback submitted successfully! It will be visible after admin approval.'
      }))

      return { success: true, comment }
    } catch (error) {
      set({ error: (error as Error).message })
      return { success: false }
    } finally {
      set({ isSubmitting: false })
    }
  }
}))
