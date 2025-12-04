import { createClient } from '@supabase/supabase-js'

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL!
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY!

export const supabase = createClient(supabaseUrl, supabaseAnonKey)

// Types for our comments table
export type CommentCategory = 'feature_request' | 'bug_report'
export type CommentStatus = 'pending' | 'approved' | 'rejected'

export interface Comment {
  id: string
  created_at: string
  category: CommentCategory
  name: string | null
  email: string | null
  title: string
  content: string
  status: CommentStatus
  admin_notes: string | null
}

export interface NewComment {
  category: CommentCategory
  name?: string
  email?: string
  title: string
  content: string
}
