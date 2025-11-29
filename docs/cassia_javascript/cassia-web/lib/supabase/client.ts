import { createBrowserClient } from '@supabase/ssr'
import { Database } from './types'
import { getSupabaseConfig } from './config'

export const createClient = () => {
  try {
    const { url, anonKey } = getSupabaseConfig()
    return createBrowserClient<Database>(url, anonKey)
  } catch (error) {
    console.error('Supabase client creation failed:', error)
    throw error
  }
}