'use client'

import { useEffect, ReactNode, useRef } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { useResultsStore } from '@/lib/stores/results-store'
import { cleanupOldStorage } from '@/lib/utils/storage-cleanup'

interface AuthProviderProps {
  children: ReactNode
}

export function AuthProvider({ children }: AuthProviderProps) {
  const { initialize, isAuthenticated } = useAuthStore()
  const { loadResults, loadSessions } = useResultsStore()
  const initialized = useRef(false)
  
  useEffect(() => {
    // Prevent re-initialization on hot reloads
    if (!initialized.current) {
      initialized.current = true

      // Clean up old storage format to prevent migration issues
      cleanupOldStorage()
      
      // Initialize auth on app start
      initialize()
    }
  }, [initialize])
  
  useEffect(() => {
    // Load user data when authentication state changes
    if (isAuthenticated) {
      console.log('User authenticated in AuthProvider, loading user data...')
      // Load user's results and sessions
      loadResults()
      loadSessions()
      // Note: API keys are now loaded manually via LoadApiKeysButton
    } else {
      console.log('User not authenticated in AuthProvider')
    }
  }, [isAuthenticated, loadResults, loadSessions])
  
  return <>{children}</>
}