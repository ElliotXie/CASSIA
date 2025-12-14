'use client'

import { useEffect, ReactNode, useRef } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { cleanupOldStorage } from '@/lib/utils/storage-cleanup'

interface AuthProviderProps {
  children: ReactNode
}

export function AuthProvider({ children }: AuthProviderProps) {
  const { initialize } = useAuthStore()
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

  return <>{children}</>
}
