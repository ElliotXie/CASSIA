'use client'

import { useEffect, useState } from 'react'
import { useSearchParams } from 'next/navigation'
import { CheckCircle } from 'lucide-react'

/**
 * Handles displaying email confirmation notification from URL search params.
 * Must be wrapped in Suspense boundary when used.
 */
export function EmailConfirmationHandler() {
  const [showNotification, setShowNotification] = useState(false)
  const searchParams = useSearchParams()

  useEffect(() => {
    if (searchParams.get('confirmed') === 'true') {
      setShowNotification(true)

      // Auto-hide after 5 seconds
      const timer = setTimeout(() => {
        setShowNotification(false)
      }, 5000)

      // Clean up URL without reloading
      window.history.replaceState({}, '', '/')

      return () => clearTimeout(timer)
    }
  }, [searchParams])

  if (!showNotification) return null

  return (
    <div className="fixed top-4 right-4 z-50 flex items-center gap-2 bg-green-600 text-white px-4 py-2 rounded-lg shadow-lg animate-in fade-in slide-in-from-top-2">
      <CheckCircle className="h-4 w-4" />
      <span>Email confirmed! You can now sign in.</span>
    </div>
  )
}
