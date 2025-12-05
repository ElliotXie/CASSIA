'use client'

import { useState, useEffect } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { Button } from '@/components/ui/button'
import { LoginForm } from './LoginForm'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from '@/components/ui/dialog'
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import { User, LogOut, CheckCircle } from 'lucide-react'

export function AuthButton() {
  const { user, profile, signOut, isAuthenticated, successMessage, clearSuccess } = useAuthStore()
  const [showLogin, setShowLogin] = useState(false)
  const [showNotification, setShowNotification] = useState(false)

  // Show notification when signed in or signed out
  useEffect(() => {
    if (successMessage) {
      setShowNotification(true)
      const timer = setTimeout(() => {
        setShowNotification(false)
        clearSuccess()
      }, 3000)
      return () => clearTimeout(timer)
    }
  }, [successMessage, clearSuccess])

  const handleSignOut = async () => {
    await signOut()
  }
  
  if (!isAuthenticated) {
    return (
      <>
        {/* Success notification */}
        {showNotification && successMessage && (
          <div className="fixed top-4 right-4 z-50 flex items-center gap-2 bg-green-600 text-white px-4 py-2 rounded-lg shadow-lg animate-in fade-in slide-in-from-top-2">
            <CheckCircle className="h-4 w-4" />
            <span>{successMessage}</span>
          </div>
        )}
        <Dialog open={showLogin} onOpenChange={setShowLogin}>
          <DialogTrigger asChild>
            <Button variant="outline" size="sm" className="glass border-white/30 hover:bg-white/20 btn-modern">
              <User className="h-4 w-4 mr-2" />
              Sign In
            </Button>
          </DialogTrigger>
          <DialogContent className="sm:max-w-md">
            <DialogHeader>
              <DialogTitle>Authentication</DialogTitle>
            </DialogHeader>
            <LoginForm onSuccess={() => setShowLogin(false)} />
          </DialogContent>
        </Dialog>
      </>
    )
  }
  
  return (
    <>
      {/* Success notification */}
      {showNotification && successMessage && (
        <div className="fixed top-4 right-4 z-50 flex items-center gap-2 bg-green-600 text-white px-4 py-2 rounded-lg shadow-lg animate-in fade-in slide-in-from-top-2">
          <CheckCircle className="h-4 w-4" />
          <span>{successMessage}</span>
        </div>
      )}
      <div className="flex items-center gap-2">
        <DropdownMenu>
        <DropdownMenuTrigger asChild>
          <Button variant="outline" size="sm" className="glass border-white/30 hover:bg-white/20 btn-modern flex items-center gap-2">
            <User className="h-4 w-4" />
            {profile?.full_name || user?.email || 'User'}
          </Button>
        </DropdownMenuTrigger>
        <DropdownMenuContent align="end" className="w-56">
          <DropdownMenuItem onClick={handleSignOut}>
            <LogOut className="h-4 w-4 mr-2" />
            Sign Out
          </DropdownMenuItem>
        </DropdownMenuContent>
      </DropdownMenu>
      </div>
    </>
  )
}