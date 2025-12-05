'use client'

import { useState, useEffect } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogHeader,
  DialogTitle,
} from '@/components/ui/dialog'
import { AlertCircle, X, CheckCircle } from 'lucide-react'

interface LoginFormProps {
  onToggleMode?: () => void
  onSuccess?: () => void
}

export function LoginForm({ onToggleMode, onSuccess }: LoginFormProps) {
  const [email, setEmail] = useState('')
  const [password, setPassword] = useState('')
  const [isSignUp, setIsSignUp] = useState(false)
  const [fullName, setFullName] = useState('')
  const [showErrorDialog, setShowErrorDialog] = useState(false)
  const [showSuccessDialog, setShowSuccessDialog] = useState(false)

  const { signIn, signUp, signInWithGoogle, isLoading, error, successMessage, clearError, clearSuccess, forceResetLoading } = useAuthStore()
  
  // Force reset loading state when component mounts
  useEffect(() => {
    console.log('LoginForm: Component mounted, forcing reset loading state')
    forceResetLoading()
  }, [forceResetLoading])
  
  // Show error dialog when error occurs
  useEffect(() => {
    if (error) {
      setShowErrorDialog(true)
    }
  }, [error])

  // Show success dialog when success message occurs
  useEffect(() => {
    if (successMessage) {
      setShowSuccessDialog(true)
      // Clear form fields after successful registration
      setEmail('')
      setPassword('')
      setFullName('')
    }
  }, [successMessage])
  
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    clearError()

    try {
      if (isSignUp) {
        await signUp(email, password, fullName)
        // Don't call onSuccess for sign-up - let user see success dialog first
        // onSuccess will be called when they dismiss the success dialog
      } else {
        const success = await signIn(email, password)
        // Only close dialog if sign-in actually succeeded
        if (success) {
          onSuccess?.()
        }
      }
    } catch (error) {
      // Error is handled by the store
      console.error('Authentication error:', error)
    }
  }
  
  const toggleMode = () => {
    setIsSignUp(!isSignUp)
    clearError()
    clearSuccess()
    setShowErrorDialog(false)
    setShowSuccessDialog(false)
    onToggleMode?.()
  }

  const handleCloseErrorDialog = () => {
    setShowErrorDialog(false)
    clearError()
  }

  const handleCloseSuccessDialog = () => {
    setShowSuccessDialog(false)
    clearSuccess()
    // If registration required email confirmation, switch to sign in mode
    if (successMessage?.includes('check your email')) {
      setIsSignUp(false)
    }
    // Close parent dialog after user dismisses success message
    onSuccess?.()
  }
  
  return (
    <Card className="w-full max-w-md">
      <CardHeader>
        <CardTitle>{isSignUp ? 'Create Account' : 'Sign In'}</CardTitle>
        <CardDescription>
          {isSignUp 
            ? 'Create a new account to save your analysis results and API keys securely'
            : 'Sign in to access your saved analysis results and API keys'
          }
        </CardDescription>
      </CardHeader>
      <CardContent>
        <form onSubmit={handleSubmit} className="space-y-4">
          {isSignUp && (
            <div className="space-y-2">
              <Label htmlFor="fullName">Full Name</Label>
              <Input
                id="fullName"
                type="text"
                placeholder="Enter your full name"
                value={fullName}
                onChange={(e) => setFullName(e.target.value)}
                required
              />
            </div>
          )}
          
          <div className="space-y-2">
            <Label htmlFor="email">Email</Label>
            <Input
              id="email"
              type="email"
              placeholder="Enter your email"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
              required
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="password">Password</Label>
            <Input
              id="password"
              type="password"
              placeholder="Enter your password"
              value={password}
              onChange={(e) => setPassword(e.target.value)}
              required
              minLength={6}
            />
          </div>
          
          
          <Button
            type="submit"
            className="w-full"
            disabled={isLoading}
          >
            {isLoading ? 'Please wait...' : (isSignUp ? 'Create Account' : 'Sign In')}
          </Button>

          {/* Divider */}
          <div className="relative my-4">
            <div className="absolute inset-0 flex items-center">
              <span className="w-full border-t" />
            </div>
            <div className="relative flex justify-center text-xs uppercase">
              <span className="bg-background px-2 text-muted-foreground">Or continue with</span>
            </div>
          </div>

          {/* Google Sign-In Button */}
          <Button
            type="button"
            variant="outline"
            onClick={() => signInWithGoogle()}
            disabled={isLoading}
            className="w-full"
          >
            <svg className="mr-2 h-4 w-4" viewBox="0 0 24 24">
              <path
                fill="currentColor"
                d="M22.56 12.25c0-.78-.07-1.53-.2-2.25H12v4.26h5.92c-.26 1.37-1.04 2.53-2.21 3.31v2.77h3.57c2.08-1.92 3.28-4.74 3.28-8.09z"
              />
              <path
                fill="currentColor"
                d="M12 23c2.97 0 5.46-.98 7.28-2.66l-3.57-2.77c-.98.66-2.23 1.06-3.71 1.06-2.86 0-5.29-1.93-6.16-4.53H2.18v2.84C3.99 20.53 7.7 23 12 23z"
              />
              <path
                fill="currentColor"
                d="M5.84 14.09c-.22-.66-.35-1.36-.35-2.09s.13-1.43.35-2.09V7.07H2.18C1.43 8.55 1 10.22 1 12s.43 3.45 1.18 4.93l2.85-2.22.81-.62z"
              />
              <path
                fill="currentColor"
                d="M12 5.38c1.62 0 3.06.56 4.21 1.64l3.15-3.15C17.45 2.09 14.97 1 12 1 7.7 1 3.99 3.47 2.18 7.07l3.66 2.84c.87-2.6 3.3-4.53 6.16-4.53z"
              />
            </svg>
            Continue with Google
          </Button>

          <div className="text-center">
            <Button 
              type="button" 
              variant="link" 
              onClick={toggleMode}
              className="text-sm"
            >
              {isSignUp 
                ? 'Already have an account? Sign in'
                : "Don't have an account? Create one"
              }
            </Button>
          </div>
        </form>
      </CardContent>
      
      {/* Error Dialog */}
      <Dialog open={showErrorDialog} onOpenChange={setShowErrorDialog}>
        <DialogContent className="sm:max-w-md">
          <DialogHeader>
            <DialogTitle className="flex items-center gap-2 text-red-600">
              <AlertCircle className="h-5 w-5" />
              Authentication Error
            </DialogTitle>
            <DialogDescription className="text-left">
              {error}
            </DialogDescription>
          </DialogHeader>
          <div className="flex justify-end gap-2 mt-4">
            <Button
              variant="outline"
              onClick={handleCloseErrorDialog}
              className="flex items-center gap-2"
            >
              <X className="h-4 w-4" />
              Close
            </Button>
          </div>
        </DialogContent>
      </Dialog>

      {/* Success Dialog */}
      <Dialog open={showSuccessDialog} onOpenChange={setShowSuccessDialog}>
        <DialogContent className="sm:max-w-md">
          <DialogHeader>
            <DialogTitle className="flex items-center gap-2 text-green-600">
              <CheckCircle className="h-5 w-5" />
              Success
            </DialogTitle>
            <DialogDescription className="text-left">
              {successMessage}
            </DialogDescription>
          </DialogHeader>
          <div className="flex justify-end gap-2 mt-4">
            <Button
              onClick={handleCloseSuccessDialog}
              className="flex items-center gap-2"
            >
              OK
            </Button>
          </div>
        </DialogContent>
      </Dialog>
    </Card>
  )
}