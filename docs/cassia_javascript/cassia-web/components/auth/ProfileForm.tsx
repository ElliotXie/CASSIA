'use client'

import { useState, useEffect } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'

export function ProfileForm() {
  const { user, profile, updateProfile, isLoading, error, clearError } = useAuthStore()
  const [fullName, setFullName] = useState('')
  const [email, setEmail] = useState('')
  const [isEditing, setIsEditing] = useState(false)
  
  useEffect(() => {
    if (profile) {
      setFullName(profile.full_name || '')
      setEmail(profile.email || '')
    } else if (user) {
      setEmail(user.email || '')
    }
  }, [profile, user])
  
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    clearError()
    
    try {
      await updateProfile({
        full_name: fullName,
        email: email
      })
      setIsEditing(false)
    } catch (error) {
      console.error('Profile update error:', error)
    }
  }
  
  const handleEdit = () => {
    setIsEditing(true)
    clearError()
  }
  
  const handleCancel = () => {
    setIsEditing(false)
    clearError()
    // Reset to original values
    if (profile) {
      setFullName(profile.full_name || '')
      setEmail(profile.email || '')
    }
  }
  
  if (!user) {
    return null
  }
  
  return (
    <Card className="w-full max-w-md">
      <CardHeader>
        <CardTitle>Profile</CardTitle>
        <CardDescription>
          Manage your account information
        </CardDescription>
      </CardHeader>
      <CardContent>
        <form onSubmit={handleSubmit} className="space-y-4">
          <div className="space-y-2">
            <Label htmlFor="email">Email</Label>
            <Input
              id="email"
              type="email"
              value={email}
              onChange={(e) => setEmail(e.target.value)}
              disabled={!isEditing}
              required
            />
          </div>
          
          <div className="space-y-2">
            <Label htmlFor="fullName">Full Name</Label>
            <Input
              id="fullName"
              type="text"
              value={fullName}
              onChange={(e) => setFullName(e.target.value)}
              disabled={!isEditing}
              placeholder="Enter your full name"
            />
          </div>
          
          <div className="space-y-2">
            <Label>User ID</Label>
            <Input
              value={user.id}
              disabled
              className="text-gray-500 text-sm"
            />
          </div>
          
          <div className="space-y-2">
            <Label>Member Since</Label>
            <Input
              value={new Date(user.created_at).toLocaleDateString()}
              disabled
              className="text-gray-500 text-sm"
            />
          </div>
          
          {error && (
            <div className="text-red-500 text-sm p-2 bg-red-50 rounded">
              {error}
            </div>
          )}
          
          <div className="flex gap-2">
            {isEditing ? (
              <>
                <Button 
                  type="submit" 
                  className="flex-1"
                  disabled={isLoading}
                >
                  {isLoading ? 'Saving...' : 'Save Changes'}
                </Button>
                <Button 
                  type="button" 
                  variant="outline" 
                  onClick={handleCancel}
                  disabled={isLoading}
                >
                  Cancel
                </Button>
              </>
            ) : (
              <Button 
                type="button" 
                onClick={handleEdit}
                className="w-full"
              >
                Edit Profile
              </Button>
            )}
          </div>
        </form>
      </CardContent>
    </Card>
  )
}