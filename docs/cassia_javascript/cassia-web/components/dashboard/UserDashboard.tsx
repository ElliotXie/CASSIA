'use client'

import { useAuthStore } from '@/lib/stores/auth-store'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { ApiKeyManager } from './ApiKeyManager'

interface UserDashboardProps {
  onNavigate?: (path: string) => void
}

export function UserDashboard({ onNavigate }: UserDashboardProps) {
  const { user, profile } = useAuthStore()

  if (!user) {
    return (
      <div className="p-6 text-center">
        <h2 className="text-2xl font-bold mb-4">Sign in to access your dashboard</h2>
        <p className="text-gray-600">Your API keys and settings will appear here.</p>
      </div>
    )
  }

  return (
    <div className="p-6 max-w-6xl mx-auto space-y-6">
      {/* Header */}
      <div>
        <h1 className="text-3xl font-bold">Welcome back, {profile?.full_name || user.email}</h1>
        <p className="text-gray-600 mt-1">Manage your API keys and settings</p>
      </div>

      {/* Stats Cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium">Member Since</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-sm font-medium">
              {new Date(user.created_at).toLocaleDateString()}
            </div>
          </CardContent>
        </Card>
      </div>

      {/* API Key Management */}
      <ApiKeyManager />
    </div>
  )
}
