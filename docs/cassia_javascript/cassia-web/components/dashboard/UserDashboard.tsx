'use client'

import { useState, useEffect } from 'react'
import { useAuthStore } from '@/lib/stores/auth-store'
import { useResultsStore } from '@/lib/stores/results-store'
import { useApiKeyStore } from '@/lib/stores/api-key-store-simple'
import { Button } from '@/components/ui/button'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Badge } from '@/components/ui/badge'
import { 
  Database, 
  Download, 
  FileText, 
  Settings, 
  Trash2, 
  Upload,
  Calendar,
  Activity
} from 'lucide-react'
import { AnalysisResult, AnalysisSession } from '@/lib/supabase/types'

interface UserDashboardProps {
  onNavigate?: (path: string) => void
}

export function UserDashboard({ onNavigate }: UserDashboardProps) {
  const { user, profile } = useAuthStore()
  const { 
    results, 
    sessions, 
    loadResults, 
    loadSessions, 
    deleteResult, 
    deleteSession,
    exportResults,
    isLoading,
    error 
  } = useResultsStore()
  const { apiKeys } = useApiKeyStore()
  
  const [selectedResults, setSelectedResults] = useState<string[]>([])
  
  useEffect(() => {
    if (user) {
      loadResults(20) // Load last 20 results
      loadSessions()
    }
  }, [user, loadResults, loadSessions])
  
  const handleExport = async () => {
    if (selectedResults.length === 0) return
    
    try {
      const blob = await exportResults(selectedResults)
      const url = URL.createObjectURL(blob)
      const a = document.createElement('a')
      a.href = url
      a.download = `cassia-results-${new Date().toISOString().split('T')[0]}.json`
      document.body.appendChild(a)
      a.click()
      document.body.removeChild(a)
      URL.revokeObjectURL(url)
    } catch (error) {
      console.error('Export error:', error)
    }
  }
  
  const handleDeleteResult = async (id: string) => {
    if (confirm('Are you sure you want to delete this result?')) {
      await deleteResult(id)
    }
  }
  
  const handleDeleteSession = async (id: string) => {
    if (confirm('Are you sure you want to delete this session?')) {
      await deleteSession(id)
    }
  }
  
  const getAnalysisTypeColor = (type: string) => {
    switch (type) {
      case 'batch': return 'bg-blue-100 text-blue-800'
      case 'symphony': return 'bg-purple-100 text-purple-800'
      case 'scoring': return 'bg-green-100 text-green-800'
      case 'subclustering': return 'bg-orange-100 text-orange-800'
      case 'annotation-boost': return 'bg-pink-100 text-pink-800'
      default: return 'bg-gray-100 text-gray-800'
    }
  }
  
  const getStatusColor = (status: string) => {
    switch (status) {
      case 'active': return 'bg-green-100 text-green-800'
      case 'completed': return 'bg-blue-100 text-blue-800'
      case 'archived': return 'bg-gray-100 text-gray-800'
      default: return 'bg-gray-100 text-gray-800'
    }
  }
  
  const hasApiKeys = Object.values(apiKeys).some(key => key.length > 0)
  
  if (!user) {
    return (
      <div className="p-6 text-center">
        <h2 className="text-2xl font-bold mb-4">Sign in to access your dashboard</h2>
        <p className="text-gray-600">Your analysis history and saved results will appear here.</p>
      </div>
    )
  }
  
  return (
    <div className="p-6 max-w-6xl mx-auto space-y-6">
      {/* Header */}
      <div className="flex justify-between items-center">
        <div>
          <h1 className="text-3xl font-bold">Welcome back, {profile?.full_name || user.email}</h1>
          <p className="text-gray-600 mt-1">Manage your analysis results and settings</p>
        </div>
        <Button onClick={() => onNavigate?.('/settings')}>
          <Settings className="h-4 w-4 mr-2" />
          Settings
        </Button>
      </div>
      
      {/* Stats Cards */}
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium">Total Results</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{results.length}</div>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium">Active Sessions</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">
              {sessions.filter(s => s.status === 'active').length}
            </div>
          </CardContent>
        </Card>
        
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-sm font-medium">API Keys</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">
              {Object.values(apiKeys).filter(key => key.length > 0).length}
            </div>
          </CardContent>
        </Card>
        
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
      
      {/* Quick Actions */}
      <Card>
        <CardHeader>
          <CardTitle>Quick Actions</CardTitle>
          <CardDescription>Get started with your analysis</CardDescription>
        </CardHeader>
        <CardContent>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <Button 
              variant="outline" 
              onClick={() => onNavigate?.('/batch')}
              className="h-20 flex-col"
            >
              <Database className="h-6 w-6 mb-2" />
              Batch Analysis
            </Button>
            <Button 
              variant="outline" 
              onClick={() => onNavigate?.('/agents/symphony-compare')}
              className="h-20 flex-col"
            >
              <Activity className="h-6 w-6 mb-2" />
              Symphony Compare
            </Button>
            <Button 
              variant="outline" 
              onClick={() => onNavigate?.('/agents/scoring')}
              className="h-20 flex-col"
            >
              <FileText className="h-6 w-6 mb-2" />
              Scoring
            </Button>
            <Button 
              variant="outline" 
              onClick={() => onNavigate?.('/subclustering')}
              className="h-20 flex-col"
            >
              <Settings className="h-6 w-6 mb-2" />
              Subclustering
            </Button>
          </div>
        </CardContent>
      </Card>
      
      {/* Recent Results */}
      <Card>
        <CardHeader>
          <div className="flex justify-between items-center">
            <div>
              <CardTitle>Recent Analysis Results</CardTitle>
              <CardDescription>Your latest analysis results</CardDescription>
            </div>
            <div className="flex gap-2">
              <Button 
                variant="outline" 
                size="sm"
                onClick={handleExport}
                disabled={selectedResults.length === 0}
              >
                <Download className="h-4 w-4 mr-2" />
                Export ({selectedResults.length})
              </Button>
            </div>
          </div>
        </CardHeader>
        <CardContent>
          {isLoading ? (
            <div className="text-center py-8">Loading results...</div>
          ) : results.length === 0 ? (
            <div className="text-center py-8">
              <FileText className="h-12 w-12 text-gray-400 mx-auto mb-4" />
              <p className="text-gray-600">No analysis results yet</p>
              <p className="text-sm text-gray-500 mt-2">
                Start an analysis to see your results here
              </p>
            </div>
          ) : (
            <div className="space-y-4">
              {results.slice(0, 10).map((result) => (
                <div key={result.id} className="flex items-center space-x-4 p-4 border rounded-lg">
                  <input
                    type="checkbox"
                    checked={selectedResults.includes(result.id)}
                    onChange={(e) => {
                      if (e.target.checked) {
                        setSelectedResults([...selectedResults, result.id])
                      } else {
                        setSelectedResults(selectedResults.filter(id => id !== result.id))
                      }
                    }}
                    className="rounded"
                  />
                  <div className="flex-1">
                    <div className="flex items-center gap-2">
                      <h3 className="font-semibold">{result.title || 'Untitled Analysis'}</h3>
                      <Badge className={getAnalysisTypeColor(result.analysis_type)}>
                        {result.analysis_type}
                      </Badge>
                    </div>
                    <p className="text-sm text-gray-600 mt-1">{result.description}</p>
                    <div className="flex items-center gap-4 mt-2 text-xs text-gray-500">
                      <span className="flex items-center gap-1">
                        <Calendar className="h-3 w-3" />
                        {new Date(result.created_at).toLocaleDateString()}
                      </span>
                    </div>
                  </div>
                  <div className="flex items-center gap-2">
                    <Button 
                      variant="ghost" 
                      size="sm"
                      onClick={() => handleDeleteResult(result.id)}
                    >
                      <Trash2 className="h-4 w-4 text-red-500" />
                    </Button>
                  </div>
                </div>
              ))}
            </div>
          )}
        </CardContent>
      </Card>
      
      {/* Active Sessions */}
      <Card>
        <CardHeader>
          <CardTitle>Active Sessions</CardTitle>
          <CardDescription>Resume your ongoing analysis work</CardDescription>
        </CardHeader>
        <CardContent>
          {sessions.filter(s => s.status === 'active').length === 0 ? (
            <div className="text-center py-8">
              <Activity className="h-12 w-12 text-gray-400 mx-auto mb-4" />
              <p className="text-gray-600">No active sessions</p>
            </div>
          ) : (
            <div className="space-y-4">
              {sessions
                .filter(s => s.status === 'active')
                .slice(0, 5)
                .map((session) => (
                  <div key={session.id} className="flex items-center space-x-4 p-4 border rounded-lg">
                    <div className="flex-1">
                      <div className="flex items-center gap-2">
                        <h3 className="font-semibold">{session.session_name || 'Untitled Session'}</h3>
                        <Badge className={getAnalysisTypeColor(session.analysis_type || 'unknown')}>
                          {session.analysis_type}
                        </Badge>
                        <Badge className={getStatusColor(session.status)}>
                          {session.status}
                        </Badge>
                      </div>
                      <div className="flex items-center gap-4 mt-2 text-xs text-gray-500">
                        <span className="flex items-center gap-1">
                          <Calendar className="h-3 w-3" />
                          {new Date(session.created_at).toLocaleDateString()}
                        </span>
                      </div>
                    </div>
                    <div className="flex items-center gap-2">
                      <Button 
                        variant="outline" 
                        size="sm"
                        onClick={() => onNavigate?.(`/${session.analysis_type}?session=${session.id}`)}
                      >
                        Resume
                      </Button>
                      <Button 
                        variant="ghost" 
                        size="sm"
                        onClick={() => handleDeleteSession(session.id)}
                      >
                        <Trash2 className="h-4 w-4 text-red-500" />
                      </Button>
                    </div>
                  </div>
                ))}
            </div>
          )}
        </CardContent>
      </Card>
      
      {/* Setup Reminder */}
      {!hasApiKeys && (
        <Card className="border-yellow-200 bg-yellow-50">
          <CardHeader>
            <CardTitle className="text-yellow-800">Setup Required</CardTitle>
            <CardDescription className="text-yellow-700">
              You haven't configured any API keys yet. Set up your API keys to start running analyses.
            </CardDescription>
          </CardHeader>
          <CardContent>
            <Button onClick={() => onNavigate?.('/settings')}>
              <Settings className="h-4 w-4 mr-2" />
              Configure API Keys
            </Button>
          </CardContent>
        </Card>
      )}
      
      {error && (
        <Card className="border-red-200 bg-red-50">
          <CardContent className="pt-6">
            <p className="text-red-800">{error}</p>
          </CardContent>
        </Card>
      )}
    </div>
  )
}