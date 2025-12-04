"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Badge } from "@/components/ui/badge"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import {
  Loader2,
  CheckCircle,
  XCircle,
  Trash2,
  Lightbulb,
  Bug,
  Clock,
  User,
  Mail,
  Lock,
} from "lucide-react"
import { type Comment, type CommentStatus } from "@/lib/supabase"

export default function AdminCommentsPage() {
  const [password, setPassword] = useState("")
  const [isAuthenticated, setIsAuthenticated] = useState(false)
  const [authError, setAuthError] = useState("")
  const [comments, setComments] = useState<Comment[]>([])
  const [loading, setLoading] = useState(false)
  const [activeTab, setActiveTab] = useState<CommentStatus | "all">("pending")
  const [actionLoading, setActionLoading] = useState<string | null>(null)

  // Check if already authenticated (sessionStorage)
  useEffect(() => {
    const savedPassword = sessionStorage.getItem("admin_password")
    if (savedPassword) {
      setPassword(savedPassword)
      setIsAuthenticated(true)
    }
  }, [])

  // Fetch comments when authenticated
  useEffect(() => {
    if (isAuthenticated) {
      fetchComments()
    }
  }, [isAuthenticated, activeTab])

  const handleLogin = async (e: React.FormEvent) => {
    e.preventDefault()
    setAuthError("")
    setLoading(true)

    try {
      const response = await fetch("/api/admin/comments?status=pending", {
        headers: {
          "x-admin-password": password,
        },
      })

      if (response.status === 401) {
        setAuthError("Invalid password")
        return
      }

      if (!response.ok) {
        throw new Error("Failed to authenticate")
      }

      sessionStorage.setItem("admin_password", password)
      setIsAuthenticated(true)
    } catch (error) {
      setAuthError("Authentication failed")
    } finally {
      setLoading(false)
    }
  }

  const fetchComments = async () => {
    setLoading(true)
    try {
      const statusParam = activeTab === "all" ? "" : `?status=${activeTab}`
      const response = await fetch(`/api/admin/comments${statusParam}`, {
        headers: {
          "x-admin-password": password,
        },
      })

      if (!response.ok) {
        throw new Error("Failed to fetch comments")
      }

      const data = await response.json()
      setComments(data.comments || [])
    } catch (error) {
      console.error("Error fetching comments:", error)
    } finally {
      setLoading(false)
    }
  }

  const updateCommentStatus = async (id: string, status: CommentStatus) => {
    setActionLoading(id)
    try {
      const response = await fetch("/api/admin/comments", {
        method: "PATCH",
        headers: {
          "Content-Type": "application/json",
          "x-admin-password": password,
        },
        body: JSON.stringify({ id, status }),
      })

      if (!response.ok) {
        throw new Error("Failed to update comment")
      }

      // Refresh comments list
      fetchComments()
    } catch (error) {
      console.error("Error updating comment:", error)
    } finally {
      setActionLoading(null)
    }
  }

  const deleteComment = async (id: string) => {
    if (!confirm("Are you sure you want to delete this comment?")) {
      return
    }

    setActionLoading(id)
    try {
      const response = await fetch(`/api/admin/comments?id=${id}`, {
        method: "DELETE",
        headers: {
          "x-admin-password": password,
        },
      })

      if (!response.ok) {
        throw new Error("Failed to delete comment")
      }

      // Refresh comments list
      fetchComments()
    } catch (error) {
      console.error("Error deleting comment:", error)
    } finally {
      setActionLoading(null)
    }
  }

  const handleLogout = () => {
    sessionStorage.removeItem("admin_password")
    setIsAuthenticated(false)
    setPassword("")
    setComments([])
  }

  // Login form
  if (!isAuthenticated) {
    return (
      <div className="container mx-auto py-12 px-4">
        <div className="max-w-md mx-auto">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Lock className="h-5 w-5" />
                Admin Login
              </CardTitle>
            </CardHeader>
            <CardContent>
              <form onSubmit={handleLogin} className="space-y-4">
                {authError && (
                  <div className="bg-red-50 dark:bg-red-900/20 p-3 rounded-md text-red-800 dark:text-red-300 text-sm">
                    {authError}
                  </div>
                )}
                <div className="space-y-2">
                  <Label htmlFor="password">Password</Label>
                  <Input
                    id="password"
                    type="password"
                    value={password}
                    onChange={(e) => setPassword(e.target.value)}
                    placeholder="Enter admin password"
                    required
                  />
                </div>
                <Button type="submit" disabled={loading} className="w-full">
                  {loading ? <Loader2 className="h-4 w-4 animate-spin mr-2" /> : null}
                  Login
                </Button>
              </form>
            </CardContent>
          </Card>
        </div>
      </div>
    )
  }

  // Admin dashboard
  const pendingCount = comments.filter((c) => c.status === "pending").length
  const approvedCount = comments.filter((c) => c.status === "approved").length
  const rejectedCount = comments.filter((c) => c.status === "rejected").length

  return (
    <div className="container mx-auto py-12 px-4">
      <div className="max-w-4xl mx-auto">
        <div className="flex items-center justify-between mb-6">
          <h1 className="text-3xl font-bold">Comment Moderation</h1>
          <Button variant="outline" size="sm" onClick={handleLogout}>
            Logout
          </Button>
        </div>

        <Tabs value={activeTab} onValueChange={(v) => setActiveTab(v as CommentStatus | "all")}>
          <TabsList className="mb-6">
            <TabsTrigger value="pending" className="flex items-center gap-2">
              <Clock className="h-4 w-4" />
              Pending
              {activeTab !== "pending" && pendingCount > 0 && (
                <Badge variant="secondary" className="ml-1">{pendingCount}</Badge>
              )}
            </TabsTrigger>
            <TabsTrigger value="approved">
              <CheckCircle className="h-4 w-4 mr-2" />
              Approved
            </TabsTrigger>
            <TabsTrigger value="rejected">
              <XCircle className="h-4 w-4 mr-2" />
              Rejected
            </TabsTrigger>
            <TabsTrigger value="all">All</TabsTrigger>
          </TabsList>

          {loading ? (
            <div className="flex items-center justify-center py-12">
              <Loader2 className="h-6 w-6 animate-spin" />
            </div>
          ) : comments.length === 0 ? (
            <div className="text-center py-12 text-muted-foreground">
              <p>No {activeTab === "all" ? "" : activeTab} comments found.</p>
            </div>
          ) : (
            <div className="space-y-4">
              {comments.map((comment) => (
                <AdminCommentCard
                  key={comment.id}
                  comment={comment}
                  onApprove={() => updateCommentStatus(comment.id, "approved")}
                  onReject={() => updateCommentStatus(comment.id, "rejected")}
                  onDelete={() => deleteComment(comment.id)}
                  isLoading={actionLoading === comment.id}
                />
              ))}
            </div>
          )}
        </Tabs>
      </div>
    </div>
  )
}

interface AdminCommentCardProps {
  comment: Comment
  onApprove: () => void
  onReject: () => void
  onDelete: () => void
  isLoading: boolean
}

function AdminCommentCard({ comment, onApprove, onReject, onDelete, isLoading }: AdminCommentCardProps) {
  const isFeature = comment.category === "feature_request"
  const Icon = isFeature ? Lightbulb : Bug

  const statusColors = {
    pending: "bg-yellow-100 text-yellow-800 dark:bg-yellow-900/30 dark:text-yellow-300",
    approved: "bg-green-100 text-green-800 dark:bg-green-900/30 dark:text-green-300",
    rejected: "bg-red-100 text-red-800 dark:bg-red-900/30 dark:text-red-300",
  }

  return (
    <Card>
      <CardHeader className="pb-2">
        <div className="flex items-start justify-between gap-4">
          <div className="flex items-center gap-2 min-w-0">
            <Icon className={`h-4 w-4 flex-shrink-0 ${isFeature ? "text-amber-500" : "text-red-500"}`} />
            <h3 className="font-medium">{comment.title}</h3>
          </div>
          <div className="flex items-center gap-2 flex-shrink-0">
            <Badge variant={isFeature ? "default" : "destructive"} className="text-xs">
              {isFeature ? "Feature" : "Bug"}
            </Badge>
            <span className={`px-2 py-0.5 rounded text-xs font-medium ${statusColors[comment.status]}`}>
              {comment.status}
            </span>
          </div>
        </div>
        <div className="flex flex-wrap items-center gap-3 text-xs text-muted-foreground mt-1">
          <span className="flex items-center gap-1">
            <User className="h-3 w-3" />
            {comment.name || "Anonymous"}
          </span>
          {comment.email && (
            <span className="flex items-center gap-1">
              <Mail className="h-3 w-3" />
              {comment.email}
            </span>
          )}
          <span className="flex items-center gap-1">
            <Clock className="h-3 w-3" />
            {new Date(comment.created_at).toLocaleString()}
          </span>
        </div>
      </CardHeader>
      <CardContent className="pt-2">
        <p className="text-sm text-muted-foreground whitespace-pre-wrap mb-4">{comment.content}</p>
        <div className="flex items-center gap-2">
          {comment.status !== "approved" && (
            <Button
              size="sm"
              variant="default"
              onClick={onApprove}
              disabled={isLoading}
              className="bg-green-600 hover:bg-green-700"
            >
              {isLoading ? <Loader2 className="h-4 w-4 animate-spin" /> : <CheckCircle className="h-4 w-4 mr-1" />}
              Approve
            </Button>
          )}
          {comment.status !== "rejected" && (
            <Button size="sm" variant="outline" onClick={onReject} disabled={isLoading}>
              {isLoading ? <Loader2 className="h-4 w-4 animate-spin" /> : <XCircle className="h-4 w-4 mr-1" />}
              Reject
            </Button>
          )}
          <Button size="sm" variant="ghost" onClick={onDelete} disabled={isLoading} className="text-red-500 hover:text-red-600 hover:bg-red-50">
            <Trash2 className="h-4 w-4" />
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}
