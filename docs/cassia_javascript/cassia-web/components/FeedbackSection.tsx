'use client'

import React, { useState, useEffect } from 'react'
import { MessageSquare, Send, CheckCircle, AlertCircle, Bug, Lightbulb, User, Mail, Loader2, Clock } from 'lucide-react'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Label } from '@/components/ui/label'
import { Badge } from '@/components/ui/badge'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { useAuthStore } from '@/lib/stores/auth-store'
import { useCommentsStore } from '@/lib/stores/comments-store'
import emailjs from '@emailjs/browser'
import { CommentType, FeedbackComment } from '@/lib/supabase/types'

interface FeedbackFormData {
  comment_type: CommentType | ''
  name: string
  email: string
  message: string
}

export function FeedbackSection() {
  const { isAuthenticated } = useAuthStore()
  const {
    comments,
    isLoading,
    isSubmitting,
    error,
    successMessage,
    loadApprovedComments,
    submitComment,
    clearError,
    clearSuccess
  } = useCommentsStore()

  const [showForm, setShowForm] = useState(false)
  const [formData, setFormData] = useState<FeedbackFormData>({
    comment_type: '',
    name: '',
    email: '',
    message: ''
  })
  const [emailSending, setEmailSending] = useState(false)

  // Load approved comments on mount
  useEffect(() => {
    loadApprovedComments()
  }, [loadApprovedComments])

  // Clear messages after timeout
  useEffect(() => {
    if (successMessage) {
      const timer = setTimeout(() => {
        clearSuccess()
        setShowForm(false)
        setFormData({ comment_type: '', name: '', email: '', message: '' })
      }, 3000)
      return () => clearTimeout(timer)
    }
  }, [successMessage, clearSuccess])

  const handleInputChange = (
    e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>
  ) => {
    const { name, value } = e.target
    setFormData(prev => ({ ...prev, [name]: value }))
  }

  const handleTypeChange = (value: CommentType) => {
    setFormData(prev => ({ ...prev, comment_type: value }))
  }

  const sendAdminNotification = async (comment: FeedbackComment) => {
    try {
      setEmailSending(true)
      const approvalUrl = `${window.location.origin}/api/comments/approve?token=${comment.approval_token}`

      await emailjs.send(
        'service_m58e32j',
        'template_po0safd',
        {
          name: 'CASSIA Feedback System',
          email: 'noreply@cassia.bio',
          title: `New ${comment.comment_type === 'feature_request' ? 'Feature Request' : 'Bug Report'} Submitted`,
          message: `
New feedback submission requires approval:

Type: ${comment.comment_type === 'feature_request' ? 'Feature Request' : 'Bug Report'}
From: ${comment.name || 'Anonymous'} (${comment.email || 'No email provided'})

Message:
${comment.message}

---
Click here to approve this comment:
${approvalUrl}
          `
        },
        'qS1BL05ocjba4E3nE'
      )
    } catch (err) {
      console.error('Failed to send admin notification:', err)
    } finally {
      setEmailSending(false)
    }
  }

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    clearError()

    if (!formData.comment_type) {
      return
    }

    const result = await submitComment({
      comment_type: formData.comment_type as CommentType,
      name: formData.name || undefined,
      email: formData.email || undefined,
      message: formData.message
    })

    if (result.success && result.comment) {
      // Send notification email to admin
      await sendAdminNotification(result.comment)
    }
  }

  const formatDate = (dateString: string) => {
    return new Date(dateString).toLocaleDateString('en-US', {
      year: 'numeric',
      month: 'short',
      day: 'numeric'
    })
  }

  return (
    <div className="space-y-4">
      {/* Feedback Toggle Button */}
      <Button
        variant="outline"
        className="w-full justify-start"
        onClick={() => setShowForm(!showForm)}
      >
        <span className="w-6 text-center flex-shrink-0">
          <MessageSquare className="h-4 w-4 inline" />
        </span>
        <span className="ml-2">Submit Feedback</span>
      </Button>

      {/* Feedback Form */}
      {showForm && (
        <Card className="border-2 border-primary/20">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center space-x-2 text-base">
              <MessageSquare className="h-4 w-4" />
              <span>Share Your Feedback</span>
            </CardTitle>
            <CardDescription className="text-xs">
              Submit feature requests or report bugs
            </CardDescription>
          </CardHeader>
          <CardContent>
            {!isAuthenticated ? (
              <div className="p-3 rounded-lg bg-amber-50 border border-amber-200 text-amber-800">
                <div className="flex items-center space-x-2 mb-1">
                  <AlertCircle className="h-4 w-4" />
                  <span className="font-medium text-sm">Login Required</span>
                </div>
                <p className="text-xs">
                  Please sign in to submit feedback.
                </p>
              </div>
            ) : (
              <form onSubmit={handleSubmit} className="space-y-3">
                {/* Comment Type */}
                <div className="space-y-1">
                  <Label htmlFor="comment_type" className="text-xs">Feedback Type *</Label>
                  <Select
                    value={formData.comment_type}
                    onValueChange={handleTypeChange}
                  >
                    <SelectTrigger className="h-9">
                      <SelectValue placeholder="Select type" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="feature_request">
                        <div className="flex items-center space-x-2">
                          <Lightbulb className="h-3 w-3 text-blue-500" />
                          <span>Feature Request</span>
                        </div>
                      </SelectItem>
                      <SelectItem value="bug_report">
                        <div className="flex items-center space-x-2">
                          <Bug className="h-3 w-3 text-red-500" />
                          <span>Bug Report</span>
                        </div>
                      </SelectItem>
                    </SelectContent>
                  </Select>
                </div>

                {/* Optional Name & Email */}
                <div className="grid grid-cols-2 gap-2">
                  <div className="space-y-1">
                    <Label htmlFor="name" className="text-xs flex items-center space-x-1">
                      <User className="h-3 w-3" />
                      <span>Name</span>
                    </Label>
                    <Input
                      id="name"
                      name="name"
                      value={formData.name}
                      onChange={handleInputChange}
                      placeholder="Optional"
                      className="h-8 text-sm"
                    />
                  </div>
                  <div className="space-y-1">
                    <Label htmlFor="email" className="text-xs flex items-center space-x-1">
                      <Mail className="h-3 w-3" />
                      <span>Email</span>
                    </Label>
                    <Input
                      id="email"
                      name="email"
                      type="email"
                      value={formData.email}
                      onChange={handleInputChange}
                      placeholder="Optional"
                      className="h-8 text-sm"
                    />
                  </div>
                </div>

                {/* Message */}
                <div className="space-y-1">
                  <Label htmlFor="message" className="text-xs">Your Feedback *</Label>
                  <Textarea
                    id="message"
                    name="message"
                    value={formData.message}
                    onChange={handleInputChange}
                    required
                    placeholder="Describe your feedback..."
                    rows={3}
                    className="text-sm"
                  />
                </div>

                {/* Error/Success Messages */}
                {error && (
                  <div className="flex items-center text-red-600 text-xs space-x-1">
                    <AlertCircle className="h-3 w-3" />
                    <span>{error}</span>
                  </div>
                )}
                {successMessage && (
                  <div className="flex items-center text-green-600 text-xs space-x-1">
                    <CheckCircle className="h-3 w-3" />
                    <span>{successMessage}</span>
                  </div>
                )}

                {/* Submit Button */}
                <Button
                  type="submit"
                  disabled={isSubmitting || emailSending || !formData.comment_type || !formData.message}
                  className="w-full h-9"
                  size="sm"
                >
                  {isSubmitting || emailSending ? (
                    <>
                      <Loader2 className="h-3 w-3 mr-2 animate-spin" />
                      <span>Submitting...</span>
                    </>
                  ) : (
                    <>
                      <Send className="h-3 w-3 mr-2" />
                      <span>Submit Feedback</span>
                    </>
                  )}
                </Button>
              </form>
            )}
          </CardContent>
        </Card>
      )}

      {/* Approved Comments List */}
      {comments.length > 0 && (
        <Card>
          <CardHeader className="pb-2">
            <CardTitle className="text-base">Community Feedback</CardTitle>
          </CardHeader>
          <CardContent className="space-y-3">
            {isLoading ? (
              <div className="flex items-center justify-center py-4">
                <Loader2 className="h-5 w-5 animate-spin text-muted-foreground" />
              </div>
            ) : (
              comments.slice(0, 5).map((comment) => (
                <div
                  key={comment.id}
                  className="p-3 border rounded-lg bg-muted/30 space-y-2"
                >
                  <div className="flex items-center justify-between">
                    <Badge
                      variant={comment.comment_type === 'feature_request' ? 'default' : 'destructive'}
                      className="text-xs h-5"
                    >
                      {comment.comment_type === 'feature_request' ? (
                        <><Lightbulb className="h-3 w-3 mr-1" />Feature</>
                      ) : (
                        <><Bug className="h-3 w-3 mr-1" />Bug</>
                      )}
                    </Badge>
                    <span className="text-xs text-muted-foreground flex items-center">
                      <Clock className="h-3 w-3 mr-1" />
                      {formatDate(comment.created_at)}
                    </span>
                  </div>
                  <p className="text-sm">{comment.message}</p>
                  {comment.name && (
                    <p className="text-xs text-muted-foreground">
                      - {comment.name}
                    </p>
                  )}
                </div>
              ))
            )}
          </CardContent>
        </Card>
      )}
    </div>
  )
}
