"use client"

import { useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Textarea } from "@/components/ui/textarea"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { AlertCircle, CheckCircle2 } from "lucide-react"
import { supabase, type CommentCategory, type NewComment } from "@/lib/supabase"

interface CommentFormProps {
  category: CommentCategory
  onSuccess?: () => void
}

export function CommentForm({ category, onSuccess }: CommentFormProps) {
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [isSubmitted, setIsSubmitted] = useState(false)
  const [error, setError] = useState("")
  const [formData, setFormData] = useState({
    name: "",
    email: "",
    title: "",
    content: "",
  })

  const categoryLabel = category === "feature_request" ? "Feature Request" : "Bug Report"

  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { id, value } = e.target
    setFormData((prev) => ({
      ...prev,
      [id]: value,
    }))
  }

  const handleSubmit = async (event: React.FormEvent) => {
    event.preventDefault()
    setIsSubmitting(true)
    setError("")

    try {
      const newComment: NewComment = {
        category,
        name: formData.name || undefined,
        email: formData.email || undefined,
        title: formData.title,
        content: formData.content,
      }

      const { error: insertError } = await supabase
        .from("comments")
        .insert([newComment])

      if (insertError) {
        throw new Error(insertError.message)
      }

      setIsSubmitted(true)
      setFormData({ name: "", email: "", title: "", content: "" })
      onSuccess?.()
    } catch (err) {
      setError(err instanceof Error ? err.message : "Something went wrong. Please try again.")
    } finally {
      setIsSubmitting(false)
    }
  }

  const resetForm = () => {
    setIsSubmitted(false)
    setError("")
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle>Submit {categoryLabel}</CardTitle>
        <CardDescription>
          {category === "feature_request"
            ? "Share your ideas for new features or improvements to CASSIA."
            : "Report any bugs or issues you've encountered while using CASSIA."}
        </CardDescription>
      </CardHeader>
      <CardContent>
        {isSubmitted ? (
          <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-md text-green-800 dark:text-green-300">
            <div className="flex items-center gap-2 mb-2">
              <CheckCircle2 className="h-5 w-5" />
              <p className="font-medium">Thank you for your {categoryLabel.toLowerCase()}!</p>
            </div>
            <p className="text-sm mb-3">
              Your submission has been received and will be reviewed by our team. Once approved, it will appear publicly.
            </p>
            <Button variant="outline" size="sm" onClick={resetForm}>
              Submit another
            </Button>
          </div>
        ) : (
          <form className="space-y-4" onSubmit={handleSubmit}>
            {error && (
              <div className="bg-red-50 dark:bg-red-900/20 p-4 rounded-md text-red-800 dark:text-red-300 flex items-start gap-3">
                <AlertCircle className="h-5 w-5 mt-0.5 flex-shrink-0" />
                <p>{error}</p>
              </div>
            )}
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
              <div className="space-y-2">
                <Label htmlFor="name">Name (optional)</Label>
                <Input
                  id="name"
                  placeholder="Your name"
                  value={formData.name}
                  onChange={handleChange}
                />
              </div>
              <div className="space-y-2">
                <Label htmlFor="email">Email (optional)</Label>
                <Input
                  id="email"
                  type="email"
                  placeholder="your@email.com"
                  value={formData.email}
                  onChange={handleChange}
                />
              </div>
            </div>
            <div className="space-y-2">
              <Label htmlFor="title">Title *</Label>
              <Input
                id="title"
                placeholder={category === "feature_request" ? "Brief description of your feature idea" : "Brief description of the bug"}
                required
                value={formData.title}
                onChange={handleChange}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="content">Description *</Label>
              <Textarea
                id="content"
                placeholder={
                  category === "feature_request"
                    ? "Describe your feature request in detail. What problem would it solve? How would it work?"
                    : "Please describe the bug in detail. Include steps to reproduce, expected behavior, and actual behavior."
                }
                className="min-h-[120px]"
                required
                value={formData.content}
                onChange={handleChange}
              />
            </div>
            <Button type="submit" disabled={isSubmitting}>
              {isSubmitting ? "Submitting..." : `Submit ${categoryLabel}`}
            </Button>
          </form>
        )}
      </CardContent>
    </Card>
  )
}
