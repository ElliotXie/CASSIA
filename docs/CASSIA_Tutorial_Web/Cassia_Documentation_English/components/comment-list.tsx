"use client"

import { useEffect, useState } from "react"
import { supabase, type Comment, type CommentCategory } from "@/lib/supabase"
import { CommentCard } from "./comment-card"
import { Loader2, MessageSquare } from "lucide-react"

interface CommentListProps {
  category: CommentCategory
  refreshTrigger?: number
}

export function CommentList({ category, refreshTrigger }: CommentListProps) {
  const [comments, setComments] = useState<Comment[]>([])
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState<string | null>(null)

  const categoryLabel = category === "feature_request" ? "feature requests" : "bug reports"

  useEffect(() => {
    async function fetchComments() {
      setLoading(true)
      setError(null)

      try {
        const { data, error: fetchError } = await supabase
          .from("comments")
          .select("*")
          .eq("category", category)
          .eq("status", "approved")
          .order("created_at", { ascending: false })

        if (fetchError) {
          throw new Error(fetchError.message)
        }

        setComments(data || [])
      } catch (err) {
        setError(err instanceof Error ? err.message : "Failed to load comments")
      } finally {
        setLoading(false)
      }
    }

    fetchComments()
  }, [category, refreshTrigger])

  if (loading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
      </div>
    )
  }

  if (error) {
    return (
      <div className="text-center py-8 text-red-500">
        <p>Error loading comments: {error}</p>
      </div>
    )
  }

  if (comments.length === 0) {
    return (
      <div className="text-center py-12 text-muted-foreground">
        <MessageSquare className="h-12 w-12 mx-auto mb-3 opacity-50" />
        <p className="text-lg font-medium">No {categoryLabel} yet</p>
        <p className="text-sm mt-1">Be the first to submit one!</p>
      </div>
    )
  }

  return (
    <div className="space-y-4">
      <h3 className="text-lg font-semibold">
        Approved {category === "feature_request" ? "Feature Requests" : "Bug Reports"} ({comments.length})
      </h3>
      <div className="space-y-3">
        {comments.map((comment) => (
          <CommentCard key={comment.id} comment={comment} />
        ))}
      </div>
    </div>
  )
}
