"use client"

import { Card, CardContent, CardHeader } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Lightbulb, Bug, User, Clock } from "lucide-react"
import { type Comment } from "@/lib/supabase"

interface CommentCardProps {
  comment: Comment
}

function formatDate(dateString: string): string {
  const date = new Date(dateString)
  const now = new Date()
  const diffMs = now.getTime() - date.getTime()
  const diffDays = Math.floor(diffMs / (1000 * 60 * 60 * 24))

  if (diffDays === 0) {
    const diffHours = Math.floor(diffMs / (1000 * 60 * 60))
    if (diffHours === 0) {
      const diffMins = Math.floor(diffMs / (1000 * 60))
      return diffMins <= 1 ? "just now" : `${diffMins} minutes ago`
    }
    return diffHours === 1 ? "1 hour ago" : `${diffHours} hours ago`
  }
  if (diffDays === 1) return "yesterday"
  if (diffDays < 7) return `${diffDays} days ago`
  if (diffDays < 30) {
    const weeks = Math.floor(diffDays / 7)
    return weeks === 1 ? "1 week ago" : `${weeks} weeks ago`
  }

  return date.toLocaleDateString("en-US", {
    year: "numeric",
    month: "short",
    day: "numeric",
  })
}

export function CommentCard({ comment }: CommentCardProps) {
  const isFeature = comment.category === "feature_request"
  const Icon = isFeature ? Lightbulb : Bug

  return (
    <Card className="transition-shadow hover:shadow-md">
      <CardHeader className="pb-2">
        <div className="flex items-start justify-between gap-4">
          <div className="flex items-center gap-2 min-w-0">
            <Icon className={`h-4 w-4 flex-shrink-0 ${isFeature ? "text-amber-500" : "text-red-500"}`} />
            <h3 className="font-medium text-base leading-tight truncate">{comment.title}</h3>
          </div>
          <Badge variant={isFeature ? "default" : "destructive"} className="flex-shrink-0 text-xs">
            {isFeature ? "Feature" : "Bug"}
          </Badge>
        </div>
        <div className="flex items-center gap-3 text-xs text-muted-foreground mt-1">
          <span className="flex items-center gap-1">
            <User className="h-3 w-3" />
            {comment.name || "Anonymous"}
          </span>
          <span className="flex items-center gap-1">
            <Clock className="h-3 w-3" />
            {formatDate(comment.created_at)}
          </span>
        </div>
      </CardHeader>
      <CardContent className="pt-2">
        <p className="text-sm text-muted-foreground whitespace-pre-wrap">{comment.content}</p>
      </CardContent>
    </Card>
  )
}
