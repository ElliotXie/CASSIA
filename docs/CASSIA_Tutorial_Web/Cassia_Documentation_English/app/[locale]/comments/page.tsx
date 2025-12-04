"use client"

import { useState } from "react"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { CommentForm } from "@/components/comment-form"
import { CommentList } from "@/components/comment-list"
import { Lightbulb, Bug } from "lucide-react"
import { type CommentCategory } from "@/lib/supabase"

export default function CommentsPage() {
  const [activeTab, setActiveTab] = useState<CommentCategory>("feature_request")
  const [refreshTrigger, setRefreshTrigger] = useState(0)

  const handleCommentSubmitted = () => {
    // Trigger a refresh of the comment list
    setRefreshTrigger((prev) => prev + 1)
  }

  return (
    <div className="container mx-auto py-12 px-4">
      <div className="max-w-3xl mx-auto">
        <h1 className="text-3xl font-bold mb-2">Feature Requests & Bug Reports</h1>
        <p className="text-muted-foreground mb-8">
          Help us improve CASSIA by sharing your ideas or reporting issues. All submissions are reviewed before being published.
        </p>

        <Tabs value={activeTab} onValueChange={(v) => setActiveTab(v as CommentCategory)}>
          <TabsList className="grid w-full grid-cols-2 mb-6">
            <TabsTrigger value="feature_request" className="flex items-center gap-2">
              <Lightbulb className="h-4 w-4" />
              Feature Requests
            </TabsTrigger>
            <TabsTrigger value="bug_report" className="flex items-center gap-2">
              <Bug className="h-4 w-4" />
              Bug Reports
            </TabsTrigger>
          </TabsList>

          <TabsContent value="feature_request" className="space-y-8">
            <CommentForm category="feature_request" onSuccess={handleCommentSubmitted} />
            <CommentList category="feature_request" refreshTrigger={refreshTrigger} />
          </TabsContent>

          <TabsContent value="bug_report" className="space-y-8">
            <CommentForm category="bug_report" onSuccess={handleCommentSubmitted} />
            <CommentList category="bug_report" refreshTrigger={refreshTrigger} />
          </TabsContent>
        </Tabs>
      </div>
    </div>
  )
}
