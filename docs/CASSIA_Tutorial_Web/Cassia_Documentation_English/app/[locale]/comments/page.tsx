"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Textarea } from "@/components/ui/textarea"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { useState } from "react"
import { AlertCircle } from "lucide-react"

export default function CommentsPage() {
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [isSubmitted, setIsSubmitted] = useState(false)
  const [error, setError] = useState("")
  const [formData, setFormData] = useState({
    name: "",
    subject: "",
    comment: "",
  })

  const handleChange = (e) => {
    const { id, value } = e.target
    setFormData((prev) => ({
      ...prev,
      [id]: value,
    }))
  }

  const handleSubmit = async (event) => {
    event.preventDefault()
    setIsSubmitting(true)
    setError("")

    try {
      // Using Formspree endpoint
      const response = await fetch("https://formspree.io/f/xnndqlvb", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify(formData),
      })

      if (!response.ok) {
        throw new Error("Failed to send comment. Please try again.")
      }

      setIsSubmitted(true)
    } catch (err) {
      setError(err.message || "Something went wrong. Please try again.")
    } finally {
      setIsSubmitting(false)
    }
  }

  return (
    <div className="container mx-auto py-12 px-4">
      <h1 className="text-3xl font-bold mb-8">Comments</h1>

      <div className="max-w-2xl mx-auto">
        {/* Comment form */}
        <Card>
          <CardHeader>
            <CardTitle>Leave a comment</CardTitle>
            <CardDescription>
              Share your thoughts, questions, or feedback about CASSIA. All comments will be sent directly to the CASSIA
              team.
            </CardDescription>
          </CardHeader>
          <CardContent>
            {isSubmitted ? (
              <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-md text-green-800 dark:text-green-300">
                <p className="font-medium">Thank you for your comment!</p>
                <p className="mt-1">Your message has been sent to the CASSIA team. We appreciate your feedback.</p>
              </div>
            ) : (
              <form className="space-y-4" onSubmit={handleSubmit} method="POST">
                {error && (
                  <div className="bg-red-50 dark:bg-red-900/20 p-4 rounded-md text-red-800 dark:text-red-300 flex items-start gap-3">
                    <AlertCircle className="h-5 w-5 mt-0.5" />
                    <p>{error}</p>
                  </div>
                )}
                <div className="space-y-2">
                  <Label htmlFor="name">Name (optional)</Label>
                  <Input id="name" placeholder="Your name" value={formData.name} onChange={handleChange} />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="subject">Subject</Label>
                  <Input
                    id="subject"
                    placeholder="Comment subject"
                    required
                    value={formData.subject}
                    onChange={handleChange}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="comment">Comment</Label>
                  <Textarea
                    id="comment"
                    placeholder="Write your comment here..."
                    className="min-h-[120px]"
                    required
                    value={formData.comment}
                    onChange={handleChange}
                  />
                </div>
                <Button type="submit" disabled={isSubmitting}>
                  {isSubmitting ? "Sending..." : "Submit Comment"}
                </Button>
              </form>
            )}
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
