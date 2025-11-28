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
      <h1 className="text-3xl font-bold mb-8">评论</h1>

      <div className="max-w-2xl mx-auto">
        {/* Comment form */}
        <Card>
          <CardHeader>
            <CardTitle>留下评论</CardTitle>
            <CardDescription>
              分享您对CASSIA的想法、问题或反馈。所有评论将直接发送给CASSIA团队。
            </CardDescription>
          </CardHeader>
          <CardContent>
            {isSubmitted ? (
              <div className="bg-green-50 dark:bg-green-900/20 p-4 rounded-md text-green-800 dark:text-green-300">
                <p className="font-medium">感谢您的评论！</p>
                <p className="mt-1">您的信息已发送给CASSIA团队。我们感谢您的反馈。</p>
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
                  <Label htmlFor="name">姓名（可选）</Label>
                  <Input id="name" placeholder="您的姓名" value={formData.name} onChange={handleChange} />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="subject">主题</Label>
                  <Input
                    id="subject"
                    placeholder="评论主题"
                    required
                    value={formData.subject}
                    onChange={handleChange}
                  />
                </div>
                <div className="space-y-2">
                  <Label htmlFor="comment">评论</Label>
                  <Textarea
                    id="comment"
                    placeholder="在此处写下您的评论..."
                    className="min-h-[120px]"
                    required
                    value={formData.comment}
                    onChange={handleChange}
                  />
                </div>
                <Button type="submit" disabled={isSubmitting}>
                  {isSubmitting ? "发送中..." : "提交评论"}
                </Button>
              </form>
            )}
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
