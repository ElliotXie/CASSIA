'use client'

import React, { useState } from 'react'
import { Mail, Github, ExternalLink, HelpCircle, Send, CheckCircle, AlertCircle } from 'lucide-react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Textarea } from '@/components/ui/textarea'
import { Label } from '@/components/ui/label'
import emailjs from '@emailjs/browser'

export function ContactDialog() {
  const [formData, setFormData] = useState({
    name: '',
    email: '',
    subject: '',
    message: ''
  })
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [submitStatus, setSubmitStatus] = useState<'idle' | 'success' | 'error'>('idle')
  const [showEmailForm, setShowEmailForm] = useState(false)

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target
    setFormData(prev => ({ ...prev, [name]: value }))
  }

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault()
    setIsSubmitting(true)
    setSubmitStatus('idle')

    try {
      // Configure with your EmailJS credentials
      const result = await emailjs.send(
        'service_m58e32j',
        'template_po0safd',
        {
          name: formData.name,
          email: formData.email,
          title: formData.subject,
          message: formData.message
        },
        'qS1BL05ocjba4E3nE'
      )

      if (result.text === 'OK') {
        setSubmitStatus('success')
        setFormData({ name: '', email: '', subject: '', message: '' })
        setTimeout(() => setShowEmailForm(false), 2000)
      }
    } catch (error) {
      console.error('Email send error:', error)
      setSubmitStatus('error')
    } finally {
      setIsSubmitting(false)
    }
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center space-x-2">
          <HelpCircle className="h-5 w-5" />
          <span>Get Help & Support</span>
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Email Contact */}
        <div className="space-y-2">
          <label className="text-sm font-medium">Email Support</label>
          <div className="p-3 border rounded-lg bg-muted/30">
            <div className="flex items-center justify-between">
              <div className="flex items-center space-x-2">
                <Mail className="h-4 w-4 text-primary" />
                <span className="text-sm">xie227@wisc.edu</span>
              </div>
              <Button
                variant="outline"
                size="sm"
                onClick={() => setShowEmailForm(!showEmailForm)}
              >
                <Send className="h-3 w-3 mr-1" />
                Send Email
              </Button>
            </div>
            <p className="text-xs text-muted-foreground mt-2">
              Direct email support for technical questions and assistance
            </p>
            
            {/* Contact Form */}
            {showEmailForm && (
              <div className="mt-4 p-4 border rounded-lg bg-background">
                <form onSubmit={handleSubmit} className="space-y-4">
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <Label htmlFor="name">Name</Label>
                      <Input
                        id="name"
                        name="name"
                        value={formData.name}
                        onChange={handleInputChange}
                        required
                        placeholder="Your name"
                      />
                    </div>
                    <div>
                      <Label htmlFor="email">Email</Label>
                      <Input
                        id="email"
                        name="email"
                        type="email"
                        value={formData.email}
                        onChange={handleInputChange}
                        required
                        placeholder="your.email@example.com"
                      />
                    </div>
                  </div>
                  <div>
                    <Label htmlFor="subject">Subject</Label>
                    <Input
                      id="subject"
                      name="subject"
                      value={formData.subject}
                      onChange={handleInputChange}
                      required
                      placeholder="Brief description of your issue"
                    />
                  </div>
                  <div>
                    <Label htmlFor="message">Message</Label>
                    <Textarea
                      id="message"
                      name="message"
                      value={formData.message}
                      onChange={handleInputChange}
                      required
                      placeholder="Describe your issue in detail..."
                      rows={4}
                    />
                  </div>
                  <div className="flex items-center justify-between">
                    <Button
                      type="submit"
                      disabled={isSubmitting}
                      className="flex items-center space-x-2"
                    >
                      {isSubmitting ? (
                        <>
                          <div className="animate-spin h-4 w-4 border-2 border-white border-t-transparent rounded-full" />
                          <span>Sending...</span>
                        </>
                      ) : (
                        <>
                          <Send className="h-4 w-4" />
                          <span>Send Message</span>
                        </>
                      )}
                    </Button>
                    {submitStatus === 'success' && (
                      <div className="flex items-center text-green-600 text-sm">
                        <CheckCircle className="h-4 w-4 mr-1" />
                        Email sent successfully!
                      </div>
                    )}
                    {submitStatus === 'error' && (
                      <div className="flex items-center text-red-600 text-sm">
                        <AlertCircle className="h-4 w-4 mr-1" />
                        Failed to send. Please try email client.
                      </div>
                    )}
                  </div>
                </form>
              </div>
            )}
          </div>
        </div>

        {/* GitHub Issues */}
        <div className="space-y-2">
          <label className="text-sm font-medium">GitHub Issues</label>
          <div className="p-3 border rounded-lg bg-muted/30">
            <div className="flex items-center justify-between">
              <div className="flex items-center space-x-2">
                <Github className="h-4 w-4 text-primary" />
                <span className="text-sm">Report Issues & Bugs</span>
              </div>
              <Button
                variant="outline"
                size="sm"
                onClick={() => window.open('https://github.com/your-repo/cassia/issues', '_blank')}
              >
                <ExternalLink className="h-3 w-3 mr-1" />
                Open Issues
              </Button>
            </div>
            <p className="text-xs text-muted-foreground mt-2">
              Report bugs, request features, or browse existing issues
            </p>
          </div>
        </div>

        {/* Help Information */}
        <div className="p-3 rounded-lg bg-blue-50 border border-blue-200">
          <div className="text-sm font-medium text-blue-800 mb-2">💡 Before contacting support:</div>
          <ul className="text-xs text-blue-700 space-y-1">
            <li>• Check the documentation and examples</li>
            <li>• Verify your API key is correctly configured</li>
            <li>• Try refreshing the page or clearing browser cache</li>
            <li>• Include error messages and steps to reproduce issues</li>
          </ul>
        </div>
      </CardContent>
    </Card>
  )
}