import { NextRequest, NextResponse } from 'next/server'
import { createClient } from '@supabase/supabase-js'

export async function GET(request: NextRequest) {
  try {
    const { searchParams } = new URL(request.url)
    const token = searchParams.get('token')

    if (!token) {
      return NextResponse.json(
        { error: 'Approval token is required' },
        { status: 400 }
      )
    }

    // Use service role key for admin operations
    const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL
    const serviceRoleKey = process.env.SUPABASE_SERVICE_ROLE_KEY

    if (!supabaseUrl || !serviceRoleKey) {
      console.error('[Comment Approval] Missing Supabase configuration')
      return NextResponse.json(
        { error: 'Server configuration error' },
        { status: 500 }
      )
    }

    const supabaseAdmin = createClient(supabaseUrl, serviceRoleKey)

    // Find comment by approval token
    const { data: comment, error: findError } = await supabaseAdmin
      .from('feedback_comments')
      .select('*')
      .eq('approval_token', token)
      .single()

    if (findError || !comment) {
      return NextResponse.json(
        { error: 'Comment not found or invalid token' },
        { status: 404 }
      )
    }

    if (comment.is_approved) {
      // Already approved - show success page
      return new NextResponse(getSuccessHtml(comment, true), {
        status: 200,
        headers: { 'Content-Type': 'text/html' }
      })
    }

    // Approve the comment
    const { error: updateError } = await supabaseAdmin
      .from('feedback_comments')
      .update({ is_approved: true })
      .eq('id', comment.id)

    if (updateError) throw updateError

    // Return success HTML page
    return new NextResponse(getSuccessHtml(comment, false), {
      status: 200,
      headers: { 'Content-Type': 'text/html' }
    })

  } catch (error) {
    console.error('[Comment Approval] Error:', error)
    return NextResponse.json(
      { error: error instanceof Error ? error.message : 'Approval failed' },
      { status: 500 }
    )
  }
}

function getSuccessHtml(comment: { comment_type: string; message: string }, alreadyApproved: boolean): string {
  const typeLabel = comment.comment_type === 'feature_request' ? 'Feature Request' : 'Bug Report'
  const typeClass = comment.comment_type === 'feature_request' ? 'feature' : 'bug'
  const message = alreadyApproved ? 'This comment was already approved.' : 'The feedback has been approved and is now visible on the website.'

  return `
    <!DOCTYPE html>
    <html>
      <head>
        <title>Comment Approved - CASSIA</title>
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <style>
          body {
            font-family: system-ui, -apple-system, sans-serif;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
            margin: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          }
          .card {
            background: white;
            padding: 2rem;
            border-radius: 12px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            text-align: center;
            max-width: 400px;
            margin: 1rem;
          }
          .success-icon {
            font-size: 48px;
            margin-bottom: 1rem;
          }
          h1 { color: #10b981; margin: 0 0 1rem; }
          p { color: #666; margin: 0; }
          .comment-preview {
            margin-top: 1rem;
            padding: 1rem;
            background: #f3f4f6;
            border-radius: 8px;
            text-align: left;
          }
          .type-badge {
            display: inline-block;
            padding: 0.25rem 0.5rem;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 600;
            margin-bottom: 0.5rem;
          }
          .feature { background: #dbeafe; color: #1d4ed8; }
          .bug { background: #fee2e2; color: #dc2626; }
          .message-text {
            margin-top: 0.5rem;
            color: #374151;
            word-break: break-word;
          }
        </style>
      </head>
      <body>
        <div class="card">
          <div class="success-icon">${alreadyApproved ? 'ℹ️' : '✅'}</div>
          <h1>${alreadyApproved ? 'Already Approved' : 'Comment Approved!'}</h1>
          <p>${message}</p>
          <div class="comment-preview">
            <span class="type-badge ${typeClass}">
              ${typeLabel}
            </span>
            <p class="message-text">${escapeHtml(comment.message.substring(0, 300))}${comment.message.length > 300 ? '...' : ''}</p>
          </div>
        </div>
      </body>
    </html>
  `
}

function escapeHtml(text: string): string {
  const div = { innerHTML: '' }
  return text
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#039;')
}
