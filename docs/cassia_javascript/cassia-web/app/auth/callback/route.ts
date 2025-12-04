import { createClient } from '@/utils/supabase/server'
import { NextResponse } from 'next/server'
import { cookies } from 'next/headers'

export async function GET(request: Request) {
  const requestUrl = new URL(request.url)
  const code = requestUrl.searchParams.get('code')
  const origin = requestUrl.origin

  if (code) {
    const cookieStore = await cookies()
    const supabase = await createClient()

    try {
      // Exchange the code for a session
      const { error } = await supabase.auth.exchangeCodeForSession(code)

      if (error) {
        console.error('Error exchanging code for session:', error)
        // Redirect to home with error
        return NextResponse.redirect(`${origin}/?error=auth_callback_error`)
      }

      // Successful email confirmation - redirect to home
      return NextResponse.redirect(`${origin}/?confirmed=true`)
    } catch (error) {
      console.error('Exception in auth callback:', error)
      return NextResponse.redirect(`${origin}/?error=auth_callback_exception`)
    }
  }

  // No code present - redirect to home
  return NextResponse.redirect(`${origin}/`)
}
