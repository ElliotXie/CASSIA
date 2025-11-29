import { createServerClient } from '@supabase/ssr'
import { NextResponse, type NextRequest } from 'next/server'
import { getSupabaseConfig, isSupabaseConfigured } from './config'

export async function updateSession(request: NextRequest) {
  try {
    // Check if environment variables are available
    if (!isSupabaseConfigured()) {
      console.error('Missing Supabase environment variables')
      return NextResponse.next()
    }

    const { url, anonKey } = getSupabaseConfig()

    let supabaseResponse = NextResponse.next({
      request,
    })

    const supabase = createServerClient(
      url,
      anonKey,
      {
        cookies: {
          getAll() {
            return request.cookies.getAll()
          },
          setAll(cookiesToSet) {
            cookiesToSet.forEach(({ name, value, options }) => {
              try {
                request.cookies.set(name, value)
              } catch (error) {
                console.error(`Error setting cookie ${name}:`, error)
              }
            })
            supabaseResponse = NextResponse.next({
              request,
            })
            cookiesToSet.forEach(({ name, value, options }) => {
              try {
                supabaseResponse.cookies.set(name, value, options)
              } catch (error) {
                console.error(`Error setting response cookie ${name}:`, error)
              }
            })
          },
        },
      }
    )

    // This will refresh session if expired - required for Server Components
    try {
      await supabase.auth.getUser()
    } catch (error) {
      console.error('Error getting user in middleware:', error)
    }

    return supabaseResponse
  } catch (error) {
    console.error('Middleware updateSession error:', error)
    return NextResponse.next()
  }
}