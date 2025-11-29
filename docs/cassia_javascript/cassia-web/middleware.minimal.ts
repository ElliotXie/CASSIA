import { type NextRequest, NextResponse } from 'next/server'

// Minimal middleware - use this if the full Supabase middleware causes issues
export async function middleware(request: NextRequest) {
  // Just pass through all requests without Supabase session handling
  return NextResponse.next()
}

export const config = {
  matcher: [
    // Only match specific paths that need middleware
    '/((?!api|_next/static|_next/image|favicon.ico|.*\\.(?:svg|png|jpg|jpeg|gif|webp|ico|css|js)$).*)',
  ],
}

// Note: This minimal version disables automatic session refresh
// Authentication will still work, but users might need to refresh
// the page manually after their session expires