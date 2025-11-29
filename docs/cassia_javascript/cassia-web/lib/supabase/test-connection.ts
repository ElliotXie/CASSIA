import { createClient } from '@/utils/supabase/client'

export async function testSupabaseConnection() {
  const supabase = createClient()
  const results = {
    connection: false,
    auth: false,
    tableAccess: false,
    user: null as any,
    error: null as any
  }

  try {
    // Test 1: Basic connection
    console.log('Testing Supabase connection...')
    const { data: { session } } = await supabase.auth.getSession()
    results.connection = true
    console.log('✓ Connection successful')

    // Test 2: Authentication
    if (session) {
      results.auth = true
      results.user = session.user
      console.log('✓ Authenticated as:', session.user.email)
    } else {
      console.log('✗ Not authenticated')
      return results
    }

    // Test 3: Table access
    console.log('Testing table access...')
    const { data, error } = await supabase
      .from('user_api_keys')
      .select('count')
      .eq('user_id', session.user.id)
      .single()

    if (error) {
      console.error('✗ Table access error:', error)
      results.error = error
    } else {
      results.tableAccess = true
      console.log('✓ Table access successful')
    }

  } catch (error) {
    console.error('Connection test error:', error)
    results.error = error
  }

  return results
}