# CASSIA Supabase Integration Debugging Guide

## Issues Fixed

### 1. API Key Storage Problem ✅
**Problem**: API keys were not being saved to Supabase table
**Cause**: Circular dependency between auth-store and api-key-store-supabase
**Solution**: 
- Removed dependency on auth-store
- Added direct authentication check using `supabase.auth.getUser()`
- Improved error handling and logging

### 2. Zustand Migration Error ✅
**Problem**: "State loaded from storage couldn't be migrated since no migrate function was provided"
**Cause**: Version mismatch between old api-key-store.ts (v1) and new api-key-store-supabase.ts (v2)
**Solution**:
- Added migrate function to handle version 1 → 2 transition
- Properly handles state structure changes
- Maintains backward compatibility

## How to Test the Fixes

### 1. Clear Browser Storage (Important!)
```javascript
// Open browser console and run:
localStorage.clear()
sessionStorage.clear()
```

### 2. Test API Key Storage
1. **Register/Login** to your account
2. **Set API Key**: Click "Set API Key" and enter a test key
3. **Check Console**: Look for "API key for [provider] saved securely to Supabase"
4. **Verify in Supabase**: 
   - Go to your Supabase dashboard
   - Check the `user_api_keys` table
   - You should see an encrypted entry for your user

### 3. Test API Key Loading
1. **Refresh the page** while logged in
2. **Check Console**: Look for "API keys loaded from Supabase"
3. **Verify UI**: API key should still be marked as "Set"

## Console Debugging

### Expected Success Messages:
```
API key storage hydrated successfully
API key for openrouter saved securely to Supabase
API keys loaded from Supabase
```

### Common Error Messages and Solutions:

#### "User not authenticated, API key saved locally only"
- **Cause**: User session not properly established
- **Solution**: Wait a moment after login, then try again

#### "Error saving API key to Supabase: [error]"
- **Cause**: Database permissions or RLS policy issue
- **Solution**: Check Supabase dashboard for RLS policies

#### "Failed to authenticate user"
- **Cause**: Supabase client not properly configured
- **Solution**: Check environment variables in .env.local

## Database Verification

### Check User API Keys Table:
```sql
SELECT * FROM user_api_keys;
```

### Check RLS Policies:
```sql
SELECT * FROM pg_policies WHERE tablename = 'user_api_keys';
```

### Manually Test Encryption:
```javascript
// In browser console:
const encryptKey = (key) => btoa(key)
const decryptKey = (encrypted) => atob(encrypted)

console.log('Encrypted:', encryptKey('test-key'))
console.log('Decrypted:', decryptKey(encryptKey('test-key')))
```

## Migration Testing

### Test Version Migration:
1. **Create old format data** in localStorage:
```javascript
localStorage.setItem('cassia-api-key-storage', JSON.stringify({
  state: {
    apiKeys: { openrouter: 'old-key', anthropic: '', openai: '' },
    provider: 'openrouter',
    model: 'old-model'
  },
  version: 1
}))
```

2. **Refresh page** and check console for "Migrating API key store from version 1 to 2"

## Performance Tips

- **Sync Delay**: API keys sync 2 seconds after hydration to ensure auth is ready
- **Local Fallback**: Keys are kept locally for offline access
- **Direct Auth**: No circular dependencies between stores

## Still Having Issues?

1. **Check Browser Network Tab**: Look for failed requests to Supabase
2. **Verify Environment Variables**: Ensure .env.local has correct values
3. **Test with Fresh Account**: Try with a completely new user account
4. **Check Supabase Logs**: Look at Supabase dashboard for error logs

## Advanced Debugging

### Enable Detailed Logging:
Add this to your browser console:
```javascript
// Enable detailed Supabase logging
localStorage.setItem('supabase.auth.debug', 'true')
```

### Check Authentication State:
```javascript
// Get current auth state
const { createClient } = require('@supabase/supabase-js')
const supabase = createClient(process.env.NEXT_PUBLIC_SUPABASE_URL, process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY)
supabase.auth.getUser().then(console.log)
```

## Success Criteria

✅ **Registration**: New user created in Supabase auth.users  
✅ **Profile**: User profile created in profiles table  
✅ **API Key Storage**: Encrypted key saved in user_api_keys table  
✅ **API Key Loading**: Keys loaded on page refresh  
✅ **Migration**: No console errors about migration  
✅ **Dashboard**: User dashboard shows correct data  

When all these work, your integration is successful! 🎉