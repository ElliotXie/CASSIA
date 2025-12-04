# Supabase Email Confirmation Configuration Guide

## Problem Fixed
Email confirmation links were redirecting to `localhost:3000` instead of your production domain, causing "otp_expired" errors.

## Solution Implemented

### 1. Code Changes (Already Done ✅)

**Updated `lib/stores/auth-store.ts`:**
- Added `emailRedirectTo` option to signUp call (line 106)
- Now redirects to `${window.location.origin}/auth/callback`
- Works for both local dev (localhost:3000) and production (Vercel)

**Created `app/auth/callback/route.ts`:**
- Handles email confirmation token exchange
- Exchanges code from email link for a session
- Redirects users to home page after successful confirmation

### 2. Supabase Dashboard Configuration (REQUIRED - Manual Steps)

You **MUST** configure these settings in your Supabase dashboard for email confirmation to work properly.

#### Step 1: Access Authentication Settings

1. Go to: https://supabase.com/dashboard/project/fauimshrishydpnactfc/auth/url-configuration
2. Log in to your Supabase account
3. Navigate to **Authentication** → **URL Configuration**

#### Step 2: Configure Site URL

**Important:** Set this to your PRIMARY production domain

**For Vercel Deployment:**
```
Site URL: https://your-vercel-app-name.vercel.app
```

**To find your Vercel URL:**
1. Go to https://vercel.com/dashboard
2. Find your CASSIA project
3. Copy the deployment URL (e.g., `https://cassia-web-xyz123.vercel.app`)
4. Use that as your Site URL in Supabase

**For local development only:**
```
Site URL: http://localhost:3000
```
(But you'll need to change this when deploying to production)

#### Step 3: Add Redirect URLs

Add these URLs to the **Redirect URLs** allowlist (all of them):

**For Production (Vercel):**
```
https://your-vercel-app-name.vercel.app/auth/callback
https://your-vercel-app-name.vercel.app/**
```

**For Local Development:**
```
http://localhost:3000/auth/callback
http://localhost:3000/**
```

**Example (replace with your actual Vercel URL):**
```
https://cassia-web-abc123.vercel.app/auth/callback
https://cassia-web-abc123.vercel.app/**
http://localhost:3000/auth/callback
http://localhost:3000/**
```

The `/**` wildcard allows all routes under your domain.

#### Step 4: Save Configuration

Click **Save** at the bottom of the URL Configuration page.

### 3. Verify Vercel Environment Variables

Ensure these are set in your Vercel project settings:

1. Go to https://vercel.com/dashboard
2. Select your CASSIA project
3. Go to **Settings** → **Environment Variables**
4. Verify these exist:
   ```
   NEXT_PUBLIC_SUPABASE_URL=https://fauimshrishydpnactfc.supabase.co
   NEXT_PUBLIC_SUPABASE_ANON_KEY=[your-anon-key]
   ```

If missing, add them and redeploy.

## Testing the Fix

### Test Email Confirmation Flow

1. **Register a new account** with a fresh email address
2. **Check your email** - you should receive a confirmation email
3. **Click the "Confirm your mail" link** in the email
4. **Expected behavior:**
   - Link should redirect to your Vercel URL (not localhost)
   - URL should be: `https://your-app.vercel.app/auth/callback?code=...`
   - After successful confirmation: redirects to `https://your-app.vercel.app/?confirmed=true`
   - You can now sign in with your email and password

### If You Still Get Errors

**Error: "otp_expired" or "Email link is invalid"**
- The email link may have actually expired (they expire after 1 hour by default)
- Try registering with a new email and clicking the link within a few minutes

**Error: "access_denied"**
- Check that your Vercel URL is correctly added to Redirect URLs in Supabase
- Ensure Site URL matches your primary domain
- Make sure you clicked "Save" in Supabase dashboard

**Error: Link still goes to localhost**
- Clear your browser cache
- Make sure the latest code is deployed to Vercel
- Check that `emailRedirectTo` is in auth-store.ts (line 106)

## Configuration Checklist

- [ ] Found your Vercel production URL
- [ ] Set Supabase Site URL to Vercel domain
- [ ] Added production callback URL to Redirect URLs (`https://your-app.vercel.app/auth/callback`)
- [ ] Added production wildcard to Redirect URLs (`https://your-app.vercel.app/**`)
- [ ] Added local callback URL for development (`http://localhost:3000/auth/callback`)
- [ ] Added local wildcard for development (`http://localhost:3000/**`)
- [ ] Saved configuration in Supabase dashboard
- [ ] Verified environment variables in Vercel
- [ ] Deployed latest code to Vercel
- [ ] Tested registration with a fresh email

## Important Notes

### Email Link Expiration
- Supabase email confirmation links expire after 1 hour by default
- If testing, register and click the link within a few minutes
- To change expiration time: Supabase Dashboard → Authentication → Email Settings

### Development vs Production
- **Local development:** Use `http://localhost:3000` as Site URL
- **Production:** Use your Vercel URL as Site URL
- You can have both URLs in the Redirect URLs list simultaneously
- The `emailRedirectTo` uses `window.location.origin` so it adapts automatically

### Custom Domains
If you add a custom domain to Vercel (e.g., `cassia.yourdomain.com`):
1. Update Supabase Site URL to your custom domain
2. Add custom domain callback to Redirect URLs
3. Keep Vercel default URL in Redirect URLs as backup

## Quick Reference

**Your Supabase Project:**
- Project URL: https://fauimshrishydpnactfc.supabase.co
- Dashboard: https://supabase.com/dashboard/project/fauimshrishydpnactfc

**Required Dashboard Pages:**
- URL Configuration: https://supabase.com/dashboard/project/fauimshrishydpnactfc/auth/url-configuration
- Email Templates: https://supabase.com/dashboard/project/fauimshrishydpnactfc/auth/templates

**Your Deployment:**
- Platform: Vercel
- Find URL at: https://vercel.com/dashboard
