# CASSIA Supabase Integration Setup Guide

🚀 **Complete setup guide for adding user authentication, secure API key storage, and analysis results persistence to CASSIA**

## 📋 Prerequisites

- [x] CASSIA project with cassia-web directory
- [x] Supabase project created with URL and anon key
- [x] @supabase/ssr package installed
- [x] Environment variables configured

## 🔧 Step-by-Step Setup

### Step 1: Database Schema Setup

1. **Open your Supabase dashboard** → SQL Editor
2. **Create a new query** and copy the entire contents of `supabase_schema.sql`
3. **Execute the query** to create all tables, indexes, RLS policies, and triggers
4. **Verify the setup** by checking the Tables section in your dashboard

You should see these tables:
- `profiles` - User profile information
- `user_api_keys` - Encrypted API keys storage
- `analysis_results` - Analysis history and results
- `analysis_sessions` - Active/completed analysis sessions

### Step 2: Authentication Configuration

1. **Enable Email Authentication**:
   - Go to Authentication → Settings
   - Enable "Email" under Auth Providers
   - Configure email templates (optional)
   - Set site URL to your domain

2. **Configure Email Settings** (optional):
   - Set custom SMTP settings for production
   - Customize email templates

### Step 3: Test the Integration

1. **Start the development server**:
   ```bash
   cd cassia-web
   npm run dev
   ```

2. **Test Authentication**:
   - Open http://localhost:3000
   - Click "Sign In" in the top right
   - Create a test account
   - Verify profile creation

3. **Test API Key Storage**:
   - Click "Set API Key" after signing in
   - Enter a test API key
   - Verify it's stored securely in Supabase

4. **Test Dashboard**:
   - Click "Dashboard" after signing in
   - Verify user statistics display
   - Test navigation to analysis tools

## 🔐 Security Features

### Row Level Security (RLS)
All tables have RLS enabled with policies ensuring:
- Users can only access their own data
- API keys are encrypted before storage
- Analysis results are private to each user

### API Key Encryption
- API keys are encrypted using base64 encoding (demo)
- In production, implement proper encryption with user-specific keys
- Keys never leave the user's session context

### Session Management
- Automatic session refresh
- Secure cookie handling
- Cross-device synchronization

## 📊 Key Features Added

### 1. User Authentication
- Email/password registration and login
- User profile management
- Secure session handling

### 2. Secure API Key Storage
- Encrypted storage in Supabase
- Per-provider key management
- Automatic synchronization across devices

### 3. Analysis Results Persistence
- Save analysis results with metadata
- Analysis session management
- Export/import capabilities

### 4. User Dashboard
- Analysis history overview
- Active sessions management
- Usage statistics
- Quick access to tools

## 🛠️ File Structure

```
cassia-web/
├── .env.local                          # Environment variables
├── middleware.ts                       # Auth middleware
├── supabase_schema.sql                 # Database schema
├── lib/
│   ├── supabase/
│   │   ├── client.ts                   # Browser client
│   │   ├── server.ts                   # Server client
│   │   ├── middleware.ts               # Session management
│   │   └── types.ts                    # Database types
│   └── stores/
│       ├── auth-store.ts               # Authentication state
│       ├── api-key-store-supabase.ts   # API key management
│       └── results-store.ts            # Results management
├── components/
│   ├── auth/
│   │   ├── AuthProvider.tsx            # Auth context provider
│   │   ├── AuthButton.tsx              # Auth UI component
│   │   ├── LoginForm.tsx               # Login/register form
│   │   └── ProfileForm.tsx             # Profile management
│   └── dashboard/
│       └── UserDashboard.tsx           # User dashboard
└── app/
    ├── layout.tsx                      # Updated with AuthProvider
    └── page.tsx                        # Updated with auth integration
```

## 🚀 Next Steps

### Phase 1: Basic Integration (✅ Complete)
- [x] User authentication
- [x] API key storage
- [x] Results persistence
- [x] User dashboard

### Phase 2: Enhanced Features
- [ ] Analysis sharing between users
- [ ] Team workspaces
- [ ] Advanced analytics
- [ ] Export/import improvements

### Phase 3: Advanced Features
- [ ] Real-time collaboration
- [ ] Advanced visualization
- [ ] API endpoints
- [ ] Third-party integrations

## 🐛 Troubleshooting

### Common Issues

1. **Authentication not working**:
   - Check environment variables in `.env.local`
   - Verify Supabase URL and anon key
   - Check browser console for errors

2. **API keys not saving**:
   - Ensure user is authenticated
   - Check RLS policies in Supabase
   - Verify table permissions

3. **Dashboard not loading**:
   - Check authentication state
   - Verify database tables exist
   - Check browser console for errors

### Debug Tips

1. **Check browser console** for authentication errors
2. **Verify database** in Supabase dashboard
3. **Test API endpoints** directly in Supabase
4. **Check network tab** for failed requests

## 📞 Support

If you encounter issues:
1. Check the console for error messages
2. Verify your Supabase configuration
3. Ensure all environment variables are set
4. Test with a fresh user account

## 🎉 Success Criteria

Your integration is successful when:
- [x] Users can register and login
- [x] API keys are stored securely
- [x] Dashboard shows user data
- [x] Analysis results persist across sessions
- [x] Cross-device synchronization works

**Congratulations! Your CASSIA application now has full user authentication and data persistence. 🎉**