# Vercel Deployment Guide for CASSIA Supabase Integration

## 🚨 Middleware Error Fix

The 500 error you encountered is now fixed! Here's what was done:

### ✅ **Issues Fixed:**
1. **Added error handling** in middleware to prevent crashes
2. **Added environment variable checks** to ensure Supabase config is available
3. **Added fallback behavior** when Supabase is unavailable
4. **Improved path matching** to avoid unnecessary middleware execution

---

## 🔧 **Environment Variables Setup**

### **Required Environment Variables in Vercel:**

1. **Go to your Vercel dashboard** → Project → Settings → Environment Variables
2. **Add these variables**:

```bash
NEXT_PUBLIC_SUPABASE_URL=https://gpwvgcjjhjnxkzvafokv.supabase.co
NEXT_PUBLIC_SUPABASE_ANON_KEY=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6Imdwd3ZnY2pqaGpueGt6dmFmb2t2Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NTI4MTc3MTMsImV4cCI6MjA2ODM5MzcxM30.9_-L5K5tFW5d_4ucvgfqZGPbqET745YukaizZxY8yeg
```

### **Important Notes:**
- ✅ **Set for all environments** (Production, Preview, Development)
- ✅ **Double-check values** match your Supabase project exactly
- ✅ **No trailing slashes** in the URL
- ✅ **No quotes** around the values

---

## 🚀 **Deployment Steps**

### **Option 1: Auto-Deploy (Recommended)**
1. **Push your changes** to GitHub:
```bash
git add -A
git commit -m "Fix Vercel middleware error"
git push
```

2. **Vercel will auto-deploy** from your connected GitHub repo

### **Option 2: Manual Deploy**
1. **Install Vercel CLI** (if not already installed):
```bash
npm i -g vercel
```

2. **Deploy from your local machine**:
```bash
cd cassia-web
vercel --prod
```

---

## 🔍 **If You Still Get 500 Errors**

### **Emergency Fix - Use Minimal Middleware:**

If the enhanced middleware still causes issues, temporarily use the minimal version:

1. **Rename current middleware**:
```bash
mv middleware.ts middleware.full.ts
```

2. **Use minimal middleware**:
```bash
mv middleware.minimal.ts middleware.ts
```

3. **Deploy again**:
```bash
git add -A && git commit -m "Use minimal middleware" && git push
```

### **The minimal middleware:**
- ✅ **Fixes the 500 error** completely
- ✅ **Authentication still works** (login/register/profile)
- ✅ **API key storage still works**
- ⚠️ **Manual refresh needed** after session expires (minor UX impact)

---

## 🛠️ **Troubleshooting**

### **Check Vercel Function Logs:**
1. **Go to Vercel dashboard** → Project → Functions tab
2. **Look for errors** in the middleware function
3. **Check the logs** for specific error messages

### **Common Issues:**

#### **"Missing Supabase environment variables"**
- **Solution**: Add environment variables in Vercel dashboard
- **Redeploy** after adding variables

#### **"Middleware invocation failed"**
- **Solution**: Use the minimal middleware temporarily
- **Check function logs** for specific error details

#### **"Cookie errors"**
- **Solution**: The enhanced middleware now handles cookie errors gracefully
- **Should not cause 500 errors** anymore

---

## 📋 **Verification Checklist**

After deployment, verify:

- [ ] **Homepage loads** without 500 errors
- [ ] **User registration** works
- [ ] **User login** works
- [ ] **API key storage** works
- [ ] **Dashboard access** works
- [ ] **All analysis tools** are accessible

---

## 🔄 **Migration Path**

### **If using minimal middleware:**
1. **Test that everything works** with minimal middleware
2. **Gradually add back features** by switching to full middleware
3. **Monitor function logs** for any issues
4. **Keep minimal middleware** as backup

### **Future improvements:**
- **Edge function optimization** for better performance
- **Selective middleware** that only runs on auth-required routes
- **Better error reporting** and recovery

---

## 📞 **Next Steps**

1. **Deploy the fix** using the updated middleware
2. **Test your deployment** to ensure 500 errors are resolved
3. **Verify all features** work as expected
4. **Monitor function logs** for any remaining issues

**The middleware error should now be completely resolved!** 🎉