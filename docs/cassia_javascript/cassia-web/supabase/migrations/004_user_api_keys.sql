-- User API keys table for encrypted API key storage
CREATE TABLE IF NOT EXISTS public.user_api_keys (
    id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    provider VARCHAR(50) NOT NULL CHECK (provider IN ('openrouter', 'anthropic', 'openai')),
    encrypted_key TEXT NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    CONSTRAINT unique_user_provider UNIQUE(user_id, provider)
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_user_api_keys_user_id ON public.user_api_keys(user_id);
CREATE INDEX IF NOT EXISTS idx_user_api_keys_provider ON public.user_api_keys(provider);

-- RLS
ALTER TABLE public.user_api_keys ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Users can view own API keys" ON public.user_api_keys
    FOR SELECT USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own API keys" ON public.user_api_keys
    FOR INSERT WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own API keys" ON public.user_api_keys
    FOR UPDATE USING (auth.uid() = user_id);

CREATE POLICY "Users can delete own API keys" ON public.user_api_keys
    FOR DELETE USING (auth.uid() = user_id);

-- updated_at trigger
CREATE OR REPLACE TRIGGER handle_updated_at_user_api_keys
    BEFORE UPDATE ON public.user_api_keys
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();
