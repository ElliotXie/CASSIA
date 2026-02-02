-- Free API key pool for R/Python package users without their own API keys
CREATE TABLE IF NOT EXISTS public.free_api_keys (
    id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
    provider VARCHAR(50) NOT NULL CHECK (provider IN ('google', 'together', 'openrouter')),
    encrypted_key TEXT NOT NULL,
    is_active BOOLEAN DEFAULT true,
    priority INTEGER DEFAULT 0,  -- Higher = preferred
    daily_limit INTEGER DEFAULT 1000,  -- Max calls per day per key
    current_daily_usage INTEGER DEFAULT 0,
    last_reset_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_free_api_keys_provider ON public.free_api_keys(provider, is_active);

-- RLS (only service role can access)
ALTER TABLE public.free_api_keys ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Service role can manage free API keys" ON public.free_api_keys
    FOR ALL USING (true);

-- updated_at trigger
CREATE OR REPLACE TRIGGER handle_updated_at_free_api_keys
    BEFORE UPDATE ON public.free_api_keys
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();
