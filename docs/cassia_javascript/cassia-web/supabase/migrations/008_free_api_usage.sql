-- Usage tracking per machine_id (hashed fingerprint) for free API rate limiting
CREATE TABLE IF NOT EXISTS public.free_api_usage (
    id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
    machine_id VARCHAR(64) NOT NULL,
    provider VARCHAR(50) NOT NULL,
    job_count INTEGER DEFAULT 0,  -- Jobs today
    cluster_count INTEGER DEFAULT 0,  -- Total clusters today
    last_job_at TIMESTAMP WITH TIME ZONE,
    usage_date DATE DEFAULT CURRENT_DATE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    CONSTRAINT unique_machine_provider_date UNIQUE(machine_id, provider, usage_date)
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_free_api_usage_machine ON public.free_api_usage(machine_id, usage_date);

-- RLS (allow access for usage tracking)
ALTER TABLE public.free_api_usage ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Allow usage tracking" ON public.free_api_usage
    FOR ALL USING (true);

-- updated_at trigger
CREATE OR REPLACE TRIGGER handle_updated_at_free_api_usage
    BEFORE UPDATE ON public.free_api_usage
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();
