-- Migration: Convert free API from daily rate limits to lifetime cluster limit
-- Old: 2 jobs/day, 30 clusters/job, per-provider, per-date
-- New: 2 cluster annotations per machine, lifetime, across all providers

-- Drop the old table (daily tracking is no longer needed)
DROP TABLE IF EXISTS public.free_api_usage;

-- Create new lifetime usage table
CREATE TABLE IF NOT EXISTS public.free_api_usage (
    id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
    machine_id VARCHAR(64) NOT NULL UNIQUE,
    clusters_used INTEGER DEFAULT 0,
    last_used_at TIMESTAMP WITH TIME ZONE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Index
CREATE INDEX IF NOT EXISTS idx_free_api_usage_machine ON public.free_api_usage(machine_id);

-- RLS
ALTER TABLE public.free_api_usage ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Allow usage tracking" ON public.free_api_usage
    FOR ALL USING (true);

-- updated_at trigger
CREATE OR REPLACE TRIGGER handle_updated_at_free_api_usage
    BEFORE UPDATE ON public.free_api_usage
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();
