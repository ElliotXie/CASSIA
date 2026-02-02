-- Analysis results table
CREATE TABLE IF NOT EXISTS public.analysis_results (
    id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
    user_id UUID REFERENCES auth.users(id) ON DELETE CASCADE,
    analysis_type VARCHAR(50) NOT NULL CHECK (analysis_type IN ('batch', 'symphony', 'scoring', 'subclustering', 'annotation-boost')),
    title VARCHAR(255),
    description TEXT,
    input_data JSONB,
    results JSONB,
    settings JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_analysis_results_user_id ON public.analysis_results(user_id);
CREATE INDEX IF NOT EXISTS idx_analysis_results_type ON public.analysis_results(analysis_type);
CREATE INDEX IF NOT EXISTS idx_analysis_results_created_at ON public.analysis_results(created_at DESC);

-- RLS
ALTER TABLE public.analysis_results ENABLE ROW LEVEL SECURITY;

CREATE POLICY "Users can view own analysis results" ON public.analysis_results
    FOR SELECT USING (auth.uid() = user_id);

CREATE POLICY "Users can insert own analysis results" ON public.analysis_results
    FOR INSERT WITH CHECK (auth.uid() = user_id);

CREATE POLICY "Users can update own analysis results" ON public.analysis_results
    FOR UPDATE USING (auth.uid() = user_id);

CREATE POLICY "Users can delete own analysis results" ON public.analysis_results
    FOR DELETE USING (auth.uid() = user_id);

-- updated_at trigger
CREATE OR REPLACE TRIGGER handle_updated_at_analysis_results
    BEFORE UPDATE ON public.analysis_results
    FOR EACH ROW EXECUTE PROCEDURE public.handle_updated_at();
