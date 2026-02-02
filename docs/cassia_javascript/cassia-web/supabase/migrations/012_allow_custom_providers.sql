-- Allow custom provider keys in user_api_keys table
-- The original CHECK constraint only allowed: 'openrouter', 'anthropic', 'openai'
-- Custom providers use the format 'custom:<preset>' (e.g., 'custom:deepseek', 'custom:qwen')
--
-- This migration dynamically finds and drops ALL check constraints on the table
-- regardless of their auto-generated name, then adds the correct permissive one.

DO $$
DECLARE
    r RECORD;
BEGIN
    FOR r IN
        SELECT con.conname
        FROM pg_constraint con
        JOIN pg_class rel ON rel.oid = con.conrelid
        JOIN pg_namespace nsp ON nsp.oid = rel.relnamespace
        WHERE rel.relname = 'user_api_keys'
        AND nsp.nspname = 'public'
        AND con.contype = 'c'
    LOOP
        EXECUTE 'ALTER TABLE public.user_api_keys DROP CONSTRAINT ' || quote_ident(r.conname);
    END LOOP;
END $$;

ALTER TABLE public.user_api_keys
  ADD CONSTRAINT user_api_keys_provider_check
  CHECK (
    provider IN ('openrouter', 'anthropic', 'openai', 'custom')
    OR provider LIKE 'custom:%'
  );
