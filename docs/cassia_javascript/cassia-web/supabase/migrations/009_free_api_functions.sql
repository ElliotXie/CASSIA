-- Reset daily usage counters for free API keys
CREATE OR REPLACE FUNCTION public.reset_daily_key_usage()
RETURNS void AS $$
BEGIN
    UPDATE public.free_api_keys
    SET current_daily_usage = 0,
        last_reset_at = NOW()
    WHERE DATE(last_reset_at) < CURRENT_DATE;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Get next available key with round-robin and auto-reset
-- Returns the least-used active key that hasn't hit its daily limit
CREATE OR REPLACE FUNCTION public.get_next_free_key(p_provider VARCHAR)
RETURNS TABLE(key_id UUID, decrypted_key TEXT) AS $$
DECLARE
    v_key RECORD;
BEGIN
    -- Reset daily usage if needed
    PERFORM public.reset_daily_key_usage();

    -- Get the least-used active key that hasn't hit its limit
    SELECT k.id, public.decrypt_api_key(k.encrypted_key) as dkey
    INTO v_key
    FROM public.free_api_keys k
    WHERE k.provider = p_provider
      AND k.is_active = true
      AND k.current_daily_usage < k.daily_limit
    ORDER BY k.current_daily_usage ASC, k.priority DESC
    LIMIT 1;

    IF v_key.id IS NOT NULL THEN
        -- Increment usage counter
        UPDATE public.free_api_keys
        SET current_daily_usage = current_daily_usage + 1
        WHERE id = v_key.id;

        key_id := v_key.id;
        decrypted_key := v_key.dkey;
        RETURN NEXT;
    END IF;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
