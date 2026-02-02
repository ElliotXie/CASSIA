-- Utility function: auto-update updated_at on row changes
CREATE OR REPLACE FUNCTION public.handle_updated_at()
RETURNS trigger AS $$
BEGIN
    NEW.updated_at = now();
    RETURN NEW;
END;
$$ language plpgsql;

-- Utility function: encrypt API key (base64 encoding)
-- In production, use proper encryption with user-specific keys
CREATE OR REPLACE FUNCTION public.encrypt_api_key(api_key TEXT)
RETURNS TEXT AS $$
BEGIN
    RETURN encode(convert_to(api_key, 'UTF8'), 'base64');
END;
$$ language plpgsql;

-- Utility function: decrypt API key (base64 decoding)
CREATE OR REPLACE FUNCTION public.decrypt_api_key(encrypted_key TEXT)
RETURNS TEXT AS $$
BEGIN
    RETURN convert_from(decode(encrypted_key, 'base64'), 'UTF8');
END;
$$ language plpgsql;
