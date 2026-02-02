import { NextRequest, NextResponse } from 'next/server';
import { createClient } from '@supabase/supabase-js';

// Lifetime limit: 2 free cluster annotations per machine, ever
const MAX_FREE_CLUSTERS = 2;

// Initialize Supabase client with service role (server-side only)
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL;
const supabaseServiceKey = process.env.SUPABASE_SERVICE_ROLE_KEY;

interface GetKeyRequest {
    machine_id: string;
    provider: 'google' | 'together' | 'openrouter';
    num_clusters: number;
}

export async function POST(request: NextRequest) {
    // Check Supabase configuration
    if (!supabaseUrl || !supabaseServiceKey) {
        console.error('[get-key] Supabase not configured');
        return NextResponse.json(
            { error: 'Service not configured' },
            { status: 503 }
        );
    }

    const supabaseAdmin = createClient(supabaseUrl, supabaseServiceKey);

    try {
        const body: GetKeyRequest = await request.json();
        const { machine_id, provider, num_clusters } = body;

        // Validate request
        if (!machine_id || !provider || num_clusters === undefined) {
            return NextResponse.json(
                { error: 'Missing required fields: machine_id, provider, num_clusters' },
                { status: 400 }
            );
        }

        // Validate machine_id format (should be 32 char hex)
        if (!/^[a-f0-9]{32}$/i.test(machine_id)) {
            return NextResponse.json(
                { error: 'Invalid machine_id format' },
                { status: 400 }
            );
        }

        if (!['google', 'together', 'openrouter'].includes(provider)) {
            return NextResponse.json(
                { error: 'Invalid provider. Must be: google, together, or openrouter' },
                { status: 400 }
            );
        }

        if (num_clusters < 1) {
            return NextResponse.json(
                { error: 'num_clusters must be at least 1' },
                { status: 400 }
            );
        }

        // Check lifetime usage
        const { data: usageData } = await supabaseAdmin
            .from('free_api_usage')
            .select('clusters_used')
            .eq('machine_id', machine_id)
            .single();

        const clustersUsed = usageData?.clusters_used || 0;
        const remaining = MAX_FREE_CLUSTERS - clustersUsed;

        if (remaining <= 0) {
            return NextResponse.json(
                {
                    error: `Free API limit reached. You have used all ${MAX_FREE_CLUSTERS} free cluster annotations. Set your own API key: CASSIA.set_api_key('provider', 'your-key')`,
                    limit: MAX_FREE_CLUSTERS,
                    used: clustersUsed,
                    remaining: 0
                },
                { status: 429 }
            );
        }

        if (num_clusters > remaining) {
            return NextResponse.json(
                {
                    error: `Only ${remaining} free cluster annotation(s) remaining. Requested ${num_clusters}.`,
                    limit: MAX_FREE_CLUSTERS,
                    used: clustersUsed,
                    remaining: remaining
                },
                { status: 429 }
            );
        }

        // Get next available key using round-robin
        const { data: keyData, error: keyError } = await supabaseAdmin
            .rpc('get_next_free_key', { p_provider: provider });

        if (keyError || !keyData || keyData.length === 0) {
            console.error('[get-key] Failed to get free key:', keyError);
            return NextResponse.json(
                { error: 'No API keys available for this provider. Please try again later or use a different provider.' },
                { status: 503 }
            );
        }

        const { decrypted_key } = keyData[0];

        // Update lifetime usage tracking (upsert)
        const newClustersUsed = clustersUsed + num_clusters;
        const { error: upsertError } = await supabaseAdmin
            .from('free_api_usage')
            .upsert({
                machine_id,
                clusters_used: newClustersUsed,
                last_used_at: new Date().toISOString()
            }, {
                onConflict: 'machine_id'
            });

        if (upsertError) {
            console.error('[get-key] Failed to update usage:', upsertError);
            // Continue anyway - key distribution is more important
        }

        // Log for monitoring (no sensitive data)
        console.log(`[get-key] Distributed ${provider} key to ${machine_id.substring(0, 8)}... (${num_clusters} clusters, ${newClustersUsed}/${MAX_FREE_CLUSTERS} lifetime used)`);

        // Return the API key
        return NextResponse.json({
            api_key: decrypted_key,
            provider,
            remaining_clusters: remaining - num_clusters,
            limit: MAX_FREE_CLUSTERS
        });

    } catch (error) {
        console.error('[get-key] Error:', error);
        return NextResponse.json(
            { error: 'Internal server error' },
            { status: 500 }
        );
    }
}
