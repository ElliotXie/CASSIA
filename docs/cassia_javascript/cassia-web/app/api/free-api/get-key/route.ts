import { NextRequest, NextResponse } from 'next/server';
import { createClient } from '@supabase/supabase-js';

// Constants
const MAX_JOBS_PER_DAY = 2;
const MAX_CLUSTERS_PER_JOB = 30;

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

        if (num_clusters > MAX_CLUSTERS_PER_JOB) {
            return NextResponse.json(
                {
                    error: `Cluster limit exceeded. Maximum ${MAX_CLUSTERS_PER_JOB} clusters per job.`,
                    limit: MAX_CLUSTERS_PER_JOB,
                    requested: num_clusters
                },
                { status: 400 }
            );
        }

        if (num_clusters < 1) {
            return NextResponse.json(
                { error: 'num_clusters must be at least 1' },
                { status: 400 }
            );
        }

        // Check rate limit
        const today = new Date().toISOString().split('T')[0];
        const { data: usageData } = await supabaseAdmin
            .from('free_api_usage')
            .select('job_count, cluster_count')
            .eq('machine_id', machine_id)
            .eq('provider', provider)
            .eq('usage_date', today)
            .single();

        const currentJobCount = usageData?.job_count || 0;

        if (currentJobCount >= MAX_JOBS_PER_DAY) {
            // Calculate time until reset (midnight UTC)
            const now = new Date();
            const tomorrow = new Date(now);
            tomorrow.setUTCDate(tomorrow.getUTCDate() + 1);
            tomorrow.setUTCHours(0, 0, 0, 0);
            const hoursUntilReset = Math.ceil((tomorrow.getTime() - now.getTime()) / (1000 * 60 * 60));

            return NextResponse.json(
                {
                    error: `Daily limit reached. Maximum ${MAX_JOBS_PER_DAY} jobs per day. Resets in ~${hoursUntilReset} hours.`,
                    limit: MAX_JOBS_PER_DAY,
                    used: currentJobCount,
                    reset_at: tomorrow.toISOString()
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

        // Update usage tracking (upsert)
        const { error: upsertError } = await supabaseAdmin
            .from('free_api_usage')
            .upsert({
                machine_id,
                provider,
                usage_date: today,
                job_count: currentJobCount + 1,
                cluster_count: (usageData?.cluster_count || 0) + num_clusters,
                last_job_at: new Date().toISOString()
            }, {
                onConflict: 'machine_id,provider,usage_date'
            });

        if (upsertError) {
            console.error('[get-key] Failed to update usage:', upsertError);
            // Continue anyway - key distribution is more important
        }

        // Log for monitoring (no sensitive data)
        console.log(`[get-key] Distributed ${provider} key to ${machine_id.substring(0, 8)}... (${num_clusters} clusters, job ${currentJobCount + 1}/${MAX_JOBS_PER_DAY})`);

        // Return the API key
        return NextResponse.json({
            api_key: decrypted_key,
            provider,
            remaining_jobs: MAX_JOBS_PER_DAY - currentJobCount - 1,
            cluster_limit: MAX_CLUSTERS_PER_JOB
        });

    } catch (error) {
        console.error('[get-key] Error:', error);
        return NextResponse.json(
            { error: 'Internal server error' },
            { status: 500 }
        );
    }
}
