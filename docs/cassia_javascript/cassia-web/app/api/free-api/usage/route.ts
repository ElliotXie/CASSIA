import { NextRequest, NextResponse } from 'next/server';
import { createClient } from '@supabase/supabase-js';

const MAX_JOBS_PER_DAY = 2;
const MAX_CLUSTERS_PER_JOB = 30;

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL;
const supabaseServiceKey = process.env.SUPABASE_SERVICE_ROLE_KEY;

/**
 * Get usage statistics for a machine
 * Used by R/Python packages to check remaining quota
 */
export async function GET(request: NextRequest) {
    // Check Supabase configuration
    if (!supabaseUrl || !supabaseServiceKey) {
        return NextResponse.json(
            { error: 'Service not configured' },
            { status: 503 }
        );
    }

    const supabaseAdmin = createClient(supabaseUrl, supabaseServiceKey);

    try {
        const { searchParams } = new URL(request.url);
        const machine_id = searchParams.get('machine_id');

        if (!machine_id) {
            return NextResponse.json(
                { error: 'machine_id query parameter is required' },
                { status: 400 }
            );
        }

        // Validate machine_id format
        if (!/^[a-f0-9]{32}$/i.test(machine_id)) {
            return NextResponse.json(
                { error: 'Invalid machine_id format' },
                { status: 400 }
            );
        }

        const today = new Date().toISOString().split('T')[0];

        const { data: usage, error } = await supabaseAdmin
            .from('free_api_usage')
            .select('provider, job_count, cluster_count, last_job_at')
            .eq('machine_id', machine_id)
            .eq('usage_date', today);

        if (error) {
            console.error('[usage] Failed to get usage:', error);
            return NextResponse.json({ providers: {} });
        }

        // Build response with usage per provider
        const providers: Record<string, {
            jobs_used: number;
            jobs_remaining: number;
            clusters_used: number;
            last_job_at: string | null;
        }> = {};

        for (const row of usage || []) {
            providers[row.provider] = {
                jobs_used: row.job_count,
                jobs_remaining: MAX_JOBS_PER_DAY - row.job_count,
                clusters_used: row.cluster_count,
                last_job_at: row.last_job_at
            };
        }

        // Add providers with no usage
        for (const provider of ['google', 'together', 'openrouter']) {
            if (!providers[provider]) {
                providers[provider] = {
                    jobs_used: 0,
                    jobs_remaining: MAX_JOBS_PER_DAY,
                    clusters_used: 0,
                    last_job_at: null
                };
            }
        }

        return NextResponse.json({
            machine_id: machine_id.substring(0, 8) + '...',  // Truncate for privacy
            date: today,
            providers,
            limits: {
                max_jobs_per_day: MAX_JOBS_PER_DAY,
                max_clusters_per_job: MAX_CLUSTERS_PER_JOB
            }
        });

    } catch (error) {
        console.error('[usage] Error:', error);
        return NextResponse.json(
            { error: 'Internal server error' },
            { status: 500 }
        );
    }
}
