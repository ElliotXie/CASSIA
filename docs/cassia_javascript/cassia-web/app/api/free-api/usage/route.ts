import { NextRequest, NextResponse } from 'next/server';
import { createClient } from '@supabase/supabase-js';

const MAX_FREE_CLUSTERS = 2;

const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL;
const supabaseServiceKey = process.env.SUPABASE_SERVICE_ROLE_KEY;

/**
 * Get lifetime usage statistics for a machine
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

        const { data: usage, error } = await supabaseAdmin
            .from('free_api_usage')
            .select('clusters_used, last_used_at')
            .eq('machine_id', machine_id)
            .single();

        if (error && error.code !== 'PGRST116') {
            // PGRST116 = no rows found (new user), which is fine
            console.error('[usage] Failed to get usage:', error);
        }

        const clustersUsed = usage?.clusters_used || 0;

        return NextResponse.json({
            machine_id: machine_id.substring(0, 8) + '...',
            clusters_used: clustersUsed,
            clusters_remaining: Math.max(0, MAX_FREE_CLUSTERS - clustersUsed),
            last_used_at: usage?.last_used_at || null,
            limits: {
                max_free_clusters: MAX_FREE_CLUSTERS
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
