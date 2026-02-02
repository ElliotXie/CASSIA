import { NextResponse } from 'next/server';

/**
 * Health check endpoint for CASSIA Free API service
 * Used by R/Python packages to verify server availability
 */
export async function GET() {
    return NextResponse.json({
        status: 'ok',
        service: 'cassia-free-api',
        version: '1.0.0',
        timestamp: new Date().toISOString(),
        limits: {
            max_jobs_per_day: 2,
            max_clusters_per_job: 30,
            supported_providers: ['google', 'together', 'openrouter']
        }
    });
}
