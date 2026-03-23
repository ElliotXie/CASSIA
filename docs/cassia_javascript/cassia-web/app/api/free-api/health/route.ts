import { NextResponse } from 'next/server';

/**
 * Health check endpoint for CASSIA Free API service
 * Used by R/Python packages to verify server availability
 */
export async function GET() {
    // Service temporarily paused
    return NextResponse.json({
        status: 'paused',
        message: 'Free API service is temporarily paused',
    }, { status: 503 });

    // Original response kept below for when service resumes
    return NextResponse.json({
        status: 'ok',
        service: 'cassia-free-api',
        version: '1.0.0',
        timestamp: new Date().toISOString(),
        limits: {
            max_free_clusters: 2,
            supported_providers: ['google', 'together', 'openrouter']
        }
    });
}
