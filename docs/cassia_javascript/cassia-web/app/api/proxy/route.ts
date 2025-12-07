import { NextRequest, NextResponse } from 'next/server';

/**
 * Proxy API route for custom OpenAI-compatible endpoints
 * This bypasses CORS restrictions by making server-side requests
 */
export async function POST(request: NextRequest) {
    try {
        const body = await request.json();
        const { baseUrl, apiKey, ...openaiRequest } = body;

        if (!baseUrl || !apiKey) {
            return NextResponse.json(
                { error: 'baseUrl and apiKey are required' },
                { status: 400 }
            );
        }

        // Construct the full URL for chat completions
        const url = baseUrl.endsWith('/v1')
            ? `${baseUrl}/chat/completions`
            : `${baseUrl}/v1/chat/completions`;

        console.log(`[Proxy] Forwarding request to: ${url}`);

        const response = await fetch(url, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'Authorization': `Bearer ${apiKey}`,
            },
            body: JSON.stringify(openaiRequest),
        });

        if (!response.ok) {
            const errorText = await response.text();
            console.error(`[Proxy] Error from upstream: ${response.status} - ${errorText}`);
            return NextResponse.json(
                { error: `Upstream error: ${response.status}`, details: errorText },
                { status: response.status }
            );
        }

        const data = await response.json();
        return NextResponse.json(data);

    } catch (error) {
        console.error('[Proxy] Error:', error);
        return NextResponse.json(
            { error: error instanceof Error ? error.message : 'Proxy request failed' },
            { status: 500 }
        );
    }
}
