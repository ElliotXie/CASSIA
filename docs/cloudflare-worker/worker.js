/**
 * CASSIA API Proxy - Cloudflare Worker
 *
 * Reverse proxy for users in regions where US API providers are blocked.
 * Routes requests based on path prefix to the appropriate upstream API.
 *
 * Deploy: wrangler deploy
 * Bind to: proxy.cassia.bio
 */

const UPSTREAM_MAP = {
  '/openai/':      'https://api.openai.com/',
  '/anthropic/':   'https://api.anthropic.com/',
  '/openrouter/':  'https://openrouter.ai/',
};

export default {
  async fetch(request) {
    const url = new URL(request.url);

    // Find matching upstream
    let upstream = null;
    let prefix = null;
    for (const [p, u] of Object.entries(UPSTREAM_MAP)) {
      if (url.pathname.startsWith(p)) {
        upstream = u;
        prefix = p;
        break;
      }
    }

    if (!upstream) {
      return new Response(
        JSON.stringify({
          error: 'Unknown route. Use /openai/, /anthropic/, or /openrouter/ prefix.',
          usage: {
            openai: 'proxy.cassia.bio/openai/v1/chat/completions',
            anthropic: 'proxy.cassia.bio/anthropic/v1/messages',
            openrouter: 'proxy.cassia.bio/openrouter/api/v1/chat/completions',
          }
        }),
        { status: 404, headers: { 'Content-Type': 'application/json' } }
      );
    }

    // Build upstream URL: strip the prefix and forward the rest
    const targetPath = url.pathname.slice(prefix.length);
    const targetUrl = upstream + targetPath + url.search;

    // Forward the request with original headers and body
    const headers = new Headers(request.headers);
    headers.set('Host', new URL(upstream).host);

    const response = await fetch(targetUrl, {
      method: request.method,
      headers: headers,
      body: request.body,
    });

    // Return response with CORS headers for web usage
    const responseHeaders = new Headers(response.headers);
    responseHeaders.set('Access-Control-Allow-Origin', '*');
    responseHeaders.set('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS');
    responseHeaders.set('Access-Control-Allow-Headers', '*');

    // Handle preflight
    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: responseHeaders });
    }

    return new Response(response.body, {
      status: response.status,
      statusText: response.statusText,
      headers: responseHeaders,
    });
  },
};
