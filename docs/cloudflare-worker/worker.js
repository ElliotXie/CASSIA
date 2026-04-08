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

const CORS_HEADERS = {
  'Access-Control-Allow-Origin': '*',
  'Access-Control-Allow-Methods': 'GET, POST, PUT, DELETE, OPTIONS',
  'Access-Control-Allow-Headers': '*',
  'Access-Control-Max-Age': '86400',
};

// Headers that should NOT be forwarded to upstream:
// - browser-only headers that confuse upstream APIs (Anthropic in particular
//   treats requests with `Origin` set as browser-originated and will reject
//   them unless `anthropic-dangerous-direct-browser-access: true` is set)
// - hop-by-hop headers (per RFC 7230) which must not be proxied
// - `host` is set automatically by the runtime from the target URL
const STRIP_REQUEST_HEADERS = new Set([
  'host',
  'origin',
  'referer',
  'cookie',
  'connection',
  'keep-alive',
  'proxy-authorization',
  'proxy-connection',
  'te',
  'trailer',
  'transfer-encoding',
  'upgrade',
]);

function filterRequestHeaders(srcHeaders) {
  const out = new Headers();
  for (const [key, value] of srcHeaders.entries()) {
    const k = key.toLowerCase();
    if (STRIP_REQUEST_HEADERS.has(k)) continue;
    if (k.startsWith('cf-')) continue;          // Cloudflare-internal
    if (k.startsWith('sec-')) continue;         // browser fetch metadata
    if (k.startsWith('x-forwarded-')) continue; // proxy chain metadata
    out.set(key, value);
  }
  return out;
}

export default {
  async fetch(request) {
    // Handle CORS preflight FIRST — never forward OPTIONS to upstream.
    // (Previously this check happened *after* the upstream fetch, which made
    // every preflight pay a real round-trip to api.openai.com / api.anthropic.com
    // / openrouter.ai and could fail if the upstream rejected OPTIONS.)
    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: CORS_HEADERS });
    }

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
        {
          status: 404,
          headers: { 'Content-Type': 'application/json', ...CORS_HEADERS },
        }
      );
    }

    // Build upstream URL: strip the prefix and forward the rest
    const targetPath = url.pathname.slice(prefix.length);
    const targetUrl = upstream + targetPath + url.search;

    // Forward the request, dropping browser-only / hop-by-hop headers so
    // upstream APIs see what looks like a normal server-to-server call.
    const headers = filterRequestHeaders(request.headers);

    const upstreamResponse = await fetch(targetUrl, {
      method: request.method,
      headers,
      body: request.body,
    });

    // Mirror upstream response, but force permissive CORS headers so the
    // browser will accept the response cross-origin.
    const responseHeaders = new Headers(upstreamResponse.headers);
    for (const [k, v] of Object.entries(CORS_HEADERS)) {
      responseHeaders.set(k, v);
    }

    return new Response(upstreamResponse.body, {
      status: upstreamResponse.status,
      statusText: upstreamResponse.statusText,
      headers: responseHeaders,
    });
  },
};
