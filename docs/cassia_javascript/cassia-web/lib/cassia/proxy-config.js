/**
 * Proxy configuration for CASSIA Web.
 *
 * Auto-detects China users (via timezone) and routes API calls through
 * a Cloudflare Worker reverse proxy. Users still BYOK.
 *
 * localStorage key "cassia_proxy" controls behavior:
 *   - unset           → auto-detect (timezone heuristic)
 *   - "off" / "none"  → force disabled, no auto-detect
 *   - "china"         → force enabled (use the china proxy)
 */

const PROXY_BASE = 'https://cassia-proxy.solodev66.workers.dev';

const PROXY_ENDPOINTS = {
  china: {
    openai: `${PROXY_BASE}/openai/v1`,
    anthropic: `${PROXY_BASE}/anthropic`,
    openrouter: `${PROXY_BASE}/openrouter/api/v1`,
  }
};

/**
 * Detect if user is likely in China based on browser timezone.
 */
function detectChina() {
  try {
    const tz = Intl.DateTimeFormat().resolvedOptions().timeZone || '';
    // Direct China timezone match
    if (tz.startsWith('Asia/Shanghai') || tz.startsWith('Asia/Chongqing') ||
        tz.startsWith('Asia/Harbin') || tz.startsWith('Asia/Urumqi') ||
        tz === 'PRC') {
      return true;
    }
    // Browser language as additional signal for UTC+8 zones
    const lang = navigator.language || '';
    if (lang.startsWith('zh') && tz.startsWith('Asia/')) {
      return true;
    }
  } catch (e) {
    // ignore
  }
  return false;
}

/**
 * Get the active proxy name, or null if disabled.
 *
 * EMERGENCY KILL-SWITCH (2026-04-07):
 * Auto-detection is currently DISABLED. The deployed Cloudflare worker lives at
 * `cassia-proxy.solodev66.workers.dev` and the `*.workers.dev` apex is DNS-poisoned
 * in mainland China, so routing China users through it just makes Anthropic /
 * OpenRouter (which would otherwise work directly) fail too. Until the worker is
 * rebound to a non-poisoned custom domain, this function only returns a proxy when
 * the user has *explicitly* opted in via `localStorage.cassia_proxy = "china"`.
 *
 * To re-enable auto-detection later: restore the `detectChina()` fallback below.
 */
export function getProxy() {
  if (typeof window !== 'undefined') {
    const override = localStorage.getItem('cassia_proxy');
    if (override === 'off' || override === 'none') return null;
    if (override && PROXY_ENDPOINTS[override]) return override;
  }
  // Auto-detection intentionally disabled — see comment above.
  return null;
}

/**
 * Set proxy override. Persists in localStorage.
 * @param {string|null} proxyName
 *   - null    → clear override, fall back to auto-detect
 *   - "off"   → force disabled (no proxy, even in China)
 *   - "china" → force enabled (route through china proxy)
 */
export function setProxy(proxyName) {
  if (typeof window === 'undefined') return;
  if (proxyName === null) {
    // Clear override so getProxy() falls back to auto-detect.
    localStorage.removeItem('cassia_proxy');
  } else if (proxyName === 'off') {
    localStorage.setItem('cassia_proxy', 'off');
  } else {
    localStorage.setItem('cassia_proxy', proxyName);
  }
}

/**
 * Get the base URL for a provider, using proxy if active.
 * @param {string} provider - "openai", "anthropic", or "openrouter"
 * @returns {{ baseUrl: string|null, isProxied: boolean }}
 */
export function getProviderBaseUrl(provider) {
  const proxy = getProxy();
  if (proxy && PROXY_ENDPOINTS[proxy]) {
    const url = PROXY_ENDPOINTS[proxy][provider.toLowerCase()];
    if (url) return { baseUrl: url, isProxied: true };
  }
  return { baseUrl: null, isProxied: false };
}
